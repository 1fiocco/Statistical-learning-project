# Section 3 (Paper) — Algorithm 2: Projection-based Conformal Bands

suppressPackageStartupMessages({
  library(data.table)
  library(fda)
  library(mvtnorm)
  library(mclust)
})

# ============================================================
# NEW (PLOTTING ONLY): X-axis like your screenshot
# - major ticks: 0,100,...,600 with labels "HH:MM\n(100)"
# - minor ticks + vertical guides each 60 min (1-hour bands)
# ============================================================

.make_time_labels_from_hhmm <- function(start_hhmm, ticks_min) {
  parts <- strsplit(start_hhmm, ":", fixed = TRUE)[[1]]
  if (length(parts) != 2) stop("start_hhmm must be 'HH:MM'")
  start_min <- as.integer(parts[1]) * 60 + as.integer(parts[2])
  
  labs_min <- (start_min + ticks_min) %% (24 * 60)
  sprintf("%02d:%02d", labs_min %/% 60, labs_min %% 60)
}

add_time_axis <- function(start_hhmm,
                          x0 = 0,
                          total_min = 600,
                          major_by = 100,
                          minor_by = 60,
                          cex.axis = 0.9) {
  major_ticks <- seq(0, total_min, by = major_by)
  minor_ticks <- seq(0, total_min, by = minor_by)
  
  # minor ticks (hourly) without labels
  axis(1, at = x0 + minor_ticks, labels = FALSE, tcl = -0.25)
  
  # major labels like screenshot
  times <- .make_time_labels_from_hhmm(start_hhmm, major_ticks)
  labs  <- sprintf("%s\n(%d)", times, major_ticks)
  axis(1, at = x0 + major_ticks, labels = labs, cex.axis = cex.axis)
  
  invisible(list(major = major_ticks, minor = minor_ticks))
}

add_time_axis_dual <- function(start_day_hhmm,
                               start_night_hhmm,
                               x0 = 0,
                               total_min = 600,
                               major_by = 100,
                               minor_by = 60,
                               cex.axis = 0.9) {
  major_ticks <- seq(0, total_min, by = major_by)
  minor_ticks <- seq(0, total_min, by = minor_by)
  
  axis(1, at = x0 + minor_ticks, labels = FALSE, tcl = -0.25)
  
  td <- .make_time_labels_from_hhmm(start_day_hhmm, major_ticks)
  tn <- .make_time_labels_from_hhmm(start_night_hhmm, major_ticks)
  labs <- sprintf("%s | %s\n(%d)", td, tn, major_ticks)
  axis(1, at = x0 + major_ticks, labels = labs, cex.axis = cex.axis)
  
  invisible(list(major = major_ticks, minor = minor_ticks))
}

add_hour_guides <- function(x0 = 0, total_min = 600, by = 60) {
  v <- x0 + seq(0, total_min, by = by)
  abline(v = v, lty = 3, col = "grey85")
  invisible(v)
}

# ---------------------------
# Helpers
# ---------------------------

subset_fd <- function(fdobj, idx) {
  fd(coef = fdobj$coef[, idx, drop = FALSE], basisobj = fdobj$basis)
}

fpca_scores_for_all <- function(fd_all, meanfd_train, harmfd_train) {
  # raw projection scores: ncurves x p
  scores <- inprod(fd_all, harmfd_train)
  
  # subtract projection of the mean (safe even if it's ~0)
  mean_scores <- as.numeric(inprod(meanfd_train, harmfd_train))  # length p
  scores <- sweep(scores, 2, mean_scores, "-")
  
  scores
}

regularize_sigma <- function(S, eps = 1e-6) {
  # ensure PD
  S + diag(eps, nrow(S))
}

extract_gmm_params_vvv <- function(mclust_fit) {
  if (mclust_fit$modelName != "VVV") {
    stop("Expected modelName='VVV' but got: ", mclust_fit$modelName,
         "\nTip: call Mclust(..., modelNames='VVV').")
  }
  G <- mclust_fit$G
  means <- mclust_fit$parameters$mean               # p x G
  if (is.null(dim(means))) means <- matrix(means, ncol = G)
  sigmas <- mclust_fit$parameters$variance$sigma    # p x p x G
  pis <- mclust_fit$parameters$pro                  # length G
  list(G = G, means = means, sigmas = sigmas, pis = pis)
}

max_component_score <- function(scores_mat, params) {
  # score_i = max_k pi_k * dmvnorm(score_i; mu_k, Sigma_k)
  n <- nrow(scores_mat)
  G <- params$G
  
  contrib <- matrix(NA_real_, nrow = n, ncol = G)
  for (k in 1:G) {
    mu <- params$means[, k]
    S  <- regularize_sigma(params$sigmas[, , k, drop = FALSE][,,1])
    dk <- mvtnorm::dmvnorm(scores_mat, mean = mu, sigma = S, log = FALSE)
    contrib[, k] <- params$pis[k] * dk
  }
  
  kstar <- max.col(contrib, ties.method = "first")
  score <- contrib[cbind(1:n, kstar)]
  list(score = score, kstar = kstar, contrib = contrib)
}

conformal_threshold <- function(cal_scores, alpha) {
  # Split conformal threshold for CONFORMITY scores (higher = more conforming)
  # We want Tn = { score >= lambda }, so lambda is a lower-tail cutoff.
  n2 <- length(cal_scores)
  s_sorted <- sort(cal_scores)
  
  k <- ceiling((n2 + 1) * alpha) - 1
  k <- max(1L, min(k, n2))  # R is 1-indexed
  
  lambda <- s_sorted[k]
  list(lambda = lambda, k = k, n_cal = n2)
}

ellipsoid_radii_from_lambda <- function(params, lambda) {
  # Condition per component:
  #   pi_k * phi(x; mu_k, Sigma_k) >= lambda
  # => Mahalanobis^2 <= r_k^2
  G <- params$G
  p <- nrow(params$means)
  
  r2 <- rep(NA_real_, G)
  keep <- rep(FALSE, G)
  
  for (k in 1:G) {
    pi_k <- params$pis[k]
    S    <- regularize_sigma(params$sigmas[, , k, drop = FALSE][,,1])
    
    log_det <- as.numeric(determinant(S, logarithm = TRUE)$modulus)
    log_phi_const <- -(p / 2) * log(2*pi) - 0.5 * log_det
    
    val <- 2 * (log(pi_k) + log_phi_const - log(lambda))
    
    if (is.finite(val) && val > 0) {
      r2[k] <- val
      keep[k] <- TRUE
    }
  }
  
  list(r2 = r2, r = sqrt(r2), keep = keep)
}

# ============================================================
# BUILD BANDS  (MODIFIED: also compute CENTRAL FUNCTIONS)
# ============================================================

build_component_bands <- function(pca_train, params, radii, tgrid) {
  # pca_train$harmonics: fd object of eigenfunctions (p harmonics)
  # pca_train$meanfd: mean function (fd)
  phi_eval <- eval.fd(tgrid, pca_train$harmonics)          # L x p
  mu_eval  <- as.numeric(eval.fd(tgrid, pca_train$meanfd)) # L
  
  G <- params$G
  p <- ncol(phi_eval)
  
  comp_bands   <- vector("list", G)
  comp_centers <- vector("list", G)
  
  # Mixture center in score space: sum_k pi_k * mu_k
  mu_mix_scores <- as.numeric(params$means %*% params$pis)   # length p
  center_mix <- mu_eval + as.numeric(phi_eval %*% mu_mix_scores)
  
  for (k in 1:G) {
    mu_k <- params$means[, k]
    
    # Center curve for component k: c_k(t) = mu(t) + phi(t)^T mu_k
    center_k <- mu_eval + as.numeric(phi_eval %*% mu_k)
    comp_centers[[k]] <- center_k
    
    if (!radii$keep[k]) {
      comp_bands[[k]] <- NULL
      next
    }
    
    S_k <- regularize_sigma(params$sigmas[, , k, drop = FALSE][,,1])
    r_k <- radii$r[k]
    
    # compute sqrt(phi(t)^T S phi(t)) per t
    tmp <- phi_eval %*% S_k
    at_S_a <- rowSums(tmp * phi_eval)
    at_S_a[at_S_a < 0] <- 0
    sd_line <- sqrt(at_S_a)
    
    lower <- center_k - r_k * sd_line
    upper <- center_k + r_k * sd_line
    
    comp_bands[[k]] <- list(lower = lower, upper = upper)
  }
  
  # Envelope (outer band)
  lowers <- do.call(cbind, lapply(comp_bands, function(b) if (is.null(b)) NULL else b$lower))
  uppers <- do.call(cbind, lapply(comp_bands, function(b) if (is.null(b)) NULL else b$upper))
  
  env_lower <- apply(lowers, 1, min, na.rm = TRUE)
  env_upper <- apply(uppers, 1, max, na.rm = TRUE)
  
  list(
    component = comp_bands,
    envelope  = list(lower = env_lower, upper = env_upper),
    center    = list(
      meanfd    = mu_eval,        # mu(t) (train mean)
      mixture   = center_mix,     # overall mixture center
      component = comp_centers    # list of c_k(t)
    ),
    phi_eval  = phi_eval,
    mu_eval   = mu_eval
  )
}

# ============================================================
# Plot bands (MODIFIED: time axis like screenshot)
# ============================================================

plot_bands <- function(tgrid, Xmat, band_obj,
                       main = "",
                       n_show = 20L,
                       use_envelope = TRUE,
                       draw_center = c("none", "mixture", "meanfd"),
                       draw_component_centers = TRUE,
                       lwd_band = 2,
                       lwd_center = 3,
                       col_center = "black",
                       lty_center = 3,
                       # NEW plotting-only:
                       start_hhmm = NULL,
                       total_min = 600,
                       major_by = 100,
                       minor_by = 60,
                       add_hour_lines = TRUE,
                       cex.axis = 0.9) {
  
  draw_center <- match.arg(draw_center)
  
  n <- nrow(Xmat)
  idx <- if (n <= n_show) seq_len(n) else sample.int(n, n_show)
  
  y_rng <- range(Xmat[idx, ],
                 band_obj$envelope$lower, band_obj$envelope$upper,
                 band_obj$center$mixture, band_obj$center$meanfd,
                 na.rm = TRUE)
  
  x0 <- min(tgrid)
  xlim_use <- c(x0, x0 + total_min)
  
  plot(NA, xlim = xlim_use, ylim = y_rng,
       xlab = "time of day (HH:MM) / (minute index)",
       ylab = "HR", main = main, xaxt = "n")
  grid()
  if (isTRUE(add_hour_lines)) add_hour_guides(x0 = x0, total_min = total_min, by = minor_by)
  
  # observed curves
  for (i in idx) lines(tgrid, Xmat[i, ], col = "grey70")
  
  # bands
  if (use_envelope) {
    lines(tgrid, band_obj$envelope$lower, lwd = lwd_band)
    lines(tgrid, band_obj$envelope$upper, lwd = lwd_band)
  } else {
    for (k in seq_along(band_obj$component)) {
      b <- band_obj$component[[k]]
      if (is.null(b)) next
      lines(tgrid, b$lower, lwd = lwd_band)
      lines(tgrid, b$upper, lwd = lwd_band)
    }
    if (draw_component_centers) {
      for (k in seq_along(band_obj$center$component)) {
        ck <- band_obj$center$component[[k]]
        if (is.null(ck)) next
        lines(tgrid, ck, lwd = 2, lty = 3, col = "grey30")
      }
    }
  }
  
  # central function (optional)
  if (draw_center == "mixture") {
    lines(tgrid, band_obj$center$mixture, lwd = lwd_center, col = col_center, lty = lty_center)
  } else if (draw_center == "meanfd") {
    lines(tgrid, band_obj$center$meanfd, lwd = lwd_center, col = col_center, lty = lty_center)
  }
  
  # NEW axis
  if (!is.null(start_hhmm)) {
    add_time_axis(start_hhmm, x0 = x0, total_min = total_min,
                  major_by = major_by, minor_by = minor_by, cex.axis = cex.axis)
  } else {
    axis(1)
  }
}

# ---------------------------
# Main runner for ONE period
# ---------------------------

run_section3_one_period <- function(obj_rds_path,
                                    period_name = "day",
                                    alpha = 0.10,
                                    p = 6,
                                    n1_frac = 0.70,
                                    seed = 123,
                                    G_max = 5,
                                    out_dir = ".",
                                    prefix = NULL) {
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  if (is.null(prefix)) prefix <- period_name
  
  obj <- readRDS(obj_rds_path)
  
  # Expecting: X (n x 600), continuous_functions (fd), tgrid, meta
  X <- obj$X
  fd_all <- obj$continuous_functions
  meta <- obj$meta
  
  n <- nrow(X)
  n1 <- floor(n1_frac * n)
  if (n1 < 10) stop("Too few training curves. Increase n1_frac or check dataset.")
  n2 <- n - n1
  
  set.seed(seed)
  idx_train <- sample.int(n, n1)
  idx_cal   <- setdiff(seq_len(n), idx_train)
  
  fd_train <- subset_fd(fd_all, idx_train)
  
  # FPCA basis on train set
  pca_train <- pca.fd(fd_train, nharm = p, centerfns = TRUE)
  
  # Scores for all curves w.r.t train mean/harmonics
  scores_all <- fpca_scores_for_all(fd_all, pca_train$meanfd, pca_train$harmonics)
  scores_train <- scores_all[idx_train, , drop = FALSE]
  scores_cal   <- scores_all[idx_cal, , drop = FALSE]
  
  # ------------------------------------------------------------
  # GMM selection via BIC (keep ONLY successful VVV fits)
  # ------------------------------------------------------------
  fits <- vector("list", G_max)
  bic_vals <- rep(NA_real_, G_max)
  model_name <- rep(NA_character_, G_max)
  
  for (g in 1:G_max) {
    fit_g <- tryCatch(
      mclust::Mclust(scores_train, G = g, modelNames = "VVV", verbose = FALSE),
      error = function(e) NULL
    )
    
    # keep only if it's a proper object and truly VVV
    if (!is.null(fit_g) && !is.null(fit_g$modelName) && fit_g$modelName == "VVV") {
      fits[[g]] <- fit_g
      bic_vals[g] <- fit_g$bic
      model_name[g] <- fit_g$modelName
    } else {
      fits[[g]] <- NULL
      bic_vals[g] <- NA_real_
      model_name[g] <- if (!is.null(fit_g$modelName)) fit_g$modelName else NA_character_
    }
  }
  
  bic_table <- data.frame(G = 1:G_max, model = model_name, BIC = bic_vals)
  message("[", period_name, "] BIC values (requested=VVV, kept only VVV):")
  print(bic_table)
  
  best_idx <- which.max(ifelse(is.na(bic_vals), -Inf, bic_vals))
  if (!is.finite(bic_vals[best_idx])) {
    stop(
      "No valid VVV fit succeeded for G=1..", G_max, ".\n",
      "Try one (or more) of:\n",
      "  - reduce p (nharm) (e.g., p=3 or p=4)\n",
      "  - reduce G_max (e.g., 2 or 3)\n",
      "  - increase regularization eps in regularize_sigma()\n",
      "  - check for near-constant score dimensions (var ~ 0)\n"
    )
  }
  
  gmm <- fits[[best_idx]]
  message("[", period_name, "] Selected G = ", gmm$G, " (BIC=", round(gmm$bic, 3), ")")
  params <- extract_gmm_params_vvv(gmm)
  
  # conformity scores on calibration set: max_k pi_k * N(...)
  cal_sc <- max_component_score(scores_cal, params)$score
  
  thr <- conformal_threshold(cal_sc, alpha = alpha)
  lambda <- thr$lambda
  
  radii <- ellipsoid_radii_from_lambda(params, lambda)
  
  # Build band(s) + centers
  tgrid <- obj$tgrid
  band_obj <- build_component_bands(pca_train, params, radii, tgrid)
  
  # Save objects
  out <- list(
    period = period_name,
    alpha = alpha,
    p = p,
    n = n, n1 = n1, n2 = n2,
    idx_train = idx_train,
    idx_cal = idx_cal,
    pca_train = pca_train,
    gmm = gmm,
    params = params,
    lambda = lambda,
    radii = radii,
    band = band_obj
  )
  
  saveRDS(out, file = file.path(out_dir, paste0("section3_fit_", prefix, ".rds")))
  
  inside_cal <- mean(cal_sc >= lambda)
  message("[", period_name, "] Calibration inside-set fraction: ", round(inside_cal, 3),
          " (target ≈ ", 1 - alpha, ")")
  
  # Plot to PNG (union of component bands) + centers
  # NEW (plotting-only): choose start time for axis
  start_hhmm_plot <- NULL
  pn <- tolower(period_name)
  if (pn %in% c("day", "day_log")) start_hhmm_plot <- "09:30"
  if (pn %in% c("night", "night_log")) start_hhmm_plot <- "22:00"
  
  png(file.path(out_dir, paste0("section3_band_", prefix, ".png")),
      width = 1400, height = 900, res = 150)
  
  plot_bands(
    tgrid = tgrid,
    Xmat = X,
    band_obj = band_obj,
    main = paste0("Section 3 band (", period_name, "), alpha=", alpha,
                  ", p=", p, ", G=", params$G),
    n_show = 25,
    use_envelope = FALSE,
    draw_center = "mixture",          # << central curve
    draw_component_centers = TRUE,    # << component centers
    start_hhmm = start_hhmm_plot,     # << ONLY axis
    total_min = 600,
    major_by = 100,
    minor_by = 60,
    add_hour_lines = TRUE,
    cex.axis = 0.9
  )
  
  dev.off()
  
  invisible(out)
}

# ==========================
# RUN: day + night
# ==========================
# to be modified (INPUT) -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
day_obj_rds   <- "C:\\Users\\2003l\\OneDrive\\Documenti\\Minor\\SDS_R\\Project functional\\heart_rate\\outputs_rds\\functional_object_day.rds"
night_obj_rds <- "C:\\Users\\2003l\\OneDrive\\Documenti\\Minor\\SDS_R\\Project functional\\heart_rate\\outputs_rds\\functional_object_night.rds"
out_dir       <- "C:\\Users\\2003l\\OneDrive\\Documenti\\Minor\\SDS_R\\Project functional\\heart_rate\\outputs_section3"

fit_day <- run_section3_one_period(
  obj_rds_path = day_obj_rds,
  period_name = "day",
  alpha = 0.10,
  p = 6,
  n1_frac = 0.70,
  seed = 123,
  G_max = 4,
  out_dir = out_dir,
  prefix = "day"
)

fit_night <- run_section3_one_period(
  obj_rds_path = night_obj_rds,
  period_name = "night",
  alpha = 0.10,
  p = 6,
  n1_frac = 0.70,
  seed = 123,
  G_max = 4,
  out_dir = out_dir,
  prefix = "night"
)

message("Done. Files written to: ", out_dir)

# ============================================================
# PLOTS (after your code): day + night
# 1) Separate plots (day, night)  -> component bands pairwise-colored + centers
# 2) Same plot (overlay day vs night) -> envelopes + mixture centers
# ============================================================

# Load original functional objects
day_obj   <- readRDS(day_obj_rds)
night_obj <- readRDS(night_obj_rds)

X_day     <- day_obj$X
tgrid_day <- day_obj$tgrid

X_night     <- night_obj$X
tgrid_night <- night_obj$tgrid

# ----------------------------
# Helper: plot component bands with pairwise colors + centers
# (NEW: axis only)
# ----------------------------
plot_component_bands_pairwise <- function(tgrid, Xmat, band_obj,
                                          main = "",
                                          n_show = 25L,
                                          seed_show = 123,
                                          col_curves = "grey75",
                                          lwd_band = 2,
                                          add_legend = TRUE,
                                          legend_pos = "topright",
                                          draw_mixture_center = TRUE,
                                          draw_component_centers = TRUE,
                                          # NEW plotting-only:
                                          start_hhmm = NULL,
                                          total_min = 600,
                                          major_by = 100,
                                          minor_by = 60,
                                          add_hour_lines = TRUE,
                                          cex.axis = 0.9) {
  comps <- band_obj$component
  keep_idx <- which(!vapply(comps, is.null, logical(1)))
  
  if (length(keep_idx) == 0) stop("No component bands to plot (all components are NULL).")
  
  cols <- if (exists("hcl.colors")) {
    hcl.colors(length(keep_idx), palette = "Dark 3")
  } else {
    rainbow(length(keep_idx))
  }
  
  n <- nrow(Xmat)
  set.seed(seed_show)
  idx <- if (n <= n_show) seq_len(n) else sample.int(n, n_show)
  
  band_lowers <- do.call(cbind, lapply(comps[keep_idx], `[[`, "lower"))
  band_uppers <- do.call(cbind, lapply(comps[keep_idx], `[[`, "upper"))
  y_rng <- range(Xmat[idx, ], band_lowers, band_uppers,
                 band_obj$center$mixture, band_obj$center$meanfd,
                 na.rm = TRUE)
  
  x0 <- min(tgrid)
  xlim_use <- c(x0, x0 + total_min)
  
  plot(NA, xlim = xlim_use, ylim = y_rng,
       xlab = "time of day (HH:MM) / (minute index)",
       ylab = "HR", main = main, xaxt = "n")
  grid()
  if (isTRUE(add_hour_lines)) add_hour_guides(x0 = x0, total_min = total_min, by = minor_by)
  
  for (i in idx) lines(tgrid, Xmat[i, ], col = col_curves)
  
  if (draw_mixture_center) {
    lines(tgrid, band_obj$center$mixture, col = "black", lwd = 3, lty = 3)
  }
  
  for (j in seq_along(keep_idx)) {
    k <- keep_idx[j]
    b <- comps[[k]]
    
    lines(tgrid, b$lower, col = cols[j], lwd = lwd_band)
    lines(tgrid, b$upper, col = cols[j], lwd = lwd_band)
    
    if (draw_component_centers) {
      ck <- band_obj$center$component[[k]]
      if (!is.null(ck)) {
        lines(tgrid, ck, col = cols[j], lwd = 2, lty = 3)
      }
    }
  }
  
  # NEW axis
  if (!is.null(start_hhmm)) {
    add_time_axis(start_hhmm, x0 = x0, total_min = total_min,
                  major_by = major_by, minor_by = minor_by, cex.axis = cex.axis)
  } else {
    axis(1)
  }
  
  if (add_legend) {
    leg <- c("mixture center", paste0("component ", keep_idx), "component centers")
    col_leg <- c("black", cols, "grey30")
    lwd_leg <- c(3, rep(lwd_band, length(cols)), 2)
    lty_leg <- c(3, rep(1, length(cols)), 3)
    
    legend(legend_pos,
           legend = leg,
           col = col_leg,
           lwd = lwd_leg,
           lty = lty_leg,
           bty = "n")
  }
  
  invisible(list(keep_idx = keep_idx, colors = cols))
}

# ----------------------------
# 1) SEPARATE PLOTS (side-by-side)
# ----------------------------
op <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))

plot_component_bands_pairwise(
  tgrid = tgrid_day,
  Xmat  = X_day,
  band_obj = fit_day$band,
  main = paste0("DAY band (alpha=", fit_day$alpha, ", p=", fit_day$p, ", G=", fit_day$params$G, ")"),
  n_show = 25,
  seed_show = 123,
  draw_mixture_center = TRUE,
  draw_component_centers = TRUE,
  start_hhmm = "09:30"
)

plot_component_bands_pairwise(
  tgrid = tgrid_night,
  Xmat  = X_night,
  band_obj = fit_night$band,
  main = paste0("NIGHT band (alpha=", fit_night$alpha, ", p=", fit_night$p, ", G=", fit_night$params$G, ")"),
  n_show = 25,
  seed_show = 123,
  draw_mixture_center = TRUE,
  draw_component_centers = TRUE,
  start_hhmm = "22:00"
)

par(op)

# Save separate plots
png(file.path(out_dir, "section3_bands_pairwise_colored_day_night_with_centers.png"),
    width = 1800, height = 900, res = 150)
op <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))

plot_component_bands_pairwise(
  tgrid = tgrid_day,
  Xmat  = X_day,
  band_obj = fit_day$band,
  main = paste0("DAY band alpha=", fit_day$alpha,
                ", p=", fit_day$p, ", G=", fit_day$params$G),
  n_show = 25,
  seed_show = 123,
  draw_mixture_center = TRUE,
  draw_component_centers = TRUE,
  start_hhmm = "09:30"
)

plot_component_bands_pairwise(
  tgrid = tgrid_night,
  Xmat  = X_night,
  band_obj = fit_night$band,
  main = paste0("NIGHT band alpha=", fit_night$alpha,
                ", p=", fit_night$p, ", G=", fit_night$params$G),
  n_show = 25,
  seed_show = 123,
  draw_mixture_center = TRUE,
  draw_component_centers = TRUE,
  start_hhmm = "22:00"
)

par(op)
dev.off()

# ----------------------------
# 2) OVERLAY DAY vs NIGHT (envelopes + mixture centers)
# (NEW: dual axis only)
# ----------------------------
if (length(tgrid_day) != length(tgrid_night) || any(tgrid_day != tgrid_night)) {
  warning("Day and night tgrid differ. Overlay will use day tgrid only; check alignment.")
}
tgrid <- tgrid_day

day_lower   <- fit_day$band$envelope$lower
day_upper   <- fit_day$band$envelope$upper
night_lower <- fit_night$band$envelope$lower
night_upper <- fit_night$band$envelope$upper

day_center   <- fit_day$band$center$mixture
night_center <- fit_night$band$center$mixture

y_rng <- range(day_lower, day_upper, night_lower, night_upper,
               day_center, night_center, na.rm = TRUE)

x0 <- min(tgrid)
xlim_use <- c(x0, x0 + 600)

plot(NA, xlim = xlim_use, ylim = y_rng,
     xlab = "DAY time | NIGHT time (HH:MM)  /  (minute index)", ylab = "HR",
     main = "Overlay: DAY vs NIGHT conformal bands (envelopes) + centers",
     xaxt = "n")
grid()
add_hour_guides(x0 = x0, total_min = 600, by = 60)

lines(tgrid, day_lower, lwd = 2, col = "dodgerblue3")
lines(tgrid, day_upper, lwd = 2, col = "dodgerblue3")
lines(tgrid, day_center, lwd = 3, col = "dodgerblue3", lty = 3)

lines(tgrid, night_lower, lwd = 2, col = "firebrick3", lty = 2)
lines(tgrid, night_upper, lwd = 2, col = "firebrick3", lty = 2)
lines(tgrid, night_center, lwd = 3, col = "firebrick3", lty = 3)

add_time_axis_dual("09:30", "22:00", x0 = x0, total_min = 600,
                   major_by = 100, minor_by = 60, cex.axis = 0.9)

legend("topright",
       legend = c("Day envelope", "Day center", "Night envelope", "Night center"),
       col = c("dodgerblue3", "dodgerblue3", "firebrick3", "firebrick3"),
       lty = c(1, 3, 2, 3),
       lwd = c(2, 3, 2, 3),
       bty = "n")

png(file.path(out_dir, "section3_band_overlay_day_vs_night_with_centers.png"),
    width = 1400, height = 900, res = 150)

plot(NA, xlim = xlim_use, ylim = y_rng,
     xlab = "DAY time | NIGHT time (HH:MM)  /  (minute index)", ylab = "HR",
     main = "Overlay: DAY vs NIGHT conformal bands (envelopes) + centers",
     xaxt = "n")
grid()
add_hour_guides(x0 = x0, total_min = 600, by = 60)

lines(tgrid, day_lower, lwd = 2, col = "dodgerblue3")
lines(tgrid, day_upper, lwd = 2, col = "dodgerblue3")
lines(tgrid, day_center, lwd = 3, col = "dodgerblue3", lty = 3)

lines(tgrid, night_lower, lwd = 2, col = "firebrick3", lty = 2)
lines(tgrid, night_upper, lwd = 2, col = "firebrick3", lty = 2)
lines(tgrid, night_center, lwd = 3, col = "firebrick3", lty = 3)

add_time_axis_dual("09:30", "22:00", x0 = x0, total_min = 600,
                   major_by = 100, minor_by = 60, cex.axis = 0.9)

# --- legend like "top" but shifted a bit to the RIGHT ---
usr <- par("usr")
xw  <- usr[2] - usr[1]
yw  <- usr[4] - usr[3]

legend(x = usr[1] + 0.50 * xw + 0.25 * xw,   # <-- 0.03 = spostamento a destra (aumenta/diminuisci)
       y = usr[4] - 0.02 * yw,               # vicino al bordo alto
       xjust = 0.5, yjust = 1,
       legend = c("Day envelope (HR)", "Day center (HR)",
                  "Night envelope (HR)", "Night center (HR)"),
       col = c("dodgerblue3", "dodgerblue3", "firebrick3", "firebrick3"),
       lty = c(1, 3, 2, 3),
       lwd = c(2, 3, 2, 3),
       bty = "n")


dev.off()

message("Plots done. Check: ", out_dir)

# ============================================================
# OPTIONAL DIAGNOSTICS (unchanged)
# ============================================================

inspect_clusters <- function(fit_obj, obj_rds_path, which_scores = c(1,2)) {
  obj <- readRDS(obj_rds_path)
  fd_all <- obj$continuous_functions
  
  scores_all <- fpca_scores_for_all(fd_all, fit_obj$pca_train$meanfd, fit_obj$pca_train$harmonics)
  scores_train <- scores_all[fit_obj$idx_train, , drop = FALSE]
  
  gmm <- fit_obj$gmm
  cls <- gmm$classification  # hard assignment for TRAIN curves
  
  cat("\nCluster sizes (TRAIN):\n")
  print(table(cls))
  cat("\nMixing proportions pi_k:\n")
  print(gmm$parameters$pro)
  
  j1 <- which_scores[1]; j2 <- which_scores[2]
  plot(scores_train[, j1], scores_train[, j2],
       col = cls, pch = 16,
       xlab = paste0("score ", j1), ylab = paste0("score ", j2),
       main = paste0("FPCA score space (train) - G=", gmm$G))
  legend("topright", legend = paste0("cluster ", sort(unique(cls))),
         col = sort(unique(cls)), pch = 16, bty = "n")
}

inspect_clusters(fit_day, day_obj_rds)
inspect_clusters(fit_night, night_obj_rds)

summary(apply(fit_day$gmm$z, 1, max))
summary(apply(fit_night$gmm$z, 1, max))

check_nested_and_get_outer <- function(fit_obj) {
  comps <- fit_obj$band$component
  keep <- which(!sapply(comps, is.null))
  if (length(keep) != 2) stop("Serve G=2 con entrambe le componenti keep=TRUE.")
  
  b1 <- comps[[keep[1]]]
  b2 <- comps[[keep[2]]]
  
  b1_contains_b2 <- all(b1$lower <= b2$lower & b1$upper >= b2$upper)
  b2_contains_b1 <- all(b2$lower <= b1$lower & b2$upper >= b1$upper)
  
  list(
    nested = b1_contains_b2 || b2_contains_b1,
    outer  = if (b1_contains_b2) b1 else if (b2_contains_b1) b2 else NULL,
    b1_contains_b2 = b1_contains_b2,
    b2_contains_b1 = b2_contains_b1
  )
}

res_day   <- check_nested_and_get_outer(fit_day)
res_night <- check_nested_and_get_outer(fit_night)

res_day
res_night

check_union_equals_envelope <- function(fit_obj) {
  comps <- fit_obj$band$component
  keep <- which(!sapply(comps, is.null))
  stopifnot(length(keep) == 2)
  
  b1 <- comps[[keep[1]]]
  b2 <- comps[[keep[2]]]
  
  disjoint <- (b1$upper < b2$lower) | (b2$upper < b1$lower)
  
  list(
    any_disjoint = any(disjoint),
    frac_disjoint = mean(disjoint),
    first_disjoint_t = if (any(disjoint)) which(disjoint)[1] else NA_integer_
  )
}

check_union_equals_envelope(fit_day)
check_union_equals_envelope(fit_night)

# ============================================================
# APPEND-ONLY CODE (paste at the END of your script)
# Goal: rerun the SAME conformal pipeline on a POSITIVE scale
#       by applying a LOG transform before FPCA + GMM + conformal.
# ============================================================

# ----------------------------
# 0) Settings
# ----------------------------
eps_log <- 1.0          # safety floor to avoid log(0); 1 bpm is harmless
lambda_smooth <- 1e-2   # smoothing for log-curves (you can tune)
alpha_use <- 0.10
p_use <- 6
n1_frac_use <- 0.70
seed_use <- 123
G_max_use <- 4

# ----------------------------
# 1) Build "log-functional objects" from your saved objects
# ----------------------------
make_log_functional_object <- function(obj_rds_in, obj_rds_out,
                                       eps = 1.0,
                                       lambda = 1e-2,
                                       Lfdobj = 2) {
  obj <- readRDS(obj_rds_in)
  
  X <- obj$X
  tgrid <- obj$tgrid
  
  # ensure strictly positive before log
  X_pos <- pmax(X, eps)
  X_log <- log(X_pos)
  
  # build a new fd object for log(HR) using SAME basis as before
  basisobj <- obj$continuous_functions$basis
  y <- t(X_log)  # [length(tgrid) x ncurves]
  
  fdParobj <- fdPar(basisobj, Lfdobj = Lfdobj, lambda = lambda)
  sm <- smooth.basis(argvals = tgrid, y = y, fdParobj = fdParobj)
  fd_log <- sm$fd
  
  obj_log <- obj
  obj_log$X <- X_log
  obj_log$continuous_functions <- fd_log
  
  # ✅ store transform info at TOP LEVEL (not inside meta)
  obj_log$transform_info <- list(
    type   = "log(pmax(HR, eps))",
    eps    = eps,
    lambda = lambda
  )
  
  saveRDS(obj_log, obj_rds_out)
  message("Saved LOG object: ", obj_rds_out)
  
  invisible(obj_rds_out)
}

# choose where to save the log objects (same out_dir is fine)
day_obj_rds_log   <- file.path(out_dir, "functional_object_day_LOG.rds")
night_obj_rds_log <- file.path(out_dir, "functional_object_night_LOG.rds")

make_log_functional_object(day_obj_rds, day_obj_rds_log, eps = eps_log, lambda = lambda_smooth)
make_log_functional_object(night_obj_rds, night_obj_rds_log, eps = eps_log, lambda = lambda_smooth)

# ----------------------------
# 2) Rerun conformal (same function, but on LOG objects)
# ----------------------------
fit_day_log <- run_section3_one_period(
  obj_rds_path = day_obj_rds_log,
  period_name = "day_log",
  alpha = alpha_use,
  p = p_use,
  n1_frac = n1_frac_use,
  seed = seed_use,
  G_max = G_max_use,
  out_dir = out_dir,
  prefix = "day_LOG"
)

fit_night_log <- run_section3_one_period(
  obj_rds_path = night_obj_rds_log,
  period_name = "night_log",
  alpha = alpha_use,
  p = p_use,
  n1_frac = n1_frac_use,
  seed = seed_use,
  G_max = G_max_use,
  out_dir = out_dir,
  prefix = "night_LOG"
)

message("LOG-scale conformal fits done.")

# ----------------------------
# 3) Plot OVERLAY in ORIGINAL HR SCALE (exp back-transform)
# ----------------------------
tgrid_day_log   <- readRDS(day_obj_rds_log)$tgrid
tgrid_night_log <- readRDS(night_obj_rds_log)$tgrid

if (length(tgrid_day_log) != length(tgrid_night_log) || any(tgrid_day_log != tgrid_night_log)) {
  warning("Day and night tgrid differ in LOG objects. Overlay will use day tgrid only; check alignment.")
}
tgrid <- tgrid_day_log

day_lower_hr   <- exp(fit_day_log$band$envelope$lower)
day_upper_hr   <- exp(fit_day_log$band$envelope$upper)
night_lower_hr <- exp(fit_night_log$band$envelope$lower)
night_upper_hr <- exp(fit_night_log$band$envelope$upper)

day_center_hr   <- exp(fit_day_log$band$center$mixture)
night_center_hr <- exp(fit_night_log$band$center$mixture)

y_rng <- range(day_lower_hr, day_upper_hr, night_lower_hr, night_upper_hr,
               day_center_hr, night_center_hr, na.rm = TRUE)

x0 <- min(tgrid)
xlim_use <- c(x0, x0 + 600)

# --- SAVE TO PNG (correct order: open device -> plot -> close device)
out_png <- file.path(out_dir, "section3_band_overlay_day_vs_night_LOGfit_back_to_HR.png")
png(out_png, width = 1400, height = 900, res = 150)
on.exit(dev.off(), add = TRUE)

plot(NA, xlim = xlim_use, ylim = y_rng,
     xlab = "DAY time | NIGHT time (HH:MM)  /  (minute index)", ylab = "HR",
     main = "Overlay: DAY vs NIGHT conformal bands (LOG-fit, exp back to HR)",
     xaxt = "n")
grid()
add_hour_guides(x0 = x0, total_min = 600, by = 60)

lines(tgrid, day_lower_hr, lwd = 2, col = "dodgerblue3")
lines(tgrid, day_upper_hr, lwd = 2, col = "dodgerblue3")
lines(tgrid, day_center_hr, lwd = 3, col = "dodgerblue3", lty = 3)

lines(tgrid, night_lower_hr, lwd = 2, col = "firebrick3", lty = 2)
lines(tgrid, night_upper_hr, lwd = 2, col = "firebrick3", lty = 2)
lines(tgrid, night_center_hr, lwd = 3, col = "firebrick3", lty = 3)

add_time_axis_dual("09:30", "22:00", x0 = x0, total_min = 600,
                   major_by = 100, minor_by = 60, cex.axis = 0.9)

# --- legend like "top" but shifted a bit to the RIGHT ---
usr <- par("usr")
xw  <- usr[2] - usr[1]
yw  <- usr[4] - usr[3]

legend(x = usr[1] + 0.50 * xw + 0.15 * xw,   # <-- 0.03 = spostamento a destra (aumenta/diminuisci)
       y = usr[4] - 0.02 * yw,               # vicino al bordo alto
       xjust = 0.5, yjust = 1,
       legend = c("Day envelope (HR)", "Day center (HR)",
                  "Night envelope (HR)", "Night center (HR)"),
       col = c("dodgerblue3", "dodgerblue3", "firebrick3", "firebrick3"),
       lty = c(1, 3, 2, 3),
       lwd = c(2, 3, 2, 3),
       bty = "n")

dev.off()
message("Saved: ", out_png)

cat("\nMin lower (DAY, HR scale):  ", min(day_lower_hr, na.rm = TRUE), "\n")
cat("Min lower (NIGHT, HR scale):", min(night_lower_hr, na.rm = TRUE), "\n")

