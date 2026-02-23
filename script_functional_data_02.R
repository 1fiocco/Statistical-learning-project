suppressPackageStartupMessages({
  library(data.table)
  library(funData)
  library(zoo)
  library(fda)
})

make_discrete_and_continuous <- function(csv_path,
                                         impute_na = TRUE,
                                         tz = "Europe/Rome",
                                         nbasis = 35,
                                         norder = 4) {
  
  # 0) Read data
  dt <- fread(csv_path)
  
  # 1) Identify and order minute columns (m000..m599)
  mcols <- grep("^m\\d{3}$", names(dt), value = TRUE)
  mcols <- mcols[order(as.integer(sub("^m", "", mcols)))]
  stopifnot(length(mcols) == 600)
  
  # 2) Matrix: n_curves x 600
  X <- as.matrix(dt[, ..mcols])
  mode(X) <- "numeric"
  
  # 3) Grid (minutes since window start)
  tgrid <- 0:599
  
  # 4) Optional: impute missing values within each curve
  if (impute_na) {
    impute_curve <- function(x) {
      if (all(is.na(x))) return(x)  # if fully NA, stays NA
      x <- zoo::na.approx(x, x = tgrid, na.rm = FALSE, rule = 2)
      x <- zoo::na.locf(x, na.rm = FALSE)
      x <- zoo::na.locf(x, fromLast = TRUE, na.rm = FALSE)
      x
    }
    X <- t(apply(X, 1, impute_curve))
  }
  
  # Drop curves still all-NA (if any)
  keep <- rowSums(!is.na(X)) > 0
  X <- X[keep, , drop = FALSE]
  dt <- dt[keep]
  
  # 5) Metadata (labels kept separately)
  meta <- dt[, .(date, period, start_time, end_time, coverage, max_gap, second_gap)]
  meta[, date := as.IDate(date)]
  meta[, start_time := as.POSIXct(start_time, tz = tz)]
  meta[, end_time   := as.POSIXct(end_time,   tz = tz)]
  
  # 6A) DISCRETE functional object (funData)
  discrete_functions <- funData(argvals = list(tgrid), X = X)
  
  # 6B) CONTINUOUS functional object (fda::fd via B-splines)
  # smooth.basis expects y as (length(tgrid) x ncurves)
  Y <- t(X)  # 600 x ncurves
  
  basis <- create.bspline.basis(rangeval = c(min(tgrid), max(tgrid)),
                                nbasis = nbasis, norder = norder)
  
  sm <- smooth.basis(argvals = tgrid, y = Y, fdParobj = basis)
  continuous_functions <- sm$fd
  
  list(
    discrete_functions = discrete_functions,
    continuous_functions = continuous_functions,
    X = X,
    tgrid = tgrid,
    meta = meta,
    smoothing = list(nbasis = nbasis, norder = norder, basis = basis, smooth = sm)
  )
}

# ------------------------------------------------------------
# NEW: run the pipeline on 3 datasets + save separate outputs
# ------------------------------------------------------------
build_three_datasets <- function(full_csv,
                                 day_csv,
                                 night_csv,
                                 impute_na = TRUE,
                                 tz = "Europe/Rome",
                                 nbasis = 35,
                                 norder = 4,
                                 out_dir = ".",
                                 save_full_objects = TRUE) {
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  paths <- list(full = full_csv, day = day_csv, night = night_csv)
  
  results <- lapply(names(paths), function(name) {
    p <- paths[[name]]
    
    obj <- make_discrete_and_continuous(
      csv_path = p,
      impute_na = impute_na,
      tz = tz,
      nbasis = nbasis,
      norder = norder
    )
    
    # Save RDS outputs
    saveRDS(obj$discrete_functions,
            file = file.path(out_dir, paste0("discrete_functions_", name, ".rds")))
    saveRDS(obj$continuous_functions,
            file = file.path(out_dir, paste0("continuous_functions_", name, ".rds")))
    
    if (isTRUE(save_full_objects)) {
      saveRDS(obj,
              file = file.path(out_dir, paste0("functional_object_", name, ".rds")))
    }
    
    obj
  })
  names(results) <- names(paths)
  
  results
}

# to be modified (INPUT) -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
full_csv  <- "C:\\Users\\2003l\\OneDrive\\Documenti\\Minor\\SDS_R\\Project functional\\heart_rate\\com.samsung.shealth.tracker.heart_rate\\full_clean.csv"
day_csv   <- "C:\\Users\\2003l\\OneDrive\\Documenti\\Minor\\SDS_R\\Project functional\\heart_rate\\com.samsung.shealth.tracker.heart_rate\\day_clean.csv"
night_csv <- "C:\\Users\\2003l\\OneDrive\\Documenti\\Minor\\SDS_R\\Project functional\\heart_rate\\com.samsung.shealth.tracker.heart_rate\\night_clean.csv"

# OUTPUT 
out_dir <- "C:\\Users\\2003l\\OneDrive\\Documenti\\Minor\\SDS_R\\Project functional\\heart_rate\\outputs_rds"

res <- build_three_datasets(
  full_csv  = full_csv,
  day_csv   = day_csv,
  night_csv = night_csv,
  impute_na = TRUE,
  nbasis = 35,
  norder = 4,
  out_dir = out_dir
)

# Required objects (like before), but now you have 3 versions:
discrete_full   <- res$full$discrete_functions
continuous_full <- res$full$continuous_functions

discrete_day    <- res$day$discrete_functions
continuous_day  <- res$day$continuous_functions

discrete_night  <- res$night$discrete_functions
continuous_night<- res$night$continuous_functions










# ============================================================
# REPORT-READY PLOTS (continuous lines everywhere, points only where observed)
# 4 single-window plots: discrete+continuous for DAY and NIGHT
# + 2 spaghetti plots (all windows) with clean English labels
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(fda)
  library(zoo)
})

# --- output dir (use your existing out_dir variable)
plot_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 0) Load RAW minutes (NA gaps preserved) from the clean CSVs
# ============================================================
dt_day_raw   <- fread(day_csv)
dt_night_raw <- fread(night_csv)

dt_day_raw[,   date := as.IDate(date)]
dt_night_raw[, date := as.IDate(date)]

mcols <- sprintf("m%03d", 0:599)

# --- choose ONE date present in BOTH day and night (safe strings)
common_dates_chr <- intersect(as.character(dt_day_raw$date), as.character(dt_night_raw$date))
if (length(common_dates_chr) == 0) stop("No common dates between day and night datasets.")
date_str    <- "2025-12-13"
chosen_date <- as.IDate(date_str)
# --- extract raw vectors for that date (NA preserved)
get_raw_vec <- function(dt, d) {
  r <- dt[date == d]
  if (nrow(r) == 0) stop("No row for chosen date in raw CSV.")
  as.numeric(r[1, ..mcols])
}

y_day_raw   <- get_raw_vec(dt_day_raw, chosen_date)
y_night_raw <- get_raw_vec(dt_night_raw, chosen_date)

# ============================================================
# 1) Get smoothed curves from your fda objects (continuous_day/night)
#    (line must stay continuous -> DO NOT break on missing)
# ============================================================
meta_day   <- as.data.table(res$day$meta);   meta_day[,   date := as.IDate(date)]
meta_night <- as.data.table(res$night$meta); meta_night[, date := as.IDate(date)]

i_day   <- which(meta_day$date   == chosen_date)[1]
i_night <- which(meta_night$date == chosen_date)[1]
if (is.na(i_day) || is.na(i_night)) stop("Could not match chosen_date in res$...$meta.")

tgrid <- 0:599
yC_day   <- as.numeric(eval.fd(tgrid, continuous_day)[, i_day])
yC_night <- as.numeric(eval.fd(tgrid, continuous_night)[, i_night])

# ============================================================
# 2) Axis ticks: HH:MM on top and (minute index) below
#    Use FIXED window start times as requested
# ============================================================
tz <- "Europe/Rome"
st_day   <- ymd_hms(paste0(date_str, " 09:30:00"), tz = tz)
st_night <- ymd_hms(paste0(date_str, " 22:00:00"), tz = tz)

make_ticks_labels <- function(start_time_posix) {
  ticks <- c(0, 100, 200, 300, 400, 500, 600)
  times <- start_time_posix + minutes(ticks)
  labs  <- sprintf("%s\n(%d)", format(times, "%H:%M"), ticks)
  list(ticks = ticks, labs = labs)
}

# ============================================================
# 3) Imputation ONLY for drawing the DISCRETE LINE (to keep it continuous)
#    Points will still come ONLY from raw observations.
# ============================================================
impute_for_line <- function(y) {
  x <- 0:599
  if (all(is.na(y))) return(rep(NA_real_, length(y)))

  # linear interpolation over existing points
  y_imp <- zoo::na.approx(y, x = x, na.rm = FALSE, rule = 2)

  # carry forward/backward remaining NA (start/end)
  y_imp <- zoo::na.locf(y_imp, na.rm = FALSE)
  y_imp <- zoo::na.locf(y_imp, fromLast = TRUE, na.rm = FALSE)

  y_imp
}

y_day_line   <- impute_for_line(y_day_raw)
y_night_line <- impute_for_line(y_night_raw)

# ============================================================
# 4) Plot functions (LINE always black and continuous; POINTS colored)
# ============================================================
col_line      <- "black"
col_day_pts   <- "royalblue3"
col_night_pts <- "navy"

plot_discrete_window <- function(y_line, y_raw, start_time_posix, main, pch, col_pts, out_png) {
  ticks <- make_ticks_labels(start_time_posix)
  x <- 0:599
  ok <- which(!is.na(y_raw))

  png(out_png, width = 1400, height = 800, res = 160)
  par(mar = c(5, 5, 3.5, 1))

  ylim_rng <- range(c(y_line, y_raw), na.rm = TRUE)
  plot(NA, xlim = c(0, 600), ylim = ylim_rng,
       xlab = "time of day (HH:MM) / (minute index)",
       ylab = "heart rate (bpm)",
       main = main, xaxt = "n")
  grid()

  # Continuous discrete line (always drawn)
  lines(x, y_line, col = col_line, lwd = 1)

  # Points ONLY where observed
  points(x[ok], y_raw[ok], pch = pch, col = col_pts, cex = 0.55)

  axis(1, at = ticks$ticks, labels = ticks$labs, cex.axis = 0.9)
  invisible(dev.off())
}

plot_continuous_window <- function(y_cont, y_raw, start_time_posix, main, pch, col_pts, out_png) {
  ticks <- make_ticks_labels(start_time_posix)
  x <- 0:599
  ok <- which(!is.na(y_raw))

  png(out_png, width = 1400, height = 800, res = 160)
  par(mar = c(5, 5, 3.5, 1))

  ylim_rng <- range(c(y_cont, y_raw), na.rm = TRUE)
  plot(NA, xlim = c(0, 600), ylim = ylim_rng,
       xlab = "time of day (HH:MM) / (minute index)",
       ylab = "heart rate (bpm)",
       main = main, xaxt = "n")
  grid()

  # Smooth curve always continuous
  lines(x, y_cont, col = col_line, lwd = 1.6)

  # Points ONLY where observed
  points(x[ok], y_raw[ok], pch = pch, col = col_pts, cex = 0.50)

  axis(1, at = ticks$ticks, labels = ticks$labs, cex.axis = 0.9)
  invisible(dev.off())
}

# ============================================================
# 5) NIGHT: print approx sleep drop & wake rise (computed on RAW)
# ============================================================
ys <- zoo::rollmedian(y_night_raw, k = 11, fill = NA, align = "center")
d  <- diff(ys)

d2 <- d; d2[!is.finite(d2)] <- Inf
drop_idx <- which.min(d2) + 1

d3 <- d; d3[!is.finite(d3)] <- -Inf
if (drop_idx < length(d3)) d3[1:(drop_idx-1)] <- -Inf
rise_idx <- which.max(d3) + 1

drop_time <- st_night + minutes(drop_idx - 1)
rise_time <- st_night + minutes(rise_idx - 1)

cat("\nChosen date:", date_str, "\n")
cat("NIGHT window start:", format(st_night, "%Y-%m-%d %H:%M"), "\n")
cat("Approx. sharp drop (sleep) at:", format(drop_time, "%H:%M"),
    " (minute", drop_idx - 1, ")\n")
cat("Approx. sharp rise (wake) at:", format(rise_time, "%H:%M"),
    " (minute", rise_idx - 1, ")\n\n")

# ============================================================
# 6) Save the 4 plots
# Markers: NIGHT ★ (pch=8), DAY △ (pch=2)
# Line: black; Points: colored
# ============================================================
plot_discrete_window(
  y_line = y_day_line, y_raw = y_day_raw, start_time_posix = st_day,
  main = sprintf("Discrete heart rate (DAY) — %s", date_str),
  pch = 2, col_pts = col_day_pts,
  out_png = file.path(plot_dir, sprintf("plot_discrete_day_%s.png", date_str))
)

plot_discrete_window(
  y_line = y_night_line, y_raw = y_night_raw, start_time_posix = st_night,
  main = sprintf("Discrete heart rate (NIGHT) — %s", date_str),
  pch = 8, col_pts = col_night_pts,
  out_png = file.path(plot_dir, sprintf("plot_discrete_night_%s.png", date_str))
)

plot_continuous_window(
  y_cont = yC_day, y_raw = y_day_raw, start_time_posix = st_day,
  main = sprintf("Smoothed heart rate (DAY) — %s", date_str),
  pch = 2, col_pts = col_day_pts,
  out_png = file.path(plot_dir, sprintf("plot_continuous_day_%s.png", date_str))
)

plot_continuous_window(
  y_cont = yC_night, y_raw = y_night_raw, start_time_posix = st_night,
  main = sprintf("Smoothed heart rate (NIGHT) — %s", date_str),
  pch = 8, col_pts = col_night_pts,
  out_png = file.path(plot_dir, sprintf("plot_continuous_night_%s.png", date_str))
)

cat("Saved 4 single-window plots in:\n", plot_dir, "\n\n")

# ============================================================
# 7) Spaghetti plots (ALL windows)
# Requirement: lines continuous -> use an imputed matrix for drawing lines
# We do: impute each row like above (approx + locf) and draw alpha-black lines
# ============================================================
impute_matrix_for_lines <- function(Xraw) {
  Ximp <- Xraw
  for (i in seq_len(nrow(Ximp))) {
    Ximp[i, ] <- impute_for_line(as.numeric(Ximp[i, ]))
  }
  Ximp
}

X_day_raw_mat   <- as.matrix(dt_day_raw[, ..mcols]);   mode(X_day_raw_mat) <- "numeric"
X_night_raw_mat <- as.matrix(dt_night_raw[, ..mcols]); mode(X_night_raw_mat) <- "numeric"

X_day_line_all   <- impute_matrix_for_lines(X_day_raw_mat)
X_night_line_all <- impute_matrix_for_lines(X_night_raw_mat)

plot_spaghetti <- function(X_line, start_time_posix, main, out_png) {
  ticks <- make_ticks_labels(start_time_posix)
  x <- 0:599
  col_alpha <- rgb(0, 0, 0, 0.12)

  png(out_png, width = 1600, height = 900, res = 160)
  par(mar = c(5, 5, 3.5, 1))

  ylim_rng <- range(X_line, na.rm = TRUE)
  plot(NA, xlim = c(0, 600), ylim = ylim_rng,
       xlab = "time of day (HH:MM) / (minute index)",
       ylab = "heart rate (bpm)",
       main = main, xaxt = "n")
  grid()

  matlines(x, t(X_line), lty = 1, lwd = 0.7, col = col_alpha)

  axis(1, at = ticks$ticks, labels = ticks$labs, cex.axis = 0.9)
  invisible(dev.off())
}

plot_spaghetti(
  X_day_line_all, st_day,
  main = "Heart rate curves — DAY windows (all days)",
  out_png = file.path(plot_dir, "spaghetti_day_all.png")
)

plot_spaghetti(
  X_night_line_all, st_night,
  main = "Heart rate curves — NIGHT windows (all days)",
  out_png = file.path(plot_dir, "spaghetti_night_all.png")
)

cat("Saved spaghetti plots in:\n", plot_dir, "\n")






# ============================================================
# PLOT + TIME AXIS + AUTO NIGHT WINDOW (NO WARNINGS)
# - robust: selects curve manually via eval.fd()
# ============================================================

library(fda)

# --- ensure fd: converts fd-like lists (with $coefs and $basis) into a proper 'fd' object
ensure_fd <- function(x) {
  if (inherits(x, "fd")) return(x)
  
  if (is.list(x) && !is.null(x$coefs) && !is.null(x$basis)) {
    fdnames <- if (!is.null(x$fdnames)) x$fdnames else NULL
    return(fda::fd(coef = x$coefs, basisobj = x$basis, fdnames = fdnames))
  }
  
  stop("Object is not an 'fd' and cannot be auto-converted (missing $coefs/$basis).")
}

# --- custom time axis (HH:MM labels)
add_time_axis <- function(start_hhmm, fdobj = NULL, by = 60) {
  parts <- strsplit(start_hhmm, ":", fixed = TRUE)[[1]]
  start_min <- as.integer(parts[1]) * 60 + as.integer(parts[2])
  
  if (!is.null(fdobj) && inherits(fdobj, "fd")) {
    rg <- fdobj$basis$rangeval
    x_start <- rg[1]
    n <- round(rg[2] - rg[1] + 1)
  } else {
    x_start <- 0
    n <- 600
  }
  
  at <- seq(x_start, x_start + n - 1, by = by)
  labs_min <- (start_min + (at - x_start)) %% (24 * 60)
  labs <- sprintf("%02d:%02d", labs_min %/% 60, labs_min %% 60)
  
  axis(1, at = at, labels = labs)
}

# --- helper: plot ONE curve of an fd object (idx) without using plot.fd(index=...)
plot_one_fd_curve <- function(fdobj, idx, start_hhmm, file_png, title_prefix) {
  if (!inherits(fdobj, "fd")) stop("fdobj must be an 'fd' object.")
  ncurves <- ncol(fdobj$coefs)
  if (idx < 1 || idx > ncurves) stop("idx out of range.")
  
  rng <- fdobj$basis$rangeval
  tgrid <- seq(rng[1], rng[2], by = 1)        # minute grid
  Y <- eval.fd(tgrid, fdobj)                  # length(tgrid) x ncurves
  y <- Y[, idx]
  
  png(file_png, width = 1200, height = 800, res = 160)
  plot(tgrid, y, type = "l", xaxt = "n",
       xlab = "time of day (HH:MM) / (minute index)",
       ylab = "heart rate (bpm)",
       main = sprintf("Smoothed HR — %s (curve %d)", title_prefix, idx))
  add_time_axis(start_hhmm, fdobj = fdobj)
  dev.off()
  
  invisible(TRUE)
}

# --- coerce to fd (robust)
continuous_full  <- ensure_fd(continuous_full)
continuous_day   <- ensure_fd(continuous_day)
continuous_night <- ensure_fd(continuous_night)

# --- choose last curve index
n_full  <- ncol(continuous_full$coefs)
n_day   <- ncol(continuous_day$coefs)
n_night <- ncol(continuous_night$coefs)

# --- save plots
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

plot_one_fd_curve(
  fdobj = continuous_full,
  idx = n_full,
  start_hhmm = "00:00",
  file_png = file.path(out_dir, sprintf("plot_full_curve_%d.png", n_full)),
  title_prefix = "FULL"
)

plot_one_fd_curve(
  fdobj = continuous_day,
  idx = n_day,
  start_hhmm = "09:30",
  file_png = file.path(out_dir, sprintf("plot_day_curve_%d.png", n_day)),
  title_prefix = "DAY"
)

plot_one_fd_curve(
  fdobj = continuous_night,
  idx = n_night,
  start_hhmm = "22:00",
  file_png = file.path(out_dir, sprintf("plot_night_curve_%d.png", n_night)),
  title_prefix = "NIGHT"
)

# ============================================================
# --- find the XXX-XXX minutes automatically (night)
# ============================================================

thr    <- 60    # bpm threshold
p_cut  <- 0.70  # "vast majority" below thr (change if you want)
exc_cut<- 0.10  # "few exceptions" above thr (change if you want)

rng   <- continuous_night$basis$rangeval
tgrid <- seq(rng[1], rng[2], by = 1)
Y     <- eval.fd(tgrid, continuous_night)   # length(tgrid) x N

p_low  <- rowMeans(Y <  thr, na.rm = TRUE)
p_high <- rowMeans(Y >= thr, na.rm = TRUE)

longest_run <- function(flag){
  r <- rle(flag); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  if (!any(r$values)) return(c(NA, NA))
  ok <- which(r$values)
  k <- ok[which.max(r$lengths[ok])]
  c(starts[k], ends[k])
}

sleep_idx <- longest_run(p_low >= p_cut)
sleep_min <- tgrid[sleep_idx[1]]; sleep_max <- tgrid[sleep_idx[2]]

inside <- seq(sleep_idx[1], sleep_idx[2])
exc_idx_local <- longest_run(p_high[inside] >= exc_cut)
exc_min <- if (!anyNA(exc_idx_local)) tgrid[inside[exc_idx_local[1]]] else NA
exc_max <- if (!anyNA(exc_idx_local)) tgrid[inside[exc_idx_local[2]]] else NA

cat(sprintf("SLEEP: minute %d to %d\n", sleep_min, sleep_max))
cat(sprintf("EXCEPTIONS (HR>%d): minute %s to %s\n", thr, exc_min, exc_max))






# ============================================================
# FINAL GRAPH (ONE PLOT): 6 curves together
# - smoothed (continuous) curves
# - observed points only where present (NA -> no marker)
# - DAY = reds/oranges (left legend column)
# - NIGHT = blues (right legend column)
# - x-axis: DAY HH:MM | NIGHT HH:MM  +  (minute index)
# - legend centered with symmetric margins
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
})

# ------------------------------------------------------------
# 0) Load the clean tables (if you already have day_clean/night_clean, skip fread)
# ------------------------------------------------------------
if (!exists("day_clean"))   day_clean   <- fread(day_csv)
if (!exists("night_clean")) night_clean <- fread(night_csv)

day_clean[,   date := as.IDate(date)]
night_clean[, date := as.IDate(date)]

m_cols <- sprintf("m%03d", 0:599)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
extract_curve_600 <- function(one_row_dt) {
  as.numeric(one_row_dt[, ..m_cols])
}

make_curve_label <- function(one_row_dt, prefix) {
  d <- format(as.IDate(one_row_dt$date), "%Y-%m-%d")
  paste0(prefix, " ", d)
}

pick_top3 <- function(dt_period) {
  if (nrow(dt_period) == 0) return(integer(0))
  tmp <- copy(dt_period)
  tmp[, row_id := .I]
  # Prefer high coverage, then small gaps (if columns exist)
  has_gap_cols <- all(c("coverage","max_gap","second_gap") %in% names(tmp))
  if (has_gap_cols) {
    setorder(tmp, -coverage, max_gap, second_gap)
  } else if ("coverage" %in% names(tmp)) {
    setorder(tmp, -coverage)
  }
  head(tmp$row_id, min(3L, nrow(tmp)))
}

build_curve_infos <- function(dt_period, idx, colors, pchs, prefix) {
  infos <- vector("list", length(idx))
  for (j in seq_along(idx)) {
    i <- idx[j]
    y <- extract_curve_600(dt_period[i])
    infos[[j]] <- list(
      y     = y,
      label = make_curve_label(dt_period[i], prefix),
      col   = colors[j],
      pch   = pchs[j]
    )
  }
  infos
}

# ------------------------------------------------------------
# 1) Choose curves (deterministic: top 3 by coverage/gaps)
# ------------------------------------------------------------
idx_day   <- pick_top3(day_clean)
idx_night <- pick_top3(night_clean)

# If for some reason you have <3, still works
curve_pch <- c(8, 1, 2)  # star, circle, triangle

# DAY = reds/oranges, NIGHT = blues  (as you requested)
day_cols   <- c("red3", "darkorange3", "orange")
night_cols <- c("navy", "royalblue3", "deepskyblue2")

day_infos   <- build_curve_infos(day_clean,   idx_day,   day_cols,   curve_pch, "DAY")
night_infos <- build_curve_infos(night_clean, idx_night, night_cols, curve_pch, "NIGHT")

# Plot order (draw NIGHT first, then DAY on top)
infos_plot <- c(night_infos, day_infos)


# Legend order must be: DAY (left column) then NIGHT (right column)
leg_labels <- c(vapply(day_infos, `[[`, "", "label"),
                vapply(night_infos, `[[`, "", "label"))
leg_cols   <- c(vapply(day_infos, `[[`, "", "col"),
                vapply(night_infos, `[[`, "", "col"))
leg_pch    <- c(vapply(day_infos, `[[`, 1, "pch"),
                vapply(night_infos, `[[`, 1, "pch"))

if (length(infos_plot) == 0) stop("No curves available to plot (day_clean/night_clean empty).")

# ------------------------------------------------------------
# 2) X-axis labels: DAY HH:MM | NIGHT HH:MM  + (minute index)
# ------------------------------------------------------------
tz <- "Europe/Rome"
st_day_ref   <- ymd_hms("2000-01-01 09:30:00", tz = tz)
st_night_ref <- ymd_hms("2000-01-01 22:00:00", tz = tz)

make_ticks_labels_dual <- function(st_day, st_night) {
  ticks <- c(0, 100, 200, 300, 400, 500, 600)
  t_day   <- st_day   + minutes(ticks)
  t_night <- st_night + minutes(ticks)
  labs <- sprintf("%s | %s\n(%d)",
                  format(t_day, "%H:%M"),
                  format(t_night, "%H:%M"),
                  ticks)
  list(ticks = ticks, labs = labs)
}
ticks <- make_ticks_labels_dual(st_day_ref, st_night_ref)

# ------------------------------------------------------------
# 3) Y-range from observed points (ignore NA)
# ------------------------------------------------------------
all_obs <- unlist(lapply(infos_plot, function(ci) ci$y), use.names = FALSE)
y_rng <- range(all_obs, na.rm = TRUE)
if (!all(is.finite(y_rng))) y_rng <- c(0, 1)
pad <- 0.05 * diff(y_rng); if (!is.finite(pad) || pad == 0) pad <- 1
y_rng <- y_rng + c(-pad, pad)

# ------------------------------------------------------------
# 4) Plot + save
# ------------------------------------------------------------
out_png <- file.path(plot_dir, "plot_6curves_day_reds_night_blues.png")

x_min  <- 0:599
x_fine <- seq(0, 599, by = 0.25)

png(out_png, width = 1800, height = 950, res = 170)
par(mar = c(6, 5, 4, 2))

plot(NA,
     xlim = c(0, 600), ylim = y_rng,
     xlab = "DAY time | NIGHT time (HH:MM)  /  (minute index)",
     ylab = "heart rate (bpm)",
     main = "Day & Night windows — smoothed curves with observed points",
     xaxt = "n")
grid()

# draw curves + observed points
for (ci in infos_plot) {
  y <- ci$y
  ok <- which(!is.na(y))
  if (length(ok) < 2) next
  
  x_obs <- x_min[ok]
  y_obs <- y[ok]
  
  # smoothed continuous curve (always drawn)
  if (length(x_obs) >= 4) {
    fit <- smooth.spline(x_obs, y_obs)
    y_smooth <- predict(fit, x_fine)$y
  } else {
    y_smooth <- approx(x_obs, y_obs, xout = x_fine, rule = 2)$y
  }
  
  lines(x_fine, y_smooth, col = ci$col, lwd = 2)
  points(x_obs, y_obs, pch = ci$pch, col = ci$col, cex = 0.70)
}

axis(1, at = ticks$ticks, labels = ticks$labs, cex.axis = 0.95)

# ------------------------------------------------------------
# 5) Legend centered (symmetric margins), 2 columns:
#    left = DAY, right = NIGHT
# ------------------------------------------------------------
usr <- par("usr")
y_top <- usr[4] - 0.02 * (usr[4] - usr[3])

leg_args <- list(
  legend    = leg_labels,
  col       = leg_cols,
  lwd       = 2,
  pch       = leg_pch,
  text.col  = leg_cols,
  bty       = "n",
  ncol      = 2,
  xjust     = 0,
  yjust     = 1,
  x.intersp = 1.0,
  y.intersp = 1.15,
  cex       = 0.95
)

tmp <- do.call(legend, c(list(x = usr[1], y = y_top, plot = FALSE), leg_args))
plot_w <- usr[2] - usr[1]
leg_w  <- tmp$rect$w
x_final <- usr[1] + (plot_w - leg_w) / 2
x_final <- max(usr[1], min(x_final, usr[2] - leg_w))

do.call(legend, c(list(x = x_final, y = y_top), leg_args))

dev.off()

cat("Saved combined plot here:\n", out_png, "\n")


