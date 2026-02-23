# ============================================================
# ONE-SHOT (EN): build HR day/night windows + summary + diagnostics
# ============================================================

suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
  library(lubridate)
  library(pbapply)
})

# ---------------- Helpers ----------------

read_one_binning_file <- function(path) {
  x <- tryCatch(fromJSON(path), error = function(e) NULL)
  if (is.null(x) || length(x) == 0) return(NULL)
  
  dt <- as.data.table(x)
  needed <- c("start_time", "end_time", "heart_rate")
  if (!all(needed %in% names(dt))) return(NULL)
  
  dt[, file := path]
  dt[]
}

load_all_minutes <- function(root_dir) {
  files <- list.files(
    root_dir,
    pattern = "heart_rate\\.binning_data\\.json$",
    recursive = TRUE,
    full.names = TRUE
  )
  files <- files[!grepl("_REPORT_OUT", files, fixed = TRUE)]
  
  message("Found ", length(files), " JSON files. Reading with progress bar...")
  chunks <- pblapply(files, read_one_binning_file)
  dt <- rbindlist(chunks, use.names = TRUE, fill = TRUE)
  
  if (nrow(dt) == 0) stop("No records read. Check the path / pattern.")
  dt[]
}

top_two_na_gaps <- function(v) {
  na <- is.na(v)
  if (!any(na)) return(c(gap1 = 0L, gap2 = 0L))
  r <- rle(na)
  gaps <- r$lengths[r$values]
  if (length(gaps) == 0) return(c(gap1 = 0L, gap2 = 0L))
  gaps <- sort(gaps, decreasing = TRUE)
  gap1 <- gaps[1]
  gap2 <- if (length(gaps) >= 2) gaps[2] else 0L
  c(gap1 = as.integer(gap1), gap2 = as.integer(gap2))
}

make_window_bounds <- function(date_d, start_hm, end_hm, tz) {
  # Fix Date class dropping inside loops (can become numeric)
  if (!inherits(date_d, "Date")) {
    date_d <- as.Date(date_d, origin = "1970-01-01")
  }
  
  st <- ymd_hm(paste(date_d, start_hm), tz = tz)
  en <- ymd_hm(paste(date_d, end_hm), tz = tz)
  
  if (en <= st) en <- en + days(1)
  list(start = st, end = en)
}

make_minute_grid <- function(st, en) {
  seq(from = st, to = en - minutes(1), by = "1 min")
}

# ---------------- Core (build + diagnostics together) ----------------

build_hr_windows_with_diag <- function(input_dir,
                                       tz = "Europe/Rome",
                                       keep_start = "2025-01-01 00:00:00",
                                       keep_end_excl = "2026-02-01 00:00:00",
                                       coverage_thr = 85,
                                       gap1_max = 30,
                                       gap1_trigger = 20,
                                       gap2_max_if_trigger = 20,
                                       verbose = TRUE) {
  
  hr <- load_all_minutes(input_dir)
  
  # ms epoch -> POSIXct UTC -> local tz
  hr[, start_ts_utc := as.POSIXct(start_time / 1000, origin = "1970-01-01", tz = "UTC")]
  hr[, start_ts := with_tz(start_ts_utc, tzone = tz)]
  
  keep_start_ts <- ymd_hms(keep_start, tz = tz)
  keep_end_ts   <- ymd_hms(keep_end_excl, tz = tz)
  
  hr <- hr[start_ts >= keep_start_ts & start_ts < keep_end_ts]
  if (nrow(hr) == 0) stop("No data after date filtering. Check date range / timezone.")
  
  setorder(hr, start_ts)
  hr <- hr[!duplicated(start_ts)]
  setkey(hr, start_ts)
  
  # Candidate dates (same logic as your original)
  date_min <- as.Date(keep_start_ts, tz = tz)
  date_max <- as.Date(keep_end_ts - seconds(1), tz = tz)
  all_dates <- seq.Date(date_min, date_max, by = "day")
  
  if (verbose) {
    message("Building windows for ", length(all_dates), " dates (day + night)...")
  }
  
  # diagnostics counters (same fields as before)
  diag <- data.table(
    period = c("day", "night"),
    total_candidates = 0L,
    in_range = 0L,
    dst_skipped = 0L,
    cov_failed = 0L,
    gap1_failed = 0L,
    gap2_rule_failed = 0L,
    kept = 0L
  )
  bump <- function(p, field) diag[period == p, (field) := get(field) + 1L]
  
  m_names <- sprintf("m%03d", 0:599)
  rows <- vector("list", length(all_dates) * 2L)
  k <- 0L
  
  build_one <- function(d, period, st_hm, en_hm) {
    bump(period, "total_candidates")
    
    b <- make_window_bounds(d, st_hm, en_hm, tz)
    if (!(b$start >= keep_start_ts && b$end <= keep_end_ts)) return(NULL)
    bump(period, "in_range")
    
    grid <- make_minute_grid(b$start, b$end)
    if (length(grid) != 600L) {
      bump(period, "dst_skipped")
      return(NULL)
    }
    
    vals <- hr[data.table(start_ts = grid), heart_rate, on = "start_ts"]
    
    covg <- 100 * sum(!is.na(vals)) / 600
    if (covg < coverage_thr) {
      bump(period, "cov_failed")
      return(NULL)
    }
    
    gaps <- top_two_na_gaps(vals)
    gap1 <- gaps["gap1"]
    gap2 <- gaps["gap2"]
    
    if (gap1 > gap1_max) {
      bump(period, "gap1_failed")
      return(NULL)
    }
    
    if (gap1 > gap1_trigger && gap2 > gap2_max_if_trigger) {
      bump(period, "gap2_rule_failed")
      return(NULL)
    }
    
    bump(period, "kept")
    
    r <- data.table(
      date = as.Date(d, origin = "1970-01-01"),  # safe if numeric
      period = period,
      start_time = b$start,
      end_time = b$end,
      coverage = covg,
      max_gap = gap1,
      second_gap = gap2
    )
    r[, (m_names) := as.list(vals)]
    r
  }
  
  for (d in all_dates) {
    r_day <- build_one(d, "day", "09:30", "19:30")
    if (!is.null(r_day)) { k <- k + 1L; rows[[k]] <- r_day }
    
    r_night <- build_one(d, "night", "22:00", "08:00")
    if (!is.null(r_night)) { k <- k + 1L; rows[[k]] <- r_night }
  }
  
  out <- rbindlist(rows[seq_len(k)], use.names = TRUE, fill = TRUE)
  setorder(out, date, period)
  
  if (verbose) message("Kept windows: ", nrow(out))
  
  # totals + percentages
  total <- diag[, lapply(.SD, sum), .SDcols = names(diag)[-1]]
  pct <- copy(diag)
  pct[, `:=`(
    pct_in_range         = round(100 * in_range / total_candidates, 1),
    pct_dst_skipped      = round(100 * dst_skipped / total_candidates, 1),
    pct_cov_failed       = round(100 * cov_failed / total_candidates, 1),
    pct_gap1_failed      = round(100 * gap1_failed / total_candidates, 1),
    pct_gap2_rule_failed = round(100 * gap2_rule_failed / total_candidates, 1),
    pct_kept             = round(100 * kept / total_candidates, 1)
  )]
  
  list(data = out, diag = diag, total = total, pct = pct)
}

# ---------------- Summary ----------------

print_summary <- function(ds) {
  dt <- as.data.table(ds)
  
  overall <- dt[, .(
    n_total = .N,
    n_day   = sum(period == "day"),
    n_night = sum(period == "night"),
    
    gap1_mean = round(mean(max_gap), 2),
    gap1_min  = min(max_gap),
    gap1_max  = max(max_gap),
    
    gap2_mean = round(mean(second_gap), 2),
    gap2_min  = min(second_gap),
    gap2_max  = max(second_gap),
    
    coverage_mean = round(mean(coverage), 2),
    coverage_min  = min(coverage),
    coverage_max  = max(coverage)
  )]
  
  by_period <- dt[, .(
    n = .N,
    gap1_mean = round(mean(max_gap), 2),
    gap1_min  = min(max_gap),
    gap1_max  = max(max_gap),
    gap2_mean = round(mean(second_gap), 2),
    gap2_min  = min(second_gap),
    gap2_max  = max(second_gap),
    coverage_mean = round(mean(coverage), 2),
    coverage_min  = min(coverage),
    coverage_max  = max(coverage)
  ), by = period][order(period)]
  
  cat("\n=== OVERALL SUMMARY ===\n")
  print(overall)
  cat("\n=== SUMMARY BY PERIOD ===\n")
  print(by_period)
}

# ============================================================
# RUN
# ============================================================

# to be modified (INPUT) -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
data_dir <- "C:\\Users\\2003l\\OneDrive\\Documenti\\Minor\\SDS_R\\Project functional\\heart_rate\\com.samsung.shealth.tracker.heart_rate"

params <- list(
  tz = "Europe/Rome",
  keep_start = "2025-01-01 00:00:00",
  keep_end_excl = "2026-02-01 00:00:00",
  coverage_thr = 85,
  gap1_max = 30,
  gap1_trigger = 20,
  gap2_max_if_trigger = 20
)

res <- do.call(build_hr_windows_with_diag, c(list(input_dir = data_dir), params))

full_clean  <- res$data
day_clean   <- full_clean[period == "day"]
night_clean <- full_clean[period == "night"]

# Save outputs (current working directory)
saveRDS(full_clean,  "full_clean.rds")
saveRDS(day_clean,   "day_clean.rds")
saveRDS(night_clean, "night_clean.rds")

fwrite(full_clean,  "full_clean.csv")
fwrite(day_clean,   "day_clean.csv")
fwrite(night_clean, "night_clean.csv")

# Also save CSVs inside data_dir
fwrite(full_clean,  file.path(data_dir, "full_clean.csv"))
fwrite(day_clean,   file.path(data_dir, "day_clean.csv"))
fwrite(night_clean, file.path(data_dir, "night_clean.csv"))

# Summaries
print_summary(full_clean)
cat("\n--- DAY ONLY ---\n");   print_summary(day_clean)
cat("\n--- NIGHT ONLY ---\n"); print_summary(night_clean)

# Diagnostics
cat("\n==================== DIAGNOSTICS (TOTAL) ====================\n")
print(as.data.table(res$total))

cat("\n==================== DIAGNOSTICS (BY PERIOD) ================\n")
print(res$diag)

cat("\n==================== PERCENTAGES (ON CANDIDATES) =============\n")
print(res$pct[, .(period, total_candidates,
                  pct_in_range, pct_dst_skipped, pct_cov_failed,
                  pct_gap1_failed, pct_gap2_rule_failed, pct_kept)])

