# ============================================================
# Daily coverage (% minutes present) — Samsung Health heart_rate
# From scratch: JSON -> unique minutes per day -> coverage plot
# ============================================================

suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
  library(ggplot2)
})

# to be modified (INPUT) -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
data_dir <- "C:/Users/2003l/OneDrive/Documenti/Minor/SDS_R/Project functional/heart_rate/com.samsung.shealth.tracker.heart_rate"
out_dir  <- file.path(data_dir, "_REPORT_OUT")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

tz <- "Europe/Rome"

# --------- Helpers ---------
list_json_files <- function(root) {
  files <- list.files(root, pattern = "\\.json$", recursive = TRUE, full.names = TRUE)
  files[!grepl("_REPORT_OUT", files, fixed = TRUE)]
}

extract_start_times <- function(obj) {
  if (is.data.frame(obj) && "start_time" %in% names(obj)) return(obj$start_time)
  
  if (is.list(obj)) {
    if (!is.null(obj$start_time)) return(obj$start_time)
    for (nm in names(obj)) {
      x <- obj[[nm]]
      if (is.data.frame(x) && "start_time" %in% names(x)) return(x$start_time)
      if (is.list(x) && !is.null(x$start_time)) return(x$start_time)
    }
  }
  numeric(0)
}

read_file_start_times <- function(f) {
  x <- tryCatch(jsonlite::fromJSON(f, flatten = TRUE), error = function(e) NULL)
  if (is.null(x)) return(numeric(0))
  st <- extract_start_times(x)
  st <- suppressWarnings(as.numeric(st))
  st[is.finite(st)]
}

# --------- 1) Read all timestamps ---------
files <- list_json_files(data_dir)
if (length(files) == 0) stop("No .json files found under data_dir.")

start_time_raw <- unlist(lapply(files, read_file_start_times), use.names = FALSE)
start_time_raw <- start_time_raw[is.finite(start_time_raw)]
if (length(start_time_raw) < 2) stop("Too few timestamps: start_time not found (or not readable).")

# ms vs sec auto-detect
med <- median(start_time_raw, na.rm = TRUE)
start_time_sec <- if (med > 1e11) start_time_raw / 1000 else start_time_raw

dt <- as.POSIXct(start_time_sec, origin = "1970-01-01", tz = tz)

# --------- 2) Floor to minute and count unique minutes per day ---------
minute_epoch <- floor(as.numeric(dt) / 60) * 60
minute_dt <- as.POSIXct(minute_epoch, origin = "1970-01-01", tz = tz)

# Build a proper table and group correctly
min_dt <- data.table(minute_dt = minute_dt)
min_dt[, date := as.Date(minute_dt, tz = tz)]

# Unique minutes *within day*
cov_dt <- min_dt[
  , .(minutes_present = uniqueN(minute_dt),
      expected = 1440L),
  by = date
][order(date)]

cov_dt[, coverage_pct := 100 * minutes_present / expected]

# quick sanity check
cat("Coverage % summary:\n")
print(summary(cov_dt$coverage_pct))
cat("Max minutes_present:", max(cov_dt$minutes_present), "\n")  # should be <= 1440 usually

# Save table
fwrite(cov_dt, file.path(out_dir, "daily_coverage.csv"))

# --------- 3) Plot (English title) ---------
p <- ggplot(cov_dt, aes(x = date, y = coverage_pct)) +
  geom_line(linewidth = 0.4, color = "black") +
  labs(
    title = "Daily coverage (% minutes present)",
    x = "date",
    y = "coverage %"
  ) +
  coord_cartesian(ylim = c(0, 100)) +   # IMPORTANT: doesn’t drop points like ylim()
  theme_grey(base_size = 18)

ggsave(
  filename = file.path(out_dir, "daily_coverage_en.png"),
  plot = p,
  width = 11, height = 6.5, dpi = 200
)

print(p)
message("Saved plot: ", file.path(out_dir, "daily_coverage_en.png"))

ggsave("C:/Users/2003l/Downloads/gap_plot.png", plot = p, width = 11, height = 6.5, dpi = 300)




# ============================================================
# Daily coverage plot — ONLY 2025-01-13 to 2026-01-31
# (uses your existing cov_dt)
# ============================================================

start_date <- as.Date("2025-01-13")
end_date   <- as.Date("2026-01-31")

# Filter the interval
cov_dt_sub <- cov_dt[date >= start_date & date <= end_date]

# Plot (same style)
p_sub <- ggplot(cov_dt_sub, aes(x = date, y = coverage_pct)) +
  geom_line(linewidth = 0.4, color = "black") +
  labs(
    title = "Daily coverage (% minutes present)",
    x = "date",
    y = "coverage %"
  ) +
  scale_x_date(limits = c(start_date, end_date),
               date_breaks = "2 months",
               date_labels = "%Y-%m") +
  coord_cartesian(ylim = c(0, 100)) +
  theme_grey(base_size = 18)

print(p_sub)

# Save plot + table
ggsave(
  filename = file.path(out_dir, "daily_coverage_2025_to_2026-01.png"),
  plot = p_sub,
  width = 11, height = 6.5, dpi = 300
)

fwrite(cov_dt_sub, file.path(out_dir, "daily_coverage_2025_to_2026-01.csv"))
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  
  
  
  
  suppressPackageStartupMessages({
    library(data.table)
    library(lubridate)
    library(ggplot2)
  })
tz <- "Europe/Rome"


# to be modified (INPUT) -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
data_dir <- "C:/Users/2003l/OneDrive/Documenti/Minor/SDS_R/Project functional/heart_rate/com.samsung.shealth.tracker.heart_rate"

out_dir <- file.path(data_dir, "_REPORT_OUT")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# to be modified (INPUT) -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source("C:\\Users\\2003l\\OneDrive\\Documenti\\Minor\\SDS_R\\Project functional\\script functional data1.R")
hr <- load_all_minutes(data_dir)
med <- median(hr$start_time, na.rm = TRUE)
div <- if (med > 1e11) 1000 else 1
hr[, start_ts_utc := as.POSIXct(start_time / div, origin = "1970-01-01", tz = "UTC")]
hr[, start_ts := with_tz(start_ts_utc, tzone = tz)]
hr[, heart_rate := as.numeric(heart_rate)]
hr <- hr[!is.na(start_ts) & !is.na(heart_rate), .(start_ts, heart_rate)]
setorder(hr, start_ts)
hr <- hr[!duplicated(start_ts)]






suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(ggplot2)
})

set.seed(123)
available_days <- sort(unique(as.Date(hr$start_ts, tz = tz)))

random_day <- sample(available_days, 1)
# random_day <- as.Date("2025-07-10")


day_raw <- hr[as.Date(start_ts, tz = tz) == random_day & !is.na(heart_rate)]
if (nrow(day_raw) == 0) stop("No data for the chosen day.")


day_raw[, minute_ts := floor_date(start_ts, unit = "minute")]

day_min <- day_raw[, .(heart_rate_1min = mean(heart_rate, na.rm = TRUE)), by = minute_ts]
setorder(day_min, minute_ts)


day_start <- as.POSIXct(paste0(format(random_day, "%Y-%m-%d"), " 00:00:00"), tz = tz)
grid_1440 <- data.table(minute_ts = seq(day_start, by = "1 min", length.out = 1440L))
day_1440 <- day_min[grid_1440, on = "minute_ts"]
day_1440[, time_of_day := format(minute_ts, "%H:%M")]

# ------------------------------------------------------------
# Missing minutes as red dots at y = 0
# ------------------------------------------------------------
day_1440[, missing := is.na(heart_rate_1min)]
day_1440[, y_plot := fifelse(missing, 0, heart_rate_1min)]

p_points <- ggplot() +
  # missing minutes (red at y=0)
  geom_point(
    data = day_1440[missing == TRUE],
    aes(x = minute_ts, y = y_plot),
    color = "red3",
    size = 0.55,
    alpha = 0.9,
    na.rm = TRUE
  ) +
  # observed minutes (black)
  geom_point(
    data = day_1440[missing == FALSE],
    aes(x = minute_ts, y = y_plot),
    color = "black",
    size = 0.55,
    alpha = 0.85,
    na.rm = TRUE
  ) +
  labs(
    title = sprintf("Heart rate observations (1-min grid) — %s", format(random_day, "%Y-%m-%d")),
    subtitle = "Black = observed minute, Red = missing minute",
    x = "time of day",
    y = "heart rate (bpm)"
  ) +
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "2 hours") +
  coord_cartesian(ylim = c(0, max(day_1440$heart_rate_1min, na.rm = TRUE) * 1.05)) +
  theme_grey(base_size = 18)

print(p_points)

out_png <- file.path(out_dir, sprintf("heart_rate_points_1440_with_missing_%s.png", format(random_day, "%Y-%m-%d")))
ggsave(out_png, plot = p_points, width = 11, height = 6.5, dpi = 300)

cat("Saved plot to:\n", out_png, "\n")
cat("Missing minutes:", sum(day_1440$missing), "/ 1440\n")


out_dir <- file.path(data_dir, "_REPORT_OUT")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_png <- file.path(out_dir, sprintf("heart_rate_points_1440_%s.png", format(random_day, "%Y-%m-%d")))
ggsave(out_png, plot = p_points, width = 11, height = 6.5, dpi = 300)

cat("Saved 1440-point plot to:\n", out_png, "\n")
cat("Minutes with data:", sum(!is.na(day_1440$heart_rate_1min)), " / 1440\n")

day_1440[, hour := format(minute_ts, "%H")]
max_obs_1h <- day_1440[!is.na(heart_rate_1min), .N, by = hour][, max(N)]
cat("Max osservations in 1 hour:", max_obs_1h, "\n")

