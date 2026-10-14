# optional_refit_vt.R
# Optional branch: TACs and arterial input behind the VT estimates.
# Load and QC functions (run with run_pipeline(tac_qc = TRUE)).
# A refit_vt() function (2TCM) can be added here later.

# --- Load ---------------------------------------------------------------------------

load_tacs <- function(path = file.path("data", "raw", "tacs.csv")) {
  if (!file.exists(path)) stop("File not found: ", path)
  df <- read.csv(path, stringsAsFactors = FALSE)

  keys     <- c("subject", "session", "region", "frame")
  num_cols <- c("frame_start_min", "frame_end_min", "activity_kBq_mL")

  check_columns(df, c(keys, num_cols), path)
  check_types(df, c("frame", num_cols), is.numeric, "numeric", path)
  check_unique_keys(df, keys, path)
  df
}

load_plasma <- function(path = file.path("data", "raw", "plasma_input.csv")) {
  if (!file.exists(path)) stop("File not found: ", path)
  df <- read.csv(path, stringsAsFactors = FALSE)

  keys     <- c("subject", "session", "time_min")
  num_cols <- c("time_min", "plasma_parent_kBq_mL")

  check_columns(df, c(keys, num_cols), path)
  check_types(df, num_cols, is.numeric, "numeric", path)
  check_unique_keys(df, keys, path)
  df
}

# --- QC of TACs --------------------------------------------------------------------

qc_tacs <- function(tacs, early_min = 1) {
  stage  <- "qc_tacs"
  tacs   <- tacs[order(tacs$subject, tacs$session, tacs$region, tacs$frame), ]
  tac_id <- paste(tacs$subject, tacs$session, tacs$region)

  # Number of frames per TAC (subject x session x region)
  n_frames <- ave(tacs$frame, tac_id, FUN = length)

  # Negative activity after the first minute (noise can make the very
  # first frames slightly negative; later frames should not be)
  late_neg <- tacs$activity_kBq_mL < 0 & tacs$frame_start_min >= early_min

  # Frame gaps or overlaps: each frame should start where the previous ended
  prev_end <- ave(tacs$frame_end_min, tac_id, FUN = function(x) c(NA, head(x, -1)))
  gap      <- !is.na(prev_end) & abs(tacs$frame_start_min - prev_end) > 1e-3

  first <- !duplicated(tac_id)   # one row per TAC, for TAC-level flags
  flags <- rbind(
    make_flags(tacs[first, ], n_frames[first] != max(n_frames), stage,
               "missing_frames", "warning", n_frames[first],
               "TAC has fewer frames than the others"),
    make_flags(tacs, late_neg, stage, "negative_activity", "warning",
               tacs$activity_kBq_mL, paste0("negative activity after ", early_min, " min")),
    make_flags(tacs, gap, stage, "frame_gap", "warning",
               tacs$frame_start_min - prev_end, "frame does not start where the previous ended")
  )
  list(data = tacs, flags = flags)
}

# --- QC of the arterial input ---------------------------------------------------------

qc_plasma <- function(plasma, max_peak_min = 3, max_tail_dev_pct = 50) {
  stage  <- "qc_plasma"
  plasma <- plasma[order(plasma$subject, plasma$session, plasma$time_min), ]
  curve_id <- paste(plasma$subject, plasma$session)
  conc   <- plasma$plasma_parent_kBq_mL

  # Peak time per curve
  is_peak  <- ave(conc, curve_id, FUN = function(x) seq_along(x) == which.max(x)) == 1
  peaks    <- plasma[is_peak, ]

  # Tail samples (after the peak) that deviate from their neighbours:
  # compare each sample with log-linear interpolation between the samples
  # before and after it
  after_peak <- ave(conc, curve_id, FUN = function(x) seq_along(x) > which.max(x)) == 1
  tail_dev <- ave(seq_along(conc), curve_id, FUN = function(i) {
    t <- plasma$time_min[i]; y <- conc[i]; n <- length(y); d <- rep(NA_real_, n)
    for (k in 2:(n - 1)) {
      if (y[k - 1] > 0 && y[k + 1] > 0) {
        pred <- exp(approx(t[c(k - 1, k + 1)], log(y[c(k - 1, k + 1)]), xout = t[k])$y)
        d[k] <- 100 * (y[k] / pred - 1)
      }
    }
    d
  })

  flags <- rbind(
    make_flags(plasma, is.na(conc), stage, "missing_sample", "warning",
               conc, "plasma sample missing"),
    make_flags(plasma, after_peak & conc < 0, stage, "negative_after_peak", "warning",
               conc, "negative concentration after the peak"),
    make_flags(peaks, peaks$time_min > max_peak_min, stage, "late_peak", "warning",
               peaks$time_min, paste0("peak later than ", max_peak_min, " min")),
    make_flags(plasma, after_peak & abs(tail_dev) > max_tail_dev_pct, stage,
               "tail_outlier", "warning", tail_dev,
               paste0("sample deviates > ", max_tail_dev_pct,
                      "% from its neighbours (contamination?)"))
  )

  # Mark the flagged samples in the data, for plotting
  plasma$qc_flagged <- (after_peak & conc < 0) |
    (after_peak & !is.na(tail_dev) & abs(tail_dev) > max_tail_dev_pct)
  list(data = plasma, flags = flags)
}

# --- Fitted 2TCM curves --------------------------------------------------------------

# Frame-mean 2TCM curve for one set of rate constants and one measured input
tcm2_frames <- function(K1, k2, k3, k4, plasma_t, plasma_c, frame_start, frame_end,
                        vB = 0.05, dt = 1 / 60) {
  t_grid <- seq(0, max(frame_end), by = dt)
  cp <- pmax(approx(c(0, plasma_t), c(0, plasma_c), xout = t_grid, rule = 2)$y, 0)

  # Convolution of the input with exp(-a t), trapezoidal recursion
  conv <- function(a) {
    e <- exp(-a * dt)
    x <- dt / 2 * (cp + e * c(0, head(cp, -1)))
    as.numeric(stats::filter(x, e, method = "recursive"))
  }
  s  <- k2 + k3 + k4
  r  <- max(sqrt(max(s^2 - 4 * k2 * k4, 0)), 1e-8)
  a1 <- (s - r) / 2
  a2 <- (s + r) / 2
  ct  <- K1 / (a2 - a1) * ((k3 + k4 - a1) * conv(a1) + (a2 - k3 - k4) * conv(a2))
  pet <- (1 - vB) * ct + vB * cp

  # Mean over the grid points inside each frame
  first <- ceiling(frame_start / dt - 1e-9) + 1
  last  <- ceiling(frame_end / dt - 1e-9)
  cs <- c(0, cumsum(pet))
  (cs[last + 1] - cs[first]) / (last - first + 1)
}
# Add the fitted curve (from the K1-k4 in vt) to every TAC frame
add_fitted_tacs <- function(tacs, plasma, vt) {
  tacs <- tacs[order(tacs$subject, tacs$session, tacs$region, tacs$frame), ]
  tacs$fitted <- NA_real_
  for (i in seq_len(nrow(vt))) {
    v  <- vt[i, ]
    ti <- tacs$subject == v$subject & tacs$session == v$session & tacs$region == v$region
    pi <- plasma$subject == v$subject & plasma$session == v$session
    if (!any(ti) || !any(pi)) next
    tacs$fitted[ti] <- tcm2_frames(v$K1, v$k2, v$k3, v$k4,
                                   plasma$time_min[pi], plasma$plasma_parent_kBq_mL[pi],
                                   tacs$frame_start_min[ti], tacs$frame_end_min[ti])
  }
  tacs
}
