# 07_qc_ec50.R
# Stage 7 (QC EC50): is the EC50 estimate robust?
# Per scan: distance from the Emax curve and influence (leave-one-out).
# Whole fit: CI width, number of scans, EC50 within the measured range.

qc_ec50 <- function(ec50, max_resid_z = 2.5, max_loo_change_pct = 15,
                    max_ci_ratio = 2, min_scans = 8) {
  stage   <- "qc_ec50"
  pts     <- ec50$points
  result  <- ec50$result
  occ_max <- result$occ_max
  
  # Residuals scaled by the residual SD of the fit
  resid_z <- pts$residual / result$resid_sd
  
  # Leave-one-out: refit without each scan, % change in EC50
  loo_ec50 <- vapply(seq_len(nrow(pts)), function(i) {
    p <- pts[-i, ]
    f <- nls(occupancy ~ occ_max * conc_ng_mL / (conc_ng_mL + exp(log_ec50)),
             data = p, start = list(log_ec50 = log(result$EC50_ng_mL)))
    exp(coef(f)[["log_ec50"]])
  }, numeric(1))
  loo_change <- 100 * (loo_ec50 / result$EC50_ng_mL - 1)
  
  # Fit-level checks use a one-row table without subject/session
  fit_row <- data.frame(subject = NA_character_, session = NA_character_)
  conc    <- pts$conc_ng_mL

  flags <- rbind(
    # Per scan
    make_flags(pts, abs(resid_z) > max_resid_z, stage, "ec50_outlier", "warning",
               resid_z, paste0("scan far from the Emax curve (|residual z| > ",
                               max_resid_z, ")")),
    make_flags(pts, abs(loo_change) > max_loo_change_pct, stage, "influential_scan",
               "warning", loo_change,
               paste0("leaving this scan out changes EC50 by > ",
                      max_loo_change_pct, "%")),
    # Whole fit
    make_flags(fit_row, result$ci_high / result$ci_low > max_ci_ratio, stage,
               "wide_ci", "warning", result$ci_high / result$ci_low,
               paste0("95% CI upper/lower ratio > ", max_ci_ratio)),
    make_flags(fit_row, result$n_scans < min_scans, stage, "few_scans", "warning",
               result$n_scans, paste0("fewer than ", min_scans, " scans in the fit")),
    make_flags(fit_row, result$EC50_ng_mL < min(conc) | result$EC50_ng_mL > max(conc),
               stage, "ec50_outside_data", "warning", result$EC50_ng_mL,
               "EC50 outside the measured concentration range (extrapolated)")
  )
  
  pts$loo_change_pct <- loo_change
  list(data = apply_flags(pts, flags, c("subject", "session")), flags = flags)
}