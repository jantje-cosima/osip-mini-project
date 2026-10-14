# 02_qc_input.R
# Stage 2 (QC input data): flag problems in the loaded data.
# QC never removes data; it returns flags + a qc_exclude column.

qc_vt <- function(vt, max_se_pct = 25, max_robust_z = 3.5, k4_lower = 0.005,
                  sessions = c("baseline", "block1", "block2")) {
  stage <- "qc_vt"

  # Robust z-score of VT within each region, across subjects, at baseline only.
  # (Block sessions differ by dose, so between-subject spread there is real;
  #  block outliers are checked later against the Lassen line.)
  base <- vt$session == "baseline"
  z <- rep(NA_real_, nrow(vt))
  z[base] <- ave(vt$VT[base], vt$region[base],
                 FUN = function(x) (x - median(x, na.rm = TRUE)) / mad(x, na.rm = TRUE))

  flags <- rbind(
    make_flags(vt, !vt$converged, stage, "not_converged", "exclude",
               vt$VT, "2TCM fit did not converge"),
    make_flags(vt, vt$VT_SE_pct > max_se_pct, stage, "high_se", "exclude",
               vt$VT_SE_pct, paste0("VT SE > ", max_se_pct, "%")),
    make_flags(vt, is.na(vt$VT) | vt$VT <= 0, stage, "invalid_vt", "exclude",
               vt$VT, "VT missing or <= 0"),
    make_flags(vt, vt$at_bound, stage, "at_bound", "warning",
               vt$VT, "rate constant at fit bound; SE may be underestimated"),
    make_flags(vt, vt$k4 <= k4_lower * 1.001, stage, "k4_at_lower_bound", "exclude",
               vt$k4, "k4 at its lower fit bound: 2TCM fit degenerate, VT unreliable"),           
    make_flags(vt, abs(z) > max_robust_z, stage, "regional_outlier", "warning",
               z, paste0("baseline VT deviates from other subjects (|robust z| > ",
                         max_robust_z, ")"))
  )

  # Completeness: every subject should have all sessions x all regions.
  # Missing rows cannot be flagged on vt itself, so compare with the full grid.
  expected <- expand.grid(subject = unique(vt$subject),
                          session = sessions,
                          region  = unique(vt$region),
                          stringsAsFactors = FALSE)
  key <- function(d) paste(d$subject, d$session, d$region)
  missing <- !(key(expected) %in% key(vt))
  flags <- rbind(flags,
                 make_flags(expected, missing, stage, "missing_row", "warning",
                            rep(NA_real_, nrow(expected)), "subject x session x region not in data"))

  list(data  = apply_flags(vt, flags, c("subject", "session", "region")),
       flags = flags)
}

qc_drug <- function(drug, ratio_range = c(0.5, 2)) {
  stage <- "qc_drug"

  # One row per scan: start and end sample side by side
  start <- drug[drug$sample == "scan_start",
                c("subject", "session", "dose_mg", "time_post_dose_h", "drug_ng_mL")]
  end   <- drug[drug$sample == "scan_end",
                c("subject", "session", "time_post_dose_h", "drug_ng_mL")]
  scans <- merge(start, end, by = c("subject", "session"), all = TRUE,
                 suffixes = c("_start", "_end"))

  ratio  <- scans$drug_ng_mL_end / scans$drug_ng_mL_start
  n_dose <- ave(scans$dose_mg, scans$subject,
                FUN = function(x) length(unique(x[!is.na(x)])))

  flags <- rbind(
    # Per sample
    make_flags(drug, is.na(drug$drug_ng_mL) | drug$drug_ng_mL <= 0, stage,
               "invalid_conc", "exclude",
               drug$drug_ng_mL, "drug concentration missing or <= 0"),
    # Per scan
    make_flags(scans, is.na(scans$drug_ng_mL_start) | is.na(scans$drug_ng_mL_end),
               stage, "missing_sample", "warning",
               rep(NA_real_, nrow(scans)), "scan lacks a start or end drug sample"),
    make_flags(scans, scans$time_post_dose_h_end <= scans$time_post_dose_h_start,
               stage, "time_order", "warning",
               scans$time_post_dose_h_end, "end sample not later than start sample"),
    make_flags(scans, ratio < ratio_range[1] | ratio > ratio_range[2],
               stage, "conc_ratio", "warning",
               ratio, paste0("end/start concentration ratio outside ",
                             ratio_range[1], "-", ratio_range[2])),
    make_flags(scans, n_dose > 1, stage, "dose_mismatch", "warning",
               scans$dose_mg, "subject has different doses in different sessions")
  )

  list(data  = apply_flags(drug, flags, c("subject", "session")),
       flags = flags)
}