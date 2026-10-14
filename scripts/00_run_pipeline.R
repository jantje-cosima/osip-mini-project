# 00_run_pipeline.R
# Master script: runs the whole pipeline from data/raw/ to the final report.
#
# Run from the project root:
#   Rscript scripts/00_run_pipeline.R
# or, interactively in Positron:
#   source("scripts/00_run_pipeline.R"); res <- run_pipeline()

# --- Load function definitions ---------------------------------------------------
source(file.path("scripts", "qc_utils.R"))
source(file.path("scripts", "01_load_data.R"))
source(file.path("scripts", "02_qc_input.R"))
source(file.path("scripts", "03_drug_exposure.R"))
source(file.path("scripts", "04_fit_lassen.R"))
source(file.path("scripts", "05_qc_lassen.R"))
source(file.path("scripts", "06_fit_ec50.R"))
source(file.path("scripts", "07_qc_ec50.R"))
source(file.path("scripts", "08_make_outputs.R"))
source(file.path("scripts", "optional_refit_vt.R"))
# ...

run_pipeline <- function(refit = FALSE, tac_qc = FALSE, qc_report = TRUE,
                         report = TRUE) {

  # 1. Load
  vt   <- load_vt()
  drug <- load_drug()

  # Optional: replace VT with a fresh 2TCM fit
  if (refit) {
    stop("refit = TRUE is not implemented yet")
  }

  # 2. QC input data
  qc_v <- qc_vt(vt)
  qc_d <- qc_drug(drug)
  vt   <- qc_v$data      # now with qc_exclude / qc_reason columns
  drug <- qc_d$data
  qc_flags <- rbind(qc_v$flags, qc_d$flags)
  message("QC input: ", sum(qc_flags$severity == "exclude"), " exclusion flag(s), ",
          sum(qc_flags$severity == "warning"), " warning(s)")
  
  # Optional: QC of the TACs and arterial input behind the VT estimates
  plasma <- NULL
  tacs   <- NULL
  if (tac_qc) {
    tacs   <- load_tacs()
    plasma <- load_plasma()
    qc_t <- qc_tacs(tacs)
    qc_p <- qc_plasma(plasma)
    plasma   <- qc_p$data           # now with qc_flagged per sample
    tacs     <- add_fitted_tacs(tacs, plasma, vt)   # 2TCM curve from the VT fits
    qc_flags <- rbind(qc_flags, qc_t$flags, qc_p$flags)
    message("QC TACs/input: ", nrow(qc_t$flags) + nrow(qc_p$flags), " warning(s)")
  }

  # 3. Drug exposure
  exposure <- drug_exposure(drug)

  # 4. Lassen fits
  lassen <- fit_lassen(vt)

  # 5. QC Lassen
  qc_l <- qc_lassen(lassen)
  lassen$fits <- qc_l$data        # now with qc_exclude / qc_reason
  qc_flags <- rbind(qc_flags, qc_l$flags)
  message("QC Lassen: ", sum(qc_l$data$qc_exclude), " of ", nrow(qc_l$data),
          " fits excluded")

  # 6. EC50 fit
  ec50 <- fit_ec50(lassen$fits, exposure)
  message(sprintf("EC50: %.1f ng/mL (95%% CI %.1f-%.1f, %d scans)",
                  ec50$result$EC50_ng_mL, ec50$result$ci_low,
                  ec50$result$ci_high, ec50$result$n_scans))

  # 7. QC EC50
  qc_e <- qc_ec50(ec50)
  ec50$points <- qc_e$data
  qc_flags <- rbind(qc_flags, qc_e$flags)
  message("QC EC50: ", nrow(qc_e$flags), " warning(s)")

  # Collect results
  res <- list(vt = vt, drug = drug, exposure = exposure, lassen = lassen,
              ec50 = ec50, tacs = tacs, plasma = plasma, qc_flags = qc_flags)

  # 8. Figures and tables
  make_outputs(res)
  message("Outputs written to results/")

  # 9. Save results; the reports read this file instead of re-running the analysis
  saveRDS(res, file.path("results", "pipeline_results.rds"))

  # Render the QC report
  if (qc_report) {
    quarto::quarto_render(file.path("reports", "qc_report.qmd"), quiet = TRUE)
    message("QC report: reports/qc_report.html")
  }

  # Render the final report
  if (report) {
    quarto::quarto_render(file.path("reports", "report.qmd"), quiet = TRUE)
    message("Report: reports/report.html")
  }
  
  message("Pipeline finished.")
  invisible(res)
}

# Run automatically only when called with Rscript, not when sourced
if (sys.nframe() == 0) {
  run_pipeline()
}