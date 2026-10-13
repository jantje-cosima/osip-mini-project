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
# source(file.path("scripts", "03_drug_exposure.R"))
# ...

run_pipeline <- function(refit = FALSE) {

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
  # 3. Drug exposure
  # 4. Lassen fits
  # 5. QC Lassen
  # 6. EC50 fit
  # 7. QC EC50
  # 8. Figures and tables
  # 9. Render report

  message("Pipeline finished.")
  invisible(list(vt = vt, drug = drug, qc_flags = qc_flags))
}

# Run automatically only when called with Rscript, not when sourced
if (sys.nframe() == 0) {
  run_pipeline()
}