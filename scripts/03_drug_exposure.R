# 03_drug_exposure.R
# Stage 3 (Analysis): one plasma drug concentration per PET scan,
# used as the exposure measure in the EC50 fit.

drug_exposure <- function(drug) {
  # Usable samples only (QC exclusions from qc_drug())
  ok <- drug[!drug$qc_exclude, ]

  # Mean concentration per scan (start and end sample)
  mean_conc <- aggregate(drug_ng_mL ~ subject + session + dose_mg,
                         data = ok, FUN = mean)
  names(mean_conc)[names(mean_conc) == "drug_ng_mL"] <- "conc_ng_mL"

  # Number of samples behind each mean (2 normally; 1 if one was excluded)
  n_conc <- aggregate(drug_ng_mL ~ subject + session, data = ok, FUN = length)
  names(n_conc)[names(n_conc) == "drug_ng_mL"] <- "n_samples"

  # One row per scan
  exposure <- merge(mean_conc, n_conc, by = c("subject", "session"))
  exposure[order(exposure$subject, exposure$session), ]
}
