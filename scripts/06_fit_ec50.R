# 06_fit_ec50.R
# Stage 6 (Analysis): plasma EC50 from occupancy vs plasma drug concentration.
#   Emax model: occupancy = occ_max * C / (C + EC50), occ_max fixed (default 1)
fit_ec50 <- function(lassen_fits, exposure, occ_max = 1, start_ec50 = 10) {

  # One point per scan: occupancy (y) and plasma concentration (x)
  ok  <- lassen_fits[!lassen_fits$qc_exclude & !is.na(lassen_fits$occupancy), ]
  pts <- merge(ok[, c("subject", "session", "occupancy", "occupancy_se")],
               exposure[, c("subject", "session", "dose_mg", "conc_ng_mL")],
               by = c("subject", "session"))
  
  # Emax model: occ = occ_max * C / (C + EC50), fitted on log(EC50)
  fit <- nls(occupancy ~ occ_max * conc_ng_mL / (conc_ng_mL + exp(log_ec50)),
             data = pts, start = list(log_ec50 = log(start_ec50)))
  
  cf <- coef(summary(fit))
  log_ec50 <- cf["log_ec50", "Estimate"]
  log_se   <- cf["log_ec50", "Std. Error"]
  ci <- exp(log_ec50 + c(-1, 1) * qt(0.975, df.residual(fit)) * log_se)
  
  pts$occupancy_fit <- fitted(fit)
  pts$residual      <- residuals(fit)

  result <- data.frame(EC50_ng_mL = exp(log_ec50),
                       ci_low     = ci[1],
                       ci_high    = ci[2],
                       occ_max    = occ_max,
                       n_scans    = nrow(pts),
                       resid_sd   = sigma(fit))

  list(result = result, points = pts)
}