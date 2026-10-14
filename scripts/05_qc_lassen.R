# 05_qc_lassen.R
# Stage 5 (QC Lassen): check each Lassen fit (is the occupancy usable?)
# and each point (is a region far from its Lassen line?).

qc_lassen <- function(lassen, occ_range = c(-0.2, 1.2), min_r2 = 0.5,
                      max_resid_z = 2.5) {
  stage <- "qc_lassen"
  fits  <- lassen$fits
  pairs <- lassen$pairs
  
  # Lowest baseline VT in each scan: VND must lie between 0 and this value
  min_vt <- aggregate(VT_base ~ subject + session, data = pairs, FUN = min)
  names(min_vt)[names(min_vt) == "VT_base"] <- "min_VT_base"
  fits <- merge(fits, min_vt, by = c("subject", "session"), all.x = TRUE)
  
  # Residuals scaled by their SD within each scan
  resid_z <- ave(pairs$residual, pairs$subject, pairs$session,
                 FUN = function(r) r / sd(r))
  
  flags <- rbind(
    # Per Lassen fit
    make_flags(fits, is.na(fits$occupancy), stage, "too_few_regions", "exclude",
               fits$n_regions, "too few regions for a Lassen fit"),
    make_flags(fits, fits$occupancy < occ_range[1] | fits$occupancy > occ_range[2],
               stage, "implausible_occupancy", "exclude",
               fits$occupancy, paste0("occupancy outside ", occ_range[1], " to ",
                                      occ_range[2])),
    make_flags(fits, fits$VND <= 0 | fits$VND >= fits$min_VT_base,
               stage, "implausible_vnd", "warning",
               fits$VND, "VND not between 0 and the lowest baseline VT"),
    make_flags(fits, fits$r_squared < min_r2, stage, "low_r2", "warning",
               fits$r_squared, paste0("Lassen R2 < ", min_r2)),
    # Per point (region) in a Lassen plot
    make_flags(pairs, abs(resid_z) > max_resid_z, stage, "lassen_outlier", "warning",
               resid_z, paste0("region far from the Lassen line (|residual z| > ",
                               max_resid_z, ")"))
  )
  
  fits$min_VT_base <- NULL
  list(data  = apply_flags(fits, flags[is.na(flags$region), ], c("subject", "session")),
       flags = flags)
}