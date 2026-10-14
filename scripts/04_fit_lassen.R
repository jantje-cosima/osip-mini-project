# 04_fit_lassen.R
# Stage 4 (Analysis): Lassen plot per subject x block session.
#   y = VT_baseline - VT_block,  x = VT_baseline
#   slope = occupancy,  x-intercept = VND
fit_lassen <- function(vt, min_regions = 4) {

  # Usable VT only
  ok <- vt[!vt$qc_exclude, c("subject", "session", "region", "VT")]
 
  # Pair each block VT with the baseline VT of the same subject and region
  base  <- ok[ok$session == "baseline", c("subject", "region", "VT")]
  names(base)[names(base) == "VT"] <- "VT_base"
  block <- ok[ok$session != "baseline", ]
  names(block)[names(block) == "VT"] <- "VT_block"
  pairs <- merge(base, block, by = c("subject", "region"))
  pairs$delta <- pairs$VT_base - pairs$VT_block 
  
  # One Lassen regression per subject x block session
  scans <- unique(pairs[, c("subject", "session")])
  fits <- lapply(seq_len(nrow(scans)), function(i) {
    p <- pairs[pairs$subject == scans$subject[i] & pairs$session == scans$session[i], ]
    out <- data.frame(subject = scans$subject[i], session = scans$session[i],
                      n_regions = nrow(p), occupancy = NA_real_, occupancy_se = NA_real_,
                      intercept = NA_real_, VND = NA_real_, r_squared = NA_real_)
    if (nrow(p) < min_regions) return(out)

    fit <- lm(delta ~ VT_base, data = p)
    cf  <- coef(summary(fit))
    out$occupancy    <- cf["VT_base", "Estimate"]
    out$occupancy_se <- cf["VT_base", "Std. Error"]
    out$intercept    <- cf["(Intercept)", "Estimate"]
    out$VND          <- -out$intercept / out$occupancy
    out$r_squared    <- summary(fit)$r.squared
    out
  })
  fits <- do.call(rbind, fits)
  fits <- fits[order(fits$subject, fits$session), ]
  
  # Fitted line and residual for every point (used by QC and plots)
  pairs <- merge(pairs, fits[, c("subject", "session", "occupancy", "intercept")],
                 by = c("subject", "session"))
  pairs$delta_fit <- pairs$intercept + pairs$occupancy * pairs$VT_base
  pairs$residual  <- pairs$delta - pairs$delta_fit
  pairs$occupancy <- NULL
  pairs$intercept <- NULL

  list(fits = fits, pairs = pairs)
}