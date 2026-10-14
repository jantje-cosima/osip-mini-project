# 08_make_outputs.R
# Stage 8 (Outputs): final figures and tables in results/, for the report,
# slides and README.

# --- Plot helpers ----------------------------------------------------------------

# One Lassen plot: points, fitted line, occupancy and VND in the title
plot_lassen <- function(fit, pairs, xlim, ylim) {
  p <- pairs[pairs$subject == fit$subject & pairs$session == fit$session, ]
  plot(p$VT_base, p$delta, pch = 19, xlim = xlim, ylim = ylim,
       xlab = expression(V[T] ~ "baseline (mL/cm"^3 * ")"),
       ylab = expression(V[T] ~ "baseline" - V[T] ~ "block"),
       main = sprintf("%s %s\noccupancy %.2f, VND %.1f", fit$subject, fit$session,
                      fit$occupancy, fit$VND), cex.main = 0.95)
  abline(h = 0, col = "grey80")
  abline(a = fit$intercept, b = fit$occupancy, col = "steelblue", lwd = 2)
}

# EC50 curve: occupancy vs plasma concentration with the Emax fit
plot_ec50 <- function(ec50) {
  er  <- ec50$result
  ept <- ec50$points
  conc_grid <- 10^seq(log10(min(ept$conc_ng_mL) / 2), log10(max(ept$conc_ng_mL) * 2),
                      length.out = 200)
  emax <- function(conc, ec) er$occ_max * conc / (conc + ec)

  plot(ept$conc_ng_mL, ept$occupancy, log = "x", pch = 19, ylim = c(-0.1, 1.1),
       col = ifelse(ept$session == "block1", "black", "darkorange"),
       xlab = "Plasma drug concentration (ng/mL)", ylab = "Receptor occupancy")
  polygon(c(conc_grid, rev(conc_grid)),
          c(emax(conc_grid, er$ci_low), rev(emax(conc_grid, er$ci_high))),
          col = adjustcolor("steelblue", 0.2), border = NA)
  lines(conc_grid, emax(conc_grid, er$EC50_ng_mL), col = "steelblue", lwd = 2)
  abline(h = c(0, er$occ_max), col = "grey85")
  legend("bottomright", bty = "n", cex = 0.85,
         legend = c("block 1 (~3 h post dose)", "block 2 (~24 h post dose)",
                    sprintf("Emax fit, EC50 = %.1f ng/mL (95%% CI %.1f-%.1f)",
                            er$EC50_ng_mL, er$ci_low, er$ci_high)),
         pch = c(19, 19, NA), lty = c(NA, NA, 1), lwd = c(NA, NA, 2),
         col = c("black", "darkorange", "steelblue"))
}

# --- Write figures and tables -------------------------------------------------------

make_outputs <- function(res, out_dir = "results") {
  fig_dir <- file.path(out_dir, "figures")
  tab_dir <- file.path(out_dir, "tables")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

  fits  <- res$lassen$fits
  pairs <- res$lassen$pairs
  
  # Figure 1: EC50 curve
  png(file.path(fig_dir, "ec50_curve.png"), width = 7, height = 5, units = "in", res = 300)
  plot_ec50(res$ec50)
  dev.off()
  
  # Figure 2: example Lassen plots at low, medium and high occupancy
  # (chosen among good fits, R2 >= 0.8, so they show the method clearly)
  ok  <- fits[!fits$qc_exclude & !is.na(fits$occupancy) & fits$r_squared >= 0.8, ]
  ex  <- ok[order(ok$occupancy), ][round(c(0.15, 0.5, 0.9) * nrow(ok)), ]
  png(file.path(fig_dir, "lassen_examples.png"), width = 10, height = 3.8, units = "in", res = 300)
  op <- par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3.5, 1))
  for (i in seq_len(nrow(ex))) {
    plot_lassen(ex[i, ], pairs, xlim = c(0, max(pairs$VT_base)),
                ylim = range(c(0, pairs$delta)))
  }
  par(op)
  dev.off()
  
  # Table 1: occupancy per scan, with exposure
  occ <- merge(fits, res$exposure, by = c("subject", "session"), all.x = TRUE)
  occ <- occ[order(occ$subject, occ$session),
             c("subject", "session", "dose_mg", "conc_ng_mL", "n_regions",
               "occupancy", "occupancy_se", "VND", "r_squared", "qc_exclude")]
  write.csv(occ, file.path(tab_dir, "occupancy_per_scan.csv"), row.names = FALSE)

  # Table 2: EC50 result
  write.csv(res$ec50$result, file.path(tab_dir, "ec50_result.csv"), row.names = FALSE)

  # Table 3: all QC flags
  write.csv(res$qc_flags, file.path(tab_dir, "qc_flags.csv"), row.names = FALSE)

  invisible(list.files(out_dir, recursive = TRUE))
}
