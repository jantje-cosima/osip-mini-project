# Estimate VT (with SE) for every subject x session x region by fitting a
# two-tissue compartment model (2TCM) to the simulated TACs, using the
# measured (noisy) metabolite-corrected plasma samples as input.
#
# Run from the project root:  Rscript scripts/fit_vt.R
# Reads data/raw/, writes data/processed/vt_estimates_refit.csv.

in_dir  <- file.path("data", "raw")
out_dir <- file.path("data", "processed")

tacs  <- read.csv(file.path(in_dir, "tacs.csv"))
blood <- read.csv(file.path(in_dir, "plasma_input.csv"))

# Analysis assumptions
vB_fit <- 0.05                  # blood volume fraction, fixed in the fit
tracer_half_life_min <- 20.4    # 11C, used for frame weights
dt_min <- 1 / 60                # integration step for the model (1 s)
lambda <- log(2) / tracer_half_life_min
lower <- c(K1 = 0.01, k2 = 0.005, k3 = 0,   k4 = 0.005)
upper <- c(K1 = 2,    k2 = 2,     k3 = 2,   k4 = 2)
starts <- list(c(0.4, 0.15, 0.10, 0.06),
               c(0.3, 0.10, 0.03, 0.03),
               c(0.5, 0.30, 0.30, 0.15))

# --- Model -------------------------------------------------------------------

conv_exp <- function(cp, a, dt) {
  e <- exp(-a * dt)
  x <- dt / 2 * (cp + e * c(0, head(cp, -1)))
  as.numeric(stats::filter(x, e, method = "recursive"))
}

tcm2 <- function(cp, K1, k2, k3, k4, dt) {
  s <- k2 + k3 + k4
  r <- sqrt(max(s^2 - 4 * k2 * k4, 0))
  r <- max(r, 1e-8)               # avoid a1 == a2 (k3 = 0, k2 = k4)
  a1 <- (s - r) / 2
  a2 <- (s + r) / 2
  K1 / (a2 - a1) * ((k3 + k4 - a1) * conv_exp(cp, a1, dt) +
                      (a2 - k3 - k4) * conv_exp(cp, a2, dt))
}

# Frame means of a fine-grid curve; fi holds first/last grid index per frame
frame_mean <- function(x, fi) {
  cs <- c(0, cumsum(x))
  (cs[fi$last + 1] - cs[fi$first]) / (fi$last - fi$first + 1)
}

# Frame-averaged model TAC for parameters p on the fine grid
model_tac <- function(p, cp, fi) {
  pet <- (1 - vB_fit) * tcm2(cp, p[1], p[2], p[3], p[4], dt_min) + vB_fit * cp
  frame_mean(pet, fi)
}

vt_fun <- function(p) p[1] / p[2] * (1 + p[3] / p[4])

# Numerical gradient (central differences)
num_grad <- function(f, p, h = 1e-6) {
  sapply(seq_along(p), function(i) {
    d <- replace(numeric(length(p)), i, h * max(abs(p[i]), 1e-3))
    (f(p + d) - f(p - d)) / (2 * d[i])
  })
}

fit_one <- function(y, w, cp, fi) {
  wrss <- function(p) {
    v <- sum(w * (y - model_tac(p, cp, fi))^2)
    if (is.finite(v)) v else 1e10
  }
  fits <- lapply(starts, function(s)
    tryCatch(optim(s, wrss, method = "L-BFGS-B", lower = lower, upper = upper,
                   control = list(maxit = 500, parscale = s, factr = 1e9)),
             error = function(e) NULL))
  fits <- Filter(Negate(is.null), fits)
  best <- fits[[which.min(sapply(fits, `[[`, "value"))]]
  p <- best$par

  n <- length(y); k <- length(p)
  sigma2 <- best$value / (n - k)
  # Covariance from the Jacobian of the weighted residuals: s2 * (J'WJ)^-1.
  # Parameters stuck at a bound are treated as fixed.
  at_b <- abs(p - lower) < 1e-6 | abs(p - upper) < 1e-6
  if (at_b[3]) at_b[4] <- TRUE    # k3 = 0: k4 drops out of the model
  free <- which(!at_b)
  J <- sapply(free, function(i) {
    d <- replace(numeric(k), i, 1e-6 * max(abs(p[i]), 1e-3))
    (model_tac(p + d, cp, fi) - model_tac(p - d, cp, fi)) / (2 * d[i])
  })
  cov_free <- tryCatch(sigma2 * solve(t(J) %*% (w * J)), error = function(e) NULL)

  vt <- vt_fun(p)
  vt_se <- if (is.null(cov_free)) NA else {
    g <- num_grad(vt_fun, p)[free]
    sqrt(max(as.numeric(t(g) %*% cov_free %*% g), 0))
  }
  at_bound <- any(at_b)
  data.frame(K1 = p[1], k2 = p[2], k3 = p[3], k4 = p[4],
             VT = vt, VT_SE = vt_se, VT_SE_pct = 100 * vt_se / vt,
             wrss = best$value, converged = best$convergence == 0,
             at_bound = at_bound)
}

# --- Fit all TACs ------------------------------------------------------------
keys <- unique(tacs[, c("subject", "session")])

fit_scan <- function(i) {
  s <- keys$subject[i]; ses <- keys$session[i]
  b <- blood[blood$subject == s & blood$session == ses, ]
  tac_s <- tacs[tacs$subject == s & tacs$session == ses, ]
  tac_s <- tac_s[order(tac_s$region, tac_s$frame), ]

  fr <- unique(tac_s[, c("frame", "frame_start_min", "frame_end_min")])
  fr <- fr[order(fr$frame), ]
  t_grid <- seq(0, max(fr$frame_end_min), by = dt_min)
  idx <- findInterval(t_grid, fr$frame_start_min)
  idx[t_grid >= max(fr$frame_end_min)] <- NA
  fi <- list(first = match(seq_len(nrow(fr)), idx),
             last  = length(idx) + 1 - match(seq_len(nrow(fr)), rev(idx)))
  dur <- fr$frame_end_min - fr$frame_start_min
  mid <- (fr$frame_start_min + fr$frame_end_min) / 2

  # Input: linear interpolation of the measured samples, zero before the first
  cp <- approx(c(0, b$time_min), c(0, b$plasma_parent_kBq_mL), xout = t_grid,
               rule = 2)$y
  cp <- pmax(cp, 0)

  out <- lapply(unique(tac_s$region), function(r) {
    y <- tac_s$activity_kBq_mL[tac_s$region == r]
    # Weights ~ 1 / variance of the noise model
    w <- dur * exp(-lambda * mid) / pmax(y, 0.05 * max(y))
    w <- w / mean(w)
    cbind(data.frame(subject = s, session = ses, region = r),
          fit_one(y, w, cp, fi))
  })
  do.call(rbind, out)
}

# One scan per worker (base R 'parallel', works on Windows)
n_workers <- max(1, min(parallel::detectCores() - 1, nrow(keys)))
cl <- parallel::makeCluster(n_workers)
parallel::clusterExport(cl, ls())
res <- parallel::parLapply(cl, seq_len(nrow(keys)), fit_scan)
parallel::stopCluster(cl)

res <- do.call(rbind, res)
num <- sapply(res, is.numeric)
res[num] <- lapply(res[num], signif, 5)
write.csv(res, file.path(out_dir, "vt_estimates_refit.csv"), row.names = FALSE)
cat("Wrote", nrow(res), "VT estimates to", out_dir, "\n")
