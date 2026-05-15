# =============================================================================
# main_replication.R
# -----------------------------------------------------------------------------
# Drives the SAVS TVP-VAR replication:
#   1. simulate a sparse TVP-VAR(p) with known coefficient paths and mask
#      (or REPLACE WITH YOUR DATA: load your own (T x M) panel),
#   2. estimate the model equation by equation with `sparse.SAVS`,
#   3. compute reduced-form coefficients, IRFs and forecasts,
#   4. plot with truth overlays:
#         - sparsity_pip.pdf     true mask vs posterior PIPs
#         - tvp_coefficients.pdf TVP coefficient paths (true vs estimated)
#         - irf_slices.pdf       IRFs at three time slices (true vs estimated)
#         - forecast_fan.pdf     forecast fan + realised future y
#
# Closely mirrors the real-data driver `estim.VAR.R`; main difference is that
# `Yraw` here comes from the DGP and the truth is overlaid on the figures.
# =============================================================================

rm(list = ls())

this_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)
setwd(this_dir)

suppressPackageStartupMessages({
  library(Rcpp)
  library(bvarsv)
  library(GIGrvg)
  library(mvtnorm)
  library(stochvol)
  library(MASS)
  library(ggplot2)
  library(reshape2)
})

source("auxilliary_functions.R")   # also sourceCpp's threshold_functions.cpp
source("ng_SAVS.R")
source("simulate_data.R")

dir.create("figures", showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Data
# -----------------------------------------------------------------------------
# REPLACE WITH YOUR DATA: build a (T x M) matrix `Yraw`. For the empirical
# exercise see `estim.VAR.R` which loads the small / medium / large FRED-QD
# panel from `datasetMNg.RData`.
set.seed(123)
sim       <- simulate_savs_tvp_var(T = 150, M = 3, p = 2,
                                   sparsity = 0.6, drift_sd = 0.01,
                                   seed = 42, fhorz = 8)
Yraw      <- sim$Y
Y_future  <- sim$Y_future
A_path_tr <- sim$A_path      # (T x K x M) true TVP lag coefficients
mask_true <- sim$mask        # (K x M) sparsity pattern

# Standardise (the original driver standardises before the sampler runs).
Y_sd      <- apply(Yraw, 2, sd)
Y_mu      <- apply(Yraw, 2, mean)
Yraw      <- apply(Yraw, 2, function(x) (x - mean(x)) / sd(x))

# -----------------------------------------------------------------------------
# 2. MCMC settings
# -----------------------------------------------------------------------------
prior.choice <- "horse"      # "horse", "LASSO", "SSVS", "NG", "DL"
nu           <- 2            # SAVS penalty exponent ( 1 / |beta|^nu )
nsave        <- 500
nburn        <- 500
horz.fcst    <- 8            # predictive horizon
h            <- 0            # in-sample window (no hold-out)
sl.marg      <- seq_len(3)
nhor         <- 20           # IRF horizon
laglen       <- 2

# -----------------------------------------------------------------------------
# 3. Build VAR regressors (mirrors estim.VAR.R)
# -----------------------------------------------------------------------------
Xraw <- mlag(Yraw, laglen)
Y    <- Yraw[(laglen + 1):(nrow(Yraw) - h), , drop = FALSE]
X    <- Xraw[(laglen + 1):(nrow(Xraw) - h), , drop = FALSE]

m <- ncol(Y); T_eff <- nrow(X); K <- ncol(X)

A.matrix   <- array(NA, c(nsave, T_eff, m, K))   # lag coefs per eq.
A0.matrix  <- array(0,  c(nsave, T_eff, m, m))   # contemporaneous coefs
H.matrix   <- array(NA, c(nsave, m, T_eff))      # log volatilities
PIP.matrix <- matrix(NA, K, m)                   # post. inclusion prob, last t

# -----------------------------------------------------------------------------
# 4. Equation-by-equation SAVS sampler
# -----------------------------------------------------------------------------
for (i in seq_len(m)) {
  sl.i <- if (i == 1) NULL else seq_len(i - 1)
  X.i  <- cbind(Y[, sl.i, drop = FALSE], X)
  Y.i  <- Y[, i, drop = FALSE]

  cat(sprintf("Fitting SAVS for equation %d / %d (%d regressors)...\n",
              i, m, ncol(X.i)))
  sim.1 <- sparse.SAVS(Y.i, X.i, h = 0, p = 0,
                       nsave = nsave, nburn = nburn, sv = TRUE,
                       prior.mean = c(0, 0), nu = nu,
                       shrinkage = prior.choice)

  if (length(sl.i) > 0)
    A0.matrix[, , i, sl.i] <- -1 * sim.1$A.thrsh[, , sl.i, 1]
  if (i == 1) {
    A.matrix[, , i, ] <- sim.1$A.thrsh[, , , 1]
  } else {
    A.matrix[, , i, ] <- sim.1$A.thrsh[, , -sl.i, 1]
  }
  H.matrix[, i, ] <- sim.1$hv

  # PIP at the last in-sample period
  last_t <- T_eff
  cols   <- if (length(sl.i) > 0) seq(length(sl.i) + 1, ncol(X.i)) else seq_len(ncol(X.i))
  PIP.matrix[, i] <- apply((sim.1$A.thrsh[, last_t, cols, 1] != 0) * 1, 2, mean)
}

# -----------------------------------------------------------------------------
# 5. Reduced form, IRFs, forecasts
# -----------------------------------------------------------------------------
A.red.array   <- array(NA, c(nsave, T_eff, K, m))
Sig.red.array <- array(NA, c(nsave, T_eff, m, m))
IRF.array     <- array(NA, c(nsave, m, m, nhor, T_eff))
pred.array    <- array(NA, c(nsave, length(sl.marg), horz.fcst))

for (irep in seq_len(nsave)) {
  for (t in seq_len(T_eff)) {
    A0 <- A0.matrix[irep, t, , ]; diag(A0) <- 1
    A0inv <- t(solve(A0))
    At    <- A0inv %*% A.matrix[irep, t, , ]
    A.red.array[irep, t, , ] <- t(At)

    Sig.red <- A0inv %*% diag(exp(H.matrix[irep, , t])) %*% t(A0inv)
    Sig.red.array[irep, t, , ] <- Sig.red

    PHI_array <- array(0, c(m, m, laglen))
    Att <- t(At)
    for (ss in seq_len(laglen))
      PHI_array[, , ss] <- t(Att[((ss - 1) * m + 1):(ss * m), ])
    IRF.array[irep, , , , t] <- impulsdtrf(PHI_array, t(chol(Sig.red)), nhor)
  }
  # Forecasts (predictive mean only; deliberately the same trick as estim.VAR.R)
  pred.sparse <- get.pred(
    zt = t(Xraw[nrow(Xraw), , drop = FALSE]),
    laglen = laglen, m = m, K = K, At = At, SIGMA = Sig.red,
    fhorz = horz.fcst, Yraw.full = Yraw, Xraw.full = Xraw,
    sl.joint = sl.marg
  )
  pred.array[irep, , ] <- pred.sparse$yhat
}

A.red.med  <- apply(A.red.array,  c(2, 3, 4), median)
IRF.median <- apply(IRF.array,    c(2, 3, 4, 5), median)
IRF.low    <- apply(IRF.array,    c(2, 3, 4, 5), quantile, 0.16, na.rm = TRUE)
IRF.high   <- apply(IRF.array,    c(2, 3, 4, 5), quantile, 0.84, na.rm = TRUE)
fcst.med   <- apply(pred.array,   c(2, 3), median)
fcst.lo16  <- apply(pred.array,   c(2, 3), quantile, 0.16, na.rm = TRUE)
fcst.hi84  <- apply(pred.array,   c(2, 3), quantile, 0.84, na.rm = TRUE)

# -----------------------------------------------------------------------------
# 6. TRUE objects from the DGP (on the standardised scale)
# -----------------------------------------------------------------------------
# Coefficient paths from the DGP are on the *original* scale. The model sees
# y on the standardised scale; rescale the truth accordingly: a coefficient
# on regressor x_j for response y_i transforms as A_ij * sd(x_j) / sd(y_i).
# All regressors are lags of standardised y, so sd(x_j) = 1 by construction
# of the standardisation; rescaling is then 1 / sd(y_i).
A_path_std <- A_path_tr * 0
for (t in seq_len(dim(A_path_tr)[1]))
  for (i in seq_len(m))
    A_path_std[t, , i] <- A_path_tr[t, , i] / Y_sd[i]   # K x 1 per (t, i)

# The estimator drops the first `laglen` observations; align truth to the
# same in-sample window.
A_path_eff <- A_path_std[(laglen + 1):dim(A_path_std)[1], , , drop = FALSE]

# True IRFs at the same slice times we plot below.
slice_times <- unique(round(seq(1, T_eff, length.out = 3)))
true_irf_slice <- array(NA, c(m, m, nhor, length(slice_times)))
for (ix in seq_along(slice_times)) {
  tt <- slice_times[ix]
  PHI_t <- array(0, c(m, m, laglen))
  At    <- A_path_eff[tt, , ]
  Att   <- At
  for (ss in seq_len(laglen))
    PHI_t[, , ss] <- t(Att[((ss - 1) * m + 1):(ss * m), ])
  # True residual covariance (rescaled by 1/sd_y on both sides).
  Sig_true <- diag(1 / Y_sd) %*% sim$L %*% diag(sim$sig[tt + laglen, ]^2) %*%
              t(sim$L) %*% diag(1 / Y_sd)
  true_irf_slice[, , , ix] <- impulsdtrf(PHI_t, t(chol(Sig_true)), nhor)
}

# -----------------------------------------------------------------------------
# 7. Plot: PIP vs true sparsity mask
# -----------------------------------------------------------------------------
# Reorder PIPs so the rows correspond to the lag structure used by mlag().
pip_df <- data.frame(
  regressor = rep(paste0("x", seq_len(K)), times = m),
  equation  = rep(paste0("eq", seq_len(m)),  each  = K),
  pip       = as.numeric(PIP.matrix),
  active    = as.numeric(mask_true)
)
pip_df$regressor <- factor(pip_df$regressor, levels = paste0("x", seq_len(K)))

p_pip <- ggplot(pip_df, aes(x = regressor, y = equation)) +
  geom_tile(aes(fill = pip), colour = "grey60") +
  geom_point(data = subset(pip_df, active == 1),
             colour = "red3", shape = 4, size = 3, stroke = 1.1) +
  scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1)) +
  labs(title  = "Posterior inclusion probabilities (SAVS, last in-sample t)",
       subtitle = "Red 'x' marks regressors that are active in the DGP.",
       x = "regressor (lag x variable)", y = "equation") +
  theme_minimal(base_size = 11)

ggsave("figures/sparsity_pip.pdf", p_pip, width = 8, height = 4)
ggsave("figures/sparsity_pip.png", p_pip, width = 8, height = 4, dpi = 150)

# -----------------------------------------------------------------------------
# 8. Plot: TVP coefficient paths (truth vs posterior)
# -----------------------------------------------------------------------------
# Pick a handful of (regressor, equation) pairs to display: the own-lag-1
# diagonals (always active by construction) plus the first off-diagonal of
# equation 1.
pick <- list()
for (i in seq_len(m)) pick[[length(pick) + 1]] <- c(i, i)   # own lag-1
if (m >= 2) pick[[length(pick) + 1]] <- c(2, 1)             # x2 -> eq1

coef_df <- do.call(rbind, lapply(pick, function(ij) {
  reg <- ij[1]; eq <- ij[2]
  med <- apply(A.red.array[, , reg, eq], 2, median)
  lo  <- apply(A.red.array[, , reg, eq], 2, quantile, 0.16, na.rm = TRUE)
  hi  <- apply(A.red.array[, , reg, eq], 2, quantile, 0.84, na.rm = TRUE)
  data.frame(time = seq_len(T_eff),
             panel = sprintf("eq%d <- x%d", eq, reg),
             median = med, lo = lo, hi = hi,
             true = A_path_eff[, reg, eq])
}))

p_coef <- ggplot(coef_df, aes(x = time)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.25) +
  geom_line(aes(y = median), colour = "steelblue4", linewidth = 0.7) +
  geom_line(aes(y = true), colour = "red3", linetype = "dashed",
            linewidth = 0.7) +
  facet_wrap(~ panel, scales = "free_y") +
  labs(title  = "Time-varying coefficients: posterior vs truth",
       subtitle = "Blue: posterior median + 68% band.  Red dashed: true DGP path.",
       x = "time", y = "coefficient") +
  theme_minimal(base_size = 11)

ggsave("figures/tvp_coefficients.pdf", p_coef, width = 9, height = 5)
ggsave("figures/tvp_coefficients.png", p_coef, width = 9, height = 5, dpi = 150)

# -----------------------------------------------------------------------------
# 9. Plot: IRFs at three slice times (truth vs posterior)
# -----------------------------------------------------------------------------
shock_var <- 1L
irf_slice_df <- do.call(rbind, lapply(seq_along(slice_times), function(ix) {
  tt <- slice_times[ix]
  do.call(rbind, lapply(seq_len(m), function(v) {
    data.frame(
      time_idx = paste0("t = ", tt),
      var      = paste0("y", v),
      h        = seq_len(nhor),
      median   = IRF.median[v, shock_var, , tt],
      lo       = IRF.low[   v, shock_var, , tt],
      hi       = IRF.high[  v, shock_var, , tt],
      true     = true_irf_slice[v, shock_var, , ix]
    )
  }))
}))

p_irf <- ggplot(irf_slice_df, aes(x = h, y = median)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25, fill = "steelblue") +
  geom_line(colour = "steelblue4", linewidth = 0.7) +
  geom_line(aes(y = true), colour = "red3", linetype = "dashed",
            linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey50") +
  facet_grid(var ~ time_idx, scales = "free_y") +
  labs(title  = sprintf("IRFs to a one-s.d. shock in y%d (Cholesky)", shock_var),
       subtitle = "Blue: posterior median + 68%.  Red dashed: true IRF (DGP).",
       x = "horizon", y = "response") +
  theme_minimal(base_size = 11)

ggsave("figures/irf_slices.pdf", p_irf, width = 9, height = 6)
ggsave("figures/irf_slices.png", p_irf, width = 9, height = 6, dpi = 150)

# -----------------------------------------------------------------------------
# 10. Plot: forecast fan with realised future y
# -----------------------------------------------------------------------------
# True future on the standardised scale.
Y_future_std <- t((t(Y_future) - Y_mu) / Y_sd)
fan_df <- do.call(rbind, lapply(seq_along(sl.marg), function(k) {
  data.frame(
    h        = seq_len(horz.fcst),
    variable = colnames(Yraw)[sl.marg[k]],
    median   = fcst.med[k, ],
    lo       = fcst.lo16[k, ],
    hi       = fcst.hi84[k, ],
    true     = Y_future_std[seq_len(horz.fcst), sl.marg[k]]
  )
}))

p_fan <- ggplot(fan_df, aes(x = h)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.30) +
  geom_line(aes(y = median),  colour = "steelblue4", linewidth = 0.7) +
  geom_line(aes(y = true),    colour = "red3", linetype = "dashed",
            linewidth = 0.7) +
  geom_point(aes(y = true),   colour = "red3", size = 1.3) +
  facet_wrap(~ variable, scales = "free_y") +
  labs(title  = "Recursive predictive density (standardised scale)",
       subtitle = "Blue: posterior median + 68%.  Red dashed: realised future y from the DGP.",
       x = "forecast horizon", y = "response") +
  theme_minimal(base_size = 11)

ggsave("figures/forecast_fan.pdf", p_fan, width = 9, height = 5)
ggsave("figures/forecast_fan.png", p_fan, width = 9, height = 5, dpi = 150)

cat("\nDone. Figures written to ./figures/.\n")
