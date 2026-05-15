# ============================================================================
# Main replication driver for the TVP-SVD model of
# Hauzenberger, Huber, Koop & Onorante (2021, JBES).
#
# This script
#   1. loads (synthetic) data via simulate_data();
#   2. estimates the TVP-WN-SVD model with sparse-mixture g-prior;
#   3. summarises posterior of the time-varying coefficients;
#   4. computes an h-step-ahead predictive density;
#   5. produces diagnostic plots in figures/.
#
# Run with default settings (nsave = 2000, nburn = 2000):
#   Rscript main_replication.R
#
# Or call directly:
#   source("main_replication.R")
# ============================================================================

rm(list = ls())
set.seed(20210501)

# ---- working directory (works from inside the replication folder) --------
this_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)),
                     error = function(e) getwd())
setwd(this_dir)

# ---- packages -------------------------------------------------------------
suppressPackageStartupMessages({
  require(coda); require(GIGrvg); require(MASS); require(Matrix)
  require(shrinkTVP); require(stochvol); require(zoo); require(mvtnorm)
  require(ggplot2); require(reshape2); require(bayesm)
})

dir.create("figures", showWarnings = FALSE)

# ---- data ----------------------------------------------------------------
source("simulate_data.R")
sim   <- simulate_data(T = 200, M_x = 4, sigma = 0.5, seed = 1234)
y     <- sim$y
X     <- sim$X
beta_true <- sim$beta_true
T <- nrow(y); K <- ncol(X)

# ---- estimator ------------------------------------------------------------
source("tvpsvd_estim.R")

# MCMC settings (override on the command line for quick smoke tests)
args <- commandArgs(trailingOnly = TRUE)
nsave <- if (length(args) >= 1) as.integer(args[1]) else 2000
nburn <- if (length(args) >= 2) as.integer(args[2]) else 2000
thin  <- if (length(args) >= 3) as.integer(args[3]) else 1

# Forecast horizon for h-step-ahead predictive density
h <- 4

model.setup <- list(
  model.type     = "TVP",
  tvp.type       = "WN-SVD",
  set.prior.time = "g",
  pooling        = TRUE,
  ub.sc          = 0.02,   # upper bound on innovation variance scaling
                            # (xi <= ub.sc * T / K^2).  For the synthetic
                            # data with sigma = 0.5 noise we use 0.02; for
                            # tight macro applications (e.g. the NKPC
                            # exercise in the paper) the authors use 0.005.
  shrnk.type     = "ng",
  a_tau          = 0.1,
  sv             = TRUE,
  fast           = TRUE,
  mix.model.setup = list(
    G            = 30,
    e0           = 0.1,
    sample.e0    = TRUE,
    perm.sampler = TRUE,
    common.mean  = TRUE
  )
)

cat("\nEstimating TVP-WN-SVD model ...\n")
cat("  T  =", T, "   K =", K, "\n")
cat("  nsave =", nsave, "   nburn =", nburn, "   thin =", thin, "\n\n")

Regr.obj <- TVPSVD_estim(y = y, X = X,
                         nsave = nsave, nburn = nburn, thin = thin,
                         model.setup = model.setup,
                         verbose = TRUE)

cat("\nMCMC complete.  Elapsed:",
    round(Regr.obj$time.min, 2), "minutes.\n\n")

# ---- posterior summaries of beta_t ---------------------------------------
b.full <- Regr.obj$b.full.store      # nsave x T x K
b.med  <- apply(b.full, c(2, 3), median)
b.lo16 <- apply(b.full, c(2, 3), quantile, 0.16)
b.hi84 <- apply(b.full, c(2, 3), quantile, 0.84)
b.lo05 <- apply(b.full, c(2, 3), quantile, 0.05)
b.hi95 <- apply(b.full, c(2, 3), quantile, 0.95)

# ---- h-step-ahead predictive density -------------------------------------
# Pseudo-out-of-sample: condition on the last in-sample observation and
# carry the most recent beta_t forward (WN-SVD: beta does not drift) while
# adding the posterior variance from the latent volatility process.
nsave_eff <- dim(b.full)[1]
y_pred <- matrix(NA, nsave_eff, h)
# Generate a forward draw of regressors: use bootstrap from the empirical
# distribution of X to keep things data-driven.
for (s in 1:nsave_eff){
  beta_T  <- b.full[s, T, ]
  sigma_T <- Regr.obj$eht.store[s, T]
  for (hh in 1:h){
    x_future <- X[sample.int(T, 1), ]
    y_pred[s, hh] <- sum(x_future*beta_T) + rnorm(1, 0, sigma_T)
  }
}

pred.med  <- apply(y_pred, 2, median)
pred.lo16 <- apply(y_pred, 2, quantile, 0.16)
pred.hi84 <- apply(y_pred, 2, quantile, 0.84)
pred.lo05 <- apply(y_pred, 2, quantile, 0.05)
pred.hi95 <- apply(y_pred, 2, quantile, 0.95)

# ---- plot 1: TVP coefficients with credible bands vs. true ---------------
plot_df <- do.call(rbind, lapply(1:K, function(k){
  data.frame(time = 1:T,
             coef = paste0("beta_", k - 1),
             true = beta_true[, k],
             med  = b.med[, k],
             lo16 = b.lo16[, k], hi84 = b.hi84[, k],
             lo05 = b.lo05[, k], hi95 = b.hi95[, k])
}))
plot_df$coef <- factor(plot_df$coef,
                       levels = paste0("beta_", 0:(K - 1)))

p1 <- ggplot(plot_df, aes(x = time)) +
  geom_ribbon(aes(ymin = lo05, ymax = hi95), fill = "navyblue", alpha = 0.15) +
  geom_ribbon(aes(ymin = lo16, ymax = hi84), fill = "navyblue", alpha = 0.30) +
  geom_line(aes(y = med),  colour = "navyblue", linewidth = 0.9) +
  geom_line(aes(y = true), colour = "red3",    linewidth = 0.7,
            linetype = "dashed") +
  facet_wrap(~ coef, scales = "free_y", ncol = 2) +
  labs(title = "Posterior of time-varying coefficients (TVP-WN-SVD)",
       subtitle = "Blue: median with 68%/90% bands.  Red dashed: true path.",
       x = "time", y = NULL) +
  theme_bw(base_size = 12)

ggsave("figures/tvp_coefficients.pdf", p1, width = 9, height = 6)
ggsave("figures/tvp_coefficients.png", p1, width = 9, height = 6, dpi = 120)

# ---- plot 2: predictive fan chart ----------------------------------------
# TRUE future y: extend the DGP one regime forward (beta is piecewise constant
# after the last break, so the truth holds beta_true[T,] going forward). For
# each future step bootstrap an x and add the true sigma noise.
set.seed(20210501 + 7)
n_true_paths <- 200
true_future <- matrix(NA, n_true_paths, h)
for (s in seq_len(n_true_paths)) {
  for (hh in seq_len(h)) {
    x_future <- X[sample.int(T, 1), ]
    true_future[s, hh] <- sum(x_future * beta_true[T, ]) +
                         rnorm(1, 0, sim$sigma)
  }
}
true.med  <- apply(true_future, 2, median)
true.lo16 <- apply(true_future, 2, quantile, 0.16)
true.hi84 <- apply(true_future, 2, quantile, 0.84)

fan_df <- data.frame(
  time = (T + 1):(T + h),
  med  = pred.med,
  lo16 = pred.lo16, hi84 = pred.hi84,
  lo05 = pred.lo05, hi95 = pred.hi95,
  true.med  = true.med,
  true.lo16 = true.lo16,
  true.hi84 = true.hi84)

hist_df <- data.frame(time = 1:T, y = as.numeric(y))

p2 <- ggplot() +
  geom_line(data = hist_df, aes(x = time, y = y),
            colour = "grey30") +
  geom_ribbon(data = fan_df, aes(x = time, ymin = lo05, ymax = hi95),
              fill = "navyblue", alpha = 0.20) +
  geom_ribbon(data = fan_df, aes(x = time, ymin = lo16, ymax = hi84),
              fill = "navyblue", alpha = 0.40) +
  geom_line(data = fan_df, aes(x = time, y = med),
            colour = "navyblue", linewidth = 0.9) +
  # truth from the DGP (median + 68% band)
  geom_ribbon(data = fan_df, aes(x = time, ymin = true.lo16, ymax = true.hi84),
              fill = "red3", alpha = 0.15) +
  geom_line(data = fan_df, aes(x = time, y = true.med),
            colour = "red3", linetype = "dashed", linewidth = 0.7) +
  labs(title = paste0("h = ", h,
                      "-step-ahead predictive density (TVP-WN-SVD)"),
       subtitle = "Blue: model median + 68/90% bands. Red dashed: DGP median + 68% band.",
       x = "time", y = "y") +
  theme_bw(base_size = 12)

ggsave("figures/predictive_fan.pdf", p2, width = 8, height = 5)
ggsave("figures/predictive_fan.png", p2, width = 8, height = 5, dpi = 120)

# ---- console summary ------------------------------------------------------
cat("\nPosterior summary at break points (true vs. median):\n")
brks <- c(50, 70, 100, 140, 180)
for (bk in brks){
  cat(sprintf("  t = %3d   true = [%s]   median = [%s]\n",
              bk,
              paste(sprintf("%6.2f", beta_true[bk, ]), collapse = " "),
              paste(sprintf("%6.2f", b.med[bk, ]),     collapse = " ")))
}

cat("\nFigures written to figures/.\n")
cat("Replication complete.\n")
