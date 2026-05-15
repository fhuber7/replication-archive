# ============================================================================
# main_replication.R
#
# Headline replication driver for the Subspace Shrinkage Conjugate VAR of
# Huber & Koop (Journal of Applied Econometrics).
#
# This script
#   (1) loads synthetic factor-driven data (REPLACE WITH YOUR DATA for real
#       applications),
#   (2) sets up Minnesota dummy-observation matrices via get.dum(),
#   (3) calls subVAR.0() on a grid of (omega, q, theta) hyperparameters and
#       integrates them out via marginal-likelihood weights,
#   (4) summarises the posterior of the VAR coefficients, recursive Cholesky
#       IRFs (with credible bands) and h-step ahead recursive forecasts,
#   (5) writes plots to ./figures/ :
#         - omega_q_heatmap.pdf  : posterior weights over (omega, q)
#         - irfs.pdf             : Cholesky IRFs with 68% bands
#         - fanchart.pdf         : fan chart of forecasts for variable 1
#
# Reduce / enlarge the grids and nsave/nburn below to trade off run-time vs
# Monte Carlo error. The defaults are small to be quick on a laptop.
# ============================================================================

# ----- working directory -----
this_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(this_dir) || !nzchar(this_dir)) this_dir <- getwd()
setwd(this_dir)

# ----- dependencies -----
suppressPackageStartupMessages({
  library(Matrix)
  library(MASS)
  library(mvtnorm)
  library(ggplot2)
  library(reshape2)
  library(zoo)
  library(scales)
})

source("subspacevar_estim.R")
source("simulate_data.R")

dir.create("figures", showWarnings = FALSE)

# ============================================================================
# 1) Data
# ============================================================================
# REPLACE WITH YOUR DATA: load a T x M numeric matrix in `Y`.
Y      <- make_synthetic_data(M = 8, T = 200, q_true = 2, seed = 1)
truth  <- attr(Y, "truth")
q_true <- truth$q_true
M      <- ncol(Y); T_obs <- nrow(Y)
cat(sprintf("Data: T = %d, M = %d (q_true = %d)\n", T_obs, M, q_true))

# ============================================================================
# 2) Hyperparameters / grids
# ============================================================================
p          <- 1                     # VAR lag order
nsave      <- 30                    # posterior draws (small for demo)
nburn      <- 30
fhor       <- 8                     # forecast horizon for variable 1
irfhor     <- 16                    # IRF horizon

grid.omega <- seq(0.05, 0.95, length.out = 5)         # 5 omegas
grid.q     <- c(1, 2, 4)                              # 3 subspace ranks
grid.theta <- c(0.05, 0.1, 0.5)                       # 3 Minnesota tightness

a.omega    <- 1
b.omega    <- 1
theta0     <- 0.1                                     # init / default theta

# ============================================================================
# 3) Estimation
# ============================================================================
set.seed(42)
fit <- subVAR.0(
  nsave        = nsave,
  nburn        = nburn,
  Y            = Y,
  p            = p,
  prior        = "Minn",
  sample.omega = TRUE,
  grid.omega   = grid.omega,
  grid.q       = grid.q,
  grid.theta   = grid.theta,
  theta        = theta0,
  a.omega      = a.omega,
  b.omega      = b.omega,
  a.mean       = 0,
  fhor         = fhor,
  irfhor       = irfhor,
  verbose      = FALSE
)

cat(sprintf(
  "Posterior summaries: median q = %.1f, mean omega = %.2f, median theta = %.3f\n",
  median(fit$q_posterior),
  mean(fit$omega_posterior),
  median(fit$theta_posterior)
))

# ============================================================================
# 4) Cholesky IRFs with 68% credible bands
# ============================================================================
irf_lo  <- apply(fit$IRF_posterior, c(2, 3, 4), quantile, probs = 0.16)
irf_md  <- apply(fit$IRF_posterior, c(2, 3, 4), quantile, probs = 0.50)
irf_hi  <- apply(fit$IRF_posterior, c(2, 3, 4), quantile, probs = 0.84)

irf_df <- do.call(rbind, lapply(seq_len(M), function(resp)
  do.call(rbind, lapply(seq_len(M), function(shock) {
    data.frame(
      horizon  = 0:(irfhor - 1),
      response = factor(paste0("y", resp), levels = paste0("y", 1:M)),
      shock    = factor(paste0("y", shock), levels = paste0("y", 1:M)),
      lo       = irf_lo[resp, shock, ],
      md       = irf_md[resp, shock, ],
      hi       = irf_hi[resp, shock, ]
    )
  }))
))

p_irf <- ggplot(irf_df, aes(x = horizon, y = md)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25, fill = "steelblue") +
  geom_line(colour = "steelblue4", linewidth = 0.6) +
  geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
  facet_grid(response ~ shock, scales = "free_y") +
  labs(x = "horizon", y = "response",
       title = "Cholesky IRFs (68% credible bands)") +
  theme_bw(base_size = 8)

ggsave("figures/irfs.pdf", p_irf, width = 9, height = 7)

# ============================================================================
# 5) Heatmap of posterior weights over (omega, q)
#    Marginalised over theta.
# ============================================================================
grid <- fit$grid
w    <- as.numeric(fit$weights)

heat_df <- aggregate(w, by = list(omega = grid$omega, q = grid$q), FUN = sum)
names(heat_df)[3] <- "weight"
heat_df$omega <- factor(round(heat_df$omega, 3))
heat_df$q     <- factor(heat_df$q)

p_heat <- ggplot(heat_df, aes(x = omega, y = q, fill = weight)) +
  geom_tile(colour = "white") +
  scale_fill_gradient(low = "white", high = "firebrick",
                      labels = scales::percent_format(accuracy = 1)) +
  geom_hline(aes(yintercept = which(levels(heat_df$q) == as.character(q_true))),
             colour = "blue", linetype = "dashed", linewidth = 0.8) +
  labs(x = expression(omega), y = "q (subspace dim.)",
       title = sprintf("Posterior weights over (omega, q)   |   q_true = %d (blue dashed)",
                       q_true)) +
  theme_minimal(base_size = 11)

ggsave("figures/omega_q_heatmap.pdf", p_heat, width = 6, height = 4)

# ============================================================================
# 6) Fan chart for variable 1
# ============================================================================
pred1 <- fit$pred_density[, 1, , drop = FALSE]       # nsave x 1 x fhor
pred1 <- matrix(pred1, nrow = nsave, ncol = fhor)

quants <- t(apply(pred1, 2, quantile,
                  probs = c(0.05, 0.16, 0.5, 0.84, 0.95)))
hist_last <- tail(Y[, 1], 24)
fan_df <- data.frame(
  time = c(seq_along(hist_last),
           length(hist_last) + seq_len(fhor)),
  type = c(rep("hist", length(hist_last)), rep("fcst", fhor)),
  q05  = c(hist_last, quants[, 1]),
  q16  = c(hist_last, quants[, 2]),
  q50  = c(hist_last, quants[, 3]),
  q84  = c(hist_last, quants[, 4]),
  q95  = c(hist_last, quants[, 5])
)

p_fan <- ggplot(fan_df, aes(x = time, y = q50)) +
  geom_ribbon(data = subset(fan_df, type == "fcst"),
              aes(ymin = q05, ymax = q95), fill = "steelblue", alpha = 0.15) +
  geom_ribbon(data = subset(fan_df, type == "fcst"),
              aes(ymin = q16, ymax = q84), fill = "steelblue", alpha = 0.35) +
  geom_line(colour = "black", linewidth = 0.6) +
  geom_vline(xintercept = length(hist_last) + 0.5,
             colour = "grey50", linetype = "dashed") +
  labs(x = "time", y = "y1",
       title = sprintf("Recursive forecast, variable 1 (h = 1..%d)", fhor)) +
  theme_bw(base_size = 11)

ggsave("figures/fanchart.pdf", p_fan, width = 7, height = 4)

# ============================================================================
# 7) Save posterior for downstream use
# ============================================================================
save(fit, Y, file = "subspaceVAR_fit.RData")

cat("Done. Figures written to ./figures/, posterior saved to ",
    "subspaceVAR_fit.RData\n", sep = "")
