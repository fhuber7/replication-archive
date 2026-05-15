# =============================================================================
# main_replication.R
# -----------------------------------------------------------------------------
# End-to-end driver for the mixBART / flexBART replication:
#   1. simulate (or load) a small multivariate time series,
#   2. fit a Bayesian VAR with BART-based conditional means and a factor
#      stochastic-volatility shock structure (`flexBART` in
#      `estimation_functions/flexBART.R`),
#   3. plot the one-step and full-horizon predictive densities for each
#      variable, and the posterior median volatility factors.
#
# Notes on cost:
#   * Each Gibbs sweep grows trees for each of the M equations (BART) and
#     does an exact factor-SV update per period; runtime is roughly
#     O((nburn + nsave) * (M * tree_cost + T * Q^2)).
#   * Defaults below (nsave = nburn = 500) are a minimal exploratory run.
#     Production runs in the underlying paper use several thousand draws.
# =============================================================================

rm(list = ls())

this_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)
setwd(this_dir)

suppressPackageStartupMessages({
  library(invgamma)
  library(stochvol)
  library(factorstochvol)
  library(MASS)
  library(dbarts)
  library(mvtnorm)
  library(ggplot2)
  library(reshape2)
})

source("aux_func.R")
source("estimation_functions/flexBART.R")
source("estimation_functions/fullBART.R")
source("estimation_functions/errorBART.R")
source("simulate_data.R")

dir.create("figures", showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
# REPLACE WITH YOUR DATA: assign a (T x M) numeric matrix to `Ytrn`. The paper
# applies the model to a quarterly US macro panel (GDP, UNR, CPI, FFR, ...).
set.seed(123)
sim       <- simulate_mixbart_data(T = 200, M = 5, p = 2, fhorz = 12, seed = 42)
Ytrn      <- sim$Y
Y_future  <- sim$Y_future       # truth retained for overlay below
M         <- ncol(Ytrn)
labs      <- sim$var.names

# -----------------------------------------------------------------------------
# 2. Sampler configuration
# -----------------------------------------------------------------------------
# `model`: "mixBART" -> mixture-of-BART conditional means (the headline model),
#          "BART"    -> standard BART means per equation,
#          switch to `fullBART(...)` instead of `flexBART(...)` to fit the
#          full-VAR variant in estimation_functions/fullBART.R.
# `sv`   : "homo" (homoskedastic), "SV" (factor SV), or "heteroBART".
# `fc.approx`: "exact" or an approximation flag handled in flexBART.R.
sl.mod     <- "mixBART"
sl.sv.mod  <- "homo"
sl.fc.type <- "exact"
fhorz      <- 12

nsave <- 500
nburn <- 500

# Optional prior mean for the AR(1) part of each equation. 0 by default; set
# the diagonal towards 1 for series you believe behave near unit-root.
prior.mean <- matrix(0, M, M)

# -----------------------------------------------------------------------------
# 3. Fit
# -----------------------------------------------------------------------------
fit <- flexBART(Ytrn,
                nburn     = nburn,
                nsave     = nsave,
                thinfac   = 1,
                prior     = "HS",
                prior.sig = c(3, 0.9),
                model     = sl.mod,
                sv        = sl.sv.mod,
                fc.approx = sl.fc.type,
                fhorz     = fhorz,
                quiet     = FALSE,
                pr.mean   = prior.mean)

# -----------------------------------------------------------------------------
# 4. Posterior summaries
# -----------------------------------------------------------------------------
# fit$fcst has dim (nsave, fhorz, M). Compute pointwise 16/50/84 quantiles.
fcst_q <- apply(fit$fcst, c(2, 3), quantile,
                probs = c(0.16, 0.5, 0.84), na.rm = TRUE)
# fcst_q dim: 3 x fhorz x M
fcst_df <- data.frame()
for (m in seq_len(M)) {
  fcst_df <- rbind(fcst_df, data.frame(
    variable = labs[m],
    horizon  = seq_len(fhorz),
    lo       = fcst_q[1, , m],
    med      = fcst_q[2, , m],
    hi       = fcst_q[3, , m],
    true     = Y_future[seq_len(fhorz), m]
  ))
}

p_fan <- ggplot(fcst_df, aes(x = horizon, y = med)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#1f77b4", alpha = 0.25) +
  geom_line(colour = "#1f77b4", linewidth = 0.8) +
  geom_line(aes(y = true), colour = "red3", linetype = "dashed",
            linewidth = 0.7) +
  geom_point(aes(y = true), colour = "red3", size = 1.3) +
  facet_wrap(~ variable, scales = "free_y") +
  labs(title = paste0("Predictive density (16/50/84%) - model = ", sl.mod,
                      ", sv = ", sl.sv.mod),
       subtitle = "Blue: posterior median + 68%.  Red dashed: realised future y from the DGP.",
       x = "Forecast horizon", y = "Response (unstandardised)") +
  theme_minimal(base_size = 11)

ggsave("figures/forecast_fan.pdf", p_fan, width = 9, height = 5)
ggsave("figures/forecast_fan.png", p_fan, width = 9, height = 5, dpi = 150)

# --- One-step-ahead predictive density for variable 1 -----------------------
y1_h1 <- fit$fcst[, 1, 1]
pdf("figures/onestep_density_y1.pdf", width = 6, height = 4)
plot(density(y1_h1, na.rm = TRUE),
     main = paste0("One-step-ahead predictive density: ", labs[1]),
     xlab = labs[1], lwd = 2, col = "#1f77b4")
abline(v = 0, col = "grey50", lty = 2)
dev.off()

# -----------------------------------------------------------------------------
# 5. Persist the fit
# -----------------------------------------------------------------------------
saveRDS(fit, file = "figures/mixbart_fit.rds")

cat("\nDone. Figures written to ./figures/.\n")
