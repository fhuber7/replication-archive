# =============================================================================
# main_replication.R
# -----------------------------------------------------------------------------
# End-to-end driver for the Subspace-Inflation-BNP replication:
#   1. simulate (or load) a single-equation (y, X) training set,
#   2. fit the Gaussian-process model with a Dirichlet-process mixture on the
#      shock distribution (`gp_bnp` in gpsubspace_function.R),
#   3. fit the BART counterpart (`BART_bnp` in BART_function.R),
#   4. plot the one-step-ahead predictive densities and the in-sample fits.
#
# Notes on cost:
#   * `gp_bnp` does, per iteration, one O(n^3) Cholesky on the GP kernel and a
#     slice-sampler update over `G.max` mixture components.
#   * `BART_bnp` replaces the GP kernel with a BART ensemble; cost is then
#     O(nsave * n * trees) plus the mixture update.
#   * Defaults below (nburn = 500, ntot = 1000) are minimal. Production runs
#     in the paper use nburn = 2000, ntot = 4000.
# =============================================================================

rm(list = ls())

this_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)
setwd(this_dir)

suppressPackageStartupMessages({
  library(Rcpp)
  library(Matrix)
  library(stochvol)
  library(bvarsv)
  library(flexmix)
  library(dbarts)
  library(GIGrvg)
  library(ggplot2)
  library(reshape2)
})

# `gpsubspace_function.R` expects the GaussKernel_star symbol from
# BAKRGibbs.cpp to be available in the global environment.
Rcpp::sourceCpp("BAKRGibbs.cpp")

source("gpsubspace_function.R")    # gp_bnp
source("BART_function.R")          # BART_bnp, BART_mixSV
source("ucsv.R")                   # UCSV baseline
source("simulate_data.R")

dir.create("figures", showWarnings = FALSE)
dir.create("output",  showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
# REPLACE WITH YOUR DATA: build `y` (T x 1), `X` (T x p) on your own series,
# then standardise `y` to zero mean / unit variance and `X` column-wise.
# The paper targets CPIAUCSL or CPILFESL with FRED-QD predictors (see
# `infdata_script.R`).
set.seed(123)
sim <- simulate_subspace_data(T = 200, p_cov = 15, k_active = 4, seed = 42)
y   <- sim$y
X   <- sim$X

# Truth retained for overlay: true active predictors and the standardised
# true conditional mean of y given x_t.
true_active      <- sim$active
true_active_lbl  <- sim$var.names[true_active]

# Standardise (the samplers expect already-standardised data)
y_mu <- mean(y);  y_sd <- sd(y)
y <- (y - y_mu) / y_sd
X <- apply(X, 2, function(z) (z - mean(z)) / sd(z))
y_true_mean_std <- (sim$y_true_mean - y_mu) / y_sd

# Train / holdout split: hold out the last observation
n     <- nrow(y)
y.in  <- y[1:(n - 1), , drop = FALSE]
X.in  <- X[1:(n - 1), , drop = FALSE]
y.out <- y[n, ]
X.out <- X[n, , drop = FALSE]

# Add a constant column (as in `infdata_script.R`)
X.in  <- cbind(X.in, 1)
Xho   <- c(X.out, 1)

# The samplers expect `yho` as `c(<placeholder>, true_holdout_value)`
yho <- c(0, as.numeric(y.out))

# -----------------------------------------------------------------------------
# 2. Sampler configuration
# -----------------------------------------------------------------------------
nburn    <- 500
ntot     <- 1000
G.max    <- 10              # mixture components
PCA      <- TRUE            # SVD-based subspace step for X
mix.sv   <- FALSE           # set TRUE to enable the mixture-SV variant

# -----------------------------------------------------------------------------
# 3. Fit the GP-subspace model (`gp_bnp`)
# -----------------------------------------------------------------------------
fit_gp <- gp_bnp(y           = y.in,
                 X           = X.in,
                 yho         = yho,
                 Xho         = Xho,
                 nburn       = nburn,
                 ntot        = ntot,
                 thinning    = 2,
                 PDP         = FALSE,
                 sample.omega= TRUE,
                 tau2.start  = 0.5,
                 G.max       = G.max,
                 fcst        = TRUE,
                 a0          = 1, a1 = 1,
                 PCA         = PCA,
                 mix.sv      = mix.sv)

# -----------------------------------------------------------------------------
# 4. Fit the BART counterpart (`BART_bnp`)
# -----------------------------------------------------------------------------
fit_bart <- BART_bnp(y           = y.in,
                     X           = X.in,
                     yho         = yho,
                     Xho         = Xho,
                     nburn       = nburn,
                     ntot        = ntot,
                     thinning    = 2,
                     PDP         = FALSE,
                     sample.omega= FALSE,
                     tau2.start  = 1e+3,
                     G.max       = G.max,
                     fcst        = TRUE,
                     a0          = 1, a1 = 1,
                     PCA         = PCA,
                     mix.sv      = mix.sv)

# -----------------------------------------------------------------------------
# 5. Posterior summaries
# -----------------------------------------------------------------------------

# ----- 5a. Predictive densities for the held-out observation ----------------
pred_gp   <- as.numeric(fit_gp$predictions)
pred_bart <- as.numeric(fit_bart$predictions)

dens_df <- rbind(
  data.frame(model = "GP-subspace", value = pred_gp),
  data.frame(model = "BART",        value = pred_bart)
)

true_mean_out <- y_true_mean_std[n]
p_dens <- ggplot(dens_df, aes(x = value, colour = model, fill = model)) +
  geom_density(alpha = 0.2, linewidth = 0.8) +
  geom_vline(xintercept = as.numeric(y.out), colour = "black", linetype = 2) +
  geom_vline(xintercept = true_mean_out, colour = "red3", linetype = 3,
             linewidth = 0.7) +
  scale_colour_manual(values = c("GP-subspace" = "#1f77b4",
                                 "BART"        = "#d95f02")) +
  scale_fill_manual(  values = c("GP-subspace" = "#1f77b4",
                                 "BART"        = "#d95f02")) +
  labs(title = "Out-of-sample predictive density",
       subtitle = sprintf("Black dashed: realised y = %.3f.  Red dotted: true mean = %.3f.  True active: %s",
                          as.numeric(y.out), true_mean_out,
                          paste(true_active_lbl, collapse = ", ")),
       x = "Standardised y", y = "Density") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", legend.title = element_blank())

ggsave("figures/predictive_density.pdf", p_dens, width = 7, height = 4.5)
ggsave("figures/predictive_density.png", p_dens, width = 7, height = 4.5, dpi = 150)

# ----- 5b. In-sample fit for GP-subspace ------------------------------------
fit_q <- apply(fit_gp$fit, 2, quantile,
               probs = c(0.16, 0.5, 0.84), na.rm = TRUE)
fit_df <- data.frame(t    = seq_len(ncol(fit_q)),
                     obs  = as.numeric(y.in),
                     true = y_true_mean_std[seq_len(ncol(fit_q))],
                     lo   = fit_q[1, ],
                     med  = fit_q[2, ],
                     hi   = fit_q[3, ])

p_fit <- ggplot(fit_df, aes(x = t)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#1f77b4", alpha = 0.25) +
  geom_line(aes(y = med),  colour = "#1f77b4", linewidth = 0.7) +
  geom_line(aes(y = true), colour = "red3", linetype = "dashed",
            linewidth = 0.6) +
  geom_point(aes(y = obs), colour = "grey20", size = 0.8) +
  labs(title = "GP-subspace in-sample fit (16/50/84% bands)",
       subtitle = "Blue: posterior median.  Red dashed: true E[y|x] from the DGP.  Points: observed y.",
       x = "Time index", y = "Standardised y") +
  theme_minimal(base_size = 11)

ggsave("figures/insample_fit_gp.pdf", p_fit, width = 9, height = 4)
ggsave("figures/insample_fit_gp.png", p_fit, width = 9, height = 4, dpi = 150)

# ----- 5c. Log predictive likelihoods ---------------------------------------
lpl_summary <- function(z) {
  z <- as.numeric(z); z <- z[is.finite(z)]
  c(mean = mean(z), median = median(z),
    sd = sd(z), lo16 = quantile(z, 0.16), hi84 = quantile(z, 0.84))
}
lpl_tab <- rbind(`GP-subspace` = lpl_summary(fit_gp$lpl),
                 `BART`        = lpl_summary(fit_bart$lpl))
write.csv(lpl_tab, "output/log_predictive_likelihood.csv")

cat("\nLog predictive likelihood at the held-out point:\n")
print(round(lpl_tab, 3))

# -----------------------------------------------------------------------------
# 6. Persist the fits
# -----------------------------------------------------------------------------
saveRDS(list(gp = fit_gp, bart = fit_bart), file = "output/subspace_fit.rds")

cat("\nDone. Figures written to ./figures/.\n")
