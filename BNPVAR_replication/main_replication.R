# =============================================================================
# main_replication.R
# -----------------------------------------------------------------------------
# End-to-end driver for the BNP-VAR replication: simulate (or load) data,
# fit the Bayesian VAR with a Dirichlet-process mixture on the shocks
# (`bvar.mix` in Codes_and_data/aux.R), and plot
#   * regime-specific impulse responses to the last variable in the system,
#   * the posterior probability of being in each cluster over time,
#   * the posterior of the number of active clusters G.
#
# Notes on cost:
#   * Each Gibbs sweep does, for each of G.max clusters, a Wishart draw on a
#     (M x M) scale plus, per period, slice-sampling of the regime indicator.
#     Plus a VAR-coefficient update of O(M^2 K) and an IRF computation of
#     O(M^3 irf.hor) per cluster. With T = 200, M = 5, G.max = 10,
#     p = 2, nsave = nburn = 500, the script below takes a few minutes.
#   * Production runs in the paper use nsave = nburn = 5000.
# =============================================================================

rm(list = ls())

this_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)
setwd(this_dir)

suppressPackageStartupMessages({
  library(MASS)
  library(mvtnorm)
  library(mvnfast)
  library(stochvol)
  library(GIGrvg)
  library(bayesm)
  library(ggplot2)
  library(reshape2)
})

source("Codes_and_data/aux.R")
source("simulate_data.R")

dir.create("figures", showWarnings = FALSE)
dir.create("output",  showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
# REPLACE WITH YOUR DATA: assign a (T x M) numeric matrix to `Xraw`. The paper
# uses a quarterly US macro panel pulled from FRED-MD via
# Codes_and_data/Data_QD/prepare_dat.R, e.g.
#   var.S.true <- c("GDPC1","UNRATE","CPIAUCSL","FEDFUNDS",
#                   "GS10TB3Mx","BAA10YM","S.P.500")
set.seed(123)
sim   <- simulate_bnpvar_data(T = 200, M = 5, p = 2, seed = 42)
Xraw  <- sim$Y

# Truth from the DGP, retained for overlays below.
true_regime   <- sim$regime
true_G        <- sim$G_true
true_A_list   <- sim$A_list
true_L_calm   <- sim$L1
true_L_turb   <- sim$L2

# Standardise (the sampler does not standardise internally)
X.mean <- apply(Xraw, 2, mean)
X.sd   <- apply(Xraw, 2, sd)
Xraw   <- apply(Xraw, 2, function(x) (x - mean(x)) / sd(x))

M     <- ncol(Xraw)
labs  <- colnames(Xraw)

# -----------------------------------------------------------------------------
# 2. Sampler configuration
# -----------------------------------------------------------------------------
p       <- 2
nsave   <- 500
nburn   <- 500
G.max   <- 10                # max number of mixture components
SV.idi  <- FALSE             # SV on idiosyncratic shocks; FALSE for speed
fhorz   <- 1
irf.hor <- 24

# -----------------------------------------------------------------------------
# 3. Fit
# -----------------------------------------------------------------------------
fit <- bvar.mix(Xraw,
                nsave     = nsave,
                nburn     = nburn,
                G.max     = G.max,
                p         = p,
                SV.idi    = SV.idi,
                fhorz     = fhorz,
                irf.hor   = irf.hor,
                prior.cov = "AR",
                sample.w  = "conditional")

# -----------------------------------------------------------------------------
# 4. Posterior summaries
# -----------------------------------------------------------------------------

# ----- 4a. Number of active mixture components -------------------------------
# fit$G has dim (nsave, 2) with column 2 = G.plus (active components).
G_plus <- as.numeric(fit$G[, 2])
G_tab  <- prop.table(table(factor(G_plus, levels = seq_len(G.max))))
G_df   <- data.frame(G = as.integer(names(G_tab)),
                     prob = as.numeric(G_tab))

p_G <- ggplot(G_df, aes(x = factor(G), y = prob)) +
  geom_col(fill = "#1f77b4") +
  geom_vline(xintercept = which(levels(factor(G_df$G)) == as.character(true_G)),
             colour = "red3", linetype = "dashed", linewidth = 0.9) +
  labs(title = "Posterior of the number of active mixture components",
       subtitle = sprintf("Red dashed: true G = %d", true_G),
       x = "G+ (active clusters)", y = "Posterior probability") +
  theme_minimal(base_size = 11)

ggsave("figures/cluster_count.pdf", p_G, width = 6, height = 4)
ggsave("figures/cluster_count.png", p_G, width = 6, height = 4, dpi = 150)

# ----- 4b. Regime allocation probabilities over time ------------------------
# fit$S has dim (nsave, T). Compute, for each t, the posterior probability of
# being in each cluster (after the ex-post re-ordering inside bvar.mix).
S_prob <- matrix(0, ncol(fit$S), G.max)
for (j in seq_len(G.max)) {
  S_prob[, j] <- apply(fit$S == j, 2, mean)
}
keep <- which(colMeans(S_prob) > 0.01)         # only show non-trivial clusters
S_prob <- S_prob[, keep, drop = FALSE]
colnames(S_prob) <- paste0("cluster ", keep)
S_df <- melt(S_prob); colnames(S_df) <- c("time","cluster","prob")

# True regime indicator (1 = calm, 2 = turbulent). The estimator drops the
# first p observations, so align the truth on the same indexing.
true_regime_aln <- tail(true_regime, ncol(fit$S))
truth_path_df <- data.frame(time = seq_along(true_regime_aln),
                            value = as.numeric(true_regime_aln == 2))

p_S <- ggplot(S_df, aes(x = time, y = prob, colour = cluster)) +
  geom_line(linewidth = 0.7) +
  geom_step(data = truth_path_df, aes(x = time, y = value),
            colour = "red3", linetype = "dashed", linewidth = 0.5,
            inherit.aes = FALSE) +
  ylim(0, 1) +
  labs(title = "Posterior cluster-allocation probabilities",
       subtitle = "Red dashed step: true 'turbulent' regime indicator (DGP regime 2)",
       x = "Time index", y = "Probability") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", legend.title = element_blank())

ggsave("figures/cluster_path.pdf", p_S, width = 9, height = 4)
ggsave("figures/cluster_path.png", p_S, width = 9, height = 4, dpi = 150)

# ----- 4c. Impulse responses for cluster 1 ----------------------------------
# fit$IRF has dim (nsave, M, M, irf.hor, G.max). We plot the IRFs of all
# variables to a shock in the *last* variable (column M), under cluster 1.
shock_idx <- M
irf_q <- apply(fit$IRF[, , shock_idx, , 1], c(2, 3),
               quantile, probs = c(0.16, 0.5, 0.84), na.rm = TRUE)
# irf_q dim: 3 x M x irf.hor
irf_df <- data.frame()
for (m in seq_len(M)) {
  irf_df <- rbind(irf_df, data.frame(
    variable = labs[m],
    horizon  = seq_len(irf.hor),
    lo       = irf_q[1, m, ],
    med      = irf_q[2, m, ],
    hi       = irf_q[3, m, ]
  ))
}

# True IRF for the calm regime (regime 1 = L1, the largest-mass cluster).
# Standardisation step rescales rows of L by 1/sd(y_m); apply the same here.
sd_y <- X.sd
L_std <- diag(1 / sd_y) %*% true_L_calm
PHI   <- lapply(true_A_list, function(A) diag(1 / sd_y) %*% A %*% diag(sd_y))
true_irf_arr <- array(0, c(M, M, irf.hor))
true_irf_arr[, , 1] <- L_std
for (h in 2:irf.hor) for (l in 1:min(length(PHI), h - 1)) {
  true_irf_arr[, , h] <- true_irf_arr[, , h] + PHI[[l]] %*% true_irf_arr[, , h - l]
}
truth_irf_df <- data.frame()
for (m in seq_len(M)) {
  truth_irf_df <- rbind(truth_irf_df, data.frame(
    variable = labs[m],
    horizon  = seq_len(irf.hor),
    true     = true_irf_arr[m, shock_idx, ]
  ))
}

p_irf <- ggplot(irf_df, aes(x = horizon, y = med)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#d95f02", alpha = 0.25) +
  geom_line(colour = "#d95f02", linewidth = 0.8) +
  geom_line(data = truth_irf_df, aes(x = horizon, y = true),
            colour = "red3", linetype = "dashed", linewidth = 0.7,
            inherit.aes = FALSE) +
  facet_wrap(~ variable, scales = "free_y") +
  labs(title = paste0("Cluster-1 IRFs to a one-s.d. shock in ", labs[shock_idx]),
       subtitle = "Orange: posterior median + 68%.  Red dashed: true calm-regime IRF (DGP).",
       x = "Horizon", y = "Response") +
  theme_minimal(base_size = 11)

ggsave("figures/irf_cluster1.pdf", p_irf, width = 9, height = 5)
ggsave("figures/irf_cluster1.png", p_irf, width = 9, height = 5, dpi = 150)

# -----------------------------------------------------------------------------
# 5. Persist the fit
# -----------------------------------------------------------------------------
saveRDS(fit, file = "output/bnpvar_fit.rds")

cat("\nDone. Figures written to ./figures/.\n")
