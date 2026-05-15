# =============================================================================
# Sparser TVP-VARs - replication driver
# Hauzenberger, Huber, Onorante (JAE, 2020):
#   "Combining Shrinkage and Sparsity in Conjugate Vector Autoregressive Models"
#
# Steps:
#   1. Load (synthetic) data; stub for real data is provided.
#   2. Run the conjugate Minnesota VAR with lag-wise SAVS sparsification + SV.
#   3. Compare sparsified vs. un-sparsified posterior medians.
#   4. Compute Cholesky IRFs with credible bands.
#   5. Compute recursive forecasts (fan chart) for the last variable.
#   6. Write all figures to ./figures/.
# =============================================================================

# --- working directory -------------------------------------------------------
this.dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
if (!nzchar(this.dir) || is.null(this.dir)) this.dir <- getwd()
setwd(this.dir)
dir.create("figures", showWarnings = FALSE)

suppressPackageStartupMessages({
  library(MASS)
  library(stochvol)
  library(mvtnorm)
  library(ggplot2)
  library(reshape2)
})

source("sparsevar_estim.R")
source("simulate_data.R")

# -----------------------------------------------------------------------------
# 1. Data
# -----------------------------------------------------------------------------
# REPLACE WITH YOUR DATA:
#   Yraw <- as.matrix(read.csv("your_data.csv"))   # T x M with column names
# To replicate the paper exactly use the FRED-MD data shipped in the source
# folder (see VAR/!submission/Code and Data/FREDdata/ and data.script.r).

sim <- simulate_sparse_var(M = 5, T = 250, p = 4, sparsity = 0.70, seed = 42)
Yraw <- sim$Yraw

cat("Data:", nrow(Yraw), "x", ncol(Yraw), "\n")
cat("True sparsity of slope matrix:", mean(sim$mask.true == 0), "\n")

# -----------------------------------------------------------------------------
# 2. Run the estimator
# -----------------------------------------------------------------------------
# Override defaults via command-line `Rscript main_replication.R 30 20` if you
# only want a smoke test.
args <- commandArgs(trailingOnly = TRUE)
nsave <- if (length(args) >= 1) as.integer(args[1]) else 3000
nburn <- if (length(args) >= 2) as.integer(args[2]) else 2000

cat("\nRunning sparsevar_estim with nsave =", nsave, ", nburn =", nburn, "\n")
t0 <- Sys.time()
fit <- sparsevar_estim(Yraw         = Yraw,
                       p            = 4,
                       nsave        = nsave,
                       nburn        = nburn,
                       cons         = TRUE,
                       prmean       = rep(0, ncol(Yraw)),
                       shrink.1     = 0.2,
                       kappa.base   = 2,
                       lambda.1st   = 1,
                       lambda.2nd   = 1,
                       lambda.decay = 0.5,
                       kappa.cov    = 4,
                       sigma.sparse = TRUE,
                       sv           = "pca.sv",
                       nhor         = 24,
                       fhorz        = 8,
                       verbose      = TRUE)
cat("Sampler done in", format(Sys.time() - t0), "\n")

# -----------------------------------------------------------------------------
# 3. Coefficient comparison: sparsified vs. un-sparsified
# -----------------------------------------------------------------------------
M <- fit$M; K <- fit$K; p <- fit$p
Kslope  <- M * p                      # number of slope rows (excl. intercept)
A.true  <- sim$A.true                 # Kslope x M (no intercept)
A.est   <- fit$A.median[1:Kslope, ]        # drop intercept row
A.full  <- fit$ALPHA.median[1:Kslope, ]

rmse <- function(x, y) sqrt(mean((x - y)^2))
cat("\n-- Coefficient comparison --\n")
cat(sprintf("  RMSE un-sparsified vs. truth: %.4f\n", rmse(A.full, A.true)))
cat(sprintf("  RMSE sparsified  vs. truth:   %.4f\n", rmse(A.est,  A.true)))
cat(sprintf("  True sparsity:               %.2f\n",  mean(A.true == 0)))
cat(sprintf("  Estimated sparsity:          %.2f\n",  mean(A.est  == 0)))

# Sparsity-pattern heatmap (true vs. estimated)
mat.true <- (A.true != 0) * 1
mat.est  <- (A.est  != 0) * 1
df.heat  <- rbind(
  data.frame(panel = "True",      melt(mat.true)),
  data.frame(panel = "Estimated", melt(mat.est))
)
names(df.heat) <- c("panel", "Lag", "Eq", "Active")
p.heat <- ggplot(df.heat, aes(x = factor(Eq), y = factor(Lag), fill = factor(Active))) +
  geom_tile(colour = "grey50") +
  facet_wrap(~ panel) +
  scale_fill_manual(values = c("0" = "white", "1" = "steelblue"),
                    name = "Coefficient",
                    labels = c("zero", "non-zero")) +
  scale_y_discrete(limits = rev) +
  labs(x = "Equation", y = "Row (lag x variable)",
       title = "Sparsity pattern: true vs. estimated") +
  theme_bw(base_size = 11)
ggsave("figures/sparsity_pattern.pdf", p.heat, width = 8, height = 6)

# -----------------------------------------------------------------------------
# 4. Cholesky impulse responses with credible bands (+ true IRFs from the DGP)
# -----------------------------------------------------------------------------
# Compute true IRFs: PHI_l = A.true rows ((l-1)*M+1):(l*M), impact = chol(Sigma.true)
H_irf <- dim(fit$irf.store)[4]
shk   <- 1
A_l   <- lapply(seq_len(p), function(l) sim$A.true[((l-1)*M+1):(l*M), ])
PHI_l <- lapply(A_l, t)   # response x shock convention used by impulsdtrf-style recursion
irf_true_arr <- array(0, c(M, M, H_irf))
irf_true_arr[ , , 1] <- t(chol(sim$Sigma.true))
for (h in 2:H_irf) for (l in 1:min(p, h - 1)) {
  irf_true_arr[ , , h] <- irf_true_arr[ , , h] + PHI_l[[l]] %*% irf_true_arr[ , , h - l]
}
true_irf_df <- data.frame(
  h        = rep(0:(H_irf - 1), times = M),
  response = rep(fit$var.names, each = H_irf),
  true     = as.numeric(t(irf_true_arr[ , shk, ])),
  stringsAsFactors = FALSE
)

mk_irf_df <- function(arr, label) {
  # arr: nsave x M x M x H ; shock = column 1 -> all responses
  arr <- arr[!is.na(arr[, 1, 1, 1]), , , , drop = FALSE]
  med <- apply(arr, c(2, 3, 4), median)
  lo  <- apply(arr, c(2, 3, 4), quantile, 0.16, na.rm = TRUE)
  hi  <- apply(arr, c(2, 3, 4), quantile, 0.84, na.rm = TRUE)
  H <- dim(arr)[4]
  shk <- 1                                  # shock to first variable
  data.frame(
    h        = rep(0:(H - 1), times = M),
    response = rep(fit$var.names, each = H),
    median   = as.numeric(t(med[, shk, ])),
    lo       = as.numeric(t(lo[,  shk, ])),
    hi       = as.numeric(t(hi[,  shk, ])),
    model    = label
  )
}
df.irf <- rbind(mk_irf_df(fit$irf.store,        "un-sparsified"),
                mk_irf_df(fit$irf.sparse.store, "sparsified"))
p.irf <- ggplot(df.irf, aes(x = h, y = median, colour = model, fill = model)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25, colour = NA) +
  geom_line() +
  geom_line(data = true_irf_df,
            aes(x = h, y = true), colour = "red3", linetype = "dashed",
            linewidth = 0.7, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  facet_wrap(~ response, scales = "free_y") +
  labs(x = "Horizon", y = "Response",
       title = sprintf("IRF: shock to %s (68%% band; red dashed = true)",
                       fit$var.names[1])) +
  theme_bw(base_size = 11)
ggsave("figures/irf.pdf", p.irf, width = 9, height = 6)

# -----------------------------------------------------------------------------
# 5. Predictive fan chart for the last variable
# -----------------------------------------------------------------------------
mk_pred_df <- function(arr, label) {
  arr <- arr[!is.na(arr[, 1, 1]), , , drop = FALSE]
  v <- M
  med <- apply(arr[, v, ], 2, median, na.rm = TRUE)
  q05 <- apply(arr[, v, ], 2, quantile, 0.05, na.rm = TRUE)
  q16 <- apply(arr[, v, ], 2, quantile, 0.16, na.rm = TRUE)
  q84 <- apply(arr[, v, ], 2, quantile, 0.84, na.rm = TRUE)
  q95 <- apply(arr[, v, ], 2, quantile, 0.95, na.rm = TRUE)
  H <- length(med)
  data.frame(h = 1:H, median = med, q05 = q05, q16 = q16,
             q84 = q84, q95 = q95, model = label)
}
df.pred <- rbind(mk_pred_df(fit$pred.store,        "un-sparsified"),
                 mk_pred_df(fit$pred.sparse.store, "sparsified"))
p.pred <- ggplot(df.pred, aes(x = h, y = median)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), fill = "steelblue", alpha = 0.2) +
  geom_ribbon(aes(ymin = q16, ymax = q84), fill = "steelblue", alpha = 0.4) +
  geom_line(colour = "navy") +
  facet_wrap(~ model) +
  labs(x = "Forecast horizon", y = sprintf("y_%d", M),
       title = sprintf("Predictive fan chart for %s (5/16/84/95)",
                       fit$var.names[M])) +
  theme_bw(base_size = 11)
ggsave("figures/forecast_fan.pdf", p.pred, width = 9, height = 5)

cat("\nFigures written to ./figures/\n")
cat("  - sparsity_pattern.pdf\n  - irf.pdf\n  - forecast_fan.pdf\n")
cat("\nDone.\n")
