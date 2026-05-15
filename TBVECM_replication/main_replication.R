# ============================================================================
# main_replication.R
#
# End-to-end driver for the TBVECM replication package: simulates a
# three-regime threshold VECM, estimates the model, and produces
# regime-conditional impulse responses and regime-probability plots.
# ============================================================================

# Working directory ----------------------------------------------------------
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) &&
    !is.null(rstudioapi::getSourceEditorContext()$path)) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
} else if (sys.nframe() > 0L) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA)
  if (!is.na(this_file)) setwd(dirname(this_file))
}

dir.create("figures", showWarnings = FALSE)

source("tbvecm_estim.R")
source("simulate_data.R")

# ---------------------------------------------------------------------------
# 1. Data
# ---------------------------------------------------------------------------
# REPLACE WITH YOUR DATA -----------------------------------------------------
# To use your own time series, replace the block below with code that
# loads a (T x M) matrix `Yraw` of levels of the endogenous variables.
# Example:
#   raw   <- read.csv("my_macro_data.csv")
#   Yraw  <- as.matrix(raw[, c("CPI", "IP", "INT")])
# ---------------------------------------------------------------------------

sim   <- simulate_tbvecm_data(T = 300, M = 3, seed = 42)
Yraw  <- sim$Y

# ---------------------------------------------------------------------------
# 2. MCMC settings
# ---------------------------------------------------------------------------
# Allow CLI override (useful for the smoke test). Both forms work:
#   Rscript main_replication.R 30 20
#   Rscript main_replication.R nsave=30 nburn=20
.parse_args <- function(defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  out  <- defaults
  pos  <- 1L
  for (a in args) {
    if (grepl("=", a, fixed = TRUE)) {
      kv  <- strsplit(a, "=", fixed = TRUE)[[1]]
      key <- kv[1]; val <- suppressWarnings(as.integer(kv[2]))
      if (key %in% names(out) && !is.na(val)) out[[key]] <- val
    } else {
      val <- suppressWarnings(as.integer(a))
      if (!is.na(val) && pos <= length(out)) {
        out[[names(out)[pos]]] <- val
        pos <- pos + 1L
      }
    }
  }
  out
}
.opts <- .parse_args(list(nsave = 3000L, nburn = 2000L, nhor = 24L))
nsave <- .opts$nsave
nburn <- .opts$nburn
nhor  <- .opts$nhor

message(sprintf("Running TBVECM with nsave=%d, nburn=%d, nhor=%d",
                nsave, nburn, nhor))

t0  <- Sys.time()
fit <- tbvecm_estim(
  Yraw           = Yraw,
  r              = 1,
  p              = 1,
  nsave          = nsave,
  nburn          = nburn,
  nhor           = nhor,
  beta_prior_var = 1.0,
  coint.shrink   = FALSE,
  sample.theta   = TRUE,
  cons           = TRUE,
  verbose_every  = max(100, (nsave + nburn) %/% 10)
)
t1  <- Sys.time()
message("Elapsed: ", format(round(difftime(t1, t0, units = "secs"), 1)))

# ---------------------------------------------------------------------------
# 3. Posterior summaries
# ---------------------------------------------------------------------------
beta_post  <- apply(fit$beta,  c(2, 3), median, na.rm = TRUE)
gamma_post <- apply(fit$gamma, 2,      median, na.rm = TRUE)

cat("\nPosterior median of beta (linearly normalised):\n"); print(beta_post)
cat("\nTrue beta used in the DGP:\n");                       print(sim$beta_true)
cat("\nPosterior median of thresholds (gamma1, gamma2):\n"); print(gamma_post)
cat("True simulation thresholds:\n");                        print(sim$thresholds)

# Regime probabilities at each t -------------------------------------------
M       <- fit$dims$M
T_used  <- fit$dims$T
nsv     <- nsave
reg_arr <- array(0, c(nsv, T_used, 3))
for (rg in 0:2) {
  reg_arr[, , rg + 1] <- (fit$st == rg) * 1
}
regime_prob <- apply(reg_arr, c(2, 3), mean, na.rm = TRUE)
colnames(regime_prob) <- paste0("Regime_", 1:3)

# True regime indicator (DGP). simulate_data.R retains the latent regime
# sequence in sim$states (full sample including p initial lags); the
# estimator uses the same final-T-of-the-window indexing.
true_states <- tail(sim$states, T_used)
cat("\nTrue regime counts vs. estimated (Pr > 0.5):\n")
print(rbind(true  = table(factor(true_states, levels = 1:3)),
            postPr = colSums(regime_prob > 0.5)))

# Save posterior summaries --------------------------------------------------
saveRDS(
  list(beta = fit$beta, gamma = fit$gamma, regime_prob = regime_prob,
       dims = fit$dims),
  file = "posterior_summaries.rds"
)

# ---------------------------------------------------------------------------
# 4. Plot: regime probabilities
# ---------------------------------------------------------------------------
have_gg <- requireNamespace("ggplot2",  quietly = TRUE) &&
           requireNamespace("reshape2", quietly = TRUE)

if (have_gg) {
  library(ggplot2)
  library(reshape2)
  df_reg  <- data.frame(t = seq_len(T_used), regime_prob)
  df_long <- reshape2::melt(df_reg, id.vars = "t",
                            variable.name = "regime",
                            value.name = "prob")
  # Build a "true" indicator panel: 1 in the column matching the true state.
  true_long <- do.call(rbind, lapply(1:3, function(rg) {
    data.frame(t = seq_len(T_used),
               regime = factor(paste0("Regime_", rg),
                               levels = paste0("Regime_", 1:3)),
               truth = as.integer(true_states == rg))
  }))
  df_long <- merge(df_long, true_long, by = c("t", "regime"))
  p_reg <- ggplot(df_long, aes(x = t)) +
    geom_area(aes(y = prob, fill = regime), alpha = 0.7) +
    geom_step(aes(y = truth), colour = "red3", linewidth = 0.4) +
    facet_wrap(~ regime, ncol = 1) +
    labs(title  = "Posterior regime probabilities (red step = true regime)",
         x = "Time", y = "Probability") +
    theme_minimal() +
    theme(legend.position = "none")
  ggsave("figures/regime_probabilities.pdf", p_reg, width = 7, height = 6)
} else {
  pdf("figures/regime_probabilities.pdf", width = 7, height = 6)
  par(mfrow = c(3, 1), mar = c(3, 4, 2, 1))
  for (rg in 1:3) {
    plot(regime_prob[, rg], type = "h", ylim = c(0, 1),
         ylab = paste("Pr(Regime", rg, ")"), xlab = "Time",
         main = paste("Posterior probability, regime", rg),
         col = "steelblue")
    lines(seq_len(T_used), as.integer(true_states == rg),
          type = "s", col = "red3", lwd = 1.5)
  }
  dev.off()
}

# ---------------------------------------------------------------------------
# 5. Plot: regime-conditional impulse responses
# ---------------------------------------------------------------------------
# Shock = unit standard deviation shock to variable 1 (Cholesky ordering),
# stored as fit$IRF[draw, response_var, shock_var, horizon, regime].
shock_var <- 1
irf_med <- apply(fit$IRF[, , shock_var, , , drop = FALSE],
                 c(2, 4, 5), median, na.rm = TRUE)
irf_lo  <- apply(fit$IRF[, , shock_var, , , drop = FALSE],
                 c(2, 4, 5), quantile, probs = 0.16, na.rm = TRUE)
irf_hi  <- apply(fit$IRF[, , shock_var, , , drop = FALSE],
                 c(2, 4, 5), quantile, probs = 0.84, na.rm = TRUE)

varnames <- fit$dims$varnames
dimnames(irf_med) <- list(varnames, NULL, paste0("R", 1:3))

# ---- TRUE regime-specific IRFs from the DGP --------------------------------
# Level VAR coefficients: A_1 = I + alpha beta' + Phi, A_2 = -Phi.
# Impact = chol(Sigma) for each regime. Recursion: IRF[,,h] = A_1 IRF[,,h-1]
# + A_2 IRF[,,h-2].
true_regime_irf <- function(alpha, Phi, Sigma, beta, M, H) {
  A1 <- diag(M) + as.numeric(t(alpha)) %*% t(beta) + Phi
  A2 <- -Phi
  out <- array(0, c(M, M, H))
  out[, , 1] <- t(chol(Sigma))
  if (H >= 2) out[, , 2] <- A1 %*% out[, , 1]
  if (H >= 3) for (h in 3:H) {
    out[, , h] <- A1 %*% out[, , h - 1] + A2 %*% out[, , h - 2]
  }
  out
}
irf_true_arr <- array(NA, c(M, M, nhor, 3))
for (rg in 1:3) {
  irf_true_arr[, , , rg] <- true_regime_irf(
    alpha = sim$alpha_list[[rg]],
    Phi   = sim$Phi_list[[rg]],
    Sigma = sim$Sigma_list[[rg]],
    beta  = sim$beta_true,
    M     = M, H = nhor
  )
}
# Slice to: response of variable v to a shock in shock_var, per regime.
irf_true_slice <- irf_true_arr[, shock_var, , , drop = FALSE]
dim(irf_true_slice) <- c(M, nhor, 3)

pdf("figures/irf_regimeconditional.pdf", width = 10, height = 7)
par(mfrow = c(M, 3), mar = c(3.5, 4, 2.5, 1))
for (v in 1:M) {
  for (rg in 1:3) {
    yrange <- range(c(irf_lo[v, , rg], irf_hi[v, , rg], irf_true_slice[v, , rg], 0),
                    na.rm = TRUE)
    plot(irf_med[v, , rg], type = "l", lwd = 2, col = "darkblue",
         ylim = yrange, xlab = "Horizon",
         ylab = paste("Response of", varnames[v]),
         main = paste0("Regime ", rg, ": shock to ", varnames[shock_var]))
    polygon(c(seq_len(nhor), rev(seq_len(nhor))),
            c(irf_lo[v, , rg], rev(irf_hi[v, , rg])),
            col = adjustcolor("steelblue", alpha.f = 0.3), border = NA)
    lines(irf_med[v, , rg], lwd = 2, col = "darkblue")
    lines(seq_len(nhor), irf_true_slice[v, , rg],
          lwd = 2, col = "red3", lty = 2)
    abline(h = 0, col = "grey60", lty = 3)
    grid()
    if (v == 1 && rg == 1)
      legend("topright", c("posterior median", "68% band", "true (DGP)"),
             col = c("darkblue", adjustcolor("steelblue", alpha.f = 0.3), "red3"),
             lty = c(1, NA, 2), lwd = c(2, NA, 2),
             pch = c(NA, 15, NA), bty = "n", cex = 0.7)
  }
}
dev.off()

# ---------------------------------------------------------------------------
# 6. Plot: posterior of beta
# ---------------------------------------------------------------------------
pdf("figures/beta_posterior.pdf", width = 7, height = 8)
par(mfrow = c(M, 1), mar = c(3.5, 4, 2.5, 1))
for (j in 1:M) {
  dr <- fit$beta[, j, 1]
  dr <- dr[is.finite(dr)]
  if (length(unique(dr)) <= 1) {
    plot.new(); title(main = paste0("beta[", j, "] (fixed by normalisation = ",
                                    signif(unique(dr), 3), ")"))
  } else {
    plot(density(dr), main = paste0("Posterior density of beta[", j, "]"),
         xlab = "", lwd = 2, col = "darkblue")
    abline(v = median(dr), lty = 2)
    abline(v = sim$beta_true[j, 1], col = "red", lwd = 2)
    grid()
  }
}
legend("topright", legend = c("posterior median", "true value"),
       col = c("black", "red"), lty = c(2, 1), bty = "n")
dev.off()

message("\nWrote: figures/regime_probabilities.pdf")
message("       figures/irf_regimeconditional.pdf")
message("       figures/beta_posterior.pdf")
message("       posterior_summaries.rds")
message("Done.")
