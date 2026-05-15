# =============================================================================
# main_replication.R
# -----------------------------------------------------------------------------
# End-to-end driver for the GP-VAR replication: simulate (or load) data,
# fit the model, plot generalized impulse responses and the median SV path.
#
# Notes on cost:
#   * Each Gibbs sweep performs MM Cholesky / solve operations on T x T
#     kernel matrices, i.e. O(MM * T^3) per iteration. Empirically a single
#     iteration with M = 3, T = 120 takes a couple of seconds.
#   * The defaults below (nsave = 500, nburn = 500) are therefore a
#     *minimal* run useful for exploration. Production runs would need
#     several thousand draws.
# =============================================================================

rm(list = ls())

this_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)
setwd(this_dir)

source("simulate_data.R")
source("gpvar_estim.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
})

dir.create("figures", showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
# REPLACE WITH YOUR DATA: assign a (T x M) numeric matrix to `Y`. For the
# applied paper this is a quarterly US macro panel (GDP, UNR, CPI, FFR, GS1,
# GS10, S&P500) plus an uncertainty index.
sim <- simulate_gpvar_data(T = 120, M = 3, p = 2, seed = 42)
Y   <- sim$Y
# Truth retained for overlay below.
A_true   <- sim$A_true       # list of M x M lag matrices
L_true   <- sim$L_true       # contemporaneous Cholesky factor
sig_true <- sim$sig_true     # T x M realised SV path

# -----------------------------------------------------------------------------
# 2. Standardise (zero mean, unit variance)
# -----------------------------------------------------------------------------
Y <- apply(Y, 2, function(x) (x - mean(x)) / sd(x))

# -----------------------------------------------------------------------------
# 3. Fit the GP-VAR
# -----------------------------------------------------------------------------
# Defaults: nsave = nburn = 500. Increase for production, but expect runtime
# to grow linearly in (nsave + nburn) and cubically in T.
res <- gpvar(Y       = Y,
             p       = 2,
             nsave   = 500,
             nburn   = 500,
             sv      = TRUE,
             h.own   = 1,
             h.other = 1,
             shk.var = 1,
             shk.sc  = 1,
             nhor    = 17,
             girf.on = TRUE,
             verbose = TRUE)

# -----------------------------------------------------------------------------
# 4. Posterior summaries
# -----------------------------------------------------------------------------
MM <- res$MM
T  <- res$T
nhor <- dim(res$girf)[3]
var.names <- res$var.names

# --- GIRFs ------------------------------------------------------------------
girf_summary <- function(arr) {
  # arr is (nsave, MM, nhor)
  q <- apply(arr, c(2, 3), quantile, probs = c(0.16, 0.5, 0.84), na.rm = TRUE)
  # q has dim (3, MM, nhor) -> reshape to long
  out <- data.frame()
  for (m in seq_len(MM)) {
    out <- rbind(out, data.frame(
      variable = var.names[m],
      horizon  = seq_len(nhor) - 1,
      lo       = q[1, m, ],
      med      = q[2, m, ],
      hi       = q[3, m, ]
    ))
  }
  out
}

girf_pos <- girf_summary(res$girf[, , , 1]); girf_pos$shock <- "positive"
girf_neg <- girf_summary(res$girf[, , , 2]); girf_neg$shock <- "negative"
girf_df  <- rbind(girf_pos, girf_neg)

# --- TRUE GIRFs from the DGP via Monte Carlo --------------------------------
# Simulate two parallel paths (shocked / unshocked) through the nonlinear
# recursion in simulate_data.R and average over many draws. Uses the same
# linear backbone (A_true), contemporaneous mixing (L_true) and SV draw
# (sig_true). Standardisation by sd(y_m) in step 2 means we report GIRFs on
# the standardised scale so they are directly comparable to the model.
p_lag <- length(A_true)
M_g   <- nrow(A_true[[1]])
nMC   <- 1000
true_girf <- array(0, c(nMC, M_g, nhor, 2))   # MC x M x H x {pos, neg}
sd_y <- apply(sim$Y, 2, sd)
for (sh in 1:2) {
  shock_sign <- if (sh == 1) +1 else -1
  for (mc in seq_len(nMC)) {
    init <- t(replicate(p_lag, rnorm(M_g, 0, 0.5)))    # baseline state
    sv0  <- sig_true[1, ] / sd_y
    base_path <- shock_path <- matrix(0, nhor, M_g)
    state_b <- state_s <- init
    impact  <- shock_sign * sv0   # one-s.d. structural shock in var 1
    impact[-1] <- 0
    e_impact   <- as.numeric(L_true %*% impact)
    for (h in seq_len(nhor)) {
      eps <- rnorm(M_g)
      sv_h <- sig_true[(h %% nrow(sig_true)) + 1, ] / sd_y
      e   <- as.numeric(L_true %*% (sv_h * eps))
      lin_b <- numeric(M_g); lin_s <- numeric(M_g)
      for (l in seq_len(p_lag)) {
        lin_b <- lin_b + A_true[[l]] %*% state_b[l, ]
        lin_s <- lin_s + A_true[[l]] %*% state_s[l, ]
      }
      nl_b <- 0.2 * tanh(2 * state_b[1, ]) +
              0.1 * sin(state_b[1, c(2:M_g, 1)])
      nl_s <- 0.2 * tanh(2 * state_s[1, ]) +
              0.1 * sin(state_s[1, c(2:M_g, 1)])
      yb <- as.numeric(lin_b) + nl_b + e
      ys <- as.numeric(lin_s) + nl_s + e + (h == 1) * e_impact
      base_path[h, ]  <- yb
      shock_path[h, ] <- ys
      state_b <- rbind(yb, state_b)[1:p_lag, , drop = FALSE]
      state_s <- rbind(ys, state_s)[1:p_lag, , drop = FALSE]
    }
    true_girf[mc, , , sh] <- t(shock_path - base_path)
  }
}
true_girf_mean <- apply(true_girf, c(2, 3, 4), mean)
truth_df <- do.call(rbind, lapply(seq_len(M_g), function(m) {
  rbind(
    data.frame(variable = var.names[m], horizon = 0:(nhor - 1),
               true = true_girf_mean[m, , 1], shock = "positive"),
    data.frame(variable = var.names[m], horizon = 0:(nhor - 1),
               true = true_girf_mean[m, , 2], shock = "negative")
  )
}))

p_girf <- ggplot(girf_df, aes(x = horizon, y = med,
                              colour = shock, fill = shock)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_line(data = truth_df, aes(x = horizon, y = true, colour = shock),
            linetype = "dashed", linewidth = 0.7, inherit.aes = FALSE) +
  facet_wrap(~ variable, scales = "free_y") +
  scale_colour_manual(values = c(positive = "#d95f02", negative = "#1b9e77")) +
  scale_fill_manual(  values = c(positive = "#d95f02", negative = "#1b9e77")) +
  labs(title = paste0("GIRFs to a one-s.d. shock in ", var.names[1]),
       subtitle = "Solid + band: posterior median + 68%.  Dashed: true GIRF (Monte Carlo on the DGP).",
       x = "Horizon", y = "Response") +
  theme_minimal(base_size = 11)

ggsave("figures/girfs.pdf", p_girf, width = 9, height = 5)
ggsave("figures/girfs.png", p_girf, width = 9, height = 5, dpi = 150)

# --- Posterior median SV path ----------------------------------------------
ht_med <- apply(res$ht, c(2, 3), median, na.rm = TRUE)  # T x MM
sv_df <- data.frame()
for (m in seq_len(MM)) {
  sv_df <- rbind(sv_df, data.frame(
    time     = seq_len(T),
    variable = var.names[m],
    sv       = exp(ht_med[, m] / 2)
  ))
}

# True SV path. The DGP simulates `sig = exp(h/2)` directly; standardisation
# of Y in step 2 rescales the apparent SV path by 1/sd(y_m), so we rescale
# the truth to match the standardised series the model sees. The estimator
# drops the first p observations as initial conditions, so trim sig_true.
sd_y <- apply(sim$Y, 2, sd)
sig_true_aln <- sig_true[(nrow(sig_true) - T + 1):nrow(sig_true), , drop = FALSE]
sv_df$true   <- as.numeric(t(t(sig_true_aln) / sd_y))

p_sv <- ggplot(sv_df, aes(x = time, y = sv)) +
  geom_line(colour = "#1f77b4", linewidth = 0.8) +
  geom_line(aes(y = true), colour = "red3", linetype = "dashed",
            linewidth = 0.7) +
  facet_wrap(~ variable, scales = "free_y") +
  labs(title = "Posterior median stochastic-volatility path  (exp(h_t / 2))",
       subtitle = "Blue: posterior median.  Red dashed: realised SV path from the DGP (rescaled).",
       x = "Time", y = expression(sigma[t])) +
  theme_minimal(base_size = 11)

ggsave("figures/sv_path.pdf", p_sv, width = 9, height = 4)
ggsave("figures/sv_path.png", p_sv, width = 9, height = 4, dpi = 150)

# -----------------------------------------------------------------------------
# 5. Persist the fit
# -----------------------------------------------------------------------------
saveRDS(res, file = "figures/gpvar_fit.rds")

cat("\nDone. Figures written to ./figures/.\n")
