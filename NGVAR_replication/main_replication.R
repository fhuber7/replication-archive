# main_replication.R
# ------------------------------------------------------------------------------
# End-to-end replication driver for the NG-VAR model with stochastic volatility
# (Huber and Feldkircher, 2019, JBES). Loads the synthetic data set produced by
# simulate_data.R, runs the MCMC sampler defined in ngvar_estim.R, and produces
# posterior IRFs (Cholesky) and predictive densities for the last variable.
# ------------------------------------------------------------------------------

# Set working directory to the location of this script (works from Rscript)
# -----------------------------------------------------------------------------
if (interactive()) {
  this_dir <- getwd()
} else {
  args0 <- commandArgs(trailingOnly = FALSE)
  fileArg <- sub("^--file=", "", args0[grep("^--file=", args0)])
  this_dir <- if (length(fileArg) > 0) normalizePath(dirname(fileArg)) else getwd()
}
setwd(this_dir)

dir.create("output",  showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ---- MCMC settings ----------------------------------------------------------
# Allow command-line overrides:  Rscript main_replication.R nsave=20 nburn=10
nsave <- 2000
nburn <- 1000
seed  <- 1234
cli_args <- commandArgs(trailingOnly = TRUE)
for (a in cli_args) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  if (length(kv) == 2 && kv[1] %in% c("nsave", "nburn", "seed")) {
    assign(kv[1], as.integer(kv[2]))
  }
}
set.seed(seed)
cat("main_replication.R: nsave =", nsave, ", nburn =", nburn,
    ", seed =", seed, "\n")

# ---- load estimator and data ------------------------------------------------
source("ngvar_estim.R")

# === REPLACE WITH YOUR DATA ==================================================
# If you have your own (T x M) data matrix, build/load it here directly and
# skip simulate_data.R. Example:
#   Y      <- as.matrix(read.csv("my_data.csv"))
#   labels <- colnames(Y)
# =============================================================================
if (!file.exists("data_synth.rda")) source("simulate_data.R")
load("data_synth.rda")
cat("Loaded data:", nrow(Y), "obs x", ncol(Y), "variables (",
    paste(labels, collapse = ", "), ")\n")

# ---- estimation -------------------------------------------------------------
fit <- ngvar_estim(Yraw        = Y,
                   p           = 2,
                   nsave       = nsave,
                   nburn       = nburn,
                   cons        = 1,
                   sv          = TRUE,
                   shrink.type = "lagwise",
                   verbose     = TRUE)

# ---- posterior summary of VAR coefficients ----------------------------------
A_med <- apply(fit$A, c(2, 3), median)
A_q16 <- apply(fit$A, c(2, 3), quantile, 0.16)
A_q84 <- apply(fit$A, c(2, 3), quantile, 0.84)
colnames(A_med) <- colnames(A_q16) <- colnames(A_q84) <- labels

# ---- impulse responses (Cholesky) -------------------------------------------
nhor <- 20
IRF  <- ngvar_irf(fit, nhor = nhor)
keep <- which(apply(IRF, 1, function(x) all(is.finite(x))))
cat("IRF: kept", length(keep), "of", dim(IRF)[1], "draws (stable).\n")
IRF  <- IRF[keep, , , , drop = FALSE]

irf_med <- apply(IRF, c(2, 3, 4), median)
irf_q16 <- apply(IRF, c(2, 3, 4), quantile, 0.16)
irf_q84 <- apply(IRF, c(2, 3, 4), quantile, 0.84)

# ---- TRUE IRFs from the DGP -------------------------------------------------
# y_t = sum_l A_l y_{t-l} + Uinv eps_t,  eps_t ~ N(0, diag(exp(h_t))).
# Reference shock covariance = Uinv exp(mean h_t) Uinv'. Cholesky IRFs of the
# reduced form follow the standard recursion in PHI = A_l (l = 1..p).
# irf_true[ , , h] = response (rows) at horizon h to a one-s.d. structural
# shock in each variable (cols), under the average SV state.
true_irf <- function(A_list, Sigma, nhor) {
  M <- nrow(A_list[[1]]); p <- length(A_list)
  out <- array(0, c(M, M, nhor))
  out[ , , 1] <- t(chol(Sigma))
  for (h in 2:nhor) for (l in 1:min(p, h - 1)) {
    out[ , , h] <- out[ , , h] + A_list[[l]] %*% out[ , , h - l]
  }
  out
}
if (exists("truth", inherits = TRUE)) {
  Uinv      <- solve(truth$U_true)
  Sigma_avg <- Uinv %*% diag(exp(colMeans(truth$h_true))) %*% t(Uinv)
  irf_true  <- true_irf(truth$A_true, Sigma_avg, nhor)
} else irf_true <- NULL

# ---- predictive density for the last variable -------------------------------
fhorz <- 8
M <- fit$M; p <- fit$p; cons <- fit$cons; K <- fit$K
nsave_eff <- dim(fit$A)[1]
Xraw_last <- if (cons == 1) {
  c(as.vector(t(Y[nrow(Y):(nrow(Y) - p + 1), ])), 1)
} else {
  c(as.vector(t(Y[nrow(Y):(nrow(Y) - p + 1), ])))
}
pred_store <- array(NA_real_, c(nsave_eff, M, fhorz))
for (s in 1:nsave_eff) {
  A_s <- fit$A[s, , ]
  U_s <- fit$U[s, , ]
  H_T <- fit$H[s, fit$T, ]
  pv  <- fit$sv_pars[s, , ]
  cmp <- get_companion(A_s, c(M, cons, p))
  Mm  <- cmp$MM
  if (max(abs(Re(eigen(Mm, only.values = TRUE)$values))) > 1.05) next
  st  <- Xraw_last
  Htt <- H_T
  for (h in 1:fhorz) {
    st  <- Mm %*% st
    Sig <- U_s %*% diag(exp(Htt), M) %*% t(U_s)
    chS <- try(t(chol(Sig)), silent = TRUE)
    if (inherits(chS, "try-error")) break
    yf  <- st[1:M] + chS %*% rnorm(M)
    pred_store[s, , h] <- yf
    # propagate log-volatility one step; pv[3,] is the std dev sigma
    Htt <- pv[1, ] + pv[2, ] * (Htt - pv[1, ]) +
           rnorm(M, 0, pmax(pv[3, ], 1e-8))
  }
}
pred_q <- apply(pred_store, c(2, 3), quantile,
                c(0.05, 0.16, 0.5, 0.84, 0.95), na.rm = TRUE)

# ---- save output ------------------------------------------------------------
save(fit, A_med, A_q16, A_q84,
     IRF, irf_med, irf_q16, irf_q84,
     pred_store, pred_q,
     file = file.path("output", "ngvar_results.rda"))
cat("Saved posterior output to output/ngvar_results.rda\n")

# ---- IRF plot (base R) ------------------------------------------------------
plot_irf <- function(med, lo, hi, var_names, file, irf_true = NULL) {
  M <- dim(med)[1]; H <- dim(med)[3]
  pdf(file, width = 2.2 * M, height = 2.0 * M)
  op <- par(mfrow = c(M, M), mar = c(2.4, 2.4, 1.6, 0.6),
            mgp = c(1.3, 0.4, 0))
  for (i in 1:M) for (j in 1:M) {
    ylvec <- c(lo[i, j, ], hi[i, j, ], 0)
    if (!is.null(irf_true)) ylvec <- c(ylvec, irf_true[i, j, ])
    yl <- range(ylvec)
    plot(0:(H - 1), med[i, j, ], type = "l", lwd = 1.8,
         ylim = yl, xlab = "horizon", ylab = "",
         main = paste0(var_names[i], " <- ", var_names[j]),
         cex.main = 0.9)
    polygon(c(0:(H - 1), (H - 1):0),
            c(lo[i, j, ], rev(hi[i, j, ])),
            col = rgb(0.2, 0.4, 0.8, 0.25), border = NA)
    lines(0:(H - 1), med[i, j, ], lwd = 1.8)
    if (!is.null(irf_true))
      lines(0:(H - 1), irf_true[i, j, 1:H], lwd = 1.6, lty = 2, col = "red3")
    abline(h = 0, lty = 3, col = "grey40")
    if (i == 1 && j == 1 && !is.null(irf_true))
      legend("topright", c("posterior median", "68% band", "true"),
             col = c("black", rgb(0.2, 0.4, 0.8), "red3"),
             lty = c(1, NA, 2), lwd = c(1.8, NA, 1.6),
             pch = c(NA, 15, NA), bty = "n", cex = 0.75)
  }
  par(op); dev.off()
}
plot_irf(irf_med, irf_q16, irf_q84, labels,
         file.path("figures", "irf.pdf"),
         irf_true = irf_true)
cat("Wrote figures/irf.pdf\n")

# ---- predictive density fan chart for the last variable ---------------------
plot_fan <- function(pred_q_mat, Y_last, var_name, file) {
  # pred_q_mat is a (5 x H) matrix of quantiles (0.05, 0.16, 0.5, 0.84, 0.95).
  H <- ncol(pred_q_mat)
  pdf(file, width = 6, height = 4)
  hist_len <- min(20, length(Y_last))
  hist_seg <- Y_last[(length(Y_last) - hist_len + 1):length(Y_last)]
  xs_hist  <- (-(hist_len - 1)):0
  xs_fcst  <- 1:H
  yl <- range(c(hist_seg, pred_q_mat[1, ], pred_q_mat[5, ]), na.rm = TRUE)
  plot(c(xs_hist, xs_fcst), c(hist_seg, rep(NA, H)), type = "n",
       ylim = yl, xlab = "horizon", ylab = var_name,
       main = paste("Predictive density:", var_name))
  polygon(c(xs_fcst, rev(xs_fcst)),
          c(pred_q_mat[1, ], rev(pred_q_mat[5, ])),
          col = rgb(0.2, 0.4, 0.8, 0.15), border = NA)
  polygon(c(xs_fcst, rev(xs_fcst)),
          c(pred_q_mat[2, ], rev(pred_q_mat[4, ])),
          col = rgb(0.2, 0.4, 0.8, 0.30), border = NA)
  lines(xs_hist, hist_seg, lwd = 1.6)
  lines(xs_fcst, pred_q_mat[3, ], lwd = 1.6, col = "darkblue")
  abline(v = 0, lty = 2, col = "grey50")
  legend("topleft", lwd = c(1.6, 1.6), col = c("black", "darkblue"),
         legend = c("history", "posterior median"), bty = "n", cex = 0.85)
  dev.off()
}
last_idx  <- ncol(Y)
slice_mat <- pred_q[, last_idx, ]   # 5 x H
plot_fan(slice_mat, Y[, last_idx], labels[last_idx],
         file.path("figures", "predictive_density.pdf"))
cat("Wrote figures/predictive_density.pdf\n")
cat("\nDone.\n")
