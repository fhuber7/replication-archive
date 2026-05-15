# =============================================================================
# simulate_data.R
# -----------------------------------------------------------------------------
# Build a synthetic multivariate time series suitable for the mixBART /
# flexBART / fullBART samplers in this folder. The data-generating process is
# a stable linear VAR backbone with a mild contemporaneous shock factor
# structure and a smooth nonlinear cross-term, so that BART-type models have
# *some* nonlinearity to fit, but the series remain on the same scale as
# standardised macro data.
#
# REPLACE WITH YOUR DATA: in main_replication.R replace the call to
# simulate_mixbart_data(...) with code that loads your own (T x M) numeric
# matrix into `Ytrn` *before* standardisation.
# =============================================================================

#' Simulate a small multivariate time series with mild nonlinearity
#'
#' @param T       Number of retained observations.
#' @param M       Number of variables.
#' @param p       Lag length of the linear backbone.
#' @param burnin  Burn-in observations (discarded).
#' @param fhorz   Number of additional out-of-sample observations to retain
#'                as `Y_future` for truth-overlay on the forecast plot.
#' @param seed    RNG seed.
#' @return        A list with `Y` (T x M), `Y_future` (fhorz x M) and the
#'                DGP parameters (A_true, L_true, sig_true).
simulate_mixbart_data <- function(T = 200, M = 5, p = 2,
                                  burnin = 200, fhorz = 12, seed = 1) {
  set.seed(seed)
  Tt <- T + burnin + fhorz

  # Stable VAR(p) backbone: diagonal-dominant lag-1 and decaying higher lags
  A_list <- vector("list", p)
  for (l in seq_len(p)) {
    A <- matrix(rnorm(M * M, 0, 0.05), M, M)
    diag(A) <- 0.5 / l
    A_list[[l]] <- A
  }

  # Lower-triangular contemporaneous covariance factor
  L <- diag(M)
  L[lower.tri(L)] <- rnorm(M * (M - 1) / 2, 0, 0.3)

  # Mild stochastic volatility
  ht <- matrix(0, Tt, M)
  for (j in seq_len(M)) {
    for (t in 2:Tt) {
      ht[t, j] <- 0.95 * ht[t - 1, j] + rnorm(1, 0, 0.1)
    }
  }
  sig <- exp(ht / 2)

  Y <- matrix(0, Tt, M)
  for (t in (p + 1):Tt) {
    lin <- numeric(M)
    for (l in seq_len(p)) lin <- lin + A_list[[l]] %*% Y[t - l, ]

    # Mild nonlinearity: tanh on own lag-1 plus a sin cross-term
    nl <- 0.2 * tanh(2 * Y[t - 1, ]) +
          0.1 * sin(Y[t - 1, c(2:M, 1)])

    u <- sig[t, ] * rnorm(M)
    e <- as.numeric(L %*% u)
    Y[t, ] <- as.numeric(lin) + nl + e
  }

  Y_train  <- Y[(burnin + 1):(burnin + T), , drop = FALSE]
  Y_future <- Y[(burnin + T + 1):Tt, , drop = FALSE]
  colnames(Y_train) <- colnames(Y_future) <- paste0("y", seq_len(M))
  list(Y         = Y_train,
       Y_future  = Y_future,
       var.names = colnames(Y_train),
       A_true    = A_list,
       L_true    = L,
       sig_true  = sig)
}
