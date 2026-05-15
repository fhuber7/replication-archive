# =============================================================================
# simulate_data.R
# -----------------------------------------------------------------------------
# Build a synthetic multivariate time series with mildly nonlinear dynamics
# and mild stochastic volatility. The goal is just to give the GP-VAR
# something interesting to fit at a manageable size.
#
# REPLACE WITH YOUR DATA: in main_replication.R simply read your own
# data set into the matrix `Y` (rows = time, cols = variables, standardised).
# =============================================================================

#' Simulate a small VAR-like data set with nonlinear dynamics and mild SV
#'
#' @param T       Number of observations to retain.
#' @param M       Number of variables.
#' @param p       Lag length (for the linear part).
#' @param burnin  Burn-in observations (discarded).
#' @param seed    RNG seed.
#' @return        A list with `Y` (T x M matrix) and `var.names`.
simulate_gpvar_data <- function(T = 120, M = 3, p = 2, burnin = 200, seed = 1) {
  set.seed(seed)

  Tt <- T + burnin

  # Stable linear "backbone": small VAR(p) with diagonal-dominant lag-1 and
  # decaying higher lags. This guarantees stationarity so the recursion
  # below cannot blow up regardless of the nonlinear term.
  A_list <- vector("list", p)
  for (l in seq_len(p)) {
    A <- matrix(rnorm(M * M, 0, 0.05), M, M)
    diag(A) <- 0.5 / l
    A_list[[l]] <- A
  }

  # SV log-variances: persistent AR(1)s
  ht <- matrix(0, Tt, M)
  for (j in seq_len(M)) {
    for (t in 2:Tt) {
      ht[t, j] <- 0.95 * ht[t - 1, j] + rnorm(1, 0, 0.15)
    }
  }
  sig <- exp(ht / 2)

  # Lower-triangular contemporaneous structure
  L <- diag(M)
  L[lower.tri(L)] <- rnorm(M * (M - 1) / 2, 0, 0.3)

  Y <- matrix(0, Tt, M)
  for (t in (p + 1):Tt) {
    # Linear part
    lin <- numeric(M)
    for (l in seq_len(p)) lin <- lin + A_list[[l]] %*% Y[t - l, ]

    # Mild nonlinearity: tanh(.) on own lag-1, plus a small cross-term
    nl <- 0.2 * tanh(2 * Y[t - 1, ]) +
          0.1 * sin(Y[t - 1, c(2:M, 1)])

    # Structural shocks with SV
    u <- sig[t, ] * rnorm(M)
    e <- as.numeric(L %*% u)

    Y[t, ] <- as.numeric(lin) + nl + e
  }

  Y_full <- Y
  Y <- Y[(burnin + 1):Tt, , drop = FALSE]
  colnames(Y) <- paste0("y", seq_len(M))
  list(Y = Y, var.names = colnames(Y),
       # Truth retained for overlay in main_replication.R: linear backbone,
       # contemporaneous mixing, and the realised SV path (post burn-in).
       A_true  = A_list,
       L_true  = L,
       sig_true = sig[(burnin + 1):Tt, , drop = FALSE])
}
