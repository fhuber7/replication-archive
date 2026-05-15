# =============================================================================
# simulate_data.R
# -----------------------------------------------------------------------------
# Sparse TVP-VAR data-generating process for the SAVS replication. The DGP
# has the structure that SAVS is *designed* to recover: most of the VAR
# coefficients are zero, and the surviving entries drift slowly over time as
# random walks. Innovations carry mild stochastic volatility.
#
# REPLACE WITH YOUR DATA: in `main_replication.R` replace the call to
# `simulate_savs_tvp_var()` with code that loads a (T x M) matrix `Y` and
# (optionally) a true coefficient path. The real-data driver
# `estim.VAR.R` already shows how to use this estimator on a quarterly US
# macro panel (`datasetMNg.RData`).
# =============================================================================

#' Simulate a sparse TVP-VAR(p) with mild SV
#'
#' @param T          Number of observations to retain.
#' @param M          Number of variables.
#' @param p          Lag length.
#' @param sparsity   Fraction of *off-diagonal* lag coefficients pinned at zero.
#' @param drift_sd   RW innovation s.d. on the surviving coefficient paths.
#' @param sv_sigma   Innovation s.d. of the per-equation log-volatility AR(1).
#' @param seed       RNG seed.
#' @param fhorz      Number of additional out-of-sample observations retained
#'                   as `Y_future` for truth-overlay on the forecast plot.
#' @return           A list with
#'                     Y         (T x M) in-sample data
#'                     Y_future  (fhorz x M) realised future y
#'                     A_path    (T x K x M) true time-varying lag coefficients
#'                     mask      (K x M) sparsity pattern (1 = active)
#'                     L         (M x M) contemporaneous Cholesky factor
#'                     sig       ((T + fhorz) x M) realised SV path
simulate_savs_tvp_var <- function(T = 150, M = 3, p = 2,
                                  sparsity = 0.6,
                                  drift_sd = 0.01,
                                  sv_sigma = 0.10,
                                  seed = 1, fhorz = 8) {
  set.seed(seed)
  K  <- M * p
  Tt <- T + fhorz

  # ---- Initial coefficient matrix --------------------------------------------
  A0 <- matrix(rnorm(K * M, 0, 0.05), K, M)
  diag(A0[1:M, ]) <- runif(M, 0.45, 0.65)   # own-lag-1 persistence

  # Sparsity mask. Always keep the own-lag-1 diagonal active.
  mask <- matrix(rbinom(K * M, 1, 1 - sparsity), K, M)
  diag(mask[1:M, ]) <- 1
  A0 <- A0 * mask

  # ---- Time-varying paths: RW on active entries, zero elsewhere -------------
  A_path <- array(0, c(Tt, K, M))
  A_path[1, , ] <- A0
  for (t in 2:Tt) {
    drift <- matrix(rnorm(K * M, 0, drift_sd), K, M) * mask
    A_path[t, , ] <- A_path[t - 1, , ] + drift
  }

  # ---- Contemporaneous covariance -------------------------------------------
  L <- diag(M)
  L[lower.tri(L)] <- rnorm(M * (M - 1) / 2, 0, 0.2)

  # ---- Stochastic volatility -------------------------------------------------
  ht <- matrix(0, Tt, M)
  for (j in seq_len(M)) {
    for (t in 2:Tt) {
      ht[t, j] <- 0.95 * ht[t - 1, j] + rnorm(1, 0, sv_sigma)
    }
  }
  sig <- exp(ht / 2)

  # ---- Simulate y -----------------------------------------------------------
  Y <- matrix(0, Tt + p, M)
  Y[1:p, ] <- matrix(rnorm(p * M, 0, 0.5), p, M)
  for (t in seq_len(Tt)) {
    lags  <- as.numeric(t(Y[(p + t - 1):t, ]))            # K-vector
    mu_t  <- as.numeric(crossprod(A_path[t, , ], lags))   # M-vector
    eps_t <- as.numeric(L %*% (sig[t, ] * rnorm(M)))
    Y[p + t, ] <- mu_t + eps_t
  }
  Y <- Y[(p + 1):(Tt + p), , drop = FALSE]
  colnames(Y) <- paste0("y", seq_len(M))

  Y_in     <- Y[seq_len(T), , drop = FALSE]
  Y_future <- Y[(T + 1):Tt, , drop = FALSE]

  list(Y         = Y_in,
       Y_future  = Y_future,
       A_path    = A_path[seq_len(T), , , drop = FALSE],
       mask      = mask,
       L         = L,
       sig       = sig)
}
