# ============================================================================
# simulate_data.R
#
# Synthetic DGP for which the PC-subspace prior is a priori plausible:
# the M-dimensional series is driven by q_true latent persistent factors plus
# small idiosyncratic dynamics. This is the DGP used in Appendix B of
# Huber & Koop (JAE) for the subspace shrinkage paper.
#
# REPLACE WITH YOUR DATA in main_replication.R if you want to run the
# estimator on a real macro panel; this file only provides a stub.
# ============================================================================

sim.function <- function(N, q, T = 500,
                         sig.true = 0.01, sig.fac = 0.1,
                         omega = 1) {
  # Factor loadings
  Lambda <- matrix(rnorm(N * q, 0, 1), N, q)
  Lambda[1:q, 1:q] <- diag(q)              # identification: top block = I_q

  # Persistent factors
  F <- matrix(rnorm(q * T), T, q)
  for (t in 2:T)
    F[t, ] <- 0.99 * F[t - 1, ] + rnorm(q, 0, sig.fac)

  # Idiosyncratic VAR(1) component
  A <- matrix(rnorm(N * N, 0, 0.01), N, N)
  diag(A) <- 0.8

  # Cross-equation shock correlation
  Q.mat <- diag(N)
  Q.mat[lower.tri(Q.mat)] <- rnorm(N * (N - 1) / 2, 0, 0.01)

  Y <- matrix(rnorm(N, 0, sig.fac), T, N)
  for (t in 2:T)
    Y[t, ] <- omega * (Lambda %*% F[t - 1, ]) +
              (1 - omega) * A %*% Y[t - 1, ] +
              Q.mat %*% rnorm(N, 0, sig.true)
  list(Y = Y, Lambda = Lambda, F = F, A = A, Q = Q.mat,
       q_true = q, omega = omega)
}

# Convenience wrapper used in main_replication.R
make_synthetic_data <- function(M = 8, T = 200, q_true = 2, seed = 1) {
  set.seed(seed)
  sim <- sim.function(N = M, q = q_true, T = T,
                      sig.true = 0.01, sig.fac = 0.1, omega = 1)
  Y <- sim$Y[2:nrow(sim$Y), , drop = FALSE]
  colnames(Y) <- paste0("y", seq_len(ncol(Y)))
  attr(Y, "truth") <- list(Lambda = sim$Lambda, F = sim$F[2:nrow(sim$F), ],
                           A = sim$A, Q = sim$Q,
                           q_true = sim$q_true, omega = sim$omega)
  Y
}
