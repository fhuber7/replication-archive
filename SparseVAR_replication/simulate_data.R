# =============================================================================
# Synthetic sparse VAR(p) data-generating process
#
# Generates a stable VAR(p) where (~70%) of the slope coefficients are zero.
# The first lag's diagonal is dense (own persistence); off-diagonals and
# higher-lag coefficients are sparse.
#
# REPLACE WITH YOUR DATA: to use real macro/financial data instead, build
# a (T x M) numeric matrix `Yraw` with column names = variable names.
# =============================================================================

simulate_sparse_var <- function(M       = 5,
                                T       = 250,
                                p       = 4,
                                sparsity = 0.70,
                                seed    = 1) {
  set.seed(seed)

  K <- M * p
  A.true <- matrix(0, K, M)

  # Start with random coefficients
  raw <- matrix(rnorm(K * M, 0, 0.3), K, M)

  # Sparse mask: each entry independently zero with probability `sparsity`
  mask <- matrix(rbinom(K * M, 1, 1 - sparsity), K, M)
  A.true <- raw * mask

  # Inject own-persistence on lag 1 (diagonal) -- keep dense and shrunken
  diag(A.true[1:M, 1:M]) <- runif(M, 0.4, 0.7)

  # Decay higher-lag entries to ensure stability
  for (k in 2:p) {
    rows <- ((k - 1) * M + 1):(k * M)
    A.true[rows, ] <- A.true[rows, ] / (k^1.5)
  }

  # Check / enforce stability via companion form
  companion <- function(A, M, p) {
    K <- M * p
    MM <- matrix(0, K, K)
    MM[1:M, ] <- t(A)
    if (p > 1) MM[(M + 1):K, 1:(K - M)] <- diag(K - M)
    MM
  }
  attempt <- 0
  while (attempt < 50) {
    MM <- companion(A.true, M, p)
    rho <- max(abs(eigen(MM, only.values = TRUE)$values))
    if (rho < 0.95) break
    A.true <- A.true * 0.9                # shrink toward zero
    diag(A.true[1:M, 1:M]) <- diag(A.true[1:M, 1:M]) * 0.95
    attempt <- attempt + 1
  }

  # Error covariance: mildly correlated
  Sigma.true <- 0.5 * diag(M) + 0.2
  diag(Sigma.true) <- 1
  chol.S <- chol(Sigma.true)

  # Burn-in then sample
  Tburn <- 200
  Y <- matrix(0, T + Tburn + p, M)
  for (t in (p + 1):nrow(Y)) {
    xlag <- as.numeric(t(Y[(t - 1):(t - p), , drop = FALSE]))
    # xlag is c(Y_{t-1}, Y_{t-2}, ...): rearrange to match mlag layout
    Y[t, ] <- as.numeric(xlag %*% A.true) + rnorm(M) %*% chol.S
  }
  Yraw <- Y[(Tburn + 1):nrow(Y), , drop = FALSE]
  colnames(Yraw) <- paste0("y", 1:M)

  list(Yraw       = Yraw,
       A.true     = A.true,                # K x M coefficient matrix
       Sigma.true = Sigma.true,
       mask.true  = (A.true != 0) * 1)     # K x M zero/one pattern
}
