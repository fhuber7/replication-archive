# =============================================================================
# simulate_data.R
# -----------------------------------------------------------------------------
# Generate a single-equation (y, X) training set suitable for the GP-subspace
# and BART samplers in this folder. The DGP is a smooth nonlinear regression
# in many "spurious" predictors, plus a mild AR(1) component in y. This gives
# the subspace prior something to do (most predictors are irrelevant) and the
# Gaussian-process / BART models a smooth mean function to learn.
#
# REPLACE WITH YOUR DATA: in main_replication.R replace the call to
# simulate_subspace_data(...) with code that constructs your own pair (y, X),
# train / holdout split, and standardises y to zero mean / unit variance and
# X column-wise. The paper's exercise targets CPIAUCSL or CPILFESL with the
# FRED-QD covariates listed in infdata_script.R.
# =============================================================================

#' Simulate a sparse nonlinear regression with mild persistence in y
#'
#' @param T          Number of observations (train + holdout).
#' @param p_cov      Number of candidate predictors in X.
#' @param k_active   Number of truly active predictors (k_active <= p_cov).
#' @param sigma      Noise s.d.
#' @param seed       RNG seed.
#' @return           List with `y` (T x 1 matrix), `X` (T x p_cov matrix),
#'                   `active` (indices of active predictors) and `var.names`.
simulate_subspace_data <- function(T = 200, p_cov = 15, k_active = 4,
                                   sigma = 0.5, seed = 1) {
  set.seed(seed)

  # Predictors: a low-rank persistent block + iid noise columns
  q <- max(3, k_active)
  F_t <- matrix(0, T, q)                    # latent factors
  for (j in seq_len(q)) {
    F_t[1, j] <- rnorm(1)
    for (t in 2:T) F_t[t, j] <- 0.8 * F_t[t - 1, j] + rnorm(1, 0, 0.6)
  }
  Lambda <- matrix(rnorm(p_cov * q, 0, 0.4), p_cov, q)
  X <- F_t %*% t(Lambda) + matrix(rnorm(T * p_cov, 0, 0.3), T, p_cov)
  colnames(X) <- paste0("x", seq_len(p_cov))

  # Pick the active predictors and build a smooth nonlinear mean
  active <- sort(sample(seq_len(p_cov), k_active))
  Xa <- X[, active, drop = FALSE]
  f_mean <- 0.6 * tanh(Xa[, 1]) +
            0.5 * sin(Xa[, 2]) +
            0.3 * (Xa[, 3])^2 / (1 + abs(Xa[, 3]))
  if (k_active >= 4) f_mean <- f_mean + 0.4 * (Xa[, 4] * Xa[, 1])

  # y with mild AR(1) component
  y <- numeric(T)
  y[1] <- f_mean[1] + rnorm(1, 0, sigma)
  for (t in 2:T) y[t] <- 0.3 * y[t - 1] + f_mean[t] + rnorm(1, 0, sigma)

  # True conditional mean of y given x_t (excluding the AR(1) component).
  # Useful for overlaying the in-sample fit truth.
  y_true_mean <- numeric(T)
  y_true_mean[1] <- f_mean[1]
  for (t in 2:T) y_true_mean[t] <- 0.3 * y[t - 1] + f_mean[t]

  list(y           = matrix(y, ncol = 1),
       X           = X,
       active      = active,
       var.names   = colnames(X),
       f_mean      = f_mean,
       y_true_mean = y_true_mean,
       sigma       = sigma)
}
