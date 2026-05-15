# =============================================================================
# simulate_data.R
# -----------------------------------------------------------------------------
# Generate a small two-regime mixture-of-VAR data set so that the BNP-VAR
# sampler in `Codes_and_data/aux.R` (function `bvar.mix`) has something
# non-Gaussian to fit. The DGP is a stationary VAR(p) backbone whose
# innovations are drawn from a two-component Gaussian mixture with very
# different scales — exactly the kind of shock structure the BNP prior on the
# error distribution is designed to capture.
#
# REPLACE WITH YOUR DATA: in main_replication.R replace the call to
# simulate_bnpvar_data(...) with code that assigns a numeric (T x M) matrix
# to `Xraw` *before* the apply(.., scale) standardisation.
# =============================================================================

#' Simulate a mixture-of-shocks VAR
#'
#' @param T            Number of retained observations.
#' @param M            Number of variables.
#' @param p            Lag length.
#' @param burnin       Burn-in observations (discarded).
#' @param mix.probs    Mixture weights for the two shock components.
#' @param scale.mult   Scale multiplier for the heavy-tailed component.
#' @param seed         RNG seed.
#' @return             List with `Y` (T x M matrix), `var.names`, and the
#'                     latent allocation indicator `regime`.
simulate_bnpvar_data <- function(T = 200, M = 5, p = 2,
                                 burnin = 200,
                                 mix.probs = c(0.85, 0.15),
                                 scale.mult = 3,
                                 seed = 1) {
  set.seed(seed)
  Tt <- T + burnin

  # Stable VAR(p) backbone
  A_list <- vector("list", p)
  for (l in seq_len(p)) {
    A <- matrix(rnorm(M * M, 0, 0.05), M, M)
    diag(A) <- 0.55 / l
    A_list[[l]] <- A
  }

  # Two contemporaneous covariance components: "calm" and "turbulent"
  L1 <- diag(M);  L1[lower.tri(L1)] <- rnorm(M*(M-1)/2, 0, 0.25)
  L2 <- diag(M) * scale.mult
  L2[lower.tri(L2)] <- rnorm(M*(M-1)/2, 0, 0.6) * scale.mult

  regime <- sample(1:2, Tt, replace = TRUE, prob = mix.probs)

  Y <- matrix(0, Tt, M)
  for (t in (p + 1):Tt) {
    lin <- numeric(M)
    for (l in seq_len(p)) lin <- lin + A_list[[l]] %*% Y[t - l, ]
    L <- if (regime[t] == 1) L1 else L2
    e <- as.numeric(L %*% rnorm(M))
    Y[t, ] <- as.numeric(lin) + e
  }

  Y      <- Y[(burnin + 1):Tt, , drop = FALSE]
  regime <- regime[(burnin + 1):Tt]
  colnames(Y) <- paste0("y", seq_len(M))
  # Truth retained for overlay in main_replication.R:
  #   regime: T-length latent component indicator (1 = calm, 2 = turbulent)
  #   A_list / L1 / L2: linear backbone and the two regime Cholesky factors
  #   G_true: number of mixture components in the DGP
  list(Y         = Y,
       var.names = colnames(Y),
       regime    = regime,
       A_list    = A_list,
       L1        = L1,
       L2        = L2,
       G_true    = 2L,
       mix.probs = mix.probs)
}
