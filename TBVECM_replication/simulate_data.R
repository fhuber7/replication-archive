# ============================================================================
# Simulate a three-regime threshold VECM data set.
#
# The data-generating process produces M = 3 cointegrated series of length
# T with one cointegrating relation y_{1t} - y_{2t} + 0.5 * y_{3t}. The
# adjustment speed alpha and the short-run autoregressive matrix switch
# across three regimes that are determined by the level of the
# (demeaned) cointegrating residual.
#
# REPLACE WITH YOUR DATA: in main_replication.R the call to
# simulate_tbvecm_data(...) is the entry point. Swap it for whatever
# loads your (T x M) matrix of levels of the endogenous variables.
# ============================================================================

simulate_tbvecm_data <- function(T            = 300,
                                 M            = 3,
                                 seed         = 42,
                                 burn         = 100,
                                 thresholds   = c(-0.4, 0.4),
                                 sigma_scale  = 0.30) {

  if (M != 3) {
    stop("simulate_tbvecm_data is calibrated for M = 3; extend manually.")
  }
  set.seed(seed)

  # Cointegrating vector (linearly normalised: leading element = 1).
  beta_true  <- matrix(c(1, -1, 0.5), nrow = M, ncol = 1)

  # Three sets of regime-specific adjustment loadings (alpha) and AR1
  # coefficients (Phi1). All chosen so that the implied companion matrix
  # has spectral radius < 1 in every regime.
  alpha_list <- list(
    matrix(c(-0.10, -0.05,  0.02), 1, M),   # Regime 1: weak adjustment
    matrix(c(-0.25,  0.05, -0.05), 1, M),   # Regime 2: medium adjustment
    matrix(c(-0.45,  0.10, -0.08), 1, M)    # Regime 3: strong adjustment
  )
  Phi_list <- list(
    diag(c( 0.35,  0.40,  0.30)),
    diag(c( 0.10,  0.15,  0.20)),
    diag(c(-0.10, -0.05,  0.05))
  )
  # Per-regime innovation covariances (constant correlation, varying
  # scale across regimes).
  R <- matrix(c(1.0, 0.30, -0.10,
                0.30, 1.0,  0.20,
               -0.10, 0.20, 1.0), M, M)
  scales <- c(0.8, 1.0, 1.4) * sigma_scale
  Sigma_list <- lapply(scales, function(s) (s^2) * R)

  Tfull <- T + burn
  Y     <- matrix(0, Tfull, M)
  st    <- integer(Tfull)
  # Seed with a small Gaussian innovation
  Y[1, ] <- rnorm(M, 0, sigma_scale)
  Y[2, ] <- Y[1, ] + rnorm(M, 0, sigma_scale)

  for (t in 3:Tfull) {
    # Threshold variable: lagged cointegrating residual (demeaned in the
    # limit, but for simulation we just use the raw value).
    z <- as.numeric(Y[t - 1, , drop = FALSE] %*% beta_true)
    s <- if (z < thresholds[1])      1L
         else if (z < thresholds[2]) 2L
         else                        3L
    st[t]  <- s
    alpha  <- alpha_list[[s]]
    Phi    <- Phi_list[[s]]
    Sigma  <- Sigma_list[[s]]

    eps    <- as.vector(t(chol(Sigma)) %*% rnorm(M))
    dY_lag <- Y[t - 1, ] - Y[t - 2, ]
    dY     <- as.vector(t(alpha) %*% (t(beta_true) %*% Y[t - 1, ])) +
              as.vector(Phi %*% dY_lag) + eps
    Y[t, ] <- Y[t - 1, ] + dY
  }

  Y  <- Y[(burn + 1):Tfull, , drop = FALSE]
  st <- st[(burn + 1):Tfull]
  colnames(Y) <- c("y1", "y2", "y3")

  list(
    Y          = Y,
    states     = st,
    beta_true  = beta_true,
    alpha_list = alpha_list,
    Phi_list   = Phi_list,
    Sigma_list = Sigma_list,
    thresholds = thresholds
  )
}
