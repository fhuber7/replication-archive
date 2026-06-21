# ============================================================================
# simulate_data.R
# ----------------------------------------------------------------------------
# Synthetic data generator for the TVP-BART replication package.
#
# Generates an M-variable, p = 1 reduced-form VAR
#
#   y_t = A_t y_{t-1} + eps_t,     eps_t ~ N(0, Sigma_t)
#
# in which the autoregressive matrix drifts *smoothly* over the sample as a
# single common factor:
#
#   A_t = (1 - g_t) A_start + g_t A_end,     g_t = logistic transition in t.
#
# Because the whole coefficient surface is driven by the scalar g_t, the time
# variation is rank-1 -- exactly the situation the latent-factor BART of the
# paper is designed to compress (recoverable with Q = 1 factor).  The error
# covariance is also time-varying: a constant lower-triangular contemporaneous
# matrix L scales stochastic-volatility paths sigma_{i,t} that hump up in the
# middle of the sample (a stylised "crisis").  Both ingredients make the
# implied impulse responses time-varying.
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  REPLACE WITH YOUR DATA  <<<<<<<<<<<<<<<<<<<<<
# In a real application, replace simulate_tvpbart_data() with code returning a
# (T x M) matrix Yraw with column names.  The `truth` list is only used by the
# driver to overlay the data-generating coefficients / IRFs and can be ignored.
# ============================================================================

simulate_tvpbart_data <- function(T    = 200,
                                  M    = 3,
                                  p    = 1,
                                  burn = 200,
                                  seed = 42) {
  stopifnot(M == 3, p == 1)
  set.seed(seed)
  Ttot <- T + burn

  # ---- two stationary AR(1) regimes (the drift between them is rank-1) -----
  # y1 loses persistence (0.70 -> 0.10) and its loading on y2 flips sign
  # (0.30 -> -0.30); y2 gains persistence (0.20 -> 0.70); y3 stays put.
  A_start <- matrix(c( 0.70,  0.30, 0.00,
                       0.00,  0.20, 0.00,
                       0.10,  0.00, 0.40), M, M, byrow = TRUE)
  A_end   <- matrix(c( 0.10, -0.30, 0.00,
                       0.00,  0.70, 0.00,
                       0.10,  0.00, 0.40), M, M, byrow = TRUE)

  # smooth logistic transition over the *kept* sample, padded across burn-in
  tt  <- seq_len(Ttot)
  mid <- burn + T / 2
  g   <- 1 / (1 + exp(-0.05 * (tt - mid)))          # 0 -> 1 transition

  # ---- time-varying recursive covariance -----------------------------------
  L <- matrix(c( 1.0, 0.0, 0.0,
                 0.3, 1.0, 0.0,
                -0.2, 0.4, 1.0), M, M, byrow = TRUE)   # constant contemporaneous
  base_sd <- c(0.40, 0.35, 0.30)
  hump    <- exp(0.35 * dnorm(tt, mean = mid, sd = T / 6) / dnorm(mid, mid, T / 6))
  sigma   <- outer(hump, base_sd)                       # Ttot x M structural sd

  # ---- simulate ------------------------------------------------------------
  Y <- matrix(0, Ttot, M)
  A_all <- array(0, c(Ttot, M, M))
  for (t in 2:Ttot) {
    A_t <- (1 - g[t]) * A_start + g[t] * A_end
    A_all[t, , ] <- A_t
    eps <- L %*% (sigma[t, ] * rnorm(M))
    Y[t, ] <- A_t %*% Y[t - 1, ] + eps
  }

  keep  <- (burn + 1):Ttot
  Yraw  <- Y[keep, , drop = FALSE]
  colnames(Yraw) <- paste0("y", seq_len(M))

  truth <- list(
    A      = A_all[keep, , , drop = FALSE],   # T x M x M reduced-form lag matrix
    sigma  = sigma[keep, , drop = FALSE],     # T x M structural shock sd
    L      = L,                               # contemporaneous impact (constant)
    g      = g[keep],                         # the common transition factor
    A_start = A_start, A_end = A_end
  )

  list(Yraw = Yraw, truth = truth)
}

# True time-t Cholesky impulse responses implied by `truth`, for overlaying on
# the estimated IRFs.  Returns an (T x M x nhor) array of responses to a
# `shk.sc`-sd shock in variable `shk.var`.
true_irf <- function(truth, shk.var = 1, shk.sc = 1, nhor = 17) {
  Tn <- dim(truth$A)[1]; M <- dim(truth$A)[2]
  out <- array(NA_real_, c(Tn, M, nhor))
  for (t in seq_len(Tn)) {
    A  <- truth$A[t, , ]
    P  <- truth$L %*% diag(truth$sigma[t, ], M)     # lower-tri Cholesky of Sigma_t
    Fh <- diag(M)
    for (h in seq_len(nhor)) { out[t, , h] <- (Fh %*% P)[, shk.var] * shk.sc; Fh <- Fh %*% A }
  }
  out
}
