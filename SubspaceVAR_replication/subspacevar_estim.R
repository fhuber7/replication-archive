# ============================================================================
# subspacevar_estim.R
#
# Self-contained estimator for the Subspace-Shrinkage Conjugate Bayesian VAR
# of Huber & Koop, "Subspace Shrinkage in Conjugate Bayesian VARs"
# (Journal of Applied Econometrics).
#
# Contents
#   get.dum()        : Minnesota-style dummy-observation prior matrices.
#   get_companion()  : Companion-form helper (used to check stationarity).
#   impulsdtrf()     : Sims (2004) recursive IRF routine.
#   subVAR.0()       : Headline conjugate symmetric estimator with subspace
#                      shrinkage toward a principal-component restriction.
#                      Returns the posterior of (A, Sigma), the posterior
#                      weights over (omega, q, theta), recursive Cholesky
#                      IRFs and predictive draws.
#
# Notes
#   * The estimator is conjugate matrix-Normal/Inverse-Wishart conditional on
#     the hyperparameters (omega, q, theta). Hyperparameters are integrated
#     out via marginal-likelihood weights on a fixed grid (inverse-transform
#     sampling).
#   * omega is the weight on the PC-subspace prior; (1-omega) is the weight
#     on the diffuse / Minnesota prior. q is the number of leading principal
#     components defining the subspace. theta is the standard Minnesota
#     overall tightness.
#   * Two prior modes:
#       prior = "Minn"  -> Minnesota dummy observations are part of the prior
#       prior = "Flat"  -> diffuse base prior (Minnesota dummies essentially
#                          turned off via a large theta).
# ============================================================================

# ---------------------------------------------------------------------------
# Minnesota dummy-observation matrices (Banbura, Giannone, Reichlin style)
# ---------------------------------------------------------------------------
get.dum <- function(theta, gamma.prior, delta, mysigma, M, p) {
  ydummy <- matrix(0, 2 * M + M * (p - 1) + 1, M)
  xdummy <- matrix(0, 2 * M + M * (p - 1) + 1, M * p + 1)

  ydummy[1:M, ] <- diag((as.numeric(mysigma) * delta) / theta)
  ydummy[(M * (p - 1) + M + 1):(M * (p - 1) + 2 * M), ] <-
    diag(as.numeric(mysigma))

  jp <- diag(1:p)
  xdummy[1:(M * p), 1:(M * p)] <-
    kronecker(jp, diag(as.numeric(mysigma))) / theta
  xdummy[nrow(xdummy), ncol(xdummy)] <- gamma.prior
  list(Xdum = xdummy, Ydum = ydummy)
}

# ---------------------------------------------------------------------------
# Companion form helper. varndxv = c(M, n_deterministic, p).
# ---------------------------------------------------------------------------
get_companion <- function(Beta_, varndxv) {
  nn <- varndxv[[1]]; nd <- varndxv[[2]]; nl <- varndxv[[3]]
  nkk <- nn * nl + nd
  Jm <- matrix(0, nkk, nn); Jm[1:nn, 1:nn] <- diag(nn)
  MM <- rbind(
    t(Beta_),
    cbind(diag((nl - 1) * nn), matrix(0, (nl - 1) * nn, nn + 1)),
    c(matrix(0, 1, nn * nl), 1)
  )
  list(MM = MM, Jm = Jm)
}

# ---------------------------------------------------------------------------
# Sims (2004) recursive impulse-response routine.
# B     : M x M x p array of reduced-form VAR coefficients.
# smat  : M x M shock-impact matrix (e.g. t(chol(Sigma))).
# nstep : horizon.
# Returns an M x M x nstep array of responses.
# ---------------------------------------------------------------------------
impulsdtrf <- function(B, smat, nstep) {
  neq  <- dim(B)[1]; nvar <- dim(B)[2]; lags <- dim(B)[3]
  if (dim(smat)[2] != dim(B)[2])
    stop("B and smat conflict on # of variables")
  response <- array(0, dim = c(neq, nvar, nstep + lags - 1))
  response[, , lags] <- smat
  response <- aperm(response, c(1, 3, 2))
  irhs <- 1:(lags * nvar); ilhs <- lags * nvar + (1:nvar)
  response <- matrix(response, ncol = neq)
  B <- B[, , seq(from = lags, to = 1, by = -1)]
  B <- matrix(B, nrow = nvar)
  for (it in 1:(nstep - 1)) {
    response[ilhs, ] <- B %*% response[irhs, ]
    irhs <- irhs + nvar; ilhs <- ilhs + nvar
  }
  dim(response) <- c(nvar, nstep + lags - 1, nvar)
  if (lags > 1) {
    response <- response[, -(1:(lags - 1)), , drop = FALSE]
  }
  response <- aperm(response, c(1, 3, 2))
  response
}

# ---------------------------------------------------------------------------
# subVAR.0 : symmetric conjugate VAR with PC-subspace shrinkage
#
# Y           : T x M data matrix.
# p           : lag order.
# prior       : "Minn" (Minnesota base) or anything else (diffuse base).
# grid.omega  : grid of subspace weights in (0,1).
# grid.q      : grid of integer subspace dimensions.
# grid.theta  : grid of Minnesota tightness values.
# a.omega,b.omega : Beta-prior shape parameters on omega.
# fhor        : predictive draws horizon (0 = none).
# irfhor      : IRF horizon (0 = none).
# nsave,nburn : number of posterior draws / burn-in.
#
# Returns a list with posterior draws, hyperparameter posteriors, IRFs and
# the per-grid-point log-marginal-likelihood weights.
# ---------------------------------------------------------------------------
subVAR.0 <- function(nsave, nburn, Y, p = 2, prior = "Minn",
                     set.omega = 0.01, sample.omega = TRUE,
                     grid.omega = seq(0.01, 0.99, 0.1),
                     theta = 0.1,
                     a.omega = 25, b.omega = 1.5 / ncol(Y),
                     a.mean = 0, s.prior = 0.1,
                     fhor = 4, irfhor = 20,
                     grid.q = seq(1, 8),
                     grid.theta = c(0.001, 0.01, 0.1, 0.25, 0.5, 1),
                     a.q = 10, verbose = TRUE) {
  ntot <- nsave + nburn
  M <- ncol(Y)

  # Step A: AR(p) residual SDs (Minnesota scale)
  mysigma <- matrix(NA, M, 1)
  for (ii in seq_len(M)) {
    tmpar <- arima(Y[, ii], order = c(p, 0, 0),
                   method = "CSS", include.mean = FALSE)
    mysigma[ii, 1] <- sqrt(tmpar$sigma2)
  }

  # Design and dummy-observation matrices
  X <- embed(Y, dimension = p + 1)
  X <- cbind(X[, (M + 1):ncol(X)], 1)
  y <- Y[(p + 1):nrow(Y), ]

  dums  <- get.dum(theta, 0.01, a.mean, mysigma, M, p)
  y.tilde <- rbind(y, dums$Ydum)
  X.tilde <- rbind(X, dums$Xdum)

  if (prior != "Minn") {
    grid.theta <- 10^2
  }

  T <- nrow(X); K <- ncol(X); n <- K * M
  v.prior <- nrow(y.tilde) - nrow(y)
  S.prior <- diag(M) * s.prior

  if (!sample.omega) grid.omega <- set.omega
  grid <- expand.grid(grid.omega, grid.q, grid.theta)
  colnames(grid) <- c("omega", "q", "theta")

  V.post.array <- array(NA, c(K, K, nrow(grid), 2))
  S.post.array <- array(NA, c(M, M, nrow(grid), 2))
  A.post.array <- array(NA, c(K, M, nrow(grid)))
  error.matrix <- array(NA, c(nrow(grid), 1))

  ml.omega <- matrix(NA, nrow(grid), 1)
  T.diag   <- diag(nrow(X.tilde))
  pca.X    <- prcomp(X)
  a.star   <- max(1, sum(pca.X$sdev > 1))

  # Per-grid-point posterior moments and log marginal likelihood
  for (jj in seq_len(nrow(grid))) {
    omega.j    <- grid[jj, 1]
    q.slct     <- grid[jj, 2]
    theta.slct <- grid[jj, 3]

    dums <- get.dum(theta.slct, 1e-5, a.mean, mysigma, M, p)
    X.tilde <- rbind(X, dums$Xdum)
    y.tilde <- rbind(y, dums$Ydum)

    pca <- svd(X)
    if (q.slct == 1) {
      eig.mat <- pca$d[1]
      Psi00   <- pca$u[, 1, drop = FALSE] * eig.mat
    } else {
      eig.mat <- diag(pca$d[1:q.slct])
      Psi00   <- pca$u[, 1:q.slct, drop = FALSE] %*% eig.mat
    }
    Psi0 <- matrix(0, T + v.prior, q.slct)
    Psi0[1:T, ] <- Psi00
    T.diag[(T + 1):nrow(Psi0), (T + 1):nrow(Psi0)] <- 0

    Q0 <- Psi0 %*% solve(crossprod(Psi0) + diag(q.slct) * 1e-7) %*% t(Psi0)

    Prior.part <- omega.j / (1 - omega.j) *
      (t(X.tilde) %*% (T.diag - Q0) %*% X.tilde)
    V.post.inv <- crossprod(X.tilde) + Prior.part + diag(K) * 1e-7
    V.post     <- solve(V.post.inv)
    Prior.part.1 <- Prior.part + crossprod(dums$Xdum)

    if (prior == "Minn") {
      A.prior <- matrix(0, K, M)
      A.prior[1:M, 1:M] <- diag(M) * a.mean
      S.prior <- crossprod(dums$Ydum - dums$Xdum %*% A.prior)
    } else {
      A.prior <- matrix(0, K, M)
    }

    V.post.chol <- try(t(chol(V.post)), silent = TRUE)
    if (is(V.post.chol, "try-error"))
      V.post.chol <- t(chol(V.post, pivot = TRUE))

    A.post <- V.post %*% crossprod(X.tilde, y.tilde)
    S.post <- S.prior + crossprod(y) - t(A.post) %*% V.post.inv %*% A.post

    log.det.V <- -2 * sum(log(diag(t(chol(Prior.part.1)))))
    log.det.K <- -2 * sum(log(diag(V.post.chol)))
    ml.i <- -M / 2 * (log.det.V + log.det.K) -
      (T + v.prior) / 2 * 2 * sum(log(diag(chol(S.post)))) +
      log(dbeta(omega.j, a.omega, b.omega)) +
      log(dgamma(theta.slct, 0.2 * 0.5, 0.5))

    # Stationarity check on companion form
    MM <- get_companion(A.post, c(M, 1, p))$MM
    if (max(abs(eigen(MM)$values)) > 1.05) ml.i <- -1e5
    if (is.nan(ml.i) || is.infinite(ml.i))  ml.i <- -1e4

    A.post.array[, , jj]      <- A.post
    V.post.array[, , jj, 1]   <- V.post
    V.post.array[, , jj, 2]   <- V.post.chol
    S.post.array[, , jj, 1]   <- S.post
    S.post.array[, , jj, 2]   <- solve(S.post)
    ml.omega[jj, ]            <- ml.i

    post.fit   <- X.tilde %*% V.post %*% t(X.tilde) %*% y.tilde
    err.approx <- mean((post.fit[1:T, ] -
                        (1 - omega.j) * X %*% solve(crossprod(X.tilde)) %*%
                          crossprod(X.tilde, y.tilde) -
                        omega.j * (Q0 %*% y.tilde)[1:T, ])^2)
    error.matrix[jj, ] <- err.approx

    if (verbose) cat(sprintf("grid %d/%d  ml=%.2f\n",
                             jj, nrow(grid), ml.i))
  }

  # Posterior weights via inverse-transform sampling on the grid
  normalizer <- max(ml.omega)
  weights    <- exp(ml.omega - normalizer) / sum(exp(ml.omega - normalizer))
  best       <- which(ml.omega == max(ml.omega))[1]

  V.chol      <- V.post.array[, , best, 2]
  S.post.inv  <- S.post.array[, , best, 2]
  A.post      <- A.post.array[, , best]
  Sigma.draw  <- crossprod(y - X %*% A.post) / T

  # Storage
  A.store      <- array(NA, c(nsave, K, M))
  Sigma.store  <- array(NA, c(nsave, M, M))
  omega.store  <- array(NA, c(nsave, 1))
  q.store      <- array(NA, c(nsave, 1))
  theta.store  <- array(NA, c(nsave, 1))
  pred.store   <- if (fhor   > 0) array(NA, c(nsave, M, fhor))         else NULL
  IRF.store    <- if (irfhor > 0) array(NA, c(nsave, M, M, irfhor))    else NULL

  Sig.chol <- chol(Sigma.draw, pivot = TRUE)
  for (irep in seq_len(ntot)) {
    # Matrix-Normal coefficient draw (vec trick)
    Alpha.draw <- matrix(as.vector(A.post) +
                         as.vector(V.chol %*%
                                     matrix(rnorm(n), K, M) %*% Sig.chol),
                         K, M)
    # Inverse-Wishart variance draw
    Sigma.draw.inv <- matrix(rWishart(1, nrow(y) + v.prior, S.post.inv), M, M)
    Sigma.draw     <- solve(Sigma.draw.inv)
    Sig.chol       <- chol(Sigma.draw, pivot = TRUE)

    # Resample (omega, q, theta) from posterior weights on the grid
    ind.sample <- sample(seq_len(nrow(grid)), size = 1, prob = weights)
    V.chol     <- V.post.array[, , ind.sample, 2]
    A.post     <- A.post.array[, , ind.sample]
    S.post.inv <- S.post.array[, , ind.sample, 2]

    if (irep > nburn) {
      ii <- irep - nburn
      A.store[ii, , ]      <- Alpha.draw
      Sigma.store[ii, , ]  <- Sigma.draw
      omega.store[ii, ]    <- grid[ind.sample, 1]
      q.store[ii, ]        <- grid[ind.sample, 2]
      theta.store[ii, ]    <- grid[ind.sample, 3]

      if (fhor > 0) {
        lag.idx <- seq_len(p * M - M)        # empty when p == 1
        XT <- c(Y[nrow(Y), ], X[nrow(X), lag.idx], 1)
        Y.mat <- matrix(NA, M, fhor)
        for (h in seq_len(fhor)) {
          Y.pred <- t(XT %*% Alpha.draw) + t(Sig.chol) %*% rnorm(M)
          XT <- c(as.numeric(Y.pred), XT[lag.idx], 1)
          Y.mat[, h] <- Y.pred
        }
        pred.store[ii, , ] <- Y.mat
      }

      if (irfhor > 0) {
        PHI_array <- array(0, c(M, M, p))
        for (ss in seq_len(p))
          PHI_array[, , ss] <- t(Alpha.draw[((ss - 1) * M + 1):(ss * M), ])
        IRF.store[ii, , , ] <- impulsdtrf(PHI_array, t(Sig.chol), irfhor)
      }
    }
  }

  list(
    A_posterior      = A.store,
    Sigma_posterior  = Sigma.store,
    omega_posterior  = omega.store,
    q_posterior      = q.store,
    theta_posterior  = theta.store,
    IRF_posterior    = IRF.store,
    pred_density     = pred.store,
    grid             = grid,
    weights          = weights,
    log_ml           = ml.omega,
    A_array          = A.post.array,
    error            = error.matrix
  )
}
