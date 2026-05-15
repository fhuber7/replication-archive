## =============================================================================
## pvar_irga_estim.R
##
## Self-contained estimator for the panel VAR with Integrated Rotated Gaussian
## Approximation (IRGA), per Feldkircher, Huber, Koop, Pfarrhofer:
## "Approximate Bayesian inference and forecasting in huge-dimensional
##  multi-country VARs".
##
## Strategy (per country i, per equation j):
##   Y_{ij} = X_i beta + Z_i gamma + eps,
## where X_i are the lagged own-country variables (the focal "domestic" block)
## and Z_i are everything else (lags of OTHER countries' variables plus the
## already-sampled contemporaneous own-country shocks for ordering).
##
## Step 1 (rotation):  QR-decompose X_i = Q R.  Premultiplying by Q^T splits
## the system into a p-dim "domestic" subsystem and an (n-p)-dim "nuisance"
## subsystem that depends only on gamma.
##
## Step 2 (VAMP/VB):   Approximate the posterior of gamma in the nuisance block
## by Vector Approximate Message Passing (VAMP) with a Horseshoe / Minnesota
## prior.  This integrates the cross-country block out analytically.
##
## Step 3 (MCMC):      Sample beta from the rotated p-dim regression
##   ytilde = Q^T (RtY - RtZ %*% gamma_hat),    Sigma = sigma^2 I_p + RtZ S RtZ^T
## with a Horseshoe prior on beta.  Output is a posterior sample of beta plus
## the VAMP mean of gamma.
##
## Dependencies: Matrix, MASS, mvtnorm, GIGrvg
## Avoided:      matrixcalc, fastSOM, lars (base R fallbacks used throughout)
## =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(MASS)
  library(mvtnorm)
  library(GIGrvg)
})

## ----------------------------------------------------------------------------
## Small utilities
## ----------------------------------------------------------------------------

mlag <- function(X, lag) {
  X <- as.matrix(X)
  Traw <- nrow(X); N <- ncol(X); p <- lag
  Xlag <- matrix(0, Traw, p * N)
  for (ii in 1:p)
    Xlag[(p + 1):Traw, (N * (ii - 1) + 1):(N * ii)] <-
      X[(p + 1 - ii):(Traw - ii), 1:N, drop = FALSE]
  Xlag
}

# Inverse-gamma draw (avoids dependency on the invgamma package).
rinvgamma_ <- function(n, shape, rate) 1 / rgamma(n, shape = shape, rate = rate)

## ----------------------------------------------------------------------------
## Horseshoe-style posterior updates (Makalic & Schmidt parameterisation).
## ----------------------------------------------------------------------------

get.hs <- function(bdraw, lambda.hs, nu.hs, tau.hs, zeta.hs) {
  k <- length(bdraw)
  if (is.na(tau.hs)) {
    tau.hs <- 1
  } else {
    tau.hs <- rinvgamma_(1, shape = (k + 1) / 2,
                         rate = 1 / zeta.hs + sum(bdraw^2 / lambda.hs) / 2)
  }
  zeta.hs <- rinvgamma_(1, shape = 1, rate = 1 + 1 / tau.hs)
  lambda.hs <- rinvgamma_(k, shape = 1,
                          rate = 1 / nu.hs + bdraw^2 / (2 * tau.hs))
  nu.hs <- rinvgamma_(k, shape = 1, rate = 1 + 1 / lambda.hs)
  list(psi = lambda.hs * tau.hs, lambda = lambda.hs,
       tau = tau.hs, nu = nu.hs, zeta = zeta.hs)
}

get.minn.hs <- function(bdraw, lambda.hs, tau.hs, zeta.hs) {
  k <- length(bdraw)
  tau.hs <- rinvgamma_(1, shape = (k + 1) / 2,
                       rate = 1 / zeta.hs + sum(bdraw^2 / lambda.hs) / 2)
  zeta.hs <- rinvgamma_(1, shape = 1, rate = 1 + 1 / tau.hs)
  list(psi = lambda.hs * tau.hs, lambda = lambda.hs,
       tau = tau.hs, zeta = zeta.hs)
}

## ----------------------------------------------------------------------------
## Vector Approximate Message Passing (VAMP)
##
## Approximate posterior of alpha in:  y = X alpha + eps,  eps ~ N(0, sigma^2 I)
## with a Horseshoe prior on alpha.  This is "Step 2" of the IRGA algorithm.
## Output: posterior mean and (scalar) variance of alpha plus sigma^2.
## ----------------------------------------------------------------------------

VAMP <- function(y, X, lambda = "HS", psi = NULL, sigma.sq = NULL,
                 rho = 1, max.iter = 1000, min.rho = 1e-3,
                 a.0 = 3, b.0 = 0.03,
                 hs.lb = 1e-32, hs.ub = 0.1, verbose = FALSE) {
  y <- as.vector(y); X <- as.matrix(X)
  n <- length(y); p <- ncol(X)
  s2unknown <- is.null(sigma.sq0 <- sigma.sq)

  sv <- svd(X)
  U <- sv$u; diag_D <- sv$d; V <- sv$v
  UtY  <- crossprod(U, y)
  D_Vt <- diag_D * t(V)

  r <- numeric(p)
  tau.inv.hs <- 1; xi.inv.hs <- 1
  lambda.inv.hs <- rep(1, p); nu.inv.hs <- rep(1, p)

  if (lambda == "MN") {
    psi.vec <- psi; lambda.inv.hs <- 1 / psi
  } else {
    psi.vec <- rep(1, p)
  }
  t.sq <- 1

  for (k in 0:max.iter) {
    # 3a. shrinkage update --------------------------------------------------
    if (lambda == "HS") {
      scale.i       <- r^2
      lambda.inv.hs <- 1 / (nu.inv.hs + tau.inv.hs * scale.i / 2)
      tau.inv.hs    <- (p + 1) / (2 * xi.inv.hs + sum(lambda.inv.hs * scale.i))
      nu.inv.hs     <- 1 / (1 + lambda.inv.hs)
      xi.inv.hs     <- 1 / (1 + tau.inv.hs)
      psi.vec       <- (1 / as.numeric(lambda.inv.hs)) * (1 / tau.inv.hs)
    } else if (lambda == "MN") {
      scale.i    <- r^2
      tau.inv.hs <- (p + 1) / (2 * xi.inv.hs + sum(lambda.inv.hs * scale.i))
      xi.inv.hs  <- 1 / (1 + tau.inv.hs)
      psi.vec    <- (1 / as.numeric(lambda.inv.hs)) * (1 / tau.inv.hs)
    }
    psi.vec[psi.vec > hs.ub] <- hs.ub
    psi.vec[psi.vec < hs.lb] <- hs.lb

    beta.new <- r / (t.sq / psi.vec + 1)
    s.sq.new <- mean(1 / (1 / psi.vec + 1 / t.sq))

    if (k == 0) { beta <- beta.new; s.sq <- s.sq.new }
    else {
      beta <- rho * beta.new + (1 - rho) * beta
      s.sq <- rho * s.sq.new + (1 - rho) * s.sq
    }

    # 3b. ---------------------------------------------------------------
    denom_b <- (t.sq - s.sq)
    if (abs(denom_b) < 1e-12) denom_b <- sign(denom_b + 1e-300) * 1e-12
    t.tilde.sq <- 1 / (1 / s.sq - 1 / t.sq)
    r.tilde    <- (t.sq * beta - s.sq * r) / denom_b

    if (s2unknown)
      sigma.sq <- (2 * b.0 + sum((y - X %*% beta)^2)) / (2 * a.0 + n)

    # 3c. ---------------------------------------------------------------
    tmp <- 1 / (sigma.sq / diag_D / t.tilde.sq + diag_D)
    beta.tilde.new <- r.tilde + V %*% (tmp * (UtY - D_Vt %*% r.tilde))
    s.tilde.sq.new <- t.tilde.sq * (1 - sum(diag_D * tmp) / p)

    if (k == 0) { beta.tilde <- beta.tilde.new; s.tilde.sq <- s.tilde.sq.new }
    else {
      beta.tilde <- rho * beta.tilde.new + (1 - rho) * beta.tilde
      s.tilde.sq <- rho * s.tilde.sq.new + (1 - rho) * s.tilde.sq
    }

    # 3d. ---------------------------------------------------------------
    denom_t <- (1 / s.tilde.sq - 1 / t.tilde.sq)
    if (abs(denom_t) < 1e-12) denom_t <- sign(denom_t + 1e-300) * 1e-12
    t.sq <- 1 / denom_t
    t.sq <- min(max(t.sq, 1e-12), 1e12)

    denom_r <- (t.tilde.sq - s.tilde.sq)
    if (abs(denom_r) < 1e-12) denom_r <- sign(denom_r + 1e-300) * 1e-12
    r.old <- r
    r <- (t.tilde.sq * beta.tilde - s.tilde.sq * r.tilde) / denom_r

    Delta <- r - r.old
    relchg <- sum(Delta * Delta) / max(sum(r * r), 1e-12)
    if (is.na(relchg) || is.nan(relchg)) {
      if (rho <= min.rho) {
        warning("VAMP failed to converge: NA stopping criterion."); break
      }
      return(VAMP(y, X, lambda, psi, sigma.sq = sigma.sq0,
                  rho = rho / 2, max.iter = max.iter, min.rho = min.rho,
                  a.0 = a.0, b.0 = b.0,
                  hs.lb = hs.lb, hs.ub = hs.ub, verbose = verbose))
    }
    if (relchg < 1e-8) break
  }
  list(mean = as.vector(beta), variance = s.sq,
       sigma.sq = sigma.sq,
       lambda.hs = 1 / as.numeric(lambda.inv.hs),
       tau.hs = 1 / tau.inv.hs)
}

## ----------------------------------------------------------------------------
## Variational Bayes alternative for Step 2 (uses no SVD trick).
## Kept for completeness; the main estimator uses VAMP.
## ----------------------------------------------------------------------------

vblr_fit <- function(X, y, lambda = 1, psi = 1, a_0 = 3, b_0 = 0.03,
                     max_iter = 500, epsilon_conv = 1e-3,
                     shrink.type = "HS", is_verbose = FALSE,
                     hs.lb = 1e-32, hs.ub = NULL) {
  X <- as.matrix(X)
  D <- ncol(X); I <- nrow(X)
  XX <- crossprod(X, X); Xy <- crossprod(X, y)

  V_prior <- diag(0.01, D); V_prior_inv <- diag(1 / diag(V_prior))
  L <- numeric(max_iter)

  tau.inv.hs <- 1; xi.inv.hs <- 1
  nu.inv.hs <- rep(1, D)
  if (shrink.type == "HS") lambda.inv.hs <- rep(1, D)
  else if (shrink.type == "MN") lambda.inv.hs <- 1 / psi
  prior.mean <- rep(0, D)

  for (i in 2:max_iter) {
    S <- solve(V_prior_inv + lambda * XX)
    m <- S %*% (lambda * Xy + diag(V_prior_inv) * prior.mean)

    a.lambda <- a_0 + I / 2
    b.lambda <- b_0 + sum((y - X %*% m)^2) / 2 + sum(diag(XX %*% S)) / 2
    lambda <- a.lambda / b.lambda

    if (shrink.type == "HS") {
      scale.i       <- diag(S) + (m - prior.mean)^2
      lambda.inv.hs <- 1 / (nu.inv.hs + lambda * tau.inv.hs * scale.i / 2)
      tau.inv.hs    <- (D + 1) / (2 * xi.inv.hs + lambda * sum(lambda.inv.hs * scale.i))
      nu.inv.hs     <- 1 / (1 + lambda.inv.hs)
      xi.inv.hs     <- 1 / (1 + tau.inv.hs)
    } else if (shrink.type == "MN") {
      scale.i    <- diag(S) + (m - prior.mean)^2
      tau.inv.hs <- (D + 1) / (2 * xi.inv.hs + lambda * sum(lambda.inv.hs * scale.i))
      xi.inv.hs  <- 1 / (1 + tau.inv.hs)
    }
    V_prior <- 1 / lambda * diag((1 / as.numeric(lambda.inv.hs)) * 1 / tau.inv.hs)
    cap <- 1 / D^2
    V_prior[V_prior > cap]   <- cap
    V_prior[V_prior < 1e-32] <- 1e-32
    V_prior_inv <- diag(1 / diag(V_prior))

    L[i] <- mean(m)
    if (abs(L[i] - L[i - 1]) < epsilon_conv) break
  }
  structure(list(m = m, S = S, lambda = lambda,
                 X = X, I = I, D = D, sigma = 1 / lambda),
            class = "vblr")
}

## ----------------------------------------------------------------------------
## Step 1+2 wrapper: rotate by QR(X), then approximate alpha via VAMP or VB.
## ----------------------------------------------------------------------------

BVS_IRGA_12 <- function(y, X, Z, lambda, psi, sigma.sq = NULL,
                        max.iter = 1000, rho = 1, min.rho = 1e-3,
                        estimator = "VAMP", hs.lb = 1e-32, hs.ub = 0.1) {
  p <- ncol(X); n <- length(y)

  Q   <- qr.Q(qr(X), complete = TRUE)
  QtY <- as.vector(crossprod(Q, y))
  RtY <- QtY[1:p]
  RtX <- crossprod(Q[, 1:p, drop = FALSE], X)
  QtZ <- crossprod(Q, Z)
  RtZ <- QtZ[1:p, , drop = FALSE]

  if (p == n) {
    if (is.null(sigma.sq))
      stop("sigma.sq must be supplied when p == n.")
    covariance <- lambda * psi * tcrossprod(RtZ)
    return(list(RtY = RtY, RtX = RtX, mean = numeric(p),
                covariance = covariance, sigma.sq = sigma.sq,
                beta.mean = numeric(ncol(Z)),
                cov.S = matrix(0, 0, 0),
                tau.HS = NA, lambda.HS = NA))
  }

  StY <- QtY[-(1:p)]
  StZ <- QtZ[-(1:p), , drop = FALSE]

  if (estimator == "VAMP") {
    vr <- VAMP(y = StY, X = StZ, lambda = lambda, psi = psi,
               sigma.sq = sigma.sq, rho = rho, max.iter = max.iter,
               min.rho = min.rho, hs.lb = hs.lb, hs.ub = hs.ub)
    mean <- RtZ %*% vr$mean
    beta.mean <- vr$mean
    covariance   <- vr$variance * tcrossprod(RtZ, RtZ)
    covariance.S <- vr$variance * tcrossprod(StZ, StZ)
    sigma.sq <- vr$sigma.sq
    lambda.HS <- vr$lambda.hs; tau.HS <- vr$tau.hs
  } else {
    vbr <- vblr_fit(StZ, StY, lambda = 1, psi = psi,
                    max_iter = 500, epsilon_conv = 1e-5,
                    shrink.type = lambda, is_verbose = FALSE)
    mean <- RtZ %*% vbr$m; beta.mean <- as.vector(vbr$m)
    covariance   <- RtZ %*% vbr$S %*% t(RtZ)
    covariance.S <- StZ %*% vbr$S %*% t(StZ)
    sigma.sq <- vbr$sigma^2
    lambda.HS <- NA; tau.HS <- NA
  }
  list(RtY = RtY, RtX = RtX, mean = mean, covariance = covariance,
       sigma.sq = sigma.sq, beta.mean = beta.mean,
       cov.S = covariance.S, tau.HS = tau.HS, lambda.HS = lambda.HS)
}

## ----------------------------------------------------------------------------
## Step 3: MCMC sampler for the domestic block beta, conditional on the VAMP
## approximation to the nuisance block.
##
## Returns:
##   beta      : nthinned x p posterior draws of own-country coefficients
##   gamma     : VAMP point estimate of cross-country coefficients
##   bigSigma  : implied T x T error covariance after integrating out alpha
##   lambda.HS.dom, tau.HS.dom : Horseshoe states at the last MCMC iteration
## ----------------------------------------------------------------------------

get.reg.IRGA <- function(y, X, Z, nsave = 1000, nburn = 1000, thin = 1,
                         psi = NULL, lambda = "HS", A.prior = NULL,
                         ols.var = NULL, estimator = "VAMP",
                         hs.lb = 1e-32, hs.ub = 0.1) {
  ntot <- nsave + nburn
  n <- length(y); p <- ncol(X)
  if (is.null(A.prior)) A.prior <- rep(0, p)

  irga12 <- BVS_IRGA_12(y, X, Z, lambda = lambda, psi = psi,
                        sigma.sq = NULL, estimator = estimator,
                        hs.lb = hs.lb, hs.ub = hs.ub)

  Sigma     <- irga12$sigma.sq * diag(p) + irga12$covariance
  Sigma.inv <- solve(Sigma)
  Q         <- t(chol(Sigma.inv))

  Smu.cov  <- irga12$cov.S
  bigSigma <- matrix(0, n, n)
  bigSigma[1:p, 1:p]         <- as.matrix(irga12$covariance)
  if (n > p)
    bigSigma[(p + 1):n, (p + 1):n] <- as.matrix(Smu.cov)
  bigSigma <- bigSigma + diag(n) * irga12$sigma.sq

  ytilde <- Q %*% (irga12$RtY - irga12$mean)
  Xtilde <- Q %*% irga12$RtX
  XSX    <- crossprod(Xtilde)

  if (p > 1) {
    V.prior     <- diag(p) * 1e-1
    V.prior.inv <- diag(1 / diag(V.prior))
  } else {
    V.prior     <- 1
    V.prior.inv <- 1
  }

  gamma <- irga12$beta.mean

  lambda.A <- if (lambda == "MN") as.numeric(ols.var) else 1
  nu.A     <- 1
  tau.A    <- 1
  zeta.A   <- 1

  if (thin < 1) thin <- 1
  nthinned <- floor(nsave / thin)
  beta.store <- matrix(NA, nthinned, p)
  thin.set <- seq(nburn + 1, ntot, by = thin)
  count.thin <- 0

  for (irep in seq_len(ntot)) {
    # Conditional draw of beta -----------------------------------------------
    V.post <- solve(XSX + V.prior.inv)
    A.post <- V.post %*% (crossprod(Xtilde, ytilde) + V.prior.inv %*% A.prior)
    A.draw <- tryCatch(A.post + t(chol(V.post)) %*% rnorm(p),
                       error = function(e) NULL)
    if (is.null(A.draw))
      A.draw <- t(mvtnorm::rmvnorm(1, mean = as.vector(A.post), sigma = V.post))

    # Shrinkage hyperparameters ---------------------------------------------
    if (lambda == "HS") {
      hs_draw <- get.hs(bdraw = (as.vector(A.draw) - A.prior),
                        lambda.hs = lambda.A, nu.hs = nu.A,
                        tau.hs = tau.A, zeta.hs = zeta.A)
      lambda.A <- hs_draw$lambda; nu.A <- hs_draw$nu
      tau.A    <- hs_draw$tau;    zeta.A <- hs_draw$zeta
      hs_draw$psi[hs_draw$psi < 1e-8] <- 1e-8
      V.prior <- diag(hs_draw$psi, nrow = p)
    } else if (lambda == "MN") {
      hs_draw <- get.minn.hs(bdraw = (as.vector(A.draw) - A.prior),
                             lambda.hs = lambda.A,
                             tau.hs = tau.A, zeta.hs = zeta.A)
      lambda.A <- hs_draw$lambda; tau.A <- hs_draw$tau
      zeta.A   <- hs_draw$zeta
      hs_draw$psi[hs_draw$psi < 1e-8] <- 1e-8
      V.prior <- diag(hs_draw$psi, nrow = p)
    }
    V.prior.inv <- if (p > 1) diag(1 / diag(V.prior)) else 1 / V.prior

    if (irep %in% thin.set) {
      count.thin <- count.thin + 1
      beta.store[count.thin, ] <- as.vector(A.draw)
    }
  }

  list(beta = beta.store, bigSigma = bigSigma, gamma = gamma,
       lambda.HS.wex = irga12$lambda.HS, tau.HS.wex = irga12$tau.HS,
       lambda.HS.dom = lambda.A, tau.HS.dom = tau.A)
}

## ----------------------------------------------------------------------------
## PVAR.IRGA: top-level estimator that loops over countries and equations.
##
## Yraw    : T x M data matrix with column names of the form "CC.var"
##           (first two characters = country code).
## p       : VAR lag length.
## nsave/nburn/thin.by : MCMC settings.
## shrink.type : "HS" (Horseshoe) or "MN" (Minnesota).
## estimator   : "VAMP" or "VB".
##
## Returns:
##   A         : nthinned x K x M array of reduced-form VAR coefficients
##                 (K = M*p, ordering = [v1_p1, v2_p1, ..., vM_p1, v1_p2, ...])
##   SIGMA     : nthinned x M x M array of reduced-form innovation covariance
##   Y, lags   : data and lag length used.
##   Upsilon   : structural impact matrix (lower-triangular).
##   VAR.mat   : structural VAR coefficients before back-transformation.
## ----------------------------------------------------------------------------

PVAR.IRGA <- function(Yraw, p, nsave = 1000, nburn = 500, thin.by = 1,
                      shrink.type = "HS", estimator = "VAMP",
                      hs.lb = 1e-32, hs.ub = NULL, quiet = FALSE) {
  X <- mlag(Yraw, p)
  X <- X[(p + 1):nrow(X), , drop = FALSE]
  Y <- Yraw[(p + 1):nrow(Yraw), , drop = FALSE]

  cN <- substr(colnames(Y), 1, 2)
  cN <- cN[!duplicated(cN)]
  colnames(X) <- paste0(rep(colnames(Y), p), "_p", rep(1:p, each = ncol(Y)))

  M <- ncol(Y); K <- ncol(X); T <- nrow(X); N <- length(cN)
  if (thin.by < 1) thin.by <- 1
  nthinned <- floor(nsave / thin.by)
  if (is.null(hs.ub)) hs.ub <- 1 / max(K^2, 1)

  pr.mean <- 0
  Upsilon.big <- diag(M)
  VAR.mat   <- array(0, c(nthinned, K, M))
  Sigma.VAR <- diag(M)
  VAR.store <- vector("list", M)
  count <- 0

  t0 <- Sys.time()
  for (i in seq_len(N)) {
    if (!quiet) message(sprintf("Estimating country %d / %d (%s).",
                                i, N, cN[i]))
    cN.slct <- grepl(cN[i], colnames(Y)) * 1
    G.slct  <- sum(cN.slct)
    sl.dom <- grepl(cN[i], colnames(X))
    sl.in  <- !sl.dom

    X.i <- X[, sl.dom, drop = FALSE]
    Y.i <- Y[, cN.slct == 1, drop = FALSE]
    prior.mean <- matrix(0, ncol(X.i), ncol(Y.i))
    if (pr.mean > 0) prior.mean[1:G.slct, 1:G.slct] <- diag(G.slct) * pr.mean

    for (j in seq_len(G.slct)) {
      count <- count + 1
      if (count > 1)
        Z.i <- cbind(Y[, 1:(count - 1), drop = FALSE],
                     X[, sl.in, drop = FALSE])
      else
        Z.i <- X[, sl.in, drop = FALSE]

      psi.j <- 1 / (0.5 * T * K)
      reg.ij <- get.reg.IRGA(y = Y.i[, j], X = X.i, Z = Z.i,
                             nsave = nsave, nburn = nburn, thin = thin.by,
                             psi = psi.j, lambda = shrink.type,
                             A.prior = prior.mean[, j], estimator = estimator,
                             hs.lb = hs.lb, hs.ub = hs.ub)

      if (count > 1)
        Upsilon.big[count, 1:(count - 1)] <- -reg.ij$gamma[1:(count - 1)]
      idx_int <- (count):length(reg.ij$gamma)
      VAR.mat[, sl.in, count]  <- matrix(rep(reg.ij$gamma[idx_int],
                                             each = nthinned),
                                         nthinned, length(idx_int))
      VAR.mat[, sl.dom, count] <- reg.ij$beta
      VAR.store[[count]]       <- reg.ij
      Sigma.VAR[count, count]  <- reg.ij$bigSigma[T, T]
    }
  }

  Upsilon.inv <- solve(Upsilon.big)
  red.VAR <- apply(VAR.mat, c(1, 2), function(x) Upsilon.inv %*% x)
  A_store <- aperm(red.VAR, c(2, 3, 1))   # nthinned x K x M
  red.VC  <- Upsilon.inv %*% tcrossprod(Sigma.VAR, Upsilon.inv)

  SIGMA_store <- array(NA, dim = c(nthinned, M, M))
  for (irep in seq_len(nthinned)) SIGMA_store[irep, , ] <- red.VC

  if (!quiet) message(sprintf("PVAR-IRGA: %.1f sec.",
                              as.numeric(difftime(Sys.time(), t0,
                                                  units = "secs"))))
  list(A = A_store, SIGMA = SIGMA_store, Y = Y, lags = p,
       countries = cN, Upsilon = Upsilon.big, VAR.mat = VAR.mat,
       store = VAR.store)
}

## ----------------------------------------------------------------------------
## Cholesky impulse responses from posterior VAR coefficients.
## A_post : nthinned x K x M  reduced-form coefficients
## Sigma  : M x M             reduced-form covariance
## p      : VAR lag length
## nstep  : horizon
##
## Returns: nthinned x M x M x nstep  IRF array
## ----------------------------------------------------------------------------

irf_cholesky <- function(A_post, Sigma, p, nstep = 20) {
  ndraw <- dim(A_post)[1]
  M     <- dim(A_post)[3]
  out <- array(0, c(ndraw, M, M, nstep))
  for (d in seq_len(ndraw)) {
    Ad <- A_post[d, , ]                  # K x M
    # Reshape into (M x M x p) "by-equation" coefficient array.
    B <- array(0, c(M, M, p))
    for (lg in seq_len(p))
      B[, , lg] <- t(Ad[((lg - 1) * M + 1):(lg * M), ])  # response var x driver var
    smat <- t(chol(Sigma))
    resp <- array(0, c(M, M, nstep))
    resp[, , 1] <- smat
    for (h in 2:nstep) {
      acc <- matrix(0, M, M)
      for (lg in seq_len(min(p, h - 1)))
        acc <- acc + B[, , lg] %*% resp[, , h - lg]
      resp[, , h] <- acc
    }
    out[d, , , ] <- resp
  }
  out
}
