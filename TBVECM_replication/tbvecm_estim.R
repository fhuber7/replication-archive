# ============================================================================
# Threshold Bayesian Vector Error Correction Model (TBVECM)
#
# Three-regime threshold VECM with regime-switching short-run dynamics,
# a common cointegrating vector beta (linearly normalised), Normal-Gamma
# shrinkage priors on the autoregressive / loading coefficients and a
# Griddy-Gibbs step on the two threshold parameters.
#
# Reduced form for regime s_t in {0, 1, 2}:
#   Delta y_t = (alpha_{s_t})' beta' z_{t-1} + (A_{s_t})' x_t + eps_t
#   eps_t ~ N(0, Sigma_{s_t})
#
# where z_{t-1} are the levels driving the cointegrating relation and x_t
# stacks lagged differences (and a constant if requested). Regime
# assignment is driven by a threshold variable thrsh_t = (z_{t-1}' beta)_g
# crossing the levels (gamma_1, gamma_2).
#
# The estimator returns posterior draws for beta, alpha (per regime),
# A (per regime), Sigma (per regime), gamma_1, gamma_2 and the regime
# state s_t. Companion-form reduced-form coefficients are stored to ease
# the construction of regime-conditional generalised impulse responses.
#
# Required packages: MASS, mnormt (optional), GIGrvg, Matrix, mvtnorm,
# splines, bayesm.
# ============================================================================

requireNamespace("MASS",    quietly = TRUE)
requireNamespace("GIGrvg",  quietly = TRUE)
requireNamespace("Matrix",  quietly = TRUE)
requireNamespace("mvtnorm", quietly = TRUE)
requireNamespace("splines", quietly = TRUE)
requireNamespace("bayesm",  quietly = TRUE)

# ---------------------------------------------------------------------------
# Helper functions (lifted from aux_natural_conjugate.R, trimmed to the
# pieces actually used by the main specification)
# ---------------------------------------------------------------------------

mlag <- function(X, lag) {
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X); N <- ncol(X)
  Xlag <- matrix(0, Traw, p * N)
  for (ii in 1:p) {
    Xlag[(p + 1):Traw, (N * (ii - 1) + 1):(N * ii)] <-
      X[(p + 1 - ii):(Traw - ii), 1:N, drop = FALSE]
  }
  Xlag
}

dmnorm.alt <- function(x, mean = rep(0, d), varcov, log = FALSE) {
  d <- if (is.matrix(varcov)) ncol(varcov) else 1L
  if (d == 1) return(dnorm(x, mean, sqrt(varcov), log = log))
  if (is.vector(x)) x <- t(matrix(x)) else x <- data.matrix(x)
  if (ncol(x) != d) stop("mismatch of dimensions of 'x' and 'varcov'")
  if (is.vector(mean)) {
    mean <- outer(rep(1, nrow(x)), as.vector(matrix(mean, d)))
  }
  X <- t(x - mean)
  conc <- solve(varcov)
  log.det <- as.numeric(determinant(varcov, logarithm = TRUE)$modulus)
  Q <- colSums((conc %*% X) * X)
  logPDF <- as.vector(Q + d * logb(2 * pi) + log.det) / (-2)
  if (log) logPDF else exp(logPDF)
}

# Posterior draw of the stacked (alpha, A) coefficient block for one
# regime under a Normal-Gamma prior with mean A_prior and diagonal
# prior variance diag(Vprior).
get.ALPHA <- function(SIGMAinv, Yt, Xt, Vprior, A_prior) {
  M <- ncol(Yt)
  Vpriorinv <- diag(1 / diag(Vprior))
  IXY      <- kronecker(diag(M), crossprod(Xt, Yt))
  psi_xx   <- kronecker(SIGMAinv, crossprod(Xt))
  V_post   <- try(solve(psi_xx + Vpriorinv), silent = TRUE)
  if (inherits(V_post, "try-error")) V_post <- MASS::ginv(psi_xx + Vpriorinv)
  visig <- as.vector(SIGMAinv)
  a_post <- V_post %*% (IXY %*% visig + Vpriorinv %*% as.vector(A_prior))
  alpha  <- try(a_post + t(chol(V_post)) %*% rnorm(M * ncol(Xt), 0, 1),
                silent = TRUE)
  if (inherits(alpha, "try-error")) alpha <- MASS::mvrnorm(1, a_post, V_post)
  alpha
}

# Normal-Gamma update for the local scaling factors theta (matrix
# version: one element per coefficient).
get.theta <- function(theta, A_draw, A_prior, c_tau, d_tau, a_tau) {
  k <- ncol(A_draw) * nrow(A_draw)
  lambda2_tau <- rgamma(1, c_tau + a_tau * k,
                        d_tau + a_tau / 2 * sum(theta))
  for (mm in 1:ncol(A_draw)) {
    for (nn in 1:nrow(A_draw)) {
      theta[nn, mm] <- GIGrvg::rgig(
        n      = 1,
        lambda = a_tau - 0.5,
        chi    = (A_draw[nn, mm] - A_prior[nn, mm])^2,
        psi    = a_tau * lambda2_tau
      )
    }
  }
  theta[theta < 1e-6] <- 1e-6
  list(theta = theta, lambda2 = lambda2_tau)
}

# Normal-Gamma update for the free elements of the cointegrating vector.
get.theta.beta <- function(theta, A_draw, c_tau, d_tau, a_tau, A_prior = NULL) {
  k <- length(A_draw)
  lambda2_tau <- rgamma(1, c_tau + a_tau * k,
                        d_tau + a_tau / 2 * sum(theta))
  if (is.null(A_prior)) {
    for (mm in 1:k) theta[mm] <- GIGrvg::rgig(
      n = 1, lambda = a_tau - 0.5,
      chi = A_draw[mm]^2, psi = a_tau * lambda2_tau)
  } else {
    for (mm in 1:k) theta[mm] <- GIGrvg::rgig(
      n = 1, lambda = a_tau - 0.5,
      chi = (A_draw[mm] - A_prior[mm])^2, psi = a_tau * lambda2_tau)
  }
  theta[theta < 1e-7] <- 1e-7
  list(theta = theta, lambda2 = lambda2_tau)
}

# Random-walk MH step on the NG concentration parameter.
atau_post <- function(atau, lambda2, thetas, k = length(thetas), rat = 1 / .1) {
  sum(dgamma(thetas, atau, atau * lambda2 / 2, log = TRUE)) +
    dexp(atau, rate = rat, log = TRUE)
}

get.a.theta <- function(a.theta, theta, lambda2, scale1) {
  a_prop      <- exp(rnorm(1, 0, scale1)) * a.theta
  post_a_prop <- atau_post(a_prop, lambda2, theta)
  post_a_old  <- atau_post(a.theta, lambda2, theta)
  post.diff   <- post_a_prop - post_a_old + log(a_prop) - log(a.theta)
  if (is.nan(post.diff)) post.diff <- -Inf
  if (post.diff > log(runif(1))) a_prop else a.theta
}

# Griddy-Gibbs step for one threshold gamma. test1 is a vector of grid
# points; the conditional posterior is evaluated on the grid, smoothed
# with a cubic spline and sampled from the resulting discrete
# approximation.
ggibbs <- function(test1, case, Y, Xt, pr_mean, pr_sd,
                   thrsh, gamma1, gamma2,
                   A1, A2, A3, S1, S2, S3) {
  K <- ncol(Xt); N <- ncol(Y)
  post1 <- numeric(length(test1))
  for (k in seq_along(test1)) {
    g1 <- if (case == 1) test1[[k]] else gamma1
    g2 <- if (case == 2) test1[[k]] else gamma2
    if (g2 <= g1) { post1[k] <- -Inf; next }

    ix1 <- thrsh <  g1
    ix2 <- thrsh >= g1 & thrsh < g2
    ix3 <- thrsh >= g2

    if (sum(ix1) < (N + 2) || sum(ix2) < (N + 2) || sum(ix3) < (N + 2)) {
      post1[k] <- -Inf
      next
    }

    L1 <- sum(dmnorm.alt(Y[ix1, , drop = FALSE],
                         Xt[ix1, , drop = FALSE] %*% matrix(A1, K, N),
                         S1, log = TRUE))
    L2 <- sum(dmnorm.alt(Y[ix2, , drop = FALSE],
                         Xt[ix2, , drop = FALSE] %*% matrix(A2, K, N),
                         S2, log = TRUE))
    L3 <- sum(dmnorm.alt(Y[ix3, , drop = FALSE],
                         Xt[ix3, , drop = FALSE] %*% matrix(A3, K, N),
                         S3, log = TRUE))
    prior_part <- dnorm(test1[[k]], pr_mean, pr_sd, log = TRUE)
    post1[k]   <- L1 + L2 + L3 + prior_part
  }

  finite <- is.finite(post1)
  if (sum(finite) < 4) {
    # Posterior too irregular: keep previous draw.
    return(if (case == 1) gamma1 else gamma2)
  }
  w   <- exp(post1[finite] - max(post1[finite]))
  pts <- test1[finite]
  # Spline-smooth the discrete posterior, then sample from the
  # interpolated approximation.
  sp <- try(splines::interpSpline(w ~ seq_along(w), bSpline = TRUE),
            silent = TRUE)
  if (inherits(sp, "try-error")) {
    probs <- w / sum(w)
    return(sample(pts, 1, prob = probs))
  }
  smooth      <- predict(sp)$y
  smooth[smooth < 0] <- 0
  if (sum(smooth) <= 0) {
    probs <- w / sum(w)
    return(sample(pts, 1, prob = probs))
  }
  probs       <- smooth / sum(smooth)
  grid_interp <- approx(pts, n = length(smooth))$y
  sample(grid_interp, 1, prob = probs)
}

# Reduced-form impulse responses via direct recursion (Sims 2003 style).
impulsdtrf <- function(B, smat, nstep) {
  neq  <- dim(B)[1]; nvar <- dim(B)[2]; lags <- dim(B)[3]
  if (dim(smat)[2] != dim(B)[2]) stop("B and smat conflict on # of variables")
  response <- array(0, c(neq, nvar, nstep + lags - 1))
  response[, , lags] <- smat
  response <- aperm(response, c(1, 3, 2))
  irhs <- 1:(lags * nvar)
  ilhs <- lags * nvar + (1:nvar)
  response <- matrix(response, ncol = neq)
  B <- B[, , seq(from = lags, to = 1, by = -1), drop = FALSE]
  B <- matrix(B, nrow = nvar)
  for (it in 1:(nstep - 1)) {
    response[ilhs, ] <- B %*% response[irhs, ]
    irhs <- irhs + nvar; ilhs <- ilhs + nvar
  }
  dim(response) <- c(nvar, nstep + lags - 1, nvar)
  aperm(response[, -(1:(lags - 1)), , drop = FALSE], c(1, 3, 2))
}

# ---------------------------------------------------------------------------
# Main estimator
# ---------------------------------------------------------------------------
# Yraw          : (T0 x M) matrix of *levels* of the endogenous variables
# r             : cointegrating rank
# p             : number of lags of *differences* (so the underlying VAR
#                 has p + 1 lags in levels)
# nsave         : posterior draws to keep
# nburn         : burn-in draws
# nhor          : IRF horizon
# beta_prior_var: tightness of the Normal prior on the free elements of
#                 the cointegrating vector
# sample.theta  : whether to update the NG concentration via MH
# verbose_every : print progress every k iterations (0 to silence)
#
# beta is normalised by stacking an identity in the first r rows
# (linear normalisation), so the free elements are beta[(r+1):M, ].
# ---------------------------------------------------------------------------

tbvecm_estim <- function(Yraw,
                         r              = 1,
                         p              = 1,
                         nsave          = 3000,
                         nburn          = 2000,
                         nhor           = 24,
                         beta_prior_var = 0.01,
                         coint.shrink   = FALSE,
                         sample.theta   = TRUE,
                         cons           = TRUE,
                         a.theta.init   = 0.1,
                         verbose_every  = 200) {

  Yraw <- as.matrix(Yraw)
  if (is.null(colnames(Yraw))) colnames(Yraw) <- paste0("y", 1:ncol(Yraw))

  # Demean the levels (the constant in the cointegration relation is
  # absorbed by the demeaning, in line with the empirical application).
  Yraw <- apply(Yraw, 2, function(x) x - mean(x))

  # Build regressors -------------------------------------------------------
  # Z: levels lagged once (loads the cointegrating relation)
  Z <- mlag(Yraw, 1)
  Z <- Z[(p + 2):nrow(Z), , drop = FALSE]
  Z <- Z[1:(nrow(Z) - 1), , drop = FALSE]  # h = 1 forecast lag

  Y <- diff(Yraw)
  Xlag <- mlag(Y, p)
  if (cons) {
    X <- cbind(Xlag[(p + 1):(nrow(Xlag) - 1), , drop = FALSE], 1)
  } else {
    X <- Xlag[(p + 1):(nrow(Xlag) - 1), , drop = FALSE]
  }
  Y <- Y[(p + 1):(nrow(Y) - 1), , drop = FALSE]

  M      <- ncol(Y); K <- ncol(X); T <- nrow(Y)
  M.beta <- M

  # Linear normalisation of beta: first r rows = I_r ----------------------
  beta.draw <- matrix(0, M.beta, r)
  beta.draw[1:r, 1:r] <- diag(r)

  # Starting values for alpha, A, Sigma -----------------------------------
  D       <- Z %*% beta.draw
  Xtilde  <- cbind(D, X)
  A.start <- solve(crossprod(Xtilde) +
                     diag(1e-4, ncol(Xtilde))) %*% crossprod(Xtilde, Y)
  alpha.1 <- alpha.2 <- alpha.3 <- A.start[1:r, , drop = FALSE]
  A.draw.1 <- A.draw.2 <- A.draw.3 <- A.start[(r + 1):(K + r), , drop = FALSE]
  Sig0 <- crossprod(Y - Z %*% t(t(alpha.1) %*% t(beta.draw)) -
                      X %*% A.draw.1) / (T - K)
  Sigma.draw.1 <- Sigma.draw.2 <- Sigma.draw.3 <- Sig0
  Sigma.inv.1  <- Sigma.inv.2  <- Sigma.inv.3  <- solve(Sig0)

  # Univariate AR(1) sigmas (used for the adaptive IW scale) --------------
  sigma2.scale <- matrix(0, M, 1)
  for (nn in 1:M) {
    fit <- lm(Y[-1, nn] ~ Y[-nrow(Y), nn])
    sigma2.scale[nn, ] <- summary(fit)$sigma
  }

  # Prior on beta: identity in the top r rows, Normal-Gamma on the free
  # part (M - r) x r block.
  P.mean <- beta.draw * 1; diag(P.mean) <- 1
  P.mat.M <- matrix(beta_prior_var, M.beta, r)
  P.mat.M[1:r, 1:r][upper.tri(P.mat.M[1:r, 1:r])] <- 0
  P.mat.M[1:r, 1:r][lower.tri(P.mat.M[1:r, 1:r])] <- 0
  diag(P.mat.M) <- 1e-12

  # Selection matrices that pin the identity in the top r rows and
  # leave (M - r) x r free elements --------------------------------------
  H.i <- t(cbind(matrix(0, M.beta - r, r), diag(M.beta - r)))
  bH  <- vector("list", r); for (j in 1:r) bH[[j]] <- H.i
  H   <- as.matrix(Matrix::bdiag(bH))
  s   <- as.vector(diag(M.beta)[, 1:r, drop = FALSE])

  # Threshold initialisation: use first differenced series as the
  # threshold-driver until beta has been refined; then switch to
  # z' beta (component sl.gamma). This matches the original code.
  sl.gamma <- 1
  thrsh    <- Y[, 1]
  gamma1   <- quantile(thrsh, 0.4)
  gamma2   <- quantile(thrsh, 0.7)

  # Prior on (alpha, A) ---------------------------------------------------
  A.prior <- matrix(0, K + r, M)
  pr.mean <- 0.9
  A.prior[(r + 1):(r + M), ] <- diag(M) * pr.mean

  st <- integer(T)
  st[thrsh <  gamma1]                    <- 0
  st[thrsh >= gamma1 & thrsh < gamma2]   <- 1
  st[thrsh >= gamma2]                    <- 2

  # NG prior bookkeeping --------------------------------------------------
  d.1 <- 0.001; d.2 <- 0.001
  theta1 <- theta2 <- theta3 <- matrix(0.6, K + r, M)
  a.A.1 <- a.A.2 <- a.A.3 <-
    a.theta.1 <- a.theta.2 <- a.theta.3 <- a.beta <- a.theta.init
  scaleRWMH <- 1

  # Storage ---------------------------------------------------------------
  beta.store  <- array(NA, c(nsave, M.beta, r))
  A.store     <- array(NA, c(nsave, K, M, 3))
  alpha.store <- array(NA, c(nsave, r, M, 3))
  SIGMA.store <- array(NA, c(nsave, M, M, 3))
  PHI.store   <- array(NA, c(nsave, (p + 1) * M, M, 3))
  IRF.store   <- array(NA, c(nsave, M, M, nhor, 3))
  st.store    <- matrix(NA, nsave, T)
  gamma.store <- matrix(NA, nsave, 2)
  thrsh.store <- matrix(NA, nsave, T)

  ntot <- nsave + nburn
  if (verbose_every > 0) message("Starting MCMC: ", ntot, " iterations")

  for (irep in 1:ntot) {
    # --- Step 0: split data by current regime ----------------------------
    D.draw <- Z %*% beta.draw
    idx <- list(which(st == 0), which(st == 1), which(st == 2))
    if (any(sapply(idx, length) < (K + r + M))) {
      # Regime collapsed: nudge the thresholds to safe quantiles.
      gamma1 <- as.numeric(quantile(thrsh, 0.4))
      gamma2 <- as.numeric(quantile(thrsh, 0.7))
      st <- integer(T)
      st[thrsh <  gamma1]                  <- 0
      st[thrsh >= gamma1 & thrsh < gamma2] <- 1
      st[thrsh >= gamma2]                  <- 2
      idx <- list(which(st == 0), which(st == 1), which(st == 2))
    }

    pull <- function(M_, ix) M_[ix, , drop = FALSE]
    Z.1 <- pull(Z, idx[[1]]); Z.2 <- pull(Z, idx[[2]]); Z.3 <- pull(Z, idx[[3]])
    X.1 <- pull(X, idx[[1]]); X.2 <- pull(X, idx[[2]]); X.3 <- pull(X, idx[[3]])
    Y.1 <- pull(Y, idx[[1]]); Y.2 <- pull(Y, idx[[2]]); Y.3 <- pull(Y, idx[[3]])
    D.1 <- pull(D.draw, idx[[1]]); D.2 <- pull(D.draw, idx[[2]])
    D.3 <- pull(D.draw, idx[[3]])

    Xt.1 <- cbind(D.1, X.1); Xt.2 <- cbind(D.2, X.2); Xt.3 <- cbind(D.3, X.3)

    # --- Step I: regime-specific (alpha, A) ------------------------------
    A.full.1 <- matrix(get.ALPHA(Sigma.inv.1, Y.1, Xt.1,
                                 diag(as.vector(theta1)), A.prior),
                       K + r, M)
    A.full.2 <- matrix(get.ALPHA(Sigma.inv.2, Y.2, Xt.2,
                                 diag(as.vector(theta2)), A.prior),
                       K + r, M)
    A.full.3 <- matrix(get.ALPHA(Sigma.inv.3, Y.3, Xt.3,
                                 diag(as.vector(theta3)), A.prior),
                       K + r, M)

    alpha.1 <- A.full.1[1:r, , drop = FALSE]
    alpha.2 <- A.full.2[1:r, , drop = FALSE]
    alpha.3 <- A.full.3[1:r, , drop = FALSE]
    A.draw.1 <- A.full.1[(r + 1):(K + r), , drop = FALSE]
    A.draw.2 <- A.full.2[(r + 1):(K + r), , drop = FALSE]
    A.draw.3 <- A.full.3[(r + 1):(K + r), , drop = FALSE]

    # --- Step Ib: NG hyperparameters on (alpha, A) -----------------------
    t1A <- get.theta(theta1[(r + 1):(K + r), , drop = FALSE], A.draw.1,
                     A.prior[(r + 1):(K + r), , drop = FALSE],
                     d.1, d.2, a.A.1)
    t2A <- get.theta(theta2[(r + 1):(K + r), , drop = FALSE], A.draw.2,
                     A.prior[(r + 1):(K + r), , drop = FALSE],
                     d.1, d.2, a.A.2)
    t3A <- get.theta(theta3[(r + 1):(K + r), , drop = FALSE], A.draw.3,
                     A.prior[(r + 1):(K + r), , drop = FALSE],
                     d.1, d.2, a.A.3)
    t1a <- get.theta(theta1[1:r, , drop = FALSE], alpha.1,
                     A.prior[1:r, , drop = FALSE], d.1, d.2, a.theta.1)
    t2a <- get.theta(theta2[1:r, , drop = FALSE], alpha.2,
                     A.prior[1:r, , drop = FALSE], d.1, d.2, a.theta.2)
    t3a <- get.theta(theta3[1:r, , drop = FALSE], alpha.3,
                     A.prior[1:r, , drop = FALSE], d.1, d.2, a.theta.3)
    theta1 <- rbind(t1a$theta, t1A$theta)
    theta2 <- rbind(t2a$theta, t2A$theta)
    theta3 <- rbind(t3a$theta, t3A$theta)

    # --- Step II: cointegrating vector beta ------------------------------
    Yhat.1 <- as.vector(Y.1 - X.1 %*% A.draw.1) -
              kronecker(t(alpha.1), Z.1) %*% s
    Yhat.2 <- as.vector(Y.2 - X.2 %*% A.draw.2) -
              kronecker(t(alpha.2), Z.2) %*% s
    Yhat.3 <- as.vector(Y.3 - X.3 %*% A.draw.3) -
              kronecker(t(alpha.3), Z.3) %*% s

    Bm.1 <- kronecker(alpha.1 %*% Sigma.inv.1, t(Z.1)) %*% Yhat.1
    Bm.2 <- kronecker(alpha.2 %*% Sigma.inv.2, t(Z.2)) %*% Yhat.2
    Bm.3 <- kronecker(alpha.3 %*% Sigma.inv.3, t(Z.3)) %*% Yhat.3
    B.mean <- crossprod(H, Bm.1 + Bm.2 + Bm.3)

    Bv.1 <- kronecker(alpha.1 %*% Sigma.inv.1 %*% t(alpha.1), crossprod(Z.1))
    Bv.2 <- kronecker(alpha.2 %*% Sigma.inv.2 %*% t(alpha.2), crossprod(Z.2))
    Bv.3 <- kronecker(alpha.3 %*% Sigma.inv.3 %*% t(alpha.3), crossprod(Z.3))
    B.var <- t(H) %*% (Bv.1 + Bv.2 + Bv.3) %*% H

    P.free   <- diag(as.vector(P.mat.M[(r + 1):M.beta, , drop = FALSE]),
                     nrow = (M.beta - r) * r)
    P.inv    <- diag(1 / diag(P.free), nrow = nrow(P.free))
    xi.var   <- solve(B.var + P.inv)
    xi.mean  <- xi.var %*% (B.mean + P.inv %*%
                  as.vector(P.mean[(r + 1):M.beta, , drop = FALSE]))
    xi.draw  <- MASS::mvrnorm(1, as.vector(xi.mean), xi.var)
    beta.draw <- matrix(s + H %*% xi.draw, M.beta, r)

    # NG update on the free beta elements
    if (r > 1 && !coint.shrink) {
      pmm  <- P.mat.M[(r + 1):M.beta, 2:r, drop = FALSE]
      bdm  <- beta.draw[(r + 1):M.beta, 2:r, drop = FALSE]
      Pv   <- get.theta.beta(as.vector(pmm), as.vector(bdm),
                             d.1, d.2, a.beta)
      P.mat.M[(r + 1):M.beta, 2:r] <- matrix(Pv$theta, M - r, r - 1)
    }
    if (coint.shrink) {
      pmm  <- P.mat.M[(r + 1):M.beta, , drop = FALSE]
      bdm  <- beta.draw[(r + 1):M.beta, , drop = FALSE]
      Pv   <- get.theta.beta(as.vector(pmm), as.vector(bdm),
                             d.1, d.2, a.beta)
      P.mat.M[(r + 1):M.beta, ] <- matrix(Pv$theta, M - r, r)
    } else if (r == 1) {
      Pv <- list(theta = as.vector(P.mat.M[(r + 1):M.beta, ]), lambda2 = 1)
    }

    # MH on NG concentration parameters
    if (sample.theta) {
      a.A.1     <- get.a.theta(a.A.1,     t1A$theta, t1A$lambda2, scaleRWMH)
      a.A.2     <- get.a.theta(a.A.2,     t2A$theta, t2A$lambda2, scaleRWMH)
      a.A.3     <- get.a.theta(a.A.3,     t3A$theta, t3A$lambda2, scaleRWMH)
      a.theta.1 <- get.a.theta(a.theta.1, t1a$theta, t1a$lambda2, scaleRWMH)
      a.theta.2 <- get.a.theta(a.theta.2, t2a$theta, t2a$lambda2, scaleRWMH)
      a.theta.3 <- get.a.theta(a.theta.3, t3a$theta, t3a$lambda2, scaleRWMH)
      if (r > 1) a.beta <- get.a.theta(a.beta, Pv$theta, Pv$lambda2, scaleRWMH)
    }

    # --- Step III: regime-specific Sigma (adaptive IW prior) -------------
    err.1 <- Y.1 - Z.1 %*% t(t(alpha.1) %*% t(beta.draw)) - X.1 %*% A.draw.1
    err.2 <- Y.2 - Z.2 %*% t(t(alpha.2) %*% t(beta.draw)) - X.2 %*% A.draw.2
    err.3 <- Y.3 - Z.3 %*% t(t(alpha.3) %*% t(beta.draw)) - X.3 %*% A.draw.3

    c0 <- 2.5 + (M - 1) / 2
    g0 <- 0.5 + (M - 1) / 2
    S0 <- 100 * g0 / c0 * diag(as.numeric(sigma2.scale))
    C0 <- bayesm::rwishart(
      2 * (g0 + c0 * 3),
      0.5 * chol2inv(chol(S0 + Sigma.inv.1 + Sigma.inv.2 + Sigma.inv.3))
    )$W
    scale1 <- crossprod(err.1) / 2 + C0
    scale2 <- crossprod(err.2) / 2 + C0
    scale3 <- crossprod(err.3) / 2 + C0
    v1 <- nrow(err.1) / 2 + c0
    v2 <- nrow(err.2) / 2 + c0
    v3 <- nrow(err.3) / 2 + c0
    Sigma.draw.1 <- bayesm::rwishart(2 * v1, 0.5 * chol2inv(chol(scale1)))$IW
    Sigma.draw.2 <- bayesm::rwishart(2 * v2, 0.5 * chol2inv(chol(scale2)))$IW
    Sigma.draw.3 <- bayesm::rwishart(2 * v3, 0.5 * chol2inv(chol(scale3)))$IW
    Sigma.inv.1 <- solve(Sigma.draw.1)
    Sigma.inv.2 <- solve(Sigma.draw.2)
    Sigma.inv.3 <- solve(Sigma.draw.3)

    # --- Step IV: thresholds via Griddy-Gibbs ----------------------------
    thrsh <- as.vector((Z %*% beta.draw)[, sl.gamma, drop = FALSE])
    thrsh <- thrsh - mean(thrsh)
    X.beta <- cbind(Z %*% beta.draw, X)

    grid.lo <- approx(c(quantile(thrsh, 0.01), quantile(thrsh, 0.645)),
                      n = 20)$y
    gamma1  <- ggibbs(grid.lo, 1, Y, X.beta, 0, 10, thrsh, gamma1, gamma2,
                      A.full.1, A.full.2, A.full.3,
                      Sigma.draw.1, Sigma.draw.2, Sigma.draw.3)
    grid.hi <- approx(c(quantile(thrsh, 0.65), quantile(thrsh, 0.99)),
                      n = 20)$y
    gamma2  <- ggibbs(grid.hi, 2, Y, X.beta, 0, 10, thrsh, gamma1, gamma2,
                      A.full.1, A.full.2, A.full.3,
                      Sigma.draw.1, Sigma.draw.2, Sigma.draw.3)

    st <- integer(T)
    st[thrsh <  gamma1]                  <- 0
    st[thrsh >= gamma1 & thrsh < gamma2] <- 1
    st[thrsh >= gamma2]                  <- 2

    # --- Store -----------------------------------------------------------
    if (irep > nburn) {
      k <- irep - nburn
      beta.store[k, , ]  <- beta.draw
      thrsh.store[k, ]   <- thrsh
      st.store[k, ]      <- st
      gamma.store[k, ]   <- c(gamma1, gamma2)
      A.store[k, , , 1]  <- A.draw.1
      A.store[k, , , 2]  <- A.draw.2
      A.store[k, , , 3]  <- A.draw.3
      alpha.store[k, , , 1] <- alpha.1
      alpha.store[k, , , 2] <- alpha.2
      alpha.store[k, , , 3] <- alpha.3
      SIGMA.store[k, , , 1] <- Sigma.draw.1
      SIGMA.store[k, , , 2] <- Sigma.draw.2
      SIGMA.store[k, , , 3] <- Sigma.draw.3

      # Companion-form coefficients (PHI) and Cholesky IRFs per regime.
      # A_r is (K x M) with K = M*p + cons; the last row (if cons) is
      # the intercept, which does not enter the impulse-response
      # recursion. We strip it before building the level-form VAR.
      strip_cons <- function(A_r) {
        if (cons) A_r[1:(K - 1), , drop = FALSE] else A_r
      }
      build_phi <- function(alpha_r, A_r) {
        A_lev <- t(strip_cons(A_r))                            # M x (M*p)
        Pi    <- t(alpha_r) %*% t(beta.draw)                   # M x M
        # First lag: identity + Pi + first AR block - first AR block-lag
        if (p == 1) {
          first  <- diag(M) + Pi + A_lev[, 1:M, drop = FALSE]
          others <- NULL
        } else {
          first  <- diag(M) + Pi + A_lev[, 1:M, drop = FALSE]
          others <- matrix(0, M, M * (p - 1))
          for (ll in 2:p) {
            others[, ((ll - 2) * M + 1):((ll - 1) * M)] <-
              A_lev[, ((ll - 1) * M + 1):(ll * M), drop = FALSE] -
              A_lev[, ((ll - 2) * M + 1):((ll - 1) * M), drop = FALSE]
          }
        }
        last <- -A_lev[, ((p - 1) * M + 1):(p * M), drop = FALSE]
        cbind(first, others, last)                             # M x (M*(p+1))
      }
      PHI.mat.1 <- build_phi(alpha.1, A.draw.1)
      PHI.mat.2 <- build_phi(alpha.2, A.draw.2)
      PHI.mat.3 <- build_phi(alpha.3, A.draw.3)
      PHI.store[k, , , 1] <- t(PHI.mat.1)
      PHI.store[k, , , 2] <- t(PHI.mat.2)
      PHI.store[k, , , 3] <- t(PHI.mat.3)

      reg_arr <- function(alpha_r, A_r) {
        ph <- build_phi(alpha_r, A_r)
        a  <- array(0, c(M, M, p + 1))
        for (ll in 1:(p + 1)) {
          a[, , ll] <- ph[, ((ll - 1) * M + 1):(ll * M), drop = FALSE]
        }
        a
      }
      Sig_list <- list(Sigma.draw.1, Sigma.draw.2, Sigma.draw.3)
      phi_list <- list(reg_arr(alpha.1, A.draw.1),
                       reg_arr(alpha.2, A.draw.2),
                       reg_arr(alpha.3, A.draw.3))
      for (rg in 1:3) {
        Sg <- Sig_list[[rg]]
        ch <- try(t(chol(Sg)), silent = TRUE)
        if (inherits(ch, "try-error")) next
        irf <- impulsdtrf(phi_list[[rg]], ch, nhor)
        IRF.store[k, , , , rg] <- irf
      }
    }

    if (verbose_every > 0 && irep %% verbose_every == 0) {
      message("  iter ", irep, " / ", ntot,
              "  gamma=(", sprintf("%.3f", gamma1), ", ",
              sprintf("%.3f", gamma2), ")")
    }
  }

  list(
    beta   = beta.store,
    alpha  = alpha.store,
    A      = A.store,
    Sigma  = SIGMA.store,
    PHI    = PHI.store,
    IRF    = IRF.store,
    st     = st.store,
    gamma  = gamma.store,
    thrsh  = thrsh.store,
    dims   = list(M = M, K = K, r = r, p = p, T = T, nhor = nhor,
                  varnames = colnames(Yraw))
  )
}
