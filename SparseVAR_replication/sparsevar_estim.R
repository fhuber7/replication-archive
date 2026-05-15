# =============================================================================
# Sparser TVP-VARs - self-contained estimator
# Hauzenberger, Huber, Onorante (JAE, 2020):
#   "Combining Shrinkage and Sparsity in Conjugate Vector Autoregressive Models"
#
# Headline model:
#   * Conjugate VAR(p) with Minnesota dummy-observation prior (Banbura et al.)
#   * Optional common stochastic volatility (pca.sv, via stochvol::svsample2)
#   * After each Gibbs draw of (Alpha, Sigma), coefficients are sparsified
#     via lag-wise SAVS thresholding (decoupled shrinkage and selection, DSS)
#   * Cholesky factor of Sigma can also be sparsified
#
# Stripped of: competitors (NG-SSVS, OLS, RW), hyperparameter MH, forecast
# evaluation against held-out data, sign/external-instrument identification.
#
# Dependencies: MASS, stochvol, mvtnorm
# =============================================================================

suppressPackageStartupMessages({
  library(MASS)
  library(stochvol)
  library(mvtnorm)
})

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

mlag <- function(X, lag) {
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X); N <- ncol(X)
  Xlag <- matrix(0, Traw, p * N)
  for (ii in 1:p) {
    Xlag[(p + 1):Traw, (N * (ii - 1) + 1):(N * ii)] <-
      X[(p + 1 - ii):(Traw - ii), 1:N]
  }
  Xlag
}

norm_vec <- function(x) sqrt(sum(x^2))

# SAVS sparsification of a single regression coefficient vector b given X
# kappa controls penalty exponent, lambda is the scale.
sparsify_chol <- function(X, b, kappa) {
  K <- nrow(b)
  bs <- matrix(0, K, 1)
  for (jj in seq_len(K)) {
    mu.jj <- 1 / abs(b[jj, 1])^kappa
    nj <- norm_vec(X[, jj])^2
    if ((abs(b[jj, 1]) * nj) < mu.jj) {
      bs[jj, 1] <- 0
    } else {
      bs[jj, 1] <- sign(b[jj, 1]) * (abs(b[jj, 1]) * nj - mu.jj) / nj
    }
  }
  bs
}

# Companion form
get.companion <- function(Beta_, varndxv) {
  nn <- varndxv[[1]]; nd <- varndxv[[2]]; nl <- varndxv[[3]]
  nkk <- nn * nl + nd
  Jm <- matrix(0, nkk, nn); Jm[1:nn, 1:nn] <- diag(nn)
  if (nd <= 1) {
    MM <- rbind(t(Beta_),
                cbind(diag((nl - 1) * nn), matrix(0, (nl - 1) * nn, nn + 1)),
                c(matrix(0, 1, nn * nl), 1))
  } else {
    MM <- rbind(t(Beta_),
                cbind(diag((nl - 1) * nn), matrix(0, (nl - 1) * nn, nn + nd)),
                cbind(matrix(0, nd, nn * nl), diag(1, nd)))
  }
  list(MM = MM, Jm = Jm)
}

# Cholesky/recursive IRFs (Sims-style)
impulsdtrf <- function(B, smat, nstep) {
  neq <- dim(B)[1]; nvar <- dim(B)[2]; lags <- dim(B)[3]
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
  response <- aperm(response[, -(1:(lags - 1)), ], c(1, 3, 2))
  response
}

# Minnesota dummies (Banbura/Sims/Zha style)
get.dum.min <- function(theta, gamma.prior, delta, mysigma, M, p) {
  ydummy <- matrix(0, 2 * M + M * (p - 1) + 1, M)
  xdummy <- matrix(0, 2 * M + M * (p - 1) + 1, M * p + 1)
  ydummy[1:M, ] <- diag((as.numeric(mysigma) * delta) / theta)
  ydummy[(M * (p - 1) + M + 1):(M * (p - 1) + 2 * M), ] <- diag(as.numeric(mysigma))
  jp <- diag(1:p)
  xdummy[1:(M * p), 1:(M * p)] <- kronecker(jp, diag(as.numeric(mysigma))) / theta
  xdummy[nrow(xdummy), ncol(xdummy)] <- gamma.prior
  list(Ydum = ydummy, Xdum = xdummy)
}

# -----------------------------------------------------------------------------
# Main estimator
# -----------------------------------------------------------------------------
# Yraw         : T x M data matrix
# p            : lag order
# nsave, nburn : Gibbs iterations
# cons         : include intercept (TRUE/FALSE)
# prmean       : vector length M, prior mean of own first lag (1 = RW, 0 = WN)
# shrink.1     : Minnesota overall shrinkage hyperparameter
# kappa.base   : SAVS penalty exponent (>1 with lag-wise spec)
# lambda.1st   : SAVS scale at lag 1
# lambda.2nd   : SAVS scale at lag 2 (decays geometrically afterwards)
# lambda.decay : geometric decay of lambda across lags >= 2
# kappa.decay  : exponent decay across lags
# kappa.cov    : SAVS exponent for the Cholesky-factor off-diagonals
# sigma.sparse : if TRUE, sparsify Sigma via its Cholesky factor
# sv           : "none" or "pca.sv"
# nhor         : IRF horizon
# fhorz        : forecast horizon
# -----------------------------------------------------------------------------

sparsevar_estim <- function(Yraw,
                            p            = 4,
                            nsave        = 3000,
                            nburn        = 2000,
                            cons         = TRUE,
                            prmean       = NULL,
                            shrink.1     = 0.2,
                            kappa.base   = 2,
                            lambda.1st   = 1,
                            lambda.2nd   = 1,
                            lambda.decay = 0.5,
                            kappa.decay  = 1,
                            kappa.cov    = 4,
                            sigma.sparse = TRUE,
                            sv           = c("pca.sv", "none"),
                            nhor         = 24,
                            fhorz        = 8,
                            verbose      = TRUE) {

  sv <- match.arg(sv)
  ntot <- nsave + nburn

  # --------------------------- DATA ------------------------------------------
  Yraw <- as.matrix(Yraw)
  if (is.null(colnames(Yraw))) colnames(Yraw) <- paste0("y", seq_len(ncol(Yraw)))
  var.names <- colnames(Yraw)

  if (cons) X <- cbind(mlag(Yraw, p), 1) else X <- mlag(Yraw, p)
  Y <- Yraw[(p + 1):nrow(Yraw), , drop = FALSE]
  X <- X[(p + 1):nrow(X), , drop = FALSE]

  T <- nrow(X); K <- ncol(X); M <- ncol(Y); n <- K * M
  if (is.null(prmean)) prmean <- rep(0, M)

  # --------------------------- SAVS lag-wise penalties -----------------------
  if (p > 1) {
    kappa.p <- kappa.base
    lambda.p <- c(lambda.1st, lambda.2nd)
    for (jj in 2:p) {
      kappa.p[jj] <- (kappa.p[jj - 1] * kappa.base)^kappa.decay
      if (jj < p) lambda.p[jj + 1] <- lambda.p[jj] * lambda.decay
    }
    kappa.vec <- rep(kappa.p, each = M)
    lambda.vec <- rep(lambda.p, each = M)
    if (cons) {
      kappa.vec  <- matrix(c(kappa.vec,  kappa.vec[1]), K, 1)
      lambda.vec <- matrix(c(lambda.vec, 1e-10),         K, 1)
    } else {
      kappa.vec  <- matrix(kappa.vec,  K, 1)
      lambda.vec <- matrix(lambda.vec, K, 1)
    }
  } else {
    kappa.vec  <- matrix(kappa.base, K, 1)
    lambda.vec <- matrix(lambda.1st, K, 1)
  }

  # --------------------------- OLS init --------------------------------------
  XtXinv <- MASS::ginv(crossprod(X))
  A_OLS  <- XtXinv %*% crossprod(X, Y)
  Em     <- Y - X %*% A_OLS
  SSE    <- crossprod(Em)
  Sigma.draw <- if (K < T) SSE / (T - K + 1) else SSE / (T + 1)

  # --------------------------- Prior (dummies) -------------------------------
  gamma.prior <- 0.01
  n0 <- M + 2

  # AR(p) sigmas for Minnesota scaling
  mysigma <- matrix(NA, M, 1)
  for (ii in 1:M) {
    Y.i <- Yraw[, ii, drop = FALSE]
    X.i <- mlag(Yraw[, ii], p)
    Y.i <- Y.i[(p + 1):nrow(Y.i), , drop = FALSE]
    X.i <- X.i[(p + 1):nrow(X.i), , drop = FALSE]
    rho.i <- ginv(crossprod(X.i)) %*% crossprod(X.i, Y.i)
    er.i  <- Y.i - X.i %*% rho.i
    mysigma[ii, 1] <- sqrt(crossprod(er.i) / (nrow(X.i) - p))
  }

  dummies.min <- get.dum.min(shrink.1, gamma.prior, prmean, mysigma, M, p)
  Y.min <- dummies.min$Ydum; X.min <- dummies.min$Xdum

  XX <- rbind(X, X.min); YY <- rbind(Y, Y.min)

  # Posterior moments
  V.post <- try(solve(crossprod(XX)), silent = TRUE)
  if (inherits(V.post, "try-error")) V.post <- ginv(crossprod(XX))
  V.chol <- t(chol(V.post))
  A.post <- V.post %*% crossprod(XX, YY)
  Scale.post <- crossprod(YY - XX %*% A.post)
  Sigma.draw <- Scale.post / (nrow(YY) - K - 1)

  # --------------------------- SV initialisation -----------------------------
  if (sv == "pca.sv") {
    h_0 <- 0; V_h0 <- 10                          # prior on mu (mean log-vol)
    sv_priors <- stochvol::specify_priors(
      mu    = stochvol::sv_normal(mean = h_0, sd = sqrt(V_h0)),
      phi   = stochvol::sv_beta(shape1 = 25, shape2 = 1.5),
      sigma2 = stochvol::sv_gamma(shape = 0.5, rate = 0.5)
    )
    sv_para   <- list(mu = -10, phi = 0.9, sigma = 0.2,
                      nu = Inf, rho = 0, beta = NA, latent0 = 0)
    sv_latent <- rep(0, T)
    normalizer <- 1 / exp(sv_latent / 2)
  } else {
    normalizer <- rep(1, T)
  }

  # --------------------------- Storage ---------------------------------------
  ALPHA.store <- array(0, c(nsave, K, M))   # un-sparsified VAR coefs
  A.store     <- array(0, c(nsave, K, M))   # sparsified VAR coefs
  SIGMA.store <- array(0, c(nsave, M, M))
  S.store     <- array(0, c(nsave, M, M))
  eht.store   <- matrix(0, nsave, T)
  delta.store <- matrix(0, nsave, n)        # inclusion indicator (vec form)
  irf.store        <- array(NA, c(nsave, M, M, nhor))  # IRFs (unsparsified)
  irf.sparse.store <- array(NA, c(nsave, M, M, nhor))  # IRFs (sparsified)
  pred.store        <- array(NA, c(nsave, M, fhorz))
  pred.sparse.store <- array(NA, c(nsave, M, fhorz))

  # --------------------------- Gibbs sampler ---------------------------------
  for (irep in seq_len(ntot)) {

    # ---- Step 1: Sample Alpha from matric-variate Normal -------------------
    Alpha.draw <- matrix(as.vector(A.post) +
                           as.vector(V.chol %*% matrix(rnorm(n), K, M) %*%
                                       chol(Sigma.draw)),
                         K, M)

    # ---- Step 2: SAVS sparsification of Alpha (lag-wise) -------------------
    A.draw <- matrix(0, K, M)
    norm.i <- apply(X * normalizer, 2, norm_vec)^2
    for (ii in seq_len(M)) {
      mu.i <- exp(log(lambda.vec) - kappa.vec * log(abs(Alpha.draw[, ii])))
      ind.change <- (abs(Alpha.draw[, ii]) * norm.i) > mu.i
      A.draw[ind.change, ii] <-
        (sign(Alpha.draw[, ii]) * (1 / norm.i) *
           (abs(Alpha.draw[, ii]) * norm.i - mu.i))[ind.change]
      A.draw[!ind.change, ii] <- 0
    }
    if (cons) A.draw[K, ] <- Alpha.draw[K, ]   # never zero out intercept

    delta.draw <- as.numeric(as.vector(A.draw) != 0)

    # ---- Step 3: Sample Sigma ---------------------------------------------
    Sigma.draw.inv <- matrix(rWishart(1, n0 + T, solve(Scale.post)), M, M)
    Sigma.draw <- solve(Sigma.draw.inv)

    # ---- Step 4: Sparsify Sigma via Cholesky factor (optional) -------------
    chol.sig <- chol(Sigma.draw)
    dd.sig   <- diag(chol.sig)
    L.draw   <- t(chol.sig / dd.sig)
    DD.draw  <- diag(dd.sig)
    if (sigma.sparse && M > 1) {
      eta <- (Y - X %*% A.draw) %*% t(solve(L.draw))
      L.sparse <- diag(M)
      for (ii in 2:M) {
        b.s <- sparsify_chol(eta[, 1:(ii - 1), drop = FALSE],
                             as.matrix(L.draw[ii, 1:(ii - 1)]),
                             kappa.cov)
        L.sparse[ii, 1:(ii - 1)] <- b.s
      }
      S.draw <- L.sparse %*% (DD.draw^2) %*% t(L.sparse)
    } else {
      S.draw <- Sigma.draw
    }

    # ---- Step 5: Common stochastic volatility ------------------------------
    if (sv == "pca.sv") {
      chol.Sig <- chol(Sigma.draw)
      eps <- (Y - X %*% Alpha.draw) %*% t(solve(chol.Sig))
      # collapse to a single common shock via 1st PC
      xx <- crossprod(eps); ev <- eigen(xx)
      lam <- sqrt(M) * ev$vectors[, 1]
      fac <- as.numeric(eps %*% lam / M)
      sv_out <- suppressWarnings(
        stochvol::svsample_fast_cpp(fac,
                                    draws = 1, burnin = 0,
                                    priorspec   = sv_priors,
                                    startpara   = sv_para,
                                    startlatent = sv_latent)
      )
      sv_para$mu    <- as.numeric(sv_out$para[1, "mu"])
      sv_para$phi   <- as.numeric(sv_out$para[1, "phi"])
      sv_para$sigma <- as.numeric(sv_out$para[1, "sigma"])
      sv_latent     <- as.numeric(sv_out$latent[1, ])
      ht <- sv_latent - sv_latent[1]
      Sig2.t <- exp(ht)
      normalizer <- 1 / sqrt(Sig2.t)

      XX <- rbind(X * normalizer, X.min)
      YY <- rbind(Y * normalizer, Y.min)
      V.post <- try(solve(crossprod(XX)), silent = TRUE)
      if (inherits(V.post, "try-error")) V.post <- ginv(crossprod(XX))
      V.chol <- t(chol(V.post))
      A.post <- V.post %*% crossprod(XX, YY)
      Scale.post <- crossprod(YY - XX %*% A.post)
    } else {
      Sig2.t <- rep(1, T)
    }

    # ---- Step 6: Store + post-processing (IRFs, forecasts) -----------------
    if (irep > nburn) {
      ii <- irep - nburn
      ALPHA.store[ii, , ] <- Alpha.draw
      A.store[ii, , ]     <- A.draw
      SIGMA.store[ii, , ] <- Sigma.draw
      S.store[ii, , ]     <- S.draw
      eht.store[ii, ]     <- Sig2.t
      delta.store[ii, ]   <- delta.draw

      # Stability check
      MM       <- get.companion(Alpha.draw, c(M, cons, p))$MM
      MM.sp    <- get.companion(A.draw,     c(M, cons, p))$MM
      stab     <- max(abs(Re(eigen(MM)$values)))     <= 1.02
      stab.sp  <- max(abs(Re(eigen(MM.sp)$values)))  <= 1.02

      # ---- IRFs via Cholesky decomposition of Sigma -----------------------
      A_arr     <- array(0, c(M, M, p))
      A_arr.sp  <- array(0, c(M, M, p))
      for (k in 1:p) {
        A_arr[, , k]    <- t(Alpha.draw[(((k - 1) * M + 1):(k * M)), ])
        A_arr.sp[, , k] <- t(A.draw[(((k - 1) * M + 1):(k * M)), ])
      }
      if (stab) {
        irf.draw <- impulsdtrf(A_arr, t(chol(Sigma.draw * Sig2.t[T])), nhor)
        irf.store[ii, , , ] <- irf.draw
      }
      if (stab.sp) {
        irf.sp <- impulsdtrf(A_arr.sp, t(chol(S.draw * Sig2.t[T])), nhor)
        irf.sparse.store[ii, , , ] <- irf.sp
      }

      # ---- Forecasts ------------------------------------------------------
      if (cons) {
        X0 <- if (p > 1) c(Y[T, ], X[T, 1:(K - M - 1)], 1) else c(Y[T, ], 1)
      } else {
        X0 <- if (p > 1) c(Y[T, ], X[T, 1:(K - M)]) else c(Y[T, ])
      }
      X0 <- as.matrix(X0); X0.sp <- X0

      Jm <- get.companion(Alpha.draw, c(M, cons, p))$Jm
      Sig0 <- matrix(0, K, K); Sig0.sp <- matrix(0, K, K)

      if (stab) {
        for (h in 1:fhorz) {
          Sig0 <- MM %*% Sig0 %*% t(MM) + Jm %*% (Sigma.draw * Sig2.t[T]) %*% t(Jm)
          X0   <- MM %*% X0
          ch   <- try(t(chol(Sig0[1:M, 1:M])), silent = TRUE)
          if (inherits(ch, "try-error"))
            pred.store[ii, , h] <- X0[1:M] else
              pred.store[ii, , h] <- X0[1:M] + ch %*% rnorm(M)
        }
      }
      if (stab.sp) {
        for (h in 1:fhorz) {
          Sig0.sp <- MM.sp %*% Sig0.sp %*% t(MM.sp) + Jm %*% (S.draw * Sig2.t[T]) %*% t(Jm)
          X0.sp   <- MM.sp %*% X0.sp
          ch.sp   <- try(t(chol(Sig0.sp[1:M, 1:M])), silent = TRUE)
          if (inherits(ch.sp, "try-error"))
            pred.sparse.store[ii, , h] <- X0.sp[1:M] else
              pred.sparse.store[ii, , h] <- X0.sp[1:M] + ch.sp %*% rnorm(M)
        }
      }
    }

    if (verbose && irep %% max(1, floor(ntot / 10)) == 0)
      message(sprintf("  iter %d / %d", irep, ntot))
  }

  # --------------------------- Posterior summaries ---------------------------
  out <- list(
    A.median       = apply(A.store,     c(2, 3), median),
    ALPHA.median   = apply(ALPHA.store, c(2, 3), median),
    SIGMA.median   = apply(SIGMA.store, c(2, 3), median),
    S.median       = apply(S.store,     c(2, 3), median),
    delta.mean     = apply(delta.store, 2, mean),       # PIP-style indicator
    eht.median     = apply(eht.store,   2, median),
    irf.store      = irf.store,
    irf.sparse.store = irf.sparse.store,
    pred.store        = pred.store,
    pred.sparse.store = pred.sparse.store,
    var.names      = var.names,
    p              = p, M = M, K = K, T = T, cons = cons,
    settings       = list(shrink.1 = shrink.1, kappa.base = kappa.base,
                          lambda.1st = lambda.1st, lambda.2nd = lambda.2nd,
                          lambda.decay = lambda.decay, kappa.decay = kappa.decay,
                          kappa.cov = kappa.cov, sigma.sparse = sigma.sparse,
                          sv = sv, nhor = nhor, fhorz = fhorz,
                          nsave = nsave, nburn = nburn)
  )
  class(out) <- "sparsevar"
  out
}
