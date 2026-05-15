# =============================================================================
# gpvar_estim.R
# -----------------------------------------------------------------------------
# Self-contained estimator for the Gaussian Process Vector Autoregression
# (GP-VAR) of "Gaussian Process Vector Autoregressions and Macroeconomic
# Uncertainty".
#
# Equation-by-equation MCMC. For each equation we build two Gaussian (squared-
# exponential / "RBF") kernel matrices: one over the own lags (h.own) and one
# over the other-variable lags (h.other). Conditional on the kernels we draw
# the latent GP functionals f.1 (own) and f.2 (other) via a multivariate
# normal whose covariance is the standard GP regression posterior
#     K - K (K + diag(sig2))^{-1} K.
# We also sample (i) equation-wise stochastic volatility via stochvol, (ii)
# the lower-triangular structural covariance Upsilon with a horseshoe prior,
# and (iii) kernel-shrinkage parameters psi[1], psi[2] via an adaptive
# random-walk Metropolis ("adptMinn") step.
#
# Differences relative to cluster.R:
#   * No C++ BAKRGibbs / no rsvd / no sparseMVN.
#   * No subspace approximation: we use the full T x T kernel matrices.
#   * MCMC budget passed as arguments, with small defaults; T must stay small
#     because every Gibbs step is O(T^3) per equation.
#   * Returns posterior draws of f.1, f.2, log-volatilities ht, Upsilon,
#     kernel hyperparameters psi, plus a generalized impulse response array.
# =============================================================================

suppressPackageStartupMessages({
  library(MASS)
  library(mvtnorm)
  library(Matrix)
  library(stochvol)
  library(GIGrvg)
  library(truncnorm)
  library(invgamma)
})

# -----------------------------------------------------------------------------
# Helpers (extracted from aux_funcs.R)
# -----------------------------------------------------------------------------

#' Squared-exponential ("Gaussian") kernel
#'
#' K(x*, x) = sigf * exp( -(h / p) * sum_j (x*_j - x_j)^2 )
#' where columns of Xin / Xstar are observations and rows are features (p).
GaussKernel <- function(Xin, Xstar, h, sigf) {
  nstar <- ncol(Xstar)
  nin   <- ncol(Xin)
  p     <- nrow(Xstar)
  K <- matrix(0, nstar, nin)
  for (ii in seq_len(nstar)) {
    K[ii, ] <- sigf * exp((-h / p) *
                            apply(as.matrix((Xstar[, ii] - Xin[, ])^2), 2, sum))
  }
  return(K)
}

#' Build a matrix of stacked lags
mlag <- function(X, lag) {
  X <- as.matrix(X)
  Traw <- nrow(X); N <- ncol(X)
  Xlag <- matrix(0, Traw, lag * N)
  for (ii in seq_len(lag)) {
    Xlag[(lag + 1):Traw, (N * (ii - 1) + 1):(N * ii)] <-
      X[(lag + 1 - ii):(Traw - ii), 1:N]
  }
  return(Xlag)
}

#' One Gibbs update of the horseshoe scales given a vector of regression
#' coefficients bdraw. Returns the new (lambda, tau, nu, zeta) and the
#' implied prior variances psi = lambda * tau.
get.hs <- function(bdraw, lambda.hs, nu.hs, tau.hs, zeta.hs) {
  k <- length(bdraw)
  if (is.na(tau.hs)) {
    tau.hs <- 1
  } else {
    tau.hs <- invgamma::rinvgamma(1, shape = (k + 1) / 2,
                                  rate = 1 / zeta.hs + sum(bdraw^2 / lambda.hs) / 2)
  }
  lambda.hs <- invgamma::rinvgamma(k, shape = 1,
                                   rate = 1 / nu.hs + bdraw^2 / (2 * tau.hs))
  nu.hs   <- invgamma::rinvgamma(k, shape = 1, rate = 1 + 1 / lambda.hs)
  zeta.hs <- invgamma::rinvgamma(1, shape = 1, rate = 1 + 1 / tau.hs)
  list(psi = lambda.hs * tau.hs, lambda = lambda.hs,
       tau = tau.hs, nu = nu.hs, zeta = zeta.hs)
}

# -----------------------------------------------------------------------------
# Main estimator
# -----------------------------------------------------------------------------

#' Fit the GP-VAR
#'
#' @param Y       (T x MM) data matrix. Assumed already standardised.
#' @param p       Number of lags.
#' @param nsave   Number of retained MCMC draws.
#' @param nburn   Number of burn-in draws.
#' @param sv      Logical, use stochvol for equation variances. Default TRUE.
#' @param h.own   Bandwidth on the own-lag kernel.
#' @param h.other Bandwidth on the other-variable-lag kernel.
#' @param shk.var Index (or name) of the variable that is shocked in the GIRF.
#' @param shk.sc  Size of the shock (in standard deviations).
#' @param nhor    Horizon of the GIRF.
#' @param girf.on Whether to compute generalized impulse responses.
#' @param verbose Print progress every 25 iterations.
gpvar <- function(Y,
                  p       = 2,
                  nsave   = 500,
                  nburn   = 500,
                  sv      = TRUE,
                  h.own   = 1,
                  h.other = 1,
                  shk.var = 1,
                  shk.sc  = 1,
                  nhor    = 17,
                  girf.on = TRUE,
                  verbose = TRUE) {

  Y <- as.matrix(Y)
  if (is.null(colnames(Y))) colnames(Y) <- paste0("y", seq_len(ncol(Y)))
  var.names <- colnames(Y)

  # ---- design --------------------------------------------------------------
  Xraw <- mlag(Y, p)
  y <- as.matrix(Y[(p + 1):nrow(Y), , drop = FALSE])
  X <- as.matrix(Xraw[(p + 1):nrow(Xraw), , drop = FALSE])
  MM <- ncol(y)
  T  <- nrow(y)
  K  <- ncol(X)

  # ---- shock matrix for GIRF -----------------------------------------------
  shk.mat <- matrix(0, MM, MM)
  colnames(shk.mat) <- rownames(shk.mat) <- var.names
  if (is.character(shk.var)) {
    shk.mat[, shk.var] <- shk.sc
  } else {
    shk.mat[, shk.var] <- shk.sc
  }

  # ---- kernels (precomputed once) ------------------------------------------
  Kn.list <- vector("list", MM)
  Kn.invChol.list <- vector("list", MM)
  Dinv.list <- vector("list", MM)
  for (km in seq_len(MM)) {
    own.idx   <- seq(km, K, by = MM)
    other.idx <- setdiff(seq_len(K), own.idx)
    Kn1 <- GaussKernel(t(X[, own.idx]),   t(X[, own.idx]),   h = h.own,   sigf = 1)
    Kn2 <- GaussKernel(t(X[, other.idx]), t(X[, other.idx]), h = h.other, sigf = 1)
    diag(Kn1) <- 1; diag(Kn2) <- 1
    Kn.list[[km]] <- list(Kn1, Kn2)
    Kn.invChol.list[[km]] <- list(
      solve(chol(Kn1 + diag(1e-4, T))),
      solve(chol(Kn2 + diag(1e-6, T)))
    )
  }
  Xginv <- Matrix(MASS::ginv(X), sparse = TRUE)

  # ---- priors / state --------------------------------------------------------
  Q <- 2
  shrnk.own   <- 1
  shrnk.other <- 1
  psi.mat <- array(NA, c(Q, MM, T))
  psi.mat[1, , ] <- shrnk.own
  psi.mat[2, , ] <- shrnk.other

  # adaptive Minnesota-style RW-MH for the two kernel scaling parameters
  sc     <- c(0.02, 0.02)
  a.tau  <- c(1e-6, 1e-6)
  accept <- c(0, 0)

  # OLS init
  beta.ols <- MASS::ginv(crossprod(X)) %*% crossprod(X, y)
  Em <- Em.str <- y - X %*% beta.ols
  sigma.ols <- crossprod(y - X %*% beta.ols) / (T - K)

  # Stochastic volatility
  if (sv) {
    sv_priors <- specify_priors(
      mu     = sv_normal(mean = 0, sd = 100),
      phi    = sv_beta(shape1 = 25, shape2 = 1.5),
      sigma2 = sv_gamma(shape = 0.5, rate = 5),
      nu     = sv_infinity(),
      rho    = sv_constant(0)
    )
    svdraw <- list(mu = 0, phi = 0.99, sigma = 0.01, nu = Inf, rho = 0,
                   beta = NA, latent0 = 0)
    h_sv <- replicate(MM, svdraw, simplify = FALSE)
  } else {
    a.sig <- 1 / 2; b.sig <- 1 / 2
  }
  sigma2.draw <- t(matrix(diag(sigma.ols), MM, T))
  ht <- t(matrix(log(diag(sigma.ols)), MM, T))

  # Lower-triangular covariance + horseshoe state
  Upsilon.draw <- lower.tri(diag(MM)) * 1
  diag(Upsilon.draw) <- 1
  zeta <- Upsilon.draw

  Upsdraw.hs.vec <- c()
  tau.hs.vec  <- 1
  zeta.hs.vec <- 1
  n.cov <- (MM * (MM - 1)) / 2
  nu.hs.mat  <- rep(1, n.cov)
  lam.hs.mat <- rep(1, n.cov)
  psi.hs.mat <- rep(1, n.cov)

  # Latent GP draws
  f.1.mat <- matrix(rnorm(T * MM, 0, 1), T, MM)
  f.2.mat <- matrix(rnorm(T * MM, 0, 1), T, MM)
  f.draw   <- matrix(NA, T, Q)
  f.tilde  <- matrix(NA, T, Q)

  # Storage
  ntot <- nsave + nburn
  f.1.store    <- array(NA, c(nsave, T, MM))
  f.2.store    <- array(NA, c(nsave, T, MM))
  beta.store   <- array(NA, c(nsave, K, MM))
  ht.store     <- array(NA, c(nsave, T, MM))
  Upsilon.store <- array(NA, c(nsave, MM, MM))
  psi.store    <- array(NA, c(nsave, Q, MM))
  girf_store   <- array(0, c(nsave, MM, nhor, 2))
  dimnames(girf_store) <- list(seq_len(nsave), var.names, seq_len(nhor),
                               c("pos", "neg"))

  # =========================================================================
  # Gibbs sampler
  # =========================================================================
  for (irep in seq_len(ntot)) {

    for (mm in seq_len(MM)) {

      own.idx   <- seq(mm, K, by = MM)
      other.idx <- setdiff(seq_len(K), own.idx)

      # ---- Step 1: latent GP for "other" lags -------------------------------
      Kn2_eq <- diag(sqrt(psi.mat[2, mm, ])) %*% Kn.list[[mm]][[2]] %*%
        diag(sqrt(psi.mat[2, mm, ]))
      Dinv2 <- solve(Kn2_eq + diag(sigma2.draw[, mm]) + diag(1e-8, T))
      if (mm == 1) {
        resid_for_f2 <- y[, mm] - f.1.mat[, mm]
      } else {
        resid_for_f2 <- y[, mm] - f.1.mat[, mm] -
          Em[, 1:(mm - 1), drop = FALSE] %*% Upsilon.draw[mm, 1:(mm - 1)]
      }
      fhat <- as.numeric(Kn2_eq %*% Dinv2 %*% resid_for_f2)
      Vhat <- Kn2_eq - Kn2_eq %*% Dinv2 %*% Kn2_eq
      Vhat <- (Vhat + t(Vhat)) / 2 + diag(1e-8, T)
      f2.draw <- tryCatch(t(mvtnorm::rmvnorm(1, fhat, Vhat)),
                          error = function(e)
                            t(mvtnorm::rmvnorm(1, fhat, as.matrix(forceSymmetric(Vhat)))))
      f.2.mat[, mm] <- f2.draw

      # ---- latent GP for "own" lags ----------------------------------------
      Kn1_eq <- diag(sqrt(psi.mat[1, mm, ])) %*% Kn.list[[mm]][[1]] %*%
        diag(sqrt(psi.mat[1, mm, ]))
      Dinv1 <- solve(Kn1_eq + diag(sigma2.draw[, mm]) + diag(1e-8, T))
      if (mm == 1) {
        resid_for_f1 <- y[, mm] - f.2.mat[, mm]
      } else {
        resid_for_f1 <- y[, mm] - f.2.mat[, mm] -
          Em[, 1:(mm - 1), drop = FALSE] %*% Upsilon.draw[mm, 1:(mm - 1)]
      }
      fhat <- as.numeric(Kn1_eq %*% Dinv1 %*% resid_for_f1)
      Vhat <- Kn1_eq - Kn1_eq %*% Dinv1 %*% Kn1_eq
      Vhat <- (Vhat + t(Vhat)) / 2 + diag(1e-8, T)
      f1.draw <- tryCatch(t(mvtnorm::rmvnorm(1, fhat, Vhat)),
                          error = function(e)
                            t(mvtnorm::rmvnorm(1, fhat, as.matrix(forceSymmetric(Vhat)))))
      f.1.mat[, mm] <- f1.draw

      # ---- Step 2: covariance row of Upsilon (eqs mm > 1) -------------------
      if (mm > 1) {
        Y.i <- (y[, mm] - f.1.mat[, mm] - f.2.mat[, mm]) / sqrt(sigma2.draw[, mm])
        X.i <- Em[, 1:(mm - 1), drop = FALSE] / sqrt(sigma2.draw[, mm])
        K.i <- ncol(X.i)
        if (mm == 2) {
          V.A.pr.inv <- zeta[mm, 1:(mm - 1)]
        } else {
          V.A.pr.inv <- diag(zeta[mm, 1:(mm - 1)])
        }
        V.A_post <- tryCatch(solve(crossprod(X.i) + V.A.pr.inv),
                             error = function(e) MASS::ginv(crossprod(X.i) + V.A.pr.inv))
        A.mean   <- V.A_post %*% crossprod(X.i, Y.i)
        A_draw   <- tryCatch(A.mean + t(chol(V.A_post)) %*% rnorm(K.i, 0, 1),
                             error = function(e) t(mvtnorm::rmvnorm(1, A.mean, V.A_post)))
        Upsilon.draw[mm, 1:(mm - 1)] <- A_draw[1:K.i]
        Upsdraw.hs.vec <- c(Upsdraw.hs.vec, A_draw[1:K.i])
        Em[, mm]     <- y[, mm] - f.1.mat[, mm] - f.2.mat[, mm]
        Em.str[, mm] <- Em[, mm] - Em[, 1:(mm - 1), drop = FALSE] %*%
          Upsilon.draw[mm, 1:(mm - 1)]
      } else {
        Em[, mm] <- Em.str[, mm] <- y[, mm] - f.1.mat[, mm] - f.2.mat[, mm]
      }

      Dinv.list[[mm]] <- list(Dinv1, Dinv2)

      # ---- Step 3: kernel shrinkage psi via adaptive RW-MH ----------------
      f.draw[, 1] <- f1.draw
      f.draw[, 2] <- f2.draw
      f.tilde[, 1] <- as.numeric(Kn.invChol.list[[mm]][[1]] %*% f1.draw)
      f.tilde[, 2] <- as.numeric(Kn.invChol.list[[mm]][[2]] %*% f2.draw)
      for (qq in seq_len(Q)) {
        psi.acc  <- psi.mat[qq, mm, 1]
        psi.prop <- truncnorm::rtruncnorm(1, a = 0, b = 5,
                                          mean = psi.acc, sd = sc[qq])
        post.prop <- sum(dnorm(f.tilde[, qq], 0, sqrt(psi.prop), log = TRUE)) +
          dgamma(psi.prop, 1 / 2, 1 / (2 * a.tau[qq]), log = TRUE)
        post.acc  <- sum(dnorm(f.tilde[, qq], 0, sqrt(psi.acc), log = TRUE)) +
          dgamma(psi.acc, 1 / 2, 1 / (2 * a.tau[qq]), log = TRUE)
        if ((post.prop - post.acc) > log(runif(1))) {
          psi.mat[qq, mm, ] <- psi.prop
          accept[qq] <- accept[qq] + 1
        }
        if (irep < (0.5 * nburn)) {
          if ((accept[qq] / irep) > 0.3)  sc[qq] <- 1.01 * sc[qq]
          if ((accept[qq] / irep) < 0.15) sc[qq] <- 0.99 * sc[qq]
        }
      }
      psi.mat[psi.mat > 100]  <- 100
      psi.mat[psi.mat < 1e-8] <- 1e-8

      # ---- Step 4: variances (SV or IG) ------------------------------------
      if (sv) {
        svdraw <- svsample_fast_cpp(Em.str[, mm], startpara = h_sv[[mm]],
                                    startlatent = ht[, mm], priorspec = sv_priors)
        svdraw[c("mu", "phi", "sigma", "nu", "rho")] <-
          as.list(svdraw$para[, c("mu", "phi", "sigma", "nu", "rho")])
        h_sv[[mm]] <- svdraw
        ht.mm <- as.numeric(svdraw$latent)
        ht.mm[ht.mm < -20] <- -20
        ht[, mm] <- ht.mm
        sigma2.draw[, mm] <- exp(ht.mm)
      } else {
        a.star <- a.sig + T / 2
        b.star <- b.sig + as.numeric(crossprod(Em.str[, mm])) / 2
        sig2.mm <- 1 / rgamma(1, a.star, b.star)
        sigma2.draw[, mm] <- rep(sig2.mm, T)
        ht[, mm] <- rep(log(sig2.mm), T)
      }

      # ---- Storage ----------------------------------------------------------
      if (irep > nburn) {
        idx <- irep - nburn
        f.1.store[idx, , mm] <- f.1.mat[, mm]
        f.2.store[idx, , mm] <- f.2.mat[, mm]
        f <- f.1.mat[, mm] + f.2.mat[, mm]
        beta.store[idx, , mm] <- as.matrix(Xginv %*% f)
        ht.store[idx, , mm]   <- ht[, mm]
        psi.store[idx, , mm]  <- psi.mat[, mm, 1]
      }
    }  # end loop over equations mm

    if (irep > nburn) Upsilon.store[irep - nburn, , ] <- Upsilon.draw

    # ---- Step 5: horseshoe on lower-triangular covariances -----------------
    if (length(Upsdraw.hs.vec) == n.cov) {
      hs <- get.hs(bdraw = Upsdraw.hs.vec, lambda.hs = lam.hs.mat,
                   nu.hs = nu.hs.mat, tau.hs = tau.hs.vec,
                   zeta.hs = zeta.hs.vec)
      lam.hs.mat  <- hs$lambda
      nu.hs.mat   <- hs$nu
      tau.hs.vec  <- hs$tau
      zeta.hs.vec <- hs$zeta
      psi.hs.mat  <- pmin(hs$psi + 1e-8, 10)
      # update prior precisions for next sweep
      tmp <- Upsdraw.hs.vec
      for (rv in 2:MM) {
        zeta[rv, 1:(rv - 1)] <- tmp[1:(rv - 1)]
        tmp <- tmp[-(1:(rv - 1))]
      }
    }
    Upsdraw.hs.vec <- c()

    # ---- Step 6: Generalized impulse responses -----------------------------
    if (girf.on && irep > nburn) {
      idx <- irep - nburn
      for (gc in 1:2) {
        sign <- if (gc == 1) 1 else -1
        girf <- matrix(0, MM, nhor)
        rownames(girf) <- var.names
        shock.chol <- Upsilon.draw * shk.mat * sign
        girf[, 1] <- shock.chol[, shk.var]
        X1.girf <- matrix(c(girf[, 1], rep(0, MM * (p - 1))), 1, K)
        X0.girf <- matrix(c(rep(0, MM * p)), 1, K)
        Kn1.out <- Kn0.out <- 1
        for (hh in 2:nhor) {
          fhat1.out <- fhat0.out <- numeric(MM)
          for (mm in seq_len(MM)) {
            own.ident   <- seq(mm, K, by = MM)
            other.ident <- setdiff(seq_len(K), own.ident)
            if (mm == 1) {
              cfe.1 <- y[, mm] - f.2.mat[, mm]
              cfe.2 <- y[, mm] - f.1.mat[, mm]
            } else {
              cfe.1 <- y[, mm] - f.1.mat[, mm] -
                Em[, 1:(mm - 1), drop = FALSE] %*% Upsilon.draw[mm, 1:(mm - 1)]
              cfe.2 <- y[, mm] - f.2.mat[, mm] -
                Em[, 1:(mm - 1), drop = FALSE] %*% Upsilon.draw[mm, 1:(mm - 1)]
            }
            sigf.oth <- sqrt(psi.mat[2, mm, 1])
            sigf.own <- sqrt(psi.mat[1, mm, 1])

            KKn1 <- GaussKernel(t(X[, other.ident]),
                                t(X1.girf[, other.ident, drop = FALSE]),
                                h = h.other, sigf = sigf.oth)
            other1.out <- as.numeric(KKn1 %*% Dinv.list[[mm]][[2]] %*% cfe.1)
            KKn0 <- GaussKernel(t(X[, other.ident]),
                                t(X0.girf[, other.ident, drop = FALSE]),
                                h = h.other, sigf = sigf.oth)
            other0.out <- as.numeric(KKn0 %*% Dinv.list[[mm]][[2]] %*% cfe.1)

            KKn1 <- GaussKernel(t(X[, own.ident]),
                                t(X1.girf[, own.ident, drop = FALSE]),
                                h = h.own, sigf = sigf.own)
            own1.out <- as.numeric(KKn1 %*% Dinv.list[[mm]][[1]] %*% cfe.2)
            KKn0 <- GaussKernel(t(X[, own.ident]),
                                t(X0.girf[, own.ident, drop = FALSE]),
                                h = h.own, sigf = sigf.own)
            own0.out <- as.numeric(KKn0 %*% Dinv.list[[mm]][[1]] %*% cfe.2)

            fhat1.out[mm] <- other1.out + own1.out
            fhat0.out[mm] <- other0.out + own0.out
          }
          girf[, hh] <- fhat1.out - fhat0.out
          if (p > 1) {
            X1.girf <- matrix(c(fhat1.out, X1.girf[1:(MM * (p - 1))]), 1, K)
            X0.girf <- matrix(c(fhat0.out, X0.girf[1:(MM * (p - 1))]), 1, K)
          } else {
            X1.girf <- matrix(fhat1.out, 1, K)
            X0.girf <- matrix(fhat0.out, 1, K)
          }
        }
        girf_store[idx, , , gc] <- girf
      }
    }

    if (verbose && irep %% 25 == 0) cat("MCMC iter", irep, "/", ntot, "\n")
  }

  list(f.1     = f.1.store,
       f.2     = f.2.store,
       ht      = ht.store,
       beta    = beta.store,
       Upsilon = Upsilon.store,
       psi     = psi.store,
       girf    = girf_store,
       y       = y,
       X       = X,
       p       = p,
       MM      = MM,
       T       = T,
       var.names = var.names)
}
