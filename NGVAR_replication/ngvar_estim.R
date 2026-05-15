# ngvar_estim.R
# ------------------------------------------------------------------------------
# Self-contained estimator for a Bayesian VAR with hierarchical Normal-Gamma
# (NG) shrinkage on the autoregressive coefficients and stochastic volatility
# (SV) on the residuals. The sampler uses an equation-by-equation triangular
# factorization Y = X A + E with E = U^{-1} D^{1/2} eps, where U is unit lower
# triangular (covariance parameters appear as additional regressors in the
# m-th equation) and D is a diagonal matrix of time-varying variances.
#
# Reference: Huber, F. and Feldkircher, M. (2019).
#   "Adaptive shrinkage in Bayesian vector autoregressive models",
#   Journal of Business and Economic Statistics, 37(1), 27-39.
# ------------------------------------------------------------------------------

# ---- small helpers -----------------------------------------------------------

mlag <- function(X, lag) {
  # Build a (T x p*N) matrix of stacked lags 1..p of X (zeros for first p rows).
  p     <- lag
  X     <- as.matrix(X)
  Traw  <- nrow(X)
  N     <- ncol(X)
  Xlag  <- matrix(0, Traw, p * N)
  for (ii in 1:p) {
    Xlag[(p + 1):Traw, (N * (ii - 1) + 1):(N * ii)] <-
      X[(p + 1 - ii):(Traw - ii), 1:N, drop = FALSE]
  }
  Xlag
}

atau_post <- function(atau, lambda2, thetas, k = length(thetas), rat = 1) {
  # Log-posterior of the NG shape parameter a_tau (Gamma-Exp prior).
  sum(dgamma(thetas, atau, (atau * lambda2 / 2), log = TRUE)) +
    dexp(atau, rate = rat, log = TRUE)
}

get_companion <- function(Beta_, varndxv) {
  # Companion-form matrices (MM, Jm) of a VAR(p) with constant flag nd.
  nn <- varndxv[[1]]
  nd <- varndxv[[2]]
  nl <- varndxv[[3]]
  nkk <- nn * nl + nd
  Jm <- matrix(0, nkk, nn)
  Jm[1:nn, 1:nn] <- diag(nn)
  if (nd == 1) {
    MM <- rbind(t(Beta_),
                cbind(diag((nl - 1) * nn), matrix(0, (nl - 1) * nn, nn + nd)),
                c(matrix(0, nd, nn * nl), nd))
  } else {
    MM <- rbind(t(Beta_),
                cbind(diag((nl - 1) * nn), matrix(0, (nl - 1) * nn, nn + nd)))
  }
  list(MM = MM, Jm = Jm)
}

# ---- main estimator ----------------------------------------------------------

ngvar_estim <- function(Yraw,
                        p          = 2,
                        nsave      = 2000,
                        nburn      = 1000,
                        cons       = 1,
                        sv         = TRUE,
                        shrink.type= "lagwise",
                        a_tau      = 0.1,
                        b_tau      = 0.1,
                        c_tau      = 0.01,
                        d_tau      = 0.01,
                        var.c1     = 0.01,
                        var.c2     = 0.01,
                        a_i        = 0.01,
                        b_i        = 0.01,
                        prmean     = 0,
                        sample_a   = TRUE,
                        verbose    = TRUE) {
  # Inputs
  #   Yraw        : (Traw x M) numeric matrix of endogenous variables
  #   p           : lag order
  #   nsave/nburn : posterior draws kept / burn-in
  #   cons        : include intercept (0/1)
  #   sv          : TRUE for stochastic volatility, FALSE for homoskedastic
  #   shrink.type : one of "global", "rowwise", "columnwise", "lagwise"
  #
  # Returns a list with posterior draws of A (VAR coefs), U (lower-tri cov),
  # H (log-volatilities, T x M), SIGMA at final time T (M x M x nsave), and
  # the design matrices.

  require(GIGrvg)
  require(MASS)
  require(mvtnorm)
  require(stochvol)
  require(Matrix)

  ntot <- nsave + nburn

  # ---- build design matrix ---------------------------------------------------
  Yraw <- as.matrix(Yraw)
  if (cons == 1) {
    Xraw <- cbind(mlag(Yraw, p), 1)
  } else {
    Xraw <- mlag(Yraw, p)
  }
  Y <- Yraw[(p + 1):nrow(Yraw), , drop = FALSE]
  X <- Xraw[(p + 1):nrow(Xraw), , drop = FALSE]
  T <- nrow(X); M <- ncol(Y); K <- ncol(X); v <- (M * (M - 1)) / 2

  # ---- OLS starting values ---------------------------------------------------
  A_OLS    <- ginv(crossprod(X)) %*% crossprod(X, Y)
  A_draw   <- A_OLS
  Em       <- Em.str <- Y - X %*% A_OLS
  SIGMA    <- crossprod(Em) / (T - K)

  # ---- priors and initial values --------------------------------------------
  A_prior  <- matrix(0, K, M)
  if (prmean != 0) A_prior[1:M, 1:M] <- diag(M) * prmean
  theta    <- matrix(10, K, M)              # NG local scales (coefs)
  b_draw   <- lower.tri(SIGMA) * 1           # NG local scales (covariances)
  Upsilon_draw <- b_draw
  diag(Upsilon_draw) <- 1                    # unit-diagonal lower triangle
  eta <- vector("list", M)

  # SV initial state per equation
  sv_startpara <- c(mu = -10, phi = 0.9, sigma = 0.2, nu = Inf, rho = 0)
  sv_priorspec <- stochvol::specify_priors(
    mu    = stochvol::sv_normal(mean = 0, sd = 100),
    phi   = stochvol::sv_beta(shape1 = 25, shape2 = 1.5),
    sigma2 = stochvol::sv_gamma(shape = 0.5, rate = 0.5)
  )
  svl <- replicate(M, list(para = sv_startpara, latent = rep(-3, T)),
                   simplify = FALSE)
  H   <- matrix(-3, T, M)
  pars_var <- matrix(0, 3, M)

  # MH tuning
  accept  <- 0; accept2 <- 0
  scale1  <- 0.43; scale2 <- 0.45

  # Lag-wise shrinkage state (used only if shrink.type == "lagwise")
  if (shrink.type == "lagwise") {
    kappa_j <- matrix(0.01, p, 1)
    a_tau   <- matrix(a_tau, p, 1)
  } else if (shrink.type == "rowwise") {
    a_tau   <- matrix(a_tau, M, 1)
  } else if (shrink.type == "columnwise") {
    a_tau   <- matrix(a_tau, K, 1)
  }

  # ---- storage ---------------------------------------------------------------
  A_store       <- array(0, c(nsave, K, M))
  U_store       <- array(0, c(nsave, M, M))
  H_store       <- array(0, c(nsave, T, M))
  SIGMA_T_store <- array(0, c(nsave, M, M))
  pars_store    <- array(0, c(nsave, 3, M))

  if (verbose) pb <- txtProgressBar(min = 0, max = ntot, style = 3)
  t0 <- Sys.time()

  for (irep in 1:ntot) {

    # ---- Step I: equation-by-equation draw of A and the U-rows ---------------
    for (mm in 1:M) {
      if (mm == 1) {
        Y.i <- Y[, mm] * exp(-0.5 * H[, mm])
        X.i <- X * exp(-0.5 * H[, mm])
        V_post <- try(chol2inv(chol(crossprod(X.i) + diag(1 / theta[, mm]))),
                      silent = TRUE)
        if (inherits(V_post, "try-error"))
          V_post <- ginv(crossprod(X.i) + diag(1 / theta[, mm]))
        A_post <- V_post %*% (crossprod(X.i, Y.i) +
                                diag(1 / theta[, mm]) %*% A_prior[, mm])
        chV <- try(t(chol(V_post)), silent = TRUE)
        A.draw.i <- if (inherits(chV, "try-error"))
          mvrnorm(1, A_post, V_post) else A_post + chV %*% rnorm(K)
        A_draw[, mm] <- A.draw.i
        Em[, mm]     <- Em.str[, mm] <- Y[, mm] - X %*% A.draw.i
      } else {
        Y.i <- Y[, mm] * exp(-0.5 * H[, mm])
        X.i <- cbind(X, Em[, 1:(mm - 1), drop = FALSE]) * exp(-0.5 * H[, mm])
        prior_var <- c(theta[, mm], b_draw[mm, 1:(mm - 1)])
        prior_mean <- c(A_prior[, mm], rep(0, mm - 1))
        V_post <- try(chol2inv(chol(crossprod(X.i) + diag(1 / prior_var))),
                      silent = TRUE)
        if (inherits(V_post, "try-error"))
          V_post <- ginv(crossprod(X.i) + diag(1 / prior_var))
        A_post <- V_post %*% (crossprod(X.i, Y.i) +
                                diag(1 / prior_var) %*% prior_mean)
        chV <- try(t(chol(V_post)), silent = TRUE)
        A.draw.i <- if (inherits(chV, "try-error"))
          mvrnorm(1, A_post, V_post) else A_post + chV %*% rnorm(ncol(X.i))
        A_draw[, mm] <- A.draw.i[1:K]
        Em[, mm]     <- Y[, mm] - X %*% A.draw.i[1:K]
        Em.str[, mm] <- Em[, mm] -
          Em[, 1:(mm - 1), drop = FALSE] %*% A.draw.i[(K + 1):ncol(X.i)]
        Upsilon_draw[mm, 1:(mm - 1)] <- eta[[mm]] <- A.draw.i[(K + 1):ncol(X.i)]
      }
    }

    # ---- Step II: NG shrinkage on coefficients -------------------------------
    if (shrink.type == "global") {

      lambda2_tau <- rgamma(1, c_tau + a_tau * K * M,
                            d_tau + a_tau / 2 * sum(theta))
      for (mm in 1:M) for (nn in 1:K) {
        theta[nn, mm] <- rgig(n = 1, lambda = a_tau - 0.5,
                              (A_draw[nn, mm] - A_prior[nn, mm])^2,
                              a_tau * lambda2_tau)
      }
      theta[theta < 1e-9] <- 1e-9
      if (sample_a) {
        a_prop <- exp(rnorm(1, 0, scale1)) * a_tau
        d_post <- atau_post(a_prop, lambda2_tau, as.vector(theta)) -
                  atau_post(a_tau,  lambda2_tau, as.vector(theta))
        if (!is.nan(d_post) && d_post > log(runif(1))) {
          a_tau <- a_prop; accept <- accept + 1
        }
        if (irep < 0.5 * nburn) {
          if (accept / irep > 0.30) scale1 <- 1.01 * scale1
          if (accept / irep < 0.15) scale1 <- 0.99 * scale1
        }
      }

    } else if (shrink.type == "rowwise") {

      for (jj in 1:M) {
        lambda2_tau_j <- rgamma(1, c_tau + a_tau[jj] * K,
                                d_tau + a_tau[jj] / 2 * sum(theta[, jj]))
        for (ii in 1:K) {
          theta[ii, jj] <- rgig(n = 1, lambda = a_tau[jj] - 0.5,
                                (A_draw[ii, jj] - A_prior[ii, jj])^2,
                                a_tau[jj] * lambda2_tau_j)
        }
        if (sample_a) {
          a_prop <- exp(rnorm(1, 0, scale1)) * a_tau[jj]
          d_post <- atau_post(a_prop, lambda2_tau_j, theta[, jj]) -
                    atau_post(a_tau[jj], lambda2_tau_j, theta[, jj])
          if (!is.nan(d_post) && d_post > log(runif(1))) a_tau[jj] <- a_prop
        }
      }
      theta[theta < 1e-8] <- 1e-8

    } else if (shrink.type == "columnwise") {

      for (jj in 1:K) {
        lambda2_tau_j <- rgamma(1, c_tau + a_tau[jj] * M,
                                d_tau + a_tau[jj] / 2 * sum(theta[jj, ]))
        for (ii in 1:M) {
          theta[jj, ii] <- rgig(n = 1, lambda = a_tau[jj] - 0.5,
                                (A_draw[jj, ii] - A_prior[jj, ii])^2,
                                a_tau[jj] * lambda2_tau_j)
        }
        if (sample_a) {
          a_prop <- exp(rnorm(1, 0, scale1)) * a_tau[jj]
          d_post <- atau_post(a_prop, lambda2_tau_j, theta[jj, ]) -
                    atau_post(a_tau[jj], lambda2_tau_j, theta[jj, ])
          if (!is.nan(d_post) && d_post > log(runif(1))) a_tau[jj] <- a_prop
        }
      }
      theta[theta < 1e-7] <- 1e-7

    } else if (shrink.type == "lagwise") {

      for (ss in 1:p) {
        idx       <- ((ss - 1) * M + 1):(ss * M)
        A.lag     <- A_draw[idx, , drop = FALSE]
        A.prior.l <- A_prior[idx, , drop = FALSE]
        theta.lag <- theta[idx, , drop = FALSE]
        if (ss == 1) {
          kappa_j[ss, 1] <- rgamma(1, var.c1 + a_tau[ss] * M^2,
                                   var.c2 + a_tau[ss] / 2 * sum(theta.lag))
        } else {
          kappa_j[ss, 1] <- rgamma(1, var.c1 + a_tau[ss] * M^2,
                                   var.c2 + a_tau[ss] / 2 *
                                     prod(kappa_j[1:(ss - 1), 1]) *
                                     sum(theta.lag))
        }
        for (jj in 1:M) for (ii in 1:M) {
          sc <- prod(kappa_j[1:ss, 1])
          if (sc == 0) sc <- 1e-8
          theta.lag[jj, ii] <- rgig(n = 1, lambda = a_tau[ss] - 0.5,
                                    (A.lag[jj, ii] - A.prior.l[jj, ii])^2,
                                    a_tau[ss] * sc)
        }
        theta[idx, ] <- theta.lag
        theta[theta < 1e-7] <- 1e-7
        if (sample_a) {
          a_prop <- exp(rnorm(1, 0, scale1)) * a_tau[ss]
          lam_s  <- prod(kappa_j[1:ss, 1])
          d_post <- atau_post(a_prop, lam_s, as.vector(theta.lag)) -
                    atau_post(a_tau[ss], lam_s, as.vector(theta.lag))
          if (!is.nan(d_post) && d_post > log(runif(1))) a_tau[ss] <- a_prop
        }
      }
    }

    # ---- Step III: SV / homoskedastic update ---------------------------------
    if (sv) {
      for (jj in 1:M) {
        svdraw <- stochvol::svsample_fast_cpp(
          Em.str[, jj],
          draws       = 1, burnin = 0,
          priorspec   = sv_priorspec,
          startpara   = svl[[jj]]$para,
          startlatent = svl[[jj]]$latent)
        svl[[jj]]$para   <- svdraw$para[1, ]
        svl[[jj]]$latent <- svdraw$latent[1, ]
        H[, jj]          <- svl[[jj]]$latent
        # pars_var rows are (mu, phi, sigma) -- sigma is the standard deviation
        # of the AR(1) innovation in log-volatility.
        pars_var[, jj]   <- svl[[jj]]$para[c("mu", "phi", "sigma")]
      }
    } else {
      for (jj in 1:M) {
        S_1 <- a_i / 2 + T / 2
        S_2 <- b_i / 2 + crossprod(Em.str[, jj]) / 2
        sig_eta <- 1 / rgamma(1, S_1, S_2)
        H[, jj] <- log(sig_eta)
      }
    }
    SIGMA <- Upsilon_draw %*% diag(exp(H[T, ]), M) %*% t(Upsilon_draw)

    # ---- Step IV: NG shrinkage on the off-diagonal covariance entries --------
    lambda_omega <- rgamma(1, c_tau + b_tau * v,
                           d_tau + b_tau / 2 *
                             sum(b_draw[lower.tri(b_draw)]))
    for (nn in 2:M) {
      et_0 <- eta[[nn]]
      for (jj in 1:(nn - 1)) {
        b_draw[nn, jj] <- rgig(n = 1, lambda = b_tau - 0.5,
                               et_0[jj]^2, b_tau * lambda_omega)
      }
    }
    if (sample_a) {
      b_prop <- exp(rnorm(1, 0, scale2)) * b_tau
      tv     <- b_draw[lower.tri(b_draw)]
      d_post <- atau_post(b_prop, lambda_omega, tv, k = length(tv)) +
        dlnorm(b_tau, log(b_prop), scale2, log = TRUE) -
        atau_post(b_tau, lambda_omega, tv, k = length(tv)) -
        dlnorm(b_prop, log(b_tau), scale2, log = TRUE)
      if (!is.nan(d_post) && d_post > log(runif(1))) {
        b_tau <- b_prop; accept2 <- accept2 + 1
      }
      if (irep < 0.5 * nburn) {
        if (accept2 / irep > 0.30) scale2 <- 1.01 * scale2
        if (accept2 / irep < 0.15) scale2 <- 0.99 * scale2
      }
    }

    # ---- store ---------------------------------------------------------------
    if (irep > nburn) {
      ss <- irep - nburn
      A_store[ss, , ]       <- A_draw
      U_store[ss, , ]       <- Upsilon_draw
      H_store[ss, , ]       <- H
      SIGMA_T_store[ss, , ] <- SIGMA
      pars_store[ss, , ]    <- pars_var
    }

    if (verbose) setTxtProgressBar(pb, irep)
  }
  if (verbose) {
    close(pb); cat("\nMCMC finished in ",
                   format(Sys.time() - t0), "\n", sep = "")
  }

  list(A         = A_store,
       U         = U_store,
       H         = H_store,
       SIGMA_T   = SIGMA_T_store,
       sv_pars   = pars_store,
       Y         = Y,
       X         = X,
       Yraw      = Yraw,
       p         = p,
       cons      = cons,
       M         = M,
       K         = K,
       T         = T,
       shrink.type = shrink.type)
}

# ---- companion-form impulse responses (Cholesky identification) --------------

ngvar_irf <- function(fit, nhor = 20) {
  # Compute orthogonalised IRFs for every retained draw using a recursive
  # Cholesky factor of the contemporaneous covariance matrix.
  #   fit  : output of ngvar_estim()
  #   nhor : horizon (including impact period h = 0)
  # Returns an array of dimension (nsave x M x M x nhor).
  M     <- fit$M; p <- fit$p; cons <- fit$cons; K <- fit$K
  nsave <- dim(fit$A)[1]
  IRF   <- array(NA_real_, dim = c(nsave, M, M, nhor))
  for (s in 1:nsave) {
    A_s     <- fit$A[s, , ]
    Sigma_s <- fit$SIGMA_T[s, , ]
    cmp     <- get_companion(A_s, c(M, cons, p))
    Mm      <- cmp$MM
    if (max(abs(Re(eigen(Mm, only.values = TRUE)$values))) > 1.05) next
    cholS <- try(t(chol(Sigma_s)), silent = TRUE)
    if (inherits(cholS, "try-error")) next
    # state-space simulation: state = (y_t', y_{t-1}', ..., y_{t-p+1}', 1)
    nkk <- M * p + cons
    for (sh in 1:M) {
      shock <- rep(0, nkk)
      shock[1:M] <- cholS[, sh]
      IRF[s, , sh, 1] <- shock[1:M]
      st <- shock
      for (h in 2:nhor) {
        st <- Mm %*% st
        IRF[s, , sh, h] <- st[1:M]
      }
    }
  }
  IRF
}
