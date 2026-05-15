# ============================================================================
# Time-Varying Parameter regression via Singular Value Decomposition (SVD)
# ----------------------------------------------------------------------------
# Hauzenberger, Huber, Koop & Onorante (2021, JBES):
#   "Fast and Flexible Bayesian Inference in Time-Varying Parameter
#    Regression Models"
#
# This file extracts the headline SVD-based estimator TVPSVD_estim() from
# the original replication package and packages it as a self-contained
# function. It supports
#   * "WN-SVD": white-noise state evolution with a sparse-mixture / g-prior
#               (the headline contribution of the paper)
#   * "RW-SVD": random-walk state evolution (the alternative SVD variant)
# Stochastic volatility is included via the stochvol package.
#
# The original NKPC-specific regularisation of pr.sc has been replaced
# by a generic uniform scaling so that the routine can be applied to any
# y = X beta_t + e_t setup.
# ============================================================================

require(coda)
require(GIGrvg)
require(MASS)
require(Matrix)
require(stochvol)
require(bayesm)

# ----- helpers -------------------------------------------------------------

# create lags of a matrix
mlag <- function(X, lag){
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X); N <- ncol(X)
  Xlag <- matrix(0, Traw, p*N)
  for (ii in 1:p){
    Xlag[(p+1):Traw, (N*(ii-1)+1):(N*ii)] <- X[(p+1-ii):(Traw-ii), (1:N)]
  }
  return(Xlag)
}

# Normal-gamma shrinkage update
get.ng <- function(bdraw, psi, pr = 0, e0 = 0.01, e1 = 0.01,
                   sample_a = FALSE, a_tau = 0.1,
                   acc = 0, pscale = 0.1, irep = 1, nburn = 1){
  k <- NROW(bdraw)
  if (pr == 0) pr <- rep(0, k)
  e0_po <- e0 + a_tau*k
  e1_po <- e1 + (a_tau * sum(psi))/2
  lambda <- rgamma(1, e0_po, e1_po)
  psi <- matrix(0, k, 1)
  for (kk in 1:k){
    scale <- ifelse((bdraw[kk]-pr[kk])^2 <= 1e-12, 1e-12, (bdraw[kk]-pr[kk])^2)
    psi[kk] <- GIGrvg::rgig(n = 1, lambda = a_tau-0.5, chi = scale,
                            psi = lambda*a_tau)
  }
  return(list(psi = psi, a_tau = a_tau, acc = acc, pscale = pscale))
}

# ----- main estimator ------------------------------------------------------

TVPSVD_estim <- function(y, X,
                         nsave = 2000, nburn = 2000, thin = 1,
                         model.setup = NULL,
                         verbose = TRUE){
  # y          ... T x 1 dependent variable
  # X          ... T x K design matrix
  # model.setup .. list with elements
  #   model.type     "TVP" or "TIV"
  #   tvp.type       "WN-SVD" or "RW-SVD"
  #   set.prior.time "g"
  #   pooling        TRUE/FALSE (sparse-mixture clustering of TVPs)
  #   ub.sc          upper-bound scaling of innovation variance
  #   shrnk.type     "ng" (normal-gamma) on constant part
  #   a_tau          NG hyperparameter
  #   sv             TRUE/FALSE stochastic volatility on observation eq.
  #   fast           TRUE -> Bhattacharya-Chakraborty-Mallick sampler
  #   mix.model.setup list with G, e0, sample.e0, perm.sampler, common.mean
  if (is.null(model.setup))
    stop("model.setup must be provided")

  ms <- model.setup
  mms <- model.setup$mix.model.setup
  model.type     <- ms$model.type
  tvp.type       <- ms$tvp.type
  pooling        <- isTRUE(ms$pooling)
  ub.sc          <- ms$ub.sc
  shrnk.type     <- ms$shrnk.type
  a_tau          <- ms$a_tau
  sv             <- isTRUE(ms$sv)
  fast           <- isTRUE(ms$fast)

  G            <- mms$G
  e0           <- mms$e0
  sample.e0    <- isTRUE(mms$sample.e0)
  perm.sampler <- isTRUE(mms$perm.sampler)
  common.mean  <- isTRUE(mms$common.mean)
  e0_j <- e0

  K <- ncol(X)
  T <- nrow(y)
  ntot <- nsave*thin + nburn
  save.set <- seq(thin, nsave*thin, thin) + nburn
  save.len <- length(save.set)
  save.ind <- 0

  # -- construct Z and SVD ----------------------------------------------
  TK <- K*T
  Z <- matrix(0, T, TK)
  if (tvp.type == "WN-SVD"){
    for (i in 1:T) Z[i, (1+K*(i-1)):(K*i)] <- X[i, ]
  } else if (tvp.type == "RW-SVD"){
    for (i in 1:T) Z[i, 1:(i*K)] <- rep(X[i, ], i)
  }
  M <- T
  svd.part <- svd(t(Z), nu = M, nv = M)
  if (tvp.type == "WN-SVD"){
    U <- Matrix(svd.part$u, sparse = TRUE)
    V <- Matrix(svd.part$v, sparse = TRUE)
  } else {
    U <- svd.part$u
    V <- svd.part$v
  }
  lambda.vec <- svd.part$d[1:M]
  lambda.inv.vec <- 1/lambda.vec
  Lambda2.inv <- diag(lambda.inv.vec^2)

  # -- MH set-up --------------------------------------------------------
  mh.sc <- ub.sc/10
  xi.low  <- 1e-10
  xi.high <- ub.sc*T/K^2
  if (verbose)
    message("Upper bound xi_ub = ub.sc * T / K^2 = ", signif(xi.high, 4))

  # -- prior variances of innovations -----------------------------------
  # Generic version: equal prior scaling across coefficients.  When the
  # data are standardised (or generated synthetically with unit-variance
  # regressors) pr.sc = 1 works well.  In the original NKPC application
  # the authors rescale by AR(p) residual variances.
  pr.sc <- rep(1, K)
  xi <- (xi.high + xi.low)/2
  theta <- xi*pr.sc
  if (tvp.type == "WN-SVD" & K > 1){
    psi_draw.tvp <- rep(theta, T)
  } else {
    psi_draw.tvp <- 1/T*rep(theta, T)
  }
  Vmat <- if (K > 1) diag(theta) else xi/T

  if (tvp.type == "WN-SVD"){
    Sigma.beta.inv <- Matrix::Diagonal(TK, 1/psi_draw.tvp)
    U.hat <- U*psi_draw.tvp
    Omega.inv <- 1/diag(t(U.hat) %*% U + Lambda2.inv)
    U.Omega.inv <- t(Omega.inv*t(U.hat))
    Var.1 <- Matrix(psi_draw.tvp*t(sqrt(Omega.inv)*t(U)), sparse = TRUE)
    Var.prod <- -1*tcrossprod(Var.1)
    m2 <- as(Var.prod, "TsparseMatrix")
    if (M == T){
      chg <- m2@i == m2@j
      m2@x[chg] <- m2@x[chg] + psi_draw.tvp
    } else if (M < T){
      diag(m2) <- diag(m2) + psi_draw.tvp
    } else stop("M must be smaller or equal to T!")
  } else if (tvp.type == "RW-SVD"){
    Omega.inv <- lambda.vec^2/(1/xi + lambda.vec^2)
    Omega.tilde.inv <- lambda.vec/(1/xi + lambda.vec^2)
    U.hat   <- t(t(U)*Omega.inv)
    U.tilde <- t(t(U)*Omega.tilde.inv)
  }

  # Initialise b.tilde and prior.part with a rough rolling-OLS estimate so
  # that the sparse-mixture step starts with sensible regime means and
  # does not lock the chain into a high-volatility, badly-fit mode.
  b.tilde <- matrix(0, TK, 1)
  prior.part <- Matrix(c(rep(0, TK)), sparse = TRUE)
  if (T >= 20 && model.type == "TVP"){
    win <- max(round(T/10), 10)
    b_roll <- matrix(0, T, K)
    yvec <- as.numeric(y)
    for (tt in 1:T){
      lo <- max(1, tt - win)
      hi <- min(T, tt + win)
      Xw <- X[lo:hi, , drop = FALSE]
      yw <- yvec[lo:hi]
      bw <- try(solve(crossprod(Xw) + 0.01*diag(K), crossprod(Xw, yw)),
                silent = TRUE)
      if (!inherits(bw, "try-error")) b_roll[tt, ] <- as.numeric(bw)
    }
    # store rolling-OLS values as starting b.tilde (deviation from 0
    # constant part).
    b.tilde <- matrix(as.vector(t(b_roll)), ncol = 1)
    prior.part <- as.numeric(b.tilde)
  }

  # -- priors for constant part -----------------------------------------
  V.pr <- 1/10
  b0 <- matrix(0, K, 1)
  b.draw <- matrix(0, K, 1)
  if (K > 1){
    V0 <- V.pr*diag(K)
    V0inv <- diag(1/diag(V0))
  } else {
    V0 <- V.pr
    V0inv <- 1/V0
  }
  psi_draw <- matrix(0.1, K, 1)

  # -- stochastic volatility --------------------------------------------
  Bvarsigma.h <- 1
  sv_priors <- specify_priors(
    mu = sv_normal(mean = 0, sd = 100),
    phi = sv_beta(shape1 = 25, shape2 = 1.5),
    sigma2 = sv_gamma(shape = 0.5, rate = 1/(2*Bvarsigma.h)),
    nu = sv_infinity(),
    rho = sv_constant(0))
  # Initialise stochastic volatility from preliminary OLS residuals so
  # that the chain does not get attached to a low- or high-volatility mode.
  init_resid <- try({
    lm_init <- lm.fit(X, as.numeric(y))
    as.numeric(lm_init$residuals)
  }, silent = TRUE)
  init_sd <- if (!inherits(init_resid, "try-error"))
    max(sd(init_resid), 1e-2) else 1
  temp_sv <- list(mu = 2*log(init_sd), phi = .9, sigma = .01,
                  nu = Inf, rho = 0, beta = NA, latent0 = 2*log(init_sd))
  sv_latent <- rep(2*log(init_sd), T)
  Sig.t <- rep(init_sd, T)

  # -- mixture prior on TVPs --------------------------------------------
  si <- matrix(0, T, G)
  for (jj in 1:T) si[jj, sample(1:G, 1)] <- 1
  if (pooling){
    c_proposal <- 0.5
    d_kappa <- 0.01
    kappa <- matrix(d_kappa, K, 1)
    zeta <- matrix(1/G, G, 1)
    V0_j <- diag(K)*10^2
    b0_mat <- matrix(0, K, 1)
    b0_j <- as.vector(b0_mat)
    B0 <- diag(K)*1
    B0inv <- solve(B0)
  }
  mu_regimes <- matrix(NA, K, G)
  prior.mat <- matrix(0, K, T)

  # -- storage ----------------------------------------------------------
  b.store <- matrix(NA, save.len, K)
  xi.store <- matrix(NA, save.len, 1)
  mu.store <- array(NA, c(save.len, K, G))
  mu_t.store <- array(NA, c(save.len, T, K))
  st.store <- array(NA, c(save.len, T, G))
  b.tilde.store <- b.full.store <- array(NA, c(save.len, T, K))
  fit.store <- matrix(0, save.len, T)
  G.store <- matrix(0, save.len, 1)
  eht.store <- matrix(NA, save.len, T)
  mh.acc <- pool.acc <- 0

  prog.points <- unique(round(seq(1, ntot, length.out = 11)))

  start <- Sys.time()
  for (irep in 1:ntot){

    # ---- Step 1: constant coefficients -------------------------------
    ytilde <- y - Z %*% b.tilde
    if (K > 1){
      Xnorm <- X/Sig.t
      ynorm <- ytilde/Sig.t
    } else {
      Xnorm <- as.numeric(X/Sig.t)
      ynorm <- as.numeric(ytilde/Sig.t)
    }
    V1 <- try(solve(crossprod(Xnorm) + V0inv), silent = TRUE)
    if (inherits(V1, "try-error")) V1 <- ginv(crossprod(Xnorm) + V0inv)
    b1 <- V1 %*% (crossprod(Xnorm, ynorm) + V0inv %*% b0)
    b.draw <- try(b1 + t(chol(V1)) %*% rnorm(K), silent = TRUE)
    if (inherits(b.draw, "try-error")) b.draw <- mvrnorm(1, b1, V1)
    b.draw <- matrix(b.draw, K, 1)

    # ---- Step 2: time-varying part via SVD ---------------------------
    y.hat <- (y - X %*% b.draw)

    if (model.type == "TVP"){
      # MH on xi
      xi.prop <- exp(rnorm(1, 0, mh.sc))*xi + 1e-20
      theta.prop <- xi.prop*pr.sc
      if (tvp.type == "WN-SVD" & K > 1)
        psi_draw.tvp.prop <- rep(theta.prop, T) else
        psi_draw.tvp.prop <- 1/T*rep(theta.prop, T)

      log.pr.prop <- log(xi)
      log.pr.acc  <- log(xi.prop)
      if (xi.prop > xi.high || xi.prop < xi.low){
        post.prop <- -Inf
      } else {
        post.prop <- sum(dnorm(b.tilde, as.numeric(prior.part),
                               Sig.t*sqrt(psi_draw.tvp.prop), log = TRUE)) +
                     log.pr.prop
      }
      post.acc <- sum(dnorm(b.tilde, as.numeric(prior.part),
                            Sig.t*sqrt(psi_draw.tvp), log = TRUE)) +
                  log.pr.acc

      if ((post.prop - post.acc) > log(runif(1))){
        xi <- xi.prop
        theta <- theta.prop
        psi_draw.tvp <- psi_draw.tvp.prop
        mh.acc <- mh.acc + 1
        if (tvp.type == "WN-SVD"){
          Sigma.beta.inv <- Matrix::Diagonal(TK, 1/psi_draw.tvp)
          U.hat <- U*psi_draw.tvp
          Omega.inv <- 1/diag(t(U.hat) %*% U + Lambda2.inv)
          U.Omega.inv <- t(Omega.inv*t(U.hat))
          Var.1 <- Matrix(psi_draw.tvp*t(sqrt(Omega.inv)*t(U)),
                          sparse = TRUE)
          Var.prod <- -1*tcrossprod(Var.1)
          m2 <- as(Var.prod, "TsparseMatrix")
          if (M == T){
            chg <- m2@i == m2@j
            m2@x[chg] <- m2@x[chg] + psi_draw.tvp
          } else if (M < T){
            diag(m2) <- diag(m2) + psi_draw.tvp
          }
        } else if (tvp.type == "RW-SVD"){
          Omega.inv <- lambda.vec^2/(1/xi + lambda.vec^2)
          Omega.tilde.inv <- lambda.vec/(1/xi + lambda.vec^2)
          U.hat   <- t(t(U)*Omega.inv)
          U.tilde <- t(t(U)*Omega.tilde.inv)
        }
        if (pooling){
          if (K > 1) Vmat <- diag(theta) else Vmat <- xi/T
        }
      }
      if (irep < (0.7*nburn)){
        if (mh.acc/irep < 0.15) mh.sc <- 0.99*mh.sc
        if (mh.acc/irep > 0.25) mh.sc <- 1.01*mh.sc
      }

      scale.mat <- as.vector(t(matrix(Sig.t, T, K)))

      if (tvp.type == "WN-SVD"){
        # Fast Bhattacharya-Chakraborty-Mallick algorithm
        pr.shk <- rnorm(TK, 0, sqrt(psi_draw.tvp))
        lam.shk <- rnorm(T, 0, lambda.inv.vec)
        b.tilde.shk <- pr.shk - U.Omega.inv %*% (t(U) %*% pr.shk + lam.shk)
        mu.K <- m2 %*% (t(Z) %*% y.hat + Sigma.beta.inv %*% prior.part)
      } else if (tvp.type == "RW-SVD"){
        pr.shk <- rnorm(TK, 0, sqrt(xi))
        lam.shk <- rnorm(T, 0, lambda.inv.vec)
        U.vec <- t(U) %*% pr.shk
        b.tilde.shk <- pr.shk - U.hat %*% U.vec - sqrt(xi)*U.hat %*% lam.shk
        v.hat <- t(V) %*% y.hat
        mu.K  <- U.tilde %*% v.hat
      }

      b.tilde.hat <- mu.K + scale.mat*b.tilde.shk
      # generous stability cutoff (original code used 10; this is too tight
      # in the early MCMC phase when stochastic volatility has not yet
      # adapted).  We use 50 instead.
      if (!max(abs(b.tilde.hat)) > 50)
        b.tilde <- matrix(as.vector(t(b.tilde.hat)))
    } else {
      b.tilde <- matrix(0, TK, 1)
    }

    # ---- Step 3: variances -------------------------------------------
    if (tvp.type == "WN-SVD"){
      eta <- y - X %*% b.draw - Z %*% prior.part
      if (K > 1){
        eta <- eta/sqrt(diag(X %*% diag(theta) %*% t(X)) + 1)
      } else {
        eta <- eta/sqrt(theta*X^2 + 1)
      }
    } else {
      eta <- y - X %*% b.draw - Z %*% b.tilde
    }

    if (sv){
      temp_sv <- svsample_fast_cpp(as.numeric(eta), startpara = temp_sv,
                                    startlatent = sv_latent,
                                    priorspec = sv_priors)
      temp_sv[c("mu", "phi", "sigma", "nu", "rho")] <-
        as.list(temp_sv$para[, c("mu", "phi", "sigma", "nu", "rho")])
      sv_latent <- as.numeric(temp_sv$latent)
      sv_latent[sv_latent < -20] <- -20
      Sig.t <- as.numeric(exp(sv_latent/2))
      Sig.t[Sig.t < 1e-8] <- 1e-8
    } else {
      t1 <- (0.01 + T)/2
      S1 <- (0.01 + sum(eta^2))/2
      sigma2 <- as.numeric(1/rgamma(1, t1, as.numeric(S1)))
      Sig.t <- rep(sqrt(sigma2), T)
    }

    # ---- Step 4: shrinkage on constant part --------------------------
    if (shrnk.type == "ng"){
      ng_draw <- get.ng((b.draw - b0[1:K]), psi_draw,
                        a_tau = a_tau, sample_a = FALSE)
      psi_draw <- as.numeric(ng_draw$psi)
      if (K > 1){
        V0 <- diag(as.numeric(psi_draw))
        V0inv <- diag(1/as.numeric(psi_draw))
      } else {
        V0 <- psi_draw
        V0inv <- 1/psi_draw
      }
    }

    bb.tilde <- matrix(b.tilde, K, T)
    if (tvp.type == "RW-SVD") bb.tilde <- t(apply(t(bb.tilde), 2, cumsum))
    bb.full <- matrix(b.draw, K, T) + bb.tilde

    # ---- Steps 5-10: sparse finite mixture pooling -------------------
    if (model.type == "TVP" && pooling){
      Xhat <- Matrix(si/as.numeric(Sig.t), sparse = TRUE)
      ysi <- Matrix((bb.tilde)/as.numeric(Sig.t), sparse = TRUE)
      XsiXsi <- crossprod(Xhat)
      Xsiysi <- crossprod(Xhat, t(ysi))
      Vpost.mix <- 1/(matrix(diag(XsiXsi), K*G, 1) +
                        matrix(diag(B0inv), K*G, 1))
      Amean.mix <- (as.vector(Xsiysi) + rep(B0inv %*% b0_j, G))*(Vpost.mix)
      Vchol <- if (K == 1) matrix(sqrt(Vmat), K, 1) else
                            matrix(sqrt(diag(Vmat)), G*K, 1)
      Amix.draw <- Amean.mix + rnorm(G*K, 0, Vchol*sqrt(diag(Vpost.mix)))
      mu_regimes <- matrix(Amix.draw, K, G)

      ek <- apply(si, 2, sum) + e0_j
      zeta <- bayesm::rdirichlet(ek)
      omega <- zeta/sum(zeta)

      if (sample.e0){
        a_gam <- 0.1; b_gam <- a_gam
        Qmax <- G; const <- Qmax*b_gam
        e0_j <- e0
        le0_p <- log(e0_j) + rnorm(1, 0, c_proposal)
        e0_p <- exp(le0_p)
        omega[omega == 0] <- 10^(-50)
        lalpha1 <- (e0_p - e0_j)*sum(log(omega)) +
                   lgamma(G*e0_p) - lgamma(G*e0_j) -
                   G*(lgamma(e0_p) - lgamma(e0_j)) +
                   (a_gam - 1)*(log(e0_p) - log(e0_j)) -
                   (e0_p - e0_j)*const +
                   log(e0_p) - log(e0_j)
        alpha1 <- min(exp(lalpha1), 1)
        if (runif(1) <= alpha1){
          e0_j <- e0_p; e0 <- e0_j; pool.acc <- pool.acc + 1
        }
        if (irep < (0.5*nburn)){
          if (pool.acc/irep < 0.2) c_proposal <- 0.99*c_proposal
          if (pool.acc/irep > 0.4) c_proposal <- 1.01*c_proposal
        }
      }

      # sample mixture indicators
      si <- matrix(0, T, G)
      for (jj in 1:T){
        phi.vec <- bb.tilde[, jj]
        if (K == 1) var.vec <- Vmat*Sig.t[jj]^2 else
                    var.vec <- diag(Vmat)*Sig.t[jj]^2
        get.probs <- apply(dnorm(phi.vec, mu_regimes,
                                  sqrt(var.vec), log = TRUE), 2, sum) +
                     log(omega)
        pjj <- exp(get.probs - max(get.probs))
        pjj <- pjj/sum(pjj)
        sample.j <- sample(seq(1, G), 1, prob = pjj)
        si[jj, sample.j] <- 1
      }

      # common-mean update for mixture means
      if (common.mean){
        R <- apply(t(bb.tilde), 2, function(x) diff(range(x)))
        scale0 <- rep(1, K)
        for (jj in 1:K)
          scale0[[jj]] <- sum((mu_regimes[jj, ] - b0_j[[jj]])^2)/R[[jj]]^2
        v0 <- v1 <- 10
        aj <- 2*v1
        pk <- v0 - G/2
        Lambda <- matrix(0, K, 1)
        for (ss in 1:K)
          Lambda[ss, 1] <- GIGrvg::rgig(1, pk, scale0[[ss]] + 1e-25, aj)
        if (K > 1){
          B0 <- diag((R^2) * as.numeric(Lambda))
          B0inv <- diag(1/diag(B0))
        } else {
          B0 <- R^2*Lambda
          B0inv <- 1/B0
        }
        b0_j <- MASS::mvrnorm(1, matrix(0, K, 1), 1/G*B0)
      }

      if (perm.sampler){
        perm <- sort(mu_regimes[1, ], index.return = TRUE)$ix
        mu_regimes <- mu_regimes[, perm, drop = FALSE]
        omega <- omega[perm, 1]
        si <- si[, perm]
      }

      for (tt in seq_len(T)){
        sl.ind <- which(si[tt, ] == 1)
        prior.mat[, tt] <- mu_regimes[, sl.ind]
      }
      prior.part <- as.numeric(as.vector(prior.mat))
    }

    if (verbose && irep %in% prog.points)
      message("  iteration ", irep, " / ", ntot,
              " (MH acc. = ", round(mh.acc/irep, 2), ")")

    # ---- store -------------------------------------------------------
    if (irep %in% save.set){
      save.ind <- save.ind + 1
      b.store[save.ind, ] <- as.numeric(b.draw)
      b.tilde.store[save.ind, , ] <- t(bb.tilde)
      b.full.store[save.ind, , ] <- t(bb.full)
      xi.store[save.ind, 1] <- xi
      mu.store[save.ind, , ] <- mu_regimes
      mu_t.store[save.ind, , ] <- t(matrix(prior.part, K, T))
      G.store[save.ind, ] <- sum(apply(si, 2, mean) > 0)
      st.store[save.ind, , ] <- si
      eht.store[save.ind, ] <- Sig.t
      fit.store[save.ind, ] <- as.numeric(X %*% b.draw + Z %*% b.tilde +
                                            rnorm(T, 0, Sig.t))
    }
  }
  end <- Sys.time()
  time.min <- as.numeric(difftime(end, start, units = "mins"))

  return(list(b.store = b.store,
              b.tilde.store = b.tilde.store,
              b.full.store = b.full.store,
              fit.store = fit.store,
              eht.store = eht.store,
              G.store = G.store,
              mu.store = mu.store,
              mu_t.store = mu_t.store,
              st.store = st.store,
              xi.store = xi.store,
              K = K, T = T,
              X = X, y = y,
              tvp.type = tvp.type,
              time.min = time.min))
}
