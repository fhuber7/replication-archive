# ============================================================================
# tvpbart_estim.R
# ----------------------------------------------------------------------------
# Self-contained estimator for the time-varying parameter VAR with
# regression-tree (BART) driven coefficients of
#
#   Hauzenberger, Huber, Koop & Mitchell,
#   "Bayesian Modeling of Time-Varying Parameter Vector Autoregressions Using
#    Regression Trees", Annals of Applied Statistics.
#
# This is a compact, teaching-oriented version of the sampler shipped in
# estim_file.R / aux_funcs/tvp_softbart_function.R.  It keeps the headline
# idea of the paper -- the coefficients of a VAR drift over the sample as
# nonparametric functions of a low-dimensional state, modelled with Bayesian
# Additive Regression Trees -- but strips out the empirical-data machinery
# (mixed-frequency transforms, EBP instrument, factor SV on the covariance,
# SoftBart C++ back-end, identification variants, ...).
#
# Model (estimated equation by equation in a recursive / triangular order):
#
#   y_{i,t} = sum_k d_{ik,t} * a_{ik,t} + u_{i,t},   u_{i,t} ~ N(0, sigma^2_{i,t})
#   a_{ik,t} = BART_{ik}(z_t)
#
# where d_{i,t} stacks the contemporaneous variables ordered before i
# (recursive identification), the p lags of all variables and a constant, and
# every coefficient a_{ik,t} is its own sum-of-trees function of a
# low-dimensional state z_t (a normalised time trend plus the first lags).
# sigma^2_{i,t} follows a stochastic-volatility process (stochvol), so both the
# conditional mean and the covariance are time-varying.
#
# Relative to the published sampler this routine models each coefficient with
# a separate BART, whereas the paper compresses the K coefficients of each
# equation into a small number Q of common latent BART factors with loadings.
# The per-coefficient version is the Lambda = I special case: less
# parsimonious, but cleanly identified (each tree ensemble's design weight is
# the observed regressor d_{ik,t}, not a sampled loading) and therefore robust
# in a short, self-contained demonstration.  Further simplifications:
#   * Hand-rolled BART (grow/prune Metropolis-Hastings, Gaussian leaves) with
#     observation-specific design weights, instead of SoftBart.
#   * Recursive (Cholesky) structural form only; the maximum-variance and
#     external-instrument identifications of the paper are not reproduced.
#   * Time-t local-linear Cholesky impulse responses (the system is treated
#     as locally linear at each date's coefficient matrix and covariance).
#     The paper computes fully nonlinear generalised impulse responses.
#   * Small MCMC / model defaults so the routine runs in minutes in pure R.
# ============================================================================

suppressPackageStartupMessages({
  library(MASS)
  library(stochvol)
})

# ----------------------------------------------------------------------------
# Generic helpers
# ----------------------------------------------------------------------------

# stacked lags of a (T x N) matrix
mlag <- function(X, lag) {
  X <- as.matrix(X)
  Traw <- nrow(X); N <- ncol(X)
  Xlag <- matrix(0, Traw, lag * N)
  for (ii in seq_len(lag))
    Xlag[(lag + 1):Traw, (N * (ii - 1) + 1):(N * ii)] <- X[(lag + 1 - ii):(Traw - ii), ]
  Xlag
}

.sample1 <- function(x) x[sample.int(length(x), 1L)]

# ----------------------------------------------------------------------------
# Minimal weighted BART
# ----------------------------------------------------------------------------
# A tree is a set of growable node vectors.  Each observation t enters the
# leaf likelihood through a design weight w_t, i.e. the leaf value mu explains
# w_t * mu of the (partial) residual r_t with observation variance s2_t.  The
# weighted sufficient statistics A = sum w^2 / s2 and B = sum w r / s2 never
# divide by an individual w_t, so the multiplicative coupling between a
# regressor and its time-varying coefficient is stable even when some weights
# are close to zero.
# ----------------------------------------------------------------------------

.bt_newtree <- function() {
  list(var = NA_integer_, cut = NA_real_, left = 0L, right = 0L,
       parent = 0L, mu = 0, leaf = TRUE, active = TRUE)
}

.bt_leaves  <- function(tr) which(tr$active & tr$leaf)
.bt_nleaves <- function(tr) sum(tr$active & tr$leaf)

# internal nodes whose two children are both leaves (the prunable nodes)
.bt_W2 <- function(tr) {
  internal <- which(tr$active & !tr$leaf)
  if (!length(internal)) return(integer(0))
  internal[tr$leaf[tr$left[internal]] & tr$leaf[tr$right[internal]]]
}

.bt_depth <- function(tr, id) {
  d <- 0L
  while (tr$parent[id] > 0L) { id <- tr$parent[id]; d <- d + 1L }
  d
}

# leaf membership of every row of Z (vectorised traversal)
.bt_assign <- function(tr, Z) {
  n <- nrow(Z)
  node <- rep(1L, n)
  repeat {
    internal <- which(!tr$leaf[node])
    if (!length(internal)) break
    nd <- node[internal]
    zv <- Z[cbind(internal, tr$var[nd])]
    node[internal] <- ifelse(zv < tr$cut[nd], tr$left[nd], tr$right[nd])
  }
  node
}

.bt_predict <- function(tr, Z) tr$mu[.bt_assign(tr, Z)]

# weighted leaf statistics for a set of observation indices
.AB <- function(obs, w, r, s2) {
  if (!length(obs)) return(c(A = 0, B = 0))
  c(A = sum(w[obs]^2 / s2[obs]), B = sum(w[obs] * r[obs] / s2[obs]))
}

# log marginal likelihood of one leaf (shared data constant dropped)
.bt_ll <- function(AB, tau2) {
  A <- AB[1]; B <- AB[2]
  -0.5 * log(2 * pi) - 0.5 * log(1 + tau2 * A) + 0.5 * tau2 * B^2 / (1 + tau2 * A)
}

.do_grow <- function(tr, L, j, cut) {
  idL <- length(tr$var) + 1L; idR <- idL + 1L
  tr$leaf[L] <- FALSE; tr$var[L] <- j; tr$cut[L] <- cut
  tr$left[L] <- idL;   tr$right[L] <- idR
  for (id in c(idL, idR)) {
    tr$var[id] <- NA_integer_; tr$cut[id] <- NA_real_
    tr$left[id] <- 0L; tr$right[id] <- 0L
    tr$parent[id] <- L; tr$mu[id] <- 0; tr$leaf[id] <- TRUE; tr$active[id] <- TRUE
  }
  tr
}

.do_prune <- function(tr, N) {
  cL <- tr$left[N]; cR <- tr$right[N]
  tr$active[c(cL, cR)] <- FALSE; tr$leaf[c(cL, cR)] <- FALSE
  tr$leaf[N] <- TRUE; tr$var[N] <- NA_integer_; tr$cut[N] <- NA_real_
  tr$left[N] <- 0L; tr$right[N] <- 0L
  tr
}

.bt_grow <- function(tr, Z, r, w, s2, hyp, assign) {
  leaves <- .bt_leaves(tr); b <- length(leaves)
  L <- .sample1(leaves)
  obsL <- which(assign == L)
  if (length(obsL) < 2L) return(list(tree = tr, ok = FALSE))
  j <- sample.int(ncol(Z), 1L)
  uv <- sort(unique(Z[obsL, j]))
  if (length(uv) < 2L) return(list(tree = tr, ok = FALSE))
  cuts <- uv[-1]; cut <- .sample1(cuts); ncut <- length(cuts)
  padj <- sum(apply(Z[obsL, , drop = FALSE], 2, function(z) length(unique(z)) >= 2))
  leftobs  <- obsL[Z[obsL, j] <  cut]
  rightobs <- obsL[Z[obsL, j] >= cut]

  llr <- .bt_ll(.AB(leftobs, w, r, s2),  hyp$tau2) +
         .bt_ll(.AB(rightobs, w, r, s2), hyp$tau2) -
         .bt_ll(.AB(obsL, w, r, s2),     hyp$tau2)

  dL  <- .bt_depth(tr, L)
  psL <- hyp$alpha * (1 + dL)^(-hyp$beta)
  psC <- hyp$alpha * (2 + dL)^(-hyp$beta)
  logprior <- log(psL) + 2 * log(1 - psC) - log(1 - psL)

  sib <- if (tr$parent[L] > 0L) { P <- tr$parent[L]; ifelse(tr$left[P] == L, tr$right[P], tr$left[P]) } else 0L
  W2after  <- length(.bt_W2(tr)) - (sib > 0L && tr$leaf[sib]) + 1L
  logtrans <- log(b) + log(padj) + log(ncut) - log(W2after)

  if (log(runif(1)) < llr + logprior + logtrans)
    return(list(tree = .do_grow(tr, L, j, cut), ok = TRUE))
  list(tree = tr, ok = FALSE)
}

.bt_prune <- function(tr, Z, r, w, s2, hyp, assign) {
  W2 <- .bt_W2(tr); w2 <- length(W2)
  if (!w2) return(list(tree = tr, ok = FALSE))
  N <- .sample1(W2); cL <- tr$left[N]; cR <- tr$right[N]
  obsL <- which(assign == cL); obsR <- which(assign == cR); obsP <- c(obsL, obsR)

  llr <- .bt_ll(.AB(obsP, w, r, s2), hyp$tau2) -
         .bt_ll(.AB(obsL, w, r, s2), hyp$tau2) -
         .bt_ll(.AB(obsR, w, r, s2), hyp$tau2)

  dN  <- .bt_depth(tr, N)
  psN <- hyp$alpha * (1 + dN)^(-hyp$beta)
  psC <- hyp$alpha * (2 + dN)^(-hyp$beta)
  logprior <- log(1 - psN) - log(psN) - 2 * log(1 - psC)

  padjN <- sum(apply(Z[obsP, , drop = FALSE], 2, function(z) length(unique(z)) >= 2))
  uv    <- sort(unique(Z[obsP, tr$var[N]])); ncutN <- length(uv) - 1L
  bafter <- .bt_nleaves(tr) - 1L
  logtrans <- log(w2) - log(bafter) - log(padjN) - log(ncutN)

  if (log(runif(1)) < llr + logprior + logtrans)
    return(list(tree = .do_prune(tr, N), ok = TRUE))
  list(tree = tr, ok = FALSE)
}

.bt_drawleaves <- function(tr, Z, r, w, s2, tau2) {
  assign <- .bt_assign(tr, Z)
  for (L in .bt_leaves(tr)) {
    ab <- .AB(which(assign == L), w, r, s2)
    prec <- ab[1] + 1 / tau2
    tr$mu[L] <- rnorm(1, ab[2] / prec, sqrt(1 / prec))
  }
  tr
}

.bart_init <- function(num.trees, T) {
  list(trees    = replicate(num.trees, .bt_newtree(), simplify = FALSE),
       tree_fit = replicate(num.trees, rep(0, T), simplify = FALSE),
       fit      = rep(0, T))
}

# one Gibbs sweep over all trees of one ensemble; returns the updated object.
# Model: r_t = w_t * sum_m T_m(z_t) + N(0, s2_t).
.bart_sweep <- function(bt, Z, r, w, s2, hyp) {
  for (m in seq_along(bt$trees)) {
    gm <- bt$tree_fit[[m]]
    Rm <- r - w * (bt$fit - gm)                    # partial residual
    tr <- bt$trees[[m]]
    assign <- .bt_assign(tr, Z)
    res <- if (runif(1) < 0.5) .bt_grow(tr, Z, Rm, w, s2, hyp, assign)
           else                .bt_prune(tr, Z, Rm, w, s2, hyp, assign)
    tr <- .bt_drawleaves(res$tree, Z, Rm, w, s2, hyp$tau2)
    newgm <- .bt_predict(tr, Z)
    bt$fit <- bt$fit + (newgm - gm)
    bt$tree_fit[[m]] <- newgm
    bt$trees[[m]] <- tr
  }
  bt
}

# ----------------------------------------------------------------------------
# Impulse responses: local-linear Cholesky IRF from a time-t coefficient set
# ----------------------------------------------------------------------------
# avec[[i]] is the coefficient vector of equation i at the requested date:
#   [ contemporaneous y_1..y_{i-1} | lag coefs (Mp, lag-major) | const ].
# Map the structural (recursive) coefficients to the reduced form.
.struct_to_rf <- function(avec, M, p, cons) {
  A0  <- diag(M); cst <- numeric(M)
  B   <- vector("list", p); for (l in seq_len(p)) B[[l]] <- matrix(0, M, M)
  for (i in seq_len(M)) {
    a <- avec[[i]]; pos <- 0L
    if (i > 1L) { A0[i, 1:(i - 1)] <- -a[1:(i - 1)]; pos <- i - 1L }
    for (l in seq_len(p)) B[[l]][i, ] <- a[pos + (l - 1L) * M + 1:M]
    if (cons) cst[i] <- a[pos + M * p + 1L]
  }
  A0inv <- solve(A0)
  list(A0inv = A0inv,
       Arf   = lapply(B, function(Bl) A0inv %*% Bl),
       crf   = as.numeric(A0inv %*% cst))
}

# Local-linear Cholesky IRF from a time-t coefficient set and structural sds.
# irf[, j, h] = response of all variables at horizon h-1 to a shock in var j.
.irf_from_coef <- function(avec, sig, M, p, nhor, cons) {
  rf  <- .struct_to_rf(avec, M, p, cons)
  Arf <- rf$Arf
  P   <- rf$A0inv %*% diag(sig, M)                       # Cholesky impact
  Mp  <- M * p
  Comp <- matrix(0, Mp, Mp)
  for (l in seq_len(p)) Comp[1:M, ((l - 1) * M + 1):(l * M)] <- Arf[[l]]
  if (p > 1) Comp[(M + 1):Mp, 1:(M * (p - 1))] <- diag(M * (p - 1))
  maxeig <- max(Mod(eigen(Comp, only.values = TRUE)$values))

  J <- cbind(diag(M), matrix(0, M, Mp - M))
  irf <- array(0, c(M, M, nhor)); Fh <- diag(Mp)
  for (h in seq_len(nhor)) { irf[, , h] <- (J %*% Fh %*% t(J)) %*% P; Fh <- Fh %*% Comp }
  list(irf = irf, maxeig = maxeig, Arf = Arf, crf = rf$crf)
}

# ----------------------------------------------------------------------------
# Main estimator
# ----------------------------------------------------------------------------

#' Fit the TVP-BART VAR
#'
#' @param Yraw      (T x M) data matrix (column names become variable names).
#' @param p         number of lags.
#' @param nsave,nburn,thin  MCMC budget.
#' @param num.trees number of trees in each coefficient's BART ensemble.
#' @param cons      include an intercept.
#' @param sv        stochastic volatility (TRUE) or homoskedastic errors.
#' @param c.bart    prior scale of the (mean-zero) time variation in each
#'                  coefficient -- i.e. how far coefficients may drift from
#'                  their constant baseline.  Small values shrink toward a
#'                  constant-coefficient VAR.
#' @param state     BART state z_t: "trend" (normalised time only) or
#'                  "trend+lags" (time plus the scaled first lags; more
#'                  flexible but prone to over-fitting in short samples).
#' @param shk.var   index/name of the shocked variable for the IRFs.
#' @param shk.sc    shock size in structural standard deviations.
#' @param nhor      IRF horizon (number of steps incl. impact).
#' @param n.irf.dates number of evenly spaced dates at which time-varying IRFs
#'                  are stored.
#' @param eig.crit  draws/dates whose companion spectral radius exceeds this
#'                  are flagged unstable and set to NA in the IRFs.
#' @param verbose   print progress.
tvpbart_estim <- function(Yraw,
                          p           = 1,
                          nsave       = 1000,
                          nburn       = 1000,
                          thin        = 1,
                          num.trees   = 15,
                          cons        = TRUE,
                          sv          = TRUE,
                          c.bart      = 0.5,
                          state       = "trend",
                          shk.var     = 1,
                          shk.sc      = 1,
                          nhor        = 17,
                          n.irf.dates = 40,
                          eig.crit    = 1.05,
                          verbose     = TRUE) {

  start <- Sys.time()
  Y <- as.matrix(Yraw)
  if (is.null(colnames(Y))) colnames(Y) <- paste0("y", seq_len(ncol(Y)))
  var.names <- colnames(Y)
  M <- ncol(Y)
  if (is.character(shk.var)) shk.var <- match(shk.var, var.names)

  # ---- design ---------------------------------------------------------------
  Xlag <- mlag(Y, p)
  y  <- Y[(p + 1):nrow(Y), , drop = FALSE]
  XL <- Xlag[(p + 1):nrow(Xlag), , drop = FALSE]            # T x Mp (lag-major)
  T  <- nrow(y)
  Klag <- M * p

  # per-equation design D_i = [contemp y_1..y_{i-1} | lags | const]
  D <- vector("list", M); Ki <- integer(M)
  for (i in seq_len(M)) {
    Di <- cbind(if (i > 1) y[, 1:(i - 1), drop = FALSE] else NULL, XL)
    if (cons) Di <- cbind(Di, 1)
    D[[i]] <- Di; Ki[i] <- ncol(Di)
  }

  # ---- BART state z_t -------------------------------------------------------
  sc01 <- function(x) { rg <- range(x); if (diff(rg) < 1e-12) rep(0.5, length(x)) else (x - rg[1]) / diff(rg) }
  Z <- matrix((seq_len(T) - 1) / (T - 1), T, 1, dimnames = list(NULL, "time"))
  if (state == "trend+lags")
    Z <- cbind(Z, apply(XL[, 1:M, drop = FALSE], 2, sc01))
  colnames(Z)[-1] <- paste0(var.names, ".l1")

  hyp <- list(alpha = 0.95, beta = 2, tau2 = c.bart^2 / num.trees)

  # ---- state: constant baseline + centered BART deviation per coefficient ----
  # a_{ik,t} = abar_{ik} + h_{ik}(z_t).  The baseline (a conjugate Gaussian
  # draw) carries the level; each BART models a shrunk, mean-zero deviation, so
  # the time variation is a departure from a constant-coefficient VAR.
  abar  <- vector("list", M)   # K_i constant baseline (OLS warm start)
  barts <- vector("list", M)   # K_i BART ensembles (deviations)
  hdev  <- vector("list", M)   # T x K_i centered deviations
  for (i in seq_len(M)) {
    Di <- D[[i]]
    abar[[i]]  <- as.numeric(MASS::ginv(crossprod(Di)) %*% crossprod(Di, y[, i]))
    barts[[i]] <- lapply(seq_len(Ki[i]), function(k) .bart_init(num.trees, T))
    hdev[[i]]  <- matrix(0, T, Ki[i])
  }
  abar.prec <- 1 / 100         # diffuse normal prior precision on the baseline

  # ---- stochastic volatility ------------------------------------------------
  s2 <- matrix(0, M, T)
  if (sv) {
    sv_priors <- specify_priors(mu = sv_normal(0, 100), phi = sv_beta(25, 1.5),
                                sigma2 = sv_gamma(0.5, 5), nu = sv_infinity(),
                                rho = sv_constant(0))
    h_sv <- vector("list", M); ht <- matrix(0, M, T)
    for (i in seq_len(M)) {
      v0 <- var(y[, i] - D[[i]] %*% abar[[i]])
      h_sv[[i]] <- list(mu = log(v0), phi = 0.95, sigma = 0.1, nu = Inf, rho = 0,
                        beta = NA, latent0 = log(v0))
      ht[i, ] <- log(v0); s2[i, ] <- v0
    }
  } else {
    a.sig <- b.sig <- 0.01
    for (i in seq_len(M)) s2[i, ] <- var(y[, i] - D[[i]] %*% abar[[i]])
  }

  # ---- IRF storage grid -----------------------------------------------------
  irf.dates <- unique(round(seq(1, T, length.out = min(n.irf.dates, T))))
  nD <- length(irf.dates)

  # ---- storage --------------------------------------------------------------
  save.len  <- floor(nsave / thin)
  A.store   <- array(NA_real_, c(save.len, T, M, Klag + cons))  # reduced-form coefs
  s2.store  <- array(NA_real_, c(save.len, M, T))
  irf.store <- array(NA_real_, c(save.len, nD, M, nhor))        # responses to shk.var
  irf.avg   <- array(NA_real_, c(save.len, M, nhor))            # sample-averaged
  dimnames(irf.store) <- list(NULL, paste0("t", irf.dates), var.names, 0:(nhor - 1))

  ntot <- nburn + nsave
  for (irep in seq_len(ntot)) {

    for (i in seq_len(M)) {
      Di <- D[[i]]; yi <- y[, i]; s2i <- s2[i, ]; ki <- Ki[i]
      ab <- abar[[i]]; Hd <- hdev[[i]]
      base <- as.numeric(Di %*% ab)
      # ---- (a) mean-zero BART deviation of each coefficient --------------
      for (k in seq_len(ki)) {
        r <- yi - base - rowSums((Di * Hd)[, -k, drop = FALSE])   # = Di[,k]*h_k + u
        barts[[i]][[k]] <- .bart_sweep(barts[[i]][[k]], Z, r, Di[, k], s2i, hyp)
        h <- barts[[i]][[k]]$fit
        Hd[, k] <- h - mean(h)                              # baseline owns the level
      }
      hdev[[i]] <- Hd
      # ---- (b) constant baseline via conjugate normal -------------------
      Gw   <- Di / sqrt(s2i)
      rw   <- (yi - rowSums(Di * Hd)) / sqrt(s2i)
      Vp   <- chol2inv(chol(crossprod(Gw) + diag(abar.prec, ki)))
      ab   <- as.numeric(Vp %*% crossprod(Gw, rw) + t(chol(Vp)) %*% rnorm(ki))
      abar[[i]] <- ab
      # ---- (c) structural residual -> stochastic volatility -------------
      u <- yi - as.numeric(Di %*% ab) - rowSums(Di * Hd)
      if (sv) {
        sd <- svsample_fast_cpp(u, startpara = h_sv[[i]], startlatent = ht[i, ],
                                priorspec = sv_priors)
        sd[c("mu","phi","sigma","nu","rho")] <- as.list(sd$para[, c("mu","phi","sigma","nu","rho")])
        h_sv[[i]] <- sd
        hh <- as.numeric(sd$latent); hh[hh < -20] <- -20
        ht[i, ] <- hh; s2[i, ] <- exp(hh)
      } else {
        s2[i, ] <- 1 / rgamma(1, a.sig + T / 2, b.sig + sum(u^2) / 2)
      }
    }

    # ---- store ---------------------------------------------------------------
    if (irep > nburn && ((irep - nburn) %% thin == 0)) {
      sidx <- (irep - nburn) / thin
      s2.store[sidx, , ] <- s2

      for (tt in seq_len(T)) {                              # reduced-form coefs, all t
        avec <- lapply(seq_len(M), function(i) abar[[i]] + hdev[[i]][tt, ])
        rf <- .struct_to_rf(avec, M, p, cons)
        for (l in seq_len(p))
          A.store[sidx, tt, , ((l - 1) * M + 1):(l * M)] <- rf$Arf[[l]]
        if (cons) A.store[sidx, tt, , Klag + 1] <- rf$crf
      }

      avg <- array(0, c(M, M, nhor)); navg <- 0                # time-varying IRFs
      for (dd in seq_len(nD)) {
        tt   <- irf.dates[dd]
        avec <- lapply(seq_len(M), function(i) abar[[i]] + hdev[[i]][tt, ])
        out  <- .irf_from_coef(avec, sqrt(s2[, tt]), M, p, nhor, cons)
        if (out$maxeig <= eig.crit) {
          irf.store[sidx, dd, , ] <- out$irf[, shk.var, ] * shk.sc
          avg <- avg + out$irf * shk.sc; navg <- navg + 1
        }
      }
      if (navg > 0) irf.avg[sidx, , ] <- (avg / navg)[, shk.var, ]
    }

    if (verbose && irep %% 50 == 0) cat(sprintf("  MCMC %5d / %d\n", irep, ntot))
  }

  time.min <- as.numeric(difftime(Sys.time(), start, units = "mins"))

  list(A.store   = A.store,      # nsave x T x M x (Mp+cons)  reduced-form coefs
       s2.store  = s2.store,     # nsave x M x T              SV variances
       irf.store = irf.store,    # nsave x nD x M x nhor       IRF to shk.var by date
       irf.avg   = irf.avg,      # nsave x M x nhor            sample-averaged IRF
       irf.dates = irf.dates,
       y = y, Z = Z, var.names = var.names,
       dims = list(M = M, T = T, p = p, Klag = Klag, cons = cons,
                   nhor = nhor, shk.var = shk.var, shk.sc = shk.sc),
       time.min = time.min)
}
