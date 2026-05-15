## =============================================================================
## simulate_data.R
##
## Simulate a synthetic multi-country panel that has the same column-naming
## convention expected by PVAR.IRGA():
##   column name = "<CC>.<var>"   (first two characters = country code).
##
## Defaults: N = 5 countries, M_i = 3 variables per country, T = 150.
## Cross-country dynamic links are weak (controlled by sig.betas).
## Each country has its own structural shocks with small cross-country
## correlation.
##
## REPLACE WITH YOUR DATA: at the bottom we save the simulated matrix to
##   data/synthetic_panel.rda  --  point main_replication.R at your own .rda
##   with the SAME naming convention and the rest of the pipeline will run.
## =============================================================================

simulate_pvar <- function(N = 5, G = 3, T = 150,
                          sparse = 0.7,
                          sig.shocks = 0.4,
                          sig.betas = 0.02,
                          burn = 200,
                          seed = 42) {
  set.seed(seed)

  # Country codes:  c1, c2, ..., cN (two characters).
  cN <- sprintf("c%d", seq_len(N))
  if (any(nchar(cN) != 2)) cN <- sprintf("%02d", seq_len(N))
  varnames <- sprintf("v%d", seq_len(G))

  # Stable own-country block (eigenvalues < 1).
  A1 <- diag(G) * 0.6
  for (a in 1:G) for (b in 1:G) if (a != b)
    A1[a, b] <- 0.1 * (-1)^(a + b)
  V1 <- diag(G) * sig.betas

  Upsilon.true <- diag(G * N) * sig.shocks
  Upsilon.true[lower.tri(Upsilon.true)] <-
    rnorm(G * N * (G * N - 1) / 2, 0, sig.shocks / 10)
  Sigma.true <- crossprod(Upsilon.true)

  A.true <- matrix(0, G * N, G * N)
  ident.country <- t(matrix(seq_len(G * N), G, N))

  Tfull <- T + burn
  Xraw <- matrix(0, Tfull, G * N)
  Xraw[1, ] <- rnorm(G * N, 0, 0.01)

  for (tt in 2:Tfull) {
    shock <- t(chol(Sigma.true)) %*% rnorm(G * N)
    for (c in 1:N) {
      sl.c <- ident.country[c, ]
      sl.r <- setdiff(seq_len(G * N), sl.c)

      A.rest <- matrix(rnorm(G * (G * N - G), 0, sig.betas),
                       G, (G * N - G))
      ind.sparse <- sample(seq_len(ncol(A.rest)),
                           round(sparse * ncol(A.rest)))
      A.rest[, ind.sparse] <- 0
      # Each column of A1 receives an independent G-dim shock from N(0, V1).
      pert <- matrix(t(chol(V1)) %*% matrix(rnorm(G * G), G, G), G, G)
      A1.draw <- A1 + pert

      A.true[sl.c, sl.c] <- A1.draw
      A.true[sl.c, sl.r] <- A.rest

      Xraw[tt, sl.c] <- Xraw[tt - 1, sl.c] %*% A1.draw +
        Xraw[tt - 1, sl.r] %*% t(A.rest) + as.numeric(shock[sl.c])
    }
  }
  Xraw <- Xraw[(burn + 1):Tfull, , drop = FALSE]

  cn <- character(G * N)
  k <- 0
  for (c in seq_len(N))
    for (v in seq_len(G)) {
      k <- k + 1
      cn[k] <- paste0(cN[c], ".", varnames[v])
    }
  colnames(Xraw) <- cn

  # A.mean = block-diagonal lag-1 matrix with A1 on each own-country block and
  # zeros off-diagonal. This is the *center* of the time-varying DGP and is
  # the natural reference for computing true IRFs (the per-period A.true is a
  # noisy realisation around A.mean, with the cross-country block sparse).
  A.mean <- matrix(0, G * N, G * N)
  for (c in seq_len(N)) {
    sl.c <- ident.country[c, ]
    A.mean[sl.c, sl.c] <- A1
  }

  list(Xraw = Xraw, A.true = A.true, A.mean = A.mean, A1 = A1,
       Sigma.true = Sigma.true, ident.country = ident.country,
       cN = cN, varnames = varnames)
}

if (sys.nframe() == 0) {
  out <- simulate_pvar()
  dir.create("data", showWarnings = FALSE)
  Xraw <- out$Xraw
  meta <- list(A.true = out$A.true, A.mean = out$A.mean, A1 = out$A1,
               Sigma.true = out$Sigma.true,
               ident.country = out$ident.country,
               cN = out$cN, varnames = out$varnames)
  save(Xraw, meta, file = "data/synthetic_panel.rda")
  message(sprintf("Wrote data/synthetic_panel.rda  (T = %d, M = %d, N = %d)",
                  nrow(Xraw), ncol(Xraw), length(out$cN)))
}
