# simulate_data.R
# ------------------------------------------------------------------------------
# Generate a synthetic stationary VAR(p) data set with mild stochastic
# volatility. The output is saved as data_synth.rda containing:
#   Y      : (T x M) matrix of simulated data
#   labels : variable names (used by main_replication.R for IRF plots)
#
# REPLACE WITH YOUR DATA:
#   To run the replication on real data, comment out the call to
#   simulate_var_sv() and instead read your own Y matrix (T x M numeric),
#   then save it under the same name:
#       Y <- as.matrix(read.csv("my_data.csv"))
#       labels <- colnames(Y)
#       save(Y, labels, file = "data_synth.rda")
# ------------------------------------------------------------------------------

simulate_var_sv <- function(T = 200, M = 4, p = 2, seed = 42,
                            persistence = 0.6, cross = 0.05,
                            sv_phi = 0.95, sv_sigma = 0.15,
                            sv_mu = c(-1.2, -1.0, -0.8, -1.4)) {
  # Simulate a stable VAR(p) with stochastic volatility.
  #   persistence : own-lag coefficient at lag 1
  #   cross       : cross-lag coefficient at lag 1 (off-diagonal)
  #   sv_*        : AR(1) parameters of the log-volatility processes
  set.seed(seed)
  if (length(sv_mu) < M) sv_mu <- rep(sv_mu, length.out = M)

  # VAR coefficient matrices A_1, ..., A_p (size M x M each).
  A_list <- vector("list", p)
  A_list[[1]] <- diag(persistence, M) +
    cross * matrix(runif(M * M, -1, 1), M, M) * (1 - diag(M))
  if (p >= 2) {
    for (l in 2:p) {
      A_list[[l]] <- diag(0.15 / l, M) +
        0.5^l * cross * matrix(runif(M * M, -1, 1), M, M) * (1 - diag(M))
    }
  }

  # Check stability via the companion matrix; if unstable, shrink toward zero.
  build_comp <- function(A_list) {
    p <- length(A_list); M <- nrow(A_list[[1]])
    comp <- matrix(0, M * p, M * p)
    comp[1:M, ] <- do.call(cbind, A_list)
    if (p > 1) comp[(M + 1):(M * p), 1:(M * (p - 1))] <- diag(M * (p - 1))
    comp
  }
  shrink_factor <- 1.0
  repeat {
    comp <- build_comp(lapply(A_list, function(a) a * shrink_factor))
    if (max(abs(eigen(comp, only.values = TRUE)$values)) < 0.98) break
    shrink_factor <- shrink_factor * 0.95
  }
  A_list <- lapply(A_list, function(a) a * shrink_factor)

  # Random unit-lower-triangular contemporaneous mixing matrix (impact effects).
  U <- diag(M)
  U[lower.tri(U)] <- runif(sum(lower.tri(U)), -0.3, 0.3)
  Uinv <- solve(U)

  # Simulate log-volatilities (one AR(1) per equation) and the VAR.
  burn <- 100
  Ttot <- T + burn + p
  h    <- matrix(0, Ttot, M)
  for (j in 1:M) h[1, j] <- sv_mu[j]
  for (t in 2:Ttot) {
    h[t, ] <- sv_mu + sv_phi * (h[t - 1, ] - sv_mu) +
              rnorm(M, 0, sv_sigma)
  }

  Y <- matrix(0, Ttot, M)
  Y[1:p, ] <- matrix(rnorm(p * M, 0, 0.5), p, M)
  for (t in (p + 1):Ttot) {
    mu_t <- rep(0, M)
    for (l in 1:p) mu_t <- mu_t + A_list[[l]] %*% Y[t - l, ]
    eps  <- rnorm(M, 0, exp(0.5 * h[t, ]))
    Y[t, ] <- mu_t + Uinv %*% eps
  }
  Y <- Y[(burn + p + 1):Ttot, , drop = FALSE]

  list(Y = Y, A_true = A_list, U_true = U, h_true = h[(burn + p + 1):Ttot, ])
}

# ---- run and save ------------------------------------------------------------
sim    <- simulate_var_sv(T = 200, M = 4, p = 2)
Y      <- sim$Y
labels <- c("output", "inflation", "interest_rate", "credit")
colnames(Y) <- labels

# Truth saved alongside Y so main_replication.R can overlay it on the figures.
truth  <- list(A_true = sim$A_true,
               U_true = sim$U_true,
               h_true = sim$h_true)

save(Y, labels, truth, file = "data_synth.rda")
cat("simulate_data.R: wrote data_synth.rda (",
    nrow(Y), "x", ncol(Y), ")\n")
