## simulate_data.R
##
## Generates a synthetic stationary VAR(p) data set with two structural
## breaks in selected coefficients to illustrate the latent-threshold
## TVP-VAR model implemented in the `threshtvp` package.
##
## REPLACE WITH YOUR DATA: in main_replication.R, point the data-loading
## block to your own time-series CSV (variables in columns, time in rows).
## The script below is only for demonstration / sanity-checking.
##
## Output: an R data set with
##   Yraw : (T x M) numeric matrix of simulated observations
##   true_A_list : list of length T with the (Mp+1) x M VAR coefficient
##                 matrix that generated each period's observation
##   meta : list with M, p, T, break points used for the simulation

simulate_threshvar_data <- function(T = 200,
                                    M = 4,
                                    p = 2,
                                    breaks = c(80, 140),
                                    seed = 1) {
  set.seed(seed)

  K <- M * p + 1                # number of regressors per equation (incl. const)

  # Helper: build a stable VAR companion matrix
  make_stable_coefs <- function(M, p, scale = 0.15) {
    repeat {
      A <- matrix(rnorm(K * M, sd = scale), nrow = K, ncol = M)
      # add some persistence on the diagonal of the first own-lag block
      diag(A[1:M, ]) <- diag(A[1:M, ]) + runif(M, 0.55, 0.75)
      # very small intercept
      A[K, ] <- rnorm(M, sd = 0.02)
      # companion form (without constant) -- check eigenvalues
      MM <- rbind(t(A[1:(M * p), , drop = FALSE]),
                  cbind(diag((p - 1) * M),
                        matrix(0, (p - 1) * M, M)))
      if (max(Mod(eigen(MM, only.values = TRUE)$values)) < 0.96) {
        return(A)
      }
    }
  }

  # One regime per (break + final) segment, supports arbitrary # breaks >= 1
  breaks <- sort(unique(as.integer(breaks)))
  segs <- c(0L, breaks, T)
  n_reg <- length(segs) - 1L
  regimes <- lapply(seq_len(n_reg), function(j) {
    make_stable_coefs(M, p, scale = c(0.10, 0.18, 0.12, 0.20, 0.08)[
      ((j - 1L) %% 5L) + 1L])
  })

  # Per-period coefficient matrix
  true_A_list <- vector("list", T)
  for (tt in seq_len(T)) {
    reg <- findInterval(tt, breaks, left.open = TRUE) + 1L
    true_A_list[[tt]] <- regimes[[reg]]
  }

  # Modest stochastic volatility on the innovations
  Sigma_chol <- 0.4 * diag(M)
  h_path <- matrix(0, T, M)
  for (j in seq_len(M)) {
    h_path[, j] <- as.numeric(arima.sim(list(ar = 0.95), n = T, sd = 0.10))
  }

  # Simulate the data
  Yraw <- matrix(0, T + p, M)
  Yraw[1:p, ] <- matrix(rnorm(p * M, sd = 0.3), p, M)
  for (tt in seq_len(T)) {
    lags <- as.numeric(t(Yraw[(p + tt - 1):tt, ]))   # stack [y_{t-1}; y_{t-2}; ...]
    x_t <- c(lags, 1)                                # add constant last
    mu_t <- as.numeric(t(true_A_list[[tt]]) %*% x_t)
    eps_t <- diag(exp(0.5 * h_path[tt, ])) %*% Sigma_chol %*% rnorm(M)
    Yraw[p + tt, ] <- mu_t + as.numeric(eps_t)
  }
  Yraw <- Yraw[(p + 1):(T + p), , drop = FALSE]
  colnames(Yraw) <- paste0("y", seq_len(M))

  list(Yraw = Yraw,
       true_A_list = true_A_list,
       meta = list(M = M, p = p, T = T, breaks = breaks, h_path = h_path))
}

## ---------------------------------------------------------------------------
## When sourced standalone, also write a CSV copy of the simulated data
## so that `main_replication.R` can read it like a real data set.
## ---------------------------------------------------------------------------
if (sys.nframe() == 0) {
  sim <- simulate_threshvar_data(T = 200, M = 4, p = 2,
                                 breaks = c(80, 140), seed = 1)
  out_dir <- "data"
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(sim$Yraw, file = file.path(out_dir, "simulated_data.csv"),
            row.names = FALSE)
  saveRDS(sim, file = file.path(out_dir, "simulated_data.rds"))
  cat("Wrote simulated data to ", file.path(out_dir, "simulated_data.csv"),
      " and .rds\n", sep = "")
}
