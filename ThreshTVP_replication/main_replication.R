## main_replication.R
##
## Threshold TVP-VAR with SV -- minimal replication driver.
##
## Estimates the latent-threshold time-varying parameter VAR of
## Huber, Kastner and Feldkircher (JAE, "Should I stay or should I go?
## A latent threshold approach to large-scale mixture innovation models"),
## using the `threshtvp` R package.
##
## Pipeline:
##   1. Load (synthetic by default) data and build VAR regressors.
##   2. Call estimate_tvp() with reasonable defaults.
##   3. Summarise time-varying coefficients (posterior medians/bands).
##   4. Compute time-varying impulse responses via posterior Monte Carlo.
##   5. Save plots to ./figures/ and posterior summaries to ./results/.

## --------------------------------------------------------------------------
## 0. Setup
## --------------------------------------------------------------------------

# Switch the working directory to this script's location if launched
# with Rscript (handy for "Rscript main_replication.R").
.this_file <- tryCatch(
  normalizePath(sys.frame(1)$ofile, mustWork = FALSE),
  error = function(e) NA_character_
)
if (!is.na(.this_file) && nzchar(.this_file) && file.exists(.this_file)) {
  setwd(dirname(.this_file))
}

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

source("auxilliary_functions.R")     # mlag, impulsdtrf, get_companion, ...
source("simulate_data.R")            # simulate_threshvar_data()

suppressPackageStartupMessages({
  library(threshtvp)                 # main estimator
  library(ggplot2)
  library(reshape2)
})

## --------------------------------------------------------------------------
## 1. Data
## --------------------------------------------------------------------------
##
## REPLACE WITH YOUR DATA:
##   - Set use_synthetic <- FALSE
##   - Point my_csv to your file (rows = time, cols = variables, numeric).
##   - Apply your transformations (logs, differences) before assigning Yraw.
##
## The default synthetic block produces a 4-variable VAR(2) with two
## structural breaks in the coefficients (T = 200).

use_synthetic <- TRUE
my_csv        <- "data/my_data.csv"   # placeholder

if (use_synthetic) {
  sim   <- simulate_threshvar_data(T = 200, M = 4, p = 2,
                                   breaks = c(80, 140), seed = 1)
  Yraw  <- sim$Yraw
  truth <- sim
  cat("Loaded synthetic data: T =", nrow(Yraw),
      " M =", ncol(Yraw), "\n")
} else {
  Yraw  <- as.matrix(read.csv(my_csv))
  class(Yraw) <- "numeric"
  truth <- NULL
  cat("Loaded user data from ", my_csv, "\n", sep = "")
}

## --------------------------------------------------------------------------
## 2. Model / MCMC settings
## --------------------------------------------------------------------------
## Defaults can be overridden from the command line, either positionally
##   Rscript main_replication.R 30 20
## or via key=value pairs (recommended for clarity)
##   Rscript main_replication.R save=30 burn=20 horizon=8
.parse_args <- function(defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  out  <- defaults
  pos  <- 1L
  for (a in args) {
    if (grepl("=", a, fixed = TRUE)) {
      kv  <- strsplit(a, "=", fixed = TRUE)[[1]]
      key <- kv[1]; val <- suppressWarnings(as.integer(kv[2]))
      if (key %in% names(out) && !is.na(val)) out[[key]] <- val
    } else {
      val <- suppressWarnings(as.integer(a))
      if (!is.na(val) && pos <= length(out)) {
        out[[names(out)[pos]]] <- val
        pos <- pos + 1L
      }
    }
  }
  out
}
.opts <- .parse_args(list(save = 2000L, burn = 2000L, horizon = 16L, p = 2L))
p        <- .opts$p          # VAR lag length
nsave    <- .opts$save       # posterior draws kept (reduced for replication speed)
nburn    <- .opts$burn       # burn-in
thin     <- 1
horizon  <- .opts$horizon    # IRF horizon
n_cpu    <- 1                # 1 = sequential, >1 = snowfall parallel by equation

# Prior specification (matches the JAE replication script defaults)
prior_btheta <- list(B_1 = 3, B_2 = 0.03, kappa0 = -0.1 / 6)
prior_mu     <- c(0, 10)
pi_low       <- 0.1     # lower threshold pct
pi_high      <- 1.0     # upper threshold pct

## --------------------------------------------------------------------------
## 3. Build regressors and estimate the model
## --------------------------------------------------------------------------
Y    <- Yraw
M    <- ncol(Y)
Xlag <- mlag(Y, p)
Y    <- Y[(p + 1):nrow(Y), , drop = FALSE]
X    <- cbind(Xlag[(p + 1):nrow(Xlag), , drop = FALSE], 1)
Tobs <- nrow(Y)

cat(sprintf("Estimating Threshold TVP-VAR: T=%d  M=%d  p=%d  K=%d\n",
            Tobs, M, p, ncol(X)))
cat(sprintf("MCMC: save=%d burn=%d thin=%g CPU=%d\n",
            nsave, nburn, thin, n_cpu))

t0 <- Sys.time()
model_VAR <- estimate_tvp(
  Y, X,
  save           = nsave,
  burn           = nburn,
  p              = p,
  sv_on          = TRUE,
  thin           = thin,
  priorbtheta    = prior_btheta,
  priormu        = prior_mu,
  h0prior        = "stationary",
  grid.length    = 150,
  thrsh.pct      = pi_low,
  thrsh.pct.high = pi_high,
  TVS            = TRUE,
  CPU            = n_cpu,
  sim.kappa0     = FALSE
)
elapsed <- Sys.time() - t0
cat("Estimation finished in ", format(elapsed), "\n", sep = "")

saveRDS(model_VAR, file = "results/model_VAR.rds")

## --------------------------------------------------------------------------
## 4. Posterior summaries of the time-varying coefficients
## --------------------------------------------------------------------------
##
## $VAR_coeff$A_post has dimension (T, K, M, ndraws):
##   T draws of the time-varying reduced-form coefficient matrix per draw.
A_post <- model_VAR$VAR_coeff$A_post
ndraws <- dim(A_post)[4]

A_median <- apply(A_post, c(1, 2, 3), median)
A_low    <- apply(A_post, c(1, 2, 3), quantile, 0.16)
A_high   <- apply(A_post, c(1, 2, 3), quantile, 0.84)

## Plot a small panel of time-varying coefficients (own-lag-1 diagonals).
# Truth: for the DGP regimes from simulate_threshvar_data the (eq, eq) entry
# of true_A_list[[tt]] is the own-first-lag coefficient at time tt.
has_truth <- !is.null(truth) && !is.null(truth$true_A_list)
coef_df <- do.call(rbind, lapply(seq_len(M), function(eq) {
  median_t <- A_median[, eq, eq]
  if (has_truth) {
    # The estimator drops the first p observations from Yraw, so map tt -> tt+p.
    true_t <- sapply(seq_len(Tobs), function(t)
      truth$true_A_list[[min(length(truth$true_A_list), t + p)]][eq, eq])
  } else true_t <- NA_real_
  data.frame(
    time     = seq_len(Tobs),
    equation = paste0("eq", eq),
    median   = median_t,
    low      = A_low[,    eq, eq],
    high     = A_high[,   eq, eq],
    true     = true_t,
    stringsAsFactors = FALSE
  )
}))

g_coef <- ggplot(coef_df, aes(x = time, y = median)) +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.25) +
  geom_line() +
  facet_wrap(~ equation, scales = "free_y") +
  labs(title = "Time-varying own first-lag coefficient (posterior median, 68% band)",
       subtitle = if (has_truth) "Red dashed: true DGP coefficient path" else NULL,
       x = "time", y = "coefficient") +
  theme_minimal(base_size = 11)
if (has_truth) {
  g_coef <- g_coef + geom_line(aes(y = true), colour = "red3",
                               linetype = "dashed", linewidth = 0.7)
}
ggsave("figures/tvp_coefficients.pdf", g_coef, width = 10, height = 6)
ggsave("figures/tvp_coefficients.png", g_coef, width = 10, height = 6, dpi = 120)

## --------------------------------------------------------------------------
## 5. Time-varying impulse responses (Cholesky, shock to last variable)
## --------------------------------------------------------------------------
shock_var <- M    # default: shock to the LAST variable in the system
IRF.store <- array(NA_real_, c(ndraws, Tobs, M, horizon))
eigs.store <- array(NA_real_, c(ndraws, Tobs))

cat("Computing time-varying IRFs...\n")
for (irep in seq_len(ndraws)) {
  for (tt in seq_len(Tobs)) {
    A_draw   <- A_post[tt, 1:(M * p), , irep]
    SIG_draw <- model_VAR$VAR_coeff$S_post[tt, , , irep]

    comp     <- get_companion(A_draw, c(M, 0, p))
    eigs.store[irep, tt] <- max(abs(Re(eigen(comp$MM)$values)))

    PHI_array <- array(0, c(M, M, p))
    for (k in seq_len(p)) {
      PHI_array[, , k] <- t(A_draw[((k - 1) * M + 1):(k * M), , drop = FALSE])
    }
    cholSIG <- try(t(chol(SIG_draw)), silent = TRUE)
    if (inherits(cholSIG, "try-error")) next
    IRF.t <- impulsdtrf(B = PHI_array, smat = cholSIG, nstep = horizon)
    # response of all variables to shock to shock_var, normalised by its own impact
    IRF.store[irep, tt, , ] <- IRF.t[, shock_var, ] / IRF.t[shock_var, shock_var, 1]
  }
}

# Drop explosive draws
keep <- apply(eigs.store, 1, max, na.rm = TRUE) <= 1
if (any(!keep)) cat("Dropping ", sum(!keep), " explosive draws out of ",
                    ndraws, "\n", sep = "")
IRF.store <- IRF.store[keep, , , , drop = FALSE]

IRF.median <- apply(IRF.store, c(2, 3, 4), median, na.rm = TRUE)
IRF.low    <- apply(IRF.store, c(2, 3, 4), quantile, 0.16, na.rm = TRUE)
IRF.high   <- apply(IRF.store, c(2, 3, 4), quantile, 0.84, na.rm = TRUE)

saveRDS(list(IRF.median = IRF.median, IRF.low = IRF.low, IRF.high = IRF.high,
             shock_var = shock_var),
        file = "results/irfs.rds")

## ---- Plot 1: IRFs at three selected sample dates -------------------------
slice_times <- unique(round(seq(1, Tobs, length.out = 3)))

# True IRF at each slice time: same Cholesky structure as the model. We use
# the (constant) DGP Sigma_chol = 0.4 * diag(M).
true_slice_irf <- function(A_K_M, M, p, H, Sigma_chol) {
  # A_K_M is (M*p+1) x M; rows 1..M = lag-1, rows (M+1)..2M = lag-2, etc.
  PHI <- lapply(seq_len(p), function(l) t(A_K_M[((l-1)*M+1):(l*M), , drop = FALSE]))
  out <- array(0, c(M, M, H))
  out[ , , 1] <- Sigma_chol
  for (h in 2:H) for (l in 1:min(p, h - 1)) {
    out[ , , h] <- out[ , , h] + PHI[[l]] %*% out[ , , h - l]
  }
  out
}
Sigma_chol_true <- 0.4 * diag(M)
irf_slice_df <- do.call(rbind, lapply(slice_times, function(tt) {
  truth_tt <- if (has_truth)
    true_slice_irf(truth$true_A_list[[min(length(truth$true_A_list), tt + p)]],
                   M, p, horizon, Sigma_chol_true)
  else NULL
  # Model normalises by IRF[shock, shock, 1] so we apply the same here.
  if (!is.null(truth_tt)) {
    norm <- truth_tt[shock_var, shock_var, 1]
    truth_tt <- truth_tt / norm
  }
  do.call(rbind, lapply(seq_len(M), function(v) {
    data.frame(
      time_idx = tt,
      var      = paste0("y", v),
      h        = 0:(horizon - 1),
      median   = IRF.median[tt, v, ],
      low      = IRF.low[   tt, v, ],
      high     = IRF.high[  tt, v, ],
      true     = if (!is.null(truth_tt)) truth_tt[v, shock_var, ] else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}))
irf_slice_df$time_idx <- factor(irf_slice_df$time_idx,
                                levels = slice_times,
                                labels = paste0("t=", slice_times))

g_irf <- ggplot(irf_slice_df, aes(x = h, y = median)) +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.25) +
  geom_line() +
  geom_line(aes(y = true), colour = "red3", linetype = "dashed",
            linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey60") +
  facet_grid(var ~ time_idx, scales = "free_y") +
  labs(title = sprintf("Time-varying IRFs (shock to y%d) - red dashed = true",
                       shock_var),
       x = "horizon", y = "response") +
  theme_minimal(base_size = 11)
ggsave("figures/irf_slices.pdf", g_irf, width = 10, height = 8)
ggsave("figures/irf_slices.png", g_irf, width = 10, height = 8, dpi = 120)

## ---- Plot 2: heat-map of the median IRF for one response variable --------
focal <- 1   # show response of y1
heat_df <- melt(IRF.median[, focal, ], varnames = c("time", "horizon"),
                value.name = "irf")
heat_df$horizon <- heat_df$horizon - 1L
g_heat <- ggplot(heat_df, aes(x = time, y = horizon, fill = irf)) +
  geom_tile() +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0) +
  labs(title = sprintf("Time-varying IRF of y%d to shock in y%d (posterior median)",
                       focal, shock_var),
       x = "time", y = "horizon") +
  theme_minimal(base_size = 11)
ggsave("figures/irf_heatmap.pdf", g_heat, width = 10, height = 5)
ggsave("figures/irf_heatmap.png", g_heat, width = 10, height = 5, dpi = 120)

## ---- Plot 3: stochastic volatility per equation --------------------------
sv_df <- do.call(rbind, lapply(seq_len(M), function(eq) {
  Hmat <- model_VAR$posterior[[eq]]$H
  data.frame(
    time     = seq_len(Tobs),
    equation = paste0("eq", eq),
    median   = exp(0.5 * apply(Hmat, 2, median, na.rm = TRUE)),
    low      = exp(0.5 * apply(Hmat, 2, quantile, 0.16, na.rm = TRUE)),
    high     = exp(0.5 * apply(Hmat, 2, quantile, 0.84, na.rm = TRUE)),
    true     = if (has_truth)
                 exp(0.5 * tail(truth$meta$h_path[, eq], Tobs)) else NA_real_,
    stringsAsFactors = FALSE
  )
}))
g_sv <- ggplot(sv_df, aes(x = time, y = median)) +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.25) +
  geom_line() +
  facet_wrap(~ equation, scales = "free_y") +
  labs(title = "Posterior stochastic volatility by equation",
       subtitle = if (has_truth) "Red dashed: true exp(h_t/2) from the DGP" else NULL,
       x = "time", y = "exp(h/2)") +
  theme_minimal(base_size = 11)
if (has_truth) g_sv <- g_sv +
  geom_line(aes(y = true), colour = "red3", linetype = "dashed",
            linewidth = 0.7)
ggsave("figures/stochastic_volatility.pdf", g_sv, width = 10, height = 6)
ggsave("figures/stochastic_volatility.png", g_sv, width = 10, height = 6, dpi = 120)

cat("Done. Figures in ./figures, posterior objects in ./results.\n")
