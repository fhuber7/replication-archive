## =============================================================================
## main_replication.R
##
## Headline replication driver for the panel VAR with Integrated Rotated
## Gaussian Approximation (IRGA).  Loads (synthetic) data, runs the IRGA
## estimator, computes Cholesky impulse responses to a shock to the focal
## variable in country 1, and saves figures to ./figures.
##
## Replace data/synthetic_panel.rda with your own .rda containing a matrix
## Xraw whose columns follow the "<CC>.<var>" naming convention.
## =============================================================================

source("pvar_irga_estim.R")

## ---- Configuration ---------------------------------------------------------
N_country <- 5        # number of countries used in this run
laglen    <- 2        # VAR lag length
nsave     <- 1000     # MCMC draws after burn-in (per equation)
nburn     <- 500
thin.by   <- 1
nhor      <- 24       # IRF horizon
shrink.type <- "HS"   # "HS" Horseshoe or "MN" Minnesota
focal_var <- 1        # which own-country variable receives the shock
seed      <- 1729

set.seed(seed)
dir.create("figures", showWarnings = FALSE)
dir.create("data",    showWarnings = FALSE)

## ---- Load (or simulate) data ----------------------------------------------
dat_path <- "data/synthetic_panel.rda"
source("simulate_data.R")  # defines simulate_pvar()
if (!file.exists(dat_path)) {
  message("No synthetic data found, generating...")
  out  <- simulate_pvar()
  Xraw <- out$Xraw
  meta <- list(A.true = out$A.true, A.mean = out$A.mean, A1 = out$A1,
               Sigma.true = out$Sigma.true,
               ident.country = out$ident.country,
               cN = out$cN, varnames = out$varnames)
  save(Xraw, meta, file = dat_path)
}
load(dat_path)  # provides Xraw + meta
stopifnot(exists("Xraw"))

# Restrict to the first N_country countries (stable two-character prefix).
cN_all <- substr(colnames(Xraw), 1, 2)
cN_uq  <- unique(cN_all)
keep_cN <- cN_uq[seq_len(min(N_country, length(cN_uq)))]
Xraw <- Xraw[, cN_all %in% keep_cN, drop = FALSE]
message(sprintf("Data: T = %d, M = %d, N = %d", nrow(Xraw), ncol(Xraw),
                length(keep_cN)))

# Standardise (the original paper estimates on standardised stationary data).
Yraw <- apply(Xraw, 2, function(x) (x - mean(x)) / sd(x))

## ---- Estimate -------------------------------------------------------------
fit <- PVAR.IRGA(Yraw = Yraw, p = laglen,
                 nsave = nsave, nburn = nburn, thin.by = thin.by,
                 shrink.type = shrink.type, estimator = "VAMP",
                 quiet = FALSE)

A_post     <- fit$A                # nthinned x K x M
SIGMA_post <- fit$SIGMA            # nthinned x M x M
Y          <- fit$Y
cN         <- fit$countries
M          <- ncol(Y)
K          <- dim(A_post)[2]
nthinned   <- dim(A_post)[1]

## ---- Impulse responses ----------------------------------------------------
message("Computing Cholesky impulse responses...")
# Use the posterior mean of SIGMA (it is constant across draws in this set-up).
Sigma_mean <- apply(SIGMA_post, c(2, 3), mean)

# Compute IRFs from each posterior draw.
irf_arr <- irf_cholesky(A_post = A_post, Sigma = Sigma_mean,
                        p = laglen, nstep = nhor)
# irf_arr: ndraw x M (response) x M (shock) x nhor

# Index of the focal shock = focal_var-th variable in country 1.
country1 <- cN[1]
vars_c1  <- which(grepl(paste0("^", country1, "\\."), colnames(Y)))
shock_idx <- vars_c1[focal_var]
message(sprintf("Focal shock: %s", colnames(Y)[shock_idx]))

# Aggregate own-country IRFs of each country's first variable to the shock.
own_irfs <- array(NA, c(length(cN), nhor, 3))   # country x horizon x q(.05,.50,.95)
dimnames(own_irfs) <- list(cN, paste0("h", seq_len(nhor)),
                           c("q05", "q50", "q95"))
for (k in seq_along(cN)) {
  vars_k  <- which(grepl(paste0("^", cN[k], "\\."), colnames(Y)))
  resp_id <- vars_k[focal_var]
  resp <- irf_arr[, resp_id, shock_idx, ]
  own_irfs[k, , 1] <- apply(resp, 2, quantile, 0.05)
  own_irfs[k, , 2] <- apply(resp, 2, quantile, 0.50)
  own_irfs[k, , 3] <- apply(resp, 2, quantile, 0.95)
}

## ---- TRUE own-country IRFs from the DGP -----------------------------------
# The DGP has lag 1 only; the estimator uses laglen lags, so pad with zeros.
# Cholesky impact = chol(Sigma.true). Then IRF[ , , h] = A^{h-1} * impact,
# computed via the standard companion-form recursion.
true_own <- NULL
if (exists("meta") && !is.null(meta$A.mean)) {
  # Restrict A.mean / Sigma.true to the columns retained by N_country.
  keep_idx <- which(substr(colnames(Xraw_full <- {
    e <- new.env(); load(dat_path, envir = e); e$Xraw
  }), 1, 2) %in% keep_cN)
  A.mean_red    <- meta$A.mean[keep_idx, keep_idx]
  Sigma.true_red <- meta$Sigma.true[keep_idx, keep_idx]
  M_red <- length(keep_idx)
  irf_true_full <- array(0, c(M_red, M_red, nhor))
  irf_true_full[ , , 1] <- t(chol(Sigma.true_red))
  for (h in 2:nhor) {
    irf_true_full[ , , h] <- A.mean_red %*% irf_true_full[ , , h - 1]
  }
  # Slice to the own-country first-variable response to the focal shock.
  vars_c1_true <- which(grepl(paste0("^", country1, "\\."),
                              colnames(Xraw_full)[keep_idx]))
  shock_idx_t  <- vars_c1_true[focal_var]
  true_own <- matrix(NA, length(cN), nhor)
  for (k in seq_along(cN)) {
    vars_k_t  <- which(grepl(paste0("^", cN[k], "\\."),
                             colnames(Xraw_full)[keep_idx]))
    resp_id_t <- vars_k_t[focal_var]
    true_own[k, ] <- irf_true_full[resp_id_t, shock_idx_t, ]
  }
}

## ---- Figure 1: own-country IRFs ------------------------------------------
pdf("figures/own_country_irfs.pdf", width = 10, height = 6.5)
nc <- length(cN)
ncols <- ceiling(sqrt(nc))
nrows <- ceiling(nc / ncols)
old.par <- par(mfrow = c(nrows, ncols), mar = c(3.3, 3.3, 2.0, 0.6),
               mgp = c(2.0, 0.6, 0))
for (k in seq_along(cN)) {
  ylvec <- c(own_irfs[k, , ], 0)
  if (!is.null(true_own)) ylvec <- c(ylvec, true_own[k, ])
  yl <- range(ylvec, na.rm = TRUE)
  plot(seq_len(nhor), own_irfs[k, , 2], type = "l", lwd = 2,
       ylim = yl, xlab = "horizon", ylab = "response",
       main = sprintf("%s response to %s shock",
                      cN[k], colnames(Y)[shock_idx]))
  polygon(c(seq_len(nhor), rev(seq_len(nhor))),
          c(own_irfs[k, , 1], rev(own_irfs[k, , 3])),
          col = rgb(0.2, 0.3, 0.8, 0.25), border = NA)
  lines(seq_len(nhor), own_irfs[k, , 2], lwd = 2, col = "navy")
  if (!is.null(true_own))
    lines(seq_len(nhor), true_own[k, ], lwd = 1.8, col = "red3", lty = 2)
  abline(h = 0, col = "grey50", lty = 2)
  if (k == 1 && !is.null(true_own))
    legend("topright", c("posterior median", "90% band", "true (DGP)"),
           col = c("navy", rgb(0.2, 0.3, 0.8), "red3"),
           lty = c(1, NA, 2), lwd = c(2, NA, 1.8),
           pch = c(NA, 15, NA), bty = "n", cex = 0.75)
}
par(old.par)
dev.off()
message("Saved: figures/own_country_irfs.pdf")

## ---- Figure 2: cross-country shrinkage heatmap ----------------------------
# For each equation (M of them) we summarise the magnitude of the implied
# posterior shrinkage on coefficients linking it to other countries.  We use
# the absolute posterior median of the international coefficients, aggregated
# at the country pair level.
shrink_mat <- matrix(0, length(cN), length(cN))
rownames(shrink_mat) <- colnames(shrink_mat) <- cN
for (eqn in seq_len(M)) {
  cn_eqn <- substr(colnames(Y)[eqn], 1, 2)
  # K = M*p; lag-l coefficients sit in columns ((l-1)*M+1):(l*M).  We
  # average absolute medians across all lags.
  med <- apply(A_post[, , eqn], 2, median)
  abs_by_var <- numeric(M)
  for (lg in seq_len(laglen))
    abs_by_var <- abs_by_var + abs(med[((lg - 1) * M + 1):(lg * M)])
  abs_by_var <- abs_by_var / laglen
  for (cc in cN) {
    sl <- grepl(paste0("^", cc, "\\."), colnames(Y))
    shrink_mat[cn_eqn, cc] <- shrink_mat[cn_eqn, cc] + sum(abs_by_var[sl])
  }
}

pdf("figures/shrinkage_heatmap.pdf", width = 6.5, height = 5.5)
par(mar = c(4.5, 4.5, 2.5, 4.5))
zlim <- range(shrink_mat)
image(seq_along(cN), seq_along(cN), t(shrink_mat[length(cN):1, ]),
      axes = FALSE, xlab = "driver country", ylab = "response country",
      main = "Posterior |median| of cross-country coefficients (avg over lags)",
      col = hcl.colors(64, palette = "YlOrRd", rev = TRUE),
      zlim = zlim)
axis(1, at = seq_along(cN), labels = cN)
axis(2, at = seq_along(cN), labels = rev(cN), las = 1)
box()
# Simple colour scale legend.
lvls <- pretty(zlim, 5)
mtext(sprintf("range: [%.3f, %.3f]", zlim[1], zlim[2]),
      side = 4, line = 1, cex = 0.85)
dev.off()
message("Saved: figures/shrinkage_heatmap.pdf")

## ---- Persist results -------------------------------------------------------
save(fit, irf_arr, own_irfs, shrink_mat,
     file = "data/replication_output.rda")
message("Done. Posterior + IRFs in data/replication_output.rda")
