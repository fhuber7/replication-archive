# ============================================================================
# main_replication.R
#
# End-to-end driver for the TVP-BART replication package.  It
#   1. simulates a VAR whose coefficients drift smoothly over the sample and
#      whose volatility humps in the middle (see simulate_data.R);
#   2. estimates the time-varying parameter VAR with BART-driven coefficients
#      and stochastic volatility (see tvpbart_estim.R);
#   3. computes time-varying Cholesky impulse responses;
#   4. writes diagnostic figures to figures/ and saves the posterior summaries.
#
# Run with the defaults (nsave = 1000, nburn = 1000):
#   Rscript main_replication.R
# Quick smoke test (a minute or so):
#   Rscript main_replication.R 150 150 8
# Inside R / RStudio:
#   source("main_replication.R")
# ============================================================================

rm(list = ls())
set.seed(20230401)

# ---- working directory (Rscript, source() and RStudio all supported) -------
get_script_dir <- function() {
  ca <- commandArgs(FALSE)                                   # Rscript --file=
  m  <- grep("^--file=", ca)
  if (length(m)) return(dirname(normalizePath(sub("^--file=", "", ca[m[1]]))))
  for (i in seq_len(sys.nframe())) {                          # source()
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(dirname(normalizePath(of)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&       # RStudio
      !is.null(rstudioapi::getSourceEditorContext()$path))
    return(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}
setwd(get_script_dir())
dir.create("figures", showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
})

source("tvpbart_estim.R")
source("simulate_data.R")

# ----------------------------------------------------------------------------
# 1. Data
# ----------------------------------------------------------------------------
# REPLACE WITH YOUR DATA -----------------------------------------------------
# To use your own series, replace the next two lines with code that assigns a
# (T x M) matrix `Yraw` (with column names) of the endogenous variables.  The
# `truth` overlay is only used for the diagnostic plots and can be dropped.
# ----------------------------------------------------------------------------
sim   <- simulate_tvpbart_data(T = 200, M = 3, p = 1, seed = 42)
Yraw  <- sim$Yraw
truth <- sim$truth
cat("Data:", nrow(Yraw), "x", ncol(Yraw), "\n")

# ----------------------------------------------------------------------------
# 2. MCMC / model settings (override on the command line for a smoke test)
# ----------------------------------------------------------------------------
args      <- commandArgs(trailingOnly = TRUE)
nsave     <- if (length(args) >= 1) as.integer(args[1]) else 1000
nburn     <- if (length(args) >= 2) as.integer(args[2]) else 1000
num.trees <- if (length(args) >= 3) as.integer(args[3]) else 15

p        <- 1
nhor     <- 17
shk.var  <- 1          # shock the first variable (Cholesky ordering)
shk.sc   <- 1          # one-standard-deviation shock

cat(sprintf("\nEstimating TVP-BART VAR (nsave=%d, nburn=%d, num.trees=%d) ...\n",
            nsave, nburn, num.trees))
fit <- tvpbart_estim(Yraw, p = p, nsave = nsave, nburn = nburn,
                     num.trees = num.trees, sv = TRUE,
                     shk.var = shk.var, shk.sc = shk.sc, nhor = nhor,
                     n.irf.dates = 30, verbose = TRUE)
cat(sprintf("MCMC complete. Elapsed: %.2f minutes.\n", fit$time.min))

M <- fit$dims$M; T <- fit$dims$T; var.names <- fit$var.names

# ----------------------------------------------------------------------------
# 3. Posterior summaries
# ----------------------------------------------------------------------------
# reduced-form coefficients: nsave x T x M x (Mp+cons)
A.med <- apply(fit$A.store, c(2, 3, 4), median,   na.rm = TRUE)
A.lo  <- apply(fit$A.store, c(2, 3, 4), quantile, 0.16, na.rm = TRUE)
A.hi  <- apply(fit$A.store, c(2, 3, 4), quantile, 0.84, na.rm = TRUE)
A.tru <- truth$A[(p + 1):nrow(truth$A), , ]                  # aligned to sample

# time-varying IRFs and the truth at the stored dates
irf.med <- apply(fit$irf.store, c(2, 3, 4), median,   na.rm = TRUE)   # nD x M x nhor
irf.lo  <- apply(fit$irf.store, c(2, 3, 4), quantile, 0.16, na.rm = TRUE)
irf.hi  <- apply(fit$irf.store, c(2, 3, 4), quantile, 0.84, na.rm = TRUE)
irf.true.full <- true_irf(truth, shk.var, shk.sc, nhor)[(p + 1):nrow(truth$A), , ]
avg.med <- apply(fit$irf.avg, c(2, 3), median,   na.rm = TRUE)         # M x nhor
avg.lo  <- apply(fit$irf.avg, c(2, 3), quantile, 0.16, na.rm = TRUE)
avg.hi  <- apply(fit$irf.avg, c(2, 3), quantile, 0.84, na.rm = TRUE)
avg.tru <- apply(irf.true.full, c(2, 3), mean)

# volatilities
sd.med <- sqrt(apply(fit$s2.store, c(2, 3), median, na.rm = TRUE))     # M x T
sd.tru <- t(truth$sigma[(p + 1):nrow(truth$sigma), ])                  # M x T

# ----------------------------------------------------------------------------
# 4. Figure: selected time-varying coefficients vs. truth
# ----------------------------------------------------------------------------
key <- list(c(1, 1), c(1, 2), c(2, 2))            # the coefficients that drift
lab <- sapply(key, function(ij) sprintf("A[%d,%d]  (resp y%d, lag y%d)",
                                        ij[1], ij[2], ij[1], ij[2]))
df.coef <- do.call(rbind, lapply(seq_along(key), function(m) {
  ij <- key[[m]]
  data.frame(t = 1:T, panel = lab[m],
             med = A.med[, ij[1], ij[2]], lo = A.lo[, ij[1], ij[2]],
             hi  = A.hi[, ij[1], ij[2]],  truth = A.tru[, ij[1], ij[2]])
}))
p.coef <- ggplot(df.coef, aes(t)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "navyblue", alpha = 0.25) +
  geom_line(aes(y = med),   colour = "navyblue", linewidth = 0.8) +
  geom_line(aes(y = truth), colour = "red3", linetype = "dashed", linewidth = 0.7) +
  facet_wrap(~ panel, scales = "free_y") +
  labs(title = "Time-varying reduced-form coefficients (TVP-BART)",
       subtitle = "Blue: posterior median + 68% band.  Red dashed: true path.",
       x = "time", y = NULL) +
  theme_bw(base_size = 11)
ggsave("figures/tvp_coefficients.pdf", p.coef, width = 9, height = 3.4)

# ----------------------------------------------------------------------------
# 5. Figure: time-varying impulse responses at early / mid / late dates
# ----------------------------------------------------------------------------
d.idx  <- fit$irf.dates
sel    <- c(early = 3, mid = ceiling(length(d.idx) / 2), late = length(d.idx) - 2)
df.irf <- do.call(rbind, lapply(names(sel), function(nm) {
  dd <- sel[nm]; tt <- d.idx[dd]
  do.call(rbind, lapply(seq_len(M), function(v) {
    data.frame(h = 0:(nhor - 1), date = paste0(nm, " (t=", tt, ")"),
               response = var.names[v],
               med = irf.med[dd, v, ], lo = irf.lo[dd, v, ], hi = irf.hi[dd, v, ],
               truth = irf.true.full[tt, v, ])
  }))
}))
df.irf$date <- factor(df.irf$date, levels = unique(df.irf$date))
p.irf <- ggplot(df.irf, aes(h)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "navyblue", alpha = 0.22) +
  geom_line(aes(y = med),   colour = "navyblue", linewidth = 0.8) +
  geom_line(aes(y = truth), colour = "red3", linetype = "dashed", linewidth = 0.7) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  facet_grid(response ~ date, scales = "free_y") +
  labs(title = sprintf("Time-varying impulse responses to a %s shock", var.names[shk.var]),
       subtitle = "Blue: posterior median + 68% band.  Red dashed: true response.",
       x = "horizon", y = NULL) +
  theme_bw(base_size = 11)
ggsave("figures/irf_timevarying.pdf", p.irf, width = 9, height = 6)

# ----------------------------------------------------------------------------
# 6. Figure: sample-averaged IRF vs. truth
# ----------------------------------------------------------------------------
df.avg <- do.call(rbind, lapply(seq_len(M), function(v)
  data.frame(h = 0:(nhor - 1), response = var.names[v],
             med = avg.med[v, ], lo = avg.lo[v, ], hi = avg.hi[v, ],
             truth = avg.tru[v, ])))
p.avg <- ggplot(df.avg, aes(h)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "darkgreen", alpha = 0.20) +
  geom_line(aes(y = med),   colour = "darkgreen", linewidth = 0.8) +
  geom_line(aes(y = truth), colour = "red3", linetype = "dashed", linewidth = 0.7) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  facet_wrap(~ response, scales = "free_y") +
  labs(title = sprintf("Sample-averaged impulse response to a %s shock", var.names[shk.var]),
       subtitle = "Green: posterior median + 68% band.  Red dashed: true average response.",
       x = "horizon", y = NULL) +
  theme_bw(base_size = 11)
ggsave("figures/irf_average.pdf", p.avg, width = 9, height = 3.2)

# ----------------------------------------------------------------------------
# 7. Figure: stochastic volatilities vs. truth
# ----------------------------------------------------------------------------
df.vol <- do.call(rbind, lapply(seq_len(M), function(v)
  data.frame(t = 1:T, variable = var.names[v],
             med = sd.med[v, ], truth = sd.tru[v, ])))
p.vol <- ggplot(df.vol, aes(t)) +
  geom_line(aes(y = med),   colour = "navyblue", linewidth = 0.8) +
  geom_line(aes(y = truth), colour = "red3", linetype = "dashed", linewidth = 0.7) +
  facet_wrap(~ variable, scales = "free_y") +
  labs(title = "Time-varying volatilities (posterior median vs. truth)",
       x = "time", y = "residual sd") +
  theme_bw(base_size = 11)
ggsave("figures/volatilities.pdf", p.vol, width = 9, height = 3.0)

# ----------------------------------------------------------------------------
# 8. Console diagnostics + save posterior summaries
# ----------------------------------------------------------------------------
cat("\n-- Recovery diagnostics (correlation of posterior median with truth) --\n")
for (m in seq_along(key)) {
  ij <- key[[m]]
  cat(sprintf("  %-28s cor = % .2f\n", lab[m],
              cor(A.med[, ij[1], ij[2]], A.tru[, ij[1], ij[2]])))
}
cat(sprintf("  sample-averaged IRF        cor = % .2f\n",
            cor(as.numeric(avg.med), as.numeric(avg.tru))))
for (v in seq_len(M))
  cat(sprintf("  volatility of %-4s         cor = % .2f\n",
              var.names[v], cor(sd.med[v, ], sd.tru[v, ])))

saveRDS(list(A.med = A.med, A.lo = A.lo, A.hi = A.hi, A.true = A.tru,
             irf.med = irf.med, irf.lo = irf.lo, irf.hi = irf.hi,
             irf.dates = fit$irf.dates, irf.avg = avg.med, sd.med = sd.med,
             var.names = var.names, dims = fit$dims),
        file = "posterior_summaries.rds")

cat("\nFigures written to figures/:\n",
    "  tvp_coefficients.pdf  irf_timevarying.pdf  irf_average.pdf  volatilities.pdf\n",
    "Posterior summaries saved to posterior_summaries.rds\n", sep = "")
cat("Done.\n")
