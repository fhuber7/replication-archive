# TBVECM Replication Package

Replication code for the threshold Bayesian vector error correction
model (TBVECM) developed in

> Huber, F. and Zoerner, T. O. (2019). *Threshold cointegration in
> international exchange rates: A Bayesian approach.* International
> Journal of Forecasting 35, 458--473.

This package provides a clean, self-contained R implementation of the
estimator together with a synthetic data-generating process so the
sampler can be exercised without access to the proprietary monetary
data set used in the paper.

## Model overview

For an `M`-dimensional vector `y_t` of cointegrated levels, the
three-regime threshold VECM is

```
Delta y_t = (alpha_{s_t})' beta' z_{t-1}
          + (A_{s_t})' x_t
          + eps_t,        eps_t ~ N(0, Sigma_{s_t})
```

where

- `z_{t-1}` is the lagged level vector that enters the cointegrating
  relation,
- `x_t` stacks lagged differences and (optionally) a constant,
- `beta` is the cointegrating matrix of dimension `M x r` (rank `r`),
  linearly normalised by stacking an `r x r` identity in its top
  block,
- the regime indicator `s_t in {1, 2, 3}` is determined by whether the
  threshold variable
  `thrsh_t = (z_{t-1}' beta)_{g}` (the `g`-th cointegrating residual)
  falls below `gamma_1`, between `gamma_1` and `gamma_2`, or above
  `gamma_2`,
- the alpha-loadings, the short-run matrix `A` and the innovation
  covariance `Sigma` are regime specific.

The MCMC sampler cycles through the following blocks:

1. Conditional posteriors of `(alpha_{s}, A_{s})` per regime under a
   Normal-Gamma shrinkage prior.
2. Normal-Gamma hyperparameter updates (local scales `theta`, global
   scale `lambda^2`, concentration `a` via a random-walk Metropolis
   step).
3. Cointegrating vector `beta`: Gaussian conditional posterior on the
   free `(M-r) x r` block.
4. Regime-specific innovation covariance `Sigma_{s}` from an adaptive
   inverse-Wishart prior.
5. Threshold parameters `gamma_1`, `gamma_2` via Griddy-Gibbs with
   cubic-spline smoothing of the conditional posterior.

## Files

- `tbvecm_estim.R` -- core estimator `tbvecm_estim()` plus the helper
  functions actually invoked by the main specification (lifted from
  the original `aux_natural_conjugate.R`).
- `simulate_data.R` -- `simulate_tbvecm_data()` produces synthetic
  data with a known three-regime threshold structure.
- `main_replication.R` -- driver: simulates data, runs the sampler,
  saves draws, and produces the figures in `figures/`.
- `figures/` -- generated PDFs (regime probabilities, regime-
  conditional impulse responses, posterior of `beta`).
- `posterior_summaries.rds` -- list with the posterior draws for the
  cointegrating vector, the thresholds, and the regime probabilities.

## How to run

From the package directory:

```
Rscript main_replication.R
```

The defaults are `nsave = 3000` posterior draws after `nburn = 2000`
burn-in iterations. Two command-line arguments override these (e.g.
for a quick sanity check):

```
Rscript main_replication.R 30 20    # smoke test (a few seconds)
Rscript main_replication.R 5000 5000
```

Inside R or RStudio, source the driver instead:

```r
setwd("path/to/TBVECM_replication")
source("main_replication.R")
```

To use your own data, replace the block flagged `REPLACE WITH YOUR
DATA` in `main_replication.R` with code that assigns a `(T x M)`
matrix of levels to `Yraw`.

## Dependencies

- `MASS`, `GIGrvg`, `Matrix`, `mvtnorm`, `splines`, `bayesm`
  (required).
- `ggplot2`, `reshape2` (optional; the regime-probability plot falls
  back to base graphics if these are missing).

Install with:

```r
install.packages(c("MASS", "GIGrvg", "Matrix", "mvtnorm",
                   "bayesm", "ggplot2", "reshape2"))
```

## Notes and caveats

- The synthetic DGP in `simulate_data.R` is calibrated for `M = 3`
  variables, three regimes, rank `r = 1`, and stationary differenced
  dynamics in each regime. Extend it manually if you need a different
  dimension or rank.
- The Griddy-Gibbs step requires each regime to retain at least
  `M + 2` observations along the grid; the threshold update is
  rejected for grid points that violate this. If the threshold
  collapses the sampler nudges the regimes back to their initial
  quantile positions and continues.
- The prior on the free elements of the cointegrating vector defaults
  to `beta_prior_var = 1.0`, which is deliberately diffuse for the
  synthetic illustration. Tighten this (the paper uses `0.01`) when
  the cointegrating relation is theory-driven and known up to sign,
  e.g. for the exchange-rate application in the published paper.
- The estimator only implements the main specification of the paper
  (linearly normalised `beta`, NG shrinkage on `alpha` and `A`,
  Griddy-Gibbs thresholds, adaptive inverse-Wishart for `Sigma`).
  Identification variants and the forecasting/IRF code paths from
  the original repository are not reproduced here.
