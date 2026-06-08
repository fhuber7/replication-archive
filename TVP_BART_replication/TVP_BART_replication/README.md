# TVP-BART Replication Package

Replication code for the time-varying parameter VAR with regression-tree
driven coefficients developed in

> Hauzenberger, N., Huber, F., Koop, G. and Mitchell, J. *Bayesian Modeling
> of Time-Varying Parameter Vector Autoregressions Using Regression Trees.*
> Annals of Applied Statistics.

The model is a semiparametric TVP-VAR in which the time variation in the VAR
coefficients and in the error variance-covariance matrix is driven by
nonlinear functions of a low-dimensional state, each modelled
nonparametrically with Bayesian Additive Regression Trees (BART).

This package ships **two** things:

1. a clean, self-contained R implementation of the estimator together with a
   synthetic data-generating process, so the sampler can be exercised end to
   end (simulate -> estimate -> impulse responses) without any external data
   (`simulate_data.R`, `tvpbart_estim.R`, `main_replication.R`); and
2. the original empirical application from the paper (`estim_file.R`,
   `data/`, `aux_funcs/`), which reproduces the US macro results but requires
   the proprietary SoftBart back-end and the FRED-QD / EBP data.

The notes below describe the self-contained replication; the empirical files
are documented at the bottom.

## Model overview

For an `M`-dimensional vector `y_t`, the estimator works equation by equation
in a recursive (Cholesky) order.  Equation `i` is

```
y_{i,t} = sum_k d_{ik,t} * a_{ik,t} + u_{i,t},   u_{i,t} ~ N(0, sigma^2_{i,t})
a_{ik,t} = BART_{ik}(z_t)
```

where

- `d_{i,t}` stacks the contemporaneous variables ordered before `i` (this is
  what delivers the recursive identification), the `p` lags of all variables
  and a constant;
- every coefficient `a_{ik,t}` is its own sum-of-trees function of a
  low-dimensional state `z_t` (a normalised time trend plus the scaled first
  lags), so the coefficients drift nonparametrically over the sample;
- `sigma^2_{i,t}` follows a stochastic-volatility process, so the conditional
  mean **and** the covariance are time-varying.

The MCMC sampler cycles through, per equation:

1. a Gibbs sweep over each coefficient's BART ensemble (grow/prune
   Metropolis-Hastings tree moves with conjugate Gaussian leaves, using the
   observed regressor `d_{ik,t}` as the leaf-likelihood design weight);
2. a stochastic-volatility update of the residual variance (`stochvol`).

Impulse responses are local-linear Cholesky responses: at each date the system
is treated as the linear VAR implied by that date's coefficient matrix and
covariance, so both the coefficient drift and the volatility produce
**time-varying** impulse responses.

## Files

- `tvpbart_estim.R` -- core estimator `tvpbart_estim()` plus a hand-rolled
  weighted BART (the `.bt_*` / `.bart_*` helpers) and the impulse-response
  routine.
- `simulate_data.R` -- `simulate_tvpbart_data()` produces a VAR whose
  coefficients drift smoothly between two stationary regimes and whose
  volatility humps mid-sample; `true_irf()` returns the implied true impulse
  responses for the diagnostic overlays.
- `main_replication.R` -- driver: simulates data, runs the sampler, computes
  impulse responses, writes the figures and saves the posterior summaries.
- `figures/` -- generated PDFs (time-varying coefficients, time-varying and
  sample-averaged impulse responses, volatilities).
- `posterior_summaries.rds` -- list with the posterior medians/bands of the
  coefficients, the impulse responses and the volatilities.

## How to run

From the package directory:

```
Rscript main_replication.R
```

The defaults are `nsave = 1000` posterior draws after `nburn = 1000` burn-in
iterations with `num.trees = 15` trees per coefficient (a few minutes in pure
R).  Up to three command-line arguments override `nsave`, `nburn` and
`num.trees`:

```
Rscript main_replication.R 150 150 8      # smoke test (about a minute)
Rscript main_replication.R 5000 5000 50   # closer to the paper's resolution
```

Inside R or RStudio, source the driver instead:

```r
setwd("path/to/TVP_BART_replication")
source("main_replication.R")
```

To use your own data, replace the block flagged `REPLACE WITH YOUR DATA` in
`main_replication.R` with code that assigns a `(T x M)` matrix `Yraw` (with
column names) of the endogenous variables, and drop the `truth` overlays.

## Dependencies

- `MASS`, `stochvol` (estimator); `ggplot2`, `reshape2` (figures).

Install with:

```r
install.packages(c("MASS", "stochvol", "ggplot2", "reshape2"))
```

## Notes and caveats

- **Per-coefficient vs. latent-factor BART.** The published model compresses
  the coefficients of each equation into a small number `Q` of common latent
  BART factors with loadings.  This self-contained version models each
  coefficient with its own BART -- the `Lambda = I` special case.  It is less
  parsimonious but cleanly identified (every tree ensemble's design weight is
  an observed regressor rather than a sampled loading), which makes it robust
  in a short, synthetic demonstration.  As a result individual *reduced-form
  coefficient* paths are noisier than the impulse responses, which are the
  well-identified functionals of the coefficients and the covariance and are
  recovered accurately.
- **Identification.** Only the recursive (Cholesky) structural form is
  implemented; the maximum-variance and external-instrument identifications
  of the paper are not reproduced.
- **Impulse responses.** These are local-linear (the system is frozen at each
  date's coefficients); the paper computes fully nonlinear generalised
  impulse responses to account for the BART nonlinearity.
- **MCMC / model size.** Defaults are kept small so the routine runs in
  minutes.  Increase `nsave`, `nburn`, `num.trees` and the sample size for
  smoother estimates.
- The synthetic DGP in `simulate_data.R` is calibrated for `M = 3`, `p = 1`.
  Extend it manually for other dimensions.

## Empirical application (original paper code)

The files below reproduce the empirical application and are **not** part of the
self-contained demo; they require the SoftBart package and the data sets used
in the paper.

- `estim_file.R` -- main script for data set-up, estimation and
  post-processing of the US macro application.
- `aux_funcs/tvp_softbart_function.R` -- the full MCMC sampler.
- `aux_funcs/aux_files.R` -- helper functions.
- `data/makroUS_Q_16.rda` -- quarterly US macro data (FRED-QD).
- `data/EBPUS_Q.rda` -- excess bond premium (Gilchrist & Zakrajsek, 2012).
- `data/uncUS_Q.rda` -- auxiliary series incl. the NBER recession indicator.
