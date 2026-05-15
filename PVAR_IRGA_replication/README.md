# PVAR-IRGA replication package

This directory contains a clean, self-contained replication of the headline
estimator from

> Feldkircher, M., Huber, F., Koop, G. and Pfarrhofer, M.
> *Approximate Bayesian inference and forecasting in huge-dimensional
>  multi-country VARs.*

The estimator implements the **Integrated Rotated Gaussian Approximation
(IRGA)** strategy of Loaiza-Maya, Smith and Maneesoonthorn for very large
multi-country VARs. The full posterior is approximated by splitting it into:

1. an **own-country** block that is sampled exactly by MCMC under a
   Horseshoe (or Minnesota) prior; and
2. a **cross-country** block whose parameters are integrated out analytically
   using Vector Approximate Message Passing (VAMP) with a global-local
   shrinkage prior.

The procedure is **approximate**: the cross-country block is summarised by a
posterior mean and a Gaussian uncertainty rather than full draws. In exchange
the estimator scales to dimensions for which standard Bayesian VAR samplers
are infeasible.

## Model

For each country i and each equation j the IRGA regression is

```
Y_{i,j,t} = X_{i,t}' beta + Z_{i,t}' gamma + eps_{i,j,t},  eps ~ N(0, sigma^2)
```

where `X_{i,t}` collects the lags of country i's own variables (the
**domestic** block) and `Z_{i,t}` collects the lags of every other country's
variables together with the contemporaneous own-country variables that have
already been sampled (Cholesky ordering of the impact matrix).

After a QR rotation of `X_i`, the system splits into a `p`-dimensional
domestic subsystem and an `(n-p)`-dimensional nuisance subsystem that depends
only on `gamma`. VAMP returns the posterior mean and (scalar) variance of
`gamma`. Conditional on that, `beta` is sampled by an MCMC over the rotated
domestic subsystem under a Horseshoe prior.

Stacking countries gives the reduced-form panel VAR
`y_t = sum_l A_l y_{t-l} + u_t`, `u_t ~ N(0, Sigma)`, with posterior draws of
`(A_1, ..., A_p, Sigma)`. Structural identification uses a recursive (lower
triangular) Cholesky scheme on `Sigma`.

## Files

| file                   | purpose                                                |
| ---------------------- | ------------------------------------------------------ |
| `pvar_irga_estim.R`    | self-contained IRGA estimator (`PVAR.IRGA`, `get.reg.IRGA`, `VAMP`, `BVS_IRGA_12`, `vblr_fit`, Horseshoe helpers, `irf_cholesky`) |
| `simulate_data.R`      | synthetic multi-country panel generator (defaults: N=5, G=3, T=150) |
| `main_replication.R`   | end-to-end driver: load data, estimate, IRFs, figures  |
| `figures/`             | output PDF figures                                     |
| `data/`                | synthetic panel + posterior output (`.rda`)            |

## How to run

```bash
cd PVAR_IRGA_replication
Rscript simulate_data.R       # writes data/synthetic_panel.rda
Rscript main_replication.R    # estimates + plots
```

The default settings are `nsave = 1000`, `nburn = 500`, `laglen = 2`, `N = 5`
countries. Reduce `nsave/nburn` and `N_country` near the top of
`main_replication.R` for quick smoke tests.

To use **your own data**, save a matrix `Xraw` with column names
`"<CC>.<var>"` (two-character country code, dot, variable name) to
`data/synthetic_panel.rda`. The script auto-discovers countries from the
column prefixes. **REPLACE WITH YOUR DATA** at this point in the pipeline.

## Outputs

`main_replication.R` produces:

- `figures/own_country_irfs.pdf` -- responses of the focal variable across all
  countries to a Cholesky shock to the same variable in country 1, with 90%
  posterior bands.
- `figures/shrinkage_heatmap.pdf` -- posterior absolute-median magnitude of
  the cross-country VAR coefficients, aggregated to country-pair level. Darker
  cells = stronger estimated dynamic linkages.
- `data/replication_output.rda` -- the full fitted object, the IRF array, and
  the aggregated tables.

## Dependencies

R packages used at runtime:

- `Matrix`
- `MASS`
- `mvtnorm`
- `GIGrvg`

The original code base depended on `matrixcalc`, `fastSOM` and `lars`; those
have been replaced with base-R equivalents (`sum(diag(...))` for the matrix
trace, no spillover decomposition, no debiased lasso).

## Caveat

This is an **approximate** Bayesian procedure. The cross-country posterior is
summarised by its VAMP mean and a Gaussian uncertainty rather than full draws,
so reported uncertainty for the impulse responses comes only from the
own-country block. For full uncertainty quantification (e.g. for forecast
evaluation in production work) extend the loop in `PVAR.IRGA` to draw
international coefficients from the VAMP-implied Gaussian at each MCMC sweep.
