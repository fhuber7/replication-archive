# SubspaceVAR replication

Replication package for

> Huber, F. and Koop, G. (forthcoming) "Subspace Shrinkage in Conjugate
> Bayesian Vector Autoregressions", *Journal of Applied Econometrics*.

The package implements the headline **subVAR-Minn** estimator: a conjugate
matrix-Normal / Inverse-Wishart VAR whose prior covariance is a convex
combination (weight `omega`) of a Minnesota-style diffuse base prior and a
principal-component subspace projection of rank `q`. The hyperparameters
`(omega, q, theta)` are integrated out by computing the marginal likelihood
on a fixed grid and resampling grid points with the resulting posterior
weights (inverse-transform sampling).

## How to run

```
cd SubspaceVAR_replication
Rscript main_replication.R
```

By default the script

1. simulates a small T = 200, M = 8 factor-driven panel (`q_true = 2`),
2. estimates the subspace-shrinkage VAR on a 5 x 3 x 3 grid of
   `(omega, q, theta)`,
3. draws `nsave = 30` posterior samples after `nburn = 30` burn-in,
4. computes recursive Cholesky impulse responses (16 horizons) and h = 8
   predictive draws for variable 1,
5. writes three PDFs to `./figures/` and saves the posterior list to
   `subspaceVAR_fit.RData`.

Bumping `nsave`, `nburn`, `grid.omega`, `grid.q`, `grid.theta` and `irfhor`
in `main_replication.R` reproduces the production-scale runs of the paper.

## Files

| file | purpose |
|------|---------|
| `subspacevar_estim.R` | `get.dum`, `get_companion`, `impulsdtrf`, and the headline estimator `subVAR.0`. |
| `simulate_data.R`     | Factor-DGP `sim.function` and `make_synthetic_data` stub. Replace with your data for real runs. |
| `main_replication.R`  | End-to-end driver: data, grids, estimation, IRFs, fan chart, plots. |
| `figures/`            | Output PDFs (created on first run). |

## Bringing in real data

In `main_replication.R`, locate the line marked `REPLACE WITH YOUR DATA` and
substitute a `T x M` numeric matrix `Y` of stationary series with column
names. Everything downstream is data-agnostic.

## Output figures

- `figures/omega_q_heatmap.pdf` - posterior weights on the `(omega, q)`
  grid, marginalising over `theta`.
- `figures/irfs.pdf` - Cholesky IRFs with 68% credible bands, all
  shock/response pairs.
- `figures/fanchart.pdf` - fan chart of the recursive forecast for the
  first variable.

## Dependencies

R (>= 4.0) and the CRAN packages: `Matrix`, `MASS`, `mvtnorm`, `ggplot2`,
`reshape2`, `zoo`, `scales`.

Install with

```r
install.packages(c("Matrix","MASS","mvtnorm","ggplot2",
                   "reshape2","zoo","scales"))
```

## Notes on the prior

- `omega -> 0`: prior collapses to a standard diffuse / Minnesota VAR.
- `omega -> 1`: prior concentrates all mass on the q-dimensional
  principal-component subspace of the lagged regressors, i.e. on a
  factor-restricted VAR.
- `q`: rank of the subspace projection (number of leading PCs of `X`).
- `theta`: overall Minnesota tightness on the dummy-observation prior.

See Sections 2 and 3 of the paper for the closed-form posterior and the
marginal-likelihood derivation, and Appendix B for the factor DGP used in
the Monte-Carlo experiments mirrored in `simulate_data.R`.
