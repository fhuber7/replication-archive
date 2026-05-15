# SAVS TVP-VAR Replication Package

Replication code for the **time-varying parameter VAR with Signal Adaptive
Variable Selector (SAVS)** sparsification. Each equation of the VAR is
estimated as a single TVP regression with one of five global-local
shrinkage priors (Horseshoe, LASSO, NG, SSVS / MNIG, Dirichlet-Laplace),
followed by an *exact* SAVS thresholding step that zeroes out every
coefficient whose `|beta|^nu * ||x||^2 < 1`. The result is a TVP-VAR with a
data-driven sparsity pattern at every point in time.

## Model in one paragraph

For each equation `i = 1, ..., M`, the regression at time `t` is

    y_{i,t} = sum_{j < i} (-A^0_{ij,t}) y_{j,t}
            +  sum_{r=1}^K  beta_{ir,t} x_{r,t}
            +  exp(h_{i,t} / 2) eps_{i,t},

with `x_{r,t}` the stacked `K = M*p` lags of `y_t`, the contemporaneous
loadings `A^0_{ij,t}` recovered from the recursive Cholesky ordering, and
the log-volatility `h_{i,t}` an AR(1) (drawn via `stochvol`). The TVP
coefficients `beta_{ir,t}` follow a Gaussian random walk; their innovation
variances receive a global-local shrinkage prior (`shrinkage` argument).
After each Gibbs sweep the SAVS step solves the per-coefficient
penalised-least-squares problem with penalty
`|beta|^{-nu}` (`nu` argument), producing an exactly-sparse draw.

Reduced-form coefficients, Cholesky impulse responses, and recursive
forecasts are recomposed from the equation-wise posteriors.

## Files

| File | Purpose |
|---|---|
| `ng_SAVS.R`             | The headline sampler. `sparse.SAVS(Y, X, ...)` returns `A.thrsh` (sparsified TVP coefficients) and `hv` (log-volatility). |
| `auxilliary_functions.R`| `mlag`, `impulsdtrf`, `get.pred`, `sparsify`, `get_companion`. Also `sourceCpp`s `threshold_functions.cpp`. |
| `threshold_functions.cpp` | Rcpp helpers (per-coefficient SAVS thresholding inside the inner loop). |
| `estim.VAR.R`           | The original empirical driver (FRED-QD via `datasetMNg.RData`). Kept for reference. |
| `datasetMNg.RData`      | Pre-computed FRED-QD slice used by `estim.VAR.R` (small / medium / large information sets). |
| `simulate_data.R`       | `simulate_savs_tvp_var(T, M, p, sparsity, ...)` -- sparse TVP-VAR DGP with random-walk drift on active entries. |
| `main_replication.R`    | End-to-end synthetic-data driver: simulate -> fit -> plot with truth overlays. |
| `figures/`              | Generated PDFs / PNGs. |

## Truth overlays

`main_replication.R` builds four figures, every one with the DGP truth
overlaid in red:

* `sparsity_pip.pdf`     -- posterior inclusion probabilities at the last
  in-sample period with the true active set marked.
* `tvp_coefficients.pdf` -- posterior median + 68% band of selected
  time-varying coefficients vs. the true random-walk paths.
* `irf_slices.pdf`       -- Cholesky IRFs at three slice times vs. the IRFs
  computed from the true coefficient paths.
* `forecast_fan.pdf`     -- recursive predictive density vs. the realised
  future y simulated from the same DGP.

## How to run

```sh
cd SAVS_TVP_VARs_replication
Rscript main_replication.R
```

To plug in your own data: open `main_replication.R`, locate the
`REPLACE WITH YOUR DATA` block, and assign a numeric `(T x M)` matrix to
`Yraw`. The script standardises internally. To reproduce the original
empirical exercise, use `estim.VAR.R` instead -- it loads the FRED-QD
information sets directly from `datasetMNg.RData`.

## Choosing the prior and the SAVS exponent

```r
prior.choice <- "horse"   # "horse", "LASSO", "SSVS", "NG", "DL"
nu           <- 2         # SAVS penalty exponent; larger nu -> sparser fits
```

`nu = 0` is no thresholding; `nu = 2` is a reasonable default. The
estimator does the full Bayesian update for the shrinkage prior; SAVS is a
post-processing step that picks a *single* sparse representative per draw.

## Runtime warning

Each Gibbs sweep does, **per equation**, `T` Kalman-filter style updates
of `(K + i - 1)` coefficients plus the SAVS thresholding loop. Total cost
is roughly `O(M * (nsave + nburn) * T * K)`. With `M = 3`, `T = 150`,
`p = 2`, `K = M*p = 6`, `nsave = nburn = 500`, `main_replication.R`
finishes in a few minutes on a laptop. The empirical driver
(`estim.VAR.R`) uses larger information sets and `nsave = nburn = 5000`;
budget several hours per chain.

## Dependencies

```
Rcpp        bvarsv      GIGrvg      mvtnorm
stochvol    MASS        ggplot2     reshape2
```

Install with

```r
install.packages(c("Rcpp","bvarsv","GIGrvg","mvtnorm","stochvol","MASS",
                   "ggplot2","reshape2"))
```

`threshold_functions.cpp` is compiled on first run via
`Rcpp::sourceCpp`; a working C++ toolchain (Xcode CLT on macOS, Rtools on
Windows) is required.

## What is *not* included

* The cluster scaffolding from the original paper (the sampler is
  embarrassingly parallel across equations and the authors run each
  equation on its own SGE node).
* The grid loop over `(prior, nu)` combinations -- `main_replication.R`
  fits a single configuration; see the top of `estim.VAR.R` for the grid.
