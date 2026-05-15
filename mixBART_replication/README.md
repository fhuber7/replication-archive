# mixBART / flexBART Replication Package

Replication code for the family of Bayesian VARs with **BART-based
conditional means and a (mixture-of-) factor stochastic-volatility shock
structure** estimated in this project. The headline specification is
`mixBART` — a finite mixture of BART predictors on each equation — but the
same driver also fits standard BART (`model = "BART"`) and the full-VAR
variant (`fullBART`).

## Model in one paragraph

For each equation `m = 1, ..., M` we model

    y_{m,t} = f_m( x_t )  +  Lambda_{m, .} F_t  +  exp( h_{m,t} / 2 ) eps_{m,t}

where `x_t` collects the `p` lags of `y_t`, `f_m` is a Bayesian Additive
Regression Tree (BART) ensemble or a finite mixture of such ensembles
(`mixBART`), `F_t` is a vector of `Q = ledermann(M)` latent factors with
log-variance `g_{q,t}` following an AR(1) (factor stochastic volatility), and
the idiosyncratic log-variances `h_{m,t}` are either constant (`sv = "homo"`),
AR(1) (`sv = "SV"`), or themselves BART-driven (`sv = "heteroBART"`). The
factor loadings receive a horseshoe prior with column-wise global shrinkage.

Forecasts are produced by recursively simulating the predictive density
`fhorz` steps ahead, either with the exact integration (`fc.approx = "exact"`)
or a cheap closed-form approximation that conditions on the posterior-mean
BART fit.

## Files

| File | Purpose |
|---|---|
| `estimation_functions/flexBART.R` | Headline sampler. `flexBART(Yraw, ...)` fits BART / mixBART means with the factor-SV shock structure and returns predictive draws. |
| `estimation_functions/fullBART.R` | Variant that grows a single BART ensemble over the full stacked VAR — slower but pools information across equations. |
| `estimation_functions/errorBART.R` | BART on the *residuals* (heteroskedastic-BART variant invoked when `sv = "heteroBART"`). |
| `estimation_functions/qr.R`       | QR helpers used by the factor block. |
| `aux_func.R`            | Utilities: `mlag`, `get.hs`, factor-loadings draws, KL / qwCRPS scoring. |
| `fcst_script.R`         | Minimal pre-existing forecasting demo (kept for reference; superseded by `main_replication.R`). |
| `simulate_data.R`       | `simulate_mixbart_data(T, M, p)`: a small VAR with `tanh`/`sin` nonlinearity and mild SV. |
| `main_replication.R`    | Driver: simulate -> fit -> plot forecast fan and one-step density -> save. |
| `figures/`              | Output PDFs/PNGs and `mixbart_fit.rds`. |

## Runtime warning

Each Gibbs sweep grows BART trees for each of the `M` equations and runs an
exact factor-SV update; runtime grows linearly in `(nburn + nsave)` and in
`T`, and is sensitive to the number of trees and the variable count. With
`T = 200`, `M = 5`, default tree count and `nsave = nburn = 500`,
`main_replication.R` finishes in a few minutes on a laptop. Production runs
in the paper use several thousand draws and a larger `M`.

For exploration:

* keep `M` moderate (<= 8) and reduce `nsave + nburn` to `~200 + 200`;
* set `sv = "homo"` to avoid the factor-SV cost;
* keep `fhorz` short (<= 12).

## How to run

```sh
cd mixBART_replication
Rscript main_replication.R
```

To plug in your own data: open `main_replication.R`, locate the line marked
`REPLACE WITH YOUR DATA`, and replace the `simulate_mixbart_data(...)` call
with code that assigns a numeric `(T x M)` matrix to `Ytrn`. Standardisation
is performed inside `flexBART`.

To switch model variant:

* `sl.mod <- "BART"` -> standard per-equation BART;
* `sl.mod <- "mixBART"` -> finite-mixture BART (default);
* call `fullBART(Ytrn, ...)` instead of `flexBART(Ytrn, ...)` for the full
  stacked variant.

## Dependencies

```
invgamma          stochvol          factorstochvol
MASS              dbarts            mvtnorm
ggplot2           reshape2          GIGrvg
Matrix
```

All are on CRAN. Install with

```r
install.packages(c("invgamma","stochvol","factorstochvol","MASS","dbarts",
                   "mvtnorm","ggplot2","reshape2","GIGrvg","Matrix"))
```

`dbarts` provides the BART backend; `factorstochvol` provides the factor-SV
update.

## Outputs

After `Rscript main_replication.R`:

* `figures/forecast_fan.pdf` / `.png` — pointwise 16/50/84% predictive bands
  for each variable over `fhorz` periods.
* `figures/onestep_density_y1.pdf` — one-step-ahead predictive density for
  the first variable.
* `figures/mixbart_fit.rds` — full posterior object returned by `flexBART`
  (`fcst`, `Hfcst`, plus, in some variants, `PHI`).

## What is *not* included

* The cluster-job scaffolding (`SGE_TASK_ID`, snowfall) used in the original
  paper. The replication runs locally and serially.
* Real macro data: the headline empirical exercise uses FRED-MD / FRED-QD.
  Wire in your own panel as described above.
* The full out-of-sample loop with rolling re-estimation. The driver fits
  one window and reports the corresponding predictive density.
