# Subspace-Inflation-BNP Replication Package

Replication code for the single-equation **inflation-forecasting model with
a Bayesian nonparametric (Dirichlet-process) mixture on the shocks** and a
choice of three conditional-mean specifications:

* **GP-subspace** — Gaussian process on a low-rank SVD subspace of the
  predictors (`gp_bnp` in `gpsubspace_function.R`),
* **BART** — Bayesian Additive Regression Trees with the same shock mixture
  (`BART_bnp` / `BART_mixSV` in `BART_function.R`),
* **UCSV** — unobserved-components stochastic-volatility baseline
  (`ucsv` in `ucsv.R`).

The empirical exercise targets US headline (CPIAUCSL) or core (CPILFESL)
inflation; see `infdata_script.R`.

## Model in one paragraph

For a single target series `y_t`,

    y_t = f( x_t )  +  e_t,
    e_t ~ sum_{g=1}^{G} eta_g  N( mu_g, sd_g^2 ),

where `f` is either a Gaussian-process functional (`gp_bnp`) or a BART
ensemble (`BART_bnp`). The mixture is sampled with a slice sampler that
truncates the stick-breaking representation; the concentration parameter
follows an Escobar-and-West Gamma prior and is updated by random-walk
Metropolis. The GP variant additionally projects `x_t` onto a low-rank
SVD subspace (`PCA = TRUE`); the projection rank is set by the optimal
hard-threshold rule of Gavish-Donoho. Optionally the idiosyncratic
component carries factor-style stochastic volatility (`mix.sv = TRUE`,
which loads `gp_svdlp.R`).

The reduced "UCSV" model is the standard Stock-Watson unobserved-components
stochastic-volatility benchmark; it is provided for comparison.

## Files

```
SubspaceInflation_replication/
  README.md
  simulate_data.R                 sparse nonlinear regression generator
  main_replication.R              driver: simulate -> fit GP & BART -> figures -> save
  BAKRGibbs.cpp                   Rcpp helpers (GaussKernel_star, ...)
  gpsubspace_function.R           the GP / GP-subspace sampler `gp_bnp`
  gp_svdlp.R                      mixture-SV add-on used when mix.sv = TRUE
  BART_function.R                 `BART_bnp` and `BART_mixSV`
  ucsv.R                          stochastic-volatility UCSV baseline
  infdata_script.R                pre-existing empirical driver (FRED-QD;
                                  uses Data/makroUS_Q_h_*_<series>.rda which
                                  is not shipped with this package)
  collect.R, collect_*.R          post-processing scripts for the cluster
                                  output of the empirical exercise
  figures/                        generated plots
  output/                         generated RDS + CSV
```

## Dependencies

```
Rcpp        Matrix      stochvol    bvarsv      flexmix
dbarts      GIGrvg      ggplot2     reshape2
```

Install with

```r
install.packages(c("Rcpp","Matrix","stochvol","bvarsv","flexmix",
                   "dbarts","GIGrvg","ggplot2","reshape2"))
```

`BAKRGibbs.cpp` is compiled on first run via `Rcpp::sourceCpp`; a working
C++ toolchain (e.g. Xcode CLT on macOS, Rtools on Windows) is required.

## How to run

```sh
cd SubspaceInflation_replication
Rscript main_replication.R
```

The driver fits both the GP-subspace model and the BART model on a small
simulated `(y, X)` set, holds out the last observation, and writes the
predictive density, the in-sample fit, and the log predictive likelihoods
into `figures/` and `output/`.

To plug in your own data: open `main_replication.R`, find the line marked
`REPLACE WITH YOUR DATA`, and supply your own `y` (T x 1) and `X` (T x p).
The samplers expect both standardised; the snippet just above the marker
does the standardisation.

## Runtime warning

`gp_bnp` performs an `O(n^3)` Cholesky on the kernel matrix per Gibbs
sweep, plus the slice-sampler update over `G.max` mixture components.
Concretely: with `n = 200`, `p = 15`, `G.max = 10`, `nburn = 500`,
`ntot = 1000`, `main_replication.R` runs in a couple of minutes. The
paper's exercise uses `nburn = 2000`, `ntot = 4000`, a full hold-out
length of 185 quarters, and several `(h, model, G, data)` combinations
distributed over an SGE cluster — budget many hours per cell.

For exploration, keep `n <= 250` and `G.max <= 10`.

## Reproducing the empirical exercise

To run the empirical exercise (single hold-out cell):

```r
# inside R, with the FRED-QD slices in ./Data/
Sys.setenv(SGE_TASK_ID = "1")   # pick a row of the (h, time, model, G, data) grid
source("infdata_script.R")
```

`infdata_script.R` expects `Data/makroUS_Q_h_<h>_<series>.rda` files
containing the matrix `Xmat`. These are not shipped with the package.
The post-processing scripts (`collect.R`, `collect_densities.R`,
`collect_features.R`, `collect_ucsv.R`) aggregate the per-cell output
files into the figures and tables reported in the paper.

## Outputs

After `Rscript main_replication.R`:

* `figures/predictive_density.{pdf,png}` — overlaid kernel-density
  estimates of the one-step-ahead predictive distribution under GP-subspace
  and BART, with the realised `y` marked.
* `figures/insample_fit_gp.{pdf,png}` — posterior 16/50/84% bands of the
  GP-subspace in-sample fit, with observed `y` overlaid.
* `output/log_predictive_likelihood.csv` — summary statistics of the log
  predictive likelihood at the held-out point.
* `output/subspace_fit.rds` — the full posterior objects returned by
  `gp_bnp` and `BART_bnp`.

## What is *not* included

* The macro database snapshots used in the paper (`Data/makroUS_Q_h_*.rda`).
  They are FRED-QD slices keyed by forecast horizon and target series; the
  preparation pipeline is documented in `infdata_script.R`.
* The cluster-job scaffolding (`SGE_TASK_ID`) that distributes the rolling
  hold-out exercise.
* The frozen output cache (`Results_CPIAUCSL/` etc.) from the published
  runs; the aggregation scripts in `collect_*.R` reproduce it from the
  per-cell output files.
