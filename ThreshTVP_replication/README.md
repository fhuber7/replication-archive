# Threshold TVP-VAR replication package

Should I stay or should I go? A latent threshold approach to large-scale
mixture innovation models (Journal of Applied Econometrics).

This folder is a self-contained, minimal replication kit for the
latent-threshold time-varying-parameter VAR with stochastic volatility
proposed in the paper. It bundles a small synthetic-data simulator, a
driver script that estimates the model and produces the main diagnostic
plots, and the source of the `threshtvp` R package on which the
replication depends.

## What the model does

Each VAR coefficient follows a quasi-constant path until the absolute
size of its proposed innovation exceeds a coefficient-specific, data-driven
threshold; at that point the coefficient jumps. This combines the smooth,
random-walk-style time variation typical of TVP-VARs with the sparser,
regime-switching behaviour typical of Markov-switching VARs. Volatility
is modelled with a standard log-stochastic-volatility process estimated
equation-by-equation. Identification of impulse responses is recursive
(Cholesky), and time variation in the IRFs is obtained by drawing
posterior coefficient matrices at every sample date `t` and using a
Sims-style recursion (`impulsdtrf`).

The estimator is implemented as the `threshtvp` R package; the main
user-facing function is `estimate_tvp()`.

## Folder contents

| File / folder                      | Purpose                                                                 |
|------------------------------------|-------------------------------------------------------------------------|
| `simulate_data.R`                  | Builds a synthetic VAR(p) with regime shifts; marked REPLACE WITH YOUR DATA. |
| `main_replication.R`               | End-to-end driver: loads data, runs MCMC, computes IRFs, writes figures.|
| `auxilliary_functions.R`           | Helpers from the original JAE replication (`mlag`, `impulsdtrf`, `get_companion`, ...). |
| `threshtvp_source/threshtvp/`      | Patched package source (R + C++).                                       |
| `threshtvp_source/threshtvp_0.8.tar.gz` | Pre-built source tarball of the patched package.                   |
| `figures/`                         | Output PDFs and PNGs.                                                   |
| `results/`                         | Saved MCMC objects (`model_VAR.rds`, `irfs.rds`).                       |
| `data/`                            | Where the synthetic CSV is dropped and where user data should live.     |

## Where `threshtvp` comes from

`threshtvp` was written by Florian Huber to accompany the paper. The
canonical sources are:

- `threshtvp_source/threshtvp/`            -- the package source tree
- `threshtvp_source/threshtvp_0.8.tar.gz`  -- a freshly built tarball

The version here carries two minor patches relative to the original
release so that it builds and runs against a current R / RcppArmadillo /
stochvol stack (R >= 4.x):

1. `src/threshold_functions.cpp` -- removed an obsolete `"std"` string
   argument from one `arma::solve()` call (current Armadillo only
   accepts `solve_opts::...` flags here, and the default is the right
   thing anyway).
2. `R/MCMC_function.R` -- replaced the deprecated `stochvol::svsample2`
   single-step sampler with a thin wrapper around
   `stochvol::svsample_fast_cpp` (`.sv_step`), which is the modern
   equivalent of the original one-iteration Gibbs step.

These patches do not alter the model -- they just keep the original
algorithm compilable / runnable.

## Dependencies

R packages (installed automatically by the package, but you may have
to install some of them by hand on a fresh R):

- `stochvol`
- `GIGrvg`
- `MASS`
- `dlm`
- `Rcpp` (>= 0.11.6) and `RcppArmadillo` (for compilation)
- `snowfall` (used internally by `estimate_tvp` for parallel
  equation-by-equation estimation)
- `ggplot2` and `reshape2` (used only by the plotting code in
  `main_replication.R`)

Install them, if missing, with

```r
install.packages(c("stochvol", "GIGrvg", "MASS", "dlm",
                   "Rcpp", "RcppArmadillo",
                   "snowfall", "ggplot2", "reshape2"))
```

A working C++ toolchain is required (Xcode CLT on macOS, Rtools on
Windows, `r-base-dev` on Debian/Ubuntu).

## Installing `threshtvp`

From the shell, in this folder:

```sh
R CMD INSTALL threshtvp_source/threshtvp_0.8.tar.gz
```

or, equivalently, from R:

```r
install.packages("threshtvp_source/threshtvp_0.8.tar.gz",
                 repos = NULL, type = "source")
```

If you prefer to install directly from the source tree (e.g. after
further patching):

```sh
R CMD INSTALL threshtvp_source/threshtvp
```

`threshtvp` is not on CRAN. If compilation fails on your system, please
inspect the build log for the failing translation unit -- the most
common cause on modern toolchains is an Armadillo signature change in
`arma::solve()`, which has already been patched here.

## How to run

From the shell:

```sh
cd ThreshTVP_replication
Rscript main_replication.R
```

or, interactively in R, set the working directory to this folder and
`source("main_replication.R")`.

The default configuration is

- synthetic data, T = 200, M = 4, p = 2, two regime shifts,
- 2000 burn-in + 2000 saved MCMC draws,
- single-CPU (set `n_cpu <- 4` for snowfall-based parallelism over
  equations),
- impulse responses to a shock to the last variable, recursively
  identified (Cholesky).

The script writes PDF and PNG versions of four diagnostic figures into
`figures/`:

- `tvp_coefficients.{pdf,png}` -- posterior median and 68% bands of the
  time-varying own-first-lag coefficient in each equation.
- `irf_slices.{pdf,png}` -- IRFs at three selected sample dates.
- `irf_heatmap.{pdf,png}` -- heat-map of the median IRF of variable 1
  across (time, horizon).
- `stochastic_volatility.{pdf,png}` -- posterior SV path per equation.

Two posterior objects are saved to `results/`:

- `model_VAR.rds` -- the full `estimate_tvp()` return value.
- `irfs.rds` -- the time-varying IRF posterior medians and bands.

## Using your own data

In `main_replication.R`:

1. Set `use_synthetic <- FALSE`.
2. Set `my_csv <- "data/<your file>.csv"`. The CSV should have one
   column per time series and one row per time period; pre-transform
   variables (logs, differences, ...) as appropriate. The script reads
   it via `read.csv()` and coerces to numeric.

The rest of the pipeline (lag construction, MCMC, IRFs, plotting) is
unchanged.

## Reproducibility note

The MCMC sampler uses `set.seed` only inside `simulate_data.R`; the
sampler itself does not pin a seed by default. Run-to-run differences in
the posterior summaries are therefore expected at the order of Monte
Carlo error.

## Citation

If you use this code, please cite the JAE paper for the model and
algorithm. The package and the patches in this folder were authored by
Florian Huber.
