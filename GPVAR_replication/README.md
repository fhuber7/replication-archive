# GP-VAR Replication Package

Replication code for the model in **Gaussian Process Vector Autoregressions
and Macroeconomic Uncertainty**: a nonparametric VAR in which each equation's
conditional mean is the sum of two Gaussian-process (GP) functionals — one
indexed by the *own* lags, one by the *other* variables' lags — with
equation-wise stochastic volatility and a horseshoe-shrunk lower-triangular
structural covariance.

## Model in one paragraph

For each equation `m = 1, ..., M` we model

    y_{m,t} = f_1(x^{own}_{m,t})  +  f_2(x^{other}_{m,t})
              +  sum_{j<m}  upsilon_{m,j} * e_{j,t}
              +  exp(h_{m,t} / 2) * eps_{m,t}

with `f_1, f_2` a priori independent Gaussian processes equipped with
squared-exponential ("Gaussian") kernels `K_own` and `K_other`. Their
prior bandwidths `h.own` and `h.other` differ, so the model can apply more
or less smoothing to the own-lag block vs. the cross-variable block. Two
global scalings `psi_1, psi_2` (one per kernel) are sampled with an
adaptive random-walk Metropolis step (the "adptMinn" prior in the paper),
the contemporaneous structural covariance is sampled with a horseshoe
prior, and the equation-specific log-variances `h_{m,t}` follow an AR(1)
stochastic-volatility process drawn via `stochvol`.

Identification of `f_1`, `f_2` and the structural row of upsilon is
recursive: equation 1 has only `f_1 + f_2`; equation `m > 1` additionally
regresses the residual on the lower-triangular columns
`e_{1,t}, ..., e_{m-1,t}`.

Generalized impulse responses to a positive and a negative one-standard-
deviation shock in variable 1 are computed by simulating two parallel paths
(with and without the shock) through the *posterior-mean* GP predictive,
and differencing them.

## Files

| File | Purpose |
|---|---|
| `gpvar_estim.R`        | Self-contained estimator. Defines `gpvar(Y, p, nsave, nburn, ...)`, plus `GaussKernel`, `mlag`, `get.hs`. |
| `simulate_data.R`      | `simulate_gpvar_data(T, M, p)`: a small VAR with `tanh` / `sin` nonlinearity and mild SV. |
| `main_replication.R`   | Driver: simulate -> standardise -> fit -> plot GIRFs and SV path. |
| `figures/`             | Output PDFs/PNGs and `gpvar_fit.rds`. |

## Runtime warning: this is a slow model

Each Gibbs sweep performs, **per equation**, two Cholesky / `solve` calls
on `T x T` kernel matrices and one stochvol update. That is

    O( M  *  T^3 )

per iteration. Concretely, with `M = 3`, `p = 2`, `T = 118`,
`nsave = nburn = 500` the script `main_replication.R` runs in well under a
minute on a laptop. Doubling `T` makes it roughly **eight times** slower;
doubling the variables doubles it again. The published paper uses tighter
specifications and is typically run on a cluster; the defaults here are
deliberately small so that the replication terminates in interactive time.

For production runs:
* keep `T` modest (the paper's quarterly sample is 200-ish observations);
* increase `nsave + nburn` to several thousand;
* do not switch on more variables than you actually need.

## How to run

```sh
cd GPVAR_replication
Rscript main_replication.R
```

To plug in your own data: open `main_replication.R` and replace the
`simulate_gpvar_data(...)` call with code that loads your `Y` as a numeric
matrix (rows = time, columns = variables) **before** the `apply(... scale)`
line. The look for the comment `REPLACE WITH YOUR DATA`.

## Dependencies

```
MASS        mvtnorm      Matrix       stochvol
GIGrvg      ggplot2      reshape2     truncnorm
invgamma
```

All are on CRAN. Install with

```r
install.packages(c("MASS","mvtnorm","Matrix","stochvol","GIGrvg",
                   "ggplot2","reshape2","truncnorm","invgamma"))
```

## What is *not* included

* The C++ Gibbs helper `BAKRGibbs.cpp` and the subspace approximation in
  `gpsubspace.R` (only used for accelerated runs). The estimator here is
  the headline kernel version from `cluster.R`.
* `RATE.R` (relative-amplitude / variable-importance diagnostics).
* The cluster-job scaffolding (`snow`, `snowfall`, `SGE_TASK_ID`) used in
  the original paper. The replication runs locally.

## Outputs

After `Rscript main_replication.R`:

* `figures/girfs.pdf` / `.png` — posterior 68% bands for the GIRFs to a
  positive and a negative shock to variable 1.
* `figures/sv_path.pdf` / `.png` — posterior median `exp(h_t / 2)` for
  each equation.
* `figures/gpvar_fit.rds` — the full posterior object returned by `gpvar`.
