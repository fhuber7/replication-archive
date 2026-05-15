# TVP-SVD Replication Package

Self-contained replication of the headline **TVP-SVD** estimator of

> Hauzenberger, N., Huber, F., Koop, G. and Onorante, L. (2021).
> *Fast and Flexible Bayesian Inference in Time-Varying Parameter
> Regression Models.* Journal of Business and Economic Statistics.

The paper proposes a **singular value decomposition (SVD)** of the
expanded design matrix that turns a time-varying parameter (TVP)
regression into an equivalent constant-coefficient regression of much
lower effective rank.  Combined with a sparse-mixture g-prior on the
coefficients, this delivers a fast Bayesian sampler whose per-iteration
cost is linear in the sample size.

## Model

The single-equation TVP regression is

    y_t = x_t' beta_t + e_t,        e_t ~ N(0, sigma_t^2),
    beta_t = b + b_tilde_t,         b_tilde_t ~ N(mu_{S_t}, sigma_t^2 * Theta),

where `b` is the time-invariant component, `b_tilde_t` is the
time-varying deviation, `S_t in {1, ..., G}` is a latent regime
indicator drawn from a sparse Dirichlet mixture, and stochastic
volatility (`stochvol` package) is placed on the observation noise.

The SVD trick rewrites `y = Z * vec(b_tilde) + ...` where `Z` is a
sparse block-diagonal `T x (TK)` matrix.  Computing the thin SVD of
`Z` once at the start lets the sampler draw `b_tilde` in `O(T)` per
MCMC iteration via the Bhattacharya-Chakraborty-Mallick algorithm.

Two variants of the state evolution are implemented (set
`tvp.type` in `main_replication.R`):

* **`WN-SVD`** -- white-noise innovations with sparse-mixture
  clustering of the time-varying components (the headline
  specification of the paper).
* **`RW-SVD`** -- random-walk innovations with a ridge prior.

## How to run

From within `TVPSVD_replication/`:

```r
# default MCMC settings (nsave = 2000, nburn = 2000):
Rscript main_replication.R

# quick smoke test (nsave = 30, nburn = 20):
Rscript main_replication.R 30 20 1
```

or interactively

```r
setwd("TVPSVD_replication")
source("main_replication.R")
```

The script

1. simulates a synthetic TVP regression (T = 200, M_x = 4 regressors,
   three regimes with structural breaks at t = 70 and t = 140);
2. estimates the TVP-WN-SVD model with sparse-mixture g-prior and
   stochastic volatility;
3. summarises the posterior of `beta_t` and computes an `h = 4`-step
   ahead predictive density;
4. writes figures to `figures/`:
   * `tvp_coefficients.pdf/png` -- posterior median (blue), 68 % and
     90 % credible bands of `beta_t,k` plotted against the true path
     (red dashed);
   * `predictive_fan.pdf/png` -- predictive fan chart for
     `y_{T+1}, ..., y_{T+h}`.

## Files

| File | Contents |
| ---- | -------- |
| `main_replication.R` | Driver: data, estimation, plots. |
| `tvpsvd_estim.R` | `TVPSVD_estim()` -- self-contained SVD sampler (WN-SVD and RW-SVD variants), plus helper functions. |
| `simulate_data.R` | `simulate_data()` -- synthetic regime-switching TVP regression.  **REPLACE WITH YOUR DATA** for a real application. |
| `figures/` | Output PDFs/PNGs created on each run. |
| `README.md` | This file. |

## Plugging in real data

Replace the call to `simulate_data()` in `main_replication.R` with code
that returns a `T x 1` matrix `y` and a `T x K` design matrix `X`
(include a column of ones for an intercept if desired).  The estimator
itself is fully generic.

The hyperparameter that most strongly governs the amount of allowed
time variation is `ub.sc` (upper-bound scaling of the innovation
variance).  Recommended starting values:

* tightly standardised macro data (e.g. the NKPC inflation exercise
  in the paper): `ub.sc = 0.005`;
* synthetic / moderately noisy data with sizeable regime shifts (as
  in the supplied example): `ub.sc = 0.02`-`0.05`.

## Dependencies

Estimation:

* `coda`, `GIGrvg`, `MASS`, `Matrix`, `shrinkTVP`, `stochvol`,
  `bayesm`, `zoo`, `mvtnorm`.

Plotting:

* `ggplot2`, `reshape2`.

The original replication package used an `Rcpp` forward-filtering
backward-smoothing C++ routine for the random-walk FFBS baseline.
**That C++ dependency is not required here** -- the SVD path is pure
R and uses sparse matrix algebra from the `Matrix` package.

Install missing packages with

```r
install.packages(c("coda","GIGrvg","MASS","Matrix","shrinkTVP",
                   "stochvol","bayesm","zoo","mvtnorm",
                   "ggplot2","reshape2"))
```

## Citation

If you use this code please cite

```
@article{hauzenberger2021tvpsvd,
  author  = {Hauzenberger, Niko and Huber, Florian and
             Koop, Gary and Onorante, Luca},
  title   = {Fast and Flexible {B}ayesian Inference in
             Time-Varying Parameter Regression Models},
  journal = {Journal of Business and Economic Statistics},
  year    = {2021}
}
```
