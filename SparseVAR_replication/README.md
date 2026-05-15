# Sparser TVP-VARs replication package

Self-contained replication of the conjugate-Bayesian-VAR-with-sparsification
estimator from

  Hauzenberger, N., Huber, F., and Onorante, L. (2020),
  "Combining Shrinkage and Sparsity in Conjugate Vector Autoregressive Models",
  *Journal of Applied Econometrics*.

## Model

The headline specification is a Bayesian VAR(p)

    y_t = c + A_1 y_{t-1} + ... + A_p y_{t-p} + u_t,    u_t ~ N(0, Sigma * sigma2_t)

estimated with

* a Minnesota dummy-observation prior (Banbura/Sims/Zha style) that puts
  shrinkage toward white-noise (or random-walk) own dynamics and stronger
  cross-variable shrinkage;
* an optional common stochastic volatility factor sigma2_t modelled via
  `stochvol::svsample2` and extracted from the first principal component of
  the rotated residuals;
* a *post-processing* sparsification step that, at each Gibbs iteration,
  thresholds the freshly drawn coefficient matrix `Alpha` with a SAVS
  (signal adaptive variable selector) penalty whose scale and exponent
  decay across lags. This implements the "decoupled shrinkage and selection"
  (DSS) idea: shrinkage is done by the prior, then sparsity is imposed
  ex post.
* the same SAVS step is applied to the off-diagonal entries of the Cholesky
  factor of `Sigma`, yielding an equation-by-equation sparsified covariance
  matrix.

The estimator returns posterior summaries for both the un-sparsified and the
sparsified coefficient matrices, together with Cholesky impulse responses
and recursive density forecasts.

## Files

| File                  | Purpose                                                |
| --------------------- | ------------------------------------------------------ |
| `sparsevar_estim.R`   | Self-contained estimator function `sparsevar_estim()`. |
| `simulate_data.R`     | Synthetic sparse VAR(p) data-generating process.       |
| `main_replication.R`  | Driver: simulate -> estimate -> plot.                  |
| `figures/`            | Output figures (created at run time).                  |

## How to run

A full run with the defaults (`nsave = 3000`, `nburn = 2000`) takes a few
minutes on a laptop:

```
Rscript main_replication.R
```

A smoke test (very short MCMC):

```
Rscript main_replication.R 30 20
```

## Using your own data

In `main_replication.R`, replace the call to `simulate_sparse_var()` with
your own `T x M` numeric matrix `Yraw` (with column names = variable names),
then leave the rest of the script unchanged. The Minnesota prior mean
`prmean` should be set to `rep(1, M)` for I(1) data and `rep(0, M)` for
stationary data.

## Dependencies

```
MASS          # ginv
stochvol      # svsample2 for the common SV
mvtnorm       # dmvnorm
Matrix        # (only loaded transitively)
GIGrvg        # (only loaded transitively)
ggplot2       # plotting
reshape2      # melt() for the heatmap
```

All packages are on CRAN.

## Output

Running the driver produces three PDF figures under `figures/`:

* `sparsity_pattern.pdf` - heatmap comparing the true coefficient mask with
  the estimated sparsity pattern (a cell is "non-zero" iff the posterior
  median of that coefficient is non-zero after SAVS thresholding).
* `irf.pdf` - Cholesky impulse responses to a shock in the first variable,
  un-sparsified vs. sparsified, with 68 percent credible bands.
* `forecast_fan.pdf` - fan chart of recursive forecasts for the last
  variable, both versions of the model, with 68 and 90 percent bands.

The console also prints RMSE of the posterior median against the true
coefficient matrix for the simulated experiment.
