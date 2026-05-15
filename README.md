# Replication Archive — Bayesian VARs and Nonparametric Macro Models

Self-contained R replication packages for eleven Bayesian time-series models.
Each subfolder is independent: it ships its own estimator, a synthetic
data-generating process in `simulate_data.R`, a driver in
`main_replication.R`, and a per-folder `README.md` with the model, file
table, dependencies, runtime guidance, and outputs.

Every driver runs end-to-end on simulated data and overlays the **true** DGP
quantities (IRFs, parameter paths, regime indicators, factor dimension,
realised future y) on every relevant figure, so each package can be sanity-
checked without external data.

## Packages

| Folder | Model |
|---|---|
| [`GPVAR_replication`](GPVAR_replication/) | Gaussian-process VAR with stochastic volatility and horseshoe-shrunk structural covariance |
| [`NGVAR_replication`](NGVAR_replication/) | Bayesian VAR with hierarchical Normal-Gamma shrinkage and SV (Huber & Feldkircher, 2019, JBES) |
| [`PVAR_IRGA_replication`](PVAR_IRGA_replication/) | Panel VAR with the Integrated Rotated Gaussian Approximation |
| [`SparseVAR_replication`](SparseVAR_replication/) | Conjugate Minnesota VAR with lag-wise SAVS sparsification + SV |
| [`SubspaceVAR_replication`](SubspaceVAR_replication/) | Subspace-shrinkage conjugate VAR (Huber & Koop, JAE) |
| [`TBVECM_replication`](TBVECM_replication/) | Three-regime threshold VECM |
| [`ThreshTVP_replication`](ThreshTVP_replication/) | Latent-threshold TVP-VAR with SV (Huber, Kastner & Feldkircher, JAE) |
| [`TVPSVD_replication`](TVPSVD_replication/) | TVP-WN-SVD regression with sparse-mixture g-prior (Hauzenberger, Huber, Koop & Onorante, 2021) |
| [`BNPVAR_replication`](BNPVAR_replication/) | Bayesian VAR with Dirichlet-process mixture on the structural shocks |
| [`mixBART_replication`](mixBART_replication/) | Bayesian VAR with mixture-of-BART conditional means and factor SV |
| [`SubspaceInflation_replication`](SubspaceInflation_replication/) | Single-equation inflation forecasting with GP-subspace / BART / UCSV |

## How to run a single package

```sh
cd <folder>_replication
Rscript main_replication.R
```

Figures and the saved posterior object are written to the folder's own
`figures/` and `output/` (or top-level for the older packages).

## Truth-overlay convention

Every plot that shows IRFs, time-varying parameters, regime probabilities,
factor dimension, or out-of-sample forecasts overlays the true value from
the data-generating process — typically in **red dashed**. This makes each
package self-checking on simulated data.

## Dependencies

The R-package dependency list is per-folder (see each `README.md`). Across
the archive the recurring CRAN packages are `MASS`, `mvtnorm`, `Matrix`,
`stochvol`, `GIGrvg`, `ggplot2`, `reshape2`, `dbarts`, `bayesm`, and
`factorstochvol`. One package (`ThreshTVP_replication`) additionally
requires the local `threshtvp` package shipped in its folder.

## Using your own data

Each `main_replication.R` contains a clearly marked `REPLACE WITH YOUR
DATA` block. Replace the call to `simulate_*()` with code that reads your
own `(T x M)` matrix and the rest of the pipeline runs unchanged.
