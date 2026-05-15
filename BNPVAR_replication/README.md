# BNP-VAR Replication Package

Replication code for the Bayesian VAR with a **Bayesian nonparametric
(Dirichlet-process / sticky-breaking) mixture on the structural shocks**.
The shock distribution is left unspecified and learned along with the
coefficients; identification of the structural rotation is recursive
(Cholesky / eigen) *within* each mixture component, so the IRFs naturally
vary across regimes.

## Model in one paragraph

For a VAR(p) of dimension M,

    y_t = c + A_1 y_{t-1} + ... + A_p y_{t-p} + e_t,
    e_t ~ sum_{g=1}^{G} eta_g  N( mu_g, Sigma_g ),

with the mixture weights `eta_g` drawn from a stick-breaking process, the
component means `mu_g` and covariances `Sigma_g` updated jointly with a
Normal-Wishart prior, and the regime indicators `S_t in {1,...,G}` sampled
by a slice sampler that truncates the stick at a data-driven `G+`. The
concentration parameter `alpha` is updated by a random-walk Metropolis step
on the log scale. Optionally each equation may carry idiosyncratic
stochastic volatility (`SV.idi = TRUE`).

The structural rotation is computed by Cholesky-factoring each posterior
draw of `Sigma_g`, so impulse responses are *cluster-specific*; the
"composite" IRF reported in the paper is the `eta_g`-weighted average of
the within-cluster IRFs.

## Files

```
BNPVAR_replication/
  README.md
  simulate_data.R                       synthetic mixture-shock VAR generator
  main_replication.R                    driver: simulate -> fit -> figures -> save
  Codes_and_data/
    aux.R                               the sampler: `bvar.mix(Xraw, ...)`
    IRF.script.R                        pre-existing empirical IRF driver
                                        (kept; uses real FRED-MD data)
    Data_QD/
      prepare_dat.R                     downloads the latest FRED-QD vintage
      US_QDvariables.csv                variable info / info sets
      dat.RData                         frozen cached snapshot
  figures/                              generated plots
  output/                               generated `bnpvar_fit.rds`
```

## Dependencies

```
MASS        mvtnorm     mvnfast     stochvol
GIGrvg      bayesm      ggplot2     reshape2
```

Install with

```r
install.packages(c("MASS","mvtnorm","mvnfast","stochvol","GIGrvg","bayesm",
                   "ggplot2","reshape2"))
```

The empirical driver `Codes_and_data/IRF.script.R` additionally needs
`httr`, `ustyc`, `dplyr`, `quantmod`, `lubridate`, `readxl`, `zoo`, and
`gridExtra` (for the FRED-MD pull and the multi-panel figures).

## How to run

```sh
cd BNPVAR_replication
Rscript main_replication.R
```

This fits the BNP-VAR on a small simulated mixture-of-shocks VAR and writes
posterior summaries plus three figures into `figures/`. To use your own
data instead, open `main_replication.R`, find the line marked
`REPLACE WITH YOUR DATA`, and assign a numeric `(T x M)` matrix to `Xraw`.

To reproduce the empirical exercise from the paper (FRED-QD US panel),
use `Codes_and_data/IRF.script.R` instead — it pulls the latest vintage of
the FRED-QD data, fits the model with `G.max = 20`, and writes the full set
of regime-specific IRFs to `Insample_results/`.

## Runtime warning

Each Gibbs sweep does, **per cluster**, a Wishart draw on an `M x M` scale
matrix and, **per period**, a slice sampler over the regime indicator.
Adding the per-cluster IRF computation (`O(M^3 * irf.hor)`), total cost is

    O( (nburn + nsave) * (G.max * irf.hor * M^3 + T * G.max) )

per chain. With `T = 200`, `M = 5`, `G.max = 10`, `irf.hor = 24`,
`nsave = nburn = 500`, `main_replication.R` finishes in a few minutes on a
laptop. The paper uses `nsave = nburn = 5000`; budget several hours.

For exploration:

* keep `G.max <= 10` and `irf.hor <= 24`;
* set `SV.idi = FALSE`;
* set `sample.w = "conditional"` (the cheaper branch).

## Outputs

After `Rscript main_replication.R`:

* `figures/cluster_count.{pdf,png}` — posterior of `G+`, the number of
  occupied mixture components.
* `figures/cluster_path.{pdf,png}` — `Pr(S_t = g | data)` over `t`, for
  every cluster with non-trivial mass.
* `figures/irf_cluster1.{pdf,png}` — Cholesky impulse responses to a
  one-s.d. shock in the last variable, computed under cluster 1
  (re-ordered by mass).
* `output/bnpvar_fit.rds` — the full posterior list returned by `bvar.mix`:
  `A_coefs`, `mu_G`, `Sigma_G`, `G`, `S` (re-ordered), `IRF`, `eta`, `fit`,
  `vola`, and the predictive draws.

## What is *not* included

* The cluster-job scaffolding (`SGE_TASK_ID`) used in the paper.
* The forecasting / out-of-sample loop (`bvar.mix` does produce predictive
  draws via the `fhorz` argument; only the in-sample plots are wired up).
* The macro database snapshot used in the published runs — the empirical
  driver re-downloads the latest FRED-QD vintage at runtime.
