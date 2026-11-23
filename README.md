# Generalized Importance Sampling for Nested Simulation: code

R code for the numerical experiments (Section 5) of

> Chen, Q., & Feng, M. B. (2023). Generalized Importance Sampling for Nested Simulation.
> In *Proceedings of the 2023 Winter Simulation Conference* (pp. 409–420). IEEE.
> https://doi.org/10.1109/WSC60868.2023.10408062

The paper estimates the conditional expectations mu(S_tau,i) of a reverse iron butterfly
portfolio's P&L for m = 1,000 outer scenarios. It compares three nested estimators:
standard nested simulation (STD), generalized importance sampling with the equal mixture
(MLR), and generalized importance sampling with the NNLS weighted mixture (OIS,
Algorithm 1).

## Quick start

```
Rscript reproduce_paper.R          # from this folder; a few seconds
```

This redraws every figure of the paper into `figures/` and Table 1 into `tables/`, using
the saved summaries of the runs behind the paper (`results/published_runs/`). It then
checks Table 1 against the printed values (exact mode: every entry must round to the
printed digits).

To rerun the experiment with new random numbers, set `RERUN_SIMULATIONS <- TRUE` in
`reproduce_paper.R` (about 2 minutes on one core). Table 1 is then checked within Monte
Carlo error. The rerun alone:

```
Rscript scripts/run_experiments.R --seed 2023 --output results/rerun_seed2023
```

Requirements: R (tested with 4.3.3) and the packages `nnls` (1.5), `ggplot2` (3.4.4) and
`latex2exp` (0.9.6).

## Folder structure

| Path | Contents |
|---|---|
| `reproduce_paper.R` | Entry point: all figures and Table 1, then verification |
| `verify_against_paper.R` | Printed Table 1 values and the checks (exact / Monte Carlo) |
| `R/parameters.R` | Model and experiment settings; the 1,000 quantile scenarios |
| `R/option_portfolio.R` | Black-Scholes prices, simulation model g, true values mu(S_tau) |
| `R/estimators.R` | STD, MLR and OIS estimators; NNLS weights; the densities of Figure 3 |
| `R/utilities.R` | Summaries of the macro replications, reading/writing them, helpers |
| `R/figures_tables.R` | Figures 1–3 and Table 1 |
| `scripts/run_experiments.R` | The experiment: 200 macro replications of all estimators |
| `tools/extract_published_runs.R` | One-off: reduces the saved workspace of the published runs to `results/published_runs/` |
| `results/published_runs/` | Summaries of the runs behind the paper |
| `results/rerun_seed2023/`, `results/rerun_seed2024/` | Summaries of two reruns with new random numbers |
| `tests/test_against_original.R` | Rebuilt functions vs. the original code on the same random numbers |
| `original_2023/` | The original code of the published runs (reference for the tests) |

## Paper result -> code

| Paper | File written | Code |
|---|---|---|
| Figure 1a, simulation model g(S_T) | `figures/setting1.pdf` | `plot_figure1a()` |
| Figure 1b, mu(S_tau), outer pdf, STD bands (n = 100, 1,000) | `figures/setting2.pdf` | `plot_figure1b()` |
| Figure 2a, 95% bands of STD (Gamma = 10^6) and OIS (10^3) | `figures/STDOIS.pdf` | `plot_two_bands()` |
| Figure 2b, 95% bands of MLR (10^3) and OIS (10^3) | `figures/MLROIS.pdf` | `plot_two_bands()` |
| Figure 3, sampling densities q*, p-bar, q_beta* | `figures/densities.pdf` | `plot_figure3()`, `figure3_densities()` |
| Table 1, AMSEs and AMSE reductions | `tables/table1_amse.csv`, `.tex` | `table1_amse()` |

## How the simulations work

| Item | Setting (paper, Section 5) | Code |
|---|---|---|
| Stock | Black-Scholes, S_0 = 100, mu = 10%, r = 5%, sigma = 30% | `gis_parameters()` |
| Portfolio | Short put 125, long put 145, long call 145, short call 165; T = 1, tau = 0.5; P_0 = 17.32 | `iron_butterfly_price()` |
| Simulation model | g(S_T) = exp(-r(T - tau)) P(S_T) - P_0 exp(r tau) (P_0 carried forward to tau) | `simulation_output()` |
| Scenarios | m = 1,000 quantiles k/(m+1) of log S_tau | `outer_scenarios()` |
| True values | Black-Scholes price at tau minus P_0 exp(r tau) | `true_conditional_mean()` |
| STD | n = 1, 10, 100, 1,000 inner samples per scenario (Gamma = 10^3 to 10^6) | `std_estimates()` |
| MLR | Gamma = 1,000 samples of the equal mixture, components chosen at random | `mlr_estimates()` |
| OIS | Stage I: 100 samples of the equal mixture (components at random), NNLS weights; Stage II: round(900 beta_i) samples from component i (stratified) | `nnls_weights()`, `ois_estimates()` |
| Figure 3 weights | NNLS with 1,000 stratified samples (one per scenario) | `nnls_weights(stratified = TRUE)` |
| Macro replications | 200; 95% bands = 2.5% and 97.5% sample quantiles per scenario | `summarize_runs()` |
| Table 1 "AMSE" | average over scenarios of the sample variance of the 200 estimates | `table1_amse()` |

Because the Stage II counts are rounded, the Stage II sample size can differ from 900 by
one or two; the sampling pdf is the mixture with the rounded weights. Table 1 is also
written with the mean squared error including squared bias (`mse_with_bias`); bias is
negligible here.

## Notation

| Code | Paper |
|---|---|
| `log_s_tau`, `log_s_T` | log S_tau, log S_T (all densities are of log prices) |
| `inner_drift`, `inner_sd` | (r - sigma^2/2)(T - tau), sigma sqrt(T - tau) |
| `budget` | Gamma |
| `stage2_share` | Gamma_2 / Gamma |
| `weights` | beta* |
| `equal_mixture` | p-bar, eq. (8) |
| `mixture_density()` | q_beta, eq. (10) |
| `optimal_unnormalized()` | q* up to a constant, eq. (9) |
| `n_macro` | number of macro replications (200) |

## Reproducibility

The reruns use one L'Ecuyer-CMRG stream per macro replication (`--seed`), so a rerun
with the same seed gives identical summaries. Within a replication the estimators are
called in a fixed order (OIS, MLR, STD n = 1, 10, 100, 1,000). Every estimator draws its
random numbers in the same order as the original code: with the same seed, the rebuilt and
original functions give the same estimates (`tests/test_against_original.R`; differences
≤ 2e-16, from the order of floating-point operations only).

## Verification

* **Tests** (`Rscript tests/test_against_original.R`): 17 checks of the rebuilt functions
  against `original_2023/IronButterfly.R` on the same random numbers, and of the Figure 1
  and Figure 3 curves against those saved with the published runs. All pass.
* **Exact mode** (default `reproduce_paper.R`): Table 1 from the summaries of the
  published runs reproduces all 11 printed entries to the printed digits.
* **Figures**: the regenerated figures match the published ones; the 95% bands of
  Figure 2 extracted from the published figure files match the published-run summaries
  to within 2e-4.
* **Monte Carlo mode** (rerun, seeds 2023 and 2024): see `VERIFICATION.md`.

## Provenance

The runs behind the paper were made in April 2023 with `original_2023/IronButterfly.R`
(functions) and a plotting script (`original_2023/plots.R`, a later version of the one
that drew the figures). Their workspace (`WSC2023.RData`, 11 MB, kept with the archive
but not in git) holds the 1,000 x 200 matrices of estimates. The paper's figures and
Table 1 use the runs `result_FIT_1` (OIS), `result_MLR_1` (MLR) and
`result_STD_1/10/100/1000` (STD), and the weights `beta_fit` (Figure 3).
`tools/extract_published_runs.R` reduced them to `results/published_runs/`. The workspace
also holds runs with Gamma = 10^4 for OIS and MLR that the paper does not use. No seed was
recorded for the published runs, so they cannot be redrawn exactly; the reruns use new
random numbers.

`original_2023/` also keeps `IronButterfly2.R` (a variant with another portfolio, not used
in the paper) and `IronButterfly_2023-04-24.R` (a later copy that changes the
non-stratified branch of Stage II, which the published runs did not use).

## Run times

Laptop (Windows 11, R 4.3.3, one core): `reproduce_paper.R` from the summaries, about 5 s;
200 macro replications (`scripts/run_experiments.R`), about 2 minutes. The paper reports
no run times.

## Change log

* 2026-11: Rebuilt in plain R from the original code (Ben Feng, with Claude Code):
  vectorized estimators that draw the same random numbers as the original; one entry
  point; summaries of the published runs; tests against the original; reruns with two new
  seeds.
