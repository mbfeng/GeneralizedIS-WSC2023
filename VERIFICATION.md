# Verification against the published paper

Checked November 2026 (R 4.3.3, Windows laptop, one core).

## 1. Rebuilt functions vs. the original code

`Rscript tests/test_against_original.R`. The tests compare the rebuilt functions with
`original_2023/IronButterfly.R` on the same random numbers, and the Figure 1 and 3 curves
with those saved alongside the published runs. All 17 checks pass. The largest difference
is 2.2e-16, which comes only from the order of floating-point operations.

| Quantity | Max abs difference |
|---|---|
| P_0, scenarios, g(S_T), mu(S_tau) | 0 |
| STD (n = 1, 10, 100) | ≤ 8.7e-19 |
| MLR (Gamma = 10^3) | 3.5e-18 |
| NNLS weights (random and stratified Stage I) | 0 |
| OIS (Gamma = 10^3, three seeds) | ≤ 2.2e-16 |
| Figure 3 densities q*, p-bar, q_beta* | 0 |

## 2. Exact mode: Table 1 from the published runs

`Rscript reproduce_paper.R` (default). Every entry rounds to the printed digits.

| Estimator | Printed AMSE | Computed | Printed reduction | Computed |
|---|---|---|---|---|
| OIS (10^3) | 0.0167 | 0.016695 | 1 | 1 |
| MLR (10^3) | 0.0339 | 0.033913 | 2.03 | 2.0314 |
| STD (10^3) | 18.59 | 18.594 | 1,114 | 1113.8 |
| STD (10^4) | 1.84 | 1.8382 | 110 | 110.11 |
| STD (10^5) | 0.18 | 0.18471 | 11.06 | 11.064 |
| STD (10^6) | 0.018 | 0.018427 | 1.10 | 1.1038 |

## 3. Figures

The figures redrawn from the published-run summaries match the published figures when
rendered side by side. The Figure 2 band polygons, extracted from the published figure
files, match the summaries' 2.5% and 97.5% quantiles to within 2e-4.

## 4. Monte Carlo mode: reruns with new random numbers

`scripts/run_experiments.R` with seeds 2023 and 2024 (200 macro replications each, about
2 minutes each). Each Table 1 entry is checked with
z = (x - v) / sqrt(2 se^2 + h^2/3). In each rerun 0 of 11 entries have |z| > 2.5; about
0.14 would be expected by chance.

| Entry | Printed | Seed 2023 (z) | Seed 2024 (z) |
|---|---|---|---|
| OIS AMSE | 0.0167 | 0.0199 (1.04) | 0.0205 (1.58) |
| MLR AMSE | 0.0339 | 0.0348 (0.18) | 0.0370 (0.64) |
| MLR reduction | 2.03 | 1.75 (-0.77) | 1.81 (-0.69) |
| STD 10^3 AMSE | 18.59 | 18.65 (0.40) | 18.40 (-1.40) |
| STD 10^3 reduction | 1,114 | 936 (-1.22) | 900 (-2.05) |
| STD 10^4 AMSE | 1.84 | 1.842 (0.20) | 1.848 (0.72) |
| STD 10^4 reduction | 110 | 92.5 (-1.22) | 90.3 (-1.87) |
| STD 10^5 AMSE | 0.18 | 0.1849 (1.59) | 0.1852 (1.68) |
| STD 10^5 reduction | 11.06 | 9.28 (-1.23) | 9.05 (-1.91) |
| STD 10^6 AMSE | 0.018 | 0.01844 (1.41) | 0.01826 (0.86) |
| STD 10^6 reduction | 1.10 | 0.926 (-1.21) | 0.893 (-2.00) |

The OIS and MLR squared errors have mildly heavy tails. The worst 1% of macro
replications contribute 7-12% of their AMSE, against about 1% for STD. Their standard
errors are therefore less precise than those of STD. The reductions all share the OIS
AMSE as denominator, so their z-values move together.

A rerun with the same seed gives byte-identical summaries (checked for seed 2023).
