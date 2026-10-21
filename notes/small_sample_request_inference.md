# Small-sample options for the request dimension: quick-check results

Branch `small-sample-request-inference`. These are quick checks (200 studies per condition), not full runs. All results describe the simulation's assumed model, not people.

## What was added

Two options in `analyze_study()`, both off by default. With the defaults, every verdict, estimate, SE and p-value is bit-identical to `main` (checked on 9 datasets × 3 cluster settings). The only change is two new columns, `lo95` and `hi95`.

- **`satterthwaite = TRUE`**: the analytic and jackknifed two-way tests use *t* reference distributions, with degrees of freedom estimated for each contrast (verdicts `two_way_sw` and `two_way_jk_sw`). The request part's df follow Pustejovsky and Tipton's model-assisted Satterthwaite formula, ν = (Σ p_i′Φp_i)² / ΣΣ (p_i′Φp_j)² (eq. 11, p. 13). The working model has independent reviews and assumes the fitted model is right, so each request's summed scores have its expected information as covariance. The jackknife version uses the CR3 form, which "closely approximates the jackknife" (p. 8). The two-way df are ν_req (V / V_req)², because only the request part is noisy. Page numbers refer to the arXiv v2 (2022) copy in the library.
- **`boot_draws > 0`**: a wild score bootstrap over requests, with Rademacher weights, for the analytic two-way *t* statistic (verdict `two_way_boot`; percentile-*t* intervals and p-values). Roodman et al. (2019) describe the score bootstrap for ML estimators (pp. 31–33). Their discussion of two-way clustering reports that bootstrapping in the dimension with fewer clusters works best (p. 20). They also warn that the plain score bootstrap uses bootstrap errors rather than residuals in the variance (pp. 32–33). So each draw here takes one Newton step and recentres the scores before the two-way variance is computed. That recentring is this code's own construction, not a published method. Draws with a non-positive two-way variance are dropped, as boottest does (p. 19).

Tests (`tests/testthat/test-small-sample.R`, 49 expectations; the suite now passes 106):

- The df formula matches the exact moments of the sandwich and jackknife variances in a Gaussian linear model to 1e-10.
- Expected information matches the expected negative Hessian, computed by finite differences of the analytic scores.
- The linearized request jackknife that the df rely on correlates above .99 with the refits, and its variance is within 10% of theirs.
- The bootstrap draws equal a direct computation from recentred review-level scores.
- `run_tests()` handles both t and bootstrap references.
- The options leave the default verdicts unchanged.

## Does coverage reach 95%?

Not for every contrast, with either option.

**Defaults, script 02** (`notes/small_sample_quick_checks/coverage_95.csv`). The pooled-by-contrast numbers are in `paired_comparisons.txt` (4,000 studies per contrast; Monte Carlo SE 0.34 points):

| Contrast | Jackknife (current) | Jackknife + Satterthwaite *t* | Score bootstrap |
|---|---|---|---|
| T1 | 92.6 | 97.7 | 95.5 |
| T2 | 91.9 | 94.0 | 93.4 |
| T3 plain | 91.6 | 93.6 | 94.2 |
| T3 tactic | 91.8 | 94.3 | 95.0 |
| T4 | 92.4 | 93.6 | 96.8 |
| T5 | 92.5 | 93.8 | 97.0 |

- **Ranges over true state × contrast** (1,000 studies each, Monte Carlo SE 0.69): jackknife 91.0–93.5 (the full run gave 91.4–93.0), jackknife + *t* 92.6–98.0, bootstrap 92.8–97.6.
- **Analytic two-way + *t*** reached only 89.7–97.5; *t* critical values can't fix standard errors that run 3–16% small (`se_ratio_two_way` 0.84–0.97 in the full run).
- **Median df** (jackknife version): T1 about 8 (12 unsafe requests), T3 about 16, T4 and T5 about 21, T2 about 23.
- **The T3 contrasts decide the falsifier.** At every sample size from 50 to 150 per cell, they covered 92.6–94.9 with jackknife + *t* and 93.6–95.3 with the bootstrap, against 90.3–92.7 with the jackknife (1,600 per cell; Monte Carlo SE 0.54).

**Sensitivity scenarios, script 06** (`sensitivity_coverage_by_test.csv`; 600 studies per scenario × contrast; Monte Carlo SE about 0.9):

- **Worst contrast per scenario:** jackknife 88.1–92.5, jackknife + *t* 91.5–94.3, bootstrap 92.0–95.2.
- **Mean over the six contrasts:** jackknife 90.0–93.9, jackknife + *t* 94.3–95.6, bootstrap 94.8–96.5.
- **Request-slope sweep, worst contrast** (`sensitivity_request_slopes_sweep.csv`): 90.8–92.5, 92.8–94.0, and 92.8–95.2.
- **τ_c = 0.8**, the jackknife's worst scenario: 88.1, 91.5, and 92.2.

**Why some contrasts still miss.**

- **The df themselves are about right.** In `df_calibration.csv`, the predicted request-part df match 2·mean(V)²/var(V) across 250 studies per scenario for T1 (11 vs 9.5–12) and T3–T5 (14–16 vs 12.5–17).
- **T2's df are too high once random effects are present** (35 vs 19–22). That leaves T2 at 94%.
- **The jackknife SEs still run about 5% small** for T3–T5 (`se_ratio_two_way_jk` in the full run). *t* critical values fix heavy tails, not that bias.
- **The bootstrap over-covers T4 and T5.**

In short: the bootstrap brings T1 and T3 to within a point of 95%. Satterthwaite *t* closes half to three-quarters of the gap for T2–T5 and overshoots on T1.

## What it costs in correct verdicts

- **Defaults, 100–150 per cell** (`correct_verdict_rates.csv`, paired changes against the jackknife in `paired_comparisons.txt`): no detectable cost. Changes were −1.0 to +2.0 points, with paired SEs of 1.5 or less.
  - At 125 per cell (cue / neither / vigilance / both): jackknife .915 / .910 / .930 / 1.00; + *t* .915 / .930 / .935 / 1.00; bootstrap .925 / .930 / .935 / 1.00.
  - The smallest size with at least .80 correct in every state stays 100 per cell for all three.
- **Small samples cost something.** At 75 per cell, −3.0 to +1.5 points. At 50 per cell, cue −8.5 (+ *t*) and −10.5 (bootstrap); neither −5.0 and −7.0.
- **Request slopes, 100–200 per cell** (script 06, `sensitivity_summary.csv`): cue and neither verdicts fall 2.5–7.0 points with + *t* and 1.0–6.5 with the bootstrap. At 125 per cell, cue .655 becomes .615 and .610; neither .665 becomes .635 and .655.
- **Other scenarios:** changes ran from −5.5 to +5.5 points. The bootstrap lost 5.5 for neither with sensitivity SDs of 0.5, and both options gained about 5 for neither at τ_c = 0.8. Vigilance verdicts never dropped.

**Equivalence-test size at the SESOI** (`sensitivity_equivalence_at_sesoi.csv`; 200 studies, Monte Carlo SE 1.5 points) is not resolved by a quick check:

| Errors | T3 plain | T3 tactic |
|---|---|---|
| Jackknife | 4.0 | 5.5 |
| Jackknife + *t* | 3.5 | 4.5 |
| Bootstrap | 2.0 | 4.0 |

- The full run gave 5.8 and 7.1 for the jackknife.
- *t* changed the call in only 3 of 400 tests.
- With the wider bound from a pilot that puts τ_c at 0.3, the jackknife's 9.0 and 8.5 became 7.0 and 7.0.

## Proposed full runs

Run scripts 02 and 06 from a private copy, then copy `output/tables` back:

- `SIM_FULL=1 SIM_WORKERS=16 Rscript scripts/02_power_and_adjudication.R`
- `SIM_FULL=1 SIM_WORKERS=16 Rscript scripts/06_sensitivity.R`

These use the default 999 draws (`SIM_BOOT=999`). Seeds and default computations are unchanged, so the existing columns should reproduce the committed tables exactly, with the new columns added.

**Runtime.**

- **Overhead per study** (single process, `analyze_study` with both options vs default): +7% at 50 per cell, +27% at 125, +58% at 200, and +37% with 96 requests. The bootstrap is almost all of it; Satterthwaite costs nothing measurable.
- **Estimated full runs:** from the README's last full runs (139 and 219 minutes on 16 workers), about **2.8 h for script 02 and 4.7 h for script 06** on 16 workers.
- **Measured quick runs:** script 02, 4,000 studies in 28.9 minutes; script 06, 8,400 in 107 minutes. Both ran on 20 workers while 16 brms chains from another session shared the machine. At that throughput the full runs would take about 2.4 h and 9 h.
- **`SIM_BOOT=399`** would cut the bootstrap's share by about 60%.

## Caveats

- **This is not a CR2 correction.** Satterthwaite df are applied to the existing jackknife and analytic sandwiches.
- **The working-model df ignore random effects.** T2 shows the consequence.
- **The recentred score bootstrap is unpublished as such.** Roodman et al. (p. 33) advise against relying on the score bootstrap without evidence that it works in the case at hand. These simulations are that evidence for this design only.
- **Verdict rates at R = 200 have Monte Carlo SEs of 1.5–3.5 points.** Paired differences are more precise but still about 1–2 points.
- **The README's line that jackknife-type corrections should "err on the conservative side"** doesn't match its own coverage results (91–93%). It is left as is here.

## References

Pustejovsky, J. E., & Tipton, E. (2018). Small-sample methods for cluster-robust variance estimation and hypothesis testing in fixed effects models. *Journal of Business & Economic Statistics, 36*(4), 672–683. https://doi.org/10.1080/07350015.2016.1247004. Pages cited from arXiv:1601.01981v2.

Roodman, D., MacKinnon, J. G., Nielsen, M. Ø., & Webb, M. D. (2019). Fast and wild: Bootstrap inference in Stata using boottest. *The Stata Journal, 19*(1), 4–60. https://doi.org/10.1177/1536867X19830877
