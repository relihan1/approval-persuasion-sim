# Small-sample options for the request dimension: full-run results

Branch `small-sample-request-inference`. The numbers come from full runs of scripts 02 and 06 (`SIM_FULL=1`: 1,000 studies per condition). They are in `output/tables`, with extra summaries in `notes/small_sample_checks/`. All results describe the simulation's assumed model, not people.

## What was added

Two options in `analyze_study()`, both off by default. With the defaults, every verdict, estimate, SE and p-value is unchanged. The full runs reproduced every existing column of every table that scripts 02 and 06 write, and the full simulation files, exactly. The only change is new columns.

- **`satterthwaite = TRUE`**: the analytic and jackknifed two-way tests use *t* reference distributions, with degrees of freedom estimated for each contrast (verdicts `two_way_sw` and `two_way_jk_sw`). The request part's df follow Pustejovsky and Tipton's model-assisted Satterthwaite formula, ν = (Σ p_i′Φp_i)² / ΣΣ (p_i′Φp_j)² (eq. 11, p. 13). The working model has independent reviews and assumes the fitted model is right, so each request's summed scores have its expected information as covariance. The jackknife version uses the CR3 form, which "closely approximates the jackknife" (p. 8). The two-way df are ν_req (V / V_req)², because only the request part is noisy. Page numbers refer to the arXiv v2 (2022) copy in the library.
- **`boot_draws > 0`**: a wild score bootstrap over requests, with Rademacher weights, for the analytic two-way *t* statistic (verdict `two_way_boot`; percentile-*t* intervals and p-values). Roodman et al. (2019) describe the score bootstrap for ML estimators (pp. 31–33). Their discussion of two-way clustering reports that bootstrapping in the dimension with fewer clusters works best (p. 20). They also warn that the plain score bootstrap uses bootstrap errors rather than residuals in the variance (pp. 32–33). So each draw here takes one Newton step and recentres the scores before the two-way variance is computed. In a linear model the recentred scores are exactly the scores of the bootstrap residuals, so the procedure reduces to the unrestricted wild cluster bootstrap, WCU (pp. 10–11). Carrying it to an ordinal model through each review's expected information is this code's own step, not a published method. Roodman et al. prefer the restricted variant, WCR, in general (p. 13); WCU was used here because one set of draws gives intervals for all six contrasts. Draws with a non-positive two-way variance are dropped, as boottest does (p. 19).

Tests (`tests/testthat/test-small-sample.R`, 49 expectations; the full suite now passes 134):

- The df formula matches the exact moments of the sandwich and jackknife variances in a Gaussian linear model to 1e-10.
- Expected information matches the expected negative Hessian, computed by finite differences of the analytic scores.
- The linearized request jackknife that the df rely on correlates above .99 with the refits, and its variance is within 10% of theirs.
- The bootstrap draws equal a direct computation from recentred review-level scores.
- `run_tests()` handles both t and bootstrap references.
- The options leave the default verdicts unchanged.

## Does coverage reach 95%?

**With the bootstrap, it does for T1 and for the two T3 contrasts that decide the falsifier, under the defaults.** It doesn't for every contrast or scenario, with either option.

**Defaults, script 02** (`coverage_95.csv`, 5,000 studies per true state × contrast, Monte Carlo SE 0.31). The table pools over the four true states (20,000 studies per contrast; SE 0.15; `paired_comparisons.txt`):

| Contrast | Jackknife (current) | Jackknife + Satterthwaite *t* | Score bootstrap |
|---|---|---|---|
| T1 | 91.8 | 97.4 | 94.8 |
| T2 | 92.5 | 94.2 | 94.0 |
| T3 plain | 91.9 | 94.0 | 94.7 |
| T3 tactic | 92.5 | 94.5 | 95.2 |
| T4 | 92.7 | 94.0 | 97.0 |
| T5 | 92.2 | 93.8 | 97.2 |

- **Ranges over true state × contrast:** jackknife 91.4–93.0, jackknife + *t* 93.3–97.5, bootstrap 93.8–97.3.
- **Analytic two-way + *t*** reached only 90.4–96.7; *t* critical values can't fix standard errors that run 3–16% small (`se_ratio_two_way` 0.84–0.97).
- **T3 by sample size** (8,000 studies per size; SE 0.24): at every size from 50 to 150 per cell, the T3 contrasts covered 94.8–95.1 with the bootstrap, 94.1–94.4 with jackknife + *t*, and 91.9–92.6 with the jackknife. Across the four true states, the bootstrap's T3 coverage was 94.5–95.4.
- **Median df** (jackknife version): T1 8, T2 23, T3 16 to 17, T4 and T5 21.

**Sensitivity scenarios, script 06** (`sensitivity_coverage_by_test.csv`; 3,000 studies per scenario × contrast; Monte Carlo SE 0.4):

- **Worst contrast per scenario:** jackknife 88.3–92.9, jackknife + *t* 92.0–94.5, bootstrap 92.9–94.8.
  - With the jackknife, the worst contrast was usually T1.
  - With jackknife + *t*, it was T4 or T5 in nine scenarios and a T3 contrast in five.
  - With the bootstrap, it was T2 in seven scenarios and T1 in five.
- **Mean over the six contrasts:** jackknife 89.9–93.6, jackknife + *t* 93.9–95.1, bootstrap 95.2–96.3.
- **Request-slope sweep, worst contrast** (`sensitivity_request_slopes_sweep.csv`): 91.8–92.9, 93.7–94.3, and 93.7–94.8.
- **The hardest scenarios:**
  - τ_c = 0.8: 88.3, 92.0, and 92.9.
  - σ = 1.0: 89.3, 92.8, and 93.0.

**Why some contrasts still miss.**

- **The df themselves are about right.** In `df_calibration.csv`, the predicted request-part df match 2·mean(V)²/var(V) across 250 studies per scenario for T1 (11 vs 9.5–12) and T3–T5 (14–16 vs 12.5–17).
- **T2's df are too high once random effects are present.** The test uses the total df: 23.7 predicted against 12.9 empirical under the defaults, and 28.1 against 18.0 with request slopes. For the request part alone it is 35 against 19 to 22. With independent reviews the two agree (34.5 against 36.8). That leaves T2 at 94%.
- **The jackknife SEs still run about 5% small** for T3–T5 (`se_ratio_two_way_jk`). *t* critical values fix heavy tails, not that bias.
- **The bootstrap over-covers T4 and T5.**

## What it costs in correct verdicts

- **Defaults, 125–150 per cell** (`correct_verdict_rates.csv`, 1,000 studies per cell; paired changes against the jackknife in `paired_comparisons.txt`): nothing. Both options gained 0.5–2.4 points (paired SE 0.6 or less).
  - At 125 per cell (cue / neither / vigilance / both): jackknife .896 / .904 / .940 / 1.00; + *t* .903 / .916 / .945 / 1.00; bootstrap .905 / .917 / .947 / 1.00.
- **At 100 per cell:** −1.2 to +0.7 points (SE 0.7 or less). Cue is .860 with the jackknife, .848 with + *t* and .849 with the bootstrap.
- **Smaller samples cost more.** At 75 per cell, cue and neither lose 3.1–4.6 points. At 50 per cell they lose 5.1–8.0. The smallest size with at least .80 correct in every state stays 100 per cell for all three.
- **Request slopes, 100–200 per cell** (`sensitivity_summary.csv`): cue and neither verdicts fall 2.2–5.9 points with + *t* and 1.2–5.1 with the bootstrap. At 125 per cell, cue .631 becomes .589 and .599; neither .650 becomes .606 and .617.
- **Larger sensitivity SDs (0.5):** + *t* loses 2.2 and 0.7 points (cue, neither); the bootstrap loses 3.7 and 3.9.
- **σ = 1.5:** losses of 1.1–1.9 points.
- **Every other scenario gained:** 0.0 to 4.0 points with + *t*, 0.3 to 5.7 with the bootstrap. The largest gain was for neither at τ_c = 0.8.
- **Vigilance verdicts never dropped** (+0.1 to +2.2).

**Equivalence-test size at the SESOI** (`sensitivity_equivalence_at_sesoi.csv`; 1,000 studies; Monte Carlo SE 0.7 points; nominal at most 5%):

| Errors | T3 plain | T3 tactic |
|---|---|---|
| Jackknife | 5.8 | 7.1 |
| Jackknife + *t* | 4.6 | 5.4 |
| Bootstrap | 3.8 | 4.9 |
| Participant-clustered | 3.9 | 6.5 |
| Analytic two-way | 7.3 | 8.2 |
| Analytic two-way + *t* | 5.8 | 6.8 |

- *t* turned 29 of the jackknife's 129 equivalence calls into non-equivalence.
- With the wider bound from a pilot that puts τ_c at 0.3, the jackknife's 9.1 and 10.6 became 8.2 and 9.5 with *t*. With the 0.20 unconverted, 16.5 and 17.4 became 14.0 and 15.5. The bootstrap wasn't run with these bounds, because its draws aren't kept.

**If one option is chosen, the bootstrap fits the preregistration's needs best.**

- It is the only version whose T3 intervals cover at the nominal 95% under the defaults (94.5–95.4). Its T3 equivalence tests are at or below 5% (3.8 and 4.9); jackknife + *t* comes close (4.6 and 5.4, the latter within Monte Carlo error of 5%), but its T3 intervals cover 94.0–94.5.
- It costs no correct verdicts at 125 per cell.
- Its weak spots are T2 (94.0%) and over-wide T4 and T5 intervals.

Jackknife + Satterthwaite *t* is the simpler choice, built on a published formula. It is about a point short on T2–T5 and leaves the tactic-laden T3 equivalence test at 5.4%.

## Runtime and checks

- **Full runs on 16 workers**, while brms fits from another session shared the machine (12 chains when the runs started): script 02 took 138.7 minutes and script 06 took 204.0. That is no slower than the README's previous full runs (139 and 219).
- **Overhead per study** (single process): +7% at 50 per cell, +27% at 125, +58% at 200, and +37% with 96 requests. The bootstrap is almost all of it; `SIM_BOOT=399` would cut its share by about 60%.
- **Missing intervals:** a negative two-way variance leaves no interval in 222 of 120,000 bootstrap rows and 142 of the jackknife + *t* rows in script 02. They are excluded from coverage, as before.
- **Earlier quick checks** (200 studies per condition; commit 4fded76) agreed with these results within their Monte Carlo error.

## Caveats

- **This is not a CR2 correction.** Satterthwaite df are applied to the existing jackknife and analytic sandwiches.
- **The working-model df ignore random effects.** T2 shows the consequence.
- **The recentred score bootstrap is unpublished as such.** It reduces to WCU in a linear model, but not exactly in this one. Roodman et al. (p. 33) advise against relying on the score bootstrap without evidence that it works in the case at hand. These simulations are that evidence for this design and its assumed model only.
- **Rademacher weights.** Roodman et al. prefer Webb's six-point weights only below about 10 to 12 clusters (p. 15). With 48 requests, Rademacher weights allow 2^48 distinct draws. The df above put the request dimension's information at about 8 to 20 clusters, though, so six-point weights may be worth a check. They weren't tried.
- **The README's line that jackknife-type corrections should "err on the conservative side"** didn't match its own coverage results (91–93%) and has since been removed.

## References

Pustejovsky, J. E., & Tipton, E. (2018). Small-sample methods for cluster-robust variance estimation and hypothesis testing in fixed effects models. *Journal of Business & Economic Statistics, 36*(4), 672–683. https://doi.org/10.1080/07350015.2016.1247004. Pages cited from arXiv:1601.01981v2.

Roodman, D., MacKinnon, J. G., Nielsen, M. Ø., & Webb, M. D. (2019). Fast and wild: Bootstrap inference in Stata using boottest. *The Stata Journal, 19*(1), 4–60. https://doi.org/10.1177/1536867X19830877
