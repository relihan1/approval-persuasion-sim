# Unsafe-request prevalence: 12.5% against 25%

The design shows each participant 12 unsafe requests out of 48 (25%), far above the rate at which unsafe requests reach approvers in deployment. This check reruns the cue, vigilance, and neither states with 6 unsafe out of 48 (12.5%; 42 safe), at 125 participants per cell, with every other value at its default (`scripts/07_prevalence_check.R`). The 25% comparison is script 06's default scenario: the same sample size, the same errors, and 1,000 studies per state. All results describe the simulation's assumed model, not people.

The run used 1,100 studies per state, so the 95% Monte Carlo half-width of a verdict rate is at most 3.0 points (1.8 at a rate near .90). It took 18.3 minutes on 12 workers.

With 6 unsafe requests, each paltering stratum has 3, so the two counterbalancing lists can't split them evenly. List A shows 4 of the 6 with tactics and list B shows 2. Across participants each request still appears equally often in each version.

## Correct verdicts

`prevalence_correct_verdicts.csv`; the 25% rows are from `sensitivity_summary.csv` (scenario `default`).

| True state | Unsafe share | Participant-clustered | Two-way, jackknife | Two-way, jackknife + *t* | Two-way, score bootstrap |
|---|---|---|---|---|---|
| Cue | 25% | .93 | .89 | .90 | .91 |
| Cue | 12.5% | .69 | .62 | .47 | .48 |
| Vigilance | 25% | .95 | .93 | .94 | .94 |
| Vigilance | 12.5% | .94 | .94 | .94 | .93 |
| Neither | 25% | .95 | .89 | .90 | .91 |
| Neither | 12.5% | .72 | .62 | .47 | .48 |

Halving the unsafe requests left the vigilance verdicts where they were and cost the cue and neither verdicts 23 to 43 points. Under the recommended errors (the score bootstrap), cue fell from .91 to .48 and neither from .91 to .48.

## Why: the equivalence tests lose power

The cue and neither verdicts need both T3 contrasts to be declared equivalent to zero. Each T3 contrast is a change in *d*, and *d* is estimated from the unsafe reviews, so halving them makes T3 less precise while the equivalence bound stays at 0.173. With 12.5% unsafe, each T3 contrast was declared equivalent in 83% to 85% of cue and neither studies under participant-clustered errors, and in 68% to 71% under the score bootstrap (`prevalence_test_rates.csv`). Requiring both brings the verdict rates down to .69 to .72 and to .48.

The vigilance verdict needs T3 gains, and the gain (0.26 on the estimator's scale; `prevalence_pseudo_true.csv`) is large enough to detect at either prevalence. T3 rejected in more than 99% of vigilance studies with participant-clustered errors, and in 92% with the bootstrap.

The two-way errors lose more than the participant-clustered ones. Under the defaults, requests don't differ in how tactics and training affect them. Both training groups review the same requests, so request effects cancel in T3, and the participant-clustered intervals are valid and the narrowest. The two-way intervals also count variation between requests, which on the unsafe side now rests on 6 requests. At 12.5% the bootstrap's T3 intervals over-covered (96.3% and 95.8%, against 95.5% and 95.1% for participant clustering), and wider intervals mean fewer equivalence calls.

The threat check, T1, compares tactic-laden with plain unsafe requests, so it also rests on the unsafe reviews. With the bootstrap it rejected in 84% to 94% of studies, against essentially all studies with participant-clustered errors. T1 doesn't enter the verdict, but a missed threat check would be reported as a limit of prompted attacks.

## Coverage of 95% intervals, T1 to T3

Pooled over the three states (`prevalence_coverage.csv`; Monte Carlo SE about 0.4 points).

| Contrast | Unsafe share | Participant-clustered | Two-way, jackknife | Two-way, jackknife + *t* | Two-way, score bootstrap |
|---|---|---|---|---|---|
| T1 | 25% | 95.5 | 91.6 | 96.8 | 94.3 |
| T1 | 12.5% | 95.8 | 90.8 | 98.1 | 94.7 |
| T2 | 25% | 95.6 | 93.0 | 94.5 | 94.4 |
| T2 | 12.5% | 95.5 | 92.9 | 94.3 | 94.2 |
| T3 plain | 25% | 95.3 | 92.5 | 94.4 | 95.2 |
| T3 plain | 12.5% | 95.5 | 91.0 | 96.4 | 96.3 |
| T3 tactic | 25% | 95.2 | 92.4 | 94.7 | 95.4 |
| T3 tactic | 12.5% | 95.1 | 90.8 | 95.8 | 95.8 |

Coverage held up. Participant-clustered and bootstrap intervals stayed at 94.2% to 96.3%, but the bootstrap and jackknife + *t* intervals for T1 and T3 grew more conservative. The jackknife without *t* undercovered at either prevalence, and by 0.8 to 1.6 points more at 12.5% for T1 and T3 (90.8% to 91.0%).

## What this means for the design

- **Lower prevalence costs precision on the contrasts that decide the falsifier, not validity.** The intervals stay honest, but the study can no longer show that training left *d* unchanged as often as it should.
- **Not simulated: keep 12 unsafe requests and add safe ones.** For example, 12 unsafe out of 96 is also 12.5% prevalence, keeps the unsafe reviews behind T1 and T3, and adds safe reviews. An inference not yet simulated: this would recover most of the 25% design's verdict rates at the cost of a longer session. Script 06's 96-request scenario uses 24 unsafe (25%) and is a different design.
- **Also not simulated: more participants.** If the T3 standard errors scale with one over the square root of the number of unsafe reviews, as they roughly should when the unsafe side dominates, matching the 25% design's precision with 6 unsafe requests would take about twice as many participants. This is also an inference.
- **Prevalence may change behavior.** Prevalence could also change reviewers' criterion and attention. The simulation holds the response model fixed across prevalence, so it says nothing about that.

## Files

- `scripts/07_prevalence_check.R`: the simulation; its per-study analysis mirrors script 06's `one_rep()`.
- `output/tables/prevalence_simulations.rds`: one row per study and contrast.
- `output/tables/prevalence_correct_verdicts.csv`, `prevalence_coverage.csv`, `prevalence_test_rates.csv`, `prevalence_pseudo_true.csv`: the tables above.
