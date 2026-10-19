# approval-persuasion-sim

**A simulation-based design analysis for a study of whether AI agents can argue their way past human approval**

Dr. Daniel P. Relihan · R · MIT license

Coding agents ask people to approve risky actions, such as reading a credential, merging a change, or silencing an alert. Anthropic keeps its agent's own messages away from the classifier that screens Claude Code's riskier actions, "so the agent can't talk the classifier into making a bad call," and reports that "Claude Code users approve 93% of permission prompts" (Hughes, 2026). Levy et al. (2026) ask the human version of that question. They propose training overseers and paraphrasing AI messages as defenses, and suggest human experiments as the next step.

This repository simulates such an experiment before anyone runs it. It writes down the data-generating process I expect, fits the analysis I would preregister, and asks four questions:

1. Can the analysis tell apart two accounts of what prebunking training does?
2. Are its error rates and interval coverage what they claim to be?
3. How many participants does it need, and how sensitive is that number to the assumptions?
4. Can simulated reviewers stand in for some of the people, and can a classifier check that the stimuli don't give the answer away?

Every number below describes the assumed model, not people.

## The study being simulated

| Element | Levels |
|---|---|
| Request (within participants) | 36 safe and 12 unsafe, in random order |
| Justification (within, counterbalanced) | Plain, or with persuasion tactics. Each participant sees each request in one version, with tactics in half of each stratum |
| Training (between) | A prebunking module on the tactics, or a time-matched control |
| Screen (between) | The agent's argument, a neutral paraphrase, or the argument plus independently retrieved evidence |
| Response | Approve or reject, with confidence on a 6-point scale (1 = definitely approve, 6 = definitely reject) |

## Two accounts of what training does

Prebunking warns people about a manipulation tactic and shows them a weakened example. Simchon et al. (2026) reanalyzed 33 experiments with a preregistered hierarchical ordinal model. They found that prebunking raised discrimination, with a median gain in *d*′ of 0.31 between participants (95% ETI [0.21, 0.42]), and near-zero bias differences. A reanalysis of prebunking games, using comparable true and fake news, found more "false" responses to all items instead (Modirrousta-Galian & Higham, 2023).

My premise is that agent requests are a hard case for the first finding: an honest request can also sound urgent or cite authority. That leaves two accounts:

- **Cue account (κ).** Training makes reviewers reject tactic-laden requests more often, whether safe or unsafe. Their criterion moves, but their discrimination doesn't.
- **Vigilance account (γ).** Training raises discrimination. Unsafe requests move toward rejection, and safe ones don't.

The design crosses tactics with safety, so the two accounts make different predictions in cells that both exist. Under the cue account, trained reviewers move along the control group's ROC curve. Under the vigilance account, the curve itself rises (Figure 1).

![Figure 1](output/figures/fig1_roc_cue_vs_vigilance.png)

*Figure 1.* Population ROC curves for tactic-laden requests on the argument screen. Points mark the five rating thresholds. Left: the cue account (κ = 0.5, γ = 0). Right: the vigilance account (κ = 0, γ = 0.3). The curves are asymmetric because evidence from unsafe requests varies more than evidence from safe ones.

## Data-generating model

Reviewer *i* sees request *j* and forms latent evidence toward rejecting it:

$$x_{ij} = \mu_{ij} + \sigma_j \varepsilon_{ij}, \qquad \varepsilon_{ij} \sim N(0, 1), \qquad \sigma_j = 1.25 \text{ if unsafe, else } 1$$

$$\mu_{ij} = U_j\left[d_0 + u^{d}_i + v^{d}_j + \beta_E E_i - \pi P_j (1 - E_i) + \gamma T_i\right] - A_{ij}\left(\phi + \phi_U U_j + u^{t}_i + v^{t}_j\right)(1 - s R_i) + A_{ij} T_i (\kappa + v^{k}_j)(1 - s R_i) + u^{c}_i + v^{c}_j - \delta (t_{ij} - 1)$$

The indicator variables:

| Symbol | Equals 1 when |
|---|---|
| *U* | the request is unsafe |
| *P* | the justification palters (is true but misleading) |
| *A* | the justification is tactic-laden |
| *T* | the reviewer was trained |
| *E* | the reviewer sees the evidence screen |
| *R* | the reviewer sees the paraphrase screen |

*t* is the trial position. The terms *u* are participant effects and the terms *v* are request effects; all are independent normals. The rating is 1 plus the number of thresholds (−0.6, 0.2, 1.0, 1.6, 2.2) that *x* exceeds.

| Parameter | Value | Basis |
|---|---|---|
| *d*₀, discrimination | 1.5 | Assumption |
| σ, SD of unsafe evidence | 1.25 (zROC slope 0.8) | Assumption. Confidence ratings identify it; a single approve/reject threshold can't |
| γ, vigilance gain in *d* | 0.3 | Matches the prebunking gains in Simchon et al. (2026), which came from news headlines, a different task |
| κ, cue shift toward rejecting | 0.5 | Assumption |
| φ, φ_U, pull of tactics toward approval | 0.4, plus 0.2 on unsafe requests | Assumption |
| *s*, share of tactic cues a paraphrase removes | 0.5 | Assumption |
| β_E, evidence gain in *d* | 0.5 | Assumption |
| π, cost of paltering without evidence | 0.6 (half of unsafe requests palter) | Assumption |
| δ, drift toward approval per trial | 0.006 | Assumption (approval fatigue) |
| Participant SDs: criterion, sensitivity, susceptibility | 0.5, 0.3, 0 (0.3 in script 03) | Assumption |
| Request SDs: criterion, sensitivity | 0.3, 0.3 | Assumption |
| Request SDs: tactic pull, trained reaction | 0, 0 (0.3, 0.2 in script 06) | Assumption |

Each simulated study draws new participants and new requests. Every value is set in `dgm_defaults()` and can be changed with `utils::modifyList()`.

## Analysis

### Preregistered primary model

`scripts/05_brms_primary_model.R` fits a hierarchical Bayesian cumulative-probit model with:

- unequal variances (`disc ~ 0 + unsafe`; brms's `disc` is the inverse SD of the latent variable; Bürkner & Vuorre, 2019);
- participant intercepts and sensitivities;
- request intercepts, tactic effects, and trained reviewers' reaction to each request's tactics.

It estimates participant-level effects, so the smallest effect of interest (SESOI) is a change of 0.20 in *d*. A full-size fit is expected to take hours; this hasn't been timed here.

### Fast estimator for the power analysis

`R/analyze.R` fits the same unequal-variance ordinal probit model without random effects (`ordinal::clm` with a scale term). Four points about it:

- It fits frequency-weighted response patterns, which reproduces the row-level fit to within 7 × 10⁻¹⁵ and takes a fraction of a second.
- Standard errors come from a participant-clustered sandwich built from analytic score contributions, which the tests check against numerical derivatives.
- It targets population-averaged parameters, which are smaller than participant-level ones.
- To compare like with like, the code computes its exact pseudo-true values: it integrates the random effects analytically and fits the estimator to the expected rating frequencies.

On the fast estimator's scale, a participant-level change of 0.20 in *d* is 0.173. On the scale of the binary equal-variance model used for comparison, it is 0.137. Both conversions use the generating variance components; a real study would take them from the pilot.

### Preregistered contrasts

Here *d* is the latent distance between unsafe and safe requests, in safe-request SD units.

| Test | Contrast | Prediction |
|---|---|---|
| T1, threat check | Tactic minus plain, unsafe requests, control group, argument screen | Lower rejection |
| T2, cue effect | On safe requests, the training effect for tactic-laden minus plain justifications | Higher rejection |
| T3, falsifier | Change in *d* from training, within plain and within tactic-laden requests | Greater |
| T4, T5 | Change in *d*, evidence minus argument and evidence minus paraphrase | Greater |

The six directional tests are Holm-corrected. T3 also gets two one-sided equivalence tests against ±SESOI (Lakens et al., 2018), and the falsifier counts as unmet only if both T3 contrasts are equivalent.

The verdicts are:

- **Cue:** T2 significant, both T3 contrasts equivalent, and no T3 gain.
- **Vigilance:** a T3 gain without T2.
- **Both:** a T3 gain together with T2.
- **Neither:** no T2, no T3 gain, and both T3 contrasts equivalent.
- **Inconclusive:** anything else.

The verdict doesn't depend on T1. A failed threat check would be reported as a limit of prompted attacks.

## What the simulation shows

Results come from `SIM_R=300 Rscript run_all.R`: 300 simulated studies per true state and sample size, with quick settings elsewhere. Every number is in `output/tables/`. With 300 studies, a rate near .90 has a Monte Carlo standard error of about 2 percentage points.

### 1. Without confidence ratings, a criterion shift looks like lost discrimination

Under the cue account, training moves only the criterion, so the true change in *d* is zero.

The unequal-variance model fitted to the confidence ratings recovers that zero; its pseudo-true value is below 0.001. It concluded that training lowered *d* for tactic-laden requests (one-sided *p* < .05, uncorrected) in 3% to 6% of studies, close to the nominal 5%.

An equal-variance probit fitted to the approve/reject decisions is signal detection theory written as a probit model (DeCarlo, 1998). With a single response threshold it can't separate the variance of unsafe-request evidence from its mean, so it reads part of the criterion shift as a loss of discrimination. Its pseudo-true change is −0.075. At 100 to 150 participants per cell:

- It reached the same wrong conclusion in 43% to 52% of studies (Figure 2).
- Its equivalence test, judged against its own SESOI, passed in 31% to 42% of studies, against 95% to 99% for the unequal-variance model.
- As a result, it returned the correct cue verdict in 21% to 37% of studies, against 91% to 94%.

In the one simulated study in script 01, a binary probit GLMM with crossed participant and request intercepts estimated the change at −0.16 (*SE* = 0.05), where the unequal-variance estimate on the same data was −0.05. Simchon et al. (2026) model unequal variances too, so this concerns analyses of binary decisions, not their meta-analysis.

![Figure 2](output/figures/fig2_equal_vs_unequal_variance.png)

*Figure 2.* Estimated change in *d* from training for tactic-laden requests when training only moved the criterion (cue account, 150 participants per cell, 300 simulated studies). The shaded band is ±SESOI on each model's own scale. The true change is zero.

### 2. Coverage of the estimator's intervals

Coverage was computed over 1,500 studies for each true state and contrast, pooling the five sample sizes. The Monte Carlo band around 95% is 93.9% to 96.1%.

- **Participant-clustered intervals** covered the pseudo-true value 93.1% to 96.5% of the time, and the SD of the standardized errors ranged from 0.94 to 1.07. The four lowest rates (93.1% to 94.0%) were all for T4 and T5 under the vigilance and both states, where the SD of the standardized errors was 1.04 to 1.07, so the standard errors for the evidence contrasts run slightly small there.
- **Model-based errors** covered 91.4% to 98.1%.
- **Errors clustered on requests as well** (Cameron et al., 2011) covered:
  - 87.1% to 92.6% for T1 to T3;
  - 99.6% to 99.9% for T4 and T5;
  - and in 26 to 36 of the 1,500 studies, T1 had a negative variance.

These results hold only while requests don't differ in how tactics and training affect them, which is what the defaults assume. Section 4 drops that assumption.

### 3. Sample size

Detecting a vigilance gain the size Simchon et al. report (γ = 0.3) took few participants: the T3 tests rejected in at least 96% of studies, even at 50 per cell. The harder claim is that training did *not* change discrimination, because both equivalence tests have to pass. Each passed in 58% to 60% of cue-account studies at 50 per cell, and in 95% to 98% at 100 to 125.

Across the true states:

- At 100 per cell (600 participants, 28,800 reviews), the rule reached the correct verdict in at least 86% of studies in every true state.
- At 125 per cell (750 participants), it did so in at least 94%.
- Under the vigilance account, the correct-verdict rate levels off near 95%. Once the other tests reject, Holm's procedure tests T2 at α = .05, so about 5% of those studies falsely detect a cue effect and get the verdict "both".

![Figure 3](output/figures/fig3_adjudication.png)

*Figure 3.* (A) Share of simulated studies in which the preregistered rule reaches the correct verdict, by sample size; series are offset slightly so overlapping points stay visible, and the grey line marks .80. (B) Verdicts by true state at 100 participants per cell.

### 4. Sensitivity

Script 06 changes one assumption at a time at 125 participants per cell, with 150 simulated studies per state. Monte Carlo standard errors are 2 to 3.5 points.

| Scenario | Correct verdict: cue / vigilance / neither | Lowest participant-clustered coverage |
|---|---|---|
| Defaults | .89 / .96 / .96 | 94% |
| Requests differ in tactic pull (SD 0.3) and in trained reviewers' reaction (SD 0.2) | .84 / .89 / .83 | 79% |
| The same, at 200 per cell | .85 / .78 / .75 | 72% |
| Smaller cue effect, κ = 0.25 | .92 / – / – | 93% |
| Vigilance gain at the SESOI, γ = 0.20 | – / .96 / – | 93% |
| Sensitivity SDs of 0.5 for participants and requests | .83 / .95 / .90 | 94% |
| σ = 1.5 | .94 / .95 / .89 | 94% |
| σ = 1.0 (equal variances) | .92 / .95 / .98 | 93% |

Dashes mark states the scenario doesn't change. The full table is `output/tables/sensitivity_summary.csv`.

Most single changes moved the correct-verdict rates by less than about 10 points. Larger sensitivity variances and σ = 1.5 hurt the cue and neither verdicts most.

The exception is variation between requests in how tactics and training affect them:

- **Participant-clustered intervals** undercovered: 79% for the worst contrast at 125 per cell.
- **Adding participants made coverage worse, not better** (72% at 200 per cell). Request-level variance doesn't shrink as participants are added, while the standard errors do.
- **Clustering on requests as well** kept coverage at 90% or higher. But the wider intervals cost equivalence power, and the cue and neither verdicts fell to .68 to .79.

So the remedy lies in the design and the primary model, not in the sample of people:

- more requests per participant, or a larger bank sampled per participant;
- random slopes by request for tactics and for training's effect on them, which the brms formula in script 05 includes.

Under the default assumptions, 125 per cell is enough. Whether it's enough in practice depends on how much requests vary, and the pilot has to estimate that before the sample size is preregistered.

### 5. Simulated reviewers

On its own, a simulator that leans toward rejecting and barely responds to tactics is badly biased. Its estimate of how often people approve tactic-laden unsafe requests was 21 to 22 percentage points too low, and its intervals never covered the true value.

PPI and PPI++ (Angelopoulos et al., 2023; Angelopoulos, Duchi, & Zrnic, 2023) combine 2,000 simulated reviewers with 200 people. They removed that bias and kept coverage between 93% and 96% at every level of profile accuracy. What they save depends on how well the simulator predicts each person:

- **Approval rate.** Effective sample size rose 10% to 11% when the profile predicted each person at *r* = .50, where the person-level correlation with the human outcome was ρ = .32. It rose 59% when the profile was perfect (ρ = .64).
- **Threat effect.** Each person's threat effect is a difference between two six-request averages, which is mostly noise. Even a perfect profile reached only ρ = .25, and the gain was 3% to 7%.

For comparison, Krsteski et al. (2025) report gains of up to 14% for rectified survey simulations. Under these assumptions, simulated reviewers trim the human sample a little, mostly for level outcomes.

![Figure 4](output/figures/fig4_simulated_reviewers_ppi.png)

*Figure 4.* (A) Coverage of 95% intervals and (B) gain in effective sample size from adding 2,000 simulated reviewers to 200 people. Both are plotted against how accurately each person's profile predicts their criterion, sensitivity, and susceptibility. The simulator leans toward rejecting and responds to tactics at 20% of human strength. In panel A, series are offset slightly.

### 6. Tell audit

The audit asks whether a ridge classifier, cross-validated by request, can tell unsafe from safe requests by the justification text alone. In two example banks of 48 requests (12 unsafe):

- the leaky bank reached an AUC of .72 (permutation *p* = .034);
- the clean bank reached .48 (*p* = .547).

Over 20 audits per condition, the audit flagged the leak (unsafe justifications minimizing 70% of the time, against 30%) in 35% of 48-request banks, 60% of 96-request banks, and 75% of 192-request banks. It flagged clean banks 0% to 10% of the time.

With only 20 audits per point, these rates are rough (Figure 5 shows Wilson intervals). For this leak and this classifier, a bank of about 50 requests often misses a real tell. The bank should be several times larger than the set any one participant sees, and the audit should be reported with its power. Blind human raters remain a separate check.

![Figure 5](output/figures/fig5_tell_audit_power.png)

*Figure 5.* Share of tell audits (grouped cross-validated ridge logistic regression, permutation *p* < .05) that flag a tell, by bank size, with Wilson 95% intervals.

## Running it

```bash
# R 4.3 or later. Packages: ordinal, brms (with rstan), lme4, future.apply,
# ggplot2, patchwork, here, testthat.
Rscript tests/testthat.R                     # the test suite, a few seconds
Rscript run_all.R                            # tests plus scripts 01-06, quick mode; about 30 minutes on 2 cores
SIM_FULL=1 SIM_WORKERS=7 Rscript run_all.R   # 1,000 studies per condition; one to two hours on 8 cores (estimate)
RUN_BRMS_FULL=1 BRMS_TRUTH=cue Rscript scripts/05_brms_primary_model.R   # full-size primary model
```

Environment variables:

- `SIM_R` and `SIM_R6` set the number of studies per condition in scripts 02 and 06.
- `SIM_REUSE=1` summarizes saved simulations without rerunning them.
- `BRMS_TRUTH` sets the true state for script 05.

For a pinned environment, run `renv::init()` once and commit the lockfile.

| File | What it does |
|---|---|
| `R/design.R` | Parameters, the four true states, the latent mean and marginal SD, cells, and the contrast matrix |
| `R/simulate.R` | Simulates one study; computes the exact population for pseudo-true values |
| `R/analyze.R` | Pattern-level fit, analytic scores, cluster-robust covariance, tests, verdict, SESOI conversion |
| `R/ppi.R` | Prediction-powered inference, PPI++, and simulated reviewers |
| `R/audit.R` | Toy request bank, ridge classifier, grouped folds, permutation test |
| `R/plots.R` | APA-style figures |
| `scripts/01`–`06` | Recovery; power and adjudication; simulated reviewers; tell audit; primary model; sensitivity |
| `tests/testthat/` | Counterbalancing, marginal SD against simulation, scores against numerical derivatives (summed and per row), the sandwich's bread and its behavior under independence, pseudo-true values, decision rule, PPI, AUC, folds |

## What this does not show

- **Most parameter values are guesses.** γ is matched to prebunking results from a different task. Section 4 shows how the answer moves when single values change, but not when several change at once.
- **The power analysis describes the fast estimator, not the preregistered Bayesian model.** Script 05 gives that model a decision rule and compares its verdict with the fast estimator's on the same data, but only once per run. Its error rates need repeated full-size fits.
- **The SESOI conversion uses the true variance components.** In practice they would be estimated, which adds uncertainty to the equivalence bounds.
- **The true states sit at the extremes** (κ of 0 or 0.5, γ of 0 or 0.3), apart from two sensitivity scenarios.
- **The request-level problem isn't solved here.** Section 4 shows that it exists and where the fix lies, but not how many requests the fix needs. That takes simulations with a request bank larger than any one participant sees, which the code doesn't yet support.
- **The simulator's failure modes are settings, not measurements.** Their sizes would come from the registered simulated-reviewer predictions in the real study.
- **The request bank in the audit is a toy.**

## References

Angelopoulos, A. N., Bates, S., Fannjiang, C., Jordan, M. I., & Zrnic, T. (2023). Prediction-powered inference. *Science, 382*(6671), 669–674. https://doi.org/10.1126/science.adi6000

Angelopoulos, A. N., Duchi, J. C., & Zrnic, T. (2023). *PPI++: Efficient prediction-powered inference* [Preprint]. arXiv. https://arxiv.org/abs/2311.01453

Bürkner, P.-C., & Vuorre, M. (2019). Ordinal regression models in psychology: A tutorial. *Advances in Methods and Practices in Psychological Science, 2*(1), 77–101. https://doi.org/10.1177/2515245918823199

Cameron, A. C., Gelbach, J. B., & Miller, D. L. (2011). Robust inference with multiway clustering. *Journal of Business & Economic Statistics, 29*(2), 238–249. https://doi.org/10.1198/jbes.2010.07136

DeCarlo, L. T. (1998). Signal detection theory and generalized linear models. *Psychological Methods, 3*(2), 186–205. https://doi.org/10.1037/1082-989X.3.2.186

Hughes, J. (2026, March 25). *How we built Claude Code auto mode: A safer way to skip permissions*. Anthropic. https://www.anthropic.com/engineering/claude-code-auto-mode

Krsteski, S., Russo, G., Chang, S., West, R., & Gligorić, K. (2025). *Valid survey simulations with limited human data: The roles of prompting, fine-tuning, and rectification* [Preprint]. arXiv. https://arxiv.org/abs/2510.11408

Lakens, D., Scheel, A. M., & Isager, P. M. (2018). Equivalence testing for psychological research: A tutorial. *Advances in Methods and Practices in Psychological Science, 1*(2), 259–269. https://doi.org/10.1177/2515245918770963

Levy, J., Yang, M., & Pelrine, K. (2026). *AI persuasion as a threat to human control* [Preprint]. arXiv. https://arxiv.org/abs/2609.14796

Modirrousta-Galian, A., & Higham, P. A. (2023). Gamified inoculation interventions do not improve discrimination between true and fake news: Reanalyzing existing research with receiver operating characteristic analysis. *Journal of Experimental Psychology: General, 152*(9), 2411–2437. https://doi.org/10.1037/xge0001395

Simchon, A., Zipori, T., Teitelbaum, L., Lewandowsky, S., & van der Linden, S. (2026). A signal detection theory meta-analysis of psychological inoculation against misinformation. *Current Opinion in Psychology, 67*, Article 102194. https://doi.org/10.1016/j.copsyc.2025.102194

## Software

The code was developed with R 4.3.3, using:

- ordinal 2023.12.4
- brms 2.20.4
- rstan 2.32.5
- lme4 1.1.35.1
- future.apply 1.11.1
- ggplot2 3.4.4
- patchwork 1.2.0
- testthat 3.2.1
- here 1.0.1

Each script writes its `sessionInfo()` to `output/`.
