# 05b_brms_rule_comparison.R
# Script 05's Bayesian rule calls a directional effect when the 95% credible
# interval excludes zero, about a one-sided test at .025 with no multiplicity
# correction. The fast estimator uses one-sided tests at .05 with Holm's
# correction. A disagreement between the two verdicts could therefore come from
# the rule rather than the model. This script re-reads the saved posterior
# summaries of the full-size fits and adds a verdict under a rule matched to the
# fast estimator's: the posterior probability of the wrong direction stands in
# for the one-sided p-value, and Holm's step-down applies across the six tests
# at .05. Equivalence is unchanged (90% interval inside +/- 0.20). The matched
# rule is a diagnostic for the comparison, not a replacement for the
# preregistered Bayesian rule.
#
# It also re-derives the fast estimator's verdicts on the same simulated study
# (script 05 uses seed 5) with the current estimator, under participant-
# clustered errors and under two-way errors with the request jackknife, so the
# comparison doesn't depend on which version of the estimator script 05 ran.

source(here::here("R", "load.R"))

files <- Sys.glob(here::here("output", "tables", "brms_contrasts_full_*.csv"))
if (!length(files)) stop("No full-size fits found; run RUN_BRMS_FULL=1 scripts/05_brms_primary_model.R first.")

dirs <- attr(contrast_L(), "direction")
rows <- lapply(files, function(f) {
  truth <- sub("^brms_contrasts_full_(.*)\\.csv$", "\\1", basename(f))
  pc <- read.csv(f)
  p_wrong <- ifelse(dirs[pc$test] == "greater", 1 - pc$p_positive, pc$p_positive)
  pc$sig <- stats::p.adjust(p_wrong, method = "holm", n = nrow(pc)) < 0.05
  cmp <- read.csv(here::here("output", "tables", sprintf("brms_vs_fast_full_%s.csv", truth)))
  pars <- with_hypothesis(dgm_defaults(), truth)
  fast <- analyze_study(simulate_study(pars, seed = 5), sesoi_on_estimator_scale(pars), with_ev = FALSE)
  data.frame(truth = truth, verdict_primary = cmp$verdict_primary,
             verdict_primary_holm = verdict(pc), verdict_fast_pid = fast$verdict_pid,
             verdict_fast_two_way_jk = fast$verdict_two_way_jk,
             correct = correct_verdict[[truth]], max_rhat = cmp$max_rhat,
             min_ess_bulk = cmp$min_ess_bulk, min_ess_tail = cmp$min_ess_tail,
             # NA for fits made before script 05 checked every parameter class,
             # whose R-hat and ESS cover the location coefficients only.
             converged = if (is.null(cmp$converged)) NA else cmp$converged)
})
out <- do.call(rbind, rows)
print(out, row.names = FALSE)
write.csv(out, out_path("tables", "brms_rule_comparison.csv"), row.names = FALSE)
