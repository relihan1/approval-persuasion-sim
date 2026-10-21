# 05_brms_primary_model.R
# The preregistered primary model: a hierarchical Bayesian unequal-variance
# ordinal signal detection model with crossed participant and request effects.
#
#   rating ~ cumulative probit
#   location: unsafe x tactic x trained x screen + trial position
#             + paltering and the evidence screen's removal of its cost,
#               palter + palter:screen (the fast estimator's terms)
#             + participant criterion, sensitivity, and susceptibility to
#               tactics, (1 + unsafe + tactic | pid)
#             + request location, tactic effect, and trained reaction to
#               tactics, (1 + tactic + tactic:trained | item); script 06 shows
#               why the request slopes matter
#             + the extra spread of unsafe requests, (0 + unsafe | item),
#               uncorrelated with the request location as in the simulation
#   disc:     ~ 0 + unsafe. brms's disc is the inverse SD of the latent
#             variable, fixed at 1 for safe requests, so exp(-b_disc_unsafe)
#             is the SD of unsafe-request evidence (Bürkner & Vuorre, 2019).
# The terms match latent_mean() in R/design.R one for one. Under the default
# parameters three of the random-slope SDs are 0 (tau_t, omega_t, omega_k).
#
# Its estimates are on the participant-level scale, so the smallest effect of
# interest is the preregistered 0.20 change in d, with no conversion.
#
# A full-size fit takes hours (README), which is why the power analysis uses
# the fast population-averaged estimator. By default this
# script fits a tiny study (72 participants, 24 requests) only to check that the
# model compiles, samples, and returns the contrasts; it is far too small to
# recover them.
# Set RUN_BRMS_FULL=1 to fit a full-size study, and BRMS_TRUTH to choose the
# true state (default "cue").

source(here::here("R", "load.R"))
suppressPackageStartupMessages(library(brms))

primary_formula <- bf(
  rating ~ unsafe * tactic * trained * screen + palter + palter:screen + trial_c +
    (1 + unsafe + tactic | pid) + (1 + tactic + tactic:trained | item) + (0 + unsafe | item)
) + lf(disc ~ 0 + unsafe, cmc = FALSE)

primary_priors <- c(
  prior(normal(0, 1), class = "b"),
  prior(normal(0, 1.5), class = "Intercept"),
  prior(normal(0, 0.5), class = "sd"),
  prior(lkj(2), class = "cor"),
  prior(normal(0, 0.5), class = "b", dpar = "disc")
)

# Posterior summaries of the preregistered contrasts and the decision rule.
# A directional effect counts when its 95% credible interval excludes zero in
# the predicted direction. "No change in d" (T3) counts when the 90% interval
# lies inside +/- 0.20, the Bayesian counterpart of two one-sided tests at
# alpha = .05. The same verdict() then applies. This rule has no multiplicity
# correction; its error rates have to be checked by repeated full-size fits.
posterior_contrasts <- function(fit, L = contrast_L(), sesoi = 0.20) {
  # The formula has the palter terms, so every contrast weight is used and the
  # contrasts match the fast estimator's. drop_absent_terms() drops nothing
  # here; it would stop if a term other than a palter term were missing.
  L <- drop_absent_terms(L, sub("^b_", "", brms::variables(fit)))
  cd <- contrast_draws(fit, L)  # keeps the chains, for R-hat and ESS of each contrast
  est <- posterior::as_draws_matrix(cd)
  q <- function(p) apply(est, 2, stats::quantile, p)
  conv <- posterior::summarise_draws(cd, "rhat", "ess_bulk", "ess_tail")
  out <- data.frame(test = colnames(est), median = apply(est, 2, stats::median),
                    lo95 = q(0.025), hi95 = q(0.975), lo90 = q(0.05), hi90 = q(0.95),
                    p_positive = colMeans(est > 0), p_inside_sesoi = colMeans(abs(est) < sesoi),
                    rhat = conv$rhat, ess_bulk = conv$ess_bulk, ess_tail = conv$ess_tail,
                    row.names = NULL)
  dir <- attr(L, "direction")[out$test]
  out$sig <- ifelse(dir == "greater", out$lo95 > 0, out$hi95 < 0)
  out$equiv <- out$lo90 > -sesoi & out$hi90 < sesoi
  out
}

full  <- identical(Sys.getenv("RUN_BRMS_FULL"), "1")
truth <- Sys.getenv("BRMS_TRUTH", "cue")  # cue, vigilance, both, or null
pars  <- with_hypothesis(dgm_defaults(), truth)
if (!full) pars <- utils::modifyList(pars, list(n_per_cell = 12, n_safe = 18, n_unsafe = 6))

dat <- simulate_study(pars, seed = 5)
cat(sprintf("Fitting %d reviews from %d participants\n", nrow(dat), length(unique(dat$pid))))
tag <- sprintf("%s_%s", if (full) "full" else "tiny", truth)

# Starting values: rstan's default draws every unconstrained parameter from
# (-2, 2). brms writes the cumulative-probit likelihood as log(Phi - Phi) for
# each review, which underflows to log(0) when a linear predictor lands far
# from the thresholds, so only about 1% of those draws have a finite log
# density at full size, and a chain that finds none in 100 tries is dropped.
# From (-0.5, 0.5), about 70% do. Starting values don't change the posterior.
# Iterations: at a third of full size, a warmup iteration cost about 135
# gradient evaluations and a sampling iteration about 95, so 2,000 sampling
# iterations per chain instead of 1,000 add about 40% to the run time. With
# 1,000, the SD of participants' susceptibility to tactics (0 in the
# simulation) reached bulk ESS of only 237 from four chains, and the earlier
# full-size fits of the smaller formula reached only 139 to 276 for the
# location coefficients.
# BRMS_THREADS > 1 splits each chain's likelihood across that many threads
# (within-chain parallelism, rstan >= 2.26); it changes run time, not the
# model, but the draws differ from an unthreaded run with the same seed.
threads <- as.integer(Sys.getenv("BRMS_THREADS", "1"))
fit <- brm(primary_formula, data = dat, family = cumulative(link = "probit"),
           prior = primary_priors, backend = "rstan",
           chains = if (full) 4 else 2, iter = if (full) 3000 else 1000,
           warmup = if (full) 1000 else 500, init_r = 0.5,
           cores = min(4, future::availableCores()), seed = 20261019,
           threads = if (threads > 1) threading(threads) else NULL,
           control = list(adapt_delta = 0.9), refresh = 0)
# Saved so that its diagnostics can be revisited without refitting
# (output/fits/ is not under version control).
saveRDS(fit, out_path("fits", sprintf("brms_%s.rds", tag)))

pc <- posterior_contrasts(fit)
print(pc, digits = 3)

# Convergence diagnostics for every parameter class and the six contrasts
# (R/brms_check.R), and the sampler's own warnings.
L_fit <- drop_absent_terms(contrast_L(), sub("^b_", "", brms::variables(fit)))
by_class <- convergence_by_class(posterior::as_draws_array(fit), contrast_draws(fit, L_fit))
sampler <- sampler_summary(brms::nuts_params(fit), max_treedepth = 10)
print(by_class, digits = 3)
print(sampler, digits = 3)
checked <- by_class[!by_class$class %in% c("r_pid", "r_item"), ]
ess_checked <- checked[checked$class %in% ess_classes, ]
diagnostics <- data.frame(
  n_chains = posterior::nchains(fit),
  max_rhat = max(checked$max_rhat), worst_rhat = checked$worst_rhat[which.max(checked$max_rhat)],
  min_ess_bulk = min(ess_checked$min_ess_bulk), min_ess_tail = min(ess_checked$min_ess_tail),
  min_ess_bulk_other = min(checked$min_ess_bulk[!checked$class %in% ess_classes]),
  divergent = sum(sampler$divergent), treedepth_hits = sum(sampler$treedepth_hits),
  min_ebfmi = min(sampler$ebfmi),
  converged = convergence_ok(by_class, sampler, n_chains = posterior::nchains(fit)))
print(diagnostics)

# The fast estimator's verdict on the same data, for comparison.
fast <- analyze_study(dat, sesoi_on_estimator_scale(pars), with_ev = FALSE)
comparison <- data.frame(truth = truth, full_size = full, verdict_primary = verdict(pc),
                         verdict_fast = fast$verdict, diagnostics)
print(comparison)
sd_unsafe <- exp(-as.matrix(fit, variable = "b_disc_unsafe"))
cat(sprintf("SD of unsafe-request evidence: median %.2f (generating value %.2f)\n",
            median(sd_unsafe), pars$sigma_unsafe))
write.csv(pc, out_path("tables", sprintf("brms_contrasts_%s.csv", tag)), row.names = FALSE)
write.csv(comparison, out_path("tables", sprintf("brms_vs_fast_%s.csv", tag)), row.names = FALSE)
write.csv(by_class, out_path("tables", sprintf("brms_convergence_%s.csv", tag)), row.names = FALSE)
writeLines(capture.output(sessionInfo()), out_path("sessionInfo_05.txt"))
