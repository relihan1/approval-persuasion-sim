# 05_brms_primary_model.R
# The preregistered primary model: a hierarchical Bayesian unequal-variance
# ordinal signal detection model with crossed participant and request effects.
#
#   rating ~ cumulative probit
#   location: unsafe x tactic x trained x screen + trial position
#             + participant criterion and sensitivity, (1 + unsafe | pid)
#             + request location, tactic effect, and trained reaction to
#               tactics, (1 + tactic + tactic:trained | item); script 06 shows
#               why the request slopes matter
#   disc:     ~ 0 + unsafe. brms's disc is the inverse SD of the latent
#             variable, fixed at 1 for safe requests, so exp(-b_disc_unsafe)
#             is the SD of unsafe-request evidence (Bürkner & Vuorre, 2019).
#
# Its estimates are on the participant-level scale, so the smallest effect of
# interest is the preregistered 0.20 change in d, with no conversion.
#
# A full-size fit is expected to take hours (not timed here), which is why the
# power analysis uses the fast population-averaged estimator. By default this
# script fits a tiny study (72 participants, 24 requests) only to check that the
# model compiles, samples, and returns the contrasts; it is far too small to
# recover them.
# Set RUN_BRMS_FULL=1 to fit a full-size study, and BRMS_TRUTH to choose the
# true state (default "cue").

source(here::here("R", "load.R"))
suppressPackageStartupMessages(library(brms))

primary_formula <- bf(
  rating ~ unsafe * tactic * trained * screen + trial_c + (1 + unsafe | pid) +
    (1 + tactic + tactic:trained | item)
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
  # This formula has no palter terms (the fast estimator's does); their contrast
  # weights drop out, and the request intercepts absorb paltering's average cost.
  L <- drop_absent_terms(L, sub("^b_", "", brms::variables(fit)))
  draws <- as.matrix(fit, variable = paste0("b_", colnames(L)))[, paste0("b_", colnames(L))]
  est <- draws %*% t(L)
  q <- function(p) apply(est, 2, stats::quantile, p)
  out <- data.frame(test = colnames(est), median = apply(est, 2, stats::median),
                    lo95 = q(0.025), hi95 = q(0.975), lo90 = q(0.05), hi90 = q(0.95),
                    p_positive = colMeans(est > 0), p_inside_sesoi = colMeans(abs(est) < sesoi),
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

fit <- brm(primary_formula, data = dat, family = cumulative(link = "probit"),
           prior = primary_priors, backend = "rstan",
           chains = if (full) 4 else 2, iter = if (full) 2000 else 1000,
           cores = min(4, future::availableCores()), seed = 20261019,
           control = list(adapt_delta = 0.9), refresh = 0)

pc <- posterior_contrasts(fit)
print(pc, digits = 3)

# Convergence diagnostics for the location and disc coefficients.
dg <- posterior::summarise_draws(posterior::subset_draws(posterior::as_draws_df(fit), variable = "^b_",
                                                         regex = TRUE), "rhat", "ess_bulk", "ess_tail")
diagnostics <- data.frame(max_rhat = max(dg$rhat), min_ess_bulk = min(dg$ess_bulk),
                          min_ess_tail = min(dg$ess_tail))
print(diagnostics)

# The fast estimator's verdict on the same data, for comparison.
fast <- analyze_study(dat, sesoi_on_estimator_scale(pars), with_ev = FALSE)
comparison <- data.frame(truth = truth, full_size = full, verdict_primary = verdict(pc),
                         verdict_fast = fast$verdict, diagnostics)
print(comparison)
sd_unsafe <- exp(-as.matrix(fit, variable = "b_disc_unsafe"))
cat(sprintf("SD of unsafe-request evidence: median %.2f (generating value %.2f)\n",
            median(sd_unsafe), pars$sigma_unsafe))
tag <- sprintf("%s_%s", if (full) "full" else "tiny", truth)
write.csv(pc, out_path("tables", sprintf("brms_contrasts_%s.csv", tag)), row.names = FALSE)
write.csv(comparison, out_path("tables", sprintf("brms_vs_fast_%s.csv", tag)), row.names = FALSE)
writeLines(capture.output(sessionInfo()), out_path("sessionInfo_05.txt"))
