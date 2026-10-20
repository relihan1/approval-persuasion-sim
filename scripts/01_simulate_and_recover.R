# 01_simulate_and_recover.R
# One simulated study, end to end:
#   1. the exact values the estimator targets under each account (pseudo-true values)
#   2. the smallest effect of interest on the estimator's scale
#   3. a simulated study under the cue account, fitted and compared with the truth
#   4. why the confidence ratings matter: an equal-variance analysis of the same
#      data finds a change in discrimination that isn't there
#   5. the applied outcome: unsafe approvals at a fixed rate of safe approvals
#   6. Figure 1, ROC curves under the two accounts

source(here::here("R", "load.R"))
set.seed(20261019)
pars <- dgm_defaults()

# ---- 1. Pseudo-true values ---------------------------------------------------
pt <- lapply(names(hypotheses()), function(h) pseudo_true(with_hypothesis(pars, h)))
names(pt) <- names(hypotheses())
pt_table <- do.call(rbind, lapply(names(pt), function(h) {
  data.frame(truth = h, test = names(pt[[h]]$uv),
             unequal_variance = round(pt[[h]]$uv, 4),
             equal_variance_binary = round(pt[[h]]$ev, 4), row.names = NULL)
}))
write.csv(pt_table, out_path("tables", "pseudo_true_values.csv"), row.names = FALSE)
print(pt_table)

# ---- 2. Smallest effect of interest -----------------------------------------
# Preregistered on the participant level as a change of 0.20 in d.
sesoi    <- sesoi_on_estimator_scale(pars, delta = 0.20)
sesoi_ev <- sesoi_on_estimator_scale(pars, delta = 0.20, model = "ev")
cat(sprintf("\nSESOI: 0.20 per participant = %.3f on the unequal-variance estimator's scale, %.3f on the binary model's\n",
            sesoi, sesoi_ev))

# ---- 3. One study under the cue account -------------------------------------
p_cue <- with_hypothesis(pars, "cue")
dat <- simulate_study(p_cue, seed = 1)
cat(sprintf("\n%d participants x %d requests = %d reviews\n",
            length(unique(dat$pid)), length(unique(dat$item)), nrow(dat)))

# The pattern-level fit is the row-level fit, computed faster.
a <- aggregate_reviews(dat)
fit_agg <- fit_uv(a$agg)
dat_rows <- transform(dat, y = factor(rating, levels = 1:6, ordered = TRUE),
                      unsafe_s = unsafe, trial_c = block_c)
fit_rows <- clm(location_formula, scale = ~ unsafe_s, data = dat_rows, link = "probit")
max_diff <- max(abs(coef(fit_rows) - coef(fit_agg)))
cat(sprintf("Largest coefficient difference, row-level vs pattern-level fit: %.1e\n", max_diff))
write.csv(data.frame(check = c("sesoi_unequal_variance", "sesoi_equal_variance_binary",
                               "max_coef_diff_row_vs_pattern_fit"),
                     value = c(sesoi, sesoi_ev, max_diff)),
          out_path("tables", "estimator_checks.csv"), row.names = FALSE)

res <- analyze_study(dat, sesoi, sesoi_ev)
recovery <- transform(res$tests, truth = pt$cue$uv[test])
recovery <- recovery[, c("test", "truth", "est", "se", "se_model", "se_two_way",
                         "p_holm", "sig", "p_tost", "equiv")]
write.csv(recovery, out_path("tables", "recovery_one_study_cue.csv"), row.names = FALSE)
print(recovery, digits = 3)
cat(sprintf("Verdict: %s | zROC slope: %.2f (pseudo-true %.2f; 0.80 per participant)\n",
            res$verdict, 1 / res$sigma_ratio, 1 / pt$cue$sigma_ratio))

# ---- 4. The equal-variance analysis of the same data ------------------------
# A binary probit GLMM with crossed participant and request intercepts
# (DeCarlo, 1998). With one approve/reject threshold, each design cell gives a
# single ROC point (ratings trace the ROC within a cell; DeCarlo, 1998, p. 189),
# so with a free mean in every cell the model must assume equal variances, and
# the criterion shift that training produces on tactic-laden requests shows up
# as a change in d.
glmm <- lme4::glmer(reject ~ unsafe * tactic * trained * screen + palter + palter:screen + trial_c +
                      (1 | pid) + (1 | item),
                    data = dat, family = binomial(link = "probit"), nAGQ = 0)
L <- contrast_L()
b <- lme4::fixef(glmm)[colnames(L)]
V <- as.matrix(vcov(glmm))[colnames(L), colnames(L)]
ev_glmm <- data.frame(test = rownames(L), est = drop(L %*% b), se = sqrt(diag(L %*% V %*% t(L))))
ev_glmm$uv_est <- res$tests$est  # population-averaged scale; the GLMM is participant-level
# Under the cue account the true change in d (T3) is zero on any scale.
print(ev_glmm, digits = 3)
write.csv(ev_glmm, out_path("tables", "equal_variance_glmm_one_study_cue.csv"), row.names = FALSE)

# ---- 5. Applied outcome ------------------------------------------------------
# Unsafe approvals when the threshold approves 90% of safe requests.
g <- cell_eta(fit_agg)
d_cells <- merge(g[g$unsafe == 1, c("tactic", "trained", "screen", "eta")],
                 g[g$unsafe == 0, c("tactic", "trained", "screen", "eta")],
                 by = c("tactic", "trained", "screen"), suffixes = c("_unsafe", "_safe"))
d_cells$d <- d_cells$eta_unsafe - d_cells$eta_safe
d_cells$unsafe_approved_at_90 <- unsafe_approval_at(d_cells$d, exp(fit_agg$zeta), 0.90)
d_cells <- d_cells[order(d_cells$screen, d_cells$tactic, d_cells$trained), ]
print(d_cells[, c("screen", "tactic", "trained", "d", "unsafe_approved_at_90")], digits = 3)
write.csv(d_cells, out_path("tables", "applied_outcome_one_study_cue.csv"), row.names = FALSE)

# ---- 6. Figure 1 -----------------------------------------------------------------
fig1 <- plot_roc(list(cue = pt$cue$uv_fit, vigilance = pt$vigilance$uv_fit))
save_fig(fig1, "fig1_roc_cue_vs_vigilance", width = 7, height = 4.2)

writeLines(capture.output(sessionInfo()), out_path("sessionInfo_01.txt"))
