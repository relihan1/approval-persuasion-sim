# 02_power_and_adjudication.R
# How many participants does the study need to tell the accounts apart?
#
# For each true state (cue, vigilance, both, neither) and sample size, simulate
# many studies, run the preregistered analysis on each, and record:
#   - power for each directional test (Holm-corrected) and for equivalence
#   - which verdict the decision rule reaches (the adjudication matrix)
#   - whether 95% intervals cover the pseudo-true values (estimator validity),
#     with participant-clustered, two-way clustered, and model-based errors
#   - what the equal-variance (binary) analysis would have concluded instead

source(here::here("R", "load.R"))
library(future.apply)

pars  <- dgm_defaults()
sesoi    <- sesoi_on_estimator_scale(pars, delta = 0.20)
sesoi_ev <- sesoi_on_estimator_scale(pars, delta = 0.20, model = "ev")  # the binary model's own bound
truth <- lapply(names(hypotheses()), function(h) pseudo_true(with_hypothesis(pars, h))$uv)
names(truth) <- names(hypotheses())

n_grid <- c(50, 75, 100, 125, 150)
R <- as.integer(Sys.getenv("SIM_R", if (QUICK) 100 else 1000))  # studies per condition
workers <- as.integer(Sys.getenv("SIM_WORKERS", future::availableCores()))
plan(multisession, workers = workers)

conditions <- expand.grid(truth = names(hypotheses()), n_per_cell = n_grid, rep = seq_len(R),
                          stringsAsFactors = FALSE)
conditions$seed <- seq_len(nrow(conditions)) + 1e6  # every study has its own seed

one_rep <- function(i) {
  cnd <- conditions[i, ]
  p <- with_hypothesis(pars, cnd$truth)
  p$n_per_cell <- cnd$n_per_cell
  res <- analyze_study(simulate_study(p, seed = cnd$seed), sesoi, sesoi_ev)
  uv <- res$tests; ev <- res$tests_ev
  data.frame(truth = cnd$truth, n_per_cell = cnd$n_per_cell, rep = cnd$rep, test = uv$test,
             est = uv$est, se = uv$se, se_two_way = uv$se_two_way, se_model = uv$se_model,
             sig = uv$sig, equiv = uv$equiv, verdict = res$verdict,
             ev_est = ev$est, ev_se = ev$se, ev_sig = ev$sig, ev_equiv = ev$equiv,
             verdict_ev = res$verdict_ev)
}

# SIM_REUSE=1 summarizes the saved simulations instead of rerunning them.
rds <- out_path("tables", "power_simulations.rds")
if (identical(Sys.getenv("SIM_REUSE"), "1") && file.exists(rds)) {
  sims <- readRDS(rds)
  R <- max(sims$rep)
} else {
  t0 <- Sys.time()
  sims <- do.call(rbind, future_lapply(seq_len(nrow(conditions)), one_rep,
                                       future.seed = TRUE, future.packages = "ordinal"))
  cat(sprintf("%d simulated studies in %.1f minutes\n", nrow(conditions),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  saveRDS(sims, rds)
}

# ---- Power ---------------------------------------------------------------------
power <- aggregate(cbind(sig, equiv) ~ truth + n_per_cell + test, data = sims, FUN = mean)
power$mcse_sig <- sqrt(power$sig * (1 - power$sig) / R)  # Monte Carlo standard errors
power$mcse_equiv <- sqrt(power$equiv * (1 - power$equiv) / R)
power <- power[order(power$truth, power$test, power$n_per_cell), ]
write.csv(power, out_path("tables", "power_by_test.csv"), row.names = FALSE)

# ---- Adjudication ---------------------------------------------------------------
studies <- unique(sims[, c("truth", "n_per_cell", "rep", "verdict", "verdict_ev")])
verdicts <- c("cue", "vigilance", "both", "null", "inconclusive")
tab_verdicts <- function(col) {
  out <- do.call(rbind, lapply(split(studies, list(studies$truth, studies$n_per_cell)), function(s) {
    data.frame(truth = s$truth[1], n_per_cell = s$n_per_cell[1], verdict = verdicts,
               share = as.numeric(table(factor(s[[col]], levels = verdicts))) / nrow(s))
  }))
  rownames(out) <- NULL
  out
}
confusion_uv <- tab_verdicts("verdict")
confusion_ev <- tab_verdicts("verdict_ev")
write.csv(confusion_uv, out_path("tables", "verdicts_unequal_variance.csv"), row.names = FALSE)
write.csv(confusion_ev, out_path("tables", "verdicts_equal_variance.csv"), row.names = FALSE)

correct <- function(conf) {
  ok <- conf[conf$verdict == correct_verdict[conf$truth], c("truth", "n_per_cell", "share")]
  names(ok)[3] <- "p_correct"
  ok
}
rates <- merge(correct(confusion_uv), setNames(correct(confusion_ev), c("truth", "n_per_cell", "p_correct_ev")))
rates$mcse <- sqrt(rates$p_correct * (1 - rates$p_correct) / R)
print(rates[order(rates$truth, rates$n_per_cell), ], digits = 3)
write.csv(rates, out_path("tables", "correct_verdict_rates.csv"), row.names = FALSE)

# Smallest sample size at which every true state gets the right verdict at least 80% of the time.
ok_n <- sort(unique(rates$n_per_cell[ave(rates$p_correct, rates$n_per_cell, FUN = min) >= 0.80]))
n_show <- if (length(ok_n)) min(ok_n) else max(n_grid)
cat(sprintf("Smallest n per cell with >= .80 correct verdicts in every state: %s\n",
            if (length(ok_n)) n_show else "none in the grid"))

# ---- Estimator validity: coverage of 95% intervals --------------------------------
sims$truth_value <- mapply(function(h, t) truth[[h]][[t]], sims$truth, sims$test)
cover <- function(se) abs(sims$est - sims$truth_value) <= qnorm(0.975) * se
coverage <- aggregate(cbind(pid = cover(sims$se), two_way = cover(sims$se_two_way),
                            model = cover(sims$se_model)) ~ truth + test, data = sims,
                      FUN = mean, na.rm = TRUE, na.action = na.pass)
coverage$two_way_na <- aggregate(cbind(na = is.na(sims$se_two_way)) ~ truth + test, sims, sum)$na
coverage$studies <- aggregate(k ~ truth + test, transform(sims, k = 1), sum)$k
coverage$mcse_at_95 <- sqrt(0.95 * 0.05 / coverage$studies)  # Monte Carlo SE of a 95% coverage rate
# SD of standardized errors should be near 1 when the standard errors are right.
sims$z <- (sims$est - sims$truth_value) / sims$se
coverage <- merge(coverage, aggregate(cbind(sd_z = z) ~ truth + test, sims, sd))
# Two-way errors with t critical values on (requests - 1) df, the small-sample
# adjustment Cameron et al. (2011) use for few clusters.
k_req <- pars$n_safe + pars$n_unsafe
coverage <- merge(coverage, aggregate(cbind(two_way_t = abs(sims$est - sims$truth_value) <=
                                              qt(0.975, k_req - 1) * sims$se_two_way) ~ truth + test,
                                      data = sims, FUN = mean, na.rm = TRUE, na.action = na.pass))
# Calibration: mean standard error over the SD of the estimates, within each
# true state and sample size, then averaged over sample sizes (1 = right size).
cal <- do.call(rbind, lapply(split(sims, list(sims$truth, sims$n_per_cell, sims$test), drop = TRUE), function(d)
  data.frame(truth = d$truth[1], test = d$test[1], se_ratio_pid = mean(d$se) / sd(d$est),
             se_ratio_two_way = mean(d$se_two_way, na.rm = TRUE) / sd(d$est))))
coverage <- merge(coverage, aggregate(cbind(se_ratio_pid, se_ratio_two_way) ~ truth + test, cal, mean))
print(coverage, digits = 3)
write.csv(coverage, out_path("tables", "coverage_95.csv"), row.names = FALSE)

# ---- If the pilot misestimates the criterion variance --------------------------------
# The SESOI conversion divides by the marginal SD of safe-request evidence, which
# depends on tau_c and omega_c. Recompute every study's verdict with the bound a
# pilot would give if it misestimated tau_c (true value 0.5), or if the 0.20 were
# used without conversion.
bounds <- c(setNames(sapply(c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8), function(tc)
  sesoi_on_estimator_scale(utils::modifyList(pars, list(tau_c = tc)))), paste0("tau_c = ", c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8))),
  "no conversion" = 0.20)
by_study <- split(sims[, c("truth", "n_per_cell", "rep", "test", "est", "se", "sig")],
                  list(sims$truth, sims$n_per_cell, sims$rep), drop = TRUE)
misest <- do.call(rbind, lapply(names(bounds), function(b) {
  v <- vapply(by_study, function(d) {
    d$equiv <- run_tests(transform(d, direction = "greater"), bounds[[b]])$equiv  # TOST ignores direction
    verdict(d)
  }, character(1))
  key <- do.call(rbind, lapply(by_study, function(d) d[1, c("truth", "n_per_cell")]))
  out <- aggregate(cbind(p_correct = v == correct_verdict[key$truth]) ~ truth + n_per_cell, key, mean)
  cbind(assumed = b, sesoi = bounds[[b]], out)
}))
misest <- misest[order(misest$assumed, misest$truth, misest$n_per_cell), ]
print(misest[misest$n_per_cell %in% c(100, 125), ], digits = 3, row.names = FALSE)
write.csv(misest, out_path("tables", "sesoi_misestimation.csv"), row.names = FALSE)

# ---- What an equal-variance analysis concludes under the cue account ---------------
cue_tac <- sims[sims$truth == "cue" & sims$test == "T3_tactic", ]
ev_harm <- aggregate(cbind(
  uv_says_training_lowered_d = pnorm(est / se) < 0.05,
  ev_says_training_lowered_d = pnorm(ev_est / ev_se) < 0.05,
  uv_equivalence = equiv, ev_equivalence = ev_equiv) ~ n_per_cell, data = cue_tac, FUN = mean)
print(ev_harm, digits = 3)
write.csv(ev_harm, out_path("tables", "equal_variance_false_conclusions.csv"), row.names = FALSE)

# ---- Figures ------------------------------------------------------------------------
n_hist <- max(n_grid)
h <- cue_tac[cue_tac$n_per_cell == n_hist, ]
est_long <- rbind(data.frame(model = "Unequal-variance ordinal model (confidence ratings)", est = h$est),
                  data.frame(model = "Equal-variance binary model (approve or reject)", est = h$ev_est))
est_long$model <- factor(est_long$model, levels = unique(est_long$model))
bands <- data.frame(model = levels(est_long$model), sesoi = c(sesoi, sesoi_ev))
bands$model <- factor(bands$model, levels = levels(est_long$model))
save_fig(plot_ev_vs_uv(est_long, bands), "fig2_equal_vs_unequal_variance", width = 7, height = 4.5)

conf_show <- confusion_uv[confusion_uv$n_per_cell == n_show, ]
save_fig(plot_adjudication(rates, conf_show, n_show), "fig3_adjudication", width = 7, height = 7)

plan(sequential)
writeLines(capture.output(sessionInfo()), out_path("sessionInfo_02.txt"))
