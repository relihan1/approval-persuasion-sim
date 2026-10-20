# 06_sensitivity.R
# How much do the conclusions depend on assumptions the defaults don't vary?
# At 125 participants per cell (the size script 02 points to), each scenario
# changes one thing and reruns the cue, vigilance, and neither states:
#   - requests differ in how much their tactics pull and in how strongly
#     trained reviewers react to them (request-level random slopes), swept
#     over 100, 125, 150, and 200 participants per cell
#   - a smaller cue effect (kappa = 0.25)
#   - a vigilance gain equal to the smallest effect of interest (gamma = 0.2)
#   - larger sensitivity variance across participants and requests (0.5)
#   - unequal variance stronger (sigma = 1.5) or absent (sigma = 1.0)
#   - tactics pull unsafe requests no harder than safe ones (pull_extra = 0)
#   - larger criterion variance across participants (tau_c = 0.8) or requests
#     (omega_c = 0.6), the components the SESOI conversion depends on
#   - request slopes with twice as many requests per participant (72 safe, 24 unsafe)
# The verdict is recorded under participant-clustered and two-way clustered
# errors, and coverage is checked against each scenario's own pseudo-true values.

source(here::here("R", "load.R"))
library(future.apply)

base <- dgm_defaults()
base$n_per_cell <- 125
R <- as.integer(Sys.getenv("SIM_R6", if (QUICK) 150 else 1000))
plan(multisession, workers = as.integer(Sys.getenv("SIM_WORKERS", future::availableCores())))

slopes <- list(omega_t = 0.3, omega_k = 0.2)
scenarios <- list(
  default            = list(),
  request_slopes_100 = c(slopes, n_per_cell = 100),
  request_slopes     = slopes,  # 125 per cell
  request_slopes_150 = c(slopes, n_per_cell = 150),
  request_slopes_200 = c(slopes, n_per_cell = 200),
  request_slopes_96_requests = c(slopes, n_safe = 72, n_unsafe = 24),  # 125 per cell
  kappa_0.25         = list(kappa_override = 0.25),
  gamma_at_sesoi     = list(gamma_override = 0.2),
  sensitivity_sd_0.5 = list(tau_d = 0.5, omega_d = 0.5),
  sigma_1.5          = list(sigma_unsafe = 1.5),
  sigma_1.0          = list(sigma_unsafe = 1.0),
  pull_extra_0       = list(pull_extra = 0),
  tau_c_0.8          = list(tau_c = 0.8),
  omega_c_0.6        = list(omega_c = 0.6)
)
states <- c("cue", "vigilance", "null")

# A scenario's hypothesis values replace the defaults only where the state uses them.
scenario_pars <- function(sc, state) {
  p <- with_hypothesis(utils::modifyList(base, sc[!grepl("_override$", names(sc))]), state)
  if (!is.null(sc$kappa_override) && p$kappa > 0) p$kappa <- sc$kappa_override
  if (!is.null(sc$gamma_override) && p$gamma > 0) p$gamma <- sc$gamma_override
  p
}

grid <- expand.grid(scenario = names(scenarios), state = states, rep = seq_len(R),
                    stringsAsFactors = FALSE)
grid$seed <- seq_len(nrow(grid)) + 5e6

# SESOI and pseudo-true values for each scenario (the conversion uses that
# scenario's own variance components).
setup <- list()
for (sc in names(scenarios)) for (st in states) {
  p <- scenario_pars(scenarios[[sc]], st)
  setup[[paste(sc, st)]] <- list(pars = p, sesoi = sesoi_on_estimator_scale(p),
                                 truth = pseudo_true(p)$uv)
}

one_rep <- function(i) {
  g <- grid[i, ]
  su <- setup[[paste(g$scenario, g$state)]]
  res <- analyze_study(simulate_study(su$pars, seed = g$seed), su$sesoi, with_ev = FALSE)
  tt <- res$tests
  z <- qnorm(0.975)
  tt_two <- run_tests(transform(tt[, c("test", "est", "direction")], se = tt$se_two_way), su$sesoi)
  data.frame(scenario = g$scenario, state = g$state, rep = g$rep,
             verdict_pid = res$verdict_pid, verdict_two_way = res$verdict_two_way,
             test = tt$test, truth = unname(su$truth[tt$test]), est = tt$est,
             se_pid = tt$se_pid, se_two_way = tt$se_two_way,
             equiv_pid = tt$equiv, equiv_two_way = tt_two$equiv,
             cover_pid = abs(tt$est - su$truth[tt$test]) <= z * tt$se_pid,
             cover_two_way = abs(tt$est - su$truth[tt$test]) <= z * tt$se_two_way)
}

# SIM_REUSE=1 summarizes the saved simulations instead of rerunning them.
rds <- out_path("tables", "sensitivity_simulations.rds")
if (identical(Sys.getenv("SIM_REUSE"), "1") && file.exists(rds)) {
  sims <- readRDS(rds)
  R <- max(sims$rep)
} else {
  t0 <- Sys.time()
  sims <- do.call(rbind, future_lapply(seq_len(nrow(grid)), one_rep,
                                       future.seed = TRUE, future.packages = "ordinal"))
  cat(sprintf("%d studies in %.1f minutes\n", nrow(grid), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  saveRDS(sims, rds)
}

studies <- unique(sims[, c("scenario", "state", "rep", "verdict_pid", "verdict_two_way")])
studies$correct_pid <- studies$verdict_pid == correct_verdict[studies$state]
studies$correct_two_way <- studies$verdict_two_way == correct_verdict[studies$state]
verdicts <- aggregate(cbind(correct_pid, correct_two_way) ~ scenario + state, studies, mean)

cov <- aggregate(cbind(cover_pid, cover_two_way) ~ scenario + state, sims, mean,
                 na.action = na.pass, na.rm = TRUE)
by_test <- aggregate(cbind(cover_pid, cover_two_way) ~ scenario + test,
                     sims, mean, na.action = na.pass, na.rm = TRUE)
by_test$two_way_na <- aggregate(cbind(na = is.na(cover_two_way)) ~ scenario + test, sims, sum,
                                na.action = na.pass)$na
by_test <- by_test[order(match(by_test$scenario, names(scenarios)), by_test$test), ]
write.csv(by_test, out_path("tables", "sensitivity_coverage_by_test.csv"), row.names = FALSE)
worst <- aggregate(cbind(min_cover_pid = cover_pid, min_cover_two_way = cover_two_way) ~ scenario,
                   by_test, min)

summary_tab <- merge(merge(verdicts, cov), worst)
summary_tab$n_per_cell <- sapply(paste(summary_tab$scenario, summary_tab$state), function(k) setup[[k]]$pars$n_per_cell)
summary_tab$sesoi <- sapply(paste(summary_tab$scenario, summary_tab$state), function(k) setup[[k]]$sesoi)
summary_tab$mcse_correct <- sqrt(summary_tab$correct_pid * (1 - summary_tab$correct_pid) / R)
summary_tab$mcse_correct_two_way <- sqrt(summary_tab$correct_two_way * (1 - summary_tab$correct_two_way) / R)
summary_tab <- summary_tab[order(match(summary_tab$scenario, names(scenarios)), summary_tab$state), ]
print(summary_tab, digits = 3, row.names = FALSE)
write.csv(summary_tab, out_path("tables", "sensitivity_summary.csv"), row.names = FALSE)

# Request-slope sweep: participant-clustered against two-way clustered errors,
# by sample size.
sweep <- summary_tab[grepl("^request_slopes", summary_tab$scenario),
                     c("n_per_cell", "state", "correct_pid", "correct_two_way",
                       "cover_pid", "cover_two_way", "min_cover_pid", "min_cover_two_way")]
sweep$requests <- sapply(paste(summary_tab$scenario, summary_tab$state)[grepl("^request_slopes", summary_tab$scenario)],
                         function(k) setup[[k]]$pars$n_safe + setup[[k]]$pars$n_unsafe)
sweep <- sweep[order(sweep$state, sweep$requests, sweep$n_per_cell), ]
print(sweep, digits = 3, row.names = FALSE)
write.csv(sweep, out_path("tables", "sensitivity_request_slopes_sweep.csv"), row.names = FALSE)

# Size of the equivalence tests: when the true gain in d equals the SESOI
# (gamma_at_sesoi, vigilance state), each T3 contrast should be declared
# equivalent in at most 5% of studies.
at_sesoi <- sims[sims$scenario == "gamma_at_sesoi" & sims$state == "vigilance" &
                   sims$test %in% c("T3_plain", "T3_tactic"), ]
size_tab <- aggregate(cbind(equiv_pid, equiv_two_way) ~ test, at_sesoi, mean)
# The same, with the wider bound a pilot gives if it underestimates tau_c
# (0.3 instead of 0.5), or with the 0.20 used unconverted.
tost_rate <- function(se, bound) {
  tapply(run_tests(data.frame(test = at_sesoi$test, est = at_sesoi$est, se = se,
                              direction = "greater"), bound)$equiv, at_sesoi$test, mean)[size_tab$test]
}
wide <- c(pilot_tau_c_0.3 = sesoi_on_estimator_scale(utils::modifyList(base, list(tau_c = 0.3))),
          no_conversion = 0.20)
for (b in names(wide)) {
  size_tab[[paste0("equiv_pid_", b)]] <- unname(tost_rate(at_sesoi$se_pid, wide[[b]]))
  size_tab[[paste0("equiv_two_way_", b)]] <- unname(tost_rate(at_sesoi$se_two_way, wide[[b]]))
}
size_tab$mcse_at_5pct <- sqrt(0.05 * 0.95 / R)
print(size_tab, digits = 3, row.names = FALSE)
write.csv(size_tab, out_path("tables", "sensitivity_equivalence_at_sesoi.csv"), row.names = FALSE)

plan(sequential)
writeLines(capture.output(sessionInfo()), out_path("sessionInfo_06.txt"))
