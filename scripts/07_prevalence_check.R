# 07_prevalence_check.R
# The design shows each participant 12 unsafe requests out of 48 (25%), far
# above the rate at which unsafe requests reach approvers in deployment. This
# script asks how much the analysis depends on that: it reruns the cue,
# vigilance, and neither states with 6 unsafe out of 48 (12.5%; 42 safe), at
# 125 participants per cell, with every other value at its default. The 25%
# comparison is script 06's default scenario (same sample size, same errors,
# 1,000 studies per state).
#
# With 6 unsafe requests, each paltering stratum has 3, an odd number, so the
# two counterbalancing lists can't split them evenly: list A shows 4 of the 6
# with tactics and list B shows 2. Across participants each request still
# appears equally often in each version.
#
# The per-study analysis mirrors one_rep() in script 06, so the two sets of
# results are directly comparable. Results go to output/tables/prevalence_*.

source(here::here("R", "load.R"))
library(future.apply)

base <- utils::modifyList(dgm_defaults(), list(n_per_cell = 125, n_safe = 42, n_unsafe = 6))
# 1,100 studies per state keep the 95% Monte Carlo half-width of a verdict
# rate under 3 points at any rate (1.96 * sqrt(.25 / 1100) = .030).
R <- as.integer(Sys.getenv("SIM_R7", if (QUICK) 150 else 1100))
plan(multisession, workers = as.integer(Sys.getenv("SIM_WORKERS", future::availableCores())))

states <- c("cue", "vigilance", "null")
grid <- expand.grid(state = states, rep = seq_len(R), stringsAsFactors = FALSE)
grid$seed <- seq_len(nrow(grid)) + 7e6

setup <- lapply(stats::setNames(states, states), function(st) {
  p <- with_hypothesis(base, st)
  list(pars = p, sesoi = sesoi_on_estimator_scale(p), truth = pseudo_true(p)$uv)
})

boot_draws <- as.integer(Sys.getenv("SIM_BOOT", 999))

one_rep <- function(i) {
  g <- grid[i, ]
  su <- setup[[g$state]]
  res <- analyze_study(simulate_study(su$pars, seed = g$seed), su$sesoi, with_ev = FALSE,
                       satterthwaite = TRUE, boot_draws = boot_draws)
  tt <- res$tests
  z <- qnorm(0.975)
  truth <- unname(su$truth[tt$test])
  inside <- function(t) if (is.null(t)) NA else t$lo95 <= truth & truth <= t$hi95
  data.frame(state = g$state, rep = g$rep, test = tt$test, truth = truth, est = tt$est,
             verdict_pid = res$verdict_pid, verdict_two_way_jk = res$verdict_two_way_jk,
             verdict_two_way_jk_sw = res$verdict_two_way_jk_sw,
             verdict_two_way_boot = if (boot_draws > 0) res$verdict_two_way_boot else NA_character_,
             sig_pid = tt$sig, equiv_pid = tt$equiv,
             sig_two_way_boot = if (boot_draws > 0) res$tests_two_way_boot$sig else NA,
             equiv_two_way_boot = if (boot_draws > 0) res$tests_two_way_boot$equiv else NA,
             cover_pid = abs(tt$est - truth) <= z * tt$se_pid,
             cover_two_way_jk = abs(tt$est - truth) <= z * tt$se_two_way_jk,
             cover_two_way_jk_sw = inside(res$tests_two_way_jk_sw),
             cover_two_way_boot = inside(res$tests_two_way_boot))
}

rds <- out_path("tables", "prevalence_simulations.rds")
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

kinds <- c("pid", "two_way_jk", "two_way_jk_sw", "two_way_boot")

# Correct verdicts by state, beside script 06's default scenario (25% unsafe).
studies <- unique(sims[, c("state", "rep", paste0("verdict_", kinds))])
for (k in kinds) studies[[paste0("correct_", k)]] <- studies[[paste0("verdict_", k)]] == correct_verdict[studies$state]
verdicts <- aggregate(studies[paste0("correct_", kinds)], studies["state"], mean)
verdicts$unsafe_share <- "12.5%"
s06 <- read.csv(out_path("tables", "sensitivity_summary.csv"))
s06 <- s06[s06$scenario == "default", c("state", paste0("correct_", kinds))]
s06$unsafe_share <- "25%"
verdict_tab <- rbind(s06, verdicts[names(s06)])
verdict_tab$studies <- ifelse(verdict_tab$unsafe_share == "25%", 1000, R)
verdict_tab <- verdict_tab[order(match(verdict_tab$state, states), verdict_tab$unsafe_share), ]
print(verdict_tab, digits = 3, row.names = FALSE)
write.csv(verdict_tab, out_path("tables", "prevalence_correct_verdicts.csv"), row.names = FALSE)

# Coverage of 95% intervals for T1 to T3, pooled over the three states.
t13 <- c("T1", "T2", "T3_plain", "T3_tactic")
cov12 <- aggregate(sims[sims$test %in% t13, paste0("cover_", kinds)],
                   sims[sims$test %in% t13, "test", drop = FALSE], mean, na.rm = TRUE)
cov12$unsafe_share <- "12.5%"
c06 <- read.csv(out_path("tables", "sensitivity_coverage_by_test.csv"))
c06 <- c06[c06$scenario == "default" & c06$test %in% t13, c("test", paste0("cover_", kinds))]
c06$unsafe_share <- "25%"
cover_tab <- rbind(c06, cov12[names(c06)])
cover_tab <- cover_tab[order(cover_tab$test, cover_tab$unsafe_share), ]
cover_tab$mcse_at_95 <- sqrt(0.95 * 0.05 / (length(states) * ifelse(cover_tab$unsafe_share == "25%", 1000, R)))
print(cover_tab, digits = 3, row.names = FALSE)
write.csv(cover_tab, out_path("tables", "prevalence_coverage.csv"), row.names = FALSE)

# Which tests move the verdicts: rejection and equivalence rates for T1 to T3
# by state, participant-clustered and bootstrap.
rates <- aggregate(sims[sims$test %in% t13, c("sig_pid", "equiv_pid", "sig_two_way_boot", "equiv_two_way_boot")],
                   sims[sims$test %in% t13, c("state", "test")], mean, na.rm = TRUE)
rates <- rates[order(match(rates$state, states), rates$test), ]
print(rates, digits = 3, row.names = FALSE)
write.csv(rates, out_path("tables", "prevalence_test_rates.csv"), row.names = FALSE)

truth_tab <- do.call(rbind, lapply(states, function(st) data.frame(state = st, t(round(setup[[st]]$truth, 4)))))
write.csv(truth_tab, out_path("tables", "prevalence_pseudo_true.csv"), row.names = FALSE)

plan(sequential)
writeLines(capture.output(sessionInfo()), out_path("sessionInfo_07.txt"))
