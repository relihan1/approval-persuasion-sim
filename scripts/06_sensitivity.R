# 06_sensitivity.R
# How much do the conclusions depend on assumptions the defaults don't vary?
# At 125 participants per cell (the size script 02 points to), each scenario
# changes one thing and reruns the cue, vigilance, and neither states:
#   - requests differ in how much their tactics pull and in how strongly
#     trained reviewers react to them (request-level random slopes), also
#     at 200 participants per cell
#   - a smaller cue effect (kappa = 0.25)
#   - a vigilance gain equal to the smallest effect of interest (gamma = 0.2)
#   - larger sensitivity variance across participants and requests (0.5)
#   - unequal variance stronger (sigma = 1.5) or absent (sigma = 1.0)
# The verdict is recorded under participant-clustered and two-way clustered
# errors, and coverage is checked against each scenario's own pseudo-true values.

source(here::here("R", "load.R"))
library(future.apply)

base <- dgm_defaults()
base$n_per_cell <- 125
R <- as.integer(Sys.getenv("SIM_R6", if (QUICK) 150 else 1000))
plan(multisession, workers = as.integer(Sys.getenv("SIM_WORKERS", future::availableCores())))

scenarios <- list(
  default            = list(),
  request_slopes     = list(omega_t = 0.3, omega_k = 0.2),
  request_slopes_200 = list(omega_t = 0.3, omega_k = 0.2, n_per_cell = 200),
  kappa_0.25         = list(kappa_override = 0.25),
  gamma_at_sesoi     = list(gamma_override = 0.2),
  sensitivity_sd_0.5 = list(tau_d = 0.5, omega_d = 0.5),
  sigma_1.5          = list(sigma_unsafe = 1.5),
  sigma_1.0          = list(sigma_unsafe = 1.0)
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
  data.frame(scenario = g$scenario, state = g$state, rep = g$rep,
             verdict_pid = res$verdict_pid, verdict_two_way = res$verdict_two_way,
             test = tt$test,
             cover_pid = abs(tt$est - su$truth[tt$test]) <= z * tt$se_pid,
             cover_two_way = abs(tt$est - su$truth[tt$test]) <= z * tt$se_two_way)
}

t0 <- Sys.time()
sims <- do.call(rbind, future_lapply(seq_len(nrow(grid)), one_rep,
                                     future.seed = TRUE, future.packages = "ordinal"))
cat(sprintf("%d studies in %.1f minutes\n", nrow(grid), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
saveRDS(sims, out_path("tables", "sensitivity_simulations.rds"))

studies <- unique(sims[, c("scenario", "state", "rep", "verdict_pid", "verdict_two_way")])
studies$correct_pid <- studies$verdict_pid == correct_verdict[studies$state]
studies$correct_two_way <- studies$verdict_two_way == correct_verdict[studies$state]
verdicts <- aggregate(cbind(correct_pid, correct_two_way) ~ scenario + state, studies, mean)

cov <- aggregate(cbind(cover_pid, cover_two_way) ~ scenario + state, sims, mean,
                 na.action = na.pass, na.rm = TRUE)
worst <- aggregate(cbind(min_cover_pid = cover_pid, min_cover_two_way = cover_two_way) ~ scenario + test,
                   sims, mean, na.action = na.pass, na.rm = TRUE)
worst <- aggregate(cbind(min_cover_pid, min_cover_two_way) ~ scenario, worst, min)

summary_tab <- merge(merge(verdicts, cov), worst)
summary_tab$sesoi <- sapply(paste(summary_tab$scenario, summary_tab$state), function(k) setup[[k]]$sesoi)
summary_tab$mcse_correct <- sqrt(summary_tab$correct_pid * (1 - summary_tab$correct_pid) / R)
summary_tab <- summary_tab[order(match(summary_tab$scenario, names(scenarios)), summary_tab$state), ]
print(summary_tab, digits = 3, row.names = FALSE)
write.csv(summary_tab, out_path("tables", "sensitivity_summary.csv"), row.names = FALSE)

plan(sequential)
writeLines(capture.output(sessionInfo()), out_path("sessionInfo_06.txt"))
