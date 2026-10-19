# 03_simulated_reviewers_ppi.R
# When can simulated reviewers cut the number of people a rerun needs?
#
# 200 people do the task; 2,000 more are only profiled and simulated. The
# simulator leans toward rejecting, responds weakly to wording, and predicts
# each person only as well as their profile allows (r). We estimate two
# control-group quantities for tactic-laden unsafe requests:
#   level:  the share a person approves
#   threat: that share minus the share of plain unsafe requests approved
# and compare four estimators: people only, simulator only, PPI, and PPI++.

source(here::here("R", "load.R"))
set.seed(20261019)

pars <- utils::modifyList(dgm_defaults(), list(tau_t = 0.3))  # people differ in susceptibility
items <- ppi_items(pars)  # one fixed set of 12 unsafe requests
truth <- true_targets(items, pars)
cat("Population values:", sprintf("%s = %.3f", names(truth), truth), "\n")

r_grid <- c(0, 0.25, 0.5, 0.75, 1)
R <- if (QUICK) 500 else 5000
n_lab <- 200; N_unl <- 2000

reps <- do.call(rbind, lapply(r_grid, function(r) {
  sim <- utils::modifyList(simulator_defaults(), list(r_profile = r))
  out <- do.call(rbind, replicate(R, ppi_replicate(items, pars, sim, n_lab, N_unl, truth),
                                  simplify = FALSE))
  out$r_profile <- r
  out
}))
saveRDS(reps, out_path("tables", "ppi_replicates.rds"))

reps$covered <- reps$lo <= reps$truth & reps$truth <= reps$hi
reps$width <- reps$hi - reps$lo
summ <- aggregate(cbind(coverage = covered, width = width, est = est, var_est = se^2, rho = rho) ~
                    method + target + r_profile, data = reps, FUN = mean)
summ$bias <- summ$est - truth[summ$target]

# Effective sample size: the number of people alone that would give the same
# variance, relative to the 200 actually tested.
base <- summ[summ$method == "People only", c("target", "r_profile", "var_est")]
summ <- merge(summ, setNames(base, c("target", "r_profile", "var_people")))
summ$effective_n <- n_lab * summ$var_people / summ$var_est
summ$gain_pct <- 100 * (summ$effective_n / n_lab - 1)
summ$gain_pct[summ$method == "Simulator only"] <- NA  # not valid, so no gain to report
summ$target_label <- ifelse(summ$target == "level",
                            "Approval of tactic-laden unsafe requests",
                            "Threat effect (tactic minus plain)")
summ <- summ[order(summ$target, summ$r_profile, summ$method), ]
print(summ[, c("target", "r_profile", "method", "rho", "bias", "coverage", "width", "gain_pct")],
      digits = 3, row.names = FALSE)
write.csv(summ, out_path("tables", "ppi_summary.csv"), row.names = FALSE)

save_fig(plot_ppi(summ), "fig4_simulated_reviewers_ppi", width = 7.5, height = 7)
writeLines(capture.output(sessionInfo()), out_path("sessionInfo_03.txt"))
