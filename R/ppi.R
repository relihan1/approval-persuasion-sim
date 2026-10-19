# ppi.R
# Can simulated reviewers stand in for some of the people?
#
# Prediction-powered inference (PPI; Angelopoulos et al., 2023) combines a
# small labeled sample (people who did the task and also have a simulator
# prediction) with a large unlabeled sample (simulator predictions only). The
# labeled sample measures the simulator's error and corrects for it, so the
# interval stays valid however biased the simulator is. PPI++ (Angelopoulos,
# Duchi, & Zrnic, 2023, preprint) scales the correction by lambda and so never
# does worse than the people-only interval asymptotically.

# ---- Estimators of a population mean ---------------------------------------

classical_mean <- function(y, level = 0.95) {
  est <- mean(y); se <- stats::sd(y) / sqrt(length(y))
  ci_row(est, se, level)
}

simulator_only_mean <- function(f_unl, level = 0.95) {
  est <- mean(f_unl); se <- stats::sd(f_unl) / sqrt(length(f_unl))
  ci_row(est, se, level)
}

# lambda = 1 is PPI; lambda = 0 returns the classical estimate.
ppi_mean <- function(y_lab, f_lab, f_unl, lambda = 1, level = 0.95) {
  n <- length(y_lab); N <- length(f_unl)
  est <- mean(y_lab) + lambda * (mean(f_unl) - mean(f_lab))
  se  <- sqrt(stats::var(y_lab - lambda * f_lab) / n + lambda^2 * stats::var(f_unl) / N)
  out <- ci_row(est, se, level)
  out$lambda <- lambda
  out
}

# Plug-in lambda for the mean (Angelopoulos, Duchi, & Zrnic, 2023, Example 6.1):
# Cov_n(Y, f) / ((1 + n/N) Var_{n+N}(f)).
ppi_lambda_hat <- function(y_lab, f_lab, f_unl) {
  n <- length(y_lab); N <- length(f_unl)
  stats::cov(y_lab, f_lab) / ((1 + n / N) * stats::var(c(f_lab, f_unl)))
}

ci_row <- function(est, se, level) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  data.frame(est = est, se = se, lo = est - z * se, hi = est + z * se)
}

# ---- People and simulated reviewers ----------------------------------------
# Target: in the control group on the argument screen, (a) the share of
# tactic-laden unsafe requests a person approves, and (b) the threat effect,
# that share minus the share of plain unsafe requests approved.
#
# A simulated reviewer is built from a profile of the person. The profile
# predicts the person's criterion, sensitivity, and tactic susceptibility with
# accuracy r.
# Three failure modes reported for LLM simulations are built in: a cautious
# lean toward rejecting, weak sensitivity to wording, and person-level
# predictions limited by what a profile can reveal.

simulator_defaults <- function() {
  list(caution = 0.3,       # simulated reviewers lean toward reject (latent units)
       wording = 0.2,       # they respond to tactics at 20% of human strength
       r_profile = 0.5)     # profile-to-person accuracy
}

ppi_items <- function(pars) {
  items <- make_items(pars)
  items$v_c <- stats::rnorm(nrow(items), 0, pars$omega_c)
  items$v_d <- stats::rnorm(nrow(items), 0, pars$omega_d)
  items[items$unsafe == 1L, ]
}

draw_people <- function(n, pars) {
  data.frame(u_c = stats::rnorm(n, 0, pars$tau_c),
             u_d = stats::rnorm(n, 0, pars$tau_d),
             u_t = stats::rnorm(n, 0, pars$tau_t),
             list = sample(c("A", "B"), n, replace = TRUE))
}

# Probability of approving each unsafe request (rows: people, cols: requests),
# control group, argument screen, mid-session.
approve_prob <- function(people, items, pars, caution = 0, wording = 1,
                         u_c = people$u_c, u_d = people$u_d, u_t = people$u_t) {
  K <- pars$n_safe + pars$n_unsafe
  tac <- outer(people$list == "A", items$tactic_A == 1L, "==") * 1  # 1 = tactic version
  m <- outer(u_c + u_d + caution, items$v_c + items$v_d + pars$d0 -
               pars$palter * items$palter, "+") -
    tac * wording * (pars$pull + pars$pull_extra + u_t) -
    pars$drift * ((K + 1) / 2 - 1)
  list(p = stats::pnorm((pars$cuts[3] - m) / pars$sigma_unsafe), tactic = tac)
}

person_targets <- function(p, tac) {
  level  <- rowSums(p * tac) / rowSums(tac)
  threat <- level - rowSums(p * (1 - tac)) / rowSums(1 - tac)
  cbind(level = level, threat = threat)
}

simulate_people <- function(people, items, pars) {
  ap <- approve_prob(people, items, pars)
  y <- matrix(stats::rbinom(length(ap$p), 1, ap$p), nrow = nrow(ap$p))
  person_targets(y, ap$tactic)
}

simulate_reviewers <- function(people, items, pars, sim = simulator_defaults()) {
  r <- sim$r_profile
  guess <- function(u, tau) r^2 * u + tau * r * sqrt(1 - r^2) * stats::rnorm(length(u))
  ap <- approve_prob(people, items, pars, caution = sim$caution, wording = sim$wording,
                     u_c = guess(people$u_c, pars$tau_c),
                     u_d = guess(people$u_d, pars$tau_d),
                     u_t = guess(people$u_t, pars$tau_t))
  person_targets(ap$p, ap$tactic)  # expected value: the average of many simulator runs
}

# The population value each estimator is trying to recover (for these requests).
true_targets <- function(items, pars, n = 2e5) {
  ppl <- draw_people(n, pars)
  ap <- approve_prob(ppl, items, pars)
  colMeans(person_targets(ap$p, ap$tactic))
}

# One replicate: n labeled people with simulator predictions, N profiles with
# simulator predictions only.
ppi_replicate <- function(items, pars, sim, n = 200, N = 2000, truth) {
  lab <- draw_people(n, pars); unl <- draw_people(N, pars)
  y  <- simulate_people(lab, items, pars)
  fl <- simulate_reviewers(lab, items, pars, sim)
  fu <- simulate_reviewers(unl, items, pars, sim)
  out <- lapply(c("level", "threat"), function(target) {
    yy <- y[, target]; a <- fl[, target]; b <- fu[, target]
    lam <- ppi_lambda_hat(yy, a, b)
    res <- rbind(cbind(method = "People only", classical_mean(yy)),
                 cbind(method = "Simulator only", simulator_only_mean(b)),
                 cbind(method = "PPI", ppi_mean(yy, a, b, 1)[, 1:4]),
                 cbind(method = "PPI++", ppi_mean(yy, a, b, lam)[, 1:4]))
    res$target <- target
    res$truth <- truth[[target]]
    res$rho <- stats::cor(yy, a)
    res
  })
  do.call(rbind, out)
}
