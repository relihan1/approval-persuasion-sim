test_that("each participant reviews each request once, in a full random order", {
  p <- utils::modifyList(dgm_defaults(), list(n_per_cell = 4))
  d <- simulate_study(p, seed = 1)
  K <- p$n_safe + p$n_unsafe
  expect_equal(nrow(d), 6 * 4 * K)
  expect_true(all(table(d$pid, d$item) == 1))
  expect_true(all(tapply(d$trial, d$pid, function(t) identical(sort(t), seq_len(K)))))
  expect_true(all(d$rating %in% 1:6))
  expect_identical(d$reject, as.integer(d$rating >= 4))
})

test_that("tactic versions are counterbalanced within participants and within requests", {
  p <- utils::modifyList(dgm_defaults(), list(n_per_cell = 4))
  d <- simulate_study(p, seed = 2)
  by_person <- aggregate(tactic ~ pid + unsafe + palter, d, mean)
  expect_true(all(by_person$tactic == 0.5))
  by_item <- aggregate(tactic ~ item + trained + screen, d, mean)
  expect_true(all(by_item$tactic == 0.5))
})

test_that("the population grid is a probability distribution", {
  pop <- population_grid(dgm_defaults())
  expect_equal(sum(pop$w), 1, tolerance = 1e-10)
  expect_true(all(pop$w >= 0))
})

test_that("the analytic marginal SD matches simulation, with every random slope switched on", {
  p <- utils::modifyList(dgm_defaults(), list(tau_t = 0.3, omega_t = 0.3, omega_k = 0.2))
  set.seed(7); n <- 4e5
  draw <- function(sd) rnorm(n, 0, sd)
  x <- latent_mean(1, 0, 1, 1, "paraphrase", 1, p,
                   u_c = draw(p$tau_c), u_d = draw(p$tau_d), u_t = draw(p$tau_t),
                   v_c = draw(p$omega_c), v_d = draw(p$omega_d), v_t = draw(p$omega_t),
                   v_k = draw(p$omega_k)) + latent_sd(1, p) * rnorm(n)
  expect_equal(sd(x), marginal_sd(1, 1, 1, "paraphrase", p), tolerance = 0.005)
})
