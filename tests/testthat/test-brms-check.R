# Synthetic draws in brms's naming scheme: 4 chains of 500 iterations.
fake_draws <- function(stuck = character(), seed = 1) {
  set.seed(seed)
  vars <- c("b_Intercept[1]", "b_unsafe", "b_tactic", "b_disc_unsafe", "sd_pid__Intercept",
            "sd_item__tactic", "cor_item__Intercept__tactic", "r_pid[1,Intercept]", "lp__")
  a <- array(rnorm(500 * 4 * length(vars)), c(500, 4, length(vars)),
             dimnames = list(NULL, NULL, vars))
  for (v in stuck) a[, , v] <- a[, , v] + rep(c(0, 0, 0, 3), each = 500)  # one chain elsewhere
  posterior::as_draws_array(a)
}

test_that("parameters are classed by brms name, and internals are left out", {
  v <- c("b_Intercept[2]", "b_unsafe:tactic", "b_disc_unsafe", "sd_item__tactic:trained",
         "cor_pid__Intercept__unsafe", "r_item[3,tactic]", "z_1[1,1]", "lp__", "Intercept[1]")
  expect_identical(brms_par_class(v), c("threshold", "b", "disc", "sd_item", "cor_pid",
                                        "r_item", NA, NA, NA))
})

test_that("a request SD that has not mixed fails the check, though every coefficient passes", {
  ok <- fake_draws()
  bad <- fake_draws(stuck = "sd_item__tactic")
  smp <- data.frame(chain = 1:4, divergent = 0, treedepth_hits = 0, ebfmi = 1)
  cb_ok <- convergence_by_class(ok)
  cb_bad <- convergence_by_class(bad)
  expect_true(convergence_ok(cb_ok, smp, n_chains = 4))
  expect_false(convergence_ok(cb_bad, smp, n_chains = 4))
  expect_lt(cb_bad$max_rhat[cb_bad$class == "b"], 1.01)  # the old b_-only check passes
  expect_identical(cb_bad$worst_rhat[cb_bad$class == "sd_item"], "sd_item__tactic")
  # Request and participant effects are summarised but not part of the rule.
  expect_true(convergence_ok(convergence_by_class(fake_draws(stuck = "r_pid[1,Intercept]")),
                             smp, n_chains = 4))
  smp$divergent[2] <- 1
  expect_false(convergence_ok(cb_ok, smp, n_chains = 4))
})

test_that("the ESS floor covers what the verdict uses; R-hat covers every class", {
  cb <- data.frame(class = c("b", "threshold", "disc", "sd_item", "cor_item", "contrast", "r_pid"),
                   max_rhat = 1.005, min_ess_bulk = 2000, min_ess_tail = 2000)
  smp <- data.frame(chain = 1:4, divergent = 0, treedepth_hits = 0, ebfmi = 1)
  with_row <- function(cls, col, value) { cb[cb$class == cls, col] <- value; cb }
  expect_true(convergence_ok(cb, smp, 4))
  expect_true(convergence_ok(with_row("threshold", "min_ess_bulk", 150), smp, 4))
  expect_true(convergence_ok(with_row("b", "min_ess_tail", 150), smp, 4))
  expect_false(convergence_ok(with_row("threshold", "max_rhat", 1.02), smp, 4))
  expect_false(convergence_ok(with_row("sd_item", "min_ess_bulk", 150), smp, 4))
  expect_false(convergence_ok(with_row("contrast", "min_ess_tail", 399), smp, 4))
  expect_true(convergence_ok(with_row("r_pid", "max_rhat", 1.05), smp, 4))
  smp$ebfmi[3] <- 0.2
  expect_false(convergence_ok(cb, smp, 4))
})

test_that("coefficients that vary only between participants or only between requests have no weight in any contrast", {
  L <- contrast_L()
  between <- c("unsafe", "palter", "trained", "screenparaphrase", "screenevidence",
               "trained:screenparaphrase", "trained:screenevidence")
  expect_true(all(between %in% colnames(L)))
  expect_true(all(L[, between] == 0))
})

test_that("contrast draws keep the chains and match L %*% b draw by draw", {
  L <- rbind(T1 = c(unsafe = 1, tactic = -1), T2 = c(unsafe = 0.5, tactic = 0.5))
  d <- fake_draws(stuck = "b_tactic")
  cd <- contrast_draws(d, L)
  expect_identical(posterior::variables(cd), c("T1", "T2"))
  expect_identical(posterior::nchains(cd), 4L)
  m <- posterior::as_draws_matrix(d)
  expect_equal(as.vector(posterior::as_draws_matrix(cd)[, "T1"]),
               as.vector(m[, "b_unsafe"] - m[, "b_tactic"]))
  cb <- convergence_by_class(d, cd)
  expect_gt(cb$max_rhat[cb$class == "contrast"], 1.01)  # the stuck chain carries over
})

test_that("sampler summary counts divergences and tree-depth hits and computes E-BFMI", {
  set.seed(3)
  n <- 200
  np <- rbind(
    data.frame(Chain = rep(1:2, each = n), Iteration = rep(1:n, 2), Parameter = "divergent__",
               Value = c(rep(0, n), rep(0:1, length.out = n))),
    data.frame(Chain = rep(1:2, each = n), Iteration = rep(1:n, 2), Parameter = "treedepth__",
               Value = c(rep(10, 5), rep(6, n - 5), rep(6, n))),
    data.frame(Chain = rep(1:2, each = n), Iteration = rep(1:n, 2), Parameter = "energy__",
               Value = c(rnorm(n), cumsum(rnorm(n, 0, 0.1)))))
  s <- sampler_summary(np, max_treedepth = 10)
  expect_identical(s$divergent, c(0, n / 2))
  expect_identical(s$treedepth_hits, c(5L, 0L))
  expect_gt(s$ebfmi[1], 1)    # independent energies: E-BFMI near 2
  expect_lt(s$ebfmi[2], 0.3)  # a random walk in energy: low E-BFMI
})
