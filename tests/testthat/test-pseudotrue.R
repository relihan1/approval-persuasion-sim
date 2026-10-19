pars <- dgm_defaults()

test_that("with no training effect, training contrasts are exactly zero", {
  pt <- pseudo_true(with_hypothesis(pars, "null"))
  expect_lt(max(abs(pt$uv[c("T2", "T3_plain", "T3_tactic")])), 1e-3)
})

test_that("a pure criterion shift leaves d unchanged in the unequal-variance model only", {
  pt <- pseudo_true(with_hypothesis(pars, "cue"))
  expect_lt(max(abs(pt$uv[c("T3_plain", "T3_tactic")])), 0.005)
  expect_lt(pt$ev[["T3_tactic"]], -0.03)  # the equal-variance artifact
})

test_that("a participant-level gain in d is attenuated on the marginal scale as expected", {
  pt <- pseudo_true(with_hypothesis(pars, "vigilance"))
  expected <- 0.3 / sqrt(1 + pars$tau_c^2 + pars$omega_c^2)
  expect_equal(unname(pt$uv[["T3_plain"]]), expected, tolerance = 0.01)
})
