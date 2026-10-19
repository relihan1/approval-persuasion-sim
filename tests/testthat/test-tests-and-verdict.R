mk <- function(t2, t3p, t3t, eq = c(TRUE, TRUE)) {
  data.frame(test = c("T1", "T2", "T3_plain", "T3_tactic", "T4", "T5"),
             sig = c(TRUE, t2, t3p, t3t, TRUE, TRUE),
             equiv = c(FALSE, FALSE, eq, FALSE, FALSE))
}

test_that("the decision rule maps patterns to verdicts", {
  expect_equal(verdict(mk(TRUE, FALSE, FALSE)), "cue")
  expect_equal(verdict(mk(FALSE, TRUE, FALSE, c(FALSE, TRUE))), "vigilance")
  expect_equal(verdict(mk(TRUE, FALSE, TRUE, c(TRUE, FALSE))), "both")
  expect_equal(verdict(mk(FALSE, FALSE, FALSE)), "null")
  expect_equal(verdict(mk(TRUE, FALSE, FALSE, c(TRUE, FALSE))), "inconclusive")
})

test_that("TOST and directional p-values behave", {
  ce <- data.frame(test = "x", est = 0, se = 0.05, direction = "greater")
  r <- run_tests(ce, sesoi = 0.2)
  expect_equal(r$p_dir, 0.5)
  expect_lt(r$p_tost, 0.001)
  ce$est <- 0.2
  expect_equal(run_tests(ce, sesoi = 0.2)$p_tost, 0.5)
})
