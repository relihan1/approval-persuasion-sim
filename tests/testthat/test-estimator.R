loglik_at <- function(par, fit, data) {
  J <- length(fit$alpha); nb <- length(fit$beta)
  alpha <- par[seq_len(J)]; beta <- par[J + seq_len(nb)]; zeta <- par[-seq_len(J + nb)]
  X <- stats::model.matrix(location_terms, data)[, names(fit$beta)]
  eta <- drop(X %*% beta); sig <- exp(data$unsafe_s * zeta)
  th <- c(-Inf, alpha, Inf); y <- as.integer(data$y)
  sum(data$n * log(pnorm((th[y + 1] - eta) / sig) - pnorm((th[y] - eta) / sig)))
}

small <- function(seed = 3, ...) {
  p <- utils::modifyList(dgm_defaults(), list(n_per_cell = 40, ...))
  simulate_study(p, seed = seed)
}

test_that("pattern-level fit reproduces the row-level fit", {
  d <- small()
  a <- aggregate_reviews(d)
  f_agg <- fit_uv(a$agg)
  rows <- transform(d, y = factor(rating, levels = 1:6, ordered = TRUE), unsafe_s = unsafe, trial_c = block_c)
  f_row <- ordinal::clm(y ~ unsafe * tactic * trained * screen + trial_c, scale = ~ unsafe_s,
                        data = rows, link = "probit")
  expect_equal(unname(coef(f_agg)), unname(coef(f_row)), tolerance = 1e-6)
  expect_equal(as.numeric(logLik(f_agg)), as.numeric(logLik(f_row)), tolerance = 1e-6)
  # The bread of the sandwich: the weighted fit's covariance is the row-level one.
  expect_lt(max(abs(vcov(f_agg) - vcov(f_row)) / pmax(abs(vcov(f_row)), 1e-8)), 1e-5)
})

test_that("analytic scores match numerical derivatives and sum to zero at the MLE", {
  d <- small(seed = 4)
  a <- aggregate_reviews(d)
  f <- fit_uv(a$agg)
  par <- coef(f)
  expect_equal(loglik_at(par, f, a$agg), as.numeric(logLik(f)), tolerance = 1e-8)
  expect_lt(max(abs(colSums(clm_scores(f, a$agg) * a$agg$n))), 1e-4)

  set.seed(9)
  par2 <- par + rnorm(length(par), 0, 0.02)
  f2 <- f
  J <- length(f$alpha); nb <- length(f$beta)
  f2$alpha[] <- par2[seq_len(J)]; f2$beta[] <- par2[J + seq_len(nb)]; f2$zeta[] <- par2[-seq_len(J + nb)]
  analytic <- colSums(clm_scores(f2, a$agg) * a$agg$n)
  h <- 1e-5
  numeric <- sapply(seq_along(par2), function(k) {
    e <- replace(numeric(length(par2)), k, h)
    (loglik_at(par2 + e, f, a$agg) - loglik_at(par2 - e, f, a$agg)) / (2 * h)
  })
  # Component by component, so a sign error in one small score can't hide.
  expect_lt(max(abs(unname(analytic) - numeric) / pmax(abs(numeric), 1)), 1e-4)

  # Row-level scores (what the sandwich sums by cluster), for a few patterns.
  rows <- c(1, 50, 200, nrow(a$agg))
  for (r in rows) {
    one <- a$agg[r, ]; one$n <- 1
    num_r <- sapply(seq_along(par2), function(k) {
      e <- replace(numeric(length(par2)), k, h)
      (loglik_at(par2 + e, f, one) - loglik_at(par2 - e, f, one)) / (2 * h)
    })
    expect_lt(max(abs(unname(clm_scores(f2, one)[1, ]) - num_r)), 1e-5)
  }
})

test_that("cluster-robust and model-based errors agree when reviews are independent", {
  d <- small(seed = 5, tau_c = 0, tau_d = 0, omega_c = 0, omega_d = 0, n_per_cell = 150)
  a <- aggregate_reviews(d)
  f <- fit_uv(a$agg)
  se_model  <- sqrt(diag(vcov(f)))
  se_robust <- sqrt(diag(vcov_cluster(f, a$agg, a$key, d$pid)))
  ratio <- se_robust / se_model
  expect_true(all(ratio > 0.8 & ratio < 1.25))
})
