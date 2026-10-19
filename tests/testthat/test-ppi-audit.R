test_that("PPI with lambda = 0 is the people-only estimate", {
  set.seed(1)
  y <- rnorm(50); f <- y + rnorm(50); fu <- rnorm(500)
  expect_equal(ppi_mean(y, f, fu, lambda = 0)$est, classical_mean(y)$est)
  expect_equal(ppi_mean(y, f, fu, lambda = 0)$se, classical_mean(y)$se)
})

test_that("the PPI++ lambda follows Cov(Y, f) / ((1 + n/N) Var(f))", {
  set.seed(2)
  y <- rnorm(100); fu <- rnorm(1000)
  expect_equal(ppi_lambda_hat(y, y, fu), var(y) / ((1 + 100 / 1000) * var(c(y, fu))))
})

test_that("AUC, folds, and the ridge classifier behave", {
  expect_equal(auc(1:4, c(0, 0, 1, 1)), 1)
  expect_equal(auc(4:1, c(0, 0, 1, 1)), 0)
  expect_equal(auc(rep(1, 4), c(0, 0, 1, 1)), 0.5)
  bank <- make_text_bank(48, seed = 1)
  folds <- grouped_folds(bank$item, bank$unsafe, k = 4)
  expect_true(all(tapply(folds, bank$item, function(f) length(unique(f))) == 1))
  set.seed(3)
  X <- matrix(rnorm(300), 100); y <- rbinom(100, 1, plogis(X %*% c(1, -1, 0.5)))
  expect_equal(predict_ridge(ridge_logit(X, y, lambda = 0), X),
               unname(predict(glm(y ~ X, family = binomial))), tolerance = 1e-5)
})
