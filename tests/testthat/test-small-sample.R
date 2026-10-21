small_ss <- function(seed = 3, ...) {
  p <- utils::modifyList(dgm_defaults(), list(n_per_cell = 40, ...))
  simulate_study(p, seed = seed)
}

# Set the parameters of a copy of fit (thresholds, location, scale).
with_par <- function(fit, par) {
  J <- length(fit$alpha); nb <- length(fit$beta)
  fit$alpha[] <- par[seq_len(J)]; fit$beta[] <- par[J + seq_len(nb)]; fit$zeta[] <- par[-seq_len(J + nb)]
  fit
}

test_that("the Satterthwaite formula matches the exact moments of cluster-robust variances in a linear model", {
  # Linear model with N(0, 1) errors, so the working model is true and the
  # sandwich and jackknife variances are exact quadratic forms e'Ae:
  # E = tr(A), Var = 2 tr(A^2), df = 2 E^2 / Var.
  set.seed(1)
  G <- 9; n_g <- 6; N <- G * n_g
  cl <- rep(seq_len(G), each = n_g)
  X <- cbind(1, rnorm(N), rep(rbinom(G, 1, 0.4), each = n_g), rnorm(N) * (cl <= 3))
  cvec <- c(0, 0, 1, 0.5)
  H <- vapply(seq_len(G), function(g) crossprod(X[cl == g, , drop = FALSE]), matrix(0, 4, 4))
  XtXi <- solve(crossprod(X))
  resid <- diag(N) - X %*% XtXi %*% t(X)
  exact_df <- function(rows) {  # rows: one per cluster, c'(...) as a linear map of the errors
    A <- crossprod(rows)
    sum(diag(A))^2 / sum(A * A)
  }
  # Sandwich: c'B X_g' e_g, with e = (I - P) errors.
  sandwich_rows <- t(vapply(seq_len(G), function(g)
    drop(t(cvec) %*% XtXi %*% t(X[cl == g, , drop = FALSE]) %*% resid[cl == g, , drop = FALSE]), numeric(N)))
  expect_equal(satterthwaite_df(H, cvec), exact_df(sandwich_rows), tolerance = 1e-10)
  # Jackknife: refit without each cluster; the change in c'beta is linear in the errors.
  jack_rows <- t(vapply(seq_len(G), function(g) {
    keep <- cl != g
    fit_g <- solve(crossprod(X[keep, ]), t(X[keep, ]))  # beta_(g) = fit_g y[keep]
    full_map <- matrix(0, 4, N); full_map[, keep] <- fit_g
    drop(t(cvec) %*% (full_map - XtXi %*% t(X)))
  }, numeric(N)))
  expect_equal(satterthwaite_df(H, cvec, jackknife = TRUE), exact_df(jack_rows), tolerance = 1e-10)
  # One cluster carrying the fourth covariate's leverage: fewer df than clusters.
  expect_lt(satterthwaite_df(H, c(0, 0, 0, 1)), G - 1)
})

test_that("expected information matches the expected negative Hessian", {
  d <- small_ss(seed = 7)
  a <- aggregate_reviews(d)
  f <- fit_uv(a$agg)
  ei <- expected_info(f, a$agg, a$key)
  par <- coef(f); h <- 1e-5; K <- length(f$alpha) + 1L
  for (x in c(1, 17, dim(ei$info)[3])) {
    one <- a$agg[rep(a$key[match(x, ei$pattern)], K), ]  # a review with pattern x, once per rating
    one$y <- factor(seq_len(K), levels = seq_len(K), ordered = TRUE)
    p <- attr(clm_scores(f, one), "prob")
    # -sum_y p_y d s_y / d par, by central differences of the analytic scores.
    neg_hess <- -sapply(seq_along(par), function(k) {
      e <- replace(numeric(length(par)), k, h)
      colSums(p * (clm_scores(with_par(f, par + e), one) - clm_scores(with_par(f, par - e), one))) / (2 * h)
    })
    expect_lt(max(abs(ei$info[, , x] - neg_hess)), 1e-5)
  }
  # Summed over all reviews, close to the observed information at the MLE.
  H_all <- rowSums(cluster_info(ei, d$item), dims = 2)
  expect_equal(rowSums(cluster_info(ei, d$pid), dims = 2), H_all, tolerance = 1e-10)
  obs <- solve(vcov(f))
  expect_lt(max(abs(H_all - obs)) / max(abs(obs)), 0.02)
})

test_that("the linearized request jackknife behind the degrees of freedom tracks the refits", {
  # satterthwaite_df(jackknife = TRUE) assumes the change from leaving out
  # request g is -c'(I - B H_g)^-1 B U_g, U_g its summed scores at the estimate.
  d <- small_ss(seed = 6, n_per_cell = 30)
  a <- aggregate_reviews(d)
  f <- fit_uv(a$agg)
  L <- contrast_L()
  ei <- expected_info(f, a$agg, a$key)
  H <- cluster_info(ei, d$item)
  B <- solve(rowSums(H, dims = 2))
  U <- rowsum(clm_scores(f, a$agg)[a$key, ], d$item)  # requests x parameters, sorted like H
  C <- contrast_full(f, L)
  G <- dim(H)[3]
  lin <- sapply(seq_len(G), function(g) -drop(t(C) %*% solve(diag(nrow(B)) - B %*% H[, , g], B %*% U[g, ])))
  full <- drop(L %*% coef(f)[colnames(L)])
  items <- sort(unique(d$item))
  refit <- sapply(items, function(g) drop(L %*% coef(fit_uv(aggregate_reviews(d[d$item != g, ])$agg))[colnames(L)]) - full)
  for (j in seq_len(nrow(L))) expect_gt(cor(lin[j, ], refit[j, ]), 0.99)
  ratio <- rowSums(lin^2) / rowSums(refit^2)
  expect_true(all(ratio > 0.9 & ratio < 1.1))
})

test_that("two-way degrees of freedom scale the request part's by its share of the variance", {
  d <- small_ss(seed = 8)
  a <- aggregate_reviews(d)
  f <- fit_uv(a$agg)
  ei <- expected_info(f, a$agg, a$key)
  H <- cluster_info(ei, d$item)
  C <- contrast_full(f, contrast_L())
  nu <- apply(C, 2, function(cv) satterthwaite_df(H, cv))
  expect_true(all(nu > 1 & nu < length(unique(d$item))))
  v_req <- rep(0.01, 6)
  expect_equal(two_way_df(H, C, v_total = v_req, v_req = v_req), unname(nu))
  expect_equal(two_way_df(H, C, v_total = 2 * v_req, v_req = v_req), 4 * unname(nu))
  expect_true(all(is.na(two_way_df(H, C, v_total = c(-1, NA, rep(0.01, 4)), v_req = v_req)[1:2])))
})

test_that("score bootstrap draws equal a direct computation from recentred scores", {
  d <- small_ss(seed = 9, n_per_cell = 20)
  a <- aggregate_reviews(d)
  f <- fit_uv(a$agg)
  L <- contrast_L()
  ei <- expected_info(f, a$agg, a$key)
  set.seed(2)
  G <- length(unique(d$item))
  v <- matrix(sample(c(-1, 1), G * 3, replace = TRUE), G, 3)
  fast <- request_score_boot(f, a, d, L, ei = ei, v = v)

  B <- vcov(f)
  s <- clm_scores(f, a$agg)[a$key, ]
  g <- match(d$item, sort(unique(d$item)))
  one_way <- function(S, cl) {
    U <- rowsum(S, cl); k <- nrow(U)
    B %*% (crossprod(U) * k / (k - 1)) %*% B
  }
  C <- contrast_full(f, L)
  for (b in seq_len(ncol(v))) {
    S_star <- colSums(s * v[g, b])                                   # sum_g v_g S_g
    recentre <- t(apply(ei$info, 3, function(h) h %*% B %*% S_star))  # h_x B S*, one row per pattern
    s_star <- s * v[g, b] - recentre[ei$pattern, ]
    V <- one_way(s_star, d$pid) + one_way(s_star, d$item) -
      one_way(s_star, interaction(d$pid, d$item, drop = TRUE))
    direct <- drop(t(C) %*% B %*% S_star) / sqrt(diag(t(C) %*% V %*% C))
    expect_equal(unname(fast[, b]), unname(direct), tolerance = 1e-8)
  }
  # The recentred scores sum to about zero (exactly, if B were the inverse
  # expected rather than observed information), unlike the perturbed ones.
  expect_lt(max(abs(colSums(s_star))) / max(abs(S_star)), 0.05)
})

test_that("tests take t and bootstrap reference distributions", {
  ce <- data.frame(test = c("T1", "T2"), est = c(-0.1, 0.1), se = c(0.1, 0.05),
                   direction = c("less", "greater"))
  z <- ce$est / ce$se
  # The default is the normal reference, exactly as before the options existed.
  n <- run_tests(ce, sesoi = 0.2)
  expect_identical(n$p_dir, c(pnorm(z[1]), pnorm(z[2], lower.tail = FALSE)))
  expect_identical(n$p_tost, pmax(pnorm((ce$est + 0.2) / ce$se, lower.tail = FALSE), pnorm((ce$est - 0.2) / ce$se)))
  expect_identical(n$hi90, ce$est + qnorm(0.95) * ce$se)
  expect_equal(n$lo95, ce$est - qnorm(0.975) * ce$se)
  t5 <- run_tests(ce, sesoi = 0.2, df = c(5, 12))
  expect_equal(t5$p_dir, c(pt(z[1], 5), pt(z[2], 12, lower.tail = FALSE)))
  expect_equal(t5$hi95, ce$est + qt(0.975, c(5, 12)) * ce$se)
  expect_equal(t5$p_tost, pmax(pt((ce$est + 0.2) / ce$se, c(5, 12), lower.tail = FALSE),
                               pt((ce$est - 0.2) / ce$se, c(5, 12))))
  # Bootstrap draws: tail shares and percentile-t bounds, with NA draws dropped.
  tstar <- rbind(c(-3.5, -1, 0, 1, 2, NA), c(-2, -1, 0, 1, 2.5, 3))
  bt <- run_tests(ce, sesoi = 0.2, tstar = tstar)
  expect_equal(bt$p_dir, c(mean(c(-3.5, -1, 0, 1, 2) <= z[1]), mean(tstar[2, ] >= z[2])))
  expect_equal(bt$lo95[2], 0.1 - quantile(tstar[2, ], 0.975, names = FALSE) * 0.05)
  expect_equal(bt$hi95[2], 0.1 - quantile(tstar[2, ], 0.025, names = FALSE) * 0.05)
  # Many standard-normal draws reproduce the normal tests.
  set.seed(3)
  bz <- run_tests(ce, sesoi = 0.2, tstar = matrix(rnorm(2e5), 2))
  expect_equal(bz$p_dir, n$p_dir, tolerance = 0.05)
  expect_equal(bz$lo95, n$lo95, tolerance = 0.05)
})

test_that("the small-sample options add verdicts without changing the default ones", {
  d <- small_ss(seed = 10)
  sesoi <- 0.1
  base <- analyze_study(d, sesoi, with_ev = FALSE)
  set.seed(4)
  more <- analyze_study(d, sesoi, with_ev = FALSE, satterthwaite = TRUE, boot_draws = 199)
  for (k in c("verdict", "verdict_pid", "verdict_two_way", "verdict_two_way_jk")) expect_identical(more[[k]], base[[k]])
  expect_identical(more$tests[, names(base$tests)], base$tests)
  expect_true(all(c("verdict_two_way_sw", "verdict_two_way_jk_sw", "verdict_two_way_boot") %in% names(more)))
  expect_null(base$verdict_two_way_sw)
  # Same estimates and errors; t intervals are wider than normal ones.
  sw <- more$tests_two_way_jk_sw
  expect_equal(sw$se, base$tests_two_way_jk$se)
  ok <- !is.na(sw$se)
  expect_true(all((sw$hi95 - sw$lo95)[ok] > (base$tests_two_way_jk$hi95 - base$tests_two_way_jk$lo95)[ok]))
  expect_true(all(more$tests$df_two_way_jk[ok] > 1))
  expect_equal(more$tests_two_way_boot$se, more$tests$se_two_way)
  picked <- analyze_study(d, sesoi, with_ev = FALSE, clusters = "two_way_jk_sw", satterthwaite = TRUE)
  expect_identical(picked$verdict, more$verdict_two_way_jk_sw)
  expect_error(analyze_study(d, sesoi, clusters = "two_way_sw"))
  expect_error(analyze_study(d, sesoi, clusters = "two_way_boot"))
})
