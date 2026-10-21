# analyze.R
# The fast estimator used inside the power loop:
#   an unequal-variance ordinal probit signal detection model (ordinal::clm with
#   a scale term for unsafe requests), fitted to rating counts, with
#   cluster-robust standard errors from analytic score contributions.
#
# What it estimates: population-averaged (marginal) parameters, in units of
# the marginal safe-request SD. The preregistered primary model
# (scripts/05_brms_primary_model.R) estimates participant-level parameters
# with crossed random effects. Marginal values are smaller by a known factor,
# so the smallest effect of interest is converted to the marginal scale with
# sesoi_on_estimator_scale() rather than reused unchanged.

# ---- Collapse reviews to response patterns ---------------------------------
# With discrete covariates, the likelihood of the row-level data equals the
# frequency-weighted likelihood of the distinct patterns. Fitting at most 864
# weighted rows (safe requests never palter) instead of 43,200 reviews makes
# each fit take milliseconds. Trial position enters as the quarter of the
# session (block_c).

pattern_code <- function(d, outcome = "rating", n_out = 6L) {
  sc <- as.integer(d$screen) - 1L
  (((((d$unsafe * 2L + d$palter) * 2L + d$tactic) * 2L + d$trained) * 3L + sc) * 4L + (d$block - 1L)) * n_out +
    (d[[outcome]] - 1L) + 1L
}

pattern_table <- function(n_out = 6L) {
  g <- expand.grid(y = seq_len(n_out), block = 1:4, screen = screens, trained = 0:1,
                   tactic = 0:1, palter = 0:1, unsafe = 0:1, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  g <- g[, rev(names(g))]  # unsafe, palter, tactic, trained, screen, block, y (matches pattern_code order)
  g$screen <- factor(g$screen, levels = screens)
  g
}

aggregate_reviews <- function(dat, outcome = c("rating", "reject")) {
  outcome <- match.arg(outcome)
  n_out <- if (outcome == "rating") 6L else 2L
  d <- dat
  if (outcome == "reject") d$reject1 <- d$reject + 1L
  code <- pattern_code(d, if (outcome == "rating") "rating" else "reject1", n_out)
  tab <- pattern_table(n_out)
  tab$n <- tabulate(code, nbins = nrow(tab))
  keep <- tab$n > 0
  agg <- tab[keep, ]
  agg$trial_c  <- (agg$block - 2.5) * 0.25  # block midpoint in session units (as block_c)
  agg$unsafe_s <- agg$unsafe                # scale-term copy keeps coefficient names distinct
  agg$y <- factor(agg$y, levels = seq_len(n_out), ordered = TRUE)
  # Map each review to its row of agg, for cluster sums of score contributions.
  row_of <- integer(nrow(tab)); row_of[keep] <- seq_len(sum(keep))
  list(agg = agg, key = row_of[code])
}

# ---- Fit -------------------------------------------------------------------

location_formula <- stats::update(location_terms, y ~ .)

# start: optional starting values c(thresholds, location, scale), used by the
# request jackknife to refit from the full-sample estimates.
fit_uv <- function(agg, weights = agg$n, start = NULL) {
  agg$.w <- weights
  ordinal::clm(location_formula, scale = ~ unsafe_s, weights = .w, data = agg,
               link = "probit", start = start)
}

# Binary approve/reject probit: the equal-variance model. With one response
# threshold and a free mean in every design cell, the variance of
# unsafe-request evidence can't be separated from its mean, so binary data
# force equal variances.
fit_ev <- function(agg, weights = agg$n) {
  agg$.w <- weights
  ordinal::clm(location_formula, weights = .w, data = agg, link = "probit")
}

# ---- Score contributions ---------------------------------------------------
# clm models P(Y <= j) = Phi((theta_j - eta) / sigma), eta = x'beta,
# sigma = exp(z'zeta). For a review in category y, with
#   a = (theta_y - eta) / sigma,  b = (theta_{y-1} - eta) / sigma,
#   p = Phi(a) - Phi(b), the log-likelihood is log p, and
#   d/d theta_j = [phi(a) 1{j = y} - phi(b) 1{j = y - 1}] / (sigma p)
#   d/d beta    = -x [phi(a) - phi(b)] / (sigma p)
#   d/d zeta    =  z [-a phi(a) + b phi(b)] / p
# (theta_0 = -Inf, theta_J = Inf, where phi and a*phi vanish).

clm_scores <- function(fit, data) {
  alpha <- fit$alpha; beta <- fit$beta; zeta <- fit$zeta
  X <- stats::model.matrix(location_terms, data)[, names(beta), drop = FALSE]
  eta <- drop(X %*% beta)
  sig <- if (length(zeta)) exp(data$unsafe_s * zeta) else rep(1, nrow(data))
  y <- as.integer(data$y)
  th <- c(-Inf, alpha, Inf)
  a <- (th[y + 1] - eta) / sig
  b <- (th[y] - eta) / sig
  p <- stats::pnorm(a) - stats::pnorm(b)
  fa <- stats::dnorm(a); fb <- stats::dnorm(b)
  afa <- ifelse(is.finite(a), a * fa, 0); bfb <- ifelse(is.finite(b), b * fb, 0)

  J <- length(alpha)
  s_theta <- sapply(seq_len(J), function(j) (fa * (y == j) - fb * (y == j + 1)) / (sig * p))
  s_theta <- matrix(s_theta, ncol = J, dimnames = list(NULL, names(alpha)))
  s_beta  <- X * (-(fa - fb) / (sig * p))
  out <- cbind(s_theta, s_beta)
  if (length(zeta)) out <- cbind(out, unsafe_s = data$unsafe_s * (-afa + bfb) / p)
  attr(out, "prob") <- p  # P(y | x), for expected_info()
  out
}

# Cluster-robust (sandwich) covariance: bread = inverse information from clm,
# meat = sum over clusters of outer products of summed scores, with the
# G / (G - 1) small-sample factor. With two cluster variables (participants
# and requests) it uses the two-way formula V_1 + V_2 - V_12, where V_12
# clusters on the intersection (Cameron, Gelbach, & Miller, 2011).
vcov_cluster <- function(fit, agg, key, clusters) {
  S <- clm_scores(fit, agg)[key, , drop = FALSE]  # one row per review
  B <- stats::vcov(fit)
  one_way <- function(g) {
    U <- rowsum(S, g, reorder = FALSE)
    G <- nrow(U)
    B %*% (crossprod(U) * G / (G - 1)) %*% B
  }
  if (!is.list(clusters)) return(one_way(clusters))
  stopifnot(length(clusters) == 2)
  one_way(clusters[[1]]) + one_way(clusters[[2]]) -
    one_way(interaction(clusters[[1]], clusters[[2]], drop = TRUE))
}

# ---- Two-way errors with the request dimension jackknifed --------------------
# With 48 requests the request-clustered sandwich runs small (Cameron, Gelbach,
# & Miller, 2008). This version keeps the two-way formula but replaces its
# request-clustered part by the delete-one-request jackknife: refit leaving out
# each request in turn, and take (G - 1) / G times the sum of squared
# deviations of the leave-one-out contrasts from the full-sample ones (Cameron
# et al., 2008). The participant-clustered and observation-level
# parts stay analytic; with hundreds of participants their small-sample bias is
# small. Jackknife-type corrections tend to over-correct (Pustejovsky & Tipton,
# 2018), so this errs on the conservative side.

contrast_var <- function(fit, V, L) {
  idx <- match(colnames(L), names(stats::coef(fit)))
  diag(L %*% V[idx, idx] %*% t(L))
}

request_jackknife_var <- function(fit, a, item, L = contrast_L()) {
  items <- sort(unique(item)); G <- length(items)
  npat <- nrow(a$agg)
  # Pattern counts by request, so leaving one out only changes the weights.
  counts <- matrix(tabulate(a$key + npat * (match(item, items) - 1L), nbins = npat * G), npat, G)
  idx <- match(colnames(L), names(stats::coef(fit)))
  full <- drop(L %*% stats::coef(fit)[idx])
  loo <- vapply(seq_len(G), function(g) {
    w <- a$agg$n - counts[, g]
    keep <- w > 0
    f <- tryCatch(fit_uv(a$agg[keep, ], weights = w[keep], start = stats::coef(fit)),
                  error = function(e) NULL)
    if (is.null(f) || !all(names(stats::coef(fit)) %in% names(stats::coef(f)))) return(rep(NA_real_, nrow(L)))
    drop(L %*% stats::coef(f)[colnames(L)])
  }, numeric(nrow(L)))
  (G - 1) / G * rowSums((loo - full)^2)  # NA if any refit failed
}

# The request part is attached as attribute "request" (for two_way_df()).
vcov_two_way_jk <- function(fit, a, dat, L = contrast_L()) {
  V_pid <- vcov_cluster(fit, a$agg, a$key, dat$pid)
  V_obs <- vcov_cluster(fit, a$agg, a$key, interaction(dat$pid, dat$item, drop = TRUE))
  v_req <- request_jackknife_var(fit, a, dat$item, L)
  structure(contrast_var(fit, V_pid, L) + v_req - contrast_var(fit, V_obs, L), request = v_req)
}

# ---- Small-sample options for the request dimension ---------------------------
# Even with the jackknife, two-way intervals cover 91% to 93% instead of 95%:
# with 48 requests (12 unsafe) the request part of the variance is itself
# noisy, so (est - truth) / se has heavier tails than the normal. Two remedies,
# both off by default in analyze_study():
#   satterthwaite = TRUE: t reference distributions with degrees of freedom
#     estimated for each contrast (Pustejovsky & Tipton, 2018);
#   boot_draws > 0: a wild score bootstrap over requests for the analytic
#     two-way t statistic (Roodman, MacKinnon, Nielsen, & Webb, 2019).

# L's rows as columns over all clm parameters (thresholds, location, scale).
contrast_full <- function(fit, L) {
  C <- matrix(0, length(stats::coef(fit)), nrow(L), dimnames = list(names(stats::coef(fit)), rownames(L)))
  C[colnames(L), ] <- t(L)
  C
}

# Expected information in one review with covariate pattern x,
# I(x) = sum_y P(y | x) s(y, x) s(y, x)', for each pattern in agg (r x r x
# patterns), and each review's pattern.
expected_info <- function(fit, agg, key) {
  lab <- do.call(paste, agg[c("unsafe", "palter", "tactic", "trained", "screen", "block")])
  cid <- match(lab, unique(lab))
  cov <- agg[!duplicated(cid), ]
  K <- length(fit$alpha) + 1L
  n <- nrow(cov)
  long <- cov[rep(seq_len(n), times = K), ]
  long$y <- factor(rep(seq_len(K), each = n), levels = seq_len(K), ordered = TRUE)
  s <- clm_scores(fit, long)
  w <- attr(s, "prob")
  s[w == 0, ] <- 0  # an impossible rating carries no information
  x <- rep(seq_len(n), times = K)
  info <- vapply(seq_len(n), function(j) crossprod(s[x == j, , drop = FALSE] * w[x == j], s[x == j, , drop = FALSE]),
                 matrix(0, ncol(s), ncol(s)))
  list(info = info, pattern = cid[key])
}

# Expected information summed over the reviews in each cluster (r x r x clusters,
# clusters in sorted order).
cluster_info <- function(ei, cluster) {
  cl <- match(cluster, sort(unique(cluster)))
  r <- dim(ei$info)[1]; n <- dim(ei$info)[3]; G <- max(cl)
  counts <- matrix(tabulate(ei$pattern + n * (cl - 1L), nbins = n * G), n, G)
  array(matrix(ei$info, r * r, n) %*% counts, c(r, r, G))
}

# Satterthwaite degrees of freedom of a one-way cluster-robust variance for a
# contrast c: nu = (sum_g M_gg)^2 / sum_gk M_gk^2, with M_gk = p_g' Phi p_k
# computed under a working model (Pustejovsky & Tipton, 2018, eq. 11, p. 13).
# The working model: reviews are independent and the fitted model is right, so
# cluster g's summed scores U_g have covariance H_g, its expected information.
# Linearizing the fit, the scores at the estimate are U_g - H_g B sum_h U_h,
# B = (sum_g H_g)^-1, which gives
#   M_gk = 1{g = k} a_g' H_g a_g - a_g' H_g B H_k a_k
# with a_g = B c for the sandwich. The delete-one-cluster jackknife linearizes
# to the CR3 adjustment, which "closely approximates the jackknife" (p. 8):
# a_g = (I - B H_g)^-1 B c.
satterthwaite_df <- function(H, cvec, jackknife = FALSE) {
  B <- solve(rowSums(H, dims = 2))
  b <- drop(B %*% cvec)
  G <- dim(H)[3]
  A <- vapply(seq_len(G), function(g) if (jackknife) solve(diag(length(b)) - B %*% H[, , g], b) else b, b)
  U <- vapply(seq_len(G), function(g) drop(H[, , g] %*% A[, g]), b)
  M <- -crossprod(U, B %*% U)
  diag(M) <- diag(M) + colSums(A * U)
  sum(diag(M))^2 / sum(M^2)
}

# Degrees of freedom for a two-way variance V = V_pid + V_req - V_obs. With
# hundreds of participants, only the request part is noisy, so
# Var(V) ~ Var(V_req) = 2 V_req^2 / nu_req and nu = 2 V^2 / Var(V) = nu_req (V / V_req)^2.
two_way_df <- function(H, C, v_total, v_req, jackknife = FALSE) {
  nu_req <- apply(C, 2, function(cv) satterthwaite_df(H, cv, jackknife))
  ifelse(!is.na(v_total) & v_total > 0 & v_req > 0, nu_req * (v_total / v_req)^2, NA_real_)
}

# Wild score bootstrap over requests for the analytic two-way t statistic.
# Each draw multiplies all score contributions of request g by v_g = +/-1
# (Rademacher), the score bootstrap that Roodman et al. (2019, pp. 31-33)
# describe for ML estimators, clustering the draws in the dimension with fewer
# clusters (p. 20). The plain score bootstrap reuses the perturbed scores in
# the variance, so it "cannot account for the use of residuals instead of
# errors" (p. 32-33). Here each draw instead takes one Newton step,
# beta* - beta_hat = B S* with S* = sum_g v_g S_g, and recentres the scores there:
#   s*_k = v_g(k) s_k - h_k B S*,
# h_k being review k's expected information. Then
#   t* = c'B S* / sqrt(two-way variance of the s*_k),
# with the small-sample factors of vcov_cluster(). Everything is linear in v,
# so the draws cost only matrix products. Returns contrasts x draws; draws
# whose two-way variance is not positive are NA and are dropped, as boottest
# does (p. 19). v: optional requests x draws matrix of weights, for tests.
request_score_boot <- function(fit, a, dat, L = contrast_L(), draws = 999L, ei = NULL, v = NULL) {
  if (is.null(ei)) ei <- expected_info(fit, a$agg, a$key)
  B <- stats::vcov(fit)
  BC <- B %*% contrast_full(fit, L)
  items <- sort(unique(dat$item)); g <- match(dat$item, items); G <- length(items)
  pids <- sort(unique(dat$pid)); p <- match(dat$pid, pids); P <- length(pids)
  N <- nrow(dat); n_x <- dim(ei$info)[3]; r <- nrow(B)
  if (is.null(v)) v <- matrix(sample(c(-1, 1), G * draws, replace = TRUE), G, draws)

  s_pat <- clm_scores(fit, a$agg)
  S <- t(rowsum(s_pat[a$key, , drop = FALSE], g, reorder = TRUE))        # r x G request score sums
  d <- (s_pat %*% BC)[a$key, , drop = FALSE]                              # N x contrasts: c'B s_k
  n_req <- matrix(tabulate(ei$pattern + n_x * (g - 1L), nbins = n_x * G), n_x, G)
  n_pid <- matrix(tabulate(ei$pattern + n_x * (p - 1L), nbins = n_x * P), n_x, P)
  n_all <- rowSums(n_req)
  m <- function(k) k / (k - 1)  # as in vcov_cluster()

  out <- t(vapply(seq_len(ncol(BC)), function(j) {
    # e[, x] = B h_x B c: c'B s*_k = v_g(k) d_k - e[, x(k)]' S v.
    e <- B %*% matrix(apply(ei$info, 3, function(h) h %*% BC[, j]), r, n_x)
    D <- drop(rowsum(d[, j], g, reorder = TRUE))
    Phi_req <- diag(D, G) - crossprod(n_req, t(e)) %*% S                  # G x G
    D_pid <- matrix(0, P, G)
    cell <- rowsum(d[, j], p + P * (g - 1L), reorder = TRUE)
    D_pid[as.integer(rownames(cell))] <- cell
    Phi_pid <- D_pid - crossprod(n_pid, t(e)) %*% S                       # P x G
    # Observation level, using v_g^2 = 1:
    # sum_k (v_g d_k - e_k' S v)^2 = sum d_k^2 - 2 v'Q S v + v'S' Omega S v.
    Q <- rowsum(d[, j] * t(e)[ei$pattern, , drop = FALSE], g, reorder = TRUE)
    Omega <- e %*% (n_all * t(e))
    A <- -2 * Q %*% S + crossprod(S, Omega %*% S)
    num <- drop(crossprod(D, v))
    var_star <- m(P) * colSums((Phi_pid %*% v)^2) + m(G) * colSums((Phi_req %*% v)^2) -
      m(N) * (sum(d[, j]^2) + colSums(v * (A %*% v)))
    ifelse(var_star > 0, num / sqrt(pmax(var_star, 0)), NA_real_)
  }, numeric(ncol(v))))
  rownames(out) <- rownames(L)
  out
}

# ---- Contrasts, tests, verdict ---------------------------------------------

contrast_estimates <- function(fit, V, L = contrast_L()) {
  beta_idx <- match(colnames(L), names(stats::coef(fit)))
  stopifnot(!anyNA(beta_idx))
  est <- drop(L %*% stats::coef(fit)[beta_idx])
  v   <- diag(L %*% V[beta_idx, beta_idx] %*% t(L))
  se  <- sqrt(ifelse(v > 0, v, NA))  # the two-way formula can go negative; report NA
  data.frame(test = rownames(L), est = est, se = se,
             direction = attr(L, "direction")[rownames(L)], row.names = NULL)
}

# Directional z tests with Holm correction across all six; TOST equivalence
# against +/- sesoi for each contrast (used for T3). Both T3 contrasts must be
# equivalent to call the falsifier "not met" (intersection-union test, so no
# correction is needed across the pair).
# Reference distribution for (est - truth) / se: normal by default (df = Inf),
# t on df degrees of freedom (one per contrast), or bootstrap draws (tstar,
# one row per contrast, as from request_score_boot()).
run_tests <- function(ce, sesoi, alpha = 0.05, df = Inf, tstar = NULL) {
  upper <- function(x) {
    if (is.null(tstar)) return(stats::pt(x, df, lower.tail = FALSE))
    ifelse(is.na(x), NA_real_, rowMeans(tstar >= x, na.rm = TRUE))
  }
  lower <- function(x) {
    if (is.null(tstar)) return(stats::pt(x, df))
    ifelse(is.na(x), NA_real_, rowMeans(tstar <= x, na.rm = TRUE))
  }
  # Interval with tail probability tail on each side.
  bounds <- function(tail) {
    if (is.null(tstar)) {
      k <- stats::qt(1 - tail, df)
      return(list(ce$est - k * ce$se, ce$est + k * ce$se))
    }
    q <- function(p) apply(tstar, 1, stats::quantile, probs = p, na.rm = TRUE, names = FALSE)
    list(ce$est - q(1 - tail) * ce$se, ce$est - q(tail) * ce$se)
  }
  z <- ce$est / ce$se
  ce$p_dir  <- ifelse(ce$direction == "greater", upper(z), lower(z))
  # n = all six tests: p.adjust() otherwise drops a missing p-value (a negative
  # two-way variance) from the count and corrects over five.
  ce$p_holm <- stats::p.adjust(ce$p_dir, method = "holm", n = length(ce$p_dir))
  ce$sig    <- !is.na(ce$p_holm) & ce$p_holm < alpha  # a missing SE never counts as evidence
  ce$p_tost <- pmax(upper((ce$est + sesoi) / ce$se), lower((ce$est - sesoi) / ce$se))
  ce$equiv  <- !is.na(ce$p_tost) & ce$p_tost < alpha
  b90 <- bounds(alpha); b95 <- bounds(0.025)
  ce$lo90 <- b90[[1]]; ce$hi90 <- b90[[2]]
  ce$lo95 <- b95[[1]]; ce$hi95 <- b95[[2]]
  ce
}

verdict <- function(tt) {
  get <- function(test, col) tt[[col]][tt$test == test]
  cue_shift <- get("T2", "sig")
  d_gain    <- get("T3_plain", "sig") || get("T3_tactic", "sig")
  d_flat    <- get("T3_plain", "equiv") && get("T3_tactic", "equiv")
  if (cue_shift && d_flat && !d_gain) return("cue")
  if (!cue_shift && d_gain) return("vigilance")
  if (cue_shift && d_gain) return("both")
  if (!cue_shift && d_flat && !d_gain) return("null")
  "inconclusive"
}

# ---- Pseudo-true values and the SESOI on the estimator's scale --------------

pseudo_true <- function(pars, L = contrast_L()) {
  pop <- population_grid(pars)
  pop$y <- factor(pop$rating, levels = 1:6, ordered = TRUE)
  pop$trial_c <- pop$block_c
  pop$unsafe_s <- pop$unsafe
  uv <- fit_uv(pop, weights = pop$w * 1e5)

  pop$y2 <- as.integer(pop$rating >= 4) + 1L  # 1 = approve, 2 = reject
  bin <- stats::aggregate(w ~ unsafe + palter + tactic + trained + screen + block_c + y2, data = pop, FUN = sum)
  bin$y <- factor(bin$y2, levels = 1:2, ordered = TRUE)
  bin$trial_c <- bin$block_c
  bin$screen <- factor(bin$screen, levels = screens)
  ev <- fit_ev(bin, weights = bin$w * 1e5)

  pick <- function(fit) {
    b <- stats::coef(fit)[colnames(L)]
    stats::setNames(drop(L %*% b), rownames(L))
  }
  list(uv = pick(uv), ev = pick(ev), uv_fit = uv, ev_fit = ev,
       sigma_ratio = unname(exp(uv$zeta)))
}

# A participant-level gain of delta in d maps to this gain on the population-
# averaged scale a given model reports. "uv" is the fast estimator; "ev" is the
# binary equal-variance comparison, which has its own, smaller image of the
# same gain. The mapping uses the generating variance components; in a real
# study they would come from the pilot or from the primary model's estimates.
sesoi_on_estimator_scale <- function(pars, delta = 0.20, model = c("uv", "ev")) {
  model <- match.arg(model)
  pt <- pseudo_true(utils::modifyList(pars, list(kappa = 0, gamma = delta)))
  unname(mean(pt[[model]][c("T3_plain", "T3_tactic")]))
}

# ---- Applied outcome --------------------------------------------------------
# Share of unsafe requests approved when the approval threshold is placed so
# that a preregistered share of safe requests is approved.
unsafe_approval_at <- function(d, sigma_ratio, safe_approval = 0.90) {
  stats::pnorm((stats::qnorm(safe_approval) - d) / sigma_ratio)
}

# ---- One replicate ----------------------------------------------------------

# Tests use participant-clustered errors by default. Verdicts under two-way
# (participant and request) clustering, analytic and with the request
# jackknife, are returned too. jackknife = FALSE skips the 48 refits.
# Small-sample options for the request dimension (off by default):
#   satterthwaite = TRUE adds the analytic and jackknifed two-way tests with t
#     reference distributions on estimated degrees of freedom ("two_way_sw",
#     "two_way_jk_sw");
#   boot_draws > 0 adds the analytic two-way tests with the request score
#     bootstrap as reference distribution ("two_way_boot").
# clusters picks which set of tests becomes $tests and $verdict.
analyze_study <- function(dat, sesoi, sesoi_ev = sesoi, L = contrast_L(), with_ev = TRUE,
                          clusters = c("pid", "two_way", "two_way_jk", "two_way_sw", "two_way_jk_sw",
                                       "two_way_boot"),
                          jackknife = TRUE, satterthwaite = FALSE, boot_draws = 0L) {
  clusters <- match.arg(clusters)
  stopifnot(jackknife || !clusters %in% c("two_way_jk", "two_way_jk_sw"),
            satterthwaite || !clusters %in% c("two_way_sw", "two_way_jk_sw"),
            boot_draws > 0 || clusters != "two_way_boot")
  cl_pid <- dat$pid
  cl_two <- list(dat$pid, dat$item)

  a <- aggregate_reviews(dat, "rating")
  uv <- fit_uv(a$agg)
  t_pid <- run_tests(contrast_estimates(uv, vcov_cluster(uv, a$agg, a$key, cl_pid), L), sesoi)
  t_two <- run_tests(contrast_estimates(uv, vcov_cluster(uv, a$agg, a$key, cl_two), L), sesoi)
  t_jk <- NULL
  if (jackknife) {
    v <- vcov_two_way_jk(uv, a, dat, L)
    t_jk <- run_tests(transform(t_pid[, c("test", "est", "direction")],
                                se = sqrt(ifelse(!is.na(v) & v > 0, v, NA))), sesoi)
  }
  t_two_sw <- t_jk_sw <- t_boot <- NULL
  df_two <- df_jk <- rep(NA_real_, nrow(L))
  if (satterthwaite || boot_draws > 0) ei <- expected_info(uv, a$agg, a$key)
  if (satterthwaite) {
    H <- cluster_info(ei, dat$item)
    C <- contrast_full(uv, L)
    v_req <- contrast_var(uv, vcov_cluster(uv, a$agg, a$key, dat$item), L)
    df_two <- two_way_df(H, C, t_two$se^2, v_req)
    t_two_sw <- run_tests(t_two[, c("test", "est", "se", "direction")], sesoi, df = df_two)
    if (jackknife) {
      df_jk <- two_way_df(H, C, t_jk$se^2, attr(v, "request"), jackknife = TRUE)
      t_jk_sw <- run_tests(t_jk[, c("test", "est", "se", "direction")], sesoi, df = df_jk)
    }
  }
  if (boot_draws > 0) {
    tstar <- request_score_boot(uv, a, dat, L, draws = boot_draws, ei = ei)
    t_boot <- run_tests(t_two[, c("test", "est", "se", "direction")], sesoi, tstar = tstar)
  }
  tt <- switch(clusters, pid = t_pid, two_way = t_two, two_way_jk = t_jk, two_way_sw = t_two_sw,
               two_way_jk_sw = t_jk_sw, two_way_boot = t_boot)
  tt$se_pid     <- t_pid$se
  tt$se_two_way <- t_two$se
  tt$se_two_way_jk <- if (jackknife) t_jk$se else NA_real_
  tt$se_model   <- contrast_estimates(uv, stats::vcov(uv), L)$se
  if (satterthwaite) {
    tt$df_two_way <- df_two
    tt$df_two_way_jk <- df_jk
  }
  tt$model <- "unequal variance"
  verdict_or_na <- function(t) if (is.null(t)) NA_character_ else verdict(t)
  out <- list(tests = tt, verdict = verdict(tt), verdict_pid = verdict(t_pid),
              verdict_two_way = verdict(t_two),
              verdict_two_way_jk = if (jackknife) verdict(t_jk) else NA_character_,
              tests_two_way_jk = t_jk, sigma_ratio = unname(exp(uv$zeta)))
  if (satterthwaite) {
    out$verdict_two_way_sw <- verdict(t_two_sw)
    out$verdict_two_way_jk_sw <- verdict_or_na(t_jk_sw)
    out$tests_two_way_sw <- t_two_sw
    out$tests_two_way_jk_sw <- t_jk_sw
  }
  if (boot_draws > 0) {
    out$verdict_two_way_boot <- verdict(t_boot)
    out$tests_two_way_boot <- t_boot
  }

  if (with_ev) {
    b <- aggregate_reviews(dat, "reject")
    ev <- fit_ev(b$agg)
    cl <- if (clusters == "pid") cl_pid else cl_two
    te <- run_tests(contrast_estimates(ev, vcov_cluster(ev, b$agg, b$key, cl), L), sesoi_ev)
    te$model <- "equal variance (binary)"
    out$tests_ev <- te
    out$verdict_ev <- verdict(te)
  }
  out
}
