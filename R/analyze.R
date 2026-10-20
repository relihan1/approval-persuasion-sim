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
# frequency-weighted likelihood of the distinct patterns. Fitting 576
# weighted rows instead of 43,200 reviews makes each fit take milliseconds.
# Trial position enters as the quarter of the session (block_c).

pattern_code <- function(d, outcome = "rating", n_out = 6L) {
  sc <- as.integer(d$screen) - 1L
  ((((d$unsafe * 2L + d$tactic) * 2L + d$trained) * 3L + sc) * 4L + (d$block - 1L)) * n_out +
    (d[[outcome]] - 1L) + 1L
}

pattern_table <- function(n_out = 6L) {
  g <- expand.grid(y = seq_len(n_out), block = 1:4, screen = screens, trained = 0:1,
                   tactic = 0:1, unsafe = 0:1, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  g <- g[, rev(names(g))]  # unsafe, tactic, trained, screen, block, y (matches pattern_code order)
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

fit_uv <- function(agg, weights = agg$n) {
  agg$.w <- weights
  ordinal::clm(y ~ unsafe * tactic * trained * screen + trial_c, scale = ~ unsafe_s,
               weights = .w, data = agg, link = "probit")
}

# Binary approve/reject probit: the equal-variance model. With one response
# threshold, the variance of unsafe-request evidence can't be separated from
# its mean, so binary data force equal variances.
fit_ev <- function(agg, weights = agg$n) {
  agg$.w <- weights
  ordinal::clm(y ~ unsafe * tactic * trained * screen + trial_c,
               weights = .w, data = agg, link = "probit")
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
run_tests <- function(ce, sesoi, alpha = 0.05) {
  z <- ce$est / ce$se
  ce$p_dir  <- ifelse(ce$direction == "greater", stats::pnorm(z, lower.tail = FALSE), stats::pnorm(z))
  # n = all six tests: p.adjust() otherwise drops a missing p-value (a negative
  # two-way variance) from the count and corrects over five.
  ce$p_holm <- stats::p.adjust(ce$p_dir, method = "holm", n = length(ce$p_dir))
  ce$sig    <- !is.na(ce$p_holm) & ce$p_holm < alpha  # a missing SE never counts as evidence
  ce$p_tost <- pmax(stats::pnorm((ce$est + sesoi) / ce$se, lower.tail = FALSE),
                    stats::pnorm((ce$est - sesoi) / ce$se))
  ce$equiv  <- !is.na(ce$p_tost) & ce$p_tost < alpha
  ce$lo90 <- ce$est - stats::qnorm(1 - alpha) * ce$se
  ce$hi90 <- ce$est + stats::qnorm(1 - alpha) * ce$se
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
  bin <- stats::aggregate(w ~ unsafe + tactic + trained + screen + block_c + y2, data = pop, FUN = sum)
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

# Tests use participant-clustered errors by default. The verdict under two-way
# (participant and request) clustering is returned too, for script 06.
analyze_study <- function(dat, sesoi, sesoi_ev = sesoi, L = contrast_L(), with_ev = TRUE,
                          clusters = c("pid", "two_way")) {
  clusters <- match.arg(clusters)
  cl_pid <- dat$pid
  cl_two <- list(dat$pid, dat$item)

  a <- aggregate_reviews(dat, "rating")
  uv <- fit_uv(a$agg)
  t_pid <- run_tests(contrast_estimates(uv, vcov_cluster(uv, a$agg, a$key, cl_pid), L), sesoi)
  t_two <- run_tests(contrast_estimates(uv, vcov_cluster(uv, a$agg, a$key, cl_two), L), sesoi)
  tt <- if (clusters == "pid") t_pid else t_two
  tt$se_pid     <- t_pid$se
  tt$se_two_way <- t_two$se
  tt$se_model   <- contrast_estimates(uv, stats::vcov(uv), L)$se
  tt$model <- "unequal variance"
  out <- list(tests = tt, verdict = verdict(tt), verdict_pid = verdict(t_pid),
              verdict_two_way = verdict(t_two), sigma_ratio = unname(exp(uv$zeta)))

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
