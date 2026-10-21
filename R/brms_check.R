# brms_check.R
# Convergence checks for the brms primary model (script 05).
#
# Script 05 used to check R-hat and ESS for the location coefficients only. A
# fit can pass that check while the request-level SDs, their correlations, or
# the thresholds have not mixed, and the verdict rests on six linear
# combinations of the coefficients rather than on any one of them. These
# helpers summarise each parameter class, the six contrasts, and the sampler's
# own diagnostics. The cut-offs follow Vehtari et al. (2021): R-hat below
# 1.01, and bulk and tail ESS of at least 100 per chain. E-BFMI below 0.3 is
# the level at which Stan warns.

# Parameter class from a brms variable name; NA for variables not checked
# (lp__, lprior, the z_ and Cor_ internals, and Intercept[k], brms's thresholds
# before the predictors are un-centred, which b_Intercept[k] already covers).
brms_par_class <- function(v) {
  patterns <- c(disc = "^b_disc_", threshold = "^b_Intercept\\[", b = "^b_",
                sd_pid = "^sd_pid__", sd_item = "^sd_item__",
                cor_pid = "^cor_pid__", cor_item = "^cor_item__",
                r_pid = "^r_pid\\[", r_item = "^r_item\\[")
  out <- rep(NA_character_, length(v))
  for (k in names(patterns)) out[is.na(out) & grepl(patterns[[k]], v)] <- k
  out
}

# Draws of the contrasts L %*% b with the chains kept, so that R-hat and ESS
# can be computed for the quantities the verdict uses.
contrast_draws <- function(draws, L) {
  a <- posterior::as_draws_array(draws)
  a <- posterior::subset_draws(a, variable = paste0("b_", colnames(L)))
  dims <- dim(a)
  est <- array(NA_real_, c(dims[1], dims[2], nrow(L)),
               dimnames = list(NULL, NULL, rownames(L)))
  for (ch in seq_len(dims[2])) est[, ch, ] <- unclass(a)[, ch, ] %*% t(L)
  posterior::as_draws_array(est)
}

# R-hat and ESS by parameter class, with the worst variable in each.
convergence_by_class <- function(draws, contrasts = NULL) {
  v <- posterior::variables(draws)
  cls <- brms_par_class(v)
  s <- posterior::summarise_draws(posterior::subset_draws(draws, variable = v[!is.na(cls)]),
                                  "rhat", "ess_bulk", "ess_tail")
  s$class <- cls[!is.na(cls)]
  if (!is.null(contrasts)) {
    sc <- posterior::summarise_draws(contrasts, "rhat", "ess_bulk", "ess_tail")
    sc$class <- "contrast"
    s <- rbind(s, sc)
  }
  rows <- lapply(split(s, factor(s$class, unique(s$class))), function(x) data.frame(
    class = x$class[1], n = nrow(x),
    max_rhat = max(x$rhat), worst_rhat = x$variable[which.max(x$rhat)],
    min_ess_bulk = min(x$ess_bulk), worst_ess_bulk = x$variable[which.min(x$ess_bulk)],
    min_ess_tail = min(x$ess_tail), worst_ess_tail = x$variable[which.min(x$ess_tail)]))
  do.call(rbind, c(rows, make.row.names = FALSE))
}

# Divergences, iterations that hit the tree-depth limit, and E-BFMI per chain,
# from the post-warmup NUTS parameters (brms::nuts_params()).
sampler_summary <- function(np, max_treedepth = 10) {
  by_chain <- function(par, f) {
    x <- np[np$Parameter == par, ]
    as.vector(tapply(x$Value, x$Chain, f))
  }
  ebfmi <- function(e) sum(diff(e)^2) / sum((e - mean(e))^2)
  data.frame(chain = sort(unique(np$Chain)),
             divergent = by_chain("divergent__", sum),
             treedepth_hits = by_chain("treedepth__", function(x) sum(x >= max_treedepth)),
             ebfmi = by_chain("energy__", ebfmi))
}

# The checks a fit has to pass before its verdict is reported:
# - R-hat below 1.01 in every class except the individual participant and
#   request effects (r_pid, r_item), which are summarised but not checked;
# - bulk and tail ESS of at least 100 per chain for the six contrasts, the SD
#   of unsafe-request evidence (disc), and the random-effect SDs, which are
#   what the verdict and the SESOI conversion use. The floor doesn't apply to
#   the thresholds, the individual location coefficients (class b), or the
#   random-effect correlations; their ESS is reported but not required. The
#   contrasts are checked directly, so a coefficient's own ESS matters only
#   through them. The slowest coefficients, those that differ only between
#   participants (training, screen) or only between requests (unsafe, palter),
#   trade off against the means of the random intercepts and have no weight in
#   any contrast (see test-brms-check.R);
# - no divergent transitions, and E-BFMI of at least 0.3 in every chain.
ess_classes <- c("contrast", "disc", "sd_pid", "sd_item")

convergence_ok <- function(by_class, sampler, n_chains) {
  # A class missing from the table would otherwise pass its ESS floor unchecked.
  missing <- setdiff(ess_classes, by_class$class)
  if (length(missing)) stop("No draws for ", paste(missing, collapse = ", "), "; the ESS floor can't be checked.")
  checked <- by_class[!by_class$class %in% c("r_pid", "r_item"), ]
  ess_checked <- checked[checked$class %in% ess_classes, ]
  # isTRUE: an NA diagnostic (e.g. a constant draw) fails the check instead of returning NA.
  isTRUE(all(checked$max_rhat < 1.01) &&
    all(ess_checked$min_ess_bulk >= 100 * n_chains) &&
    all(ess_checked$min_ess_tail >= 100 * n_chains) &&
    sum(sampler$divergent) == 0 && all(sampler$ebfmi >= 0.3))
}
