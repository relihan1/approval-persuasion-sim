# simulate.R
# Simulate one study, and compute the exact population the estimator targets.

# ---- Request bank ----------------------------------------------------------
# Each request has a plain and a tactic-laden version. Two counterbalancing
# lists show each participant every request once, in one version; list B
# flips list A. Tactic versions are balanced within safe, unsafe-palter and
# unsafe-non-palter requests, so every participant sees half of each stratum
# with tactics.

make_items <- function(pars) {
  K <- pars$n_safe + pars$n_unsafe
  stopifnot(K %% 4 == 0)  # trial blocks are quarters of the session
  items <- data.frame(item = seq_len(K),
                      unsafe = rep(c(0L, 1L), c(pars$n_safe, pars$n_unsafe)))
  n_palter <- round(pars$palter_share * pars$n_unsafe)
  items$palter <- 0L
  items$palter[items$unsafe == 1L][seq_len(n_palter)] <- 1L
  stratum <- interaction(items$unsafe, items$palter, drop = TRUE)
  items$tactic_A <- as.integer(stats::ave(items$item, stratum, FUN = seq_along) %% 2 == 1)
  items
}

# ---- One simulated study ---------------------------------------------------

simulate_study <- function(pars = dgm_defaults(), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  items <- make_items(pars)
  K <- nrow(items)

  # Participants: 6 between-subject cells, lists alternate within each cell.
  cells <- expand.grid(trained = 0:1, screen = screens, stringsAsFactors = FALSE)
  P <- cells[rep(seq_len(nrow(cells)), each = pars$n_per_cell), ]
  P$pid  <- seq_len(nrow(P))
  P$list <- ifelse(stats::ave(P$pid, P$trained, P$screen, FUN = seq_along) %% 2 == 1, "A", "B")
  P$u_c  <- stats::rnorm(nrow(P), 0, pars$tau_c)
  P$u_d  <- stats::rnorm(nrow(P), 0, pars$tau_d)
  P$u_t  <- stats::rnorm(nrow(P), 0, pars$tau_t)

  # Requests: a fresh draw each study, so power generalizes over requests.
  items$v_c <- stats::rnorm(K, 0, pars$omega_c)
  items$v_d <- stats::rnorm(K, 0, pars$omega_d)
  items$v_t <- stats::rnorm(K, 0, pars$omega_t)
  items$v_k <- stats::rnorm(K, 0, pars$omega_k)

  # Every participant reviews every request once, in a random order.
  D <- cbind(P[rep(seq_len(nrow(P)), each = K), ], items[rep(seq_len(K), times = nrow(P)), ])
  D$trial  <- as.vector(replicate(nrow(P), sample.int(K)))
  D$tactic <- ifelse(D$list == "A", D$tactic_A, 1L - D$tactic_A)

  mu <- latent_mean(D$unsafe, D$palter, D$tactic, D$trained, D$screen, D$trial, pars,
                    u_c = D$u_c, u_d = D$u_d, u_t = D$u_t, v_c = D$v_c, v_d = D$v_d,
                    v_t = D$v_t, v_k = D$v_k)
  x <- mu + latent_sd(D$unsafe, pars) * stats::rnorm(nrow(D))

  out <- data.frame(
    pid = D$pid, item = D$item, trained = D$trained,
    screen = factor(D$screen, levels = screens), list = D$list,
    unsafe = D$unsafe, palter = D$palter, tactic = D$tactic, trial = D$trial,
    rating = findInterval(x, pars$cuts) + 1L
  )
  out$reject  <- as.integer(out$rating >= 4L)
  out$trial_c <- (out$trial - (K + 1) / 2) / K       # row-level trial position
  out$block   <- ceiling(out$trial / (K / 4))          # quarter of the session
  out$block_c <- (out$block - 2.5) * (K / 4) / K       # block midpoint, same units
  rownames(out) <- NULL
  out
}

# ---- Exact population ------------------------------------------------------
# Integrates the participant and request effects analytically (they are normal
# and enter the latent mean additively), so each design cell's rating
# distribution is a known mixture over palter status and trial position.
# Fitting the estimator to these expected frequencies gives its pseudo-true
# values: the numbers it converges to with unlimited participants and requests.

population_grid <- function(pars) {
  K <- pars$n_safe + pars$n_unsafe
  n_palter <- round(pars$palter_share * pars$n_unsafe)
  kinds <- data.frame(unsafe = c(0L, 1L, 1L), palter = c(0L, 1L, 0L),
                      share = c(pars$n_safe, n_palter, pars$n_unsafe - n_palter) / K)
  g <- merge(kinds, expand.grid(tactic = 0:1, trained = 0:1, screen = screens,
                                trial = seq_len(K), stringsAsFactors = FALSE))

  m <- latent_mean(g$unsafe, g$palter, g$tactic, g$trained, g$screen, g$trial, pars)
  s <- marginal_sd(g$unsafe, g$tactic, g$trained, g$screen, pars)

  cuts <- c(-Inf, pars$cuts, Inf)
  probs <- sapply(1:6, function(k) stats::pnorm((cuts[k + 1] - m) / s) - stats::pnorm((cuts[k] - m) / s))

  long <- data.frame(g[rep(seq_len(nrow(g)), 6), c("unsafe", "palter", "tactic", "trained", "screen", "trial")],
                     rating = rep(1:6, each = nrow(g)),
                     w = as.vector(probs) * rep(g$share, 6) / 2 / 6 / K)  # tactic, cell, trial weights
  long$block <- ceiling(long$trial / (K / 4))
  pop <- stats::aggregate(w ~ unsafe + tactic + trained + screen + block + rating, data = long, FUN = sum)
  pop$screen  <- factor(pop$screen, levels = screens)
  pop$block_c <- (pop$block - 2.5) * (K / 4) / K
  pop
}
