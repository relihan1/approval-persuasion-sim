# design.R
# Study design, data-generating parameters, rival hypotheses, and the
# preregistered contrasts.
#
# Scale: everything is on the latent probit scale used by signal detection
# theory. The noise SD for safe requests is 1. Higher latent values mean a
# stronger lean toward "reject". A participant rates each request on a
# 6-point scale: 1 = definitely approve ... 6 = definitely reject, and
# ratings of 4 or higher count as a rejection.

# ---- Data-generating parameters ------------------------------------------

dgm_defaults <- function() {
  list(
    # Design
    n_per_cell   = 150,   # participants in each training x screen cell (6 cells)
    n_safe       = 36,    # safe requests each participant reviews
    n_unsafe     = 12,    # unsafe requests each participant reviews (25%)

    # Signal detection structure
    d0           = 1.5,   # discrimination for plain requests on the argument screen
    sigma_unsafe = 1.25,  # SD of unsafe-request evidence (zROC slope = 1/1.25 = 0.8)
    cuts         = c(-0.6, 0.2, 1.0, 1.6, 2.2),  # rating thresholds

    # Crossed random effects
    tau_c   = 0.5,  # participants differ in overall tendency to reject
    tau_d   = 0.3,  # participants differ in discrimination
    tau_t   = 0.0,  # participants differ in susceptibility to tactics (PPI script sets 0.3)
    omega_c = 0.3,  # requests differ in how rejectable they look
    omega_d = 0.3,  # unsafe requests differ in how detectable they are
    omega_t = 0.0,  # requests differ in how much their tactics pull (script 06 sets 0.3)
    omega_k = 0.0,  # requests differ in how strongly trained reviewers react to their tactics (06: 0.2)

    # What the agent's message and the screen do
    pull         = 0.4,  # tactics move any request toward approval
    pull_extra   = 0.2,  # ...and pull harder on unsafe ones, which give the reviewer doubts to answer
    strip        = 0.5,  # share of tactic cues a neutral paraphrase removes
    evidence     = 0.5,  # independent evidence raises discrimination
    palter       = 0.6,  # true-but-misleading justifications lower d unless evidence is shown
    palter_share = 0.5,  # share of unsafe requests whose justification palters
    drift        = 0.006,  # per-trial drift toward approval (approval fatigue)

    # The two accounts of what training does
    kappa = 0.5,  # cue account: training moves tactic-laden requests toward reject
    gamma = 0.0   # vigilance account: training raises discrimination
  )
}

# The four states of the world the adjudication rule must tell apart.
# gamma = 0.3 matches the discrimination gains Simchon et al. (2026) report for
# prebunking (median 0.31 between participants, 0.29 within); kappa = 0.5 is a
# moderate criterion shift, an assumption to vary.
hypotheses <- function() {
  list(
    cue       = list(kappa = 0.5, gamma = 0.0),
    vigilance = list(kappa = 0.0, gamma = 0.3),
    both      = list(kappa = 0.5, gamma = 0.3),
    null      = list(kappa = 0.0, gamma = 0.0)
  )
}

with_hypothesis <- function(pars, h) {
  hyp <- if (is.character(h)) hypotheses()[[h]] else h
  utils::modifyList(pars, hyp)
}

screens <- c("argument", "paraphrase", "evidence")

# ---- Latent mean of one review --------------------------------------------
# Shared by the simulator and by the exact population projection, so the two
# can never disagree about the data-generating process.

latent_mean <- function(unsafe, palter, tactic, trained, screen, trial, pars,
                        u_c = 0, u_d = 0, u_t = 0, v_c = 0, v_d = 0, v_t = 0, v_k = 0) {
  ev <- as.numeric(screen == "evidence")
  pa <- as.numeric(screen == "paraphrase")
  cue_visible <- 1 - pars$strip * pa  # paraphrase hides part of the tactic cues

  signal <- unsafe * (pars$d0 + u_d + v_d +
                        pars$evidence * ev -
                        pars$palter * palter * (1 - ev) +
                        pars$gamma * trained)
  tactics <- -tactic * (pars$pull + pars$pull_extra * unsafe + u_t + v_t) * cue_visible +
    tactic * trained * (pars$kappa + v_k) * cue_visible

  signal + tactics + u_c + v_c - pars$drift * (trial - 1)
}

latent_sd <- function(unsafe, pars) ifelse(unsafe == 1, pars$sigma_unsafe, 1)

# SD of the latent variable in a design cell once participant and request
# effects are integrated out (they are independent normals entering additively).
marginal_sd <- function(unsafe, tactic, trained, screen, pars) {
  cue_visible <- 1 - pars$strip * as.numeric(screen == "paraphrase")
  sqrt(latent_sd(unsafe, pars)^2 + pars$tau_c^2 + pars$omega_c^2 +
         unsafe * (pars$tau_d^2 + pars$omega_d^2) +
         tactic * cue_visible^2 * (pars$tau_t^2 + pars$omega_t^2 + trained * pars$omega_k^2))
}

# ---- Cells and preregistered contrasts --------------------------------------

# 24 cells: request type (safe/unsafe) x justification (plain/tactic) x
# training (control/trained) x screen (argument/paraphrase/evidence).
# Unsafe cells sit at the bank's mix of paltering and non-paltering requests
# (palter = palter_share), so each cell's latent mean averages over that mix.
cell_grid <- function(palter_share = 0.5) {
  g <- expand.grid(unsafe = 0:1, tactic = 0:1, trained = 0:1, screen = screens,
                   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  g$screen <- factor(g$screen, levels = screens)
  g$palter <- g$unsafe * palter_share
  g
}

# The location model. Paltering is a known property of each unsafe request,
# and the evidence screen removes its cost (latent_mean), so it enters with its
# screen interaction. Leaving it out makes paltering a fixed difference between
# requests that request-clustered errors count as noise.
location_terms <- ~ unsafe * tactic * trained * screen + palter + palter:screen + trial_c

# Rows: cells. Columns: location coefficients (same names as the fitted models).
cell_X <- function(grid = cell_grid()) {
  grid$trial_c <- 0
  X <- stats::model.matrix(location_terms, grid)
  X[, colnames(X) != "(Intercept)", drop = FALSE]
}

# Each contrast is a weight vector over the 24 cells; L = W %*% cell_X().
# d in a cell = latent mean (unsafe) - latent mean (safe), in safe-noise SD units.
planned_contrasts <- function(grid = cell_grid()) {
  is <- function(u, tac, tr, sc) {
    as.numeric(grid$unsafe == u & grid$tactic == tac &
                 grid$trained == tr & grid$screen == sc)
  }
  d_w <- function(tac, tr, sc) is(1, tac, tr, sc) - is(0, tac, tr, sc)
  avg <- function(x) Reduce(`+`, x) / length(x)

  # T1 threat check: tactics lower rejection of unsafe requests
  # (control group, argument screen).
  T1 <- is(1, 1, 0, "argument") - is(1, 0, 0, "argument")

  # T2 cue effect: on SAFE requests, training raises rejection more for
  # tactic-laden than for plain justifications (averaged over screens).
  T2 <- avg(lapply(screens, function(s) {
    (is(0, 1, 1, s) - is(0, 1, 0, s)) - (is(0, 0, 1, s) - is(0, 0, 0, s))
  }))

  # T3 falsifier: training changes d among plain, and among tactic-laden,
  # requests (averaged over screens).
  T3_plain  <- avg(lapply(screens, function(s) d_w(0, 1, s) - d_w(0, 0, s)))
  T3_tactic <- avg(lapply(screens, function(s) d_w(1, 1, s) - d_w(1, 0, s)))

  # T4, T5: evidence improves d relative to the argument and the paraphrase.
  d_screen <- function(s) avg(list(d_w(0, 0, s), d_w(0, 1, s), d_w(1, 0, s), d_w(1, 1, s)))
  T4 <- d_screen("evidence") - d_screen("argument")
  T5 <- d_screen("evidence") - d_screen("paraphrase")

  W <- rbind(T1 = T1, T2 = T2, T3_plain = T3_plain, T3_tactic = T3_tactic,
             T4 = T4, T5 = T5)
  attr(W, "direction") <- c(T1 = "less", T2 = "greater", T3_plain = "greater",
                            T3_tactic = "greater", T4 = "greater", T5 = "greater")
  W
}

contrast_L <- function(W = planned_contrasts()) {
  L <- W %*% cell_X()
  attr(L, "direction") <- attr(W, "direction")
  L
}

# For a model without the palter terms (the brms primary model): drop those
# columns of L, keeping the test directions, which column subsetting would lose.
drop_absent_terms <- function(L, coef_names) {
  have <- colnames(L) %in% coef_names
  stopifnot(all(grepl("palter", colnames(L)[!have])))
  out <- L[, have, drop = FALSE]
  attr(out, "direction") <- attr(L, "direction")
  out
}

# The verdict each true state should produce.
correct_verdict <- c(cue = "cue", vigilance = "vigilance", both = "both", null = "null")
