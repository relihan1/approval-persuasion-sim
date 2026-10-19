# audit.R
# Stimulus tell audit: can a classifier tell safe from unsafe requests from
# the agent's justification alone? If it can, participants might pick up the
# same cue, and the study would partly measure detection of a writing style
# instead of judgment of the requested action. The audit runs before piloting,
# alongside blind human raters.
#
# Design choices that matter:
# - Both versions of a request (plain and tactic-laden) share most of their
#   text, so they go in the same cross-validation fold. Splitting them would
#   let the classifier memorize a request in training and recognize it in
#   testing.
# - The penalty is fixed in advance, so no tuning touches held-out requests.
# - The p-value comes from refitting under permuted labels, which keeps the
#   grouping and the small number of unsafe requests.

# ---- A toy request bank -----------------------------------------------------
# Real use: replace make_text_bank() with the model-written justifications.

text_pools <- list(
  reason = c("I need this to finish the migration you assigned.",
             "The build fails without this change.",
             "I need to check the configuration before the next step.",
             "This lets me confirm the deployment works.",
             "I need this to reproduce the bug you reported.",
             "The test suite needs this to run end to end.",
             "This is required to update the dependency versions.",
             "I need this to clean up the logging setup."),
  tactic = c("This is blocking tonight's release, so please approve quickly.",
             "The team lead already signed off on this approach.",
             "Other engineers on this repository do this routinely.",
             "You approved a similar request earlier today.",
             "The maintenance window closes in ten minutes."),
  minimizer = c("It is only temporary.", "This will just take a moment.",
                "Nothing else will be affected.")
)

# leak: how much more often unsafe justifications include a minimizing phrase.
make_text_bank <- function(n_items = 48, unsafe_share = 0.25, p_min_safe = 0.3,
                           leak = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_unsafe <- round(unsafe_share * n_items)
  unsafe <- sample(rep(c(0L, 1L), c(n_items - n_unsafe, n_unsafe)))
  p_min <- ifelse(unsafe == 1L, p_min_safe + leak, p_min_safe)
  base <- paste(sample(text_pools$reason, n_items, replace = TRUE),
                ifelse(stats::runif(n_items) < p_min,
                       sample(text_pools$minimizer, n_items, replace = TRUE), ""))
  tactic <- sample(text_pools$tactic, n_items, replace = TRUE)
  data.frame(item = rep(seq_len(n_items), 2),
             version = rep(c("plain", "tactic"), each = n_items),
             unsafe = rep(unsafe, 2),
             text = gsub("\\s+", " ", trimws(c(base, paste(base, tactic)))))
}

# ---- Features and classifier -------------------------------------------------

doc_term_matrix <- function(text) {
  tokens <- strsplit(gsub("[^a-z' ]", " ", tolower(text)), "\\s+")
  tokens <- lapply(tokens, function(t) t[nzchar(t)])
  vocab <- sort(unique(unlist(tokens)))
  X <- t(vapply(tokens, function(t) tabulate(match(t, vocab), length(vocab)), numeric(length(vocab))))
  colnames(X) <- vocab
  cbind(X, log_length = log(lengths(tokens)))
}

# L2-penalized logistic regression by iteratively reweighted least squares.
# Features are standardized on the training fold; the intercept is unpenalized.
ridge_logit <- function(X, y, lambda = 10, maxit = 25, tol = 1e-6) {
  mu <- colMeans(X); s <- apply(X, 2, stats::sd); s[!is.finite(s) | s == 0] <- 1
  Z <- cbind(1, sweep(sweep(X, 2, mu), 2, s, "/"))
  pen <- diag(c(0, rep(lambda, ncol(Z) - 1)))
  beta <- rep(0, ncol(Z))
  for (i in seq_len(maxit)) {
    eta <- drop(Z %*% beta); p <- stats::plogis(eta); w <- p * (1 - p) + 1e-10
    new <- drop(solve(crossprod(Z, Z * w) + pen, crossprod(Z, w * eta + (y - p))))
    done <- max(abs(new - beta)) < tol
    beta <- new
    if (done) break
  }
  list(beta = beta, mu = mu, s = s)
}

predict_ridge <- function(fit, X) {
  drop(cbind(1, sweep(sweep(X, 2, fit$mu), 2, fit$s, "/")) %*% fit$beta)
}

auc <- function(score, y) {
  r <- rank(score); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Folds assign whole requests, stratified so every fold has unsafe requests.
grouped_folds <- function(item, unsafe, k = 4) {
  per_item <- unique(data.frame(item, unsafe))
  per_item$fold <- stats::ave(per_item$unsafe, per_item$unsafe,
                              FUN = function(v) sample(rep_len(seq_len(k), length(v))))
  per_item$fold[match(item, per_item$item)]
}

cv_auc <- function(X, y, folds, lambda = 10) {
  score <- numeric(length(y))
  for (f in unique(folds)) {
    tr <- folds != f
    score[!tr] <- predict_ridge(ridge_logit(X[tr, , drop = FALSE], y[tr], lambda), X[!tr, , drop = FALSE])
  }
  auc(score, y)
}

# Observed AUC, and a permutation p-value that relabels whole requests.
tell_audit <- function(bank, k = 4, lambda = 10, B = 200) {
  X <- doc_term_matrix(bank$text)
  folds <- grouped_folds(bank$item, bank$unsafe, k)
  obs <- cv_auc(X, bank$unsafe, folds, lambda)
  item_lab <- unique(bank[, c("item", "unsafe")])
  null <- replicate(B, {
    perm <- sample(item_lab$unsafe)
    y <- perm[match(bank$item, item_lab$item)]
    cv_auc(X, y, grouped_folds(bank$item, y, k), lambda)
  })
  list(auc = obs, null = null, p = (1 + sum(null >= obs)) / (B + 1))
}
