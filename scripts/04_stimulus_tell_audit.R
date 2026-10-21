# 04_stimulus_tell_audit.R
# Before piloting: can the agent's justification alone give away which
# requests are unsafe? This runs the audit on a matched pair of toy banks, one
# clean and one where unsafe justifications more often minimize ("It is only
# temporary."), then asks how large a bank must be for the audit to catch that
# tell, and how much the classifier's AUC varies from bank to bank.
# A bank too small to catch a tell can't certify that there isn't one.

source(here::here("R", "load.R"))
library(future.apply)

B <- if (QUICK) 60 else 250     # permutations per audit (p-value resolution 1 / (B + 1))
R <- if (QUICK) 20 else 100     # audits per condition for the power curve
plan(multisession, workers = as.integer(Sys.getenv("SIM_WORKERS", future::availableCores())))

# ---- A matched pair of example banks, 48 requests each ------------------------
# Both banks come from the same seed, fixed in advance (the project seed), so
# they have the same requests, wording, and unsafe assignment and differ only in
# the leak: some unsafe justifications gain a minimizing phrase. Each audit is
# seeded too, so the examples don't depend on the permutation count or on
# anything run before them.
examples <- do.call(rbind, lapply(c(clean = 0, leaky = 0.4), function(leak) {
  bank <- make_text_bank(n_items = 48, leak = leak, seed = 20261019)
  set.seed(20261019)
  au <- tell_audit(bank, B = if (QUICK) 500 else 2000)
  plain <- bank[bank$version == "plain", ]
  minimizes <- grepl("temporary|a moment|Nothing else", plain$text)
  data.frame(leak = leak, auc = au$auc, p = au$p, null_95 = unname(quantile(au$null, 0.95)),
             minimizer_share_safe = mean(minimizes[plain$unsafe == 0L]),
             minimizer_share_unsafe = mean(minimizes[plain$unsafe == 1L]))
}))
print(examples, digits = 3)
write.csv(examples, out_path("tables", "tell_audit_examples.csv"), row.names = FALSE)

# ---- How often does the audit flag the tell, by bank size? ---------------------
# Its own seed, so the power curve doesn't depend on the examples above. Each
# audit's AUC is kept as well as its verdict, to show the spread across banks.
set.seed(20261020)
grid <- expand.grid(n_items = c(48, 96, 192), leak = c(0, 0.4))
pw <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  res <- future_sapply(seq_len(R), function(r) {
    au <- tell_audit(make_text_bank(g$n_items, leak = g$leak), B = B)
    c(flag = au$p < 0.05, auc = au$auc)
  }, future.seed = TRUE)
  auc <- res["auc", ]
  data.frame(g, p_flag = mean(res["flag", ]), audits = R,
             auc_median = stats::median(auc), auc_q25 = unname(stats::quantile(auc, 0.25)),
             auc_q75 = unname(stats::quantile(auc, 0.75)), auc_below_half = mean(auc < 0.5))
}))
print(pw, digits = 3)
write.csv(pw, out_path("tables", "tell_audit_power.csv"), row.names = FALSE)

save_fig(plot_audit(pw), "fig5_tell_audit_power", width = 6, height = 7)
plan(sequential)
writeLines(capture.output(sessionInfo()), out_path("sessionInfo_04.txt"))
