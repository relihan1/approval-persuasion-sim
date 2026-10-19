# 04_stimulus_tell_audit.R
# Before piloting: can the agent's justification alone give away which
# requests are unsafe? This runs the audit on two toy banks, one clean and one
# where unsafe justifications more often minimize ("It is only temporary."),
# then asks how large a bank must be for the audit to catch that tell.
# A bank too small to catch a tell can't certify that there isn't one.

source(here::here("R", "load.R"))
library(future.apply)
set.seed(20261019)

B <- if (QUICK) 60 else 250     # permutations per audit (p-value resolution 1 / (B + 1))
R <- if (QUICK) 20 else 100     # audits per condition for the power curve
plan(multisession, workers = as.integer(Sys.getenv("SIM_WORKERS", future::availableCores())))

# ---- One audit each: a clean bank and a leaky bank of 48 requests --------------
examples <- do.call(rbind, lapply(c(clean = 0, leaky = 0.4), function(leak) {
  bank <- make_text_bank(n_items = 48, leak = leak)
  au <- tell_audit(bank, B = if (QUICK) 500 else 2000)
  data.frame(leak = leak, auc = au$auc, p = au$p, null_95 = unname(quantile(au$null, 0.95)))
}))
print(examples, digits = 3)
write.csv(examples, out_path("tables", "tell_audit_examples.csv"), row.names = FALSE)

# ---- How often does the audit flag the tell, by bank size? ---------------------
grid <- expand.grid(n_items = c(48, 96, 192), leak = c(0, 0.4))
pw <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  flags <- future_sapply(seq_len(R), function(r) {
    tell_audit(make_text_bank(g$n_items, leak = g$leak), B = B)$p < 0.05
  }, future.seed = TRUE)
  data.frame(g, p_flag = mean(flags), audits = R)
}))
print(pw, digits = 3)
write.csv(pw, out_path("tables", "tell_audit_power.csv"), row.names = FALSE)

save_fig(plot_audit(pw), "fig5_tell_audit_power", width = 6, height = 4.2)
plan(sequential)
writeLines(capture.output(sessionInfo()), out_path("sessionInfo_04.txt"))
