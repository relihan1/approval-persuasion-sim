# paired_comparisons.R: coverage by contrast and sample size, and paired changes in
# correct verdicts against the jackknife, from the saved simulations of a run of
# scripts 02 and 06. Run from the project root: RUNQ=<root of that run> Rscript notes/small_sample_quick_checks/paired_comparisons.R
source(here::here("R", "load.R"))
tabs <- file.path(Sys.getenv("RUNQ"), "output", "tables")
s2 <- readRDS(file.path(tabs, "power_simulations.rds"))
pars <- dgm_defaults()
truth <- lapply(names(hypotheses()), function(h) pseudo_true(with_hypothesis(pars, h))$uv); names(truth) <- names(hypotheses())
s2$tv <- mapply(function(h, t) truth[[h]][[t]], s2$truth, s2$test)
z <- qnorm(0.975)
s2$c_jk <- abs(s2$est - s2$tv) <= z * s2$se_two_way_jk
s2$c_jk_sw <- s2$lo95_two_way_jk_sw <= s2$tv & s2$tv <= s2$hi95_two_way_jk_sw
s2$c_boot <- s2$lo95_two_way_boot <= s2$tv & s2$tv <= s2$hi95_two_way_boot
cat("== 02: coverage by test, pooled over truths and n (4000 studies each; MCSE 0.34) ==\n")
print(round(100 * aggregate(cbind(c_jk, c_jk_sw, c_boot) ~ test, s2, mean, na.action = na.pass, na.rm = TRUE)[, -1], 1) |> cbind(test = sort(unique(s2$test))))
cat("== 02: coverage of T3 contrasts by n (pooled over truths and both T3; 1600 per cell; MCSE 0.54) ==\n")
t3 <- s2[s2$test %in% c("T3_plain", "T3_tactic"), ]
print(aggregate(cbind(c_jk, c_jk_sw, c_boot) ~ n_per_cell, t3, function(x) round(100 * mean(x, na.rm = TRUE), 1), na.action = na.pass))
cat("== 02: coverage by test and n, jk_sw / boot ==\n")
print(reshape(aggregate(cbind(jk_sw = c_jk_sw, boot = c_boot) ~ test + n_per_cell, s2, function(x) round(100 * mean(x, na.rm = TRUE), 1), na.action = na.pass), idvar = "test", timevar = "n_per_cell", direction = "wide"))
cat("== 02: per truth, T3 ==\n")
print(aggregate(cbind(c_jk, c_jk_sw, c_boot) ~ truth + test, t3, function(x) round(100 * mean(x, na.rm = TRUE), 1), na.action = na.pass))
# Paired differences in correct verdicts (one row per study).
st <- unique(s2[, c("truth", "n_per_cell", "rep", "verdict_two_way_jk", "verdict_two_way_jk_sw", "verdict_two_way_boot")])
ok <- function(v) as.numeric(v == correct_verdict[st$truth])
st$d_sw <- ok(st$verdict_two_way_jk_sw) - ok(st$verdict_two_way_jk)
st$d_boot <- ok(st$verdict_two_way_boot) - ok(st$verdict_two_way_jk)
pd <- function(x) sprintf("%+.1f (%.1f)", 100 * mean(x), 100 * sd(x) / sqrt(length(x)))
cat("== 02: change in correct-verdict rate vs jackknife, points (paired SE) ==\n")
print(aggregate(cbind(jk_sw = d_sw, boot = d_boot) ~ truth + n_per_cell, st, pd))

s6 <- readRDS(file.path(tabs, "sensitivity_simulations.rds"))
st6 <- unique(s6[, c("scenario", "state", "rep", "verdict_two_way_jk", "verdict_two_way_jk_sw", "verdict_two_way_boot")])
ok6 <- function(v) as.numeric(v == correct_verdict[st6$state])
st6$d_sw <- ok6(st6$verdict_two_way_jk_sw) - ok6(st6$verdict_two_way_jk)
st6$d_boot <- ok6(st6$verdict_two_way_boot) - ok6(st6$verdict_two_way_jk)
cat("== 06: change in correct-verdict rate vs jackknife, points (paired SE) ==\n")
print(aggregate(cbind(jk_sw = d_sw, boot = d_boot) ~ scenario + state, st6, pd))
at <- s6[s6$scenario == "gamma_at_sesoi" & s6$state == "vigilance" & s6$test %in% c("T3_plain", "T3_tactic"), ]
cat("== 06: equivalence at SESOI, counts of 200 ==\n")
print(aggregate(cbind(jk = equiv_two_way_jk, jk_sw = equiv_two_way_jk_sw, boot = equiv_two_way_boot, pid = equiv_pid) ~ test, at, sum))
cat("jk-only (jk says equivalent, jk_sw does not):", sum(at$equiv_two_way_jk & !at$equiv_two_way_jk_sw), "\n")
cat("== 06: mean coverage over contrasts by scenario ==\n")
print(aggregate(cbind(jk = cover_two_way_jk, jk_sw = cover_two_way_jk_sw, boot = cover_two_way_boot) ~ scenario, s6,
                function(x) round(100 * mean(x, na.rm = TRUE), 1), na.action = na.pass))
