# df_calibration.R: compares the Satterthwaite degrees of freedom for the request
# part of the jackknifed two-way variance with 2 mean(V)^2 / var(V) across
# simulated studies. Run from the project root: OUTDIR=<folder> Rscript notes/small_sample_quick_checks/df_calibration.R
# Do the Satterthwaite df describe how much the request part of the variance
# varies from study to study? Empirical df = 2 mean(V)^2 / var(V) across
# simulated studies, against the mean of the per-study estimates.
source(here::here("R", "load.R"))
library(future.apply)
plan(multisession, workers = 4)
L <- contrast_L()
scen <- list(default = list(), request_slopes = list(omega_t = 0.3, omega_k = 0.2),
             independent = list(tau_c = 0, tau_d = 0, omega_c = 0, omega_d = 0))
R <- 250
grid <- expand.grid(scen = names(scen), rep = seq_len(R), stringsAsFactors = FALSE)
one <- function(i) {
  p <- utils::modifyList(with_hypothesis(dgm_defaults(), "vigilance"), c(list(n_per_cell = 125), scen[[grid$scen[i]]]))
  d <- simulate_study(p, seed = 9e6 + i)
  a <- aggregate_reviews(d); f <- fit_uv(a$agg)
  v <- vcov_two_way_jk(f, a, d, L)
  ei <- expected_info(f, a$agg, a$key); H <- cluster_info(ei, d$item); C <- contrast_full(f, L)
  nu_req <- apply(C, 2, function(cv) satterthwaite_df(H, cv, jackknife = TRUE))
  data.frame(scen = grid$scen[i], rep = grid$rep[i], test = rownames(L), v_req = attr(v, "request"),
             v_tot = as.numeric(v), nu_req = nu_req, nu_tot = two_way_df(H, C, as.numeric(v), attr(v, "request"), TRUE),
             est = drop(L %*% coef(f)[colnames(L)]))
}
t0 <- Sys.time()
res <- do.call(rbind, future_lapply(seq_len(nrow(grid)), one, future.seed = TRUE, future.packages = "ordinal"))
cat(sprintf("%.1f minutes\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
saveRDS(res, file.path(Sys.getenv("OUTDIR"), "df_calibration.rds"))
emp <- function(x) 2 * mean(x)^2 / var(x)
out <- do.call(rbind, lapply(split(res, list(res$scen, res$test)), function(d) data.frame(
  scen = d$scen[1], test = d$test[1],
  empirical_df_req = emp(d$v_req), mean_nu_req = mean(d$nu_req), median_nu_req = median(d$nu_req),
  empirical_df_tot = emp(d$v_tot[d$v_tot > 0]), median_nu_tot = median(d$nu_tot, na.rm = TRUE),
  sd_est_over_rms_se = sd(d$est) / sqrt(mean(d$v_tot[d$v_tot > 0])))))
out <- out[order(out$scen, out$test), ]
print(out, digits = 3, row.names = FALSE)
write.csv(out, file.path(Sys.getenv("OUTDIR"), "df_calibration.csv"), row.names = FALSE)
