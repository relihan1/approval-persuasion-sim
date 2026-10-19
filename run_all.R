# run_all.R: run the tests, then every script in order.
# Quick mode (default) takes about 30 minutes on two cores.
# Full mode: SIM_FULL=1 Rscript run_all.R  (set SIM_WORKERS to the cores you can spare)
# The full-size brms fit is separate: RUN_BRMS_FULL=1 Rscript scripts/05_brms_primary_model.R

testthat::test_dir(here::here("tests", "testthat"), stop_on_failure = TRUE)

for (s in c("01_simulate_and_recover.R", "02_power_and_adjudication.R",
            "03_simulated_reviewers_ppi.R", "04_stimulus_tell_audit.R",
            "05_brms_primary_model.R", "06_sensitivity.R")) {
  message("\n==== ", s, " ====")
  t0 <- Sys.time()
  source(here::here("scripts", s), echo = FALSE, local = new.env())
  message(sprintf("%s finished in %.1f minutes", s, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
