# Run from the project root: Rscript tests/testthat.R
testthat::test_dir(here::here("tests", "testthat"), stop_on_failure = TRUE)
