# load.R: attach packages, set options, and source the project functions.
# Every script starts with: source(here::here("R", "load.R"))

suppressPackageStartupMessages({
  library(ordinal)    # clm(): cumulative link models with a scale term
  library(ggplot2)
  library(patchwork)
})

options(contrasts = c("contr.treatment", "contr.poly"))  # coefficient names depend on this

for (f in c("design.R", "simulate.R", "analyze.R", "ppi.R", "audit.R", "plots.R")) {
  path <- here::here("R", f)
  if (file.exists(path)) source(path)
}

# QUICK = TRUE (the default) runs every script in minutes on a laptop.
# Set the environment variable SIM_FULL=1 for the full replication counts.
QUICK <- !identical(Sys.getenv("SIM_FULL"), "1")

out_path <- function(...) {
  p <- here::here("output", ...)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  p
}
