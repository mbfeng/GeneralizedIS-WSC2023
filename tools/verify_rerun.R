################################################################################
# verify_rerun.R                                                               #
# -- Checks the summaries of one rerun against Table 1 (Monte Carlo mode) and  #
#    reports how heavy-tailed the per-macro squared errors are                 #
#    (top1_share = share of the AMSE from the worst 1% of macro replications)  #
# -- Usage (from the Code folder):                                             #
#    Rscript tools/verify_rerun.R results/rerun_seed2024                       #
################################################################################

source("R/utilities.R")
load_library(".")
source("R/figures_tables.R")
source("verify_against_paper.R")
folder <- commandArgs(trailingOnly = TRUE)[1]
par <- gis_parameters()
summaries <- read_summaries(folder, par$n_scenarios)
table1 <- table1_amse(summaries)
cat("Summaries:", folder, "\n")
checks <- verify_table1(table1, "monte_carlo")

cat("\nShare of each estimator's AMSE from the worst 1% of macro replications\n")
for (name in estimator_names()) {
  squared_errors <- summaries$by_macro[[paste0(name, "_sq_error")]]
  n_top <- ceiling(0.01 * length(squared_errors))
  worst <- sort(squared_errors, decreasing = TRUE)[1:n_top]
  cat(sprintf("  %-10s top1_share %.3f\n", name, sum(worst) / sum(squared_errors)))
}
