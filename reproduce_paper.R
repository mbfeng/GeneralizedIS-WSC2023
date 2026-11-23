################################################################################
# reproduce_paper.R                                                            #
# -- Entry point: redraws Figures 1-3 and Table 1 of Chen & Feng (WSC 2023)    #
#    from the saved summaries in results/ (a few seconds), then checks Table 1 #
#    against the printed values (verify_against_paper.R)                       #
# -- RERUN_SIMULATIONS <- TRUE first reruns the experiment with new random     #
#    numbers (scripts/run_experiments.R, a few minutes) and uses that rerun    #
# -- Usage (from the Code folder):  Rscript reproduce_paper.R                  #
################################################################################

RERUN_SIMULATIONS <- FALSE
RERUN_SEED <- 2023
source("R/utilities.R")
load_library(".")
source("R/figures_tables.R")
source("verify_against_paper.R")
dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)
par <- gis_parameters()

if (RERUN_SIMULATIONS) {
  folder <- sprintf("results/rerun_seed%d", RERUN_SEED)
  status <- system2("Rscript", c("scripts/run_experiments.R", "--seed", RERUN_SEED,
                                 "--output", folder))
  if (status != 0) {
    stop("scripts/run_experiments.R failed")
  }
  mode <- "monte_carlo"
} else {
  folder <- "results/published_runs"                             # the runs behind the paper
  mode <- "exact"
}
summaries <- read_summaries(folder, par$n_scenarios)
cat("Summaries:", folder, "\n")

# ===== Figures (file names as in the paper's LaTeX source) ============================
plot_figure1a(par, "figures/setting1.pdf")
plot_figure1b(summaries$by_scenario, par, "figures/setting2.pdf")
plot_two_bands(summaries$by_scenario, "STD_n1000", "STD est. with $\\Gamma_{STD} = 10^6$",
               "figures/STDOIS.pdf")
plot_two_bands(summaries$by_scenario, "MLR", "MLR est. with $\\Gamma_{MLR} = 10^3$",
               "figures/MLROIS.pdf")
plot_figure3(summaries$weights, par, "figures/densities.pdf")

# ===== Table 1 ======================================================================
table1 <- table1_amse(summaries)
write.csv(table1, "tables/table1_amse.csv", row.names = FALSE)
write_table1_latex(table1, "tables/table1_amse.tex")

# ===== Verification ====================================================================
cat("\nTable 1 against the printed values (", mode, "mode)\n")
checks <- verify_table1(table1, mode)
write.csv(checks, sprintf("tables/verification_table1_%s.csv", mode), row.names = FALSE)
