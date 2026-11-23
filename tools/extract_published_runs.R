################################################################################
# extract_published_runs.R                                                     #
# -- One-off provenance step: reduces the saved workspace of the published     #
#    runs (original_2023/WSC2023.RData, April 2023, 200 macro replications)    #
#    to the summaries in results/published_runs/                               #
# -- Runs used (names in the workspace -> estimator):                          #
#    result_FIT_1 -> OIS (Gamma = 10^3: 100 + 900), result_MLR_1 -> MLR        #
#    (10^3), result_STD_1/10/100/1000 -> STD with n = 1/10/100/1000            #
#    (10^3 to 10^6); beta_fit -> NNLS weights of Figure 3                      #
# -- The workspace also holds runs with Gamma = 10^4 (result_FIT_10,           #
#    result_MLR_10) that the paper does not use; they are not extracted        #
# -- Usage (from the Code folder): Rscript tools/extract_published_runs.R      #
################################################################################

source("R/utilities.R")
load_library(".")
par <- gis_parameters()
log_s_tau <- outer_scenarios(par)
mu_true <- true_conditional_mean(log_s_tau, par)

saved <- new.env()
load("original_2023/WSC2023.RData", envir = saved)

#-- the saved scenarios and true values are the ones the rebuilt code computes
stopifnot(max(abs(saved$x_out - log_s_tau)) < 1e-12)
stopifnot(max(abs(saved$mu_true - mu_true)) < 1e-12)
stopifnot(saved$N_Macro == par$n_macro)

estimates <- list(OIS = saved$result_FIT_1, MLR = saved$result_MLR_1,
                  STD_n1 = saved$result_STD_1, STD_n10 = saved$result_STD_10,
                  STD_n100 = saved$result_STD_100, STD_n1000 = saved$result_STD_1000)
summaries <- summarize_runs(estimates, log_s_tau, mu_true)
write_summaries(summaries, saved$beta_fit, "results/published_runs")
cat("Wrote results/published_runs/ (", nrow(summaries$by_scenario), "scenarios,",
    nrow(summaries$by_macro), "macro replications )\n")
