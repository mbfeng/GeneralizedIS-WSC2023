################################################################################
# run_experiments.R                                                            #
# -- Monte Carlo rerun of the Section 5 experiment with new random numbers:    #
#    200 macro replications of the OIS, MLR and STD estimators of mu(S_tau,i)  #
#    for the m = 1,000 quantile scenarios, plus the NNLS weights of Figure 3   #
# -- One L'Ecuyer-CMRG stream per macro replication (seed --seed)              #
# -- Writes the same summaries as tools/extract_published_runs.R               #
# -- Usage (from the Code folder):                                             #
#      Rscript scripts/run_experiments.R --seed 2023 \                         #
#        --output results/rerun_seed2023                                       #
#    About 2 minutes on one core.                                              #
################################################################################

source("R/utilities.R")
load_library(".")
settings <- parse_arguments(list(seed = 2023, n_macro = 200,
                                 output = "results/rerun_seed2023"))
par <- gis_parameters()
log_s_tau <- outer_scenarios(par)
mu_true <- true_conditional_mean(log_s_tau, par)
m <- par$n_scenarios
n_macro <- settings$n_macro
names <- estimator_names()

# ===== Macro replications ==========================================================
RNGkind("L'Ecuyer-CMRG")
set.seed(settings$seed)
stream <- .Random.seed
estimates <- list()
for (name in names) {
  estimates[[name]] <- matrix(0, nrow = m, ncol = n_macro)     # [i, k]: scenario i, run k
}
stage2_sizes <- rep(0, n_macro)
start <- Sys.time()
for (k in 1:n_macro) {
  stream <- parallel::nextRNGStream(stream)
  assign(".Random.seed", stream, envir = globalenv())          # start of stream k
  #-- OIS and MLR: Gamma = 10^3 (OIS: 100 Stage I + 900 Stage II samples)
  ois <- ois_estimates(log_s_tau, par, budget = 1000, stage2_share = 0.9)
  estimates$OIS[, k] <- ois$estimates
  stage2_sizes[k] <- ois$stage2_size
  estimates$MLR[, k] <- mlr_estimates(log_s_tau, par, budget = 1000)
  #-- STD: n = 1, 10, 100, 1000 inner samples per scenario (Gamma = 10^3 to 10^6)
  estimates$STD_n1[, k] <- std_estimates(log_s_tau, par, n_inner = 1)
  estimates$STD_n10[, k] <- std_estimates(log_s_tau, par, n_inner = 10)
  estimates$STD_n100[, k] <- std_estimates(log_s_tau, par, n_inner = 100)
  estimates$STD_n1000[, k] <- std_estimates(log_s_tau, par, n_inner = 1000)
  if (k %% 20 == 0) {
    cat(sprintf("macro %d of %d, %.0f s\n", k, n_macro, seconds_since(start)))
  }
}

# ===== NNLS weights of Figure 3 (Gamma_1 = 1,000, stratified: one sample per scenario) ===
stream <- parallel::nextRNGStream(stream)
assign(".Random.seed", stream, envir = globalenv())
weights <- nnls_weights(log_s_tau, par, budget_stage1 = 1000, stratified = TRUE)

# ===== Summaries ==================================================================
summaries <- summarize_runs(estimates, log_s_tau, mu_true)
write_summaries(summaries, weights, settings$output)
cat(sprintf("Stage II sizes of OIS: %d to %d (nominal 900)\n", min(stage2_sizes),
            max(stage2_sizes)))
cat(sprintf("Done in %.0f s; wrote %s\n", seconds_since(start), settings$output))
