################################################################################
# test_against_original.R                                                      #
# -- Checks that the rebuilt functions give the same numbers as the original   #
#    code (original_2023/IronButterfly.R) for the same random numbers, and     #
#    that the deterministic inputs of Figures 1 and 3 equal the saved ones     #
#    (original_2023/WSC2023.RData)                                             #
# -- Usage (from the Code folder): Rscript tests/test_against_original.R       #
################################################################################

source("R/utilities.R")
load_library(".")
par <- gis_parameters()
log_s_tau <- outer_scenarios(par)

original <- new.env()
sys.source("original_2023/IronButterfly.R", envir = original)   # defines its functions there
payoff_original <- function(x) {
  return(original$Payoff_IB(x, original$K, original$r, original$T, original$tau, original$IB_0))
}
S_out <- original$S_out
remaining <- original$T - original$tau

n_failed <- 0
check <- function(label, rebuilt, reference, tolerance = 1e-10) {
  difference <- max(abs(rebuilt - reference))
  # Prints one line; returns 1 if the check failed, 0 otherwise.
  failed <- difference > tolerance
  status <- if (failed) "FAILED" else "ok"
  cat(sprintf("%-44s max abs difference %.2e  %s\n", label, difference, status))
  return(as.numeric(failed))
}

# ===== Deterministic quantities ===================================================
n_failed <- n_failed + check("time-0 price P_0", par$price_0, original$IB_0)
n_failed <- n_failed + check("scenarios log S_tau", log_s_tau, original$x_out)
grid <- seq(4, 6, length.out = 500)
n_failed <- n_failed + check("simulation output g(S_T)", simulation_output(grid, par),
                             payoff_original(grid))
reference <- original$Price_IB(S_out, original$r, remaining, original$vol, original$K) -
  original$IB_0 * exp(original$r * original$tau)
n_failed <- n_failed + check("true values mu(S_tau)", true_conditional_mean(log_s_tau, par),
                             reference)

# ===== Estimators on the same random numbers ========================================
for (n in c(1, 10, 100)) {
  set.seed(11)
  rebuilt <- std_estimates(log_s_tau, par, n_inner = n)
  set.seed(11)
  reference <- original$NGS_STD(S_out, original$r, remaining, original$vol, num_sample = n,
                                payoff_original)
  n_failed <- n_failed + check(sprintf("STD, n = %d", n), rebuilt, reference)
}

set.seed(12)
rebuilt <- mlr_estimates(log_s_tau, par, budget = 1000)
set.seed(12)
reference <- original$NGS_MLR(S_out, original$r, remaining, original$vol, num_sample = 1000,
                              payoff_original, strat = FALSE)$C_tau
n_failed <- n_failed + check("MLR, Gamma = 10^3", rebuilt, reference)

for (stratified in c(FALSE, TRUE)) {
  set.seed(13)
  rebuilt <- nnls_weights(log_s_tau, par, budget_stage1 = 1000, stratified = stratified)
  set.seed(13)
  reference <- original$ls.fitting(S_out, original$r, remaining, original$vol,
                                   num_sample = 1000, payoff_original, strat = stratified)
  label <- sprintf("NNLS weights, stratified = %s", stratified)
  n_failed <- n_failed + check(label, rebuilt, reference)
}

for (seed in c(14, 15, 16)) {
  set.seed(seed)
  rebuilt <- ois_estimates(log_s_tau, par, budget = 1000, stage2_share = 0.9)$estimates
  set.seed(seed)
  reference <- original$NGS_LS(S_out, original$r, remaining, original$vol, num_sample = 1000,
                               payoff_original, p = 0.9, strat = TRUE)$C_tau
  n_failed <- n_failed + check(sprintf("OIS, Gamma = 10^3, seed %d", seed), rebuilt, reference)
}

# ===== Figure 3 densities against the saved curves ====================================
# (the saved workspace is kept with the archive, not in the git repository)
workspace <- "original_2023/WSC2023.RData"
if (file.exists(workspace)) {
  saved <- new.env()
  load(workspace, envir = saved)
  densities <- figure3_densities(saved$x_plot, log_s_tau, saved$beta_fit, par)
  n_failed <- n_failed + check("Figure 3: optimal q*", densities$optimal, saved$g_star_plot,
                               1e-8)
  n_failed <- n_failed + check("Figure 3: equal mixture p-bar", densities$equal_mixture,
                               saved$f_mix_plot)
  n_failed <- n_failed + check("Figure 3: NNLS mixture q_beta*", densities$nnls_mixture,
                               saved$f_beta_plot)
  n_failed <- n_failed + check("g(S_T) on the saved plotting grid",
                               simulation_output(saved$x_plot, par), saved$gST)
} else {
  cat("Skipped the checks against", workspace, "(not present)\n")
}

cat(if (n_failed == 0) "All checks ok\n" else sprintf("%d checks FAILED\n", n_failed))
