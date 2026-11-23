################################################################################
# estimators.R                                                                 #
# -- The three nested estimators of mu(S_tau,i), i = 1..m (Section 5):         #
#    STD (standard nested simulation), MLR (GIS with the equal mixture p-bar,  #
#    eq. (8)) and OIS (GIS with the NNLS weighted mixture q_beta*, Alg. 1)     #
# -- The sampling densities of Figure 3                                        #
# -- Notation: x = log S_T; scenario i has inner pdf p(x | S_tau,i) =          #
#    N(log S_tau,i + inner_drift, inner_sd^2); Gamma = simulation budget       #
# -- Random numbers are drawn in the same order as the original code           #
#    (original_2023/IronButterfly.R), so both give the same estimates for the  #
#    same seed (tests/test_against_original.R)                                 #
################################################################################

# ===== Inner densities ==========================================================

inner_density_matrix <- function(log_s_T, inner_means, inner_sd) {
  # Conditional pdfs of every sample under every scenario.
  # Input(s): log_s_T = samples x_j (length Gamma); inner_means = means of the m inner pdfs
  # Output(s): Gamma-by-m matrix; [j, i] = p(x_j | S_tau,i)
  differences <- outer(log_s_T, inner_means, "-")                # [j, i] = x_j - mean_i
  densities <- dnorm(differences, mean = 0, sd = inner_sd)
  return(densities)
}

mixture_density <- function(log_s_T, component_means, weights, inner_sd) {
  # Weighted mixture q_beta(x) = sum_i beta_i p(x | S_tau,i), eq. (10); components with
  # zero weight are skipped.
  # Input(s): log_s_T = points x; component_means; weights beta (sum 1)
  # Output(s): vector q_beta(x)
  used <- weights > 0
  densities <- inner_density_matrix(log_s_T, component_means[used], inner_sd)
  q_beta <- as.vector(densities %*% weights[used])
  return(q_beta)
}

# ===== STD: standard nested simulation ==========================================

std_estimates <- function(log_s_tau, par, n_inner) {
  # Standard nested estimates: n inner samples per scenario, Gamma = m n.
  # Input(s): log_s_tau = m scenarios; par = gis_parameters(); n_inner = n
  # Output(s): vector of m estimates mu-hat_i
  m <- length(log_s_tau)
  normals <- matrix(rnorm(n_inner * m), nrow = n_inner)         # column i = scenario i
  inner_means <- log_s_tau + par$inner_drift
  log_s_T <- sweep(par$inner_sd * normals, 2, inner_means, "+")  # [j, i] = x_j for scenario i
  outputs <- matrix(simulation_output(log_s_T, par), nrow = n_inner)
  estimates <- colMeans(outputs)
  return(estimates)
}

# ===== MLR: GIS with the equal mixture ===========================================

mlr_estimates <- function(log_s_tau, par, budget) {
  # GIS estimates (6) with q = p-bar, the equal mixture of the m inner pdfs (eq. (8)).
  # Components are chosen at random (not stratified), as in the published runs.
  # Input(s): log_s_tau = m scenarios; par = gis_parameters(); budget = Gamma
  # Output(s): vector of m estimates mu-tilde_i
  m <- length(log_s_tau)
  inner_means <- log_s_tau + par$inner_drift
  component <- sample(1:m, budget, replace = TRUE)
  log_s_T <- inner_means[component] + par$inner_sd * rnorm(budget)
  outputs <- simulation_output(log_s_T, par)
  densities <- inner_density_matrix(log_s_T, inner_means, par$inner_sd)   # [j, i]
  equal_mixture <- rowMeans(densities)                           # p-bar(x_j)
  likelihood_ratios <- densities / equal_mixture                 # [j, i] = p(x_j | i) / p-bar(x_j)
  estimates <- colMeans(outputs * likelihood_ratios)
  return(estimates)
}

# ===== OIS: GIS with the NNLS weighted mixture (Algorithm 1) =====================

nnls_weights <- function(log_s_tau, par, budget_stage1, stratified) {
  # Stage I of Algorithm 1: mixture weights beta* by NNLS, eq. (11).
  # Input(s): log_s_tau = m scenarios; par; budget_stage1 = Gamma_1 samples from p-bar;
  #           stratified = TRUE draws Gamma_1 / m samples from each component (Figure 3),
  #           FALSE chooses components at random (the OIS runs: Gamma_1 = 100 < m)
  # Output(s): vector beta* (length m, sums to 1, mostly zeros)
  # formula: regress |g(x_j)| sqrt(sum_i p^2(x_j | i)) on p(x_j | i), i = 1..m, with
  #          non-negative coefficients; the constant 1/m inside the square root of (11)
  #          only rescales gamma and does not change the normalized weights
  m <- length(log_s_tau)
  inner_means <- log_s_tau + par$inner_drift
  if (stratified) {
    per_component <- budget_stage1 / m
    if (per_component %% 1 != 0) {
      stop("stratified Stage I needs Gamma_1 to be a multiple of m")
    }
    log_s_T <- rep(inner_means, per_component) + par$inner_sd * rnorm(budget_stage1)
  } else {
    component <- sample(1:m, budget_stage1, replace = TRUE)
    log_s_T <- inner_means[component] + par$inner_sd * rnorm(budget_stage1)
  }
  densities <- inner_density_matrix(log_s_T, inner_means, par$inner_sd)   # regressors
  response <- abs(simulation_output(log_s_T, par)) * sqrt(rowSums(densities^2))
  fit <- nnls::nnls(densities, response)
  weights <- fit$x / sum(fit$x)
  return(weights)
}

ois_estimates <- function(log_s_tau, par, budget, stage2_share) {
  # Two-stage GIS estimates of Algorithm 1 (OIS).
  # Input(s): log_s_tau = m scenarios; par; budget = Gamma;
  #           stage2_share = Gamma_2 / Gamma (0.9 in the paper: 100 + 900 samples)
  # Output(s): list
  # -- estimates: m estimates mu-tilde_i
  # -- weights: beta* of Stage I
  # -- stage2_size: number of Stage II samples actually drawn
  # Stage II is stratified: round(beta_i Gamma_2) samples from component i. The rounded
  # counts need not add up to Gamma_2 exactly; the sampling pdf is the mixture with the
  # rounded weights (as in the published runs).
  inner_means <- log_s_tau + par$inner_drift
  #-- Stage I: weights from round((1 - share) Gamma) samples of p-bar
  budget_stage1 <- round((1 - stage2_share) * budget)
  weights <- nnls_weights(log_s_tau, par, budget_stage1, stratified = FALSE)
  #-- Stage II: stratified samples from q_beta*
  used <- which(weights > 0)
  counts <- round(budget * stage2_share * weights[used])        # samples per used component
  stage2_size <- sum(counts)
  rounded_weights <- counts / stage2_size
  component <- rep(seq_along(used), counts)
  log_s_T <- inner_means[used][component] + par$inner_sd * rnorm(stage2_size)
  outputs <- simulation_output(log_s_T, par)
  sampling_pdf <- mixture_density(log_s_T, inner_means[used], rounded_weights, par$inner_sd)
  densities <- inner_density_matrix(log_s_T, inner_means, par$inner_sd)   # [j, i]
  likelihood_ratios <- densities / sampling_pdf
  estimates <- colMeans(outputs * likelihood_ratios)
  return(list(estimates = estimates, weights = weights, stage2_size = stage2_size))
}

# ===== Densities of Figure 3 =====================================================

optimal_unnormalized <- function(log_s_T, component_means, par) {
  # Unnormalized optimal sampling pdf of Proposition 2 (eq. (9)).
  # Input(s): log_s_T = points x; component_means = means of the m inner pdfs; par
  # Output(s): vector |g(x)| sqrt(sum_i p^2(x | i))
  densities <- inner_density_matrix(log_s_T, component_means, par$inner_sd)
  value <- abs(simulation_output(log_s_T, par)) * sqrt(rowSums(densities^2))
  return(value)
}

figure3_densities <- function(log_s_T, log_s_tau, weights, par) {
  # The three sampling pdfs of log S_T in Figure 3: true optimal q*, equal mixture p-bar
  # and NNLS weighted mixture q_beta*.
  # Input(s): log_s_T = grid of points; log_s_tau = m scenarios; weights = beta*
  # Output(s): data frame (log_s_T, optimal, equal_mixture, nnls_mixture)
  # As drawn for the paper, the components are centred at log S_tau,i (without the inner
  # drift (r - sigma^2/2)(T - tau) = 0.0025, which is invisible at the figure's scale).
  # q* (eq. (9)) is |g(x)| sqrt(sum_i p^2(x | i)) normalized by numerical integration over
  # outer_mean +/- 10 outer_sd.
  m <- length(log_s_tau)
  constant <- integrate(optimal_unnormalized, lower = par$outer_mean - 10 * par$outer_sd,
                        upper = par$outer_mean + 10 * par$outer_sd, subdivisions = 1e6,
                        component_means = log_s_tau, par = par)$value
  optimal <- optimal_unnormalized(log_s_T, log_s_tau, par) / constant
  equal_mixture <- mixture_density(log_s_T, log_s_tau, rep(1 / m, m), par$inner_sd)
  nnls_mixture <- mixture_density(log_s_T, log_s_tau, weights, par$inner_sd)
  return(data.frame(log_s_T = log_s_T, optimal = optimal, equal_mixture = equal_mixture,
                    nnls_mixture = nnls_mixture))
}
