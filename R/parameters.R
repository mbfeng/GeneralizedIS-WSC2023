################################################################################
# parameters.R                                                                 #
# -- Settings of the numerical experiment (Section 5 of Chen & Feng, WSC 2023) #
# -- Reverse iron butterfly on a Black-Scholes stock; risk horizon tau = 0.5   #
# -- Inner and outer variables are LOG prices: x = log S                       #
################################################################################

gis_parameters <- function() {
  # Model, portfolio and experiment settings of Section 5.
  # Output(s): list with
  # -- S0, mu, r, sigma: initial price, real-world drift, risk-free rate, volatility (annual)
  # -- maturity T, horizon tau (years); strikes K = (125, 145, 165)
  # -- n_scenarios m, n_macro (macro replications)
  # -- outer_mean, outer_sd: log S_tau ~ N(outer_mean, outer_sd^2) under P
  # -- inner_drift, inner_sd: log S_T | S_tau ~ N(log S_tau + inner_drift, inner_sd^2) under Q
  # -- price_0: time-0 price P_0 of the portfolio (17.32)
  par <- list(S0 = 100, mu = 0.10, r = 0.05, sigma = 0.30,
              maturity = 1, horizon = 0.5,
              strikes = c(125, 145, 165),
              n_scenarios = 1000, n_macro = 200)
  par$outer_mean <- log(par$S0) + (par$mu - par$sigma^2 / 2) * par$horizon
  par$outer_sd <- par$sigma * sqrt(par$horizon)
  remaining <- par$maturity - par$horizon                        # T - tau
  par$inner_drift <- (par$r - par$sigma^2 / 2) * remaining
  par$inner_sd <- par$sigma * sqrt(remaining)
  par$price_0 <- iron_butterfly_price(par$S0, par$r, par$maturity, par$sigma, par$strikes)
  return(par)
}

outer_scenarios <- function(par) {
  # The m deterministic outer scenarios: the k / (m + 1) quantiles of log S_tau, k = 1..m.
  # Output(s): vector of log S_tau (length m), increasing
  probabilities <- (1:par$n_scenarios) / (par$n_scenarios + 1)
  log_s_tau <- par$outer_mean + par$outer_sd * qnorm(probabilities)
  return(log_s_tau)
}
