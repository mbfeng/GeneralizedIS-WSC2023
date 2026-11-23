################################################################################
# option_portfolio.R                                                           #
# -- Black-Scholes prices, the simulation model g(S_T) and the true values     #
#    mu(S_tau) of Section 5                                                    #
# -- Portfolio, as run: short put K1, long put K2, long call K2, short call K3 #
#    (K = 125, 145, 165); time-0 price P_0 = 17.32                             #
################################################################################

bs_price <- function(spot, strike, r, time_left, sigma, type) {
  # Black-Scholes price of a European call or put.
  # Input(s): spot (vector), strike, rate r, time to maturity, volatility, type "call"/"put"
  # Output(s): vector of prices
  d1 <- (log(spot / strike) + (r + sigma^2 / 2) * time_left) / (sigma * sqrt(time_left))
  d2 <- d1 - sigma * sqrt(time_left)
  if (type == "call") {
    price <- spot * pnorm(d1) - strike * exp(-r * time_left) * pnorm(d2)
  } else {
    price <- strike * exp(-r * time_left) * pnorm(-d2) - spot * pnorm(-d1)
  }
  return(price)
}

iron_butterfly_price <- function(spot, r, time_left, sigma, strikes) {
  # Price of the reverse iron butterfly (short put K1, long put K2, long call K2,
  # short call K3).
  # Input(s): spot (vector), rate, time to maturity, volatility, strikes = c(K1, K2, K3)
  # Output(s): vector of prices
  short_put <- bs_price(spot, strikes[1], r, time_left, sigma, "put")
  long_put <- bs_price(spot, strikes[2], r, time_left, sigma, "put")
  long_call <- bs_price(spot, strikes[2], r, time_left, sigma, "call")
  short_call <- bs_price(spot, strikes[3], r, time_left, sigma, "call")
  price <- -short_put + long_put + long_call - short_call
  return(price)
}

simulation_output <- function(log_s_T, par) {
  # Simulation model g(S_T) = exp(-r (T - tau)) P(S_T) - P_0 (P&L at tau of one inner path).
  # Input(s): log_s_T = vector of log terminal prices; par = gis_parameters()
  # Output(s): vector g(S_T)
  # formula: P(S) = -(K1 - S)^+ + (K2 - S)^+ + (S - K2)^+ - (S - K3)^+
  s_T <- exp(log_s_T)
  K <- par$strikes
  payoff <- -pmax(K[1] - s_T, 0) + pmax(K[2] - s_T, 0) + pmax(s_T - K[2], 0) -
    pmax(s_T - K[3], 0)
  discount <- exp(-par$r * (par$maturity - par$horizon))
  g <- discount * payoff - par$price_0 * exp(par$r * par$horizon)
  return(g)
}

true_conditional_mean <- function(log_s_tau, par) {
  # True value mu(S_tau) = E^Q[g(S_T) | S_tau], by Black-Scholes.
  # Input(s): log_s_tau = vector of log scenarios; par = gis_parameters()
  # Output(s): vector mu(S_tau)
  price_tau <- iron_butterfly_price(exp(log_s_tau), par$r, par$maturity - par$horizon,
                                    par$sigma, par$strikes)
  mu_true <- price_tau - par$price_0 * exp(par$r * par$horizon)
  return(mu_true)
}
