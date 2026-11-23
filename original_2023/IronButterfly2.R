# library(functional)
require(nnls)
require(ggplot2)
require(reshape2)
require(latex2exp)
#
rm(list = ls())
######  parameters for Black-Scholes model and experiment design
S0 = 100        # current stock price
mu = 10e-2       # annual real-world rate of return
r = 5e-2        # annual risk-free rate
vol = 25e-2     # annual volatility

T = 2           # time-to-maturity from time 0
tau = 6/12      # risk horizon time-tau

mu_outer = log(S0) + (mu-0.5*vol^2)*tau
sig_outer = vol*sqrt(tau)
sig_inner = vol*sqrt(T-tau)

N_Out = 1e3     # number of outer scenarios
N_In = 1e2      # number of inner scenarios
N_Macro = 1e2   # number of macro replication (for MSE calculation)

# strikes in iron butterfly portfolio
K = c(70, 80, 90, 140, 150, 160)

OptionPrice_BS <- function(S, K, r, t, sigma, option_type = "C"){
    # Black-Scholes option pricing formula function

    d1 = (log(S/K) + (r + sigma^2/2)*t) / (sigma * sqrt(t))
    d2 = d1 - sigma * sqrt(t)

    if(option_type == "C") {
        price = S * pnorm(d1) - K * exp(-r*t) * pnorm(d2)
    } else if(option_type == "P") {
        price = K * exp(-r*t) * pnorm(-d2) - S * pnorm(-d1)
    } else {
        stop("Option type must be 'C' or 'P'")
    }
    return(price)
}

Price_IB2 <- function(S, r, t, sigma, K){
    # price of iron butterfly option portfolio
    price = (-OptionPrice_BS(S, K[1], r, t, sigma, "P") +
                 OptionPrice_BS(S, K[2], r, t, sigma, "P") +
                 OptionPrice_BS(S, K[2], r, t, sigma, "C") -
                 OptionPrice_BS(S, K[3], r, t, sigma, "C")) +
        (-OptionPrice_BS(S, K[4], r, t, sigma, "P") +
             OptionPrice_BS(S, K[5], r, t, sigma, "P") +
             OptionPrice_BS(S, K[5], r, t, sigma, "C") -
             OptionPrice_BS(S, K[6], r, t, sigma, "C"))

    return(price)
}

# time-0 price of the iron butterfly portfolio
IB2_0 = Price_IB2(S0, r, T, vol, K)

Payoff_IB2 <- function(x, K, r, T, tau, Time0Price){
    # Terminal payoff of an iron butterfly option portfolio given
    # the terminal stock return x

    return(exp(-(T-tau)*r)*(- pmax(K[1]-exp(x),0) +
                                pmax(K[2]-exp(x),0) +
                                pmax(exp(x)-K[2],0) -
                                pmax(exp(x)-K[3],0) -
                                pmax(K[4]-exp(x),0) +
                                pmax(K[5]-exp(x),0) +
                                pmax(exp(x)-K[5],0) -
                                pmax(exp(x)-K[6],0)) - Time0Price*exp(r*tau))
}

# N_Outer outer scenarios from the outer distribution's respective quantiles
# Deterministic scenarios
q_tiles = (1:N_Out)/(N_Out+1)
x_out = mu_outer + sig_outer*qnorm(q_tiles)
S_out = exp(x_out)

# Plot problem settings
# -- simulation g(S_T)
ggplot() +
    xlim(exp(mu_outer - 3*sig_outer), exp(mu_outer + 3*sig_outer)) +
    geom_function(fun = function(x) Payoff_IB2(log(x), K, r, T, tau, IB2_0),
                  linewidth=2, n = 1e4) +
    labs(x = TeX("Terminal price $S_T$"),
         y = TeX("Simulation output, $g(S_T)$")) +
    theme(plot.title = element_blank(),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 14))
# -- outer scenario distribution and true value
df1 = data.frame(x = S_out,
                 mu_true = Price_IB2(S_out, r, T-tau, vol, K) - IB2_0*exp(r*tau),
                 outer_dist = dlnorm(S_out, meanlog = mu_outer, sdlog = sig_outer))
df1_long = melt(df1, id.var='x')
ggplot(df1, aes(x=x)) +
    geom_line(aes(y=mu_true), linewidth=2, color="#F8766D") +
    geom_line(aes(y=outer_dist/2e-2-2), linewidth=2, color="#619CFF")+
    geom_ribbon(aes(ymax=outer_dist/2e-2-2), ymin=-2, alpha=0.3, fill="#619CFF") +
    scale_y_continuous(
        # Features of the first axis
        name = TeX("True cond. expectation, $\\mu(S_\\tau)$"),
        limits = c(-2,1.5),
        # Add a second axis and specify its features
        sec.axis = sec_axis(~(.+2)*2e-2, name=TeX("Outer distribution, $\\pi(S_\\tau)$"))
    ) +
    theme(plot.title = element_blank(),
          axis.title.x = element_text(size = 18),
          axis.title.y.left = element_text(color = "#F8766D", size=15),
          axis.title.y.right = element_text(color = "#619CFF", size=15),
          axis.text.y.left = element_text(size = 14,color = "#F8766D"),
          axis.text.y.right = element_text(size = 14,color = "#619CFF")) +
    labs(x = TeX("Outer scenario $S_\\tau$"))
plot(S_out,Payoff_IB2(x_out, K, r, T, tau, IB2_0))

# # Stochastic scenarios
# x_out = log(S0) + (mu-0.5*vol^2)*tau + vol*sqrt(tau)*rnorm(N_Out)
# S_out = exp(x_out)

result_true = Price_IB2(S_out, r, T-tau, vol, K) - IB2_0*exp(r*tau)
plot(S_out, result_true)

NGS_STD <- function(S_out, r, t, sigma, num_sample = 1e2, payoff_fun){
    N_Out = length(S_out)
    C_tau = rep(0, N_Out)

    inner_mean = (r-0.5*sigma^2)*t
    inner_vol = sigma*sqrt(t)

    for (i in 1:N_Out){
        xT = log(S_out[i]) + inner_mean + inner_vol*rnorm(num_sample)
        C_tau[i] = mean(payoff_fun(xT))
    }

    return(C_tau)
}

result_STD = NGS_STD(S_out, r, T-tau, vol, N_In, function(x) Payoff_IB2(x, K, r, T, tau, IB2_0))
plot(S_out, result_STD)
lines(S_out, result_true, col='blue',lwd  = 2)

g_c <- function(x, x_out = x_out, sd = sig_inner){

    result = outer(x, x_out, function(x,y) dnorm(x = x, mean = y, sd = sd))
    return(sqrt(rowSums(result^2)))
}
g_c_const <- function(x_out = x_out){

    const = integrate(f = function (x) g_c(x, x_out = x_out),
                      lower = mu_outer - 10*sig_outer,
                      upper = mu_outer + 10*sig_outer,
                      subdivisions = 1e6)$value

    return(const)
}

g_star <- function(x, x_out = x_out, h = function(x) Payoff_IB2(x, K, r, T, tau, IB2_0)){
    return(abs(h(x)) * g_c(x, x_out))
}
g_star_const <- function(x_out = x_out){

    const = integrate(f = function (x) g_star(x, x_out = x_out, h = function(x) Payoff_IB2(x, K, r, T, tau, IB2_0)),
                      lower = mu_outer - 10*sig_outer,
                      upper = mu_outer + 10*sig_outer,
                      subdivisions = 1e6)$value

    return(const)
}

x_plot = seq(from = mu_outer - 5*sig_outer, to = mu_outer + 5*sig_outer, length.out = 1e4)
S_plot = exp(x_plot)


g_c_plot = g_c(x = x_plot, x_out = x_plot)/g_c_const(x_out = x_plot)
g_star_plot = g_star(x = x_plot, x_out = x_plot)/g_star_const(x_out = x_plot)

plot(x_plot, g_c_plot)
plot(S_plot, g_c_plot)

plot(x_plot, g_star_plot)
plot(S_plot, g_star_plot)

tab_wide = data.frame(x = x_plot, S = S_plot,
                      g_c = g_c_plot, g_star = g_star_plot)
tab_melted = melt(tab_wide, id = c("x", "S"))

cols = c(1, 2)
ggplot(tab_melted) +
    geom_line(aes(x=S, y=value, color=variable, group=variable)) +
    labs(x = "One-year Stock Price", y="Density",
         title = "Comparison between different proposal densities") +
    theme(plot.title=element_text(hjust = 0.5))

f_beta <- function(x, x_out, beta, sd = sig_inner){
    ind = (beta > 0)
    x_out = x_out[ind]
    beta = beta[ind]

    f_mat = outer(x, x_out, function(x,y) dnorm(x = x, mean = y, sd = sd))

    return(f_mat %*% beta)
}

NGS_MLR <- function(S_out, r, t, sigma, num_sample = 1e2, payoff_fun, strat = FALSE){
    N_Out = length(S_out)

    n_per_out = num_sample / N_Out
    if(n_per_out %% 1 != 0 && strat){
        warning("Total number of sample is not a multiple of the number of outer scenarios.
                Stratification is turned off.")
        strat = FALSE
    }

    inner_mean = (r-0.5*sigma^2)*t
    inner_vol = sigma*sqrt(t)

    # Sampling Mixture r.v.
    if (strat){
        xT = rep(log(S_out) + inner_mean, n_per_out) + inner_vol*rnorm(num_sample)
    }
    else{
        U = sample((1:N_Out), num_sample, replace=TRUE)
        xT = log(S_out[U]) + inner_mean + inner_vol*rnorm(num_sample)
    }
    C_tau_sample = exp(-r*t)*payoff_fun(xT)

    # Calculate likelihoods
    beta = rep(1/N_Out, N_Out)
    llh.g = f_beta(x = xT, x_out = log(S_out) + inner_mean,
                   beta = beta, sd = inner_vol)

    # Calculate option payoffs
    C_tau = rep(0, N_Out)
    C_tau_sn = rep(0, N_Out)
    for (i in 1:N_Out){
        llh.f = dnorm(xT, log(S_out[i]) + inner_mean, inner_vol)
        lr = llh.f / llh.g

        C_tau[i] = mean(C_tau_sample * lr)
        C_tau_sn[i] = C_tau[i] / mean(lr)
    }

    return(list(C_tau = C_tau, C_tau_sn = C_tau_sn))
}

result_MLR = NGS_MLR(S_out, r, T-tau, vol, N_In*N_Out, function(x) Payoff_IB2(x, K, r, T, tau, IB2_0), strat = FALSE)
plot(S_out, result_MLR$C_tau)
lines(S_out, result_true, col='blue',lwd  = 2)
plot(S_out, result_MLR$C_tau_sn)
lines(S_out, result_true, col='blue',lwd  = 2)

result_MLR_strat = NGS_MLR(S_out, r, T-tau, vol, N_In*N_Out, function(x) Payoff_IB2(x, K, r, T, tau, IB2_0), strat = TRUE)
plot(S_out, result_MLR_strat$C_tau)
lines(S_out, result_true, col='blue',lwd  = 2)
plot(S_out, result_MLR_strat$C_tau_sn)
lines(S_out, result_true, col='blue',lwd  = 2)

# fit beta for OIS
ls.fitting <- function(S_out, r, t, sigma, num_sample = 1e2, payoff_fun, strat = FALSE){
    N_Out = length(S_out)

    n_per_out = num_sample / N_Out
    if(n_per_out %% 1 != 0 && strat){
        warning("Total number of sample is not a multiple of the number of outer scenarios.
                Stratification is turned off.")
        strat = FALSE
    }

    inner_mean = (r-0.5*sigma^2)*t
    inner_vol = sigma*sqrt(t)

    # Sampling Mixture r.v.
    if (strat){
        xT = rep(log(S_out) + inner_mean, n_per_out) + inner_vol*rnorm(num_sample)
    }
    else{
        U = sample((1:N_Out), num_sample, replace=TRUE)
        xT = log(S_out[U]) + inner_mean + inner_vol*rnorm(num_sample)
    }

    # Computing L2 sum
    X = outer(xT, log(S_out) + inner_mean,
              function(x,y) dnorm(x = x, mean = y, sd = inner_vol))
    Y = abs(payoff_fun(xT))*sqrt(rowSums(X^2))
    beta = nnls(X, Y)$x
    beta = beta / sum(beta)

    return(beta)
}

# LS with strat
NGS_LS <- function(S_out, r, t, sigma, num_sample = 1e2, payoff_fun, p, strat = FALSE){
    N_Out = length(S_out)

    # Estimate beta and sample
    num_sample_beta = round((1-p)*num_sample)
    sprintf("%d samples are used to estimate beta.", num_sample_beta)

    beta = ls.fitting(S_out, r, t, sigma, num_sample = num_sample_beta, payoff_fun, strat = FALSE)
    NZ_ind = (beta > 0)
    S_out_NZ = S_out[NZ_ind]
    beta_NZ = beta[NZ_ind]

    num_sample_x = num_sample*p
    N_Out_x = length(S_out_NZ)

    num_sample_rounded = round(num_sample_x * beta_NZ)
    beta_rounded = num_sample_rounded / sum(num_sample_rounded)

    # Sampling Mixture r.v.
    inner_mean = (r-0.5*sigma^2)*t
    inner_vol = sigma*sqrt(t)
    if (strat){
        U = rep((1:N_Out_x), num_sample_rounded)
    }
    else{
        U = sample((1:N_Out_x), num_sample_x, replace=TRUE)
    }
    xT = log(S_out_NZ[U]) + inner_mean + inner_vol*rnorm(num_sample_x)
    C_tau_sample = exp(-r*t)*payoff_fun(xT)

    # Calculate likelihoods
    llh.g = f_beta(xT, log(S_out_NZ) + inner_mean, beta_rounded, sd = inner_vol)

    # Calculate option payoffs
    C_tau = rep(0, N_Out)
    C_tau_sn = rep(0, N_Out)
    for (i in 1:N_Out){
        llh.f = dnorm(xT, log(S_out[i]) + inner_mean, inner_vol)
        lr = llh.f / llh.g

        C_tau[i] = mean(C_tau_sample * lr)
        C_tau_sn[i] = C_tau[i] / mean(lr)
    }

    return(list(C_tau = C_tau, C_tau_sn = C_tau_sn))
}

result_LS = NGS_LS(S_out, r, T-tau, vol, N_In*N_Out, function(x) Payoff_IB2(x, K, r, T, tau, IB2_0),
                   p = 0.99, strat = FALSE)

plot(S_out, result_LS$C_tau)
lines(S_out, result_true, col='blue',lwd  = 2)
plot(S_out, result_LS$C_tau_sn)
lines(S_out, result_true, col='blue',lwd  = 2)

result_LS_strat = NGS_LS(S_out, r, T-tau, vol, N_In*N_Out, function(x) Payoff_IB2(x, K, r, T, tau, IB2_0),
                         p = 0.99, strat = TRUE)
plot(S_out, result_LS_strat$C_tau)
lines(S_out, result_true, col='blue',lwd  = 2)
plot(S_out, result_LS_strat$C_tau_sn)
lines(S_out, result_true, col='blue',lwd  = 2)

AR <- function(S_out, r, t, sigma, num_sample = 1e2, payoff_fun, big_M, strat = FALSE){
    N_Out=length(S_out)

    n_per_out = num_sample / N_Out
    if(n_per_out %% 1 != 0 && strat){
        warning("Total number of sample is not a multiple of the number of outer scenarios.
                Stratification is turned off.")
        strat = FALSE
    }

    inner_mean = (r-0.5*sigma^2)*t
    inner_vol = sigma*sqrt(t)

    # Sampling Mixture r.v.
    if (strat){
        xT = rep(log(S_out) + inner_mean, n_per_out) + inner_vol*rnorm(num_sample)
    }
    else{
        U = sample((1:N_Out), num_sample, replace = TRUE)
        xT = log(S_out[U]) + inner_mean + inner_vol*rnorm(num_sample)
    }

    f_mat = outer(xT, log(S_out) + inner_mean,
                  function(x,y) dnorm(x = x, mean = y, sd = inner_vol))

    f_bar_denom = rowSums(f_mat)
    g_star_numer = abs(payoff_fun(xT))*sqrt(rowSums(f_mat^2))

    Z = runif(num_sample)
    accept_ind = (Z < (g_star_numer/(big_M*f_bar_denom)))

    accept_rate = sum(accept_ind)/length(accept_ind)
    c = 1/(N_Out*accept_rate*big_M)

    xT = xT[accept_ind]
    llh.g = c * g_star_numer[accept_ind]
    f_mat = f_mat[accept_ind, ,drop = F]

    return(list(xT = xT, llh.g = llh.g, f_mat = f_mat,
                accept_rate = accept_rate, c = c))
}

NGS_OIS <- function(S_out, r, t, sigma, num_sample = 1e2, payoff_fun, big_M, strat = FALSE, fixed_num_samp = FALSE){
    N_Out=length(S_out)

    AR_results = AR(S_out, r, t, sigma, num_sample, payoff_fun, big_M, strat)
    xT = AR_results$xT
    llh.g = AR_results$llh.g
    f_mat = AR_results$f_mat
    accept_rate = AR_results$accept_rate
    c = AR_results$c

    if(fixed_num_samp){
        num_total = num_sample
        xT_len = length(xT)
        while(xT_len < num_sample){
            # estimate the number of AR samples needed based on the acceptance rate
            shortage = num_sample - length(xT)
            num_sample_cur = min(shortage / accept_rate, num_sample)

            if(strat){
                num_sample_cur = ceiling(num_sample_cur/N_Out) * N_Out
            }

            AR_results = AR(S_out, r, t, sigma, num_sample_cur, payoff_fun, big_M, strat)
            xT = c(xT, AR_results$xT)
            llh.g = c(llh.g, AR_results$llh.g)
            f_mat = rbind(f_mat, AR_results$f_mat)

            num_total = num_total + num_sample_cur
            xT_len = length(xT)

            accept_rate = xT_len/num_total
            c = 1/(N_Out*accept_rate*big_M)
        }

        # trim down sample size to match num_samp
        ind = sample((1:xT_len), num_sample, replace = FALSE)
        xT = xT[ind]
        llh.g = llh.g[ind]
        f_mat = f_mat[ind, ,drop = F]
    }

    C_tau_sample = exp(-r*t)*payoff_fun(xT)

    # Calculate option payoffs
    C_tau = rep(0, N_Out)
    C_tau_sn = rep(0, N_Out)
    for (i in 1:N_Out){
        llh.f = f_mat[ , i]
        lr = llh.f / llh.g

        C_tau[i] = mean(C_tau_sample * lr)
        C_tau_sn[i] = C_tau[i] / mean(lr)
    }

    return(list(C_tau = C_tau, C_tau_sn = C_tau_sn, accept_rate = accept_rate, c_hat = c))
}

result_OIS = NGS_OIS(S_out, r, T-tau, vol, N_In*N_Out, function(x) Payoff_IB2(x, K, r, T, tau, IB2_0),
                     big_M = abs(K1-K2), strat = FALSE, fixed_num_samp = FALSE)
plot(S_out, result_OIS$C_tau)
lines(S_out, result_true, col='blue',lwd  = 2)
plot(S_out, result_OIS$C_tau_sn)
lines(S_out, result_true, col='blue',lwd  = 2)
print(result_OIS$accept_rate)

result_OIS_strat = NGS_OIS(S_out, r, T-tau, vol, N_In*N_Out, function(x) Payoff_IB2(x, K, r, T, tau, IB2_0),
                           big_M = abs(K1-K2), strat = FALSE)
plot(S_out, result_OIS_strat$C_tau)
lines(S_out, result_true, col='blue',lwd  = 2)
plot(S_out, result_OIS_strat$C_tau_sn)
lines(S_out, result_true, col='blue',lwd  = 2)
print(result_OIS_strat$accept_rate)

