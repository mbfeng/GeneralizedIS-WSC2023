source("IronButterfly.R")

N_Macro = 100

# # time-0 price of the iron butterfly portfolio
IB_0 = Price_IB(S0, r, T, vol, K)

# Plot problem settings
q_tiles_plot = (1:N_Out)/(N_Out+1)
x_plot = mu_outer + sig_outer*qnorm(q_tiles_plot)
S_plot = exp(x_plot)

# -- simulation model g(S_T)
gST = Payoff_IB(x_plot, K, r, T, tau, IB_0)
df1 = data.frame(x = S_plot,
                 gST = gST)
ggplot(df1, aes(x=x)) +
    geom_line(aes(y=gST), linewidth=2) +
    labs(x = TeX("Terminal stock price, $S_T$"),
         y = TeX("Simulation model, $g(S_T)$")) +
    theme(plot.title = element_blank(),
          axis.title = element_text(size = 18,face="bold"),
          axis.text = element_text(size = 18,face="bold"),
          plot.margin = margin(0.5, 0.5, 0, 0, "cm"))

# -- outer scenario distribution and true value
mu_true = Price_IB(S_plot, r, T-tau, vol, K) - IB_0*exp(r*tau)
outer_dist = dlnorm(S_plot, meanlog = mu_outer, sdlog = sig_outer)
result_STD_100 = matrix(rep(0,N_Macro*N_Out), ncol = N_Macro, nrow = N_Out)
result_STD_1000 = matrix(rep(0,N_Macro*N_Out), ncol = N_Macro, nrow = N_Out)
for(i in 1:N_Macro){
    print(i)
    result_STD_100[,i] = NGS_STD(S_out, r, T-tau, vol, num_sample = 1e2, function(x) Payoff_IB(x, K, r, T, tau, IB_0))
    result_STD_1000[,i] = NGS_STD(S_out, r, T-tau, vol, num_sample = 1e3, function(x) Payoff_IB(x, K, r, T, tau, IB_0))
}
Qntl_STD_100 = apply(result_STD_100 , 1 , quantile , probs = c(0.025, 0.975) , na.rm = TRUE)
Qntl_STD_1000 = apply(result_STD_1000 , 1 , quantile , probs = c(0.025, 0.975) , na.rm = TRUE)
df2 = data.frame(x = S_out,
                 mu_true = mu_true,
                 outer_dist = outer_dist,
                 lower100 = Qntl_STD_100[1,],
                 upper100 = Qntl_STD_100[2,],
                 lower1000 = Qntl_STD_1000[1,],
                 upper1000 = Qntl_STD_1000[2,])
ggplot(df2, aes(x=x)) +
    geom_ribbon(aes(ymin = lower100, ymax = lower1000, fill = "B"), alpha=0.9)+
    geom_ribbon(aes(ymin = upper1000, ymax = upper100, fill = "B"), alpha=0.9)+
    geom_ribbon(aes(ymin = lower1000, ymax = upper1000, fill = "C"), alpha=0.9)+
    geom_line(aes(y=mu_true), linewidth=1.5) +
    scale_fill_manual(values = c("#fbb4b9","#f768a1"),
                      labels = c(TeX("$n=10^2, \\Gamma_{STD} = 10^5$"),
                                 TeX("$n=10^3, \\Gamma_{STD} = 10^6$")))+
    geom_line(aes(y=outer_dist/1e-2-5), linewidth=1, color="#619CFF")+
    geom_ribbon(aes(ymax=outer_dist/1e-2-5), ymin=-5, alpha=0.3, fill="#619CFF") +
    scale_y_continuous(
        # Features of the first axis
        name = TeX("True value $\\mu(S_\\tau)$"),
        limits = c(-5,2.5),
        # Add a second axis and specify its features
        sec.axis = sec_axis(~(.+5)*1e-2, name=TeX("Outer pdf $\\pi(S_\\tau)$"))
    ) +
    theme(legend.background = element_rect(fill="transparent"),
          legend.key.width=unit(1,"cm"),
          legend.text=element_text(size = 14,face="bold"),
          legend.title = element_blank(),
          legend.text.align = 0,
          legend.position = c(0.7, 0.9),
          plot.title = element_blank(),
          axis.title.x = element_text(size = 18),
          axis.title.y.left = element_text(color = "black", size=18,face="bold"),
          axis.title.y.right = element_text(color = "#619CFF", size=18,face="bold"),
          axis.text.y.left = element_text(size = 18,color = "black",face="bold"),
          axis.text.y.right = element_text(size = 18,color = "#619CFF",face="bold"),
          axis.text.x = element_text(size = 18, face="bold"),
          plot.margin = margin(0.5, 0.5, 0, 0, "cm")) +
    labs(x = TeX("Outer scenario $S_\\tau$"))

# -- Density comparisons
g_star_plot = g_star(x = x_plot, x_out = x_out)/g_star_const(x_out = x_out)
f_mix_plot = f_beta(x = x_plot, x_out = x_out, beta = 1/length(x_out))
beta_fit = ls.fitting(S_out, r, T-tau, vol, num_sample = 1000, function(x) Payoff_IB(x, K, r, T, tau, IB_0), strat = TRUE)
f_beta_plot = f_beta(x = x_plot, x_out = x_out, beta = beta_fit)

df3 = data.frame(x = S_plot,
                 g_star = g_star_plot,
                 f_mix = f_mix_plot,
                 f_beta = f_beta_plot)
ggplot(df3, aes(x=x)) +
    geom_line(aes(y=g_star, linetype = "A"), linewidth=1.5) +
    geom_line(aes(y=f_mix, linetype = "B"), linewidth=1.5) +
    geom_line(aes(y=f_beta, linetype = "C"), linewidth=1.5) +
    scale_linetype_manual(labels = c(TeX("True optimal $q^*(S_T)$"),
                                     TeX("Equal mixture $\\bar{p}(S_T)$"),
                                     TeX("NNLS weighted mixture $q^*_{\\beta}(S_T)$")),
                          values=c("solid", "dashed", "dotted")) +
    labs(x = TeX("Terminal stock price, $S_T$"))+
    theme(legend.key.width=unit(2,"cm"),
          legend.text=element_text(size = 14,face="bold"),
          legend.title = element_blank(),
          legend.text.align = 0,
          legend.position = c(0.25, 0.7),
          plot.title = element_blank(),
          axis.title.x = element_text(size = 18,face="bold"),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 18, face="bold"),
          axis.text.y = element_text(size = 18, face="bold"),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))


#==========
result_MLR_1 = matrix(rep(0,N_Macro*N_Out), ncol = N_Macro, nrow = N_Out)
result_MLR_10 = matrix(rep(0,N_Macro*N_Out), ncol = N_Macro, nrow = N_Out)
# result_MLR_100 = matrix(rep(0,N_Macro*N_Out), ncol = N_Macro, nrow = N_Out)
for(i in 1:N_Macro){
    print(i)
    result_MLR_1[,i] = NGS_MLR(S_out, r, T-tau, vol, num_sample = 1e0*N_Out, function(x) Payoff_IB(x, K, r, T, tau, IB_0), strat = FALSE)$C_tau
    result_MLR_10[,i] = NGS_MLR(S_out, r, T-tau, vol, num_sample = 1e1*N_Out, function(x) Payoff_IB(x, K, r, T, tau, IB_0), strat = FALSE)$C_tau
    # result_MLR_100[,i] = NGS_MLR(S_out, r, T-tau, vol, num_sample = 1e2*N_Out, function(x) Payoff_IB(x, K, r, T, tau, IB_0), strat = FALSE)$C_tau
}
Qntl_MLR_1 = apply(result_MLR_1 , 1 , quantile , probs = c(0.025, 0.975) , na.rm = TRUE)
Qntl_MLR_10 = apply(result_MLR_10 , 1 , quantile , probs = c(0.025, 0.975) , na.rm = TRUE)
# Qntl_MLR_100 = apply(result_MLR_100 , 1 , quantile , probs = c(0.025, 0.975) , na.rm = TRUE)

result_FIT_1 = matrix(rep(0,N_Macro*N_Out), ncol = N_Macro, nrow = N_Out)
result_FIT_10 = matrix(rep(0,N_Macro*N_Out), ncol = N_Macro, nrow = N_Out)
# result_FIT_100 = matrix(rep(0,N_Macro*N_Out), ncol = N_Macro, nrow = N_Out)
for(i in 1:N_Macro){
    print(i)
    result_FIT_1[,i] = NGS_LS(S_out, r, T-tau, vol, num_sample = 1e0*N_Out, function(x) Payoff_IB(x, K, r, T, tau, IB_0),
                              p = 0.99, strat = TRUE)$C_tau
    result_FIT_10[,i] = NGS_LS(S_out, r, T-tau, vol, num_sample = 1e1*N_Out, function(x) Payoff_IB(x, K, r, T, tau, IB_0),
                               p = 0.99, strat = TRUE)$C_tau
    # result_FIT_100[,i] = NGS_LS(S_out, r, T-tau, vol, num_sample = 1e2*N_Out, function(x) Payoff_IB(x, K, r, T, tau, IB_0),
    #                             p = 0.99, strat = TRUE)$C_tau
}
Qntl_FIT_1 = apply(result_FIT_1 , 1 , quantile , probs = c(0.025, 0.975) , na.rm = TRUE)
Qntl_FIT_10 = apply(result_FIT_10 , 1 , quantile , probs = c(0.025, 0.975) , na.rm = TRUE)
# Qntl_FIT_100 = apply(result_FIT_100 , 1 , quantile , probs = c(0.025, 0.975) , na.rm = TRUE)

df_FITMLRSTD = data.frame(x = S_out,
                          mu_true = mu_true,
                          lower1FIT = Qntl_FIT_1[1,],
                          upper1FIT = Qntl_FIT_1[2,],
                          lower10FIT = Qntl_FIT_10[1,],
                          upper10FIT = Qntl_FIT_10[2,],
                          lower1MLR = Qntl_MLR_1[1,],
                          upper1MLR = Qntl_MLR_1[2,],
                          lower10MLR = Qntl_MLR_10[1,],
                          upper10MLR = Qntl_MLR_10[2,],
                          lower1000STD = Qntl_STD_1000[1,],
                          upper1000STD = Qntl_STD_1000[2,])
ggplot(df_FITMLRSTD, aes(x=x)) +
    geom_ribbon(aes(ymin = lower1MLR, ymax = upper1MLR, fill = "MLR"), alpha=0.7)+
    geom_ribbon(aes(ymin = lower1000STD, ymax = upper1000STD, fill = "STD"), alpha=0.7)+
    geom_line(aes(y=mu_true), linewidth=1.5) +
    scale_fill_manual(values=c("#6baed6", "#f768a1"),
                      labels = c(TeX("MLR est. with $\\Gamma_{MLR} = 10^3$"),
                                 TeX("STD est. with $\\Gamma_{STD} = 10^6$")))+
    theme(legend.background = element_rect(fill="transparent"),
          legend.key.width=unit(1,"cm"),
          legend.text=element_text(size = 14,face="bold"),
          legend.title = element_blank(),
          legend.text.align = 0,
          legend.position = c(0.7, 0.9),
          plot.title = element_blank(),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 18, face="bold"),
          axis.text.y = element_text(size = 18, face="bold"),
          plot.margin = margin(0.5, 0.5, 0, 0, "cm")) +
    scale_y_continuous(limits = c(-4, 4), breaks = c(-4,-3,-2,-1,0,1,2,3)) +
    labs(x = TeX("Outer scenario $S_\\tau$"))

ggplot(df_FITMLRSTD, aes(x=x)) +
    geom_ribbon(aes(ymin = lower1FIT, ymax = upper1FIT, fill = "OIS"), alpha=0.7)+
    geom_ribbon(aes(ymin = lower1000STD, ymax = upper1000STD, fill = "STD"), alpha=0.7)+
    geom_line(aes(y=mu_true), linewidth=1.5) +
    scale_fill_manual(values=c("#6baed6", "#f768a1"),
                      labels = c(TeX("OIS est. with $\\Gamma_{OIS} = 10^3$"),
                                 TeX("STD est. with $\\Gamma_{STD} = 10^6$")))+
    theme(legend.background = element_rect(fill="transparent"),
          legend.key.width=unit(1,"cm"),
          legend.text=element_text(size = 14,face="bold"),
          legend.title = element_blank(),
          legend.text.align = 0,
          legend.position = c(0.7, 0.9),
          plot.title = element_blank(),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 18, face="bold"),
          axis.text.y = element_text(size = 18, face="bold"),
          plot.margin = margin(0.5, 0.5, 0, 0, "cm"))+
    scale_y_continuous(limits = c(-4, 4), breaks = c(-4,-3,-2,-1,0,1,2,3)) +
    labs(x = TeX("Outer scenario $S_\\tau$"))
