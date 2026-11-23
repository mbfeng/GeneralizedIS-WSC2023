################################################################################
# figures_tables.R                                                             #
# -- Figures 1-3 and Table 1 of Chen & Feng (WSC 2023) from the summaries      #
# -- Styling follows the original plotting script (original_2023/plots.R):     #
#    ggplot2, LaTeX labels by latex2exp, 6 x 4 in (Figure 3: 10 x 4 in)        #
################################################################################

library(ggplot2)     # figures, as in the original code
library(latex2exp)   # LaTeX axis and legend labels

band_colours <- function() {
  # Fill colours of the original figures (light pink, pink, dark blue).
  return(c(light = "#fbb4b9", pink = "#f768a1", blue = "#08519c", outer = "#619CFF"))
}

bold_theme <- function(legend_position) {
  # Text sizes and legend placement shared by Figures 1b and 2.
  # Input(s): legend_position = c(x, y) inside the panel
  theme_settings <- theme(legend.background = element_rect(fill = "transparent"),
                          legend.key.width = unit(1, "cm"),
                          legend.text = element_text(size = 14, face = "bold"),
                          legend.title = element_blank(),
                          legend.text.align = 0,
                          legend.position = legend_position,
                          plot.title = element_blank(),
                          axis.title.x = element_text(size = 18),
                          axis.text.x = element_text(size = 18, face = "bold"),
                          axis.text.y = element_text(size = 18, face = "bold"),
                          plot.margin = margin(0.5, 0.5, 0, 0, "cm"))
  return(theme_settings)
}

# ===== Figure 1: problem setting ====================================================

plot_figure1a <- function(par, file) {
  # Figure 1a, drawn as published: exp(-r (T - tau)) P(S_T) - P_0 (the formula for g in
  # Section 5) on 10,000 equally spaced log prices over outer_mean +/- 3 outer_sd.
  # (simulation_output(), used by the estimators, subtracts P_0 exp(r tau) instead.)
  # Input(s): par = gis_parameters(); file = output PDF
  log_s <- seq(par$outer_mean - 3 * par$outer_sd, par$outer_mean + 3 * par$outer_sd,
               length.out = 1e4)
  carried_forward <- par$price_0 * exp(par$r * par$horizon)     # constant in simulation_output
  g_drawn <- simulation_output(log_s, par) + carried_forward - par$price_0
  data <- data.frame(S = exp(log_s), g = g_drawn)
  figure <- ggplot(data, aes(x = S)) +
    geom_line(aes(y = g), linewidth = 2) +
    labs(x = TeX("Terminal stock price, $S_T$"), y = TeX("Simulation model, $g(S_T)$")) +
    theme(plot.title = element_blank(),
          axis.title = element_text(size = 18, face = "bold"),
          axis.text = element_text(size = 18, face = "bold"),
          plot.margin = margin(0.5, 0.5, 0, 0, "cm"))
  ggsave(file, figure, width = 6, height = 4)
}

plot_figure1b <- function(summary, par, file) {
  # Figure 1b: true value mu(S_tau), 95% bands of STD with n = 100 and n = 1,000, and the
  # outer pdf pi(S_tau) (shifted and scaled onto the same axis: pdf / 0.01 - 5).
  # Input(s): summary = by_scenario summaries; par; file = output PDF
  colours <- band_colours()
  data <- summary
  data$outer_pdf <- dlnorm(summary$S_tau, meanlog = par$outer_mean, sdlog = par$outer_sd)
  data$outer_scaled <- data$outer_pdf / 1e-2 - 5
  figure <- ggplot(data, aes(x = S_tau)) +
    geom_ribbon(aes(ymin = STD_n100_q025, ymax = STD_n1000_q025, fill = "B"), alpha = 0.9) +
    geom_ribbon(aes(ymin = STD_n1000_q975, ymax = STD_n100_q975, fill = "B"), alpha = 0.9) +
    geom_ribbon(aes(ymin = STD_n1000_q025, ymax = STD_n1000_q975, fill = "C"), alpha = 0.9) +
    geom_line(aes(y = mu_true), linewidth = 1.5) +
    scale_fill_manual(values = unname(colours[c("light", "pink")]),
                      labels = c(TeX("$n=10^2, \\Gamma_{STD} = 10^5$"),
                                 TeX("$n=10^3, \\Gamma_{STD} = 10^6$"))) +
    geom_line(aes(y = outer_scaled), linewidth = 1, color = colours[["outer"]]) +
    geom_ribbon(aes(ymax = outer_scaled), ymin = -5, alpha = 0.3, fill = colours[["outer"]]) +
    scale_y_continuous(name = TeX("True value $\\mu(S_\\tau)$"), limits = c(-5, 2.5),
                       sec.axis = sec_axis(~ (. + 5) * 1e-2,
                                           name = TeX("Outer pdf $\\pi(S_\\tau)$"))) +
    bold_theme(c(0.7, 0.9)) +
    theme(axis.title.y.left = element_text(color = "black", size = 18, face = "bold"),
          axis.title.y.right = element_text(color = colours[["outer"]], size = 18,
                                            face = "bold"),
          axis.text.y.left = element_text(size = 18, color = "black", face = "bold"),
          axis.text.y.right = element_text(size = 18, color = colours[["outer"]],
                                           face = "bold")) +
    labs(x = TeX("Outer scenario $S_\\tau$"))
  ggsave(file, figure, width = 6, height = 4)
}

# ===== Figure 2: 95% confidence bands ================================================

plot_two_bands <- function(summary, first, first_label, file) {
  # Figure 2a (first = STD with n = 1,000) or 2b (first = MLR): 95% bands of the first
  # estimator and of OIS (Gamma = 10^3), and the true value.
  # Input(s): summary = by_scenario summaries; first = column prefix ("STD_n1000", "MLR");
  #           first_label = LaTeX legend label; file = output PDF
  colours <- band_colours()
  data <- data.frame(S_tau = summary$S_tau, mu_true = summary$mu_true,
                     first_lower = summary[[paste0(first, "_q025")]],
                     first_upper = summary[[paste0(first, "_q975")]],
                     ois_lower = summary$OIS_q025, ois_upper = summary$OIS_q975)
  figure <- ggplot(data, aes(x = S_tau)) +
    geom_ribbon(aes(ymin = first_lower, ymax = first_upper, fill = "A"), alpha = 0.7) +
    geom_ribbon(aes(ymin = ois_lower, ymax = ois_upper, fill = "B"), alpha = 0.7) +
    geom_line(aes(y = mu_true), linewidth = 1.5) +
    scale_fill_manual(values = unname(colours[c("pink", "blue")]),
                      labels = c(TeX(first_label),
                                 TeX("OIS est. with $\\Gamma_{OIS} = 10^3$"))) +
    scale_y_continuous(limits = c(-4, 3), breaks = -4:3) +       # as in the published figures
    bold_theme(c(0.7, 0.9)) +
    theme(axis.title.y = element_blank()) +
    labs(x = TeX("Outer scenario $S_\\tau$"))
  ggsave(file, figure, width = 6, height = 4)
}

# ===== Figure 3: sampling densities =================================================

plot_figure3 <- function(weights, par, file) {
  # Figure 3: true optimal q*, equal mixture p-bar and NNLS weighted mixture q_beta*, as
  # pdfs of log S_T plotted against S_T.
  # Input(s): weights = beta* (Figure 3 fit); par; file = output PDF
  log_s <- outer_scenarios(par)
  densities <- figure3_densities(log_s, log_s, weights, par)
  densities$S_T <- exp(densities$log_s_T)
  figure <- ggplot(densities, aes(x = S_T)) +
    geom_line(aes(y = optimal, linetype = "A"), linewidth = 1.5) +
    geom_line(aes(y = equal_mixture, linetype = "B"), linewidth = 1.5) +
    geom_line(aes(y = nnls_mixture, linetype = "C"), linewidth = 1.5) +
    scale_linetype_manual(labels = c(TeX("True optimal $q^*(S_T)$"),
                                     TeX("Equal mixture $\\bar{p}(S_T)$"),
                                     TeX("NNLS weighted mixture $q^*_{\\beta}(S_T)$")),
                          values = c("solid", "dashed", "dotted")) +
    labs(x = TeX("Terminal stock price, $S_T$")) +
    theme(legend.key.width = unit(2, "cm"),
          legend.text = element_text(size = 14, face = "bold"),
          legend.title = element_blank(),
          legend.text.align = 0,
          legend.position = c(0.25, 0.7),
          plot.title = element_blank(),
          axis.title.x = element_text(size = 18, face = "bold"),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 18, face = "bold"),
          axis.text.y = element_text(size = 18, face = "bold"),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
  ggsave(file, figure, width = 10, height = 4)
}

# ===== Table 1: AMSEs ================================================================

table1_amse <- function(summaries) {
  # Table 1: AMSE of each estimator and its ratio to the OIS AMSE.
  # Input(s): summaries = list(by_scenario, by_macro) of read_summaries()
  # Output(s): data frame, one row per estimator, with
  # -- amse: (1/m) sum_i of the sample variance (divisor K - 1) of the K estimates of
  #    mu(S_tau,i), the quantity printed in Table 1
  # -- amse_se: approximate standard error of the AMSE: sd of the per-macro average squared
  #    errors divided by sqrt(K) (exact for mse_with_bias; bias is negligible here)
  # -- mse_with_bias: (1/m) sum_i of the mean squared error (variance + squared bias)
  # -- reduction: amse / amse of OIS, with its delta-method standard error
  names <- estimator_names()
  budgets <- c(1e3, 1e3, 1e3, 1e4, 1e5, 1e6)
  n_macro <- nrow(summaries$by_macro)
  table <- data.frame(estimator = names, budget = budgets, amse = 0, amse_se = 0,
                      mse_with_bias = 0)
  for (k in seq_along(names)) {
    variances <- summaries$by_scenario[[paste0(names[k], "_var")]]
    squared_errors <- summaries$by_macro[[paste0(names[k], "_sq_error")]]
    table$amse[k] <- mean(variances)
    table$amse_se[k] <- sd(squared_errors) / sqrt(n_macro)
    table$mse_with_bias[k] <- mean(squared_errors)
  }
  table$reduction <- table$amse / table$amse[1]
  relative_se <- sqrt((table$amse_se / table$amse)^2 + (table$amse_se[1] / table$amse[1])^2)
  table$reduction_se <- table$reduction * relative_se
  table$reduction_se[1] <- 0
  return(table)
}

write_table1_latex <- function(table, file) {
  # Writes Table 1 as a LaTeX tabular fragment (same layout as the paper).
  # Input(s): table = table1_amse(); file = output .tex
  reduction_text <- c("1", sprintf("%.2f", table$reduction[2]),
                      format(round(table$reduction[3]), big.mark = "\\!,\\!"),
                      sprintf("%.0f", table$reduction[4]), sprintf("%.2f", table$reduction[5]),
                      sprintf("%.2f", table$reduction[6]))
  amse_text <- c(sprintf("%.4f", table$amse[1:2]), sprintf("%.2f", table$amse[3:5]),
                 sprintf("%.3f", table$amse[6]))
  lines <- c("\\begin{tabular}{c|cc|cccc}",
             "\\toprule",
             paste("Estimators for $\\mu(S_{\\tau,i})$& OIS & MLR &",
                   "\\multicolumn{4}{|c}{STD different simulation budgets}   \\\\"),
             paste("Simulation budget & $\\boxed{\\Gamma = 10^3}$ & $\\boxed{\\Gamma = 10^3}$ &",
                   "$\\boxed{\\Gamma = 10^3}$& $\\Gamma = 10^4$ & $\\Gamma = 10^5$ &",
                   "$\\Gamma = 10^6$   \\\\"),
             "\\midrule",
             sprintf(paste("AMSE reduction & $\\boxed{%s}$ & $\\boxed{%s}$ & $\\boxed{%s}$ &",
                           "$%s$ & $%s$ & $%s$  \\\\"),
                     reduction_text[1], reduction_text[2], reduction_text[3], reduction_text[4],
                     reduction_text[5], reduction_text[6]),
             sprintf("AMSE & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ & $%s$  \\\\", amse_text[1],
                     amse_text[2], amse_text[3], amse_text[4], amse_text[5], amse_text[6]),
             "\\bottomrule",
             "\\end{tabular}")
  writeLines(lines, file)
}
