################################################################################
# utilities.R                                                                  #
# -- Loading the library, command-line arguments, timing, and the summaries    #
#    of the macro replications that all tables and figures are drawn from      #
################################################################################

load_library <- function(code_folder) {
  # Sources all library files of the rebuilt code.
  # Input(s): code_folder = path of the Code folder
  for (file in c("parameters.R", "option_portfolio.R", "estimators.R")) {
    source(file.path(code_folder, "R", file))
  }
}

parse_arguments <- function(defaults) {
  # Reads "--name value" pairs from the command line; values keep the type of the default.
  # Input(s): defaults = named list of default values
  # Output(s): named list of settings
  args <- commandArgs(trailingOnly = TRUE)
  settings <- defaults
  k <- 1
  while (k < length(args)) {
    name <- sub("^--", "", args[k])
    if (!(name %in% names(defaults))) {
      stop(paste("unknown argument", args[k]))
    }
    if (is.numeric(defaults[[name]])) {
      settings[[name]] <- as.numeric(args[k + 1])
    } else {
      settings[[name]] <- args[k + 1]
    }
    k <- k + 2
  }
  return(settings)
}

seconds_since <- function(start) {
  # Elapsed wall-clock seconds since start (a Sys.time() value).
  return(as.numeric(difftime(Sys.time(), start, units = "secs")))
}

# ===== Summaries of the macro replications =========================================

estimator_names <- function() {
  # Short names of the six estimators of Table 1 (budget Gamma in brackets):
  # OIS (10^3), MLR (10^3), STD with n = 1, 10, 100, 1000 (10^3 to 10^6).
  return(c("OIS", "MLR", "STD_n1", "STD_n10", "STD_n100", "STD_n1000"))
}

summarize_runs <- function(estimates, log_s_tau, mu_true) {
  # Reduces the m-by-K matrices of estimates (one column per macro replication) to the
  # summaries that all tables and figures use.
  # Input(s): estimates = named list of m-by-K matrices (names = estimator_names());
  #           log_s_tau = m scenarios; mu_true = m true values
  # Output(s): list
  # -- by_scenario: one row per scenario; for each estimator the 2.5% and 97.5% quantiles
  #    of its K estimates (R's default quantile type 7, as for the paper), their mean and
  #    sample variance (divisor K - 1)
  # -- by_macro: one row per macro replication; for each estimator the average squared
  #    error over the m scenarios, mean_i (estimate_i - mu_i)^2
  by_scenario <- data.frame(scenario = seq_along(log_s_tau), S_tau = exp(log_s_tau),
                            mu_true = mu_true)
  n_macro <- ncol(estimates[[1]])
  by_macro <- data.frame(macro = 1:n_macro)
  for (name in names(estimates)) {
    values <- estimates[[name]]                                  # [i, k]: scenario i, run k
    quantiles <- apply(values, 1, quantile, probs = c(0.025, 0.975))
    by_scenario[[paste0(name, "_q025")]] <- quantiles[1, ]
    by_scenario[[paste0(name, "_q975")]] <- quantiles[2, ]
    by_scenario[[paste0(name, "_mean")]] <- rowMeans(values)
    by_scenario[[paste0(name, "_var")]] <- apply(values, 1, var)
    squared_errors <- (values - mu_true)^2                       # mu_true recycled by row
    by_macro[[paste0(name, "_sq_error")]] <- colMeans(squared_errors)
  }
  return(list(by_scenario = by_scenario, by_macro = by_macro))
}

write_summaries <- function(summaries, weights, folder) {
  # Writes the summaries of summarize_runs() and the NNLS weights of Figure 3.
  # Input(s): summaries; weights = beta* (length m); folder = output folder
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  write.csv(summaries$by_scenario, file.path(folder, "by_scenario.csv"), row.names = FALSE)
  write.csv(summaries$by_macro, file.path(folder, "by_macro.csv"), row.names = FALSE)
  used <- which(weights > 0)
  write.csv(data.frame(scenario = used, weight = weights[used]),
            file.path(folder, "figure3_nnls_weights.csv"), row.names = FALSE)
}

read_summaries <- function(folder, m) {
  # Reads the files written by write_summaries().
  # Input(s): folder; m = number of scenarios (to expand the sparse weights)
  # Output(s): list(by_scenario, by_macro, weights)
  by_scenario <- read.csv(file.path(folder, "by_scenario.csv"))
  by_macro <- read.csv(file.path(folder, "by_macro.csv"))
  sparse <- read.csv(file.path(folder, "figure3_nnls_weights.csv"))
  weights <- rep(0, m)
  weights[sparse$scenario] <- sparse$weight
  return(list(by_scenario = by_scenario, by_macro = by_macro, weights = weights))
}
