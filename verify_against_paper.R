################################################################################
# verify_against_paper.R                                                       #
# -- Compares Table 1 as computed from a set of summaries with the values      #
#    printed in the paper                                                      #
# -- Exact mode (summaries of the published runs): every entry must round to   #
#    the printed digits                                                        #
# -- Monte Carlo mode (a rerun with new random numbers): every entry must      #
#    agree within Monte Carlo error,                                           #
#      z = (x - v) / sqrt(se_x^2 + se_v^2 + h^2 / 3),                          #
#    where v is the printed value, se_v = se_x (the paper's run had the same   #
#    size) and h is half a unit in the last printed digit; |z| > 2.5 flagged   #
# -- Called by reproduce_paper.R                                               #
################################################################################

printed_table1 <- function() {
  # Table 1 as printed (text, to keep the printed precision), from the final .tex,
  # checked against the published PDF.
  return(data.frame(estimator = estimator_names(),
                    reduction = c("1", "2.03", "1114", "110", "11.06", "1.10"),
                    amse = c("0.0167", "0.0339", "18.59", "1.84", "0.18", "0.018")))
}

decimals_of <- function(text) {
  # Number of decimals of a printed number ("0.0167" -> 4, "1114" -> 0).
  parts <- strsplit(text, ".", fixed = TRUE)[[1]]
  return(if (length(parts) == 2) nchar(parts[2]) else 0)
}

verify_table1 <- function(table, mode) {
  # Input(s): table = table1_amse(); mode = "exact" or "monte_carlo"
  # Output(s): data frame of checks (also printed); stops if any exact check fails
  printed <- printed_table1()
  checks <- data.frame()
  for (k in seq_len(nrow(printed))) {
    for (quantity in c("amse", "reduction")) {
      if (quantity == "reduction" && k == 1) {
        next                                                     # OIS / OIS = 1 by definition
      }
      text <- printed[[quantity]][k]
      value <- as.numeric(text)
      decimals <- decimals_of(text)
      half_unit <- 0.5 * 10^(-decimals)
      computed <- table[[quantity]][k]
      se <- table[[paste0(quantity, "_se")]][k]
      if (mode == "exact") {
        status <- if (round(computed, decimals) == value) "ok" else "UNEXPECTED"
        z <- NA
      } else {
        z <- (computed - value) / sqrt(2 * se^2 + half_unit^2 / 3)
        status <- if (abs(z) <= 2.5) "ok" else "UNEXPECTED"
      }
      checks <- rbind(checks, data.frame(estimator = printed$estimator[k], quantity = quantity,
                                         printed = text, computed = signif(computed, 6),
                                         se = signif(se, 3), z = round(z, 2), status = status))
    }
  }
  print(checks, row.names = FALSE)
  n_flagged <- sum(checks$status != "ok")
  if (mode == "monte_carlo") {
    cat(sprintf("%d of %d checks flagged; expected by chance: %.2f\n", n_flagged, nrow(checks),
                0.0124 * nrow(checks)))
  } else {
    cat(sprintf("%d of %d entries differ from the printed digits\n", n_flagged, nrow(checks)))
    if (n_flagged > 0) {
      stop("exact verification failed")
    }
  }
  return(checks)
}
