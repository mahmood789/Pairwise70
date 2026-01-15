################################################################################
# FAST 500-ITERATION SIMULATION: V2 METHODS
# Streamlined validation of Advanced Pooling Methods V2
################################################################################

library(metafor)
library(data.table)

# Load V2 methods
source("C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/Advanced_Pooling_Methods_V2.R")

set.seed(42)

cat(strrep("=", 70), "\n")
cat("500-ITERATION SIMULATION: V2 METHODS VALIDATION\n")
cat(strrep("=", 70), "\n\n")

################################################################################
# PARAMETERS
################################################################################

N_SIM <- 500
TRUE_EFFECT <- 0.3

# Results storage
all_results <- data.table()

################################################################################
# FAST METHOD RUNNER (parallel-ready structure)
################################################################################

run_all_v2 <- function(yi, vi) {
  k <- length(yi)

  results <- list()

  # REML
  tryCatch({
    fit <- rma(yi, vi, method = "REML")
    results$REML <- c(coef(fit), fit$se, fit$ci.lb, fit$ci.ub)
  }, error = function(e) results$REML <<- c(NA, NA, NA, NA))

  # HKSJ
  tryCatch({
    fit <- rma(yi, vi, method = "REML", test = "knha")
    results$HKSJ <- c(coef(fit), fit$se, fit$ci.lb, fit$ci.ub)
  }, error = function(e) results$HKSJ <<- c(NA, NA, NA, NA))

  # MWM_v2
  tryCatch({
    res <- mafi_weighted_ma_v2(yi, vi)
    results$MWM_v2 <- c(res$estimate, res$se, res$ci_lb, res$ci_ub)
  }, error = function(e) results$MWM_v2 <<- c(NA, NA, NA, NA))

  # ARP_v2
  tryCatch({
    res <- adaptive_robust_pooling_v2(yi, vi)
    results$ARP_v2 <- c(res$estimate, res$se, res$ci_lb, res$ci_ub)
  }, error = function(e) results$ARP_v2 <<- c(NA, NA, NA, NA))

  # SIT_v2
  tryCatch({
    res <- sequential_influence_trimming_v2(yi, vi)
    results$SIT_v2 <- c(res$estimate, res$se, res$ci_lb, res$ci_ub)
  }, error = function(e) results$SIT_v2 <<- c(NA, NA, NA, NA))

  # UBSF_v2
  tryCatch({
    res <- unified_bias_stability_v2(yi, vi)
    results$UBSF_v2 <- c(res$estimate, res$se, res$ci_lb, res$ci_ub)
  }, error = function(e) results$UBSF_v2 <<- c(NA, NA, NA, NA))

  # EMA_v2
  tryCatch({
    res <- ensemble_meta_analysis_v2(yi, vi)
    results$EMA_v2 <- c(res$estimate, res$se, res$ci_lb, res$ci_ub)
  }, error = function(e) results$EMA_v2 <<- c(NA, NA, NA, NA))

  results
}

################################################################################
# SCENARIO DEFINITIONS
################################################################################

scenarios <- list(
  # Standard scenarios
  list(name = "S1: k=10, tau2=0", k = 10, tau2 = 0, outlier = FALSE, pub_bias = "none"),
  list(name = "S2: k=10, tau2=0.05", k = 10, tau2 = 0.05, outlier = FALSE, pub_bias = "none"),
  list(name = "S3: k=10, tau2=0.15", k = 10, tau2 = 0.15, outlier = FALSE, pub_bias = "none"),
  list(name = "S4: k=20, tau2=0.05", k = 20, tau2 = 0.05, outlier = FALSE, pub_bias = "none"),
  list(name = "S5: k=5, tau2=0.05", k = 5, tau2 = 0.05, outlier = FALSE, pub_bias = "none"),
  # Outlier scenarios
  list(name = "S6: k=10 + outlier", k = 10, tau2 = 0.05, outlier = TRUE, pub_bias = "none"),
  list(name = "S7: k=15 + outlier", k = 15, tau2 = 0.05, outlier = TRUE, pub_bias = "none"),
  # Publication bias scenarios
  list(name = "S8: mild pub bias", k = 15, tau2 = 0.05, outlier = FALSE, pub_bias = "mild"),
  list(name = "S9: moderate pub bias", k = 15, tau2 = 0.05, outlier = FALSE, pub_bias = "moderate"),
  list(name = "S10: severe pub bias", k = 15, tau2 = 0.05, outlier = FALSE, pub_bias = "severe")
)

################################################################################
# RUN SIMULATION
################################################################################

for (s in seq_along(scenarios)) {
  scen <- scenarios[[s]]
  cat(sprintf("%s ... ", scen$name))
  flush.console()

  for (i in 1:N_SIM) {
    # Generate base data
    n_per_arm <- sample(30:150, scen$k, replace = TRUE)
    vi <- 2 / n_per_arm
    sei <- sqrt(vi)
    theta_i <- rnorm(scen$k, TRUE_EFFECT, sqrt(scen$tau2))
    yi <- rnorm(scen$k, theta_i, sei)

    # Add outlier if needed
    if (scen$outlier) {
      idx <- sample(1:scen$k, 1)
      yi[idx] <- TRUE_EFFECT + 3 * sqrt(scen$tau2 + vi[idx])
    }

    # Add publication bias if needed
    if (scen$pub_bias != "none") {
      strength <- switch(scen$pub_bias, mild = 0.3, moderate = 0.5, severe = 0.8)
      yi <- yi + strength * (sei - mean(sei)) / sd(sei) * 0.1
    }

    # Run methods
    res <- run_all_v2(yi, vi)

    # Store results
    for (m in names(res)) {
      all_results <- rbind(all_results, data.table(
        scenario = s,
        scenario_name = scen$name,
        k = scen$k,
        tau2 = scen$tau2,
        iter = i,
        method = m,
        estimate = res[[m]][1],
        se = res[[m]][2],
        ci_lb = res[[m]][3],
        ci_ub = res[[m]][4],
        true_effect = TRUE_EFFECT
      ))
    }
  }
  cat("done\n")
}

################################################################################
# PERFORMANCE METRICS
################################################################################

cat("\n")
cat(strrep("=", 70), "\n")
cat("PERFORMANCE METRICS\n")
cat(strrep("=", 70), "\n\n")

# Overall
overall <- all_results[, .(
  Bias = mean(estimate - true_effect, na.rm = TRUE),
  RMSE = sqrt(mean((estimate - true_effect)^2, na.rm = TRUE)),
  Coverage = mean(ci_lb <= true_effect & ci_ub >= true_effect, na.rm = TRUE) * 100,
  CI_Width = mean(ci_ub - ci_lb, na.rm = TRUE)
), by = method]

cat("OVERALL (5,000 simulations per method):\n")
print(overall[order(RMSE)])

# By scenario type
cat("\n\nSTANDARD SCENARIOS (S1-S5):\n")
std <- all_results[scenario <= 5, .(
  Bias = mean(estimate - true_effect, na.rm = TRUE),
  RMSE = sqrt(mean((estimate - true_effect)^2, na.rm = TRUE)),
  Coverage = mean(ci_lb <= true_effect & ci_ub >= true_effect, na.rm = TRUE) * 100
), by = method]
print(std[order(RMSE)])

cat("\n\nOUTLIER SCENARIOS (S6-S7):\n")
out <- all_results[scenario %in% c(6, 7), .(
  Bias = mean(estimate - true_effect, na.rm = TRUE),
  RMSE = sqrt(mean((estimate - true_effect)^2, na.rm = TRUE)),
  Coverage = mean(ci_lb <= true_effect & ci_ub >= true_effect, na.rm = TRUE) * 100
), by = method]
print(out[order(RMSE)])

cat("\n\nPUBLICATION BIAS SCENARIOS (S8-S10):\n")
pb <- all_results[scenario >= 8, .(
  Bias = mean(estimate - true_effect, na.rm = TRUE),
  RMSE = sqrt(mean((estimate - true_effect)^2, na.rm = TRUE)),
  Coverage = mean(ci_lb <= true_effect & ci_ub >= true_effect, na.rm = TRUE) * 100
), by = method]
print(pb[order(abs(Bias))])

################################################################################
# KEY COMPARISONS
################################################################################

cat("\n")
cat(strrep("=", 70), "\n")
cat("KEY COMPARISONS\n")
cat(strrep("=", 70), "\n\n")

reml_rmse <- overall[method == "REML", RMSE]
hksj_cov <- overall[method == "HKSJ", Coverage]

cat("RMSE vs REML:\n")
for (m in c("HKSJ", "MWM_v2", "ARP_v2", "SIT_v2", "UBSF_v2", "EMA_v2")) {
  mr <- overall[method == m, RMSE]
  pct <- (reml_rmse - mr) / reml_rmse * 100
  cat(sprintf("  %s: %+.1f%% (%s)\n", m, pct, ifelse(pct > 0, "better", "worse")))
}

cat("\nCoverage (target=95%):\n")
for (m in overall$method) {
  mc <- overall[method == m, Coverage]
  cat(sprintf("  %s: %.1f%% (deviation: %+.1f%%)\n", m, mc, mc - 95))
}

cat("\nBest method for each scenario type:\n")
cat(sprintf("  Standard: %s (RMSE=%.4f)\n", std[which.min(RMSE), method], min(std$RMSE)))
cat(sprintf("  Outliers: %s (RMSE=%.4f)\n", out[which.min(RMSE), method], min(out$RMSE)))
cat(sprintf("  Pub Bias: %s (Bias=%.4f)\n", pb[which.min(abs(Bias)), method], pb[which.min(abs(Bias)), Bias]))

################################################################################
# SAVE RESULTS
################################################################################

output_dir <- "C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/results"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

fwrite(all_results, file.path(output_dir, "simulation_v2_500_raw.csv"))
fwrite(overall, file.path(output_dir, "simulation_v2_500_overall.csv"))

cat("\n")
cat(strrep("=", 70), "\n")
cat("COMPLETE - Results saved\n")
cat(strrep("=", 70), "\n")
