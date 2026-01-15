################################################################################
# FULL 1000-ITERATION SIMULATION: V2 METHODS
# Comprehensive validation of Advanced Pooling Methods V2
################################################################################

library(metafor)
library(data.table)

# Load V2 methods
source("C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/Advanced_Pooling_Methods_V2.R")

set.seed(42)

cat(strrep("=", 70), "\n")
cat("1000-ITERATION SIMULATION: V2 METHODS VALIDATION\n")
cat(strrep("=", 70), "\n\n")

################################################################################
# SIMULATION PARAMETERS
################################################################################

N_SIM <- 1000
K_STUDIES <- c(5, 10, 20)  # Small, medium, large
TAU2_LEVELS <- c(0, 0.05, 0.15)  # No, moderate, high heterogeneity
TRUE_EFFECT <- 0.3

# Results storage
all_results <- data.table()

################################################################################
# HELPER FUNCTIONS
################################################################################

run_methods_v2 <- function(yi, vi) {
  k <- length(yi)

  # Standard methods
  reml <- tryCatch({
    fit <- rma(yi, vi, method = "REML")
    list(est = coef(fit), se = fit$se, ci_lb = fit$ci.lb, ci_ub = fit$ci.ub)
  }, error = function(e) list(est = NA, se = NA, ci_lb = NA, ci_ub = NA))

  hksj <- tryCatch({
    fit <- rma(yi, vi, method = "REML", test = "knha")
    list(est = coef(fit), se = fit$se, ci_lb = fit$ci.lb, ci_ub = fit$ci.ub)
  }, error = function(e) list(est = NA, se = NA, ci_lb = NA, ci_ub = NA))

  # V2 methods
  mwm <- tryCatch({
    res <- mafi_weighted_ma_v2(yi, vi)
    list(est = res$estimate, se = res$se, ci_lb = res$ci_lb, ci_ub = res$ci_ub)
  }, error = function(e) list(est = NA, se = NA, ci_lb = NA, ci_ub = NA))

  arp <- tryCatch({
    res <- adaptive_robust_pooling_v2(yi, vi)
    list(est = res$estimate, se = res$se, ci_lb = res$ci_lb, ci_ub = res$ci_ub)
  }, error = function(e) list(est = NA, se = NA, ci_lb = NA, ci_ub = NA))

  sit <- tryCatch({
    res <- sequential_influence_trimming_v2(yi, vi)
    list(est = res$estimate, se = res$se, ci_lb = res$ci_lb, ci_ub = res$ci_ub)
  }, error = function(e) list(est = NA, se = NA, ci_lb = NA, ci_ub = NA))

  ubsf <- tryCatch({
    res <- unified_bias_stability_v2(yi, vi)
    list(est = res$estimate, se = res$se, ci_lb = res$ci_lb, ci_ub = res$ci_ub)
  }, error = function(e) list(est = NA, se = NA, ci_lb = NA, ci_ub = NA))

  ema <- tryCatch({
    res <- ensemble_meta_analysis_v2(yi, vi)
    list(est = res$estimate, se = res$se, ci_lb = res$ci_lb, ci_ub = res$ci_ub)
  }, error = function(e) list(est = NA, se = NA, ci_lb = NA, ci_ub = NA))

  list(
    REML = reml, HKSJ = hksj, MWM_v2 = mwm, ARP_v2 = arp,
    SIT_v2 = sit, UBSF_v2 = ubsf, EMA_v2 = ema
  )
}

################################################################################
# SCENARIO 1-3: STANDARD (varying k, tau2)
################################################################################

cat("SCENARIOS 1-9: Standard Random-Effects (k x tau2 combinations)\n")
cat(strrep("-", 70), "\n")

scenario_num <- 0
for (k in K_STUDIES) {
  for (tau2 in TAU2_LEVELS) {
    scenario_num <- scenario_num + 1

    cat(sprintf("Scenario %d: k=%d, tau2=%.2f ... ", scenario_num, k, tau2))

    for (i in 1:N_SIM) {
      # Generate data
      n_per_arm <- sample(30:150, k, replace = TRUE)
      vi <- 2 / n_per_arm
      sei <- sqrt(vi)
      theta_i <- rnorm(k, TRUE_EFFECT, sqrt(tau2))
      yi <- rnorm(k, theta_i, sei)

      # Run methods
      res <- run_methods_v2(yi, vi)

      # Store results
      for (method in names(res)) {
        all_results <- rbind(all_results, data.table(
          scenario = scenario_num,
          scenario_type = "Standard",
          k = k,
          tau2 = tau2,
          outlier = FALSE,
          pub_bias = FALSE,
          iter = i,
          method = method,
          estimate = res[[method]]$est,
          se = res[[method]]$se,
          ci_lb = res[[method]]$ci_lb,
          ci_ub = res[[method]]$ci_ub,
          true_effect = TRUE_EFFECT
        ))
      }
    }
    cat("done\n")
  }
}

################################################################################
# SCENARIO 10-12: OUTLIER SCENARIOS
################################################################################

cat("\nSCENARIOS 10-12: Outlier Contamination\n")
cat(strrep("-", 70), "\n")

for (k in c(10, 15, 20)) {
  scenario_num <- scenario_num + 1

  cat(sprintf("Scenario %d: k=%d with 1 outlier (3 SD) ... ", scenario_num, k))

  for (i in 1:N_SIM) {
    # Generate data with one outlier
    n_per_arm <- sample(30:150, k, replace = TRUE)
    vi <- 2 / n_per_arm
    sei <- sqrt(vi)
    theta_i <- rnorm(k, TRUE_EFFECT, sqrt(0.05))
    yi <- rnorm(k, theta_i, sei)

    # Replace one with outlier (3 SD away)
    outlier_idx <- sample(1:k, 1)
    yi[outlier_idx] <- TRUE_EFFECT + 3 * sqrt(0.05 + vi[outlier_idx])

    # Run methods
    res <- run_methods_v2(yi, vi)

    # Store results
    for (method in names(res)) {
      all_results <- rbind(all_results, data.table(
        scenario = scenario_num,
        scenario_type = "Outlier",
        k = k,
        tau2 = 0.05,
        outlier = TRUE,
        pub_bias = FALSE,
        iter = i,
        method = method,
        estimate = res[[method]]$est,
        se = res[[method]]$se,
        ci_lb = res[[method]]$ci_lb,
        ci_ub = res[[method]]$ci_ub,
        true_effect = TRUE_EFFECT
      ))
    }
  }
  cat("done\n")
}

################################################################################
# SCENARIO 13-15: PUBLICATION BIAS
################################################################################

cat("\nSCENARIOS 13-15: Publication Bias\n")
cat(strrep("-", 70), "\n")

for (severity in c("mild", "moderate", "severe")) {
  scenario_num <- scenario_num + 1
  k <- 15

  cat(sprintf("Scenario %d: %s publication bias ... ", scenario_num, severity))

  # Bias severity parameters
  bias_strength <- switch(severity, mild = 0.3, moderate = 0.5, severe = 0.8)

  for (i in 1:N_SIM) {
    # Generate data with publication bias (small-study effect)
    n_per_arm <- sample(30:150, k, replace = TRUE)
    vi <- 2 / n_per_arm
    sei <- sqrt(vi)
    theta_i <- rnorm(k, TRUE_EFFECT, sqrt(0.05))
    yi <- rnorm(k, theta_i, sei)

    # Induce small-study effect (smaller studies have inflated effects)
    yi <- yi + bias_strength * (sei - mean(sei)) / sd(sei) * 0.1

    # Run methods
    res <- run_methods_v2(yi, vi)

    # Store results
    for (method in names(res)) {
      all_results <- rbind(all_results, data.table(
        scenario = scenario_num,
        scenario_type = paste0("PubBias_", severity),
        k = k,
        tau2 = 0.05,
        outlier = FALSE,
        pub_bias = TRUE,
        iter = i,
        method = method,
        estimate = res[[method]]$est,
        se = res[[method]]$se,
        ci_lb = res[[method]]$ci_lb,
        ci_ub = res[[method]]$ci_ub,
        true_effect = TRUE_EFFECT
      ))
    }
  }
  cat("done\n")
}

################################################################################
# CALCULATE PERFORMANCE METRICS
################################################################################

cat("\n")
cat(strrep("=", 70), "\n")
cat("CALCULATING PERFORMANCE METRICS\n")
cat(strrep("=", 70), "\n\n")

# Calculate metrics by scenario and method
performance <- all_results[, .(
  Bias = mean(estimate - true_effect, na.rm = TRUE),
  RMSE = sqrt(mean((estimate - true_effect)^2, na.rm = TRUE)),
  Coverage = mean(ci_lb <= true_effect & ci_ub >= true_effect, na.rm = TRUE) * 100,
  CI_Width = mean(ci_ub - ci_lb, na.rm = TRUE),
  N_Valid = sum(!is.na(estimate))
), by = .(scenario, scenario_type, k, tau2, method)]

# Overall metrics
overall <- all_results[, .(
  Bias = mean(estimate - true_effect, na.rm = TRUE),
  RMSE = sqrt(mean((estimate - true_effect)^2, na.rm = TRUE)),
  Coverage = mean(ci_lb <= true_effect & ci_ub >= true_effect, na.rm = TRUE) * 100,
  CI_Width = mean(ci_ub - ci_lb, na.rm = TRUE),
  N_Valid = sum(!is.na(estimate))
), by = .(method)]

################################################################################
# DISPLAY RESULTS
################################################################################

cat("OVERALL PERFORMANCE (All 15,000 simulations per method)\n")
cat(strrep("-", 70), "\n")
print(overall[order(RMSE)])

cat("\n\nBY SCENARIO TYPE\n")
cat(strrep("-", 70), "\n")

# Standard scenarios
cat("\nStandard Random-Effects:\n")
std_perf <- performance[scenario_type == "Standard", .(
  Bias = mean(Bias),
  RMSE = mean(RMSE),
  Coverage = mean(Coverage),
  CI_Width = mean(CI_Width)
), by = method]
print(std_perf[order(RMSE)])

# Outlier scenarios
cat("\nOutlier Contamination:\n")
out_perf <- performance[scenario_type == "Outlier", .(
  Bias = mean(Bias),
  RMSE = mean(RMSE),
  Coverage = mean(Coverage),
  CI_Width = mean(CI_Width)
), by = method]
print(out_perf[order(RMSE)])

# Publication bias scenarios
cat("\nPublication Bias (all severities):\n")
pb_perf <- performance[grepl("PubBias", scenario_type), .(
  Bias = mean(Bias),
  RMSE = mean(RMSE),
  Coverage = mean(Coverage),
  CI_Width = mean(CI_Width)
), by = method]
print(pb_perf[order(abs(Bias))])

################################################################################
# COMPARISON: V2 vs STANDARD METHODS
################################################################################

cat("\n")
cat(strrep("=", 70), "\n")
cat("V2 IMPROVEMENTS vs STANDARD METHODS\n")
cat(strrep("=", 70), "\n\n")

reml_rmse <- overall[method == "REML", RMSE]
reml_cov <- overall[method == "REML", Coverage]
hksj_cov <- overall[method == "HKSJ", Coverage]

cat("RMSE Improvement over REML:\n")
for (m in c("MWM_v2", "ARP_v2", "SIT_v2", "UBSF_v2", "EMA_v2")) {
  method_rmse <- overall[method == m, RMSE]
  improvement <- (reml_rmse - method_rmse) / reml_rmse * 100
  cat(sprintf("  %s: %.1f%% %s\n", m, abs(improvement),
              ifelse(improvement > 0, "better", "worse")))
}

cat("\nCoverage vs 95% Nominal:\n")
for (m in c("REML", "HKSJ", "MWM_v2", "ARP_v2", "SIT_v2", "UBSF_v2", "EMA_v2")) {
  method_cov <- overall[method == m, Coverage]
  cat(sprintf("  %s: %.1f%% (deviation: %.1f%%)\n", m, method_cov,
              method_cov - 95))
}

################################################################################
# SCENARIO-SPECIFIC WINNERS
################################################################################

cat("\n")
cat(strrep("=", 70), "\n")
cat("BEST METHOD BY SCENARIO\n")
cat(strrep("=", 70), "\n\n")

# For each scenario, find best RMSE
for (s in unique(performance$scenario)) {
  scen_data <- performance[scenario == s]
  best_rmse <- scen_data[which.min(RMSE)]
  best_cov <- scen_data[which.min(abs(Coverage - 95))]

  cat(sprintf("Scenario %d (%s): Best RMSE=%s (%.4f), Best Coverage=%s (%.1f%%)\n",
              s, unique(scen_data$scenario_type),
              best_rmse$method, best_rmse$RMSE,
              best_cov$method, best_cov$Coverage))
}

################################################################################
# SAVE RESULTS
################################################################################

output_dir <- "C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/results"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

fwrite(all_results, file.path(output_dir, "simulation_v2_1000_raw.csv"))
fwrite(performance, file.path(output_dir, "simulation_v2_1000_performance.csv"))
fwrite(overall, file.path(output_dir, "simulation_v2_1000_overall.csv"))

cat("\n")
cat(strrep("=", 70), "\n")
cat("SIMULATION COMPLETE\n")
cat(strrep("=", 70), "\n")
cat(sprintf("Results saved to: %s\n", output_dir))
cat("  - simulation_v2_1000_raw.csv (all iterations)\n")
cat("  - simulation_v2_1000_performance.csv (by scenario)\n")
cat("  - simulation_v2_1000_overall.csv (overall summary)\n")
