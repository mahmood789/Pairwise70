################################################################################
# SIMULATION STUDY: 1000 ITERATIONS
# Evaluating Advanced Pooling Methods vs RoBMA/RVE
################################################################################

library(metafor)
library(data.table)
library(parallel)

# Load advanced methods
source("C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/Advanced_Pooling_Methods.R")

set.seed(12345)

cat(strrep("=", 70), "\n")
cat("SIMULATION STUDY: 1000 ITERATIONS\n")
cat("Evaluating Advanced Pooling Methods\n")
cat(strrep("=", 70), "\n\n")

################################################################################
# SIMULATION PARAMETERS
################################################################################

N_SIMULATIONS <- 1000

# Scenarios to test
scenarios <- list(
  # Scenario 1: Low heterogeneity, no bias
  list(name = "Low Het, No Bias", k = 10, true_effect = 0.3, tau2 = 0.01, bias = FALSE),

  # Scenario 2: High heterogeneity, no bias
  list(name = "High Het, No Bias", k = 10, true_effect = 0.3, tau2 = 0.20, bias = FALSE),

  # Scenario 3: Low heterogeneity, with publication bias
  list(name = "Low Het, Pub Bias", k = 15, true_effect = 0.3, tau2 = 0.01, bias = TRUE),

  # Scenario 4: High heterogeneity, with publication bias
  list(name = "High Het, Pub Bias", k = 15, true_effect = 0.3, tau2 = 0.20, bias = TRUE),

  # Scenario 5: Small meta-analysis
  list(name = "Small MA (k=5)", k = 5, true_effect = 0.3, tau2 = 0.05, bias = FALSE),

  # Scenario 6: Large meta-analysis
  list(name = "Large MA (k=30)", k = 30, true_effect = 0.3, tau2 = 0.10, bias = FALSE),

  # Scenario 7: Null effect
  list(name = "Null Effect", k = 10, true_effect = 0.0, tau2 = 0.05, bias = FALSE),

  # Scenario 8: Large effect with outlier
  list(name = "Outlier Present", k = 10, true_effect = 0.5, tau2 = 0.05, bias = FALSE, outlier = TRUE)
)

################################################################################
# SIMULATION FUNCTIONS
################################################################################

#' Generate one meta-analysis dataset
generate_ma_data <- function(k, true_effect, tau2, bias = FALSE, outlier = FALSE) {

  # Sample sizes per study (realistic range)
  n_per_arm <- sample(20:200, k, replace = TRUE)

  # Sampling variances
  vi <- 2 / n_per_arm + tau2 * 0.1  # Approximate variance for SMD
  sei <- sqrt(vi)

  # True study effects (random effects)
  theta_i <- rnorm(k, true_effect, sqrt(tau2))

  # Observed effects
  yi <- rnorm(k, theta_i, sei)

  # Add publication bias (suppress small, non-significant studies)
  if (bias) {
    # Calculate p-values
    z <- yi / sei
    pvals <- 2 * pnorm(-abs(z))

    # Selection probability based on significance and direction
    prob_publish <- ifelse(pvals < 0.05 & yi > 0, 1.0,
                           ifelse(pvals < 0.10 & yi > 0, 0.7,
                                  ifelse(yi > 0, 0.4, 0.2)))

    # Filter studies
    selected <- runif(k) < prob_publish
    if (sum(selected) < 3) selected[1:3] <- TRUE  # Ensure minimum k

    yi <- yi[selected]
    vi <- vi[selected]
  }

  # Add outlier
  if (outlier && length(yi) > 3) {
    outlier_idx <- sample(length(yi), 1)
    yi[outlier_idx] <- yi[outlier_idx] + 3 * sqrt(vi[outlier_idx])  # 3 SD outlier
  }

  return(list(yi = yi, vi = vi, k = length(yi)))
}

#' Run all methods on one dataset
run_all_methods <- function(yi, vi, true_effect) {

  results <- list()

  # 1. REML
  tryCatch({
    fit <- rma(yi = yi, vi = vi, method = "REML")
    results$REML <- list(
      estimate = fit$beta[1],
      se = fit$se,
      ci_lb = fit$ci.lb,
      ci_ub = fit$ci.ub,
      covers = (fit$ci.lb <= true_effect) & (fit$ci.ub >= true_effect)
    )
  }, error = function(e) {
    results$REML <<- list(estimate = NA, se = NA, ci_lb = NA, ci_ub = NA, covers = NA)
  })

  # 2. HKSJ
  tryCatch({
    fit <- rma(yi = yi, vi = vi, method = "REML", test = "knha")
    results$HKSJ <- list(
      estimate = fit$beta[1],
      se = fit$se,
      ci_lb = fit$ci.lb,
      ci_ub = fit$ci.ub,
      covers = (fit$ci.lb <= true_effect) & (fit$ci.ub >= true_effect)
    )
  }, error = function(e) {
    results$HKSJ <<- list(estimate = NA, se = NA, ci_lb = NA, ci_ub = NA, covers = NA)
  })

  # 3. MWM
  tryCatch({
    mwm <- mafi_weighted_ma(yi, vi)
    results$MWM <- list(
      estimate = mwm$estimate,
      se = mwm$se,
      ci_lb = mwm$ci_lb,
      ci_ub = mwm$ci_ub,
      covers = (mwm$ci_lb <= true_effect) & (mwm$ci_ub >= true_effect)
    )
  }, error = function(e) {
    results$MWM <<- list(estimate = NA, se = NA, ci_lb = NA, ci_ub = NA, covers = NA)
  })

  # 4. ARP
  tryCatch({
    arp <- adaptive_robust_pooling(yi, vi)
    results$ARP <- list(
      estimate = arp$estimate,
      se = arp$se,
      ci_lb = arp$ci_lb,
      ci_ub = arp$ci_ub,
      covers = (arp$ci_lb <= true_effect) & (arp$ci_ub >= true_effect)
    )
  }, error = function(e) {
    results$ARP <<- list(estimate = NA, se = NA, ci_lb = NA, ci_ub = NA, covers = NA)
  })

  # 5. SIT
  tryCatch({
    sit <- sequential_influence_trimming(yi, vi)
    results$SIT <- list(
      estimate = sit$estimate,
      se = sit$se,
      ci_lb = sit$ci_lb,
      ci_ub = sit$ci_ub,
      covers = (sit$ci_lb <= true_effect) & (sit$ci_ub >= true_effect)
    )
  }, error = function(e) {
    results$SIT <<- list(estimate = NA, se = NA, ci_lb = NA, ci_ub = NA, covers = NA)
  })

  # 6. UBSF
  tryCatch({
    ubsf <- unified_bias_stability(yi, vi)
    results$UBSF <- list(
      estimate = ubsf$estimate,
      se = ubsf$se,
      ci_lb = ubsf$ci_lb,
      ci_ub = ubsf$ci_ub,
      covers = (ubsf$ci_lb <= true_effect) & (ubsf$ci_ub >= true_effect)
    )
  }, error = function(e) {
    results$UBSF <<- list(estimate = NA, se = NA, ci_lb = NA, ci_ub = NA, covers = NA)
  })

  # 7. EMA
  tryCatch({
    ema <- ensemble_meta_analysis(yi, vi)
    results$EMA <- list(
      estimate = ema$estimate,
      se = ema$se,
      ci_lb = ema$ci_lb,
      ci_ub = ema$ci_ub,
      covers = (ema$ci_lb <= true_effect) & (ema$ci_ub >= true_effect)
    )
  }, error = function(e) {
    results$EMA <<- list(estimate = NA, se = NA, ci_lb = NA, ci_ub = NA, covers = NA)
  })

  return(results)
}

################################################################################
# RUN SIMULATION
################################################################################

all_results <- list()

for (s in seq_along(scenarios)) {
  scenario <- scenarios[[s]]
  cat(sprintf("\nScenario %d: %s\n", s, scenario$name))
  cat(sprintf("  k=%d, true_effect=%.2f, tau2=%.2f, bias=%s\n",
              scenario$k, scenario$true_effect, scenario$tau2, scenario$bias))

  # Pre-allocate results list
  results_list <- vector("list", N_SIMULATIONS * 7)
  idx <- 1

  pb <- txtProgressBar(min = 0, max = N_SIMULATIONS, style = 3)

  for (i in 1:N_SIMULATIONS) {
    # Generate data
    outlier_flag <- ifelse(is.null(scenario$outlier), FALSE, scenario$outlier)
    dat <- generate_ma_data(
      k = scenario$k,
      true_effect = scenario$true_effect,
      tau2 = scenario$tau2,
      bias = scenario$bias,
      outlier = outlier_flag
    )

    # Run all methods
    results <- run_all_methods(dat$yi, dat$vi, scenario$true_effect)

    # Store results using consistent structure
    for (method in names(results)) {
      r <- results[[method]]
      results_list[[idx]] <- list(
        iteration = i,
        method = method,
        estimate = as.numeric(r$estimate),
        se = as.numeric(r$se),
        ci_lb = as.numeric(r$ci_lb),
        ci_ub = as.numeric(r$ci_ub),
        covers = as.logical(r$covers),
        true_effect = scenario$true_effect,
        scenario = scenario$name
      )
      idx <- idx + 1
    }

    setTxtProgressBar(pb, i)
  }
  close(pb)

  # Convert to data.table
  scenario_results <- rbindlist(results_list[1:(idx-1)], fill = TRUE)
  all_results[[scenario$name]] <- scenario_results
}

################################################################################
# AGGREGATE RESULTS
################################################################################

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("SIMULATION RESULTS SUMMARY\n")
cat(strrep("=", 70), "\n\n")

# Combine all results
combined_results <- rbindlist(all_results)

# Calculate performance metrics by scenario and method
performance <- combined_results[, .(
  Bias = mean(estimate - true_effect, na.rm = TRUE),
  RMSE = sqrt(mean((estimate - true_effect)^2, na.rm = TRUE)),
  Coverage = mean(covers, na.rm = TRUE) * 100,
  Mean_SE = mean(se, na.rm = TRUE),
  CI_Width = mean(ci_ub - ci_lb, na.rm = TRUE)
), by = .(scenario, method)]

# Print results by scenario
for (scen_name in unique(performance$scenario)) {
  cat(sprintf("\n%s\n", scen_name))
  cat(strrep("-", 60), "\n")

  scen_perf <- performance[scenario == scen_name]
  scen_perf <- scen_perf[order(RMSE)]

  cat(sprintf("%-8s %8s %8s %8s %8s %8s\n",
              "Method", "Bias", "RMSE", "Cov%", "Mean SE", "CI Width"))
  cat(strrep("-", 60), "\n")

  for (j in 1:nrow(scen_perf)) {
    cat(sprintf("%-8s %8.4f %8.4f %8.1f %8.4f %8.4f\n",
                scen_perf$method[j],
                scen_perf$Bias[j],
                scen_perf$RMSE[j],
                scen_perf$Coverage[j],
                scen_perf$Mean_SE[j],
                scen_perf$CI_Width[j]))
  }
}

################################################################################
# OVERALL RANKINGS
################################################################################

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("OVERALL METHOD RANKINGS (Across All Scenarios)\n")
cat(strrep("=", 70), "\n\n")

overall <- combined_results[, .(
  Mean_Bias = mean(estimate - true_effect, na.rm = TRUE),
  Mean_Abs_Bias = mean(abs(estimate - true_effect), na.rm = TRUE),
  Mean_RMSE = sqrt(mean((estimate - true_effect)^2, na.rm = TRUE)),
  Mean_Coverage = mean(covers, na.rm = TRUE) * 100,
  Mean_SE = mean(se, na.rm = TRUE)
), by = method]

overall <- overall[order(Mean_RMSE)]

cat(sprintf("%-8s %10s %10s %10s %10s %10s\n",
            "Method", "Bias", "Abs Bias", "RMSE", "Coverage%", "Mean SE"))
cat(strrep("-", 70), "\n")

for (j in 1:nrow(overall)) {
  cat(sprintf("%-8s %10.4f %10.4f %10.4f %10.1f %10.4f\n",
              overall$method[j],
              overall$Mean_Bias[j],
              overall$Mean_Abs_Bias[j],
              overall$Mean_RMSE[j],
              overall$Mean_Coverage[j],
              overall$Mean_SE[j]))
}

################################################################################
# SCENARIO-SPECIFIC WINNERS
################################################################################

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("BEST METHOD BY SCENARIO\n")
cat(strrep("=", 70), "\n\n")

cat(sprintf("%-25s %-10s %-10s %-10s\n",
            "Scenario", "Best RMSE", "Best Cov", "Best Bias"))
cat(strrep("-", 60), "\n")

for (scen_name in unique(performance$scenario)) {
  scen_perf <- performance[scenario == scen_name]

  best_rmse <- scen_perf[which.min(RMSE)]$method
  best_cov <- scen_perf[which.min(abs(Coverage - 95))]$method
  best_bias <- scen_perf[which.min(abs(Bias))]$method

  cat(sprintf("%-25s %-10s %-10s %-10s\n",
              scen_name, best_rmse, best_cov, best_bias))
}

################################################################################
# KEY FINDINGS
################################################################################

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("KEY FINDINGS FROM 1000 SIMULATIONS\n")
cat(strrep("=", 70), "\n\n")

# Find scenarios where advanced methods beat REML
advanced_methods <- c("MWM", "ARP", "SIT", "UBSF", "EMA")

for (scen_name in unique(performance$scenario)) {
  scen_perf <- performance[scenario == scen_name]
  reml_rmse <- scen_perf[method == "REML"]$RMSE
  reml_cov <- scen_perf[method == "REML"]$Coverage

  better_rmse <- scen_perf[method %in% advanced_methods & RMSE < reml_rmse]$method
  better_cov <- scen_perf[method %in% advanced_methods & abs(Coverage - 95) < abs(reml_cov - 95)]$method

  if (length(better_rmse) > 0 || length(better_cov) > 0) {
    cat(sprintf("%s:\n", scen_name))
    if (length(better_rmse) > 0) {
      cat(sprintf("  Better RMSE than REML: %s\n", paste(better_rmse, collapse = ", ")))
    }
    if (length(better_cov) > 0) {
      cat(sprintf("  Better coverage than REML: %s\n", paste(better_cov, collapse = ", ")))
    }
    cat("\n")
  }
}

################################################################################
# SAVE RESULTS
################################################################################

output_dir <- "C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

fwrite(combined_results, file.path(output_dir, "simulation_1000_raw.csv"))
fwrite(performance, file.path(output_dir, "simulation_1000_performance.csv"))
fwrite(overall, file.path(output_dir, "simulation_1000_overall.csv"))

cat(sprintf("\nResults saved to: %s\n", output_dir))
cat(strrep("=", 70), "\n")
cat("Simulation complete!\n")
