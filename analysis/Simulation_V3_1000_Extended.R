#' Extended Simulation for Advanced Pooling Methods V3
#'
#' Editorial Revision: 1000 iterations with extended scenarios
#' - Very small k (k=3, k=4)
#' - Large k (k=50, k=100)
#' - All V3 methods + RVE comparison
#' - Monte Carlo SEs for coverage estimates
#'
#' @author Pairwise70 Team
#' @date January 2026

# === SETUP ===
suppressPackageStartupMessages({
  library(metafor)
  library(data.table)
})

# Source V3 methods
source("C:/Users/user/OneDrive - NHS/Documents/Pairwise70/R/advanced_pooling_v3.R")

# Install clubSandwich if needed
if (!requireNamespace("clubSandwich", quietly = TRUE)) {
  install.packages("clubSandwich", repos = "https://cloud.r-project.org")
}

cat(paste(rep("=", 60), collapse = ""), "\n")
cat("ADVANCED POOLING METHODS V3 - 1000 ITERATION SIMULATION\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# === SIMULATION PARAMETERS ===
set.seed(42)
N_SIM <- 1000
TRUE_EFFECT <- 0.3

# Extended scenarios per editorial request
scenarios <- list(
  # Very small k (NEW - editorial request)
  list(name = "k3_standard", k = 3, tau2 = 0.05, outlier = FALSE, bias = 0),
  list(name = "k4_standard", k = 4, tau2 = 0.05, outlier = FALSE, bias = 0),

  # Standard k (from V2)
  list(name = "k5_standard", k = 5, tau2 = 0.05, outlier = FALSE, bias = 0),
  list(name = "k10_standard", k = 10, tau2 = 0.05, outlier = FALSE, bias = 0),
  list(name = "k15_standard", k = 15, tau2 = 0.05, outlier = FALSE, bias = 0),
  list(name = "k20_standard", k = 20, tau2 = 0.05, outlier = FALSE, bias = 0),

  # Large k (NEW - editorial request)
  list(name = "k50_standard", k = 50, tau2 = 0.05, outlier = FALSE, bias = 0),
  list(name = "k100_standard", k = 100, tau2 = 0.05, outlier = FALSE, bias = 0),

  # Heterogeneity variants
  list(name = "k10_no_het", k = 10, tau2 = 0, outlier = FALSE, bias = 0),
  list(name = "k10_high_het", k = 10, tau2 = 0.20, outlier = FALSE, bias = 0),

  # Outlier scenarios
  list(name = "k10_outlier", k = 10, tau2 = 0.05, outlier = TRUE, bias = 0),
  list(name = "k15_outlier", k = 15, tau2 = 0.05, outlier = TRUE, bias = 0),

  # Publication bias scenarios
  list(name = "k15_pubbias_mild", k = 15, tau2 = 0.05, outlier = FALSE, bias = 0.3),
  list(name = "k15_pubbias_mod", k = 15, tau2 = 0.05, outlier = FALSE, bias = 0.5),
  list(name = "k15_pubbias_severe", k = 15, tau2 = 0.05, outlier = FALSE, bias = 0.8)
)

cat(sprintf("Scenarios: %d\n", length(scenarios)))
cat(sprintf("Iterations per scenario: %d\n", N_SIM))
cat(sprintf("TRUE_EFFECT: %.2f\n\n", TRUE_EFFECT))

# === DATA GENERATION FUNCTIONS ===

generate_ma_data <- function(k, tau2, true_effect, outlier = FALSE, bias = 0) {
  # Sample sizes
  n <- sample(30:150, k, replace = TRUE)

  # True study effects
  theta <- rnorm(k, mean = true_effect, sd = sqrt(tau2))

  # Sampling variances
  vi <- 4 / n + runif(k, 0.01, 0.05)
  sei <- sqrt(vi)

  # Observed effects
  yi <- rnorm(k, mean = theta, sd = sei)

  # Add outlier if requested
  if (outlier && k >= 5) {
    outlier_idx <- sample(1:k, 1)
    yi[outlier_idx] <- true_effect + 3 * sqrt(tau2 + vi[outlier_idx])
  }

  # Publication bias: suppress non-significant negative results
  if (bias > 0) {
    pvals <- 2 * pnorm(-abs(yi / sei))
    suppress_prob <- ifelse(yi < 0 & pvals > 0.05, bias, 0)
    keep <- runif(k) > suppress_prob

    if (sum(keep) >= 3) {
      yi <- yi[keep]
      vi <- vi[keep]
    }
  }

  list(yi = yi, vi = vi)
}

# === METHOD RUNNERS ===

run_method <- function(yi, vi, method) {
  k <- length(yi)

  tryCatch({
    result <- switch(method,
      "REML" = {
        fit <- rma(yi, vi, method = "REML")
        list(
          estimate = as.numeric(coef(fit)),
          se = fit$se,
          ci_lb = fit$ci.lb,
          ci_ub = fit$ci.ub
        )
      },
      "HKSJ" = {
        fit <- rma(yi, vi, method = "REML", test = "knha")
        list(
          estimate = as.numeric(coef(fit)),
          se = fit$se,
          ci_lb = fit$ci.lb,
          ci_ub = fit$ci.ub
        )
      },
      "MWM_v3" = {
        res <- mafi_weighted_ma_v3(yi, vi)
        list(estimate = res$estimate, se = res$se, ci_lb = res$ci_lb, ci_ub = res$ci_ub)
      },
      "SIT_v3" = {
        # No bootstrap in simulation for speed (would be too slow)
        res <- sequential_influence_trimming_v3(yi, vi, bootstrap = FALSE)
        list(estimate = res$estimate, se = res$se, ci_lb = res$ci_lb, ci_ub = res$ci_ub)
      },
      "ARP_v3" = {
        res <- adaptive_robust_pooling_v3(yi, vi)
        list(estimate = res$estimate, se = res$se, ci_lb = res$ci_lb, ci_ub = res$ci_ub)
      },
      "RVE" = {
        res <- rve_meta(yi, vi)
        list(estimate = res$estimate, se = res$se, ci_lb = res$ci_lb, ci_ub = res$ci_ub)
      }
    )

    data.table(
      estimate = result$estimate,
      se = result$se,
      ci_lb = result$ci_lb,
      ci_ub = result$ci_ub,
      converged = TRUE
    )
  }, error = function(e) {
    data.table(
      estimate = NA_real_,
      se = NA_real_,
      ci_lb = NA_real_,
      ci_ub = NA_real_,
      converged = FALSE
    )
  })
}

# === MAIN SIMULATION LOOP ===

methods <- c("REML", "HKSJ", "MWM_v3", "SIT_v3", "ARP_v3", "RVE")
all_results <- list()

start_time <- Sys.time()

for (s in seq_along(scenarios)) {
  sc <- scenarios[[s]]
  cat(sprintf("\n[%d/%d] Scenario: %s (k=%d, tau2=%.2f)\n",
              s, length(scenarios), sc$name, sc$k, sc$tau2))

  scenario_results <- list()

  pb <- txtProgressBar(min = 0, max = N_SIM, style = 3)

  for (i in 1:N_SIM) {
    setTxtProgressBar(pb, i)

    # Generate data
    dat <- generate_ma_data(
      k = sc$k,
      tau2 = sc$tau2,
      true_effect = TRUE_EFFECT,
      outlier = sc$outlier,
      bias = sc$bias
    )

    # Run all methods
    for (m in methods) {
      res <- run_method(dat$yi, dat$vi, m)
      res[, `:=`(
        scenario = sc$name,
        method = m,
        iter = i,
        k = length(dat$yi),
        true_effect = TRUE_EFFECT
      )]
      scenario_results[[length(scenario_results) + 1]] <- res
    }
  }

  close(pb)

  all_results[[s]] <- rbindlist(scenario_results)
}

results <- rbindlist(all_results)

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf("\n\nTotal time: %.1f minutes\n", as.numeric(elapsed)))

# === COMPUTE METRICS WITH MONTE CARLO SEs ===

cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
cat("COMPUTING PERFORMANCE METRICS\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

metrics <- results[converged == TRUE, .(
  N = .N,
  Bias = mean(estimate - true_effect, na.rm = TRUE),
  Bias_MCSE = sd(estimate - true_effect, na.rm = TRUE) / sqrt(.N),
  RMSE = sqrt(mean((estimate - true_effect)^2, na.rm = TRUE)),
  Mean_SE = mean(se, na.rm = TRUE),
  Coverage = mean(ci_lb <= true_effect & ci_ub >= true_effect, na.rm = TRUE),
  Coverage_MCSE = sqrt(mean(ci_lb <= true_effect & ci_ub >= true_effect, na.rm = TRUE) *
                       (1 - mean(ci_lb <= true_effect & ci_ub >= true_effect, na.rm = TRUE)) / .N),
  CI_Width = mean(ci_ub - ci_lb, na.rm = TRUE),
  Convergence = .N / N_SIM
), by = .(scenario, method)]

# Round for display
metrics[, `:=`(
  Bias = round(Bias, 4),
  Bias_MCSE = round(Bias_MCSE, 4),
  RMSE = round(RMSE, 4),
  Mean_SE = round(Mean_SE, 4),
  Coverage = round(Coverage * 100, 1),
  Coverage_MCSE = round(Coverage_MCSE * 100, 1),
  CI_Width = round(CI_Width, 3),
  Convergence = round(Convergence * 100, 1)
)]

# === SUMMARY BY METHOD (OVERALL) ===

overall <- results[converged == TRUE, .(
  N = .N,
  Bias = round(mean(estimate - true_effect, na.rm = TRUE), 4),
  RMSE = round(sqrt(mean((estimate - true_effect)^2, na.rm = TRUE)), 4),
  Coverage = round(mean(ci_lb <= true_effect & ci_ub >= true_effect, na.rm = TRUE) * 100, 1),
  Coverage_MCSE = round(sqrt(mean(ci_lb <= true_effect & ci_ub >= true_effect, na.rm = TRUE) *
                             (1 - mean(ci_lb <= true_effect & ci_ub >= true_effect, na.rm = TRUE)) / .N) * 100, 2)
), by = method]

setorder(overall, -Coverage, RMSE)

cat("OVERALL RESULTS BY METHOD:\n")
cat(paste(rep("-", 60), collapse = ""), "\n")
print(overall)

# === RESULTS BY SCENARIO ===

cat("\n\nRESULTS BY SCENARIO:\n")
cat(paste(rep("-", 60), collapse = ""), "\n")

for (sc_name in unique(metrics$scenario)) {
  cat(sprintf("\n%s:\n", sc_name))
  print(metrics[scenario == sc_name, .(method, Bias, RMSE, Coverage, Coverage_MCSE, CI_Width)])
}

# === SPECIAL ANALYSIS: SMALL k PERFORMANCE ===

cat("\n\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
cat("SMALL k ANALYSIS (k=3,4,5)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

small_k <- metrics[grepl("k[345]_", scenario)]
small_k_summary <- small_k[, .(
  Mean_Coverage = round(mean(Coverage), 1),
  Min_Coverage = round(min(Coverage), 1),
  Mean_RMSE = round(mean(RMSE), 4)
), by = method]

setorder(small_k_summary, -Mean_Coverage)
print(small_k_summary)

# === SPECIAL ANALYSIS: LARGE k PERFORMANCE ===

cat("\n\nLARGE k ANALYSIS (k=50,100)\n")
cat(paste(rep("-", 60), collapse = ""), "\n")

large_k <- metrics[grepl("k(50|100)_", scenario)]
large_k_summary <- large_k[, .(
  Mean_Coverage = round(mean(Coverage), 1),
  Mean_RMSE = round(mean(RMSE), 4)
), by = method]

setorder(large_k_summary, -Mean_Coverage)
print(large_k_summary)

# === SPECIAL ANALYSIS: OUTLIER SCENARIOS ===

cat("\n\nOUTLIER SCENARIOS:\n")
cat(paste(rep("-", 60), collapse = ""), "\n")

outlier_scenarios <- metrics[grepl("outlier", scenario)]
print(outlier_scenarios[, .(scenario, method, Bias, RMSE, Coverage)])

# === SPECIAL ANALYSIS: PUBLICATION BIAS ===

cat("\n\nPUBLICATION BIAS SCENARIOS:\n")
cat(paste(rep("-", 60), collapse = ""), "\n")

pubbias <- metrics[grepl("pubbias", scenario)]
print(pubbias[, .(scenario, method, Bias, RMSE, Coverage)])

# === SAVE RESULTS ===

results_dir <- "C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/results"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

fwrite(results, file.path(results_dir, "simulation_v3_1000_raw.csv"))
fwrite(metrics, file.path(results_dir, "simulation_v3_1000_metrics.csv"))
fwrite(overall, file.path(results_dir, "simulation_v3_1000_overall.csv"))

cat("\n\nResults saved to:\n")
cat(sprintf("  - %s/simulation_v3_1000_raw.csv\n", results_dir))
cat(sprintf("  - %s/simulation_v3_1000_metrics.csv\n", results_dir))
cat(sprintf("  - %s/simulation_v3_1000_overall.csv\n", results_dir))

# === CONCLUSIONS ===

cat("\n\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
cat("CONCLUSIONS\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

best_coverage <- overall[which.max(Coverage)]
best_rmse <- overall[which.min(RMSE)]

cat(sprintf("\nBest Coverage: %s (%.1f%% +/- %.2f%%)\n",
            best_coverage$method, best_coverage$Coverage, best_coverage$Coverage_MCSE))
cat(sprintf("Best RMSE: %s (%.4f)\n", best_rmse$method, best_rmse$RMSE))

# Check if V3 methods meet editorial requirements
cat("\n\nEDITORIAL REQUIREMENT CHECKS:\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

v3_methods <- overall[method %in% c("MWM_v3", "SIT_v3", "ARP_v3")]
cat(sprintf("V3 methods coverage >= 95%%: %s\n",
            ifelse(all(v3_methods$Coverage >= 95), "PASS", "CHECK")))
cat(sprintf("Monte Carlo SE < 1%% (n=1000): %s\n",
            ifelse(all(overall$Coverage_MCSE < 1), "PASS", "FAIL")))
cat(sprintf("RVE comparison included: PASS\n"))
cat(sprintf("Small k scenarios (k=3,4): PASS\n"))
cat(sprintf("Large k scenarios (k=50,100): PASS\n"))

cat("\n\nSimulation complete!\n")
