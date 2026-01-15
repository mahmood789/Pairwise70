################################################################################
# QUICK VALIDATION: V2 vs V1 METHODS
# Compare improved methods against original
################################################################################

library(metafor)
library(data.table)

# Load both versions
source("C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/Advanced_Pooling_Methods.R")
source("C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/Advanced_Pooling_Methods_V2.R")

set.seed(42)

cat(strrep("=", 70), "\n")
cat("QUICK VALIDATION: V2 vs V1 METHODS\n")
cat(strrep("=", 70), "\n\n")

################################################################################
# TEST 1: BCG Data (Real Data)
################################################################################

cat("TEST 1: BCG Vaccine Trials (k=13)\n")
cat(strrep("-", 50), "\n")

dat <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg,
              data = dat.bcg)

cat("\nV1 Methods:\n")
v1_results <- compare_all_methods(dat$yi, dat$vi)
print(v1_results)

cat("\nV2 Methods:\n")
v2_results <- compare_all_methods_v2(dat$yi, dat$vi)
print(v2_results)

################################################################################
# TEST 2: Simulated Outlier Scenario
################################################################################

cat("\n\nTEST 2: Simulated Outlier (k=10, one 3-SD outlier)\n")
cat(strrep("-", 50), "\n")

yi_outlier <- c(0.5, 0.6, 0.4, 0.7, 0.3, 0.5, 0.55, 0.45, 0.6, 2.0)  # Last one is outlier
vi_outlier <- rep(0.05, 10)

cat("\nTrue effect ~ 0.5, outlier = 2.0\n")

cat("\nV1 Methods:\n")
v1_out <- compare_all_methods(yi_outlier, vi_outlier)
print(v1_out)

cat("\nV2 Methods:\n")
v2_out <- compare_all_methods_v2(yi_outlier, vi_outlier)
print(v2_out)

# SIT specific comparison
cat("\nSIT Comparison:\n")
sit_v1 <- sequential_influence_trimming(yi_outlier, vi_outlier)
sit_v2 <- sequential_influence_trimming_v2(yi_outlier, vi_outlier)
cat(sprintf("  V1: estimate=%.4f, n_trimmed=%d\n", sit_v1$estimate, sit_v1$n_trimmed))
cat(sprintf("  V2: estimate=%.4f, n_trimmed=%d\n", sit_v2$estimate, sit_v2$n_trimmed))

################################################################################
# TEST 3: Simulated Publication Bias
################################################################################

cat("\n\nTEST 3: Simulated Publication Bias\n")
cat(strrep("-", 50), "\n")

# Small studies have larger effects
n_studies <- 15
sei <- seq(0.20, 0.05, length.out = n_studies)
yi_bias <- 0.3 + (sei - mean(sei)) * 4 + rnorm(n_studies, 0, 0.08)
vi_bias <- sei^2

cat("True effect = 0.3, small studies inflated\n")

cat("\nV1 Methods:\n")
v1_bias <- compare_all_methods(yi_bias, vi_bias)
print(v1_bias)

cat("\nV2 Methods:\n")
v2_bias <- compare_all_methods_v2(yi_bias, vi_bias)
print(v2_bias)

# UBSF specific comparison
cat("\nUBSF Comparison:\n")
ubsf_v1 <- unified_bias_stability(yi_bias, vi_bias)
ubsf_v2 <- unified_bias_stability_v2(yi_bias, vi_bias)
cat(sprintf("  V1: estimate=%.4f, bias_adj=%.4f\n", ubsf_v1$estimate, ubsf_v1$bias_adjustment))
cat(sprintf("  V2: estimate=%.4f, bias_adj=%.4f, method=%s\n",
            ubsf_v2$estimate, ubsf_v2$bias_adjustment, ubsf_v2$bias_method))

################################################################################
# TEST 4: 100-Iteration Mini Simulation
################################################################################

cat("\n\nTEST 4: 100-Iteration Mini Simulation\n")
cat(strrep("-", 50), "\n")

N_SIM <- 100
true_effect <- 0.3
tau2 <- 0.05
k <- 10

results_v1 <- data.frame()
results_v2 <- data.frame()

cat("Running 100 simulations...\n")

for (i in 1:N_SIM) {
  # Generate data
  n_per_arm <- sample(30:150, k, replace = TRUE)
  vi <- 2 / n_per_arm
  sei <- sqrt(vi)
  theta_i <- rnorm(k, true_effect, sqrt(tau2))
  yi <- rnorm(k, theta_i, sei)

  # V1 EMA
  ema_v1 <- tryCatch({
    ensemble_meta_analysis(yi, vi)
  }, error = function(e) list(estimate = NA, ci_lb = NA, ci_ub = NA))

  # V2 EMA
  ema_v2 <- tryCatch({
    ensemble_meta_analysis_v2(yi, vi)
  }, error = function(e) list(estimate = NA, ci_lb = NA, ci_ub = NA))

  results_v1 <- rbind(results_v1, data.frame(
    iter = i,
    estimate = ema_v1$estimate,
    covers = (!is.na(ema_v1$ci_lb)) && (ema_v1$ci_lb <= true_effect) && (ema_v1$ci_ub >= true_effect)
  ))

  results_v2 <- rbind(results_v2, data.frame(
    iter = i,
    estimate = ema_v2$estimate,
    covers = (!is.na(ema_v2$ci_lb)) && (ema_v2$ci_lb <= true_effect) && (ema_v2$ci_ub >= true_effect)
  ))
}

cat("\n100-Iteration Results:\n")
cat(sprintf("  EMA V1: Mean Est=%.4f, Bias=%.4f, Coverage=%.1f%%\n",
            mean(results_v1$estimate, na.rm = TRUE),
            mean(results_v1$estimate - true_effect, na.rm = TRUE),
            mean(results_v1$covers, na.rm = TRUE) * 100))
cat(sprintf("  EMA V2: Mean Est=%.4f, Bias=%.4f, Coverage=%.1f%%\n",
            mean(results_v2$estimate, na.rm = TRUE),
            mean(results_v2$estimate - true_effect, na.rm = TRUE),
            mean(results_v2$covers, na.rm = TRUE) * 100))

################################################################################
# SUMMARY
################################################################################

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("SUMMARY OF V2 IMPROVEMENTS\n")
cat(strrep("=", 70), "\n\n")

cat("KEY CHANGES MADE:\n")
cat("1. T-distribution CIs: Uses df = k-2 for coverage improvement\n")
cat("2. MWM: Reduced stability_weight from 0.5 to 0.3\n")
cat("3. SIT: Uses studentized residuals, adaptive trimming\n")
cat("4. UBSF: PET-PEESE integration for bias correction\n")
cat("5. ARP: Inverse-variance weights, Rubin's rules variance\n")
cat("6. EMA: Adaptive weights based on method precision\n\n")

cat("Run Simulation_1000_Iterations_V2.R for full validation.\n")
cat(strrep("=", 70), "\n")
