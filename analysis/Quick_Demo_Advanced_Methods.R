################################################################################
# QUICK DEMONSTRATION OF ADVANCED POOLING METHODS
# Testing on selected datasets for rapid validation
################################################################################

library(metafor)

# Load advanced methods
source("C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/Advanced_Pooling_Methods.R")

cat(strrep("=", 70), "\n")
cat("QUICK DEMONSTRATION: ADVANCED POOLING METHODS\n")
cat(strrep("=", 70), "\n\n")

################################################################################
# DEMO 1: BCG Vaccine Trial (Classic Meta-Analysis)
################################################################################

cat("DEMO 1: BCG Vaccine Trials (k=13)\n")
cat(strrep("-", 50), "\n\n")

dat <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg,
              data = dat.bcg)

# Run all methods
comparison <- compare_all_methods(dat$yi, dat$vi)
print(comparison)

cat("\n")

# Detailed MWM output
mwm <- mafi_weighted_ma(dat$yi, dat$vi)
cat(sprintf("MWM Stability Scores: %s\n",
            paste(round(mwm$study_stability_scores, 2), collapse=", ")))
cat(sprintf("MWM Adjustment from REML: %.4f\n\n", mwm$adjustment))

################################################################################
# DEMO 2: Simulated High Heterogeneity Data
################################################################################

cat("\nDEMO 2: Simulated High Heterogeneity (k=10, I2~75%%)\n")
cat(strrep("-", 50), "\n\n")

set.seed(42)
k <- 10
yi_sim <- c(0.5, 0.8, 0.2, 1.2, -0.3, 0.6, 0.4, 1.0, 0.7, 2.5)  # High variability
vi_sim <- c(0.05, 0.08, 0.04, 0.06, 0.07, 0.05, 0.09, 0.04, 0.06, 0.10)

comparison2 <- compare_all_methods(yi_sim, vi_sim)
print(comparison2)

# Highlight influence trimming
sit <- sequential_influence_trimming(yi_sim, vi_sim)
cat(sprintf("\nSIT trimmed %d influential studies\n", sit$n_trimmed))
cat(sprintf("Original estimate: %.4f, SIT estimate: %.4f\n",
            sit$estimate_original, sit$estimate))

################################################################################
# DEMO 3: Simulated Publication Bias
################################################################################

cat("\nDEMO 3: Simulated Publication Bias (small studies have larger effects)\n")
cat(strrep("-", 50), "\n\n")

set.seed(123)
# Simulate funnel plot asymmetry
n_studies <- 12
true_effect <- 0.3
sei <- seq(0.15, 0.05, length.out = n_studies)
yi_bias <- true_effect + (sei - mean(sei)) * 3 + rnorm(n_studies, 0, 0.1)
vi_bias <- sei^2

comparison3 <- compare_all_methods(yi_bias, vi_bias)
print(comparison3)

# UBSF details
ubsf <- unified_bias_stability(yi_bias, vi_bias)
cat(sprintf("\nUBSF Diagnostics:\n"))
cat(sprintf("  Bias detected: %s (Egger p=%.4f)\n", ubsf$bias_detected, ubsf$bias_pval))
cat(sprintf("  Bias adjustment: %.4f\n", ubsf$bias_adjustment))
cat(sprintf("  Stability adjustment: %.4f\n", ubsf$stability_adjustment))

################################################################################
# DEMO 4: Test on Pairwise70 Dataset
################################################################################

cat("\n\nDEMO 4: Real Cochrane Dataset (Pairwise70)\n")
cat(strrep("-", 50), "\n\n")

# Load Pairwise70
if (require("Pairwise70", quietly = TRUE) ||
    tryCatch({devtools::load_all("C:/Users/user/OneDrive - NHS/Documents/Pairwise70"); TRUE},
             error = function(e) FALSE)) {

  # Get first available dataset
  data_list <- data(package = "Pairwise70")$results[, "Item"]

  # Try to find a suitable dataset
  for (dataset_name in data_list[1:20]) {
    tryCatch({
      data(list = dataset_name, package = "Pairwise70", envir = environment())
      d <- get(dataset_name)

      # Check for binary outcome data
      if (all(c("Experimental.cases", "Experimental.N",
                "Control.cases", "Control.N") %in% names(d))) {

        # Calculate effect sizes
        es <- escalc(measure = "OR",
                     ai = d$Experimental.cases,
                     bi = d$Experimental.N - d$Experimental.cases,
                     ci = d$Control.cases,
                     di = d$Control.N - d$Control.cases)

        yi <- es$yi[!is.na(es$yi) & !is.na(es$vi) & es$vi > 0]
        vi <- es$vi[!is.na(es$yi) & !is.na(es$vi) & es$vi > 0]

        if (length(yi) >= 5) {
          cat(sprintf("Dataset: %s (k=%d)\n\n", dataset_name, length(yi)))

          comparison4 <- compare_all_methods(yi, vi)
          print(comparison4)

          # EMA details
          ema <- ensemble_meta_analysis(yi, vi)
          cat(sprintf("\nEMA Method Agreement: %.1f%%\n", ema$agreement * 100))
          cat(sprintf("Methods contributing: %d\n", ema$n_methods))

          break
        }
      }
    }, error = function(e) NULL)
  }
} else {
  cat("Pairwise70 package not available for this demo.\n")
}

################################################################################
# SUMMARY
################################################################################

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("SUMMARY: ADVANTAGES OVER RoBMA AND RVE\n")
cat(strrep("=", 70), "\n\n")

cat("1. MAFI-WEIGHTED META-ANALYSIS (MWM)\n")
cat("   - Directly addresses conclusion STABILITY (RoBMA/RVE don't)\n")
cat("   - Downweights studies that could flip results\n")
cat("   - Computationally instant (vs RoBMA's MCMC)\n\n")

cat("2. ADAPTIVE ROBUST POOLING (ARP)\n")
cat("   - Combines REML, DL, PM, HKSJ with optimal weights\n")
cat("   - Adapts to heterogeneity level automatically\n")
cat("   - Accounts for estimator uncertainty\n\n")

cat("3. SEQUENTIAL INFLUENCE TRIMMING (SIT)\n")
cat("   - Addresses influential outliers directly\n")
cat("   - Soft trimming (preserves all studies via weights)\n")
cat("   - Transparent adjustment process\n\n")

cat("4. UNIFIED BIAS-STABILITY FRAMEWORK (UBSF)\n")
cat("   - Integrates bias AND stability correction\n")
cat("   - Provides decomposed adjustments\n")
cat("   - More interpretable than black-box Bayesian\n\n")

cat("5. ENSEMBLE META-ANALYSIS (EMA)\n")
cat("   - Combines all methods with uncertainty quantification\n")
cat("   - Most robust overall performance\n")
cat("   - Reports method agreement metric\n\n")

cat("All methods are implemented and ready for use!\n")
cat(strrep("=", 70), "\n")
