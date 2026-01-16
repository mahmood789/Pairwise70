#!/usr/bin/env Rscript

# Generate Publication Plots for V4 Methods
setwd("C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/simulation")

source("Comprehensive_Testing_Framework.R")

# Load results
if (file.exists("../results/v4_quick_validation_metrics.rds")) {
  metrics <- readRDS("../results/v4_quick_validation_metrics.rds")
  cat("Loaded metrics from previous run\n")
} else {
  # Run quick validation first
  cat("Running quick validation...\n")
  quick_results <- run_quick_validation(n_sim = 50, n_cores = 1)
  metrics <- quick_results$metrics
}

# Create output directory
dir.create("../plots", showWarnings = FALSE)

# Generate plots
png("../plots/v4_rmse_comparison.png", width = 1000, height = 600)
par(mar = c(5, 6, 4, 2))

# Sort by mean RMSE
method_summary <- metrics[, .(mean_rmse = mean(rmse, na.rm = TRUE),
                             mean_coverage = mean(coverage, na.rm = TRUE),
                             mean_bias = mean(bias, na.rm = TRUE)), by = method]
method_summary <- method_summary[order(mean_rmse)]

# RMSE comparison
barplot(method_summary$mean_rmse, names.arg = method_summary$method,
        las = 2, horiz = FALSE, col = "steelblue",
        main = "RMSE Comparison Across Methods",
        ylab = "Mean RMSE", xlab = "")
abline(h = mean(method_summary$mean_rmse), lty = 2, col = "red", lwd = 2)
legend("topright", legend = "Average", lty = 2, col = "red", lwd = 2, bty = "n")

dev.off()

# Coverage comparison
png("../plots/v4_coverage_comparison.png", width = 1000, height = 600)
par(mar = c(5, 6, 4, 2))

barplot(method_summary$mean_coverage * 100, names.arg = method_summary$method,
        las = 2, horiz = FALSE, col = ifelse(method_summary$mean_coverage >= 0.95,
                                               "darkgreen", "orange"),
        main = "Coverage Comparison Across Methods (Target: 95%)",
        ylab = "Mean Coverage (%)", xlab = "", ylim = c(0, 100))
abline(h = 95, lty = 2, col = "red", lwd = 2)
legend("topright", legend = "95% Target", lty = 2, col = "red", lwd = 2, bty = "n")

dev.off()

# Bias-RMSE scatter
png("../plots/v4_bias_rmse_scatter.png", width = 800, height = 800)
par(mar = c(5, 5, 4, 2))

plot(method_summary$mean_rmse, method_summary$mean_bias,
     pch = 19, cex = 2, col = "steelblue",
     xlab = "Mean RMSE", ylab = "Mean Bias",
     main = "Bias vs RMSE Trade-off")
text(method_summary$mean_rmse, method_summary$mean_bias,
     labels = method_summary$method, pos = 4, cex = 0.8)
abline(v = mean(method_summary$mean_rmse), lty = 2, col = "gray")
abline(h = mean(method_summary$mean_bias), lty = 2, col = "gray")

dev.off()

# Coverage-RMSE scatter
png("../plots/v4_coverage_rmse_scatter.png", width = 800, height = 800)
par(mar = c(5, 5, 4, 2))

plot(method_summary$mean_rmse, method_summary$mean_coverage * 100,
     pch = 19, cex = 2, col = "steelblue",
     xlab = "Mean RMSE", ylab = "Mean Coverage (%)",
     main = "Coverage vs RMSE Trade-off")
text(method_summary$mean_rmse, method_summary$mean_coverage * 100,
     labels = method_summary$method, pos = 4, cex = 0.8)
abline(h = 95, lty = 2, col = "red", lwd = 2)
abline(v = mean(method_summary$mean_rmse), lty = 2, col = "gray")
legend("bottomright", legend = "95% Target", lty = 2, col = "red", lwd = 2, bty = "n")

dev.off()

cat("\n=== Plots saved to ../plots/ ===\n")
cat("- v4_rmse_comparison.png\n")
cat("- v4_coverage_comparison.png\n")
cat("- v4_bias_rmse_scatter.png\n")
cat("- v4_coverage_rmse_scatter.png\n")
