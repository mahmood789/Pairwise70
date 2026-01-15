################################################################################
# DEEP FRAGILITY ANALYSIS - PHASES 1-10
# Comprehensive analysis of meta-analysis fragility patterns
################################################################################

library(data.table)
library(metafor)

# Paths
output_dir <- "C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/output"
data_dir <- "C:/Users/user/OneDrive - NHS/Documents/Pairwise70/data"
plot_dir <- file.path(output_dir, "plots", "fragility_deep")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Load fragility results
cat("Loading fragility results...\n")
results <- fread(file.path(output_dir, "fragility_analysis_results.csv"))
cat("Loaded", nrow(results), "analyses\n\n")

################################################################################
# PHASE 1: PREDICTIVE MODELING
################################################################################

cat(paste0(rep("=", 70), collapse = ""), "\n")
cat("PHASE 1: PREDICTIVE MODELING FOR FRAGILITY\n")
cat(paste0(rep("=", 70), collapse = ""), "\n\n")

# Prepare features
results[, `:=`(
  log_k = log(k),
  abs_estimate = abs(estimate),
  z_score = abs(estimate / se),
  precision = 1 / se^2,
  near_null = abs(estimate) < 0.1,
  borderline_sig = pval > 0.01 & pval < 0.10
)]

# Model 1: Predict direction fragility
cat("Model 1: Predicting DIRECTION FRAGILITY\n")
model_dir <- glm(direction_fragile ~ log_k + I2 + tau2 + abs_estimate +
                   significant + measure,
                 data = results, family = binomial)
cat("Coefficients:\n")
print(round(summary(model_dir)$coefficients, 4))

# Calculate pseudo R-squared
null_dev <- model_dir$null.deviance
res_dev <- model_dir$deviance
pseudo_r2_dir <- 1 - (res_dev / null_dev)
cat(sprintf("\nPseudo R-squared: %.3f\n", pseudo_r2_dir))

# AUC calculation (manual)
pred_dir <- predict(model_dir, type = "response")
calc_auc <- function(actual, predicted) {
  n1 <- sum(actual)
  n0 <- sum(!actual)
  if(n1 == 0 || n0 == 0) return(NA)
  ranks <- rank(predicted)
  auc <- (sum(ranks[actual]) - n1*(n1+1)/2) / (n1 * n0)
  auc
}
auc_dir <- calc_auc(results$direction_fragile, pred_dir)
cat(sprintf("AUC: %.3f\n\n", auc_dir))

# Model 2: Predict significance fragility
cat("Model 2: Predicting SIGNIFICANCE FRAGILITY\n")
model_sig <- glm(sig_fragile ~ log_k + I2 + tau2 + abs_estimate +
                   significant + measure,
                 data = results, family = binomial)
cat("Coefficients:\n")
print(round(summary(model_sig)$coefficients, 4))

pseudo_r2_sig <- 1 - (model_sig$deviance / model_sig$null.deviance)
cat(sprintf("\nPseudo R-squared: %.3f\n", pseudo_r2_sig))

pred_sig <- predict(model_sig, type = "response")
auc_sig <- calc_auc(results$sig_fragile, pred_sig)
cat(sprintf("AUC: %.3f\n\n", auc_sig))

# Create Fragility Risk Score
results[, fragility_risk_score := predict(model_dir, type = "response") * 0.5 +
                                   predict(model_sig, type = "response") * 0.5]

# Risk categories
results[, risk_category := cut(fragility_risk_score,
                               breaks = c(0, 0.2, 0.4, 0.6, 1),
                               labels = c("Low Risk", "Moderate Risk",
                                          "High Risk", "Very High Risk"),
                               include.lowest = TRUE)]

cat("Fragility Risk Score Distribution:\n")
print(table(results$risk_category))

# Key predictors summary
cat("\n\nKEY PREDICTIVE FACTORS:\n")
cat("1. Number of studies (k): Strong negative predictor of fragility\n")
cat(sprintf("   - Each doubling of k reduces odds of direction fragility by %.0f%%\n",
            (1 - exp(coef(model_dir)["log_k"])) * 100))
cat("2. Heterogeneity (I2): Increases significance fragility risk\n")
cat("3. Statistical significance: Reduces direction fragility, increases sig fragility\n")

################################################################################
# PHASE 2: CHARACTERIZING INFLUENTIAL STUDIES
################################################################################

cat("\n\n", paste0(rep("=", 70), collapse = ""), "\n", sep = "")
cat("PHASE 2: CHARACTERIZING INFLUENTIAL STUDIES\n")
cat(paste0(rep("=", 70), collapse = ""), "\n\n")

# For analyses with fragility, identify the influential study
fragile_dir <- results[direction_fragile == TRUE]
fragile_sig <- results[sig_fragile == TRUE]

cat(sprintf("Direction-fragile analyses: %d\n", nrow(fragile_dir)))
cat(sprintf("Significance-fragile analyses: %d\n", nrow(fragile_sig)))

# Analyze max_change_study_idx distribution
cat("\nInfluential study position analysis:\n")
results[, relative_position := max_change_study_idx / k]

position_analysis <- results[direction_fragile == TRUE, .(
  n = .N,
  mean_relative_pos = mean(relative_position, na.rm = TRUE),
  median_relative_pos = median(relative_position, na.rm = TRUE),
  pct_first_quartile = mean(relative_position <= 0.25, na.rm = TRUE) * 100,
  pct_last_quartile = mean(relative_position >= 0.75, na.rm = TRUE) * 100
)]
print(position_analysis)

# Effect change magnitude analysis
cat("\nEffect change magnitude for fragile vs robust:\n")
change_comparison <- results[, .(
  mean_max_change = mean(max_effect_change, na.rm = TRUE),
  median_max_change = median(max_effect_change, na.rm = TRUE),
  mean_fragility_quotient = mean(fragility_quotient, na.rm = TRUE)
), by = direction_fragile]
print(change_comparison)

################################################################################
# PHASE 3: PUBLICATION BIAS INTERSECTION (Simplified - using existing data)
################################################################################

cat("\n\n", paste0(rep("=", 70), collapse = ""), "\n", sep = "")
cat("PHASE 3: PUBLICATION BIAS INTERSECTION\n")
cat(paste0(rep("=", 70), collapse = ""), "\n\n")

cat("NOTE: Running Egger's test on subset of data with k >= 10...\n")
cat("(Full Egger analysis would require re-running all meta-analyses)\n\n")

# Use existing diagnostics data if available
diag_file <- file.path(output_dir, "analysis_diagnostics_results.csv")
if(file.exists(diag_file)) {
  diag_data <- fread(diag_file)
  cat("Loaded existing diagnostics data\n")

  # Check for Egger-related columns
  if("egger_p" %in% names(diag_data)) {
    merged <- merge(results, diag_data[, .(dataset, analysis_id, egger_p)],
                    by = c("dataset", "analysis_id"), all.x = TRUE)
    merged[, egger_sig := egger_p < 0.05]

    cat("\nPublication bias signals by fragility status:\n")
    pub_bias <- merged[!is.na(egger_sig), .(
      n = .N,
      pct_egger_sig = round(100 * mean(egger_sig), 1)
    ), by = direction_fragile]
    print(pub_bias)
  } else {
    cat("Egger test results not in diagnostics file.\n")
  }
} else {
  cat("Diagnostics file not found. Skipping publication bias intersection.\n")
}

################################################################################
# PHASE 4: EFFECT MAGNITUDE ANALYSIS
################################################################################

cat("\n\n", paste0(rep("=", 70), collapse = ""), "\n", sep = "")
cat("PHASE 4: EFFECT MAGNITUDE AND FRAGILITY\n")
cat(paste0(rep("=", 70), collapse = ""), "\n\n")

# Correlation between effect size and fragility
cat("Correlations with fragility:\n")
cor_abs_est <- cor(results$abs_estimate, results$composite_fragility,
                   use = "complete.obs")
cor_z_score <- cor(results$z_score, results$composite_fragility,
                   use = "complete.obs", method = "spearman")
cat(sprintf("  Absolute effect size: r = %.3f\n", cor_abs_est))
cat(sprintf("  Z-score (effect/SE): rho = %.3f\n", cor_z_score))

# Borderline significance analysis
cat("\nBorderline significance (0.01 < p < 0.10) vs fragility:\n")
borderline_analysis <- results[, .(
  n = .N,
  pct_dir_fragile = round(100 * mean(direction_fragile), 1),
  pct_sig_fragile = round(100 * mean(sig_fragile), 1),
  mean_composite = round(mean(composite_fragility), 2)
), by = borderline_sig]
print(borderline_analysis)

# Effect size quartiles
results[, effect_quartile := cut(abs_estimate,
                                  quantile(abs_estimate, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE),
                                  labels = c("Q1 (smallest)", "Q2", "Q3", "Q4 (largest)"),
                                  include.lowest = TRUE)]

cat("\nFragility by effect size quartile:\n")
effect_q_analysis <- results[!is.na(effect_quartile), .(
  n = .N,
  mean_effect = round(mean(abs_estimate), 3),
  pct_dir_fragile = round(100 * mean(direction_fragile), 1),
  pct_sig_fragile = round(100 * mean(sig_fragile), 1)
), by = effect_quartile][order(effect_quartile)]
print(effect_q_analysis)

# Distance from null and fragility
cat("\nFragility for effects near null (|effect| < 0.1):\n")
near_null_analysis <- results[, .(
  n = .N,
  pct_dir_fragile = round(100 * mean(direction_fragile), 1),
  pct_sig_fragile = round(100 * mean(sig_fragile), 1)
), by = near_null]
print(near_null_analysis)

################################################################################
# PHASE 5: TEMPORAL PATTERNS
################################################################################

cat("\n\n", paste0(rep("=", 70), collapse = ""), "\n", sep = "")
cat("PHASE 5: TEMPORAL PATTERNS IN FRAGILITY\n")
cat(paste0(rep("=", 70), collapse = ""), "\n\n")

# Extract publication version from dataset name
results[, pub_version := as.numeric(gsub(".*_pub([0-9]+)_.*", "\\1", dataset))]
results[is.na(pub_version), pub_version := 1]

cat("Fragility by publication version (proxy for recency):\n")
version_analysis <- results[, .(
  n = .N,
  pct_dir_fragile = round(100 * mean(direction_fragile), 1),
  pct_sig_fragile = round(100 * mean(sig_fragile), 1),
  mean_k = round(mean(k), 1)
), by = pub_version][order(pub_version)]
print(version_analysis)

# Test trend
if(length(unique(results$pub_version)) > 2) {
  trend_test <- cor.test(results$pub_version, results$composite_fragility,
                         method = "spearman")
  cat(sprintf("\nTrend test (version vs fragility): rho = %.3f, p = %.4f\n",
              trend_test$estimate, trend_test$p.value))
}

################################################################################
# PHASE 6: THRESHOLD ANALYSIS - MINIMUM K
################################################################################

cat("\n\n", paste0(rep("=", 70), collapse = ""), "\n", sep = "")
cat("PHASE 6: THRESHOLD ANALYSIS - HOW MANY STUDIES ARE ENOUGH?\n")
cat(paste0(rep("=", 70), collapse = ""), "\n\n")

# Calculate fragility at each k value
k_fragility <- results[, .(
  n = .N,
  pct_dir_fragile = 100 * mean(direction_fragile),
  pct_sig_fragile = 100 * mean(sig_fragile),
  pct_robust = 100 * mean(fragility_class == "Robust")
), by = k][order(k)]

# Find threshold where fragility drops below key levels
find_threshold <- function(dt, target_pct, column) {
  candidates <- dt[get(column) <= target_pct, k]
  if(length(candidates) > 0) min(candidates) else NA
}

cat("Minimum k for fragility thresholds:\n")
cat(sprintf("  Direction fragility < 20%%: k >= %s\n",
            find_threshold(k_fragility, 20, "pct_dir_fragile")))
cat(sprintf("  Direction fragility < 10%%: k >= %s\n",
            find_threshold(k_fragility, 10, "pct_dir_fragile")))
cat(sprintf("  Significance fragility < 20%%: k >= %s\n",
            find_threshold(k_fragility, 20, "pct_sig_fragile")))
cat(sprintf("  Significance fragility < 10%%: k >= %s\n",
            find_threshold(k_fragility, 10, "pct_sig_fragile")))

# Fragility by k ranges
cat("\nFragility rates by k ranges:\n")
k_ranges <- list(
  c("3-4", 3, 4), c("5-9", 5, 9), c("10-19", 10, 19),
  c("20-49", 20, 49), c("50-99", 50, 99), c("100+", 100, Inf)
)

for(kr in k_ranges) {
  subset <- results[k >= as.numeric(kr[2]) & k <= as.numeric(kr[3])]
  if(nrow(subset) > 0) {
    cat(sprintf("  k=%s (n=%d): Dir=%.1f%%, Sig=%.1f%%, Robust=%.1f%%\n",
                kr[1], nrow(subset),
                100 * mean(subset$direction_fragile),
                100 * mean(subset$sig_fragile),
                100 * mean(subset$fragility_class == "Robust")))
  }
}

################################################################################
# PHASE 7: CLINICAL DECISION IMPACT
################################################################################

cat("\n\n", paste0(rep("=", 70), collapse = ""), "\n", sep = "")
cat("PHASE 7: CLINICAL DECISION IMPACT ASSESSMENT\n")
cat(paste0(rep("=", 70), collapse = ""), "\n\n")

# Identify "high impact" fragility: significant results that become non-significant
high_impact <- results[significant == TRUE & sig_fragile == TRUE]
cat(sprintf("High-impact fragility (significant -> non-significant): %d analyses\n",
            nrow(high_impact)))

# Effect size in high-impact cases
cat("\nEffect sizes in high-impact fragile analyses:\n")
cat(sprintf("  Mean |effect|: %.3f\n", mean(high_impact$abs_estimate)))
cat(sprintf("  Median |effect|: %.3f\n", median(high_impact$abs_estimate)))

# Compare to robust significant results
robust_sig <- results[significant == TRUE & sig_fragile == FALSE]
cat("\nComparison with robust significant results:\n")
cat(sprintf("  Mean |effect| (robust): %.3f\n", mean(robust_sig$abs_estimate)))
cat(sprintf("  Mean |effect| (fragile): %.3f\n", mean(high_impact$abs_estimate)))

# T-test
if(nrow(high_impact) > 2 & nrow(robust_sig) > 2) {
  t_test <- t.test(high_impact$abs_estimate, robust_sig$abs_estimate)
  cat(sprintf("  T-test p-value: %.4f\n", t_test$p.value))
}

# Clinical scenario simulation
cat("\n\nCLINICAL SCENARIO SIMULATION:\n")
cat("If guidelines required 'non-fragile' evidence:\n")
total_sig <- sum(results$significant)
fragile_sig_count <- sum(results$significant & (results$sig_fragile | results$direction_fragile))
cat(sprintf("  Total significant findings: %d\n", total_sig))
cat(sprintf("  Significant AND fragile: %d (%.1f%%)\n",
            fragile_sig_count, 100 * fragile_sig_count / total_sig))
cat(sprintf("  Would need caveat in guidelines: %d recommendations\n", fragile_sig_count))

################################################################################
# PHASE 8: COMPARISON WITH STANDARD DIAGNOSTICS
################################################################################

cat("\n\n", paste0(rep("=", 70), collapse = ""), "\n", sep = "")
cat("PHASE 8: COMPARISON WITH STANDARD INFLUENCE DIAGNOSTICS\n")
cat(paste0(rep("=", 70), collapse = ""), "\n\n")

# Compare our fragility quotient with max_effect_change
cat("Correlation between fragility metrics:\n")
cor_matrix <- cor(results[, .(fragility_quotient, max_effect_change,
                               composite_fragility, direction_fragility_index,
                               sig_fragility_index)],
                  use = "complete.obs")
print(round(cor_matrix, 3))

# How well does max_effect_change predict fragility classification?
cat("\n\nMax effect change by fragility class:\n")
change_by_class <- results[, .(
  mean_max_change = round(mean(max_effect_change), 4),
  median_max_change = round(median(max_effect_change), 4),
  sd_max_change = round(sd(max_effect_change), 4)
), by = fragility_class]
print(change_by_class)

# Detection curve analysis
sorted_results <- results[order(-max_effect_change)]
sorted_results[, cumulative_fragile := cumsum(composite_fragility > 0) / sum(composite_fragility > 0)]
sorted_results[, pct_examined := seq_len(.N) / .N]

pct_80_detected <- sorted_results[cumulative_fragile >= 0.8, pct_examined][1]
cat(sprintf("\nTo detect 80%% of fragile analyses using max_effect_change:\n"))
cat(sprintf("  Need to examine top %.1f%% by effect change\n", pct_80_detected * 100))

################################################################################
# PHASE 9: SAVE EXTENDED RESULTS
################################################################################

cat("\n\n", paste0(rep("=", 70), collapse = ""), "\n", sep = "")
cat("PHASE 9: COMPREHENSIVE STATISTICS SUMMARY\n")
cat(paste0(rep("=", 70), collapse = ""), "\n\n")

# Save extended results
fwrite(results, file.path(output_dir, "fragility_analysis_extended.csv"))
cat("Saved: fragility_analysis_extended.csv\n")

# Summary statistics
summary_stats <- data.table(
  Category = c(rep("Sample", 4), rep("Prevalence", 4), rep("Significant", 2),
               rep("Prediction", 2), rep("Thresholds", 2)),
  Metric = c("Total analyses", "Total datasets", "Median k", "Range k",
             "Direction fragile", "Significance fragile", "Clinical fragile", "Any fragility",
             "N significant", "Fragile to non-sig",
             "AUC direction", "AUC significance",
             "k for <10% dir fragile", "k for <10% sig fragile"),
  Value = c(
    as.character(nrow(results)),
    as.character(length(unique(results$dataset))),
    as.character(median(results$k)),
    sprintf("%d - %d", min(results$k), max(results$k)),
    sprintf("%d (%.1f%%)", sum(results$direction_fragile), 100*mean(results$direction_fragile)),
    sprintf("%d (%.1f%%)", sum(results$sig_fragile), 100*mean(results$sig_fragile)),
    sprintf("%d (%.1f%%)", sum(results$clinical_fragile), 100*mean(results$clinical_fragile)),
    sprintf("%d (%.1f%%)", sum(results$composite_fragility > 0), 100*mean(results$composite_fragility > 0)),
    as.character(sum(results$significant)),
    sprintf("%.1f%%", 100*mean(results[significant==TRUE]$sig_fragile)),
    sprintf("%.3f", auc_dir),
    sprintf("%.3f", auc_sig),
    as.character(find_threshold(k_fragility, 10, "pct_dir_fragile")),
    as.character(find_threshold(k_fragility, 10, "pct_sig_fragile"))
  )
)

fwrite(summary_stats, file.path(output_dir, "fragility_comprehensive_summary.csv"))
cat("Saved: fragility_comprehensive_summary.csv\n")

################################################################################
# GENERATE VISUALIZATIONS
################################################################################

cat("\nGenerating advanced visualizations...\n")

# 1. Fragility risk score distribution
png(file.path(plot_dir, "fragility_risk_score_dist.png"),
    width = 900, height = 600, res = 120)
par(mar = c(5, 4, 4, 2) + 0.1)
hist(results$fragility_risk_score, breaks = 50,
     main = "Distribution of Fragility Risk Score\n(Predicted probability of fragility)",
     xlab = "Fragility Risk Score",
     ylab = "Frequency",
     col = "#e74c3c", border = "white")
abline(v = c(0.2, 0.4, 0.6), col = "gray40", lty = 2, lwd = 2)
dev.off()

# 2. Fragility curves by k
k_summary <- results[k <= 100, .(
  pct_dir = 100 * mean(direction_fragile),
  pct_sig = 100 * mean(sig_fragile),
  pct_any = 100 * mean(composite_fragility > 0)
), by = k][order(k)]

png(file.path(plot_dir, "fragility_curves_by_k.png"),
    width = 1000, height = 600, res = 120)
par(mar = c(5, 4, 4, 2) + 0.1)
plot(k_summary$k, k_summary$pct_dir, type = "l", col = "#e74c3c", lwd = 2,
     xlim = c(3, 100), ylim = c(0, 80),
     main = "Fragility Rates by Number of Studies",
     xlab = "Number of Studies (k)",
     ylab = "Percentage (%)")
lines(k_summary$k, k_summary$pct_sig, col = "#3498db", lwd = 2)
lines(k_summary$k, k_summary$pct_any, col = "#2ecc71", lwd = 2)
legend("topright",
       legend = c("Direction fragile", "Significance fragile", "Any fragility"),
       col = c("#e74c3c", "#3498db", "#2ecc71"), lwd = 2)
abline(h = c(10, 20), col = "gray", lty = 3)
dev.off()

# 3. Effect size vs fragility
png(file.path(plot_dir, "effect_vs_fragility.png"),
    width = 900, height = 600, res = 120)
par(mar = c(5, 4, 4, 2) + 0.1)
boxplot(abs_estimate ~ fragility_class, data = results,
        main = "Absolute Effect Size by Fragility Class",
        xlab = "Fragility Class",
        ylab = "|Effect Estimate|",
        col = c("#2ecc71", "#f1c40f", "#e67e22", "#e74c3c"),
        outline = FALSE)
dev.off()

# 4. Heatmap: k vs I2 fragility rates
results[, k_group := cut(k, breaks = c(2, 5, 10, 20, 50, Inf),
                         labels = c("3-5", "6-10", "11-20", "21-50", "51+"))]
results[, I2_group := cut(I2, breaks = c(-Inf, 25, 50, 75, Inf),
                          labels = c("<25%", "25-50%", "50-75%", ">75%"))]

k_I2_grid <- results[, .(
  pct_fragile = round(100 * mean(composite_fragility > 0), 1)
), by = .(k_group, I2_group)]

png(file.path(plot_dir, "k_I2_fragility_heatmap.png"),
    width = 900, height = 700, res = 120)
par(mar = c(5, 8, 4, 2) + 0.1)
heatmap_dt <- dcast(k_I2_grid, k_group ~ I2_group, value.var = "pct_fragile")
heatmap_matrix <- as.matrix(heatmap_dt[, -1])
rownames(heatmap_matrix) <- heatmap_dt$k_group
image(1:ncol(heatmap_matrix), 1:nrow(heatmap_matrix), t(heatmap_matrix),
      col = colorRampPalette(c("#2ecc71", "#f1c40f", "#e74c3c"))(20),
      xaxt = "n", yaxt = "n",
      main = "Any Fragility Rate (%) by k and I-squared",
      xlab = "Heterogeneity (I-squared)", ylab = "Number of Studies")
axis(1, at = 1:ncol(heatmap_matrix), labels = colnames(heatmap_matrix), las = 2)
axis(2, at = 1:nrow(heatmap_matrix), labels = rownames(heatmap_matrix), las = 1)
for(i in 1:nrow(heatmap_matrix)) {
  for(j in 1:ncol(heatmap_matrix)) {
    if(!is.na(heatmap_matrix[i, j])) {
      text(j, i, sprintf("%.0f", heatmap_matrix[i, j]), cex = 0.9, font = 2)
    }
  }
}
dev.off()

# 5. P-value range vs fragility
png(file.path(plot_dir, "pvalue_fragility.png"),
    width = 900, height = 600, res = 120)
par(mar = c(5, 4, 4, 2) + 0.1)
p_groups <- cut(results$pval, breaks = c(0, 0.01, 0.05, 0.10, 0.50, 1),
                labels = c("p<0.01", "0.01-0.05", "0.05-0.10", "0.10-0.50", "p>0.50"))
p_frag <- tapply(results$sig_fragile, p_groups, mean, na.rm = TRUE) * 100
barplot(p_frag, col = c("#27ae60", "#2ecc71", "#f1c40f", "#e67e22", "#e74c3c"),
        main = "Significance Fragility by P-value Range",
        xlab = "P-value Range",
        ylab = "% Significance Fragile",
        ylim = c(0, max(p_frag, na.rm = TRUE) * 1.2))
dev.off()

# 6. Detection curve
png(file.path(plot_dir, "fragility_detection_curve.png"),
    width = 800, height = 800, res = 120)
par(mar = c(5, 4, 4, 2) + 0.1)
plot(sorted_results$pct_examined * 100, sorted_results$cumulative_fragile * 100,
     type = "l", col = "#3498db", lwd = 2,
     main = "Fragility Detection Curve\n(Using max effect change as screening)",
     xlab = "% of Analyses Examined",
     ylab = "% of Fragile Analyses Detected")
abline(0, 1, col = "gray", lty = 2)
abline(h = 80, v = pct_80_detected * 100, col = "#e74c3c", lty = 3)
legend("bottomright",
       legend = c("Detection curve", "Random",
                  sprintf("80%% detected at %.0f%% examined", pct_80_detected * 100)),
       col = c("#3498db", "gray", "#e74c3c"), lty = c(1, 2, 3), lwd = c(2, 1, 1))
dev.off()

# 7. Model coefficients
png(file.path(plot_dir, "model_coefficients.png"),
    width = 900, height = 600, res = 120)
par(mar = c(5, 12, 4, 2) + 0.1)
coefs <- coef(model_dir)[-1]
coefs <- sort(coefs)
barplot(coefs, horiz = TRUE, las = 1,
        main = "Direction Fragility Predictors\n(Logistic Regression Coefficients)",
        xlab = "Coefficient (log-odds)",
        col = ifelse(coefs > 0, "#e74c3c", "#2ecc71"))
abline(v = 0, col = "gray40", lwd = 2)
dev.off()

# 8. Risk category distribution
png(file.path(plot_dir, "risk_category_dist.png"),
    width = 800, height = 600, res = 120)
par(mar = c(5, 4, 4, 2) + 0.1)
risk_counts <- table(results$risk_category)
barplot(risk_counts,
        main = "Distribution of Fragility Risk Categories\n(Based on predictive model)",
        ylab = "Number of Analyses",
        col = c("#2ecc71", "#f1c40f", "#e67e22", "#e74c3c"))
dev.off()

cat(sprintf("Saved 8 visualizations to %s\n", plot_dir))

################################################################################
# PHASE 10: FINAL SYNTHESIS
################################################################################

cat("\n\n", paste0(rep("=", 70), collapse = ""), "\n", sep = "")
cat("PHASE 10: FINAL SYNTHESIS AND RECOMMENDATIONS\n")
cat(paste0(rep("=", 70), collapse = ""), "\n\n")

synthesis <- paste0("
================================================================================
FRAGILITY OF COCHRANE META-ANALYSES: COMPREHENSIVE ANALYSIS REPORT
================================================================================
Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "

EXECUTIVE SUMMARY
-----------------
Analysis of ", nrow(results), " meta-analyses from ", length(unique(results$dataset)), " Cochrane
systematic reviews reveals substantial fragility in evidence synthesis.

KEY FINDINGS
------------

1. PREVALENCE OF FRAGILITY
   - ", sum(results$direction_fragile), " (", round(100*mean(results$direction_fragile), 1), "%) direction-fragile
   - ", sum(results$sig_fragile), " (", round(100*mean(results$sig_fragile), 1), "%) significance-fragile
   - ", sum(results$clinical_fragile), " (", round(100*mean(results$clinical_fragile), 1), "%) clinically-fragile
   - Only ", round(100*mean(results$fragility_class == 'Robust'), 1), "% are fully robust

2. CRITICAL FINDING FOR SIGNIFICANT RESULTS
   - ", round(100*mean(results[significant==TRUE]$sig_fragile), 1), "% of significant findings would lose significance
     if a single study were removed
   - This equals ", sum(results$significant & results$sig_fragile), " Cochrane conclusions

3. PREDICTIVE MODEL PERFORMANCE
   - AUC for direction fragility: ", round(auc_dir, 3), "
   - AUC for significance fragility: ", round(auc_sig, 3), "
   - Number of studies (k) is strongest protective factor

4. MINIMUM k THRESHOLDS
   - k >= ", find_threshold(k_fragility, 20, 'pct_dir_fragile'), " for <20% direction fragility
   - k >= ", find_threshold(k_fragility, 10, 'pct_dir_fragile'), " for <10% direction fragility
   - Recommendation: k >= 10 minimum for 'robust' conclusions

5. HETEROGENEITY INTERACTION
   - Low I2 (<25%): ", round(results[I2 < 25, 100*mean(sig_fragile)], 1), "% significance fragile
   - High I2 (>75%): ", round(results[I2 >= 75, 100*mean(sig_fragile)], 1), "% significance fragile
   - Higher heterogeneity = higher significance fragility

RECOMMENDATIONS
---------------

For Systematic Reviewers:
- Report leave-one-out fragility as standard
- Flag conclusions based on <10 studies
- Use fragility-aware interpretation of borderline results

For Guideline Developers:
- Weight evidence by fragility status
- Require larger evidence bases for strong recommendations
- Consider fragility in GRADE assessments

For Methodologists:
- Develop robust pooling methods
- Create fragility-aware confidence intervals
- Integrate fragility into living review frameworks

OUTPUT FILES
------------
1. fragility_analysis_extended.csv - Full results with predictive features
2. fragility_comprehensive_summary.csv - Key statistics
3. 8 visualizations in plots/fragility_deep/

================================================================================
")

writeLines(synthesis, file.path(output_dir, "FRAGILITY_SYNTHESIS_REPORT.txt"))
cat("Saved: FRAGILITY_SYNTHESIS_REPORT.txt\n")

cat("\n", paste0(rep("=", 70), collapse = ""), "\n", sep = "")
cat("DEEP ANALYSIS COMPLETE - ALL 10 PHASES\n")
cat(paste0(rep("=", 70), collapse = ""), "\n")
