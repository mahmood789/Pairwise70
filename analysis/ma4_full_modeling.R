# ============================================================================
# MA4 v1.0.1 COMPREHENSIVE STATISTICAL MODELING
# Full analysis of 5,088 meta-analyses from Pairwise70 Cochrane Collection
# ============================================================================

cat("
================================================================================
   MA4 v1.0.1 COMPREHENSIVE STATISTICAL MODELING
   Pairwise70 Cochrane Collection Analysis
================================================================================
\n")

# Load results
results <- read.csv("C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/ma4_results_pairwise70.csv",
                    stringsAsFactors = FALSE)

cat("Data loaded:", nrow(results), "meta-analyses from",
    length(unique(results$review_id)), "Cochrane reviews\n\n")

# ============================================================================
# SECTION 1: DESCRIPTIVE STATISTICS
# ============================================================================

cat("
================================================================================
SECTION 1: DESCRIPTIVE STATISTICS
================================================================================
\n")

# Basic summary
cat("--- 1.1 Overall Summary ---\n")
cat("Total meta-analyses:", nrow(results), "\n")
cat("Unique Cochrane reviews:", length(unique(results$review_id)), "\n")
cat("Effect types:", paste(unique(results$effect_type), collapse = ", "), "\n")
cat("Tau estimators used:", paste(unique(results$tau_estimator), collapse = ", "), "\n")
cat("R status values:", paste(unique(results$R_status), collapse = ", "), "\n\n")

# Number of studies (k) distribution
cat("--- 1.2 Number of Studies (k) Distribution ---\n")
cat("Min k:", min(results$k), "\n")
cat("Q1 k:", quantile(results$k, 0.25), "\n")
cat("Median k:", median(results$k), "\n")
cat("Mean k:", round(mean(results$k), 2), "\n")
cat("Q3 k:", quantile(results$k, 0.75), "\n")
cat("Max k:", max(results$k), "\n")
cat("SD k:", round(sd(results$k), 2), "\n\n")

# k categories
results$k_cat <- cut(results$k,
                     breaks = c(1, 2, 3, 5, 10, 20, 50, Inf),
                     labels = c("2", "3", "4-5", "6-10", "11-20", "21-50", "51+"),
                     include.lowest = TRUE)
cat("k distribution by category:\n")
print(table(results$k_cat))
cat("\n")

# R (stability) distribution
cat("--- 1.3 R Stability Distribution ---\n")
R_valid <- results$R[!is.na(results$R)]
cat("n with valid R:", length(R_valid), "\n")
cat("Min R:", round(min(R_valid), 4), "\n")
cat("Q1 R:", round(quantile(R_valid, 0.25), 4), "\n")
cat("Median R:", round(median(R_valid), 4), "\n")
cat("Mean R:", round(mean(R_valid), 4), "\n")
cat("Q3 R:", round(quantile(R_valid, 0.75), 4), "\n")
cat("Max R:", round(max(R_valid), 4), "\n")
cat("SD R:", round(sd(R_valid), 4), "\n\n")

# R categories
results$R_cat <- cut(results$R,
                     breaks = c(0, 0.3, 0.5, 0.7, 0.8, 0.9, 1.0),
                     labels = c("Very Low (<0.3)", "Low (0.3-0.5)", "Moderate (0.5-0.7)",
                               "Good (0.7-0.8)", "High (0.8-0.9)", "Excellent (0.9-1.0)"),
                     include.lowest = TRUE)
cat("R distribution by category:\n")
print(table(results$R_cat))
cat("Percentages:\n")
print(round(prop.table(table(results$R_cat)) * 100, 1))
cat("\n")

# Theta distribution
cat("--- 1.4 Effect Size (theta) Distribution ---\n")
theta_valid <- results$theta[!is.na(results$theta) & is.finite(results$theta)]
cat("n with valid theta:", length(theta_valid), "\n")
cat("Min theta:", round(min(theta_valid), 4), "\n")
cat("Q1 theta:", round(quantile(theta_valid, 0.25), 4), "\n")
cat("Median theta:", round(median(theta_valid), 4), "\n")
cat("Mean theta:", round(mean(theta_valid), 4), "\n")
cat("Q3 theta:", round(quantile(theta_valid, 0.75), 4), "\n")
cat("Max theta:", round(max(theta_valid), 4), "\n\n")

# Tau (heterogeneity) distribution
cat("--- 1.5 Heterogeneity (tau) Distribution ---\n")
tau_valid <- results$tau[!is.na(results$tau) & is.finite(results$tau)]
cat("n with valid tau:", length(tau_valid), "\n")
cat("tau = 0 (no heterogeneity):", sum(tau_valid < 1e-6),
    "(", round(100*mean(tau_valid < 1e-6), 1), "%)\n")
cat("Min tau (>0):", round(min(tau_valid[tau_valid >= 1e-6]), 4), "\n")
cat("Median tau:", round(median(tau_valid), 4), "\n")
cat("Mean tau:", round(mean(tau_valid), 4), "\n")
cat("Q3 tau:", round(quantile(tau_valid, 0.75), 4), "\n")
cat("Max tau:", round(max(tau_valid), 4), "\n\n")

# Sigma distribution
cat("--- 1.6 Standard Error (sigma) Distribution ---\n")
sigma_valid <- results$sigma[!is.na(results$sigma) & is.finite(results$sigma) & results$sigma > 0]
cat("n with valid sigma:", length(sigma_valid), "\n")
cat("Min sigma:", round(min(sigma_valid), 4), "\n")
cat("Median sigma:", round(median(sigma_valid), 4), "\n")
cat("Mean sigma:", round(mean(sigma_valid), 4), "\n")
cat("Max sigma:", round(max(sigma_valid), 4), "\n\n")

# By effect type
cat("--- 1.7 Summary by Effect Type ---\n")
for (et in unique(results$effect_type)) {
  sub <- results[results$effect_type == et, ]
  cat("\n", et, ":\n", sep = "")
  cat("  n =", nrow(sub), "\n")
  cat("  Median k =", median(sub$k), "\n")
  cat("  Mean R =", round(mean(sub$R, na.rm = TRUE), 4), "\n")
  cat("  Median R =", round(median(sub$R, na.rm = TRUE), 4), "\n")
  cat("  Mean tau =", round(mean(sub$tau, na.rm = TRUE), 4), "\n")
  cat("  Mean |theta| =", round(mean(abs(sub$theta), na.rm = TRUE), 4), "\n")
}

# ============================================================================
# SECTION 2: CORRELATION ANALYSIS
# ============================================================================

cat("\n
================================================================================
SECTION 2: CORRELATION ANALYSIS
================================================================================
\n")

# Create numeric dataset for correlations
num_data <- data.frame(
  R = results$R,
  k = results$k,
  log_k = log(results$k),
  theta = results$theta,
  abs_theta = abs(results$theta),
  sigma = results$sigma,
  log_sigma = log(results$sigma + 1e-10),
  tau = results$tau,
  log_tau = log(results$tau + 1e-6),
  tau_over_sigma = results$tau / (results$sigma + 1e-10),
  precision = 1 / (results$sigma + 1e-10),
  I2_proxy = results$tau^2 / (results$tau^2 + results$sigma^2 + 1e-10)
)

# Remove non-finite values
complete_idx <- complete.cases(num_data) &
                is.finite(num_data$theta) &
                is.finite(num_data$log_sigma) &
                is.finite(num_data$log_tau)
num_complete <- num_data[complete_idx, ]

cat("Complete cases for correlation:", nrow(num_complete), "\n\n")

cat("--- 2.1 Pearson Correlations with R ---\n")
cor_vars <- c("k", "log_k", "abs_theta", "sigma", "log_sigma",
              "tau", "log_tau", "tau_over_sigma", "I2_proxy")
for (v in cor_vars) {
  r <- cor(num_complete$R, num_complete[[v]], use = "complete.obs")
  test <- cor.test(num_complete$R, num_complete[[v]])
  cat(sprintf("  cor(R, %s) = %.4f (p %s)\n",
              v, r,
              ifelse(test$p.value < 0.001, "< 0.001", sprintf("= %.4f", test$p.value))))
}

cat("\n--- 2.2 Spearman Correlations with R ---\n")
for (v in cor_vars) {
  rho <- cor(num_complete$R, num_complete[[v]], method = "spearman", use = "complete.obs")
  test <- cor.test(num_complete$R, num_complete[[v]], method = "spearman")
  cat(sprintf("  rho(R, %s) = %.4f (p %s)\n",
              v, rho,
              ifelse(test$p.value < 0.001, "< 0.001", sprintf("= %.4f", test$p.value))))
}

cat("\n--- 2.3 Full Correlation Matrix (key variables) ---\n")
key_vars <- c("R", "k", "abs_theta", "sigma", "tau", "I2_proxy")
cor_matrix <- cor(num_complete[, key_vars], use = "complete.obs")
print(round(cor_matrix, 3))

# ============================================================================
# SECTION 3: REGRESSION MODELING
# ============================================================================

cat("\n
================================================================================
SECTION 3: REGRESSION MODELING
================================================================================
\n")

# Prepare modeling dataset
model_data <- results[complete.cases(results[, c("R", "k", "theta", "sigma", "tau", "effect_type")]), ]
model_data <- model_data[is.finite(model_data$theta) & is.finite(model_data$sigma) & model_data$sigma > 0, ]
model_data$log_k <- log(model_data$k)
model_data$abs_theta <- abs(model_data$theta)
model_data$log_sigma <- log(model_data$sigma)
model_data$log_tau <- log(model_data$tau + 1e-6)
model_data$has_heterogeneity <- ifelse(model_data$tau > 1e-6, 1, 0)
model_data$theta_sign <- sign(model_data$theta)
model_data$theta_near_zero <- ifelse(abs(model_data$theta) < 0.1, 1, 0)
model_data$I2_proxy <- model_data$tau^2 / (model_data$tau^2 + model_data$sigma^2)

cat("Modeling dataset:", nrow(model_data), "observations\n\n")

# Model 1: Simple linear regression
cat("--- 3.1 Model 1: R ~ log(k) ---\n")
m1 <- lm(R ~ log_k, data = model_data)
cat("Coefficients:\n")
print(round(summary(m1)$coefficients, 4))
cat("R-squared:", round(summary(m1)$r.squared, 4), "\n")
cat("Adj R-squared:", round(summary(m1)$adj.r.squared, 4), "\n\n")

# Model 2: Add effect type
cat("--- 3.2 Model 2: R ~ log(k) + effect_type ---\n")
m2 <- lm(R ~ log_k + effect_type, data = model_data)
cat("Coefficients:\n")
print(round(summary(m2)$coefficients, 4))
cat("R-squared:", round(summary(m2)$r.squared, 4), "\n")
cat("Adj R-squared:", round(summary(m2)$adj.r.squared, 4), "\n\n")

# Model 3: Add heterogeneity
cat("--- 3.3 Model 3: R ~ log(k) + effect_type + log(tau+1e-6) ---\n")
m3 <- lm(R ~ log_k + effect_type + log_tau, data = model_data)
cat("Coefficients:\n")
print(round(summary(m3)$coefficients, 4))
cat("R-squared:", round(summary(m3)$r.squared, 4), "\n")
cat("Adj R-squared:", round(summary(m3)$adj.r.squared, 4), "\n\n")

# Model 4: Add absolute effect size
cat("--- 3.4 Model 4: R ~ log(k) + effect_type + log(tau) + |theta| ---\n")
m4 <- lm(R ~ log_k + effect_type + log_tau + abs_theta, data = model_data)
cat("Coefficients:\n")
print(round(summary(m4)$coefficients, 4))
cat("R-squared:", round(summary(m4)$r.squared, 4), "\n")
cat("Adj R-squared:", round(summary(m4)$adj.r.squared, 4), "\n\n")

# Model 5: Full model with interactions
cat("--- 3.5 Model 5: Full model with theta near zero indicator ---\n")
m5 <- lm(R ~ log_k + effect_type + log_tau + abs_theta + theta_near_zero +
           log_k:effect_type, data = model_data)
cat("Coefficients:\n")
print(round(summary(m5)$coefficients, 4))
cat("R-squared:", round(summary(m5)$r.squared, 4), "\n")
cat("Adj R-squared:", round(summary(m5)$adj.r.squared, 4), "\n\n")

# Model 6: I2-based model
cat("--- 3.6 Model 6: R ~ log(k) + effect_type + I2_proxy ---\n")
m6 <- lm(R ~ log_k + effect_type + I2_proxy, data = model_data)
cat("Coefficients:\n")
print(round(summary(m6)$coefficients, 4))
cat("R-squared:", round(summary(m6)$r.squared, 4), "\n")
cat("Adj R-squared:", round(summary(m6)$adj.r.squared, 4), "\n\n")

# Model comparison
cat("--- 3.7 Model Comparison (AIC/BIC) ---\n")
models <- list(m1 = m1, m2 = m2, m3 = m3, m4 = m4, m5 = m5, m6 = m6)
model_comp <- data.frame(
  Model = names(models),
  AIC = sapply(models, AIC),
  BIC = sapply(models, BIC),
  R2 = sapply(models, function(m) summary(m)$r.squared),
  Adj_R2 = sapply(models, function(m) summary(m)$adj.r.squared),
  df = sapply(models, function(m) m$df.residual)
)
print(model_comp)
cat("\nBest model by AIC:", model_comp$Model[which.min(model_comp$AIC)], "\n")
cat("Best model by BIC:", model_comp$Model[which.min(model_comp$BIC)], "\n\n")

# ANOVA for nested models
cat("--- 3.8 ANOVA: Model Comparison ---\n")
cat("m1 vs m2 (adding effect_type):\n")
print(anova(m1, m2))
cat("\nm2 vs m3 (adding log_tau):\n")
print(anova(m2, m3))
cat("\nm3 vs m4 (adding abs_theta):\n")
print(anova(m3, m4))

# ============================================================================
# SECTION 4: BEST MODEL DIAGNOSTICS
# ============================================================================

cat("\n
================================================================================
SECTION 4: BEST MODEL DIAGNOSTICS (Model 5)
================================================================================
\n")

best_model <- m5

cat("--- 4.1 Residual Statistics ---\n")
resid <- residuals(best_model)
cat("Mean residual:", round(mean(resid), 6), "(should be ~0)\n")
cat("SD residual:", round(sd(resid), 4), "\n")
cat("Min residual:", round(min(resid), 4), "\n")
cat("Max residual:", round(max(resid), 4), "\n")
cat("Skewness:", round(mean((resid - mean(resid))^3) / sd(resid)^3, 4), "\n")
cat("Kurtosis:", round(mean((resid - mean(resid))^4) / sd(resid)^4, 4), "(normal = 3)\n\n")

cat("--- 4.2 Shapiro-Wilk Test (normality, sample of 5000) ---\n")
set.seed(42)
sample_resid <- sample(resid, min(5000, length(resid)))
sw_test <- shapiro.test(sample_resid)
cat("W =", round(sw_test$statistic, 4), ", p-value =", format(sw_test$p.value, digits = 4), "\n")
cat("Note: With large n, even small deviations are significant\n\n")

cat("--- 4.3 Variance Inflation Factors ---\n")
# Manual VIF calculation
vif_calc <- function(model) {
  X <- model.matrix(model)[, -1]  # Remove intercept
  vifs <- numeric(ncol(X))
  names(vifs) <- colnames(X)
  for (j in 1:ncol(X)) {
    r2 <- summary(lm(X[, j] ~ X[, -j]))$r.squared
    vifs[j] <- 1 / (1 - r2)
  }
  vifs
}
tryCatch({
  vifs <- vif_calc(best_model)
  print(round(vifs, 2))
  cat("\nVIF > 5 indicates multicollinearity concern\n")
}, error = function(e) cat("VIF calculation error:", e$message, "\n"))

cat("\n--- 4.4 Influential Observations ---\n")
cooks_d <- cooks.distance(best_model)
influential <- which(cooks_d > 4 / nrow(model_data))
cat("Observations with Cook's D > 4/n:", length(influential), "\n")
cat("Max Cook's D:", round(max(cooks_d), 4), "\n")
if (length(influential) > 0 && length(influential) <= 10) {
  cat("Influential observation indices:", head(influential, 10), "\n")
}

# ============================================================================
# SECTION 5: SUBGROUP ANALYSES
# ============================================================================

cat("\n
================================================================================
SECTION 5: SUBGROUP ANALYSES
================================================================================
\n")

cat("--- 5.1 R by Number of Studies (k) ---\n")
k_summary <- aggregate(R ~ k_cat, data = results,
                       FUN = function(x) c(n = length(x),
                                           mean = mean(x, na.rm = TRUE),
                                           median = median(x, na.rm = TRUE),
                                           sd = sd(x, na.rm = TRUE)))
k_summary <- do.call(data.frame, k_summary)
names(k_summary) <- c("k_category", "n", "mean_R", "median_R", "sd_R")
print(k_summary)

cat("\n--- 5.2 R by Effect Type ---\n")
et_summary <- aggregate(R ~ effect_type, data = results,
                        FUN = function(x) c(n = length(x),
                                            mean = mean(x, na.rm = TRUE),
                                            median = median(x, na.rm = TRUE),
                                            sd = sd(x, na.rm = TRUE)))
et_summary <- do.call(data.frame, et_summary)
names(et_summary) <- c("effect_type", "n", "mean_R", "median_R", "sd_R")
print(et_summary)

cat("\n--- 5.3 Kruskal-Wallis Test: R by effect_type ---\n")
kw_test <- kruskal.test(R ~ effect_type, data = results)
print(kw_test)

cat("\n--- 5.4 R by Heterogeneity Level ---\n")
results$het_level <- cut(results$tau,
                         breaks = c(-Inf, 1e-6, 0.1, 0.3, 0.5, Inf),
                         labels = c("None (tau~0)", "Low (<0.1)",
                                   "Moderate (0.1-0.3)", "High (0.3-0.5)",
                                   "Very High (>0.5)"))
het_summary <- aggregate(R ~ het_level, data = results,
                         FUN = function(x) c(n = length(x),
                                             mean = mean(x, na.rm = TRUE),
                                             median = median(x, na.rm = TRUE)))
het_summary <- do.call(data.frame, het_summary)
names(het_summary) <- c("heterogeneity", "n", "mean_R", "median_R")
print(het_summary)

cat("\n--- 5.5 R by Effect Size Magnitude ---\n")
results$theta_mag <- cut(abs(results$theta),
                         breaks = c(-Inf, 0.1, 0.3, 0.5, 1.0, Inf),
                         labels = c("Negligible (<0.1)", "Small (0.1-0.3)",
                                   "Medium (0.3-0.5)", "Large (0.5-1.0)",
                                   "Very Large (>1.0)"))
theta_summary <- aggregate(R ~ theta_mag, data = results[is.finite(results$theta), ],
                           FUN = function(x) c(n = length(x),
                                               mean = mean(x, na.rm = TRUE),
                                               median = median(x, na.rm = TRUE)))
theta_summary <- do.call(data.frame, theta_summary)
names(theta_summary) <- c("effect_magnitude", "n", "mean_R", "median_R")
print(theta_summary)

# ============================================================================
# SECTION 6: SIGN FLIP ANALYSIS
# ============================================================================

cat("\n
================================================================================
SECTION 6: SIGN FLIP ANALYSIS
================================================================================
\n")

# Identify potential sign flips (R <= 0.5 often indicates sign instability)
results$potential_sign_flip <- results$R <= 0.5

cat("--- 6.1 Sign Flip Prevalence ---\n")
cat("R <= 0.5 (potential sign instability):", sum(results$potential_sign_flip, na.rm = TRUE),
    "(", round(100 * mean(results$potential_sign_flip, na.rm = TRUE), 1), "%)\n\n")

cat("--- 6.2 Sign Flip by Effect Type ---\n")
sf_by_type <- aggregate(potential_sign_flip ~ effect_type, data = results,
                        FUN = function(x) c(n = sum(x), pct = round(100 * mean(x), 1)))
sf_by_type <- do.call(data.frame, sf_by_type)
names(sf_by_type) <- c("effect_type", "n_sign_flip", "pct_sign_flip")
print(sf_by_type)

cat("\n--- 6.3 Sign Flip by theta near zero ---\n")
results$theta_near_zero_cat <- ifelse(abs(results$theta) < 0.1, "Near zero", "Away from zero")
sf_by_theta <- aggregate(potential_sign_flip ~ theta_near_zero_cat, data = results[is.finite(results$theta), ],
                         FUN = function(x) c(n = sum(x), pct = round(100 * mean(x), 1)))
sf_by_theta <- do.call(data.frame, sf_by_theta)
names(sf_by_theta) <- c("theta_category", "n_sign_flip", "pct_sign_flip")
print(sf_by_theta)

cat("\n--- 6.4 Chi-square Test: Sign flip vs theta near zero ---\n")
results$theta_near_zero_factor <- factor(results$theta_near_zero_cat)
chi_test <- chisq.test(table(results$theta_near_zero_factor[is.finite(results$theta)],
                             results$potential_sign_flip[is.finite(results$theta)]))
print(chi_test)

# ============================================================================
# SECTION 7: PREDICTION ANALYSIS
# ============================================================================

cat("\n
================================================================================
SECTION 7: PREDICTION PERFORMANCE
================================================================================
\n")

# Cross-validation (simple holdout)
set.seed(42)
n <- nrow(model_data)
train_idx <- sample(1:n, floor(0.8 * n))
test_idx <- setdiff(1:n, train_idx)

train_data <- model_data[train_idx, ]
test_data <- model_data[test_idx, ]

cat("--- 7.1 Train/Test Split ---\n")
cat("Training set:", nrow(train_data), "observations\n")
cat("Test set:", nrow(test_data), "observations\n\n")

# Fit best model on training data
train_model <- lm(R ~ log_k + effect_type + log_tau + abs_theta + theta_near_zero +
                    log_k:effect_type, data = train_data)

# Predict on test data
pred_test <- predict(train_model, newdata = test_data)

cat("--- 7.2 Prediction Metrics ---\n")
# RMSE
rmse <- sqrt(mean((test_data$R - pred_test)^2))
cat("RMSE:", round(rmse, 4), "\n")

# MAE
mae <- mean(abs(test_data$R - pred_test))
cat("MAE:", round(mae, 4), "\n")

# R-squared on test set
ss_res <- sum((test_data$R - pred_test)^2)
ss_tot <- sum((test_data$R - mean(test_data$R))^2)
r2_test <- 1 - ss_res / ss_tot
cat("Test R-squared:", round(r2_test, 4), "\n")

# Correlation
cor_test <- cor(test_data$R, pred_test)
cat("Correlation (actual vs predicted):", round(cor_test, 4), "\n\n")

cat("--- 7.3 Prediction Accuracy by R Category ---\n")
test_data$pred_R <- pred_test
test_data$pred_R_cat <- cut(test_data$pred_R,
                            breaks = c(-Inf, 0.5, 0.7, 0.8, Inf),
                            labels = c("Low", "Moderate", "Good", "High"))
test_data$actual_R_cat <- cut(test_data$R,
                              breaks = c(-Inf, 0.5, 0.7, 0.8, Inf),
                              labels = c("Low", "Moderate", "Good", "High"))

confusion <- table(Predicted = test_data$pred_R_cat, Actual = test_data$actual_R_cat)
cat("Confusion Matrix:\n")
print(confusion)
cat("\nAccuracy:", round(sum(diag(confusion)) / sum(confusion), 3), "\n")

# ============================================================================
# SECTION 8: KEY FINDINGS SUMMARY
# ============================================================================

cat("\n
================================================================================
SECTION 8: KEY FINDINGS SUMMARY
================================================================================
\n")

cat("1. OVERALL STABILITY:\n")
cat("   - Mean R = ", round(mean(results$R, na.rm = TRUE), 3), "\n")
cat("   - Median R = ", round(median(results$R, na.rm = TRUE), 3), "\n")
cat("   - ", round(100 * mean(results$R >= 0.8, na.rm = TRUE), 1), "% of meta-analyses have high stability (R >= 0.8)\n")
cat("   - ", round(100 * mean(results$R < 0.5, na.rm = TRUE), 1), "% have potential sign instability (R < 0.5)\n\n")

cat("2. STRONGEST PREDICTORS OF STABILITY:\n")
coef_m5 <- summary(m5)$coefficients
sig_coefs <- coef_m5[coef_m5[, 4] < 0.05, , drop = FALSE]
cat("   Significant predictors (p < 0.05):\n")
for (i in 1:nrow(sig_coefs)) {
  direction <- ifelse(sig_coefs[i, 1] > 0, "increases", "decreases")
  cat("   - ", rownames(sig_coefs)[i], ": ", direction, " R (coef = ",
      round(sig_coefs[i, 1], 4), ")\n", sep = "")
}

cat("\n3. EFFECT TYPE DIFFERENCES:\n")
cat("   - logRR (binary): Mean R = ", round(mean(results$R[results$effect_type == "logRR"], na.rm = TRUE), 3), "\n")
cat("   - GIV (pre-computed): Mean R = ", round(mean(results$R[results$effect_type == "GIV"], na.rm = TRUE), 3), "\n")
cat("   - MD (continuous): Mean R = ", round(mean(results$R[results$effect_type == "MD"], na.rm = TRUE), 3), "\n\n")

cat("4. HETEROGENEITY IMPACT:\n")
cat("   - Zero heterogeneity: Mean R = ",
    round(mean(results$R[results$tau < 1e-6], na.rm = TRUE), 3), "\n")
cat("   - High heterogeneity (tau > 0.5): Mean R = ",
    round(mean(results$R[results$tau > 0.5], na.rm = TRUE), 3), "\n\n")

cat("5. SAMPLE SIZE (k) IMPACT:\n")
cat("   - k = 2: Mean R = ", round(mean(results$R[results$k == 2], na.rm = TRUE), 3), "\n")
cat("   - k = 3-5: Mean R = ", round(mean(results$R[results$k >= 3 & results$k <= 5], na.rm = TRUE), 3), "\n")
cat("   - k = 6-10: Mean R = ", round(mean(results$R[results$k >= 6 & results$k <= 10], na.rm = TRUE), 3), "\n")
cat("   - k > 20: Mean R = ", round(mean(results$R[results$k > 20], na.rm = TRUE), 3), "\n\n")

cat("6. MODEL PERFORMANCE:\n")
cat("   - Best model explains ", round(100 * summary(m5)$r.squared, 1), "% of variance in R\n")
cat("   - Test set RMSE = ", round(rmse, 4), "\n")
cat("   - Test set correlation = ", round(cor_test, 3), "\n\n")

cat("7. PRACTICAL RECOMMENDATIONS:\n")
cat("   - Meta-analyses with k < 5 have notably lower stability\n")
cat("   - GIV-type analyses (pre-computed effects) show more instability\n")
cat("   - Effects near zero (|theta| < 0.1) are more prone to sign flips\n")
cat("   - High heterogeneity (tau > 0.5) reduces stability\n")
cat("   - logRR analyses (binary outcomes) tend to be most stable\n")

# Save summary to file
sink("C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/ma4_modeling_summary.txt")
cat("MA4 v1.0.1 Modeling Summary - Pairwise70\n")
cat("=========================================\n\n")
cat("Total meta-analyses:", nrow(results), "\n")
cat("Mean R:", round(mean(results$R, na.rm = TRUE), 4), "\n")
cat("Median R:", round(median(results$R, na.rm = TRUE), 4), "\n")
cat("Best model R-squared:", round(summary(m5)$r.squared, 4), "\n")
cat("Test RMSE:", round(rmse, 4), "\n")
sink()

cat("\n
================================================================================
ANALYSIS COMPLETE
Results saved to: ma4_modeling_summary.txt
================================================================================
\n")
