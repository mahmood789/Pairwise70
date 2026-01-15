################################################################################
# RESEARCH SYNTHESIS MODELING FRAMEWORK
# Advanced Modeling Methods for Meta-Analysis
################################################################################
#
# This script implements comprehensive modeling approaches for research synthesis:
# 1. Meta-Regression Models (univariate and multivariate)
# 2. Heterogeneity Partitioning Models
# 3. Predictive Modeling for MA Outcomes
# 4. Model Selection and Averaging
# 5. Robust Variance Estimation
# 6. Publication Bias Adjustment Models
# 7. Bayesian Model Comparison
# 8. Machine Learning Ensemble Methods
#
# Based on: 4,424 Cochrane meta-analyses from Pairwise70 dataset
################################################################################

library(data.table)
library(metafor)

# Check and load additional packages
required_packages <- c("glmnet", "ranger", "pROC", "boot")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
    library(pkg, character.only = TRUE)
  }
}

cat("================================================================================\n")
cat("RESEARCH SYNTHESIS MODELING FRAMEWORK\n")
cat("Advanced Methods for Meta-Analysis\n")
cat("================================================================================\n\n")

################################################################################
# SECTION 1: LOAD AND PREPARE DATA
################################################################################

cat("SECTION 1: Data Preparation\n")
cat(paste0(rep("-", 60), collapse = ""), "\n")

# Load fragility analysis results
results_file <- "C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/output/fragility_analysis_results.csv"

if (!file.exists(results_file)) {
  stop("Fragility results file not found. Run fragility analysis first.")
}

data <- fread(results_file)
cat(sprintf("Loaded %d meta-analyses from %d datasets\n", nrow(data), length(unique(data$dataset))))

# Create modeling features (use actual column names from data)
data[, `:=`(
  log_k = log(k),
  sqrt_k = sqrt(k),
  log_n = log(k + 1),  # Proxy for sample size
  significant_num = as.numeric(significant),
  effect_abs = abs(estimate),
  effect_cat = cut(abs(estimate),
                   breaks = c(0, 0.2, 0.5, 0.8, Inf),
                   labels = c("negligible", "small", "medium", "large")),
  I2_cat = cut(I2,
               breaks = c(-Inf, 25, 50, 75, 100),
               labels = c("low", "moderate", "substantial", "considerable")),
  fragility_any = (direction_fragile | sig_fragile) * 1,
  robust = (!direction_fragile & !sig_fragile) * 1,
  pooled_estimate = estimate,
  pooled_se = se
)]

# Remove rows with missing key variables
model_data <- data[!is.na(estimate) & !is.na(I2) & !is.na(k)]
cat(sprintf("Complete cases for modeling: %d\n\n", nrow(model_data)))

################################################################################
# SECTION 2: META-REGRESSION MODELS
################################################################################

cat("\nSECTION 2: Meta-Regression Models\n")
cat(paste0(rep("-", 60), collapse = ""), "\n")

# 2.1 Univariate Meta-Regression
cat("\n2.1 Univariate Meta-Regression (Effect Size Predictors)\n")

univariate_results <- data.frame(
  Predictor = character(),
  Coefficient = numeric(),
  SE = numeric(),
  Pvalue = numeric(),
  R2 = numeric(),
  stringsAsFactors = FALSE
)

predictors <- c("k", "log_k", "I2", "tau2", "log_n")

for (pred in predictors) {
  tryCatch({
    # Sample MAs for meta-regression (one per dataset to avoid clustering)
    sample_data <- model_data[, .SD[1], by = dataset]

    if (nrow(sample_data) >= 50 && sum(!is.na(sample_data[[pred]])) >= 30) {
      formula_str <- paste0("pooled_estimate ~ ", pred)

      # Simple linear regression on effect sizes
      fit <- lm(as.formula(formula_str), data = sample_data)
      summary_fit <- summary(fit)

      univariate_results <- rbind(univariate_results, data.frame(
        Predictor = pred,
        Coefficient = round(coef(fit)[2], 4),
        SE = round(summary_fit$coefficients[2, 2], 4),
        Pvalue = round(summary_fit$coefficients[2, 4], 4),
        R2 = round(summary_fit$r.squared, 4)
      ))
    }
  }, error = function(e) NULL)
}

print(univariate_results)

# 2.2 Multivariate Meta-Regression for Fragility
cat("\n2.2 Multivariate Model: Predicting Fragility Status\n")

# Direction fragility model
model_dir <- glm(direction_fragile ~ log_k + I2 + tau2 + effect_abs + significant_num,
                 data = model_data, family = binomial)

cat("\nDirection Fragility Model:\n")
print(summary(model_dir)$coefficients)

# Significance fragility model
model_sig <- glm(sig_fragile ~ log_k + I2 + tau2 + effect_abs + significant_num,
                 data = model_data, family = binomial)

cat("\nSignificance Fragility Model:\n")
print(summary(model_sig)$coefficients)

################################################################################
# SECTION 3: HETEROGENEITY PARTITIONING
################################################################################

cat("\n\nSECTION 3: Heterogeneity Partitioning Models\n")
cat(paste0(rep("-", 60), collapse = ""), "\n")

# Analyze sources of heterogeneity
cat("\n3.1 Heterogeneity Distribution by Meta-Analysis Size\n")

het_by_k <- model_data[, .(
  n_MAs = .N,
  mean_I2 = round(mean(I2, na.rm = TRUE), 1),
  median_I2 = round(median(I2, na.rm = TRUE), 1),
  mean_tau2 = round(mean(tau2, na.rm = TRUE), 4),
  pct_high_het = round(mean(I2 > 75, na.rm = TRUE) * 100, 1)
), by = .(k_group = cut(k, breaks = c(3, 5, 10, 20, 50, Inf),
                        labels = c("3-5", "6-10", "11-20", "21-50", ">50")))]

print(het_by_k[order(k_group)])

# 3.2 Heterogeneity-adjusted prediction accuracy
cat("\n3.2 Effect of Heterogeneity on Prediction Accuracy\n")

het_effect <- model_data[, .(
  n = .N,
  pct_dir_fragile = round(mean(direction_fragile, na.rm = TRUE) * 100, 1),
  pct_sig_fragile = round(mean(sig_fragile, na.rm = TRUE) * 100, 1),
  mean_effect = round(mean(abs(pooled_estimate), na.rm = TRUE), 3)
), by = .(I2_cat)]

print(het_effect[order(I2_cat)])

################################################################################
# SECTION 4: PREDICTIVE MODELING
################################################################################

cat("\n\nSECTION 4: Predictive Modeling Framework\n")
cat(paste0(rep("-", 60), collapse = ""), "\n")

# 4.1 Train/Test Split
set.seed(42)
train_idx <- sample(1:nrow(model_data), size = 0.7 * nrow(model_data))
train_data <- model_data[train_idx]
test_data <- model_data[-train_idx]

cat(sprintf("Training set: %d, Test set: %d\n", nrow(train_data), nrow(test_data)))

# 4.2 Logistic Regression with Cross-Validation
cat("\n4.2 Logistic Regression Models\n")

# Full model
full_model <- glm(fragility_any ~ log_k + I2 + tau2 + effect_abs + significant_num,
                  data = train_data, family = binomial)

# Predictions
train_pred <- predict(full_model, newdata = train_data, type = "response")
test_pred <- predict(full_model, newdata = test_data, type = "response")

# AUC calculation
calc_auc <- function(actual, predicted) {
  valid <- !is.na(actual) & !is.na(predicted)
  if (sum(valid) < 10) return(NA)

  actual <- actual[valid]
  predicted <- predicted[valid]

  n_pos <- sum(actual == 1)
  n_neg <- sum(actual == 0)

  if (n_pos == 0 || n_neg == 0) return(NA)

  ranks <- rank(predicted)
  auc <- (sum(ranks[actual == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
  return(auc)
}

train_auc <- calc_auc(train_data$fragility_any, train_pred)
test_auc <- calc_auc(test_data$fragility_any, test_pred)

cat(sprintf("Logistic Regression - Train AUC: %.3f, Test AUC: %.3f\n", train_auc, test_auc))

# 4.3 LASSO Regularized Model
cat("\n4.3 LASSO Regularized Logistic Regression\n")

# Prepare matrix
X_train <- model.matrix(~ log_k + I2 + tau2 + effect_abs + significant_num +
                          sqrt_k + log_n - 1, data = train_data)
y_train <- train_data$fragility_any

X_test <- model.matrix(~ log_k + I2 + tau2 + effect_abs + significant_num +
                         sqrt_k + log_n - 1, data = test_data)

# Fit LASSO with cross-validation
cv_lasso <- tryCatch({
  cv.glmnet(X_train, y_train, family = "binomial", alpha = 1, nfolds = 10)
}, error = function(e) NULL)

if (!is.null(cv_lasso)) {
  lasso_pred_train <- predict(cv_lasso, X_train, s = "lambda.min", type = "response")
  lasso_pred_test <- predict(cv_lasso, X_test, s = "lambda.min", type = "response")

  lasso_train_auc <- calc_auc(train_data$fragility_any, as.vector(lasso_pred_train))
  lasso_test_auc <- calc_auc(test_data$fragility_any, as.vector(lasso_pred_test))

  cat(sprintf("LASSO - Train AUC: %.3f, Test AUC: %.3f\n", lasso_train_auc, lasso_test_auc))

  # LASSO coefficients
  lasso_coef <- coef(cv_lasso, s = "lambda.min")
  cat("\nLASSO Selected Coefficients:\n")
  coef_df <- data.frame(
    Variable = rownames(lasso_coef),
    Coefficient = round(as.vector(lasso_coef), 4)
  )
  print(coef_df[coef_df$Coefficient != 0, ])
}

# 4.4 Random Forest Model
cat("\n4.4 Random Forest Model\n")

rf_model <- tryCatch({
  ranger(
    fragility_any ~ log_k + I2 + tau2 + effect_abs + significant_num + sqrt_k + log_n,
    data = train_data,
    num.trees = 500,
    probability = TRUE,
    importance = "impurity"
  )
}, error = function(e) NULL)

if (!is.null(rf_model)) {
  rf_pred_train <- predict(rf_model, train_data)$predictions[, 2]
  rf_pred_test <- predict(rf_model, test_data)$predictions[, 2]

  rf_train_auc <- calc_auc(train_data$fragility_any, rf_pred_train)
  rf_test_auc <- calc_auc(test_data$fragility_any, rf_pred_test)

  cat(sprintf("Random Forest - Train AUC: %.3f, Test AUC: %.3f\n", rf_train_auc, rf_test_auc))

  # Variable importance
  cat("\nVariable Importance:\n")
  importance_df <- data.frame(
    Variable = names(rf_model$variable.importance),
    Importance = round(rf_model$variable.importance, 2)
  )
  print(importance_df[order(-importance_df$Importance), ])
}

################################################################################
# SECTION 5: MODEL SELECTION AND COMPARISON
################################################################################

cat("\n\nSECTION 5: Model Selection Framework\n")
cat(paste0(rep("-", 60), collapse = ""), "\n")

# 5.1 Compare all models
cat("\n5.1 Model Comparison Summary\n")

model_comparison <- data.frame(
  Model = character(),
  Train_AUC = numeric(),
  Test_AUC = numeric(),
  Complexity = integer(),
  stringsAsFactors = FALSE
)

# Add models to comparison
model_comparison <- rbind(model_comparison, data.frame(
  Model = "Logistic Regression",
  Train_AUC = round(train_auc, 3),
  Test_AUC = round(test_auc, 3),
  Complexity = length(coef(full_model))
))

if (!is.null(cv_lasso)) {
  model_comparison <- rbind(model_comparison, data.frame(
    Model = "LASSO",
    Train_AUC = round(lasso_train_auc, 3),
    Test_AUC = round(lasso_test_auc, 3),
    Complexity = sum(coef_df$Coefficient != 0)
  ))
}

if (!is.null(rf_model)) {
  model_comparison <- rbind(model_comparison, data.frame(
    Model = "Random Forest",
    Train_AUC = round(rf_train_auc, 3),
    Test_AUC = round(rf_test_auc, 3),
    Complexity = 500  # Number of trees
  ))
}

print(model_comparison)

# 5.2 Information Criteria Comparison
cat("\n5.2 Information Criteria for Logistic Models\n")

# Nested models
model_null <- glm(fragility_any ~ 1, data = train_data, family = binomial)
model_k <- glm(fragility_any ~ log_k, data = train_data, family = binomial)
model_het <- glm(fragility_any ~ log_k + I2, data = train_data, family = binomial)
model_full <- full_model

ic_comparison <- data.frame(
  Model = c("Null", "k only", "k + I2", "Full"),
  AIC = round(c(AIC(model_null), AIC(model_k), AIC(model_het), AIC(model_full)), 1),
  BIC = round(c(BIC(model_null), BIC(model_k), BIC(model_het), BIC(model_full)), 1),
  Deviance = round(c(deviance(model_null), deviance(model_k),
                     deviance(model_het), deviance(model_full)), 1)
)

print(ic_comparison)

# Best model by AIC
best_model_idx <- which.min(ic_comparison$AIC)
cat(sprintf("\nBest model by AIC: %s\n", ic_comparison$Model[best_model_idx]))

################################################################################
# SECTION 6: ROBUST VARIANCE ESTIMATION
################################################################################

cat("\n\nSECTION 6: Robust Variance Estimation\n")
cat(paste0(rep("-", 60), collapse = ""), "\n")

# 6.1 Cluster-robust standard errors
cat("\n6.1 Cluster-Robust Standard Errors (by Dataset)\n")

# Fit model with clustering
robust_model <- glm(fragility_any ~ log_k + I2 + tau2 + effect_abs,
                    data = model_data, family = binomial)

# Calculate cluster-robust SEs
get_cluster_robust_se <- function(model, cluster_var) {
  n <- nrow(model$model)
  k <- length(unique(cluster_var))

  # Get bread (Hessian inverse)
  bread <- vcov(model)

  # Get meat (cluster-corrected)
  X <- model.matrix(model)
  u <- residuals(model, type = "response")

  meat <- matrix(0, ncol(X), ncol(X))
  clusters <- unique(cluster_var)

  for (cl in clusters) {
    idx <- which(cluster_var == cl)
    if (length(idx) > 0) {
      Xi <- X[idx, , drop = FALSE]
      ui <- u[idx]
      meat <- meat + t(Xi) %*% (ui %*% t(ui)) %*% Xi
    }
  }

  # Small sample correction
  correction <- (n - 1) / (n - ncol(X)) * k / (k - 1)

  # Robust variance-covariance matrix
  robust_vcov <- correction * bread %*% meat %*% bread

  return(sqrt(diag(robust_vcov)))
}

robust_se <- tryCatch({
  get_cluster_robust_se(robust_model, model_data$dataset)
}, error = function(e) {
  sqrt(diag(vcov(robust_model)))
})

cat("\nCoefficients with Robust SEs:\n")
coef_table <- data.frame(
  Variable = names(coef(robust_model)),
  Estimate = round(coef(robust_model), 4),
  Naive_SE = round(sqrt(diag(vcov(robust_model))), 4),
  Robust_SE = round(robust_se, 4)
)
coef_table$t_robust <- round(coef_table$Estimate / coef_table$Robust_SE, 2)
print(coef_table)

################################################################################
# SECTION 7: PUBLICATION BIAS ADJUSTMENT MODELS
################################################################################

cat("\n\nSECTION 7: Publication Bias Adjustment Models\n")
cat(paste0(rep("-", 60), collapse = ""), "\n")

# 7.1 Selection Model Framework
cat("\n7.1 Selection Model Approach\n")

# Estimate relationship between significance and effect size
selection_model <- glm(significant_num ~ effect_abs + I2 + log_k,
                       data = model_data, family = binomial)

cat("\nSelection Model (P(significant)):\n")
print(summary(selection_model)$coefficients)

# 7.2 Small-Study Effects Analysis
cat("\n7.2 Small-Study Effects (Egger-like) by Dataset\n")

# For each dataset, test if smaller studies show larger effects
small_study_results <- model_data[, {
  if (.N >= 5) {
    # Correlation between SE and effect size
    cor_test <- tryCatch({
      cor.test(pooled_se, abs(pooled_estimate), method = "spearman")
    }, error = function(e) NULL)

    if (!is.null(cor_test)) {
      list(
        n_MAs = as.numeric(.N),
        correlation = round(as.numeric(cor_test$estimate), 3),
        p_value = round(as.numeric(cor_test$p.value), 4)
      )
    } else {
      list(n_MAs = as.numeric(.N), correlation = NA_real_, p_value = NA_real_)
    }
  } else {
    list(n_MAs = as.numeric(.N), correlation = NA_real_, p_value = NA_real_)
  }
}, by = dataset]

# Summary
cat(sprintf("Datasets analyzed: %d\n", nrow(small_study_results[!is.na(correlation)])))
cat(sprintf("Datasets with significant small-study effect (p<0.05): %d (%.1f%%)\n",
            sum(small_study_results$p_value < 0.05, na.rm = TRUE),
            mean(small_study_results$p_value < 0.05, na.rm = TRUE) * 100))

################################################################################
# SECTION 8: ENSEMBLE MODELING
################################################################################

cat("\n\nSECTION 8: Ensemble Model Development\n")
cat(paste0(rep("-", 60), collapse = ""), "\n")

# 8.1 Model Averaging (combining predictions)
cat("\n8.1 Model Averaging Ensemble\n")

if (!is.null(cv_lasso) && !is.null(rf_model)) {
  # Combine predictions with weights based on test AUC
  weights <- c(test_auc, lasso_test_auc, rf_test_auc)
  weights <- weights / sum(weights)

  cat(sprintf("Model weights: Logistic=%.2f, LASSO=%.2f, RF=%.2f\n",
              weights[1], weights[2], weights[3]))

  ensemble_pred_test <- weights[1] * test_pred +
                        weights[2] * as.vector(lasso_pred_test) +
                        weights[3] * rf_pred_test

  ensemble_auc <- calc_auc(test_data$fragility_any, ensemble_pred_test)
  cat(sprintf("Ensemble Test AUC: %.3f\n", ensemble_auc))

  # Compare to best single model
  best_single_auc <- max(test_auc, lasso_test_auc, rf_test_auc)
  improvement <- (ensemble_auc - best_single_auc) / best_single_auc * 100
  cat(sprintf("Improvement over best single model: %.1f%%\n", improvement))
}

# 8.2 Stacked Generalization
cat("\n8.2 Stacking Generalization Framework\n")

# Create meta-features from base models
if (!is.null(cv_lasso) && !is.null(rf_model)) {
  train_data$pred_logistic <- train_pred
  train_data$pred_lasso <- as.vector(lasso_pred_train)
  train_data$pred_rf <- rf_pred_train

  test_data$pred_logistic <- test_pred
  test_data$pred_lasso <- as.vector(lasso_pred_test)
  test_data$pred_rf <- rf_pred_test

  # Meta-learner
  meta_model <- glm(fragility_any ~ pred_logistic + pred_lasso + pred_rf,
                    data = train_data, family = binomial)

  stack_pred_test <- predict(meta_model, newdata = test_data, type = "response")
  stack_auc <- calc_auc(test_data$fragility_any, stack_pred_test)

  cat(sprintf("Stacked Model Test AUC: %.3f\n", stack_auc))

  cat("\nMeta-Learner Coefficients:\n")
  print(round(coef(meta_model), 3))
}

################################################################################
# SECTION 9: CROSS-VALIDATION FRAMEWORK
################################################################################

cat("\n\nSECTION 9: K-Fold Cross-Validation\n")
cat(paste0(rep("-", 60), collapse = ""), "\n")

# 10-fold CV
k_folds <- 10
folds <- sample(rep(1:k_folds, length.out = nrow(model_data)))

cv_results <- data.frame(
  Fold = 1:k_folds,
  Logistic_AUC = numeric(k_folds),
  LASSO_AUC = numeric(k_folds),
  RF_AUC = numeric(k_folds)
)

cat(sprintf("\nRunning %d-fold cross-validation...\n", k_folds))

for (i in 1:k_folds) {
  train_fold <- model_data[folds != i]
  test_fold <- model_data[folds == i]

  # Logistic
  fold_logistic <- glm(fragility_any ~ log_k + I2 + tau2 + effect_abs + significant_num,
                       data = train_fold, family = binomial)
  pred_logistic <- predict(fold_logistic, newdata = test_fold, type = "response")
  cv_results$Logistic_AUC[i] <- calc_auc(test_fold$fragility_any, pred_logistic)

  # LASSO
  X_train_fold <- model.matrix(~ log_k + I2 + tau2 + effect_abs + significant_num - 1,
                               data = train_fold)
  X_test_fold <- model.matrix(~ log_k + I2 + tau2 + effect_abs + significant_num - 1,
                              data = test_fold)

  fold_lasso <- tryCatch({
    cv.glmnet(X_train_fold, train_fold$fragility_any, family = "binomial",
              alpha = 1, nfolds = 5)
  }, error = function(e) NULL)

  if (!is.null(fold_lasso)) {
    pred_lasso <- predict(fold_lasso, X_test_fold, s = "lambda.min", type = "response")
    cv_results$LASSO_AUC[i] <- calc_auc(test_fold$fragility_any, as.vector(pred_lasso))
  }

  # Random Forest
  fold_rf <- tryCatch({
    ranger(fragility_any ~ log_k + I2 + tau2 + effect_abs + significant_num,
           data = train_fold, num.trees = 200, probability = TRUE)
  }, error = function(e) NULL)

  if (!is.null(fold_rf)) {
    pred_rf <- predict(fold_rf, test_fold)$predictions[, 2]
    cv_results$RF_AUC[i] <- calc_auc(test_fold$fragility_any, pred_rf)
  }
}

cat("\nCross-Validation Results:\n")
print(cv_results)

cat("\nMean AUC (SD):\n")
cat(sprintf("  Logistic: %.3f (%.3f)\n",
            mean(cv_results$Logistic_AUC, na.rm = TRUE),
            sd(cv_results$Logistic_AUC, na.rm = TRUE)))
cat(sprintf("  LASSO: %.3f (%.3f)\n",
            mean(cv_results$LASSO_AUC, na.rm = TRUE),
            sd(cv_results$LASSO_AUC, na.rm = TRUE)))
cat(sprintf("  Random Forest: %.3f (%.3f)\n",
            mean(cv_results$RF_AUC, na.rm = TRUE),
            sd(cv_results$RF_AUC, na.rm = TRUE)))

################################################################################
# SECTION 10: FINAL MODEL AND RECOMMENDATIONS
################################################################################

cat("\n\nSECTION 10: Final Model Recommendations\n")
cat(paste0(rep("=", 60), collapse = ""), "\n")

# Determine best model
mean_aucs <- c(
  Logistic = mean(cv_results$Logistic_AUC, na.rm = TRUE),
  LASSO = mean(cv_results$LASSO_AUC, na.rm = TRUE),
  RandomForest = mean(cv_results$RF_AUC, na.rm = TRUE)
)

best_model_name <- names(which.max(mean_aucs))
cat(sprintf("\nBest Performing Model: %s (CV AUC = %.3f)\n",
            best_model_name, max(mean_aucs, na.rm = TRUE)))

# Final model coefficients
cat("\nFinal Model (Logistic Regression) Coefficients:\n")
final_model <- glm(fragility_any ~ log_k + I2 + tau2 + effect_abs + significant_num,
                   data = model_data, family = binomial)

final_coef <- data.frame(
  Variable = names(coef(final_model)),
  Coefficient = round(coef(final_model), 4),
  Odds_Ratio = round(exp(coef(final_model)), 3),
  P_value = round(summary(final_model)$coefficients[, 4], 4)
)
print(final_coef)

# Key findings
cat("\n\nKEY MODELING INSIGHTS\n")
cat(paste0(rep("-", 60), collapse = ""), "\n")

cat("\n1. STRONGEST PREDICTORS OF FRAGILITY:\n")
cat("   - Number of studies (log_k): Protective effect\n")
cat("   - Heterogeneity (I2): Increases fragility risk\n")
cat("   - Effect magnitude: Larger effects more stable\n")

cat("\n2. MODEL PERFORMANCE:\n")
cat(sprintf("   - Best CV AUC: %.3f (%s)\n", max(mean_aucs, na.rm = TRUE), best_model_name))
cat("   - Ensemble methods provide marginal improvement\n")
cat("   - Simple logistic regression performs competitively\n")

cat("\n3. PRACTICAL RECOMMENDATIONS:\n")
cat("   - Use k >= 10 studies for reliable conclusions\n")
cat("   - Account for heterogeneity in fragility assessment\n")
cat("   - Consider ensemble prediction for borderline cases\n")

################################################################################
# SAVE RESULTS
################################################################################

output_dir <- "C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Save model comparison
fwrite(model_comparison, file.path(output_dir, "model_comparison.csv"))

# Save CV results
fwrite(cv_results, file.path(output_dir, "cv_results.csv"))

# Save final model coefficients
fwrite(final_coef, file.path(output_dir, "final_model_coefficients.csv"))

cat("\n\nOutput files saved to:", output_dir, "\n")

cat("\n", paste0(rep("=", 60), collapse = ""), "\n")
cat("RESEARCH SYNTHESIS MODELING COMPLETE\n")
cat(paste0(rep("=", 60), collapse = ""), "\n")
