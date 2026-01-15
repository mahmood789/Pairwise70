################################################################################
#     EDITORIAL REVISION: Research Synthesis Methods (FAST VERSION)
#     ===========================================================
#     Focuses on core statistical improvements using MA4 data directly
################################################################################

cat("\n================================================================\n")
cat("EDITORIAL REVISION: Research Synthesis Methods\n")
cat("================================================================\n\n")

setwd("C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis")
dir.create("output", showWarnings = FALSE)

# Load packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lme4)
  library(pROC)
  library(caret)
  library(randomForest)
})

# Load data
ma4_data <- read.csv("ma4_results_pairwise70.csv", stringsAsFactors = FALSE)
cat(paste0("Loaded ", nrow(ma4_data), " meta-analyses\n\n"))

################################################################################
# REVISION 1: IMPROVED DOMAIN CLASSIFICATION
################################################################################

cat("=== REVISION 1: IMPROVED DOMAIN CLASSIFICATION ===\n\n")

classify_by_outcome <- function(analysis_name) {
  name_lower <- tolower(analysis_name)

  if (grepl("death|mortal|surviv|fatal", name_lower)) return("Mortality")
  if (grepl("adverse|side effect|withdraw|dropout|tolerab|safety|toxicity", name_lower)) return("Adverse Events")
  if (grepl("pain|analges|ache", name_lower)) return("Pain Outcomes")
  if (grepl("infection|sepsis|fever|pneumonia", name_lower)) return("Infection")
  if (grepl("bleed|hemorrh|transfus|haemorrhag", name_lower)) return("Bleeding")
  if (grepl("fractur|bone|osteo", name_lower)) return("Musculoskeletal")
  if (grepl("heart|cardiac|coronary|arrhythm|myocard|infarct", name_lower)) return("Cardiovascular")
  if (grepl("stroke|cerebr", name_lower)) return("Stroke")
  if (grepl("cancer|tumor|tumour|malignan|oncol|chemotherapy|carcinoma", name_lower)) return("Cancer")
  if (grepl("diabet|glucose|insulin|glyc", name_lower)) return("Diabetes/Metabolic")
  if (grepl("pregnan|birth|deliver|caesar|neonat|fetal|maternal", name_lower)) return("Pregnancy/Birth")
  if (grepl("child|infant|pediatr|paediatr", name_lower)) return("Pediatric")
  if (grepl("depress|anxiety|psychiat|mental|schizo|bipolar", name_lower)) return("Mental Health")
  if (grepl("respirat|lung|asthma|copd|pulmon|ventilat|oxygen", name_lower)) return("Respiratory")
  if (grepl("renal|kidney|dialysis|urinar", name_lower)) return("Renal")
  if (grepl("liver|hepat|gastro|bowel|intestin|colitis", name_lower)) return("GI/Hepatic")
  if (grepl("surgery|surgical|operat|incision|wound", name_lower)) return("Surgical")
  if (grepl("quality of life|qol|function|disabil", name_lower)) return("Quality of Life")
  if (grepl("length of stay|hospital|admission|discharge", name_lower)) return("Healthcare Utilization")
  if (grepl("recurrence|relapse|remission", name_lower)) return("Disease Course")
  if (grepl("response|improvement|efficacy|success", name_lower)) return("Treatment Response")

  return("Other Clinical")
}

ma4_data$outcome_domain <- sapply(ma4_data$analysis_name, classify_by_outcome)

domain_counts <- sort(table(ma4_data$outcome_domain), decreasing = TRUE)
cat("Domain Classification (22 categories):\n")
print(domain_counts)

other_pct <- round(sum(ma4_data$outcome_domain == "Other Clinical") / nrow(ma4_data) * 100, 1)
cat(paste0("\n'Other Clinical' category: ", other_pct, "%\n"))

################################################################################
# REVISION 2: EFFECT TYPE SEPARATION
################################################################################

cat("\n\n=== REVISION 2: EFFECT TYPE SEPARATION ===\n\n")

logRR_data <- ma4_data %>% filter(effect_type == "logRR")
GIV_data <- ma4_data %>% filter(effect_type == "GIV")
MD_data <- ma4_data %>% filter(effect_type == "MD")

cat("Effect Type Distribution:\n")
cat(paste0("  logRR (binary): ", nrow(logRR_data), " (", round(nrow(logRR_data)/nrow(ma4_data)*100, 1), "%)\n"))
cat(paste0("  GIV (generic): ", nrow(GIV_data), " (", round(nrow(GIV_data)/nrow(ma4_data)*100, 1), "%)\n"))
cat(paste0("  MD (continuous): ", nrow(MD_data), " (", round(nrow(MD_data)/nrow(ma4_data)*100, 1), "%)\n"))

cat("\n*** PRIMARY ANALYSIS: logRR DATA (n=", nrow(logRR_data), ") ***\n")
primary_data <- logRR_data

################################################################################
# REVISION 3: FRAGILITY METRIC INTERPRETATION
################################################################################

cat("\n\n=== REVISION 3: R STABILITY METRIC ANALYSIS ===\n\n")

# R metric is already available from MA4 analysis
# R = 1 - max(p_i) where p_i is the influence of each study

cat("Stability (R) Distribution:\n")
cat(paste0("  Mean R: ", round(mean(primary_data$R, na.rm = TRUE), 3), "\n"))
cat(paste0("  Median R: ", round(median(primary_data$R, na.rm = TRUE), 3), "\n"))
cat(paste0("  SD: ", round(sd(primary_data$R, na.rm = TRUE), 3), "\n\n"))

# Categorize
primary_data$R_category <- cut(primary_data$R,
                                breaks = c(0, 0.5, 0.8, 1.0),
                                labels = c("Fragile (R<0.5)", "Moderate (0.5-0.8)", "Stable (R>0.8)"),
                                include.lowest = TRUE)

cat("Stability Categories:\n")
r_cat_table <- table(primary_data$R_category)
print(r_cat_table)
cat(paste0("\nHighly Fragile: ", round(r_cat_table[1] / sum(r_cat_table) * 100, 1), "%\n"))
cat(paste0("Stable: ", round(r_cat_table[3] / sum(r_cat_table) * 100, 1), "%\n"))

################################################################################
# REVISION 4: PREDICTIVE MODEL WITH CV, AUC, CLASS WEIGHTS
################################################################################

cat("\n\n=== REVISION 4: IMPROVED PREDICTIVE MODEL ===\n\n")

model_data <- primary_data %>%
  filter(!is.na(R) & !is.na(tau) & !is.na(sigma) & k >= 3) %>%
  mutate(
    log_k = log(k),
    abs_theta = abs(theta),
    log_tau = log(tau + 0.001),
    log_sigma = log(sigma + 0.001),
    theta_near_zero = as.numeric(abs(theta) < 0.1),
    fragile = factor(ifelse(R < 0.5, "Fragile", "Stable"), levels = c("Stable", "Fragile"))
  ) %>%
  filter(is.finite(log_tau) & is.finite(log_sigma))

cat(paste0("Modeling sample: ", nrow(model_data), " meta-analyses\n"))
cat("Class distribution:\n")
print(table(model_data$fragile))

imbalance <- sum(model_data$fragile == "Stable") / sum(model_data$fragile == "Fragile")
cat(paste0("Class imbalance: ", round(imbalance, 1), ":1 (Stable:Fragile)\n\n"))

features <- c("log_k", "abs_theta", "log_tau", "log_sigma", "theta_near_zero")

# 10-fold stratified CV
set.seed(42)
folds <- createFolds(model_data$fragile, k = 10, list = TRUE)

cv_results <- data.frame()
all_predictions <- data.frame()

cat("Running 10-fold stratified cross-validation...\n")

for (fold_idx in 1:length(folds)) {
  test_idx <- folds[[fold_idx]]
  train_idx <- setdiff(1:nrow(model_data), test_idx)

  train <- model_data[train_idx, ]
  test <- model_data[test_idx, ]

  # Class weight of 4x for minority class
  rf <- randomForest(
    x = train[, features],
    y = train$fragile,
    ntree = 500,
    classwt = c("Stable" = 1, "Fragile" = 4),
    importance = TRUE
  )

  pred_class <- predict(rf, test[, features])
  pred_prob <- predict(rf, test[, features], type = "prob")[, "Fragile"]

  all_predictions <- rbind(all_predictions, data.frame(
    actual = test$fragile,
    predicted = pred_class,
    prob_fragile = pred_prob,
    fold = fold_idx
  ))

  cm <- confusionMatrix(pred_class, test$fragile, positive = "Fragile")

  cv_results <- rbind(cv_results, data.frame(
    fold = fold_idx,
    accuracy = cm$overall["Accuracy"],
    sensitivity = cm$byClass["Sensitivity"],
    specificity = cm$byClass["Specificity"],
    balanced_accuracy = cm$byClass["Balanced Accuracy"]
  ))
}

# AUC-ROC
roc_obj <- roc(all_predictions$actual, all_predictions$prob_fragile,
               levels = c("Stable", "Fragile"), direction = "<")
auc_value <- auc(roc_obj)
auc_ci <- ci.auc(roc_obj)

cat("\nCROSS-VALIDATION RESULTS:\n")
cat(paste(rep("-", 60), collapse = ""), "\n")
cat(paste0("AUC-ROC: ", round(auc_value, 3), " (95% CI: ", round(auc_ci[1], 3), "-", round(auc_ci[3], 3), ")\n\n"))

cv_summary <- colMeans(cv_results[, -1], na.rm = TRUE)
cat("Mean performance across 10 folds:\n")
cat(paste0("  Accuracy: ", round(cv_summary["accuracy"], 3), "\n"))
cat(paste0("  Sensitivity (Fragile detection): ", round(cv_summary["sensitivity"], 3), "\n"))
cat(paste0("  Specificity (Stable detection): ", round(cv_summary["specificity"], 3), "\n"))
cat(paste0("  Balanced Accuracy: ", round(cv_summary["balanced_accuracy"], 3), "\n"))

# Improvement over null
null_accuracy <- max(table(model_data$fragile)) / nrow(model_data)
cat(paste0("\nNull model (always Stable): ", round(null_accuracy, 3), "\n"))
cat(paste0("Improvement: +", round((cv_summary["balanced_accuracy"] - 0.5) * 100, 1), " ppts balanced acc\n"))

# Final model for feature importance
final_rf <- randomForest(
  x = model_data[, features],
  y = model_data$fragile,
  ntree = 500,
  classwt = c("Stable" = 1, "Fragile" = 4),
  importance = TRUE
)

importance_df <- data.frame(
  Feature = rownames(importance(final_rf)),
  MeanDecreaseGini = importance(final_rf)[, "MeanDecreaseGini"],
  MeanDecreaseAccuracy = importance(final_rf)[, "MeanDecreaseAccuracy"]
) %>% arrange(desc(MeanDecreaseGini))

cat("\nFEATURE IMPORTANCE (Mean Decrease Gini):\n")
print(importance_df, row.names = FALSE)

saveRDS(list(
  cv_results = cv_results,
  cv_summary = cv_summary,
  auc = auc_value,
  auc_ci = auc_ci,
  importance = importance_df
), "output/predictive_model_revised.rds")

################################################################################
# REVISION 5: MIXED-EFFECTS MODELS FOR CLUSTERING
################################################################################

cat("\n\n=== REVISION 5: MIXED-EFFECTS MODELS ===\n\n")

# Extract review base (removing _pubN suffix)
model_data$review_base <- gsub("_pub[0-9]+$", "", model_data$review_id)

n_reviews <- length(unique(model_data$review_base))
n_ma <- nrow(model_data)
cat(paste0("Unique reviews: ", n_reviews, "\n"))
cat(paste0("Meta-analyses: ", n_ma, "\n"))
cat(paste0("MAs per review: ", round(n_ma / n_reviews, 1), " (mean)\n\n"))

cat("Fitting mixed-effects model for stability (R)...\n")

me_model <- lmer(R ~ log_k + abs_theta + log_tau + log_sigma + theta_near_zero +
                   (1 | review_base),
                 data = model_data)

cat("\nMIXED-EFFECTS MODEL RESULTS:\n")
cat(paste(rep("-", 60), collapse = ""), "\n")

# Fixed effects
cat("\nFixed Effects:\n")
fe <- summary(me_model)$coefficients
print(round(fe, 4))

# Variance components
vc <- as.data.frame(VarCorr(me_model))
review_var <- vc$vcov[1]
residual_var <- vc$vcov[2]
icc <- review_var / (review_var + residual_var)

cat("\nVariance Components:\n")
cat(paste0("  Between-review variance: ", round(review_var, 4), "\n"))
cat(paste0("  Residual variance: ", round(residual_var, 4), "\n"))
cat(paste0("\nIntraclass Correlation (ICC): ", round(icc, 3), "\n"))
cat(paste0("  Interpretation: ", round(icc * 100, 1), "% of R variance is between reviews\n"))

# Model comparison
ols_model <- lm(R ~ log_k + abs_theta + log_tau + log_sigma + theta_near_zero, data = model_data)
aic_ols <- AIC(ols_model)
aic_me <- AIC(me_model)

cat("\nModel Comparison (AIC, lower = better):\n")
cat(paste0("  OLS (ignoring clustering): ", round(aic_ols, 1), "\n"))
cat(paste0("  Mixed-Effects: ", round(aic_me, 1), "\n"))
cat(paste0("  Improvement: ", round(aic_ols - aic_me, 1), " points\n"))

saveRDS(me_model, "output/mixed_effects_model.rds")

################################################################################
# REVISION 6: DOMAIN ANALYSIS WITH FDR CORRECTION
################################################################################

cat("\n\n=== REVISION 6: DOMAIN ANALYSIS WITH FDR CORRECTION ===\n\n")

# Domain-level summary
domain_analysis <- primary_data %>%
  group_by(outcome_domain) %>%
  summarise(
    n = n(),
    n_reviews = n_distinct(review_id),
    mean_R = mean(R, na.rm = TRUE),
    sd_R = sd(R, na.rm = TRUE),
    se_R = sd(R, na.rm = TRUE) / sqrt(n()),
    pct_fragile = mean(R < 0.5, na.rm = TRUE) * 100,
    mean_tau = mean(tau, na.rm = TRUE),
    mean_k = mean(k, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n >= 10)

overall_mean_R <- mean(primary_data$R, na.rm = TRUE)

# One-sample t-tests vs overall mean
domain_tests <- domain_analysis %>%
  rowwise() %>%
  mutate(
    t_stat = (mean_R - overall_mean_R) / se_R,
    df = n - 1,
    p_value = 2 * pt(-abs(t_stat), df),
    direction = ifelse(mean_R > overall_mean_R, "MORE stable", "LESS stable")
  ) %>%
  ungroup()

# FDR correction (Benjamini-Hochberg)
domain_tests$p_adjusted <- p.adjust(domain_tests$p_value, method = "BH")
domain_tests$sig_raw <- domain_tests$p_value < 0.05
domain_tests$sig_fdr <- domain_tests$p_adjusted < 0.05

cat("DOMAIN ANALYSIS WITH FDR CORRECTION:\n")
cat(paste(rep("-", 60), collapse = ""), "\n")
cat(paste0("Overall mean R: ", round(overall_mean_R, 3), "\n"))
cat(paste0("Domains tested (n >= 10): ", nrow(domain_tests), "\n\n"))

# Results table
results_table <- domain_tests %>%
  arrange(p_adjusted) %>%
  select(outcome_domain, n, mean_R, pct_fragile, p_value, p_adjusted, sig_fdr)

cat("Results sorted by adjusted p-value:\n")
print(as.data.frame(results_table), row.names = FALSE)

cat(paste0("\nSignificant before FDR correction: ", sum(domain_tests$sig_raw), "\n"))
cat(paste0("Significant after FDR correction: ", sum(domain_tests$sig_fdr), "\n"))

# Report significant domains
sig_domains <- domain_tests %>% filter(sig_fdr) %>% arrange(p_adjusted)
if (nrow(sig_domains) > 0) {
  cat("\nDOMAINS WITH SIGNIFICANTLY DIFFERENT STABILITY:\n")
  for (i in 1:nrow(sig_domains)) {
    cat(paste0("  ", sig_domains$outcome_domain[i], ": ", sig_domains$direction[i],
               " (R = ", round(sig_domains$mean_R[i], 3),
               ", p_adj = ", format.pval(sig_domains$p_adjusted[i], digits = 2), ")\n"))
  }
}

saveRDS(domain_tests, "output/domain_analysis_fdr.rds")

################################################################################
# REVISION 7: HETEROGENEITY ANALYSIS
################################################################################

cat("\n\n=== REVISION 7: HETEROGENEITY ANALYSIS ===\n\n")

# Tau distribution
cat("Heterogeneity (tau) Distribution:\n")
cat(paste0("  Mean tau: ", round(mean(primary_data$tau, na.rm = TRUE), 4), "\n"))
cat(paste0("  Median tau: ", round(median(primary_data$tau, na.rm = TRUE), 4), "\n"))
cat(paste0("  SD: ", round(sd(primary_data$tau, na.rm = TRUE), 4), "\n"))

# Large heterogeneity
large_het <- sum(primary_data$tau > 0.5, na.rm = TRUE)
cat(paste0("\nLarge heterogeneity (tau > 0.5): ", large_het, " (",
           round(large_het / nrow(primary_data) * 100, 1), "%)\n"))

# Heterogeneity by domain
het_by_domain <- primary_data %>%
  group_by(outcome_domain) %>%
  summarise(
    n = n(),
    mean_tau = mean(tau, na.rm = TRUE),
    median_tau = median(tau, na.rm = TRUE),
    pct_large_tau = mean(tau > 0.5, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  filter(n >= 10) %>%
  arrange(desc(mean_tau))

cat("\nHeterogeneity by Domain (top 10):\n")
print(head(as.data.frame(het_by_domain), 10), row.names = FALSE)

################################################################################
# FINAL SUMMARY
################################################################################

cat("\n\n================================================================\n")
cat("REVISED ANALYSIS SUMMARY - READY FOR SUBMISSION\n")
cat("================================================================\n\n")

cat("DATA SCOPE:\n")
cat(paste0("  Total meta-analyses: ", nrow(ma4_data), "\n"))
cat(paste0("  Primary analysis (logRR): ", nrow(primary_data), "\n"))
cat(paste0("  Unique Cochrane reviews: ", length(unique(primary_data$review_id)), "\n"))
cat(paste0("  Outcome domains: ", length(unique(primary_data$outcome_domain)), "\n\n"))

cat("KEY RESULTS (REVISED):\n")
cat(paste(rep("-", 60), collapse = ""), "\n")

cat("\n1. STABILITY DISTRIBUTION:\n")
cat(paste0("   Highly Fragile (R < 0.5): ", round(r_cat_table[1] / sum(r_cat_table) * 100, 1), "%\n"))
cat(paste0("   Stable (R > 0.8): ", round(r_cat_table[3] / sum(r_cat_table) * 100, 1), "%\n"))

cat("\n2. PREDICTIVE MODEL (with proper CV):\n")
cat(paste0("   AUC-ROC: ", round(auc_value, 3), " (", round(auc_ci[1], 3), "-", round(auc_ci[3], 3), ")\n"))
cat(paste0("   Balanced Accuracy: ", round(cv_summary["balanced_accuracy"], 3), "\n"))
cat(paste0("   Sensitivity: ", round(cv_summary["sensitivity"], 3), "\n"))
cat(paste0("   Specificity: ", round(cv_summary["specificity"], 3), "\n"))

cat("\n3. MIXED-EFFECTS MODEL:\n")
cat(paste0("   ICC: ", round(icc, 3), " (", round(icc * 100, 1), "% between-review variance)\n"))
cat(paste0("   AIC improvement over OLS: ", round(aic_ols - aic_me, 1), "\n"))

cat("\n4. DOMAIN DIFFERENCES (FDR-corrected):\n")
cat(paste0("   ", sum(domain_tests$sig_fdr), " of ", nrow(domain_tests),
           " domains significantly different from average\n"))

cat("\n5. TOP PREDICTORS OF FRAGILITY:\n")
for (i in 1:min(3, nrow(importance_df))) {
  cat(paste0("   ", i, ". ", importance_df$Feature[i],
             " (Gini: ", round(importance_df$MeanDecreaseGini[i], 1), ")\n"))
}

cat("\n\nFILES GENERATED:\n")
cat("  output/predictive_model_revised.rds\n")
cat("  output/mixed_effects_model.rds\n")
cat("  output/domain_analysis_fdr.rds\n")

cat("\n\nEDITORIAL CONCERNS ADDRESSED:\n")
cat("  [x] Domain classification improved (22 outcome-based categories)\n")
cat("  [x] Effect types analyzed separately (logRR primary)\n")
cat("  [x] Predictive model: 10-fold CV, AUC-ROC, class weights for imbalance\n")
cat("  [x] Mixed-effects models account for review clustering\n")
cat("  [x] FDR correction for multiple domain comparisons\n")
cat("  [x] Heterogeneity properly characterized\n")

cat("\n================================================================\n")
cat("REVISION COMPLETE\n")
cat("================================================================\n")

# Save comprehensive results
saveRDS(list(
  data = primary_data,
  domain_analysis = domain_tests,
  cv_summary = cv_summary,
  auc = auc_value,
  icc = icc,
  importance = importance_df,
  het_by_domain = het_by_domain
), "output/EDITORIAL_REVISION_RESULTS.rds")
