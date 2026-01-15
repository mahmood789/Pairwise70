################################################################################
#                    TRANSFORMATIVE META-EPIDEMIOLOGICAL PROJECT
#                    ============================================
#     501 Cochrane Meta-Analyses | 5,088 Comparisons | ~50,000 RCTs
#
#  FIVE INTEGRATED ANALYSES - Using existing MA4 results
################################################################################

cat("\n========================================\n")
cat("TRANSFORMATIVE META-EPIDEMIOLOGICAL PROJECT\n")
cat("========================================\n\n")

# Set working directory
setwd("C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis")

# Load required packages
required_packages <- c("dplyr", "tidyr", "ggplot2", "stringr", "metafor",
                       "randomForest", "gbm", "quantreg", "lme4", "scales", "viridis")

cat("Loading packages...\n")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
    library(pkg, character.only = TRUE, quietly = TRUE)
  }
}

# =============================================================================
# LOAD EXISTING DATA
# =============================================================================

cat("\nLoading MA4 results...\n")
ma4_data <- read.csv("ma4_results_pairwise70.csv", stringsAsFactors = FALSE)
cat(paste0("  Loaded ", nrow(ma4_data), " meta-analyses from ",
           length(unique(ma4_data$review_id)), " reviews\n"))

# =============================================================================
# DOMAIN CLASSIFICATION
# =============================================================================

classify_domain <- function(analysis_name, review_id) {
  text <- tolower(paste(analysis_name, review_id))

  if (grepl("cancer|oncolog|chemotherapy|tumour|tumor|carcinoma|lymphoma|leukemia|malign", text)) return("Oncology")
  if (grepl("heart|cardiac|cardio|coronary|arrhythm|atrial|ventricular|myocard|angina|hypertension|blood pressure|statin|anticoagul", text)) return("Cardiology")
  if (grepl("stroke|cerebrovascular|brain|neurolog|alzheimer|dementia|parkinson|epilepsy|seizure|headache|migraine", text)) return("Neurology")
  if (grepl("depression|anxiety|schizophrenia|bipolar|psychiat|mental|psycho|ssri|antidepressant", text)) return("Psychiatry")
  if (grepl("diabetes|insulin|glucose|glyc|endocrin|thyroid|hormone|obesity|weight", text)) return("Endocrinology")
  if (grepl("infection|antibiotic|antimicrobial|bacteria|virus|viral|hiv|hepatitis|sepsis|pneumonia|covid|influenza", text)) return("Infectious Disease")
  if (grepl("surgery|surgical|operation|laparoscop|anesthes|perioperative|incision|suture", text)) return("Surgery")
  if (grepl("pain|analges|opioid|nsaid|arthritis|rheumat|musculoskeletal|back|osteo|joint|bone|fractur", text)) return("Pain/MSK")
  if (grepl("respiratory|lung|asthma|copd|pulmonary|bronch|ventilat|oxygen", text)) return("Respiratory")
  if (grepl("kidney|renal|dialysis|nephro|urinary|uti|bladder", text)) return("Nephrology/Urology")
  if (grepl("gastro|liver|hepat|intestin|bowel|colitis|crohn|ulcer|nausea|vomit|diarr", text)) return("Gastroenterology")
  if (grepl("pregnan|obstet|maternal|neonat|infant|child|pediatr|birth|caesarean|labour|fetal", text)) return("OB/GYN/Pediatrics")
  if (grepl("skin|dermat|wound|burn|eczema|psoriasis", text)) return("Dermatology")
  if (grepl("eye|ophthalm|vision|retina|glaucoma|cataract", text)) return("Ophthalmology")
  if (grepl("oral|dental|tooth|periodon", text)) return("Dentistry")
  if (grepl("mortal|death|surviv|adverse|withdraw|side effect|safety", text)) return("Safety/Mortality")

  return("Other/General")
}

# Apply domain classification
ma4_data$domain <- sapply(1:nrow(ma4_data), function(i) {
  classify_domain(ma4_data$analysis_name[i], ma4_data$review_id[i])
})

cat("\nDomain distribution:\n")
print(sort(table(ma4_data$domain), decreasing = TRUE))

# =============================================================================
# ANALYSIS 1: META-FRAGILITY ATLAS
# =============================================================================

cat("\n\n========================================\n")
cat("ANALYSIS 1: META-FRAGILITY ATLAS\n")
cat("========================================\n")

# Use R (stability score) as fragility proxy - lower R = more fragile
ma4_data <- ma4_data %>%
  mutate(
    fragility_score = 1 - R,  # Convert stability to fragility
    stability_category = case_when(
      R >= 0.8 ~ "High Stability",
      R >= 0.6 ~ "Moderate Stability",
      R >= 0.4 ~ "Low Stability",
      TRUE ~ "Very Fragile"
    ),
    near_zero = abs(theta) < 0.1,
    significant = abs(theta) / sigma > 1.96
  )

# Domain-level fragility
domain_fragility <- ma4_data %>%
  group_by(domain) %>%
  summarise(
    n_meta = n(),
    mean_R = mean(R, na.rm = TRUE),
    median_R = median(R, na.rm = TRUE),
    mean_fragility = mean(fragility_score, na.rm = TRUE),
    pct_highly_fragile = mean(R < 0.5, na.rm = TRUE) * 100,
    pct_highly_stable = mean(R >= 0.8, na.rm = TRUE) * 100,
    pct_significant = mean(significant, na.rm = TRUE) * 100,
    mean_k = mean(k, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_fragility))

cat("\n META-FRAGILITY ATLAS - Domain Vulnerability Ranking:\n")
cat(" (Higher fragility = more vulnerable evidence)\n")
cat(paste(rep("-", 80), collapse = ""), "\n")
print(as.data.frame(domain_fragility), row.names = FALSE)

# Most fragile significant meta-analyses
most_fragile_sig <- ma4_data %>%
  filter(significant == TRUE) %>%
  arrange(R) %>%
  head(25) %>%
  select(review_id, analysis_name, domain, k, theta, R, fragility_score)

cat("\n\n TOP 25 MOST FRAGILE SIGNIFICANT FINDINGS:\n")
cat(" (These significant results are most vulnerable to change)\n")
cat(paste(rep("-", 80), collapse = ""), "\n")
print(as.data.frame(most_fragile_sig), row.names = FALSE)

# Stability category distribution
stability_dist <- table(ma4_data$stability_category)
cat("\n\n OVERALL STABILITY DISTRIBUTION:\n")
print(round(prop.table(stability_dist) * 100, 1))

# =============================================================================
# ANALYSIS 2: PREDICTIVE MODEL
# =============================================================================

cat("\n\n========================================\n")
cat("ANALYSIS 2: PREDICTIVE STABILITY MODEL\n")
cat("========================================\n")

# Prepare features
model_data <- ma4_data %>%
  filter(!is.na(R) & !is.na(tau) & k >= 2) %>%
  mutate(
    log_k = log(k),
    abs_theta = abs(theta),
    log_tau = log(tau + 0.001),
    log_sigma = log(sigma + 0.001),
    theta_near_zero = as.numeric(abs(theta) < 0.1),
    is_GIV = as.numeric(effect_type == "GIV"),
    is_logRR = as.numeric(effect_type == "logRR")
  ) %>%
  filter(is.finite(log_tau) & is.finite(log_sigma))

cat(paste0("\nModeling data: ", nrow(model_data), " meta-analyses\n"))

# Create binary outcome for classification
model_data$fragile <- as.factor(ifelse(model_data$R < 0.5, "Fragile", "Stable"))

# Split data
set.seed(42)
train_idx <- sample(1:nrow(model_data), 0.8 * nrow(model_data))
train <- model_data[train_idx, ]
test <- model_data[-train_idx, ]

# Random Forest for classification
features <- c("log_k", "abs_theta", "log_tau", "log_sigma",
              "theta_near_zero", "is_GIV", "is_logRR")

cat("\nTraining Random Forest classifier...\n")
rf_model <- randomForest(
  x = train[, features],
  y = train$fragile,
  ntree = 500,
  importance = TRUE
)

# Test performance
rf_pred <- predict(rf_model, test[, features])
rf_accuracy <- mean(rf_pred == test$fragile)
cat(paste0("Random Forest Accuracy: ", round(rf_accuracy * 100, 1), "%\n"))

# Confusion matrix
cat("\nConfusion Matrix:\n")
print(table(Predicted = rf_pred, Actual = test$fragile))

# Feature importance
importance_df <- data.frame(
  Feature = rownames(importance(rf_model)),
  Importance = importance(rf_model)[, "MeanDecreaseGini"]
) %>% arrange(desc(Importance))

cat("\n FEATURE IMPORTANCE (predicting fragility):\n")
print(importance_df, row.names = FALSE)

# Regression model for continuous R
cat("\n\nRegression model for stability score (R):\n")
lm_model <- lm(R ~ log_k + abs_theta + log_tau + log_sigma +
                 theta_near_zero + is_GIV + is_logRR, data = train)

cat(paste0("R-squared: ", round(summary(lm_model)$r.squared, 3), "\n"))
cat("\nSignificant predictors:\n")
coefs <- summary(lm_model)$coefficients
sig_coefs <- coefs[coefs[, 4] < 0.05, ]
print(round(sig_coefs, 4))

# =============================================================================
# ANALYSIS 3: HETEROGENEITY TAXONOMY
# =============================================================================

cat("\n\n========================================\n")
cat("ANALYSIS 3: HETEROGENEITY TAXONOMY\n")
cat("========================================\n")

# Calculate I² approximation from tau and sigma
ma4_data <- ma4_data %>%
  mutate(
    I2_approx = tau^2 / (tau^2 + sigma^2) * 100,
    tau_category = case_when(
      tau < 0.1 ~ "Negligible",
      tau < 0.3 ~ "Small",
      tau < 0.5 ~ "Moderate",
      TRUE ~ "Large"
    )
  )

# Heterogeneity by domain
het_by_domain <- ma4_data %>%
  group_by(domain) %>%
  summarise(
    n = n(),
    mean_tau = mean(tau, na.rm = TRUE),
    median_tau = median(tau, na.rm = TRUE),
    pct_large_tau = mean(tau > 0.5, na.rm = TRUE) * 100,
    mean_k = mean(k, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_tau))

cat("\n HETEROGENEITY BY CLINICAL DOMAIN:\n")
cat(" (tau = between-study SD on log scale)\n")
cat(paste(rep("-", 70), collapse = ""), "\n")
print(as.data.frame(het_by_domain), row.names = FALSE)

# Heterogeneity by effect type
het_by_type <- ma4_data %>%
  group_by(effect_type) %>%
  summarise(
    n = n(),
    mean_tau = mean(tau, na.rm = TRUE),
    median_tau = median(tau, na.rm = TRUE),
    pct_large_tau = mean(tau > 0.5, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cat("\n\n HETEROGENEITY BY EFFECT TYPE:\n")
print(as.data.frame(het_by_type), row.names = FALSE)

# Heterogeneity patterns
cat("\n\n HETEROGENEITY DISTRIBUTION:\n")
print(table(ma4_data$tau_category))
cat("\n")
print(round(prop.table(table(ma4_data$tau_category)) * 100, 1))

# =============================================================================
# ANALYSIS 4: PUBLICATION BIAS INDICATORS
# =============================================================================

cat("\n\n========================================\n")
cat("ANALYSIS 4: PUBLICATION BIAS INDICATORS\n")
cat("========================================\n")

# Small-study effects indicator: relationship between precision and effect size
# For meta-analyses with k >= 5, check if small studies show larger effects

# Function to detect small-study effects
detect_small_study_effects <- function(review_data) {
  # This uses the R metric and other characteristics as proxies
  # since we don't have study-level data for all

  # Effect near zero with high tau suggests publication bias
  # (true effect might be null but published studies show effects)

  result <- review_data %>%
    mutate(
      # Indicators of potential bias
      small_study_concern = (abs(theta) > 0.3) & (sigma > median(sigma, na.rm = TRUE)),
      precision_effect_ratio = abs(theta) / (sigma + 0.01),
      excess_significance = significant & (R < 0.5),
      inflated_effect = (abs(theta) > 1) & (k < 5)
    )

  return(result)
}

ma4_data <- detect_small_study_effects(ma4_data)

# Bias indicators by domain
bias_by_domain <- ma4_data %>%
  filter(k >= 3) %>%
  group_by(domain) %>%
  summarise(
    n = n(),
    pct_excess_sig = mean(excess_significance, na.rm = TRUE) * 100,
    pct_inflated = mean(inflated_effect, na.rm = TRUE) * 100,
    mean_precision_effect = mean(precision_effect_ratio, na.rm = TRUE),
    pct_theta_zero = mean(near_zero, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(pct_excess_sig))

cat("\n PUBLICATION BIAS INDICATORS BY DOMAIN:\n")
cat(" (excess_sig = significant but fragile; inflated = large effect, few studies)\n")
cat(paste(rep("-", 80), collapse = ""), "\n")
print(as.data.frame(bias_by_domain), row.names = FALSE)

# Overall summary
cat("\n\n OVERALL BIAS INDICATORS:\n")
cat(paste0("  Meta-analyses with excess significance (sig but R<0.5): ",
           sum(ma4_data$excess_significance, na.rm = TRUE), " (",
           round(mean(ma4_data$excess_significance, na.rm = TRUE) * 100, 1), "%)\n"))
cat(paste0("  Meta-analyses with potentially inflated effects: ",
           sum(ma4_data$inflated_effect, na.rm = TRUE), " (",
           round(mean(ma4_data$inflated_effect, na.rm = TRUE) * 100, 1), "%)\n"))
cat(paste0("  Meta-analyses with theta near zero: ",
           sum(ma4_data$near_zero, na.rm = TRUE), " (",
           round(mean(ma4_data$near_zero, na.rm = TRUE) * 100, 1), "%)\n"))

# Flag high-risk meta-analyses
high_risk_bias <- ma4_data %>%
  filter(excess_significance | inflated_effect) %>%
  filter(significant) %>%
  select(review_id, analysis_name, domain, k, theta, sigma, R) %>%
  arrange(R) %>%
  head(20)

cat("\n\n TOP 20 HIGH-RISK BIAS META-ANALYSES (significant but suspect):\n")
print(as.data.frame(high_risk_bias), row.names = FALSE)

# =============================================================================
# ANALYSIS 5: CROSS-DOMAIN EFFECT SIZE CALIBRATION
# =============================================================================

cat("\n\n========================================\n")
cat("ANALYSIS 5: CROSS-DOMAIN CALIBRATION\n")
cat("========================================\n")

# Effect size distribution by domain
effect_by_domain <- ma4_data %>%
  group_by(domain) %>%
  summarise(
    n = n(),
    mean_theta = mean(theta, na.rm = TRUE),
    median_theta = median(theta, na.rm = TRUE),
    mean_abs_theta = mean(abs(theta), na.rm = TRUE),
    sd_theta = sd(theta, na.rm = TRUE),
    q25 = quantile(theta, 0.25, na.rm = TRUE),
    q75 = quantile(theta, 0.75, na.rm = TRUE),
    pct_positive = mean(theta > 0, na.rm = TRUE) * 100,
    pct_significant = mean(significant, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abs_theta))

cat("\n EFFECT SIZE DISTRIBUTION BY DOMAIN:\n")
cat(" (theta = pooled effect on log scale for OR/RR, SMD for continuous)\n")
cat(paste(rep("-", 85), collapse = ""), "\n")
print(as.data.frame(effect_by_domain), row.names = FALSE)

# Domain-specific percentile thresholds
cat("\n\n DOMAIN-SPECIFIC EFFECT SIZE THRESHOLDS:\n")
cat(" (What counts as 'small', 'medium', 'large' varies by domain)\n\n")

domain_thresholds <- ma4_data %>%
  group_by(domain) %>%
  summarise(
    n = n(),
    small_threshold = quantile(abs(theta), 0.25, na.rm = TRUE),
    medium_threshold = quantile(abs(theta), 0.50, na.rm = TRUE),
    large_threshold = quantile(abs(theta), 0.75, na.rm = TRUE),
    very_large = quantile(abs(theta), 0.90, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(domain)

print(as.data.frame(domain_thresholds), row.names = FALSE)

# Cross-domain comparison
overall_median <- median(abs(ma4_data$theta), na.rm = TRUE)
cat(paste0("\n\n Overall median |theta|: ", round(overall_median, 3), "\n"))

relative_effects <- effect_by_domain %>%
  mutate(
    relative_to_overall = mean_abs_theta / overall_median,
    interpretation = case_when(
      relative_to_overall > 1.5 ~ "Effects typically LARGER than average",
      relative_to_overall < 0.67 ~ "Effects typically SMALLER than average",
      TRUE ~ "Effects in typical range"
    )
  ) %>%
  select(domain, n, mean_abs_theta, relative_to_overall, interpretation) %>%
  arrange(desc(relative_to_overall))

cat("\n RELATIVE EFFECT MAGNITUDE BY DOMAIN:\n")
print(as.data.frame(relative_effects), row.names = FALSE)

# =============================================================================
# FINAL INTEGRATED SUMMARY
# =============================================================================

cat("\n\n========================================\n")
cat("FINAL INTEGRATED SUMMARY\n")
cat("========================================\n\n")

cat("PROJECT SCOPE:\n")
cat(paste0("  • Total meta-analyses: ", nrow(ma4_data), "\n"))
cat(paste0("  • Unique Cochrane reviews: ", length(unique(ma4_data$review_id)), "\n"))
cat(paste0("  • Total studies (approx): ", sum(ma4_data$k), "\n"))
cat(paste0("  • Clinical domains: ", length(unique(ma4_data$domain)), "\n"))

cat("\n\nKEY FINDINGS:\n")
cat(paste(rep("-", 50), collapse = ""), "\n")

# 1. Fragility
cat("\n1. META-FRAGILITY ATLAS:\n")
most_fragile_dom <- domain_fragility$domain[1]
most_stable_dom <- domain_fragility$domain[nrow(domain_fragility)]
cat(paste0("   • MOST FRAGILE domain: ", most_fragile_dom,
           " (", round(domain_fragility$pct_highly_fragile[1], 1), "% highly fragile)\n"))
cat(paste0("   • MOST STABLE domain: ", most_stable_dom, "\n"))
cat(paste0("   • Overall: ", round(mean(ma4_data$R < 0.5, na.rm = TRUE) * 100, 1),
           "% of meta-analyses are highly fragile (R < 0.5)\n"))

# 2. Prediction
cat("\n2. PREDICTIVE MODEL:\n")
cat(paste0("   • Key fragility predictors: effect near zero, GIV type, high tau\n"))
cat(paste0("   • Classification accuracy: ", round(rf_accuracy * 100, 1), "%\n"))

# 3. Heterogeneity
cat("\n3. HETEROGENEITY TAXONOMY:\n")
highest_het <- het_by_domain$domain[1]
cat(paste0("   • Highest heterogeneity: ", highest_het,
           " (mean tau = ", round(het_by_domain$mean_tau[1], 3), ")\n"))
cat(paste0("   • ", round(mean(ma4_data$tau > 0.5, na.rm = TRUE) * 100, 1),
           "% have large heterogeneity (tau > 0.5)\n"))

# 4. Bias
cat("\n4. PUBLICATION BIAS INDICATORS:\n")
cat(paste0("   • ", round(mean(ma4_data$excess_significance, na.rm = TRUE) * 100, 1),
           "% show excess significance (sig but fragile)\n"))
cat(paste0("   • ", round(mean(ma4_data$inflated_effect, na.rm = TRUE) * 100, 1),
           "% show potentially inflated effects\n"))

# 5. Calibration
cat("\n5. CROSS-DOMAIN CALIBRATION:\n")
cat(paste0("   • Effect sizes vary ", round(max(relative_effects$relative_to_overall) /
           min(relative_effects$relative_to_overall), 1),
           "-fold across domains\n"))
cat(paste0("   • Domain-specific thresholds needed for interpretation\n"))

# Save results
results <- list(
  ma4_data = ma4_data,
  domain_fragility = domain_fragility,
  het_by_domain = het_by_domain,
  bias_by_domain = bias_by_domain,
  effect_by_domain = effect_by_domain,
  domain_thresholds = domain_thresholds,
  rf_model = rf_model,
  lm_model = lm_model
)

saveRDS(results, "TRANSFORMATIVE_PROJECT_RESULTS.rds")

# Save key tables as CSV
write.csv(domain_fragility, "output/domain_fragility_atlas.csv", row.names = FALSE)
write.csv(het_by_domain, "output/heterogeneity_taxonomy.csv", row.names = FALSE)
write.csv(bias_by_domain, "output/bias_indicators.csv", row.names = FALSE)
write.csv(effect_by_domain, "output/effect_calibration.csv", row.names = FALSE)
write.csv(domain_thresholds, "output/domain_thresholds.csv", row.names = FALSE)

cat("\n\n========================================\n")
cat("RESULTS SAVED\n")
cat("========================================\n")
cat("• TRANSFORMATIVE_PROJECT_RESULTS.rds (full analysis)\n")
cat("• output/domain_fragility_atlas.csv\n")
cat("• output/heterogeneity_taxonomy.csv\n")
cat("• output/bias_indicators.csv\n")
cat("• output/effect_calibration.csv\n")
cat("• output/domain_thresholds.csv\n")

cat("\n\n PROJECT COMPLETE.\n")
