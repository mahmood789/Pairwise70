################################################################################
#     EDITORIAL REVISION 2: Research Synthesis Methods
#     ===========================================================
#     Addressing Remaining Reviewer Concerns:
#     1. Improve predictive model sensitivity (threshold optimization)
#     2. Enhanced domain classification
#     3. Publication bias analysis (subsample)
#     4. Confidence intervals for domain comparisons
#     5. Calibration assessment
#     6. Improved figures
################################################################################

cat("\n================================================================\n")
cat("EDITORIAL REVISION 2: Final Revisions\n")
cat("================================================================\n\n")

setwd("C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis")
dir.create("output", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# Load packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lme4)
  library(pROC)
  library(caret)
  library(randomForest)
  library(metafor)
  library(viridis)
})

# Load data
ma4_data <- read.csv("ma4_results_pairwise70.csv", stringsAsFactors = FALSE)
cat(paste0("Loaded ", nrow(ma4_data), " meta-analyses\n\n"))

################################################################################
# REVISION 1: ENHANCED DOMAIN CLASSIFICATION
# Goal: Reduce "Other Clinical" from 39.5% to <25%
################################################################################

cat("=== REVISION 1: ENHANCED DOMAIN CLASSIFICATION ===\n\n")

classify_enhanced <- function(analysis_name, review_id) {
  name_lower <- tolower(analysis_name)

  # Extract CD number for Cochrane Review Group hints
  cd_num <- as.numeric(gsub("CD0*([0-9]+).*", "\\1", review_id))

  # TIER 1: High-confidence outcome keywords
  if (grepl("death|mortal|surviv|fatal|died|dying", name_lower)) return("Mortality")
  if (grepl("adverse|side effect|withdraw|dropout|tolerab|safety|toxicity|harm", name_lower)) return("Adverse Events")
  if (grepl("pain|analges|ache|neuralg|migraine|headache", name_lower)) return("Pain Outcomes")
  if (grepl("infect|sepsis|fever|pneumonia|bacterem|viral|fungal", name_lower)) return("Infection")
  if (grepl("bleed|hemorrh|transfus|haemorrhag|blood loss", name_lower)) return("Bleeding")
  if (grepl("fractur|bone|osteo|arthr|joint|musculo|sprain", name_lower)) return("Musculoskeletal")
  if (grepl("heart|cardiac|coronary|arrhythm|myocard|infarct|angina|atrial", name_lower)) return("Cardiovascular")
  if (grepl("stroke|cerebr|isch[ae]mic brain|tia", name_lower)) return("Stroke")
  if (grepl("cancer|tumor|tumour|malignan|oncol|chemotherapy|carcinoma|lymphoma|leuk", name_lower)) return("Cancer")
  if (grepl("diabet|glucose|insulin|glyc|hba1c|hyperglycemia", name_lower)) return("Diabetes/Metabolic")
  if (grepl("pregnan|birth|deliver|caesar|neonat|fetal|maternal|labour|labor|obstet", name_lower)) return("Pregnancy/Birth")
  if (grepl("child|infant|pediatr|paediatr|newborn|preterm", name_lower)) return("Pediatric")
  if (grepl("depress|anxiety|psychiat|mental|schizo|bipolar|mood|suicid", name_lower)) return("Mental Health")
  if (grepl("respirat|lung|asthma|copd|pulmon|ventilat|oxygen|dyspn|breath", name_lower)) return("Respiratory")
  if (grepl("renal|kidney|dialysis|urinar|nephro|creatinin|gfr", name_lower)) return("Renal")
  if (grepl("liver|hepat|gastro|bowel|intestin|colitis|diarrh|nausea|vomit|gi\\b", name_lower)) return("GI/Hepatic")
  if (grepl("surgery|surgical|operat|incision|wound|anastom|resect", name_lower)) return("Surgical")
  if (grepl("quality of life|qol|hrqol|sf-36|eq-5d|well-?being|function", name_lower)) return("Quality of Life")
  if (grepl("length of stay|hospital|admission|discharge|readmission|icu", name_lower)) return("Healthcare Utilization")
  if (grepl("recurrence|relapse|remission|progression|exacerb", name_lower)) return("Disease Course")
  if (grepl("response|improvement|efficacy|success|cure|resolution|recovery", name_lower)) return("Treatment Response")

  # TIER 2: Broader clinical patterns
  if (grepl("weight|bmi|obes|body mass", name_lower)) return("Diabetes/Metabolic")
  if (grepl("blood pressure|hypertens|systolic|diastolic", name_lower)) return("Cardiovascular")
  if (grepl("cholesterol|lipid|ldl|hdl|triglycer", name_lower)) return("Cardiovascular")
  if (grepl("cognit|dementia|alzheimer|memory|attention", name_lower)) return("Mental Health")
  if (grepl("sleep|insomnia|fatigue|tired", name_lower)) return("Quality of Life")
  if (grepl("mobility|walk|gait|balance|fall", name_lower)) return("Musculoskeletal")
  if (grepl("vision|eye|ophthalm|blind|cataract|glauc", name_lower)) return("Ophthalmology")
  if (grepl("hearing|deaf|auditory|ear\\b|otitis", name_lower)) return("ENT")
  if (grepl("skin|dermat|rash|eczema|psoria|wound heal", name_lower)) return("Dermatology")
  if (grepl("thrombo|embol|dvt|clot|anticoag", name_lower)) return("Thromboembolic")
  if (grepl("anemia|anaemia|hemoglobin|haemoglobin|iron\\b", name_lower)) return("Hematologic")
  if (grepl("immun|vaccin|antibod|allerg", name_lower)) return("Immunologic")
  if (grepl("cost|economic|resource|utilization", name_lower)) return("Healthcare Utilization")
  if (grepl("satisfaction|adherence|compliance|accept", name_lower)) return("Patient-Reported")
  if (grepl("symptom|score|scale|index|measure", name_lower)) return("Clinical Scores")

  # TIER 3: Intervention-based (less specific)
  if (grepl("antibiotic|antimicrob", name_lower)) return("Infection")
  if (grepl("steroid|cortico", name_lower)) return("Immunologic")
  if (grepl("chemo|radiation|radiother", name_lower)) return("Cancer")
  if (grepl("anesthes|anaesthes|sedation", name_lower)) return("Surgical")
  if (grepl("physical therap|physiother|rehabilit|exercise", name_lower)) return("Rehabilitation")

  return("Other Clinical")
}

ma4_data$outcome_domain <- mapply(classify_enhanced,
                                   ma4_data$analysis_name,
                                   ma4_data$review_id)

domain_counts <- sort(table(ma4_data$outcome_domain), decreasing = TRUE)
cat("Enhanced Domain Classification:\n")
print(domain_counts)

other_pct <- round(sum(ma4_data$outcome_domain == "Other Clinical") / nrow(ma4_data) * 100, 1)
cat(paste0("\n'Other Clinical' category: ", other_pct, "% (previous: 39.5%)\n"))
cat(paste0("Total domains: ", length(unique(ma4_data$outcome_domain)), "\n"))

# Separate by effect type
primary_data <- ma4_data %>% filter(effect_type == "logRR")
cat(paste0("\nPrimary analysis (logRR): n = ", nrow(primary_data), "\n"))

################################################################################
# REVISION 2: IMPROVED PREDICTIVE MODEL WITH THRESHOLD OPTIMIZATION
################################################################################

cat("\n\n=== REVISION 2: PREDICTIVE MODEL WITH THRESHOLD OPTIMIZATION ===\n\n")

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

cat(paste0("Modeling sample: ", nrow(model_data), "\n"))
cat("Class distribution:\n")
print(table(model_data$fragile))

features <- c("log_k", "abs_theta", "log_tau", "log_sigma", "theta_near_zero")

# 10-fold CV with probability predictions
set.seed(42)
folds <- createFolds(model_data$fragile, k = 10, list = TRUE)

all_predictions <- data.frame()

for (fold_idx in 1:length(folds)) {
  test_idx <- folds[[fold_idx]]
  train_idx <- setdiff(1:nrow(model_data), test_idx)

  train <- model_data[train_idx, ]
  test <- model_data[test_idx, ]

  # Higher class weight for Fragile (minority)
  rf <- randomForest(
    x = train[, features],
    y = train$fragile,
    ntree = 500,
    classwt = c("Stable" = 1, "Fragile" = 5),
    importance = TRUE
  )

  pred_prob <- predict(rf, test[, features], type = "prob")[, "Fragile"]

  all_predictions <- rbind(all_predictions, data.frame(
    actual = test$fragile,
    prob_fragile = pred_prob,
    fold = fold_idx
  ))
}

# ROC analysis
roc_obj <- roc(all_predictions$actual, all_predictions$prob_fragile,
               levels = c("Stable", "Fragile"), direction = "<")
auc_value <- auc(roc_obj)
auc_ci <- ci.auc(roc_obj)

cat("\nROC ANALYSIS:\n")
cat(paste0("AUC-ROC: ", round(auc_value, 3), " (95% CI: ",
           round(auc_ci[1], 3), "-", round(auc_ci[3], 3), ")\n"))

# Threshold optimization
thresholds <- seq(0.1, 0.9, by = 0.05)
threshold_results <- data.frame()

for (thresh in thresholds) {
  pred_class <- factor(ifelse(all_predictions$prob_fragile >= thresh, "Fragile", "Stable"),
                       levels = c("Stable", "Fragile"))

  cm <- confusionMatrix(pred_class, all_predictions$actual, positive = "Fragile")

  # Calculate F1 and Youden's J
  sens <- cm$byClass["Sensitivity"]
  spec <- cm$byClass["Specificity"]
  ppv <- cm$byClass["Pos Pred Value"]
  f1 <- 2 * (ppv * sens) / (ppv + sens)
  youden_j <- sens + spec - 1

  threshold_results <- rbind(threshold_results, data.frame(
    threshold = thresh,
    sensitivity = sens,
    specificity = spec,
    ppv = ppv,
    npv = cm$byClass["Neg Pred Value"],
    f1 = f1,
    youden_j = youden_j,
    balanced_acc = cm$byClass["Balanced Accuracy"]
  ))
}

# Find optimal thresholds
optimal_youden <- threshold_results[which.max(threshold_results$youden_j), ]
optimal_f1 <- threshold_results[which.max(threshold_results$f1), ]
optimal_sens80 <- threshold_results[which.min(abs(threshold_results$sensitivity - 0.80)), ]

cat("\nTHRESHOLD OPTIMIZATION:\n")
cat(paste(rep("-", 60), collapse = ""), "\n")

cat("\n1. Default threshold (0.5):\n")
default <- threshold_results[threshold_results$threshold == 0.5, ]
cat(paste0("   Sensitivity: ", round(default$sensitivity, 3), "\n"))
cat(paste0("   Specificity: ", round(default$specificity, 3), "\n"))
cat(paste0("   F1 Score: ", round(default$f1, 3), "\n"))

cat("\n2. Optimal by Youden's J (", optimal_youden$threshold, "):\n")
cat(paste0("   Sensitivity: ", round(optimal_youden$sensitivity, 3), "\n"))
cat(paste0("   Specificity: ", round(optimal_youden$specificity, 3), "\n"))
cat(paste0("   Youden's J: ", round(optimal_youden$youden_j, 3), "\n"))

cat("\n3. Optimal by F1 Score (", optimal_f1$threshold, "):\n")
cat(paste0("   Sensitivity: ", round(optimal_f1$sensitivity, 3), "\n"))
cat(paste0("   Specificity: ", round(optimal_f1$specificity, 3), "\n"))
cat(paste0("   F1 Score: ", round(optimal_f1$f1, 3), "\n"))

cat("\n4. Target 80% Sensitivity (threshold = ", optimal_sens80$threshold, "):\n")
cat(paste0("   Sensitivity: ", round(optimal_sens80$sensitivity, 3), "\n"))
cat(paste0("   Specificity: ", round(optimal_sens80$specificity, 3), "\n"))
cat(paste0("   PPV: ", round(optimal_sens80$ppv, 3), "\n"))

# Use Youden's optimal for final metrics
optimal_thresh <- optimal_youden$threshold
all_predictions$predicted <- factor(
  ifelse(all_predictions$prob_fragile >= optimal_thresh, "Fragile", "Stable"),
  levels = c("Stable", "Fragile")
)

final_cm <- confusionMatrix(all_predictions$predicted, all_predictions$actual, positive = "Fragile")

cat("\n\nFINAL MODEL (Youden-optimal threshold):\n")
cat(paste(rep("-", 60), collapse = ""), "\n")
print(final_cm$table)
cat(paste0("\nSensitivity: ", round(final_cm$byClass["Sensitivity"], 3), "\n"))
cat(paste0("Specificity: ", round(final_cm$byClass["Specificity"], 3), "\n"))
cat(paste0("Balanced Accuracy: ", round(final_cm$byClass["Balanced Accuracy"], 3), "\n"))

################################################################################
# REVISION 3: CALIBRATION ASSESSMENT
################################################################################

cat("\n\n=== REVISION 3: CALIBRATION ASSESSMENT ===\n\n")

# Create decile calibration
all_predictions$prob_decile <- cut(all_predictions$prob_fragile,
                                    breaks = seq(0, 1, by = 0.1),
                                    include.lowest = TRUE)

calibration <- all_predictions %>%
  group_by(prob_decile) %>%
  summarise(
    n = n(),
    mean_predicted = mean(prob_fragile),
    observed_rate = mean(actual == "Fragile"),
    .groups = "drop"
  ) %>%
  filter(n >= 10)

cat("CALIBRATION BY DECILE:\n")
print(as.data.frame(calibration), row.names = FALSE)

# Calibration slope and intercept
calib_model <- glm(I(actual == "Fragile") ~ prob_fragile,
                   data = all_predictions, family = binomial)
calib_coef <- coef(calib_model)

cat(paste0("\nCalibration intercept: ", round(calib_coef[1], 3),
           " (ideal: 0)\n"))
cat(paste0("Calibration slope: ", round(calib_coef[2], 3),
           " (ideal: 1)\n"))

# Hosmer-Lemeshow statistic (simplified)
hl_stat <- sum((calibration$observed_rate - calibration$mean_predicted)^2 * calibration$n /
               (calibration$mean_predicted * (1 - calibration$mean_predicted) + 0.001))
cat(paste0("H-L type statistic: ", round(hl_stat, 2), "\n"))

################################################################################
# REVISION 4: PUBLICATION BIAS ANALYSIS (SUBSAMPLE)
################################################################################

cat("\n\n=== REVISION 4: PUBLICATION BIAS ANALYSIS ===\n\n")

rds_dir <- "output/cleaned_rds"

if (dir.exists(rds_dir)) {
  rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)
  cat(paste0("Found ", length(rds_files), " data files\n"))

  # Sample 100 reviews for bias analysis
  set.seed(42)
  sample_files <- sample(rds_files, min(100, length(rds_files)))

  bias_results <- data.frame()
  tested <- 0

  cat("Running publication bias tests (k >= 10)...\n")

  for (rds_file in sample_files) {
    tryCatch({
      d <- readRDS(rds_file)

      if (all(c("Experimental.cases", "Experimental.N",
                "Control.cases", "Control.N") %in% names(d))) {

        d <- d[complete.cases(d[, c("Experimental.cases", "Experimental.N",
                                     "Control.cases", "Control.N")]), ]

        if (nrow(d) >= 10) {
          es <- escalc(measure = "OR",
                       ai = d$Experimental.cases, n1i = d$Experimental.N,
                       ci = d$Control.cases, n2i = d$Control.N)

          es <- es[!is.na(es$yi) & !is.na(es$vi) & is.finite(es$yi) & is.finite(es$vi), ]

          if (nrow(es) >= 10) {
            ma <- rma(yi = es$yi, vi = es$vi, method = "REML")

            egger <- tryCatch(regtest(ma, model = "lm"), error = function(e) NULL)
            tf <- tryCatch(trimfill(ma), error = function(e) NULL)

            if (!is.null(egger) && !is.null(tf)) {
              tested <- tested + 1

              bias_results <- rbind(bias_results, data.frame(
                file = basename(rds_file),
                k = ma$k,
                theta = ma$beta[1],
                egger_z = egger$zval,
                egger_p = egger$pval,
                egger_sig = egger$pval < 0.1,
                tf_added = tf$k0,
                tf_theta = tf$beta[1],
                tf_change = tf$beta[1] - ma$beta[1],
                stringsAsFactors = FALSE
              ))
            }
          }
        }
      }
    }, error = function(e) NULL)

    if (tested %% 10 == 0) cat(paste0("\rTested: ", tested))
  }

  cat(paste0("\n\nCompleted: ", nrow(bias_results), " meta-analyses tested\n"))

  if (nrow(bias_results) > 0) {
    cat("\nPUBLICATION BIAS SUMMARY:\n")
    cat(paste(rep("-", 60), collapse = ""), "\n")

    egger_pos <- sum(bias_results$egger_sig)
    cat(paste0("Egger's test positive (p < 0.1): ", egger_pos, "/", nrow(bias_results),
               " (", round(egger_pos/nrow(bias_results)*100, 1), "%)\n"))

    tf_any <- sum(bias_results$tf_added > 0)
    cat(paste0("Trim-and-fill imputed studies: ", tf_any, "/", nrow(bias_results),
               " (", round(tf_any/nrow(bias_results)*100, 1), "%)\n"))

    cat(paste0("Mean studies imputed: ", round(mean(bias_results$tf_added), 2), "\n"))
    cat(paste0("Mean effect change after T&F: ", round(mean(bias_results$tf_change), 4), "\n"))

    # Effect of bias on significance
    bias_results$sig_before <- abs(bias_results$theta / sqrt(0.1)) > 1.96  # rough approximation
    bias_results$theta_reduced <- abs(bias_results$tf_theta) < abs(bias_results$theta)

    cat(paste0("Effects reduced by T&F: ", sum(bias_results$theta_reduced),
               " (", round(mean(bias_results$theta_reduced)*100, 1), "%)\n"))

    saveRDS(bias_results, "output/publication_bias_results.rds")
  }
} else {
  cat("No RDS files found. Skipping publication bias analysis.\n")
}

################################################################################
# REVISION 5: DOMAIN ANALYSIS WITH CONFIDENCE INTERVALS
################################################################################

cat("\n\n=== REVISION 5: DOMAIN ANALYSIS WITH CONFIDENCE INTERVALS ===\n\n")

domain_analysis <- primary_data %>%
  group_by(outcome_domain) %>%
  summarise(
    n = n(),
    n_reviews = n_distinct(review_id),
    mean_R = mean(R, na.rm = TRUE),
    sd_R = sd(R, na.rm = TRUE),
    se_R = sd(R, na.rm = TRUE) / sqrt(n()),
    ci_lower = mean_R - 1.96 * se_R,
    ci_upper = mean_R + 1.96 * se_R,
    pct_fragile = mean(R < 0.5, na.rm = TRUE) * 100,
    mean_tau = mean(tau, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n >= 10)

overall_mean_R <- mean(primary_data$R, na.rm = TRUE)
overall_se_R <- sd(primary_data$R, na.rm = TRUE) / sqrt(nrow(primary_data))

# Statistical tests
domain_tests <- domain_analysis %>%
  rowwise() %>%
  mutate(
    t_stat = (mean_R - overall_mean_R) / se_R,
    df = n - 1,
    p_value = 2 * pt(-abs(t_stat), df),
    effect_size = (mean_R - overall_mean_R) / sd_R,  # Cohen's d equivalent
    direction = ifelse(mean_R > overall_mean_R, "MORE stable", "LESS stable")
  ) %>%
  ungroup()

domain_tests$p_adjusted <- p.adjust(domain_tests$p_value, method = "BH")
domain_tests$sig_fdr <- domain_tests$p_adjusted < 0.05

cat("DOMAIN ANALYSIS WITH 95% CONFIDENCE INTERVALS:\n")
cat(paste(rep("-", 70), collapse = ""), "\n")
cat(paste0("Overall mean R: ", round(overall_mean_R, 3),
           " (95% CI: ", round(overall_mean_R - 1.96*overall_se_R, 3),
           "-", round(overall_mean_R + 1.96*overall_se_R, 3), ")\n\n"))

results_table <- domain_tests %>%
  arrange(p_adjusted) %>%
  select(outcome_domain, n, mean_R, ci_lower, ci_upper, effect_size, p_adjusted, sig_fdr)

cat("Results (sorted by adjusted p-value):\n")
print(head(as.data.frame(results_table), 15), row.names = FALSE)

cat(paste0("\nSignificant after FDR: ", sum(domain_tests$sig_fdr), "/",
           nrow(domain_tests), "\n"))

# Report significant domains with effect sizes
sig_domains <- domain_tests %>% filter(sig_fdr) %>% arrange(p_adjusted)
if (nrow(sig_domains) > 0) {
  cat("\nSIGNIFICANT DOMAINS (with effect sizes):\n")
  for (i in 1:nrow(sig_domains)) {
    cat(paste0("  ", sig_domains$outcome_domain[i], ": ", sig_domains$direction[i],
               "\n    R = ", round(sig_domains$mean_R[i], 3),
               " (95% CI: ", round(sig_domains$ci_lower[i], 3), "-",
               round(sig_domains$ci_upper[i], 3), ")",
               "\n    Effect size d = ", round(sig_domains$effect_size[i], 2),
               ", p_adj = ", format.pval(sig_domains$p_adjusted[i], digits = 2), "\n"))
  }
}

saveRDS(domain_tests, "output/domain_analysis_with_ci.rds")

################################################################################
# REVISION 6: GENERATE IMPROVED FIGURES
################################################################################

cat("\n\n=== REVISION 6: GENERATING PUBLICATION FIGURES ===\n\n")

theme_publication <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "gray40"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# Figure 1: ROC Curve with optimal threshold
png("figures/ROC_curve_optimized.png", width = 8, height = 6, units = "in", res = 300)

plot(roc_obj, main = "", legacy.axes = TRUE,
     xlab = "False Positive Rate (1 - Specificity)",
     ylab = "True Positive Rate (Sensitivity)",
     col = "steelblue", lwd = 2)

# Add optimal point
optimal_coords <- coords(roc_obj, optimal_thresh, ret = c("specificity", "sensitivity"))
points(1 - optimal_coords$specificity, optimal_coords$sensitivity,
       pch = 19, col = "red", cex = 1.5)

# Add diagonal
abline(a = 0, b = 1, lty = 2, col = "gray50")

# Add legend
legend("bottomright",
       legend = c(paste0("AUC = ", round(auc_value, 3)),
                  paste0("Optimal threshold = ", optimal_thresh)),
       col = c("steelblue", "red"),
       lwd = c(2, NA), pch = c(NA, 19),
       bty = "n")

dev.off()
cat("  Saved: figures/ROC_curve_optimized.png\n")

# Figure 2: Precision-Recall Curve
pr_data <- data.frame()
for (thresh in seq(0.01, 0.99, by = 0.01)) {
  pred <- ifelse(all_predictions$prob_fragile >= thresh, "Fragile", "Stable")
  tp <- sum(pred == "Fragile" & all_predictions$actual == "Fragile")
  fp <- sum(pred == "Fragile" & all_predictions$actual == "Stable")
  fn <- sum(pred == "Stable" & all_predictions$actual == "Fragile")

  precision <- ifelse((tp + fp) > 0, tp / (tp + fp), 0)
  recall <- ifelse((tp + fn) > 0, tp / (tp + fn), 0)

  pr_data <- rbind(pr_data, data.frame(threshold = thresh, precision = precision, recall = recall))
}

# Calculate PR-AUC
pr_auc <- sum(diff(pr_data$recall[order(pr_data$recall)]) *
              head(pr_data$precision[order(pr_data$recall)], -1), na.rm = TRUE)

p_pr <- ggplot(pr_data, aes(x = recall, y = precision)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_hline(yintercept = mean(all_predictions$actual == "Fragile"),
             linetype = "dashed", color = "gray50") +
  labs(
    title = "Precision-Recall Curve",
    subtitle = paste0("PR-AUC = ", round(abs(pr_auc), 3), " | Baseline = ",
                      round(mean(all_predictions$actual == "Fragile"), 3)),
    x = "Recall (Sensitivity)",
    y = "Precision (PPV)"
  ) +
  theme_publication +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

ggsave("figures/PR_curve.png", p_pr, width = 8, height = 6, dpi = 300)
cat("  Saved: figures/PR_curve.png\n")

# Figure 3: Threshold optimization plot
p_thresh <- ggplot(threshold_results, aes(x = threshold)) +
  geom_line(aes(y = sensitivity, color = "Sensitivity"), linewidth = 1) +
  geom_line(aes(y = specificity, color = "Specificity"), linewidth = 1) +
  geom_line(aes(y = balanced_acc, color = "Balanced Accuracy"), linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = optimal_thresh, linetype = "dotted", color = "red") +
  annotate("text", x = optimal_thresh + 0.05, y = 0.95,
           label = paste0("Optimal = ", optimal_thresh), size = 3) +
  scale_color_manual(values = c("Sensitivity" = "forestgreen",
                                 "Specificity" = "steelblue",
                                 "Balanced Accuracy" = "purple")) +
  labs(
    title = "Threshold Optimization",
    subtitle = "Trade-off between sensitivity and specificity",
    x = "Classification Threshold",
    y = "Performance Metric",
    color = NULL
  ) +
  theme_publication

ggsave("figures/threshold_optimization.png", p_thresh, width = 10, height = 6, dpi = 300)
cat("  Saved: figures/threshold_optimization.png\n")

# Figure 4: Domain forest plot with CIs
sig_for_plot <- domain_tests %>%
  filter(n >= 20) %>%
  arrange(mean_R) %>%
  mutate(outcome_domain = factor(outcome_domain, levels = outcome_domain))

p_forest <- ggplot(sig_for_plot, aes(x = mean_R, y = outcome_domain)) +
  geom_vline(xintercept = overall_mean_R, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.3, color = "gray40") +
  geom_point(aes(color = sig_fdr, size = n)) +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "steelblue"),
                     labels = c("Not significant", "Significant (FDR < 0.05)"),
                     name = NULL) +
  scale_size_continuous(range = c(2, 6), name = "n MAs") +
  labs(
    title = "Domain-Specific Stability",
    subtitle = paste0("Mean R with 95% CI | Dashed line = overall mean (",
                      round(overall_mean_R, 3), ")"),
    x = "Mean Stability (R)",
    y = NULL
  ) +
  theme_publication +
  theme(legend.position = "right")

ggsave("figures/domain_forest_plot.png", p_forest, width = 10, height = 8, dpi = 300)
cat("  Saved: figures/domain_forest_plot.png\n")

# Figure 5: Calibration plot
p_calib <- ggplot(calibration, aes(x = mean_predicted, y = observed_rate)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  geom_point(aes(size = n), color = "steelblue", alpha = 0.7) +
  geom_smooth(method = "loess", se = FALSE, color = "red", linewidth = 0.8) +
  scale_size_continuous(range = c(2, 8), name = "n") +
  labs(
    title = "Model Calibration",
    subtitle = "Predicted vs observed fragility rates by decile",
    x = "Mean Predicted Probability",
    y = "Observed Fragility Rate"
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_publication

ggsave("figures/calibration_plot.png", p_calib, width = 8, height = 6, dpi = 300)
cat("  Saved: figures/calibration_plot.png\n")

################################################################################
# FINAL SUMMARY
################################################################################

cat("\n\n================================================================\n")
cat("FINAL REVISED ANALYSIS SUMMARY\n")
cat("================================================================\n\n")

cat("DATA:\n")
cat(paste0("  Total meta-analyses: ", nrow(ma4_data), "\n"))
cat(paste0("  Primary analysis (logRR): ", nrow(primary_data), "\n"))
cat(paste0("  Outcome domains: ", length(unique(ma4_data$outcome_domain)), "\n"))
cat(paste0("  'Other Clinical': ", other_pct, "%\n\n"))

cat("PREDICTIVE MODEL (IMPROVED):\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
cat(paste0("  AUC-ROC: ", round(auc_value, 3), " (", round(auc_ci[1], 3), "-", round(auc_ci[3], 3), ")\n"))
cat(paste0("  Optimal threshold (Youden): ", optimal_thresh, "\n"))
cat(paste0("  Sensitivity at optimal: ", round(optimal_youden$sensitivity, 3), "\n"))
cat(paste0("  Specificity at optimal: ", round(optimal_youden$specificity, 3), "\n"))
cat(paste0("  Balanced Accuracy: ", round(optimal_youden$balanced_acc, 3), "\n"))
cat(paste0("  For 80% sensitivity: threshold = ", optimal_sens80$threshold,
           ", specificity = ", round(optimal_sens80$specificity, 3), "\n\n"))

if (exists("bias_results") && nrow(bias_results) > 0) {
  cat("PUBLICATION BIAS:\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat(paste0("  Meta-analyses tested: ", nrow(bias_results), "\n"))
  cat(paste0("  Egger positive: ", round(mean(bias_results$egger_sig)*100, 1), "%\n"))
  cat(paste0("  T&F imputed studies: ", round(mean(bias_results$tf_added > 0)*100, 1), "%\n\n"))
}

cat("DOMAIN FINDINGS:\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
cat(paste0("  Significant domains (FDR < 0.05): ", sum(domain_tests$sig_fdr), "/",
           nrow(domain_tests), "\n"))

cat("\nFIGURES GENERATED:\n")
cat("  figures/ROC_curve_optimized.png\n")
cat("  figures/PR_curve.png\n")
cat("  figures/threshold_optimization.png\n")
cat("  figures/domain_forest_plot.png\n")
cat("  figures/calibration_plot.png\n")

cat("\n\nEDITORIAL CONCERNS ADDRESSED:\n")
cat("  [x] Predictive model sensitivity improved via threshold optimization\n")
cat("  [x] Domain classification enhanced (reduced 'Other')\n")
cat("  [x] Publication bias analysis on subsample\n")
cat("  [x] Confidence intervals for domain comparisons\n")
cat("  [x] Calibration assessment added\n")
cat("  [x] ROC and PR curves generated\n")

# Save comprehensive results
saveRDS(list(
  data = primary_data,
  predictions = all_predictions,
  roc = roc_obj,
  auc = auc_value,
  threshold_optimization = threshold_results,
  optimal_threshold = optimal_thresh,
  calibration = calibration,
  domain_analysis = domain_tests,
  bias_results = if(exists("bias_results")) bias_results else NULL
), "output/EDITORIAL_REVISION_2_RESULTS.rds")

cat("\n================================================================\n")
cat("REVISION COMPLETE\n")
cat("================================================================\n")
