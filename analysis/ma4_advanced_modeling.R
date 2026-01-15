# ============================================================================
# MA4 v1.0.1 ADVANCED STATISTICAL MODELING
# Beta regression, Quantile regression, Validation against metafor
# ============================================================================

cat("
================================================================================
   MA4 v1.0.1 ADVANCED STATISTICAL MODELING
   Beta Regression, Quantile Regression, metafor Validation
================================================================================
\n")

# Load data
results <- read.csv("C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/ma4_results_pairwise70.csv",
                    stringsAsFactors = FALSE)

cat("Loaded", nrow(results), "meta-analyses\n\n")

# Prepare modeling data
model_data <- results[complete.cases(results[, c("R", "k", "theta", "sigma", "tau", "effect_type")]), ]
model_data <- model_data[is.finite(model_data$theta) & is.finite(model_data$sigma) & model_data$sigma > 0, ]
model_data$log_k <- log(model_data$k)
model_data$abs_theta <- abs(model_data$theta)
model_data$log_tau <- log(model_data$tau + 1e-6)
model_data$theta_near_zero <- ifelse(abs(model_data$theta) < 0.1, 1, 0)
model_data$I2_proxy <- model_data$tau^2 / (model_data$tau^2 + model_data$sigma^2)

# Transform R for beta regression (avoid 0 and 1 boundaries)
eps <- 1e-6
model_data$R_beta <- pmax(eps, pmin(1 - eps, model_data$R))

# ============================================================================
# SECTION 1: BETA REGRESSION (for bounded [0,1] response)
# ============================================================================

cat("
================================================================================
SECTION 1: BETA REGRESSION
(Appropriate for bounded [0,1] response variable)
================================================================================
\n")

# Check if betareg is available
if (requireNamespace("betareg", quietly = TRUE)) {
  library(betareg)

  cat("--- 1.1 Beta Regression: Full Model ---\n")
  beta_m1 <- betareg(R_beta ~ log_k + effect_type + log_tau + abs_theta + theta_near_zero,
                     data = model_data)
  cat("\nModel Summary:\n")
  print(summary(beta_m1))

  cat("\n--- 1.2 Beta Regression: Pseudo R-squared ---\n")
  null_ll <- logLik(betareg(R_beta ~ 1, data = model_data))
  full_ll <- logLik(beta_m1)
  pseudo_r2 <- 1 - as.numeric(full_ll) / as.numeric(null_ll)
  cat("McFadden Pseudo R-squared:", round(pseudo_r2, 4), "\n")

  cat("\n--- 1.3 Beta Regression: Marginal Effects ---\n")
  # Average marginal effects (approximate)
  coefs <- coef(beta_m1)
  cat("Coefficient interpretation (log-odds scale):\n")
  for (i in 2:length(coefs)) {
    cat("  ", names(coefs)[i], ": ", round(coefs[i], 4), "\n", sep = "")
  }

  cat("\n--- 1.4 Beta vs OLS Comparison ---\n")
  ols_m <- lm(R ~ log_k + effect_type + log_tau + abs_theta + theta_near_zero, data = model_data)

  ols_pred <- predict(ols_m)
  beta_pred <- predict(beta_m1, type = "response")

  cat("OLS Predictions: min=", round(min(ols_pred), 4),
      ", max=", round(max(ols_pred), 4), "\n")
  cat("Beta Predictions: min=", round(min(beta_pred), 4),
      ", max=", round(max(beta_pred), 4), "\n")
  cat("Correlation OLS vs Beta predictions:", round(cor(ols_pred, beta_pred), 4), "\n")

} else {
  cat("betareg package not installed. Skipping beta regression.\n")
  cat("Install with: install.packages('betareg')\n")
}

# ============================================================================
# SECTION 2: QUANTILE REGRESSION
# ============================================================================

cat("\n
================================================================================
SECTION 2: QUANTILE REGRESSION
(Modeling different quantiles of R distribution)
================================================================================
\n")

if (requireNamespace("quantreg", quietly = TRUE)) {
  library(quantreg)

  cat("--- 2.1 Quantile Regression at tau = 0.25, 0.50, 0.75 ---\n")

  # 25th percentile
  qr_25 <- rq(R ~ log_k + effect_type + log_tau + theta_near_zero,
              data = model_data, tau = 0.25)
  cat("\n25th Percentile (tau=0.25):\n")
  print(round(summary(qr_25, se = "boot", R = 200)$coefficients, 4))

  # 50th percentile (median)
  qr_50 <- rq(R ~ log_k + effect_type + log_tau + theta_near_zero,
              data = model_data, tau = 0.50)
  cat("\n50th Percentile (median, tau=0.50):\n")
  print(round(summary(qr_50, se = "boot", R = 200)$coefficients, 4))

  # 75th percentile
  qr_75 <- rq(R ~ log_k + effect_type + log_tau + theta_near_zero,
              data = model_data, tau = 0.75)
  cat("\n75th Percentile (tau=0.75):\n")
  print(round(summary(qr_75, se = "boot", R = 200)$coefficients, 4))

  cat("\n--- 2.2 Coefficient Comparison Across Quantiles ---\n")
  coef_comparison <- data.frame(
    tau_25 = coef(qr_25),
    tau_50 = coef(qr_50),
    tau_75 = coef(qr_75),
    OLS = coef(lm(R ~ log_k + effect_type + log_tau + theta_near_zero, data = model_data))
  )
  print(round(coef_comparison, 4))

} else {
  cat("quantreg package not installed. Skipping quantile regression.\n")
  cat("Install with: install.packages('quantreg')\n")
}

# ============================================================================
# SECTION 3: VALIDATION AGAINST METAFOR
# ============================================================================

cat("\n
================================================================================
SECTION 3: VALIDATION AGAINST METAFOR
(Compare MA4 estimates with metafor's rma())
================================================================================
\n")

if (requireNamespace("metafor", quietly = TRUE)) {
  library(metafor)

  # Load original data to compute metafor estimates
  data_dir <- "C:/Users/user/OneDrive - NHS/Documents/Pairwise70/data"
  rda_files <- list.files(data_dir, pattern = "\\.rda$", full.names = TRUE)

  cat("Validating MA4 against metafor on sample of meta-analyses...\n\n")

  # Sample 100 meta-analyses for validation
  set.seed(42)
  sample_idx <- sample(1:nrow(results), min(200, nrow(results)))

  validation_results <- list()
  n_validated <- 0
  n_matched <- 0

  for (i in sample_idx) {
    row <- results[i, ]

    # Find and load the corresponding .rda file
    file_pattern <- paste0(row$review_id, "_data.rda")
    file_match <- grep(file_pattern, rda_files, value = TRUE)

    if (length(file_match) == 0) next

    e <- new.env()
    tryCatch({
      load(file_match[1], envir = e)
    }, error = function(err) next)

    objs <- ls(envir = e)
    if (length(objs) == 0) next

    df <- get(objs[1], envir = e)
    if (!is.data.frame(df)) next

    # Get the specific analysis
    if ("Analysis.number" %in% names(df)) {
      analysis_df <- df[df$Analysis.number == as.numeric(row$analysis_number), ]
    } else {
      analysis_df <- df
    }

    if (nrow(analysis_df) < 2) next

    # Filter for main analysis
    if ("Applicability" %in% names(analysis_df)) {
      main_df <- analysis_df[grepl("OVERALL|^$", analysis_df$Applicability, ignore.case = TRUE), ]
      if (nrow(main_df) < 2) main_df <- analysis_df
    } else {
      main_df <- analysis_df
    }

    if ("Study" %in% names(main_df)) {
      main_df <- main_df[!duplicated(main_df$Study), ]
    }

    # Extract y, se based on effect type
    y <- se <- NULL

    if (row$effect_type == "GIV" && "GIV.Mean" %in% names(main_df) && "GIV.SE" %in% names(main_df)) {
      y <- as.numeric(main_df$GIV.Mean)
      se <- as.numeric(main_df$GIV.SE)
    } else if (row$effect_type == "logRR" &&
               all(c("Experimental.cases", "Experimental.N", "Control.cases", "Control.N") %in% names(main_df))) {
      et <- as.numeric(main_df$Experimental.cases)
      nt <- as.numeric(main_df$Experimental.N)
      ec <- as.numeric(main_df$Control.cases)
      nc <- as.numeric(main_df$Control.N)

      valid <- !is.na(et) & !is.na(nt) & !is.na(ec) & !is.na(nc) &
               nt > 0 & nc > 0 & et >= 0 & ec >= 0 & et <= nt & ec <= nc
      if (sum(valid) < 2) next

      et <- et[valid]; nt <- nt[valid]; ec <- ec[valid]; nc <- nc[valid]
      a <- et; b <- nt - et; c <- ec; d <- nc - ec
      need_cc <- (a == 0) | (b == 0) | (c == 0) | (d == 0)
      cc <- ifelse(need_cc, 0.5, 0.0)
      a <- a + cc; b <- b + cc; c <- c + cc; d <- d + cc
      rr <- (a / (a + b)) / (c / (c + d))
      y <- log(rr)
      se <- sqrt(1/a - 1/(a + b) + 1/c - 1/(c + d))
    }

    if (is.null(y) || is.null(se)) next

    valid <- is.finite(y) & is.finite(se) & se > 0
    y <- y[valid]; se <- se[valid]

    if (length(y) < 2) next
    if (length(y) != row$k) next  # Ensure same k

    # Run metafor
    tryCatch({
      mf_result <- rma(yi = y, sei = se, method = "REML")

      validation_results[[length(validation_results) + 1]] <- data.frame(
        review_id = row$review_id,
        analysis_number = row$analysis_number,
        k = row$k,
        ma4_theta = row$theta,
        mf_theta = as.numeric(mf_result$beta),
        ma4_se = row$sigma,
        mf_se = mf_result$se,
        ma4_tau = row$tau,
        mf_tau = sqrt(mf_result$tau2),
        theta_diff = abs(row$theta - as.numeric(mf_result$beta)),
        se_diff = abs(row$sigma - mf_result$se),
        tau_diff = abs(row$tau - sqrt(mf_result$tau2))
      )
      n_validated <- n_validated + 1

      if (abs(row$theta - as.numeric(mf_result$beta)) < 0.01 &&
          abs(row$sigma - mf_result$se) < 0.01) {
        n_matched <- n_matched + 1
      }

    }, error = function(err) NULL)

    rm(e)
    if (n_validated >= 100) break
  }

  if (length(validation_results) > 0) {
    validation_df <- do.call(rbind, validation_results)

    cat("--- 3.1 Validation Summary ---\n")
    cat("Meta-analyses validated:", n_validated, "\n")
    cat("Exact matches (diff < 0.01):", n_matched, "(", round(100*n_matched/n_validated, 1), "%)\n\n")

    cat("--- 3.2 Theta (Effect Size) Comparison ---\n")
    cat("Correlation MA4 vs metafor:", round(cor(validation_df$ma4_theta, validation_df$mf_theta), 6), "\n")
    cat("Mean absolute difference:", round(mean(validation_df$theta_diff), 6), "\n")
    cat("Max absolute difference:", round(max(validation_df$theta_diff), 6), "\n")
    cat("RMSE:", round(sqrt(mean(validation_df$theta_diff^2)), 6), "\n\n")

    cat("--- 3.3 SE (Sigma) Comparison ---\n")
    cat("Correlation MA4 vs metafor:", round(cor(validation_df$ma4_se, validation_df$mf_se), 6), "\n")
    cat("Mean absolute difference:", round(mean(validation_df$se_diff), 6), "\n")
    cat("Max absolute difference:", round(max(validation_df$se_diff), 6), "\n\n")

    cat("--- 3.4 Tau (Heterogeneity) Comparison ---\n")
    cat("Correlation MA4 vs metafor:", round(cor(validation_df$ma4_tau, validation_df$mf_tau), 6), "\n")
    cat("Mean absolute difference:", round(mean(validation_df$tau_diff), 6), "\n")
    cat("Max absolute difference:", round(max(validation_df$tau_diff), 6), "\n\n")

    cat("--- 3.5 Sample Comparison (first 10) ---\n")
    print(head(validation_df[, c("review_id", "k", "ma4_theta", "mf_theta", "theta_diff",
                                  "ma4_tau", "mf_tau", "tau_diff")], 10))

    # Save validation results
    write.csv(validation_df,
              "C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/ma4_metafor_validation.csv",
              row.names = FALSE)
    cat("\nValidation results saved to: ma4_metafor_validation.csv\n")
  }

} else {
  cat("metafor package not installed. Skipping validation.\n")
  cat("Install with: install.packages('metafor')\n")
}

# ============================================================================
# SECTION 4: HIERARCHICAL/MIXED EFFECTS MODEL
# ============================================================================

cat("\n
================================================================================
SECTION 4: HIERARCHICAL MODEL (Random effects for reviews)
================================================================================
\n")

if (requireNamespace("lme4", quietly = TRUE)) {
  library(lme4)

  cat("--- 4.1 Mixed Effects Model: R ~ fixed effects + (1|review_id) ---\n")

  # Add review_id as factor
  model_data$review_id_factor <- factor(model_data$review_id)

  # Fit mixed model
  mixed_m <- lmer(R ~ log_k + effect_type + log_tau + theta_near_zero +
                    (1 | review_id_factor), data = model_data)

  cat("\nFixed Effects:\n")
  print(round(summary(mixed_m)$coefficients, 4))

  cat("\nRandom Effects Variance:\n")
  print(VarCorr(mixed_m))

  # ICC
  vc <- as.data.frame(VarCorr(mixed_m))
  icc <- vc$vcov[1] / sum(vc$vcov)
  cat("\nIntraclass Correlation (ICC):", round(icc, 4), "\n")
  cat("Interpretation:", round(icc * 100, 1), "% of variance in R is between-review\n")

  # Compare with OLS
  cat("\n--- 4.2 Comparison: Mixed vs OLS ---\n")
  ols_m <- lm(R ~ log_k + effect_type + log_tau + theta_near_zero, data = model_data)
  cat("OLS R-squared:", round(summary(ols_m)$r.squared, 4), "\n")

  if (requireNamespace("MuMIn", quietly = TRUE)) {
    rsq <- MuMIn::r.squaredGLMM(mixed_m)
    cat("Marginal R-squared (fixed effects only):", round(rsq[1], 4), "\n")
    cat("Conditional R-squared (fixed + random):", round(rsq[2], 4), "\n")
  }

} else {
  cat("lme4 package not installed. Skipping hierarchical model.\n")
  cat("Install with: install.packages('lme4')\n")
}

# ============================================================================
# SECTION 5: ROBUSTNESS CHECKS
# ============================================================================

cat("\n
================================================================================
SECTION 5: ROBUSTNESS CHECKS
================================================================================
\n")

# 5.1 Outlier sensitivity
cat("--- 5.1 Outlier Sensitivity Analysis ---\n")

# Identify outliers (|residual| > 3 SD)
ols_full <- lm(R ~ log_k + effect_type + log_tau + abs_theta + theta_near_zero,
               data = model_data)
resid <- rstandard(ols_full)
outliers <- which(abs(resid) > 3)
cat("Observations with |standardized residual| > 3:", length(outliers), "\n")

# Refit without outliers
if (length(outliers) > 0 && length(outliers) < nrow(model_data) * 0.05) {
  ols_no_outliers <- lm(R ~ log_k + effect_type + log_tau + abs_theta + theta_near_zero,
                        data = model_data[-outliers, ])
  cat("\nFull model R-squared:", round(summary(ols_full)$r.squared, 4), "\n")
  cat("Without outliers R-squared:", round(summary(ols_no_outliers)$r.squared, 4), "\n")

  cat("\nCoefficient comparison:\n")
  coef_comp <- data.frame(
    Full = coef(ols_full),
    No_Outliers = coef(ols_no_outliers),
    Diff = coef(ols_full) - coef(ols_no_outliers)
  )
  print(round(coef_comp, 4))
}

# 5.2 Effect type stratified analysis
cat("\n--- 5.2 Stratified Analysis by Effect Type ---\n")

for (et in unique(model_data$effect_type)) {
  sub <- model_data[model_data$effect_type == et, ]
  if (nrow(sub) > 50) {
    m <- lm(R ~ log_k + log_tau + theta_near_zero, data = sub)
    cat("\n", et, " (n=", nrow(sub), "):\n", sep = "")
    cat("  R-squared:", round(summary(m)$r.squared, 4), "\n")
    cat("  Coefficients:\n")
    print(round(coef(m), 4))
  }
}

# 5.3 Sample size sensitivity
cat("\n--- 5.3 Sensitivity by Sample Size (k) ---\n")

for (k_thresh in c(3, 5, 10)) {
  sub <- model_data[model_data$k >= k_thresh, ]
  m <- lm(R ~ log_k + effect_type + log_tau + theta_near_zero, data = sub)
  cat("k >=", k_thresh, "(n=", nrow(sub), "): R-squared =",
      round(summary(m)$r.squared, 4), "\n")
}

# ============================================================================
# SECTION 6: PREDICTIVE MODEL FOR R CATEGORIES
# ============================================================================

cat("\n
================================================================================
SECTION 6: CLASSIFICATION MODEL FOR R CATEGORIES
================================================================================
\n")

# Create R categories
model_data$R_class <- cut(model_data$R,
                          breaks = c(0, 0.5, 0.8, 1.0),
                          labels = c("Low", "Medium", "High"),
                          include.lowest = TRUE)

cat("--- 6.1 R Category Distribution ---\n")
print(table(model_data$R_class))

# Multinomial logistic regression
if (requireNamespace("nnet", quietly = TRUE)) {
  library(nnet)

  cat("\n--- 6.2 Multinomial Logistic Regression ---\n")
  multi_m <- multinom(R_class ~ log_k + effect_type + log_tau + theta_near_zero,
                      data = model_data, trace = FALSE)

  cat("\nCoefficients (vs Low baseline):\n")
  print(round(coef(multi_m), 4))

  cat("\nZ-values:\n")
  z <- summary(multi_m)$coefficients / summary(multi_m)$standard.errors
  print(round(z, 2))

  # Prediction accuracy
  pred_class <- predict(multi_m)
  accuracy <- mean(pred_class == model_data$R_class)
  cat("\nTraining accuracy:", round(accuracy, 4), "\n")

  cat("\nConfusion matrix:\n")
  print(table(Predicted = pred_class, Actual = model_data$R_class))

} else {
  cat("nnet package not installed. Skipping multinomial regression.\n")
}

# ============================================================================
# SECTION 7: FINAL SUMMARY
# ============================================================================

cat("\n
================================================================================
SECTION 7: FINAL MODEL SUMMARY
================================================================================
\n")

cat("BEST LINEAR MODEL (OLS):\n")
best_ols <- lm(R ~ log_k + effect_type + log_tau + abs_theta + theta_near_zero +
                 log_k:effect_type, data = model_data)
cat("R-squared:", round(summary(best_ols)$r.squared, 4), "\n")
cat("Adjusted R-squared:", round(summary(best_ols)$adj.r.squared, 4), "\n")
cat("Residual SE:", round(summary(best_ols)$sigma, 4), "\n")
cat("\nCoefficients:\n")
print(round(summary(best_ols)$coefficients, 4))

cat("\n
================================================================================
ADVANCED MODELING COMPLETE
================================================================================
\n")
