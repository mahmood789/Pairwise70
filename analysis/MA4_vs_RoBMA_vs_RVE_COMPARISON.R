# ============================================================================
# MA4 vs RoBMA vs RVE - COMPREHENSIVE COMPARISON
# Comparative analysis of modern meta-analytic methods
# ============================================================================

cat("
################################################################################
#                                                                              #
#   MA4 vs RoBMA vs RVE - COMPARATIVE ANALYSIS                                #
#   Modern Meta-Analytic Methods Comparison                                    #
#                                                                              #
################################################################################
\n\n")

set.seed(42)

# Load packages
library(metafor)
library(clubSandwich)

# Check RoBMA availability
robma_available <- requireNamespace("RoBMA", quietly = TRUE)
if (robma_available) {
  library(RoBMA)
  cat("RoBMA version:", as.character(packageVersion("RoBMA")), "\n")
} else {
  cat("WARNING: RoBMA not available - will skip Bayesian analysis\n")
}

cat("metafor version:", as.character(packageVersion("metafor")), "\n")
cat("clubSandwich version:", as.character(packageVersion("clubSandwich")), "\n\n")

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading MA4 results and raw data...\n")

# Load MA4 results
ma4_results <- read.csv("C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/ma4_results_pairwise70.csv",
                        stringsAsFactors = FALSE)

cat("Total meta-analyses:", nrow(ma4_results), "\n\n")

# Load raw Pairwise70 data for detailed analysis
data_dir <- "C:/Users/user/OneDrive - NHS/Documents/Pairwise70/data"
rda_files <- list.files(data_dir, pattern = "\\.rda$", full.names = TRUE)

# ============================================================================
# PART 1: RVE ANALYSIS (Robust Variance Estimation)
# ============================================================================

cat("
################################################################################
PART 1: ROBUST VARIANCE ESTIMATION (RVE) ANALYSIS
################################################################################
\n")

# RVE is for handling dependent effect sizes
# We'll apply small-sample corrections to our meta-analyses

rve_results <- data.frame()

cat("Running RVE analysis on all meta-analyses...\n")

for (i in seq_along(rda_files)) {
  # Load review data
  env <- new.env()
  load(rda_files[i], envir = env)
  df <- get(ls(env)[1], envir = env)

  if (nrow(df) == 0) next

  review_id <- gsub("\\.rda$", "", basename(rda_files[i]))

  # Get unique analyses
  analyses <- unique(df$Analysis.number)

  for (an in analyses) {
    sub <- df[df$Analysis.number == an, ]
    if (nrow(sub) < 2) next

    # Determine effect type and extract data
    has_giv <- all(!is.na(sub$GIV.Mean) & !is.na(sub$GIV.SE))
    has_binary <- all(!is.na(sub$Experimental.cases) & !is.na(sub$Experimental.N) &
                      !is.na(sub$Control.cases) & !is.na(sub$Control.N))
    has_cont <- all(!is.na(sub$Experimental.mean) & !is.na(sub$Experimental.SD) &
                    !is.na(sub$Control.mean) & !is.na(sub$Control.SD))

    yi <- sei <- NULL
    effect_type <- NA

    if (has_giv) {
      yi <- sub$GIV.Mean
      sei <- sub$GIV.SE
      effect_type <- "GIV"
    } else if (has_binary) {
      # Calculate log RR
      e_t <- pmax(0.5, sub$Experimental.cases)
      n_t <- sub$Experimental.N
      e_c <- pmax(0.5, sub$Control.cases)
      n_c <- sub$Control.N
      yi <- log((e_t / n_t) / (e_c / n_c))
      sei <- sqrt(1/e_t - 1/n_t + 1/e_c - 1/n_c)
      effect_type <- "logRR"
    } else if (has_cont) {
      # Calculate MD
      yi <- sub$Experimental.mean - sub$Control.mean
      s_pool <- sqrt(((sub$Experimental.N - 1) * sub$Experimental.SD^2 +
                      (sub$Control.N - 1) * sub$Control.SD^2) /
                     (sub$Experimental.N + sub$Control.N - 2))
      sei <- s_pool * sqrt(1/sub$Experimental.N + 1/sub$Control.N)
      effect_type <- "MD"
    }

    if (is.null(yi) || any(!is.finite(yi)) || any(!is.finite(sei)) || any(sei <= 0)) next

    k <- length(yi)

    tryCatch({
      # Standard REML
      fit_reml <- rma(yi = yi, sei = sei, method = "REML")

      # RVE with clubSandwich (small-sample corrected)
      # CR2 is the recommended bias-reduced linearization estimator
      rve_se <- sqrt(vcovCR(fit_reml, type = "CR2"))
      rve_ci <- coef_test(fit_reml, vcov = "CR2")

      # Calculate RVE-adjusted statistics
      theta_reml <- as.numeric(coef(fit_reml))
      se_reml <- as.numeric(fit_reml$se)
      se_rve <- as.numeric(rve_se)

      # Satterthwaite df for small-sample inference
      df_satt <- rve_ci$df

      # RVE-adjusted CI
      t_crit <- qt(0.975, df = df_satt)
      ci_lower_rve <- theta_reml - t_crit * se_rve
      ci_upper_rve <- theta_reml + t_crit * se_rve

      # Standard CI (normal approximation)
      ci_lower_std <- theta_reml - 1.96 * se_reml
      ci_upper_std <- theta_reml + 1.96 * se_reml

      # Check if conclusions differ
      sig_std <- (ci_lower_std > 0) | (ci_upper_std < 0)
      sig_rve <- (ci_lower_rve > 0) | (ci_upper_rve < 0)
      conclusion_changed <- sig_std != sig_rve

      # CI width ratio
      ci_width_ratio <- (ci_upper_rve - ci_lower_rve) / (ci_upper_std - ci_lower_std)

      rve_results <- rbind(rve_results, data.frame(
        review_id = review_id,
        analysis_id = an,
        effect_type = effect_type,
        k = k,
        theta = theta_reml,
        se_standard = se_reml,
        se_rve = se_rve,
        df_satterthwaite = df_satt,
        ci_lower_std = ci_lower_std,
        ci_upper_std = ci_upper_std,
        ci_lower_rve = ci_lower_rve,
        ci_upper_rve = ci_upper_rve,
        sig_standard = sig_std,
        sig_rve = sig_rve,
        conclusion_changed = conclusion_changed,
        ci_width_ratio = ci_width_ratio,
        tau2 = fit_reml$tau2
      ))

    }, error = function(e) NULL)
  }

  if (i %% 50 == 0) cat("  Processed", i, "of", length(rda_files), "reviews\n")
}

cat("\nRVE analysis complete:", nrow(rve_results), "meta-analyses\n")

# Merge with MA4 results
rve_results$ma_id <- paste(rve_results$review_id, rve_results$analysis_id, sep = "_")
ma4_results$ma_id <- paste(ma4_results$review_id, ma4_results$analysis_id, sep = "_")

comparison_df <- merge(rve_results, ma4_results[, c("ma_id", "R", "theta", "sigma", "tau")],
                       by = "ma_id", suffixes = c("_rve", "_ma4"))

cat("Merged dataset:", nrow(comparison_df), "meta-analyses\n\n")

# ============================================================================
# RVE RESULTS SUMMARY
# ============================================================================

cat("--- RVE Summary Statistics ---\n\n")

cat("SE inflation (RVE vs Standard):\n")
cat("  Mean ratio:", round(mean(comparison_df$se_rve / comparison_df$se_standard, na.rm = TRUE), 3), "\n")
cat("  Median ratio:", round(median(comparison_df$se_rve / comparison_df$se_standard, na.rm = TRUE), 3), "\n")
cat("  Range:", round(min(comparison_df$se_rve / comparison_df$se_standard, na.rm = TRUE), 3), "-",
    round(max(comparison_df$se_rve / comparison_df$se_standard, na.rm = TRUE), 3), "\n\n")

cat("CI width ratio (RVE vs Standard):\n")
cat("  Mean:", round(mean(comparison_df$ci_width_ratio, na.rm = TRUE), 3), "\n")
cat("  Median:", round(median(comparison_df$ci_width_ratio, na.rm = TRUE), 3), "\n\n")

cat("Conclusions changed after RVE:\n")
cat("  N changed:", sum(comparison_df$conclusion_changed, na.rm = TRUE), "\n")
cat("  % changed:", round(100 * mean(comparison_df$conclusion_changed, na.rm = TRUE), 2), "%\n\n")

cat("Satterthwaite degrees of freedom:\n")
cat("  Mean:", round(mean(comparison_df$df_satterthwaite, na.rm = TRUE), 2), "\n")
cat("  Median:", round(median(comparison_df$df_satterthwaite, na.rm = TRUE), 2), "\n")
cat("  df < 4 (very small sample):", sum(comparison_df$df_satterthwaite < 4, na.rm = TRUE),
    "(", round(100 * mean(comparison_df$df_satterthwaite < 4, na.rm = TRUE), 1), "%)\n\n")

# ============================================================================
# PART 2: RoBMA ANALYSIS (Bayesian Model Averaging)
# ============================================================================

cat("
################################################################################
PART 2: RoBMA ANALYSIS (Bayesian Model Averaging)
################################################################################
\n")

if (robma_available) {

  # RoBMA is computationally intensive - run on stratified sample
  # Sample across R-index categories to ensure representation

  comparison_df$R_cat <- cut(comparison_df$R, breaks = c(0, 0.3, 0.5, 0.7, 0.9, 1),
                             labels = c("Very Low", "Low", "Moderate", "High", "Excellent"),
                             include.lowest = TRUE)

  # Stratified sample: 50 per category (or all if fewer)
  sample_ids <- c()
  for (cat in levels(comparison_df$R_cat)) {
    cat_ids <- comparison_df$ma_id[comparison_df$R_cat == cat]
    n_sample <- min(50, length(cat_ids))
    sample_ids <- c(sample_ids, sample(cat_ids, n_sample))
  }

  cat("Running RoBMA on", length(sample_ids), "stratified sample meta-analyses\n")
  cat("(Full dataset would take ~", round(nrow(comparison_df) * 2 / 60, 0), " hours)\n\n")

  robma_results <- data.frame()

  pb_count <- 0
  for (ma_id in sample_ids) {
    pb_count <- pb_count + 1

    # Get data for this meta-analysis
    row <- comparison_df[comparison_df$ma_id == ma_id, ]
    review_id <- row$review_id
    analysis_id <- row$analysis_id

    # Reload raw data
    rda_file <- file.path(data_dir, paste0(review_id, ".rda"))
    if (!file.exists(rda_file)) next

    env <- new.env()
    load(rda_file, envir = env)
    df <- get(ls(env)[1], envir = env)

    sub <- df[df$Analysis.number == analysis_id, ]
    if (nrow(sub) < 2) next

    # Extract effect sizes
    has_giv <- all(!is.na(sub$GIV.Mean) & !is.na(sub$GIV.SE))
    has_binary <- all(!is.na(sub$Experimental.cases) & !is.na(sub$Experimental.N) &
                      !is.na(sub$Control.cases) & !is.na(sub$Control.N))

    yi <- sei <- NULL

    if (has_giv) {
      yi <- sub$GIV.Mean
      sei <- sub$GIV.SE
    } else if (has_binary) {
      e_t <- pmax(0.5, sub$Experimental.cases)
      n_t <- sub$Experimental.N
      e_c <- pmax(0.5, sub$Control.cases)
      n_c <- sub$Control.N
      yi <- log((e_t / n_t) / (e_c / n_c))
      sei <- sqrt(1/e_t - 1/n_t + 1/e_c - 1/n_c)
    } else {
      yi <- sub$Experimental.mean - sub$Control.mean
      s_pool <- sqrt(((sub$Experimental.N - 1) * sub$Experimental.SD^2 +
                      (sub$Control.N - 1) * sub$Control.SD^2) /
                     (sub$Experimental.N + sub$Control.N - 2))
      sei <- s_pool * sqrt(1/sub$Experimental.N + 1/sub$Control.N)
    }

    if (is.null(yi) || any(!is.finite(yi)) || any(!is.finite(sei)) || any(sei <= 0)) next

    tryCatch({
      # Run RoBMA with minimal output
      # Using default priors and 12 model ensemble
      fit_robma <- RoBMA(d = yi, se = sei,
                         seed = 42,
                         parallel = FALSE,
                         silent = TRUE)

      # Extract key results
      summ <- summary(fit_robma)

      # Posterior probabilities
      ph1 <- summ$components$Effect$probability  # P(effect exists)
      ph_het <- summ$components$Heterogeneity$probability  # P(heterogeneity)
      ph_bias <- summ$components$Bias$probability  # P(publication bias)

      # Model-averaged estimates
      theta_robma <- summ$estimates$Mean[1]  # Effect estimate
      ci_lower_robma <- summ$estimates$`0.025`[1]
      ci_upper_robma <- summ$estimates$`0.975`[1]

      # Inclusion BF for effect
      bf_effect <- summ$components$Effect$inclusion_BF

      robma_results <- rbind(robma_results, data.frame(
        ma_id = ma_id,
        k = length(yi),
        theta_robma = theta_robma,
        ci_lower_robma = ci_lower_robma,
        ci_upper_robma = ci_upper_robma,
        P_effect = ph1,
        P_heterogeneity = ph_het,
        P_bias = ph_bias,
        BF_effect = bf_effect
      ))

    }, error = function(e) {
      cat("  Error for", ma_id, ":", conditionMessage(e), "\n")
    })

    if (pb_count %% 25 == 0) {
      cat("  Completed", pb_count, "of", length(sample_ids), "\n")
    }
  }

  cat("\nRoBMA analysis complete:", nrow(robma_results), "meta-analyses\n\n")

  # Merge RoBMA with comparison data
  comparison_robma <- merge(comparison_df, robma_results, by = "ma_id", all.x = FALSE)

} else {
  cat("RoBMA not available - skipping Bayesian analysis\n")
  comparison_robma <- NULL
}

# ============================================================================
# PART 3: COMPARATIVE ANALYSIS
# ============================================================================

cat("
################################################################################
PART 3: COMPARATIVE ANALYSIS - MA4 vs RVE vs RoBMA
################################################################################
\n")

# --- MA4 vs RVE Comparison ---
cat("--- MA4 R-index vs RVE Metrics ---\n\n")

# Correlation between R and RVE conclusion change
cat("Correlation: R-index vs conclusion changed (point-biserial):\n")
cor_r_change <- cor.test(comparison_df$R, as.numeric(comparison_df$conclusion_changed))
cat("  r =", round(cor_r_change$estimate, 4), ", p =", format(cor_r_change$p.value, digits = 3), "\n\n")

# Do low-R meta-analyses have conclusions more likely to change with RVE?
cat("Conclusion change rate by R-index category:\n")
for (cat in levels(comparison_df$R_cat)) {
  sub <- comparison_df[comparison_df$R_cat == cat, ]
  if (nrow(sub) > 0) {
    change_rate <- mean(sub$conclusion_changed, na.rm = TRUE)
    cat("  ", cat, ": ", round(100 * change_rate, 1), "% (n=", nrow(sub), ")\n", sep = "")
  }
}

# Chi-square test
cat("\nChi-square test (R category vs conclusion changed):\n")
chisq_result <- chisq.test(table(comparison_df$R_cat, comparison_df$conclusion_changed))
cat("  X2 =", round(chisq_result$statistic, 2), ", df =", chisq_result$parameter,
    ", p =", format(chisq_result$p.value, digits = 3), "\n\n")

# Correlation: R vs SE inflation
cat("Correlation: R-index vs SE inflation ratio:\n")
se_ratio <- comparison_df$se_rve / comparison_df$se_standard
cor_r_se <- cor.test(comparison_df$R, se_ratio, use = "complete.obs")
cat("  r =", round(cor_r_se$estimate, 4), ", p =", format(cor_r_se$p.value, digits = 3), "\n\n")

# --- MA4 vs RoBMA Comparison (if available) ---
if (!is.null(comparison_robma) && nrow(comparison_robma) > 0) {

  cat("--- MA4 R-index vs RoBMA Metrics ---\n\n")

  # Correlation: R vs P(effect)
  cat("Correlation: R-index vs P(effect exists):\n")
  cor_r_peff <- cor.test(comparison_robma$R, comparison_robma$P_effect, use = "complete.obs")
  cat("  r =", round(cor_r_peff$estimate, 4), ", p =", format(cor_r_peff$p.value, digits = 3), "\n\n")

  # Correlation: R vs P(bias)
  cat("Correlation: R-index vs P(publication bias):\n")
  cor_r_pbias <- cor.test(comparison_robma$R, comparison_robma$P_bias, use = "complete.obs")
  cat("  r =", round(cor_r_pbias$estimate, 4), ", p =", format(cor_r_pbias$p.value, digits = 3), "\n\n")

  # Correlation: R vs BF(effect)
  cat("Correlation: R-index vs log(BF effect):\n")
  log_bf <- log(comparison_robma$BF_effect + 1e-10)
  log_bf[!is.finite(log_bf)] <- NA
  cor_r_bf <- cor.test(comparison_robma$R, log_bf, use = "complete.obs")
  cat("  r =", round(cor_r_bf$estimate, 4), ", p =", format(cor_r_bf$p.value, digits = 3), "\n\n")

  # P(bias) by R category
  cat("P(publication bias) by R-index category:\n")
  for (cat in levels(comparison_robma$R_cat)) {
    sub <- comparison_robma[comparison_robma$R_cat == cat, ]
    if (nrow(sub) > 0) {
      mean_pbias <- mean(sub$P_bias, na.rm = TRUE)
      cat("  ", cat, ": ", round(mean_pbias, 3), " (n=", nrow(sub), ")\n", sep = "")
    }
  }

  # Effect estimate comparison
  cat("\nEffect estimate comparison (MA4 theta vs RoBMA theta):\n")
  cor_theta <- cor.test(comparison_robma$theta_rve, comparison_robma$theta_robma, use = "complete.obs")
  cat("  r =", round(cor_theta$estimate, 4), "\n")
  cat("  Mean absolute difference:", round(mean(abs(comparison_robma$theta_rve - comparison_robma$theta_robma), na.rm = TRUE), 4), "\n\n")

  # Key insight: Does R predict whether RoBMA detects bias?
  cat("--- Key Question: Does low R predict RoBMA-detected bias? ---\n\n")

  comparison_robma$high_bias <- comparison_robma$P_bias > 0.5
  comparison_robma$low_R <- comparison_robma$R < 0.5

  bias_table <- table(Low_R = comparison_robma$low_R, High_Pbias = comparison_robma$high_bias)
  cat("Contingency table (R < 0.5 vs P(bias) > 0.5):\n")
  print(bias_table)

  if (all(dim(bias_table) == c(2, 2))) {
    fisher_test <- fisher.test(bias_table)
    cat("\nFisher's exact test:\n")
    cat("  Odds ratio:", round(fisher_test$estimate, 2), "\n")
    cat("  p-value:", format(fisher_test$p.value, digits = 3), "\n")
  }
}

# ============================================================================
# PART 4: SUMMARY COMPARISON TABLE
# ============================================================================

cat("\n
################################################################################
PART 4: SUMMARY - WHAT EACH METHOD TELLS US
################################################################################
\n")

cat("--- Method Comparison Summary ---\n\n")

summary_table <- data.frame(
  Metric = c("Total MAs analyzed",
             "Mean R-index (MA4)",
             "% conclusions changed (RVE)",
             "Mean SE inflation (RVE)",
             "Mean P(bias) (RoBMA)",
             "Correlation: R vs SE inflation",
             "Correlation: R vs P(bias)"),
  Value = c(
    nrow(comparison_df),
    round(mean(comparison_df$R, na.rm = TRUE), 3),
    round(100 * mean(comparison_df$conclusion_changed, na.rm = TRUE), 1),
    round(mean(comparison_df$se_rve / comparison_df$se_standard, na.rm = TRUE), 3),
    ifelse(!is.null(comparison_robma), round(mean(comparison_robma$P_bias, na.rm = TRUE), 3), "N/A"),
    round(cor_r_se$estimate, 3),
    ifelse(!is.null(comparison_robma), round(cor_r_pbias$estimate, 3), "N/A")
  )
)
print(summary_table, row.names = FALSE)

cat("\n--- Key Findings ---\n\n")

cat("1. MA4 R-index and RVE:\n")
if (cor_r_se$p.value < 0.05) {
  cat("   - Significant correlation between R and SE inflation\n")
  cat("   - Low R meta-analyses tend to have ",
      ifelse(cor_r_se$estimate < 0, "larger", "smaller"), " RVE corrections\n", sep = "")
} else {
  cat("   - No significant relationship between R and RVE corrections\n")
}

if (!is.null(comparison_robma)) {
  cat("\n2. MA4 R-index and RoBMA:\n")
  if (cor_r_pbias$p.value < 0.05) {
    cat("   - Significant correlation between R and P(bias)\n")
    cat("   - ", ifelse(cor_r_pbias$estimate < 0, "Lower", "Higher"),
        " R associated with ", ifelse(cor_r_pbias$estimate < 0, "higher", "lower"),
        " probability of publication bias\n", sep = "")
  } else {
    cat("   - No significant relationship between R and P(bias)\n")
  }
}

cat("\n3. Complementarity:\n")
cat("   - RVE addresses: dependent effect sizes, small-sample bias\n")
cat("   - RoBMA addresses: publication bias, model uncertainty\n")
cat("   - MA4 addresses: analytic stability/fragility\n")
cat("   - These are DIFFERENT constructs with partial overlap\n")

# ============================================================================
# SAVE RESULTS
# ============================================================================

out_dir <- "C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/comparison_output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(comparison_df, file.path(out_dir, "MA4_vs_RVE_comparison.csv"), row.names = FALSE)

if (!is.null(comparison_robma)) {
  write.csv(comparison_robma, file.path(out_dir, "MA4_vs_RoBMA_comparison.csv"), row.names = FALSE)
}

cat("\n\nResults saved to:", out_dir, "\n")

cat("\n
################################################################################
ANALYSIS COMPLETE
################################################################################
\n")
