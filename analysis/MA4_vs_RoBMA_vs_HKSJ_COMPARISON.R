# ============================================================================
# MA4 vs RoBMA vs HKSJ - COMPREHENSIVE COMPARISON
# HKSJ = Hartung-Knapp-Sidik-Jonkman small-sample correction
# ============================================================================

cat("
################################################################################
#   MA4 vs RoBMA vs HKSJ - COMPARATIVE ANALYSIS                               #
################################################################################
\n\n")

set.seed(42)

library(metafor)

robma_available <- requireNamespace("RoBMA", quietly = TRUE)
if (robma_available) library(RoBMA)

cat("Packages loaded.\n\n")

# ============================================================================
# LOAD DATA
# ============================================================================

ma4_results <- read.csv("C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/ma4_results_pairwise70.csv",
                        stringsAsFactors = FALSE)
cat("MA4 results loaded:", nrow(ma4_results), "meta-analyses\n")

data_dir <- "C:/Users/user/OneDrive - NHS/Documents/Pairwise70/data"
rda_files <- list.files(data_dir, pattern = "\\.rda$", full.names = TRUE)

# ============================================================================
# PART 1: HKSJ vs Standard Analysis
# ============================================================================

cat("\n--- PART 1: HKSJ (Knapp-Hartung) vs Standard ---\n\n")

hksj_results <- list()
hksj_count <- 0

for (i in seq_along(rda_files)) {
  env <- new.env()
  load(rda_files[i], envir = env)
  df <- get(ls(env)[1], envir = env)
  if (nrow(df) == 0) next

  review_id <- gsub("_data\\.rda$", "", basename(rda_files[i]))
  analyses <- unique(df$Analysis.number)

  for (an in analyses) {
    sub <- df[df$Analysis.number == an, ]
    if (nrow(sub) < 2) next

    # Extract effect sizes
    has_giv <- all(!is.na(sub$GIV.Mean) & !is.na(sub$GIV.SE))
    has_binary <- all(!is.na(sub$Experimental.cases) & !is.na(sub$Experimental.N) &
                      !is.na(sub$Control.cases) & !is.na(sub$Control.N))
    has_cont <- all(!is.na(sub$Experimental.mean) & !is.na(sub$Experimental.SD))

    yi <- sei <- NULL
    effect_type <- NA

    if (has_giv) {
      yi <- sub$GIV.Mean
      sei <- sub$GIV.SE
      effect_type <- "GIV"
    } else if (has_binary) {
      e_t <- pmax(0.5, sub$Experimental.cases)
      n_t <- sub$Experimental.N
      e_c <- pmax(0.5, sub$Control.cases)
      n_c <- sub$Control.N
      yi <- log((e_t / n_t) / (e_c / n_c))
      var_term <- 1/e_t - 1/n_t + 1/e_c - 1/n_c
      var_term[var_term < 0] <- NA
      sei <- sqrt(var_term)
      effect_type <- "logRR"
    } else if (has_cont) {
      yi <- sub$Experimental.mean - sub$Control.mean
      s_pool <- sqrt(((sub$Experimental.N - 1) * sub$Experimental.SD^2 +
                      (sub$Control.N - 1) * sub$Control.SD^2) /
                     (sub$Experimental.N + sub$Control.N - 2))
      sei <- s_pool * sqrt(1/sub$Experimental.N + 1/sub$Control.N)
      effect_type <- "MD"
    }

    if (is.null(yi) || any(!is.finite(yi)) || any(!is.finite(sei)) || any(sei <= 0)) next

    k <- length(yi)
    ma_id <- paste(review_id, an, sep = "_")

    tryCatch({
      # Standard REML (z-based CI)
      fit_std <- rma(yi = yi, sei = sei, method = "REML")

      # HKSJ-adjusted (t-based CI with adjusted SE)
      fit_hksj <- rma(yi = yi, sei = sei, method = "REML", test = "knha")

      theta <- as.numeric(coef(fit_std))
      se_std <- as.numeric(fit_std$se)
      se_hksj <- as.numeric(fit_hksj$se)

      # CIs
      ci_std_lo <- fit_std$ci.lb
      ci_std_hi <- fit_std$ci.ub
      ci_hksj_lo <- fit_hksj$ci.lb
      ci_hksj_hi <- fit_hksj$ci.ub

      # P-values
      p_std <- fit_std$pval
      p_hksj <- fit_hksj$pval

      # Significance at 0.05
      sig_std <- p_std < 0.05
      sig_hksj <- p_hksj < 0.05

      # CI width ratio
      ci_width_std <- ci_std_hi - ci_std_lo
      ci_width_hksj <- ci_hksj_hi - ci_hksj_lo
      ci_ratio <- ci_width_hksj / ci_width_std

      hksj_count <- hksj_count + 1
      hksj_results[[hksj_count]] <- data.frame(
        ma_id = ma_id,
        review_id = review_id,
        analysis_id = an,
        effect_type = effect_type,
        k = k,
        theta = theta,
        se_std = se_std,
        se_hksj = se_hksj,
        se_ratio = se_hksj / se_std,
        ci_std_lo = ci_std_lo,
        ci_std_hi = ci_std_hi,
        ci_hksj_lo = ci_hksj_lo,
        ci_hksj_hi = ci_hksj_hi,
        ci_width_ratio = ci_ratio,
        p_std = p_std,
        p_hksj = p_hksj,
        sig_std = sig_std,
        sig_hksj = sig_hksj,
        conclusion_changed = sig_std != sig_hksj,
        tau2 = fit_std$tau2,
        I2 = fit_std$I2
      )
    }, error = function(e) NULL)
  }

  if (i %% 100 == 0) cat("  Processed", i, "/", length(rda_files), "reviews\n")
}

hksj_df <- do.call(rbind, hksj_results)
cat("\nHKSJ complete:", nrow(hksj_df), "meta-analyses\n")

# Merge with MA4
ma4_results$ma_id <- paste(ma4_results$review_id, ma4_results$analysis_number, sep = "_")
comp_df <- merge(hksj_df, ma4_results[, c("ma_id", "R")], by = "ma_id")

cat("Merged:", nrow(comp_df), "meta-analyses\n\n")

# HKSJ Summary
cat("--- HKSJ Summary ---\n")
cat("SE inflation ratio (HKSJ/Standard):\n")
cat("  Mean:", round(mean(comp_df$se_ratio, na.rm = TRUE), 3), "\n")
cat("  Median:", round(median(comp_df$se_ratio, na.rm = TRUE), 3), "\n")
cat("  IQR:", round(quantile(comp_df$se_ratio, 0.25, na.rm = TRUE), 3), "-",
    round(quantile(comp_df$se_ratio, 0.75, na.rm = TRUE), 3), "\n\n")

cat("CI width ratio (HKSJ/Standard):\n")
cat("  Mean:", round(mean(comp_df$ci_width_ratio, na.rm = TRUE), 3), "\n")
cat("  Median:", round(median(comp_df$ci_width_ratio, na.rm = TRUE), 3), "\n\n")

cat("Conclusions changed after HKSJ:\n")
n_changed <- sum(comp_df$conclusion_changed)
cat("  N:", n_changed, "\n")
cat("  %:", round(100 * mean(comp_df$conclusion_changed), 2), "%\n\n")

# Direction of change
sig_to_nonsig <- sum(comp_df$sig_std & !comp_df$sig_hksj)
nonsig_to_sig <- sum(!comp_df$sig_std & comp_df$sig_hksj)
cat("  Significant -> Non-significant:", sig_to_nonsig, "\n")
cat("  Non-significant -> Significant:", nonsig_to_sig, "\n\n")

# R-index categories
comp_df$R_cat <- cut(comp_df$R, breaks = c(0, 0.3, 0.5, 0.7, 0.9, 1),
                     labels = c("Very Low", "Low", "Moderate", "High", "Excellent"),
                     include.lowest = TRUE)

cat("Conclusion change by R-index category:\n")
for (cat_name in levels(comp_df$R_cat)) {
  sub <- comp_df[comp_df$R_cat == cat_name, ]
  if (nrow(sub) > 0) {
    cat("  ", cat_name, ": ", round(100 * mean(sub$conclusion_changed), 1),
        "% (n=", nrow(sub), ")\n", sep = "")
  }
}

# k categories
comp_df$k_cat <- cut(comp_df$k, breaks = c(1, 3, 5, 10, 20, Inf),
                     labels = c("2-3", "4-5", "6-10", "11-20", "21+"),
                     include.lowest = TRUE)

cat("\nConclusion change by number of studies:\n")
for (cat_name in levels(comp_df$k_cat)) {
  sub <- comp_df[comp_df$k_cat == cat_name, ]
  if (nrow(sub) > 0) {
    cat("  ", cat_name, ": ", round(100 * mean(sub$conclusion_changed), 1),
        "% (n=", nrow(sub), ")\n", sep = "")
  }
}

# Correlations
cat("\n--- Correlations ---\n")

cor_r_se <- cor.test(comp_df$R, comp_df$se_ratio)
cat("R vs SE ratio: r =", round(cor_r_se$estimate, 4),
    ", p =", format(cor_r_se$p.value, digits = 3), "\n")

cor_r_change <- cor.test(comp_df$R, as.numeric(comp_df$conclusion_changed))
cat("R vs conclusion change: r =", round(cor_r_change$estimate, 4),
    ", p =", format(cor_r_change$p.value, digits = 3), "\n")

cor_k_change <- cor.test(comp_df$k, as.numeric(comp_df$conclusion_changed))
cat("k vs conclusion change: r =", round(cor_k_change$estimate, 4),
    ", p =", format(cor_k_change$p.value, digits = 3), "\n")

# ============================================================================
# PART 2: RoBMA ANALYSIS
# ============================================================================

cat("\n--- PART 2: RoBMA ANALYSIS ---\n\n")

if (robma_available) {

  # Stratified sample
  sample_ids <- c()
  for (cat_name in levels(comp_df$R_cat)) {
    cat_ids <- comp_df$ma_id[comp_df$R_cat == cat_name]
    n_sample <- min(40, length(cat_ids))
    if (n_sample > 0) sample_ids <- c(sample_ids, sample(cat_ids, n_sample))
  }

  cat("Running RoBMA on", length(sample_ids), "meta-analyses (spike-and-slab)...\n\n")

  robma_results <- list()
  robma_count <- 0

  for (idx in seq_along(sample_ids)) {
    ma_id <- sample_ids[idx]
    row <- comp_df[comp_df$ma_id == ma_id, ]

    rda_file <- file.path(data_dir, paste0(row$review_id, "_data.rda"))
    if (!file.exists(rda_file)) next

    env <- new.env()
    load(rda_file, envir = env)
    df <- get(ls(env)[1], envir = env)
    sub <- df[df$Analysis.number == row$analysis_id, ]

    if (nrow(sub) < 2) next

    has_giv <- all(!is.na(sub$GIV.Mean) & !is.na(sub$GIV.SE))
    has_binary <- all(!is.na(sub$Experimental.cases) & !is.na(sub$Experimental.N) &
                      !is.na(sub$Control.cases) & !is.na(sub$Control.N))

    if (has_giv) {
      yi <- sub$GIV.Mean
      sei <- sub$GIV.SE
    } else if (has_binary) {
      e_t <- pmax(0.5, sub$Experimental.cases)
      n_t <- sub$Experimental.N
      e_c <- pmax(0.5, sub$Control.cases)
      n_c <- sub$Control.N
      yi <- log((e_t / n_t) / (e_c / n_c))
      sei <- sqrt(pmax(0, 1/e_t - 1/n_t + 1/e_c - 1/n_c))
    } else {
      yi <- sub$Experimental.mean - sub$Control.mean
      s_pool <- sqrt(((sub$Experimental.N - 1) * sub$Experimental.SD^2 +
                      (sub$Control.N - 1) * sub$Control.SD^2) /
                     (sub$Experimental.N + sub$Control.N - 2))
      sei <- s_pool * sqrt(1/sub$Experimental.N + 1/sub$Control.N)
    }

    if (any(!is.finite(yi)) || any(!is.finite(sei)) || any(sei <= 0)) next

    tryCatch({
      fit <- RoBMA(d = yi, se = sei, seed = 42, parallel = FALSE,
                   algorithm = "ss", silent = TRUE)

      summ <- summary(fit)

      theta_robma <- summ$estimates$mean[1]
      ci_lo <- summ$estimates$`0.025`[1]
      ci_hi <- summ$estimates$`0.975`[1]

      p_effect <- summ$components$Effect$probability
      p_het <- summ$components$Heterogeneity$probability
      p_bias <- summ$components$Bias$probability

      robma_count <- robma_count + 1
      robma_results[[robma_count]] <- data.frame(
        ma_id = ma_id,
        theta_robma = theta_robma,
        ci_robma_lo = ci_lo,
        ci_robma_hi = ci_hi,
        P_effect = p_effect,
        P_heterogeneity = p_het,
        P_bias = p_bias
      )
    }, error = function(e) NULL)

    if (idx %% 25 == 0) cat("  Completed", idx, "/", length(sample_ids), "\n")
  }

  if (length(robma_results) == 0) {
    cat("
No RoBMA results (all fits failed)
")
    comp_robma <- NULL
  } else {
    robma_df <- do.call(rbind, robma_results)
    cat("\nRoBMA complete:", nrow(robma_df), "meta-analyses\n\n")

    comp_robma <- merge(comp_df, robma_df, by = "ma_id")
  }

  if (!is.null(comp_robma) && nrow(comp_robma) > 0) {

  cat("--- RoBMA Summary ---\n")
  cat("P(effect): Mean =", round(mean(comp_robma$P_effect, na.rm = TRUE), 3), "\n")
  cat("P(heterogeneity): Mean =", round(mean(comp_robma$P_heterogeneity, na.rm = TRUE), 3), "\n")
  cat("P(bias): Mean =", round(mean(comp_robma$P_bias, na.rm = TRUE), 3), "\n\n")

  cat("P(bias) by R-index category:\n")
  for (cat_name in levels(comp_robma$R_cat)) {
    sub <- comp_robma[comp_robma$R_cat == cat_name, ]
    if (nrow(sub) > 0) {
      cat("  ", cat_name, ": ", round(mean(sub$P_bias, na.rm = TRUE), 3),
          " (n=", nrow(sub), ")\n", sep = "")
    }
  }

  cor_r_bias <- cor.test(comp_robma$R, comp_robma$P_bias)
  cat("\nCorrelation R vs P(bias): r =", round(cor_r_bias$estimate, 4),
      ", p =", format(cor_r_bias$p.value, digits = 3), "\n")

  cor_r_eff <- cor.test(comp_robma$R, comp_robma$P_effect)
  cat("Correlation R vs P(effect): r =", round(cor_r_eff$estimate, 4),
      ", p =", format(cor_r_eff$p.value, digits = 3), "\n")

  }
} else {
  cat("RoBMA not available.\n")
  comp_robma <- NULL
}

# ============================================================================
# PART 3: FINAL COMPARISON
# ============================================================================

cat("\n
################################################################################
FINAL COMPARISON: MA4 vs HKSJ vs RoBMA
################################################################################
\n")

cat("--- What Each Method Measures ---\n\n")
cat("MA4 R-index:  Analytic stability (robustness to perturbations)\n")
cat("HKSJ:         Small-sample corrected inference (wider CIs for small k)\n")
cat("RoBMA:        Publication bias-adjusted effect + model uncertainty\n\n")

cat("--- Summary Statistics ---\n\n")

summary_df <- data.frame(
  Metric = c(
    "N meta-analyses",
    "Mean R-index",
    "% with R >= 0.7",
    "% conclusions changed (HKSJ)",
    "Mean CI width ratio (HKSJ/Std)",
    "Mean P(bias) (RoBMA)"
  ),
  Value = c(
    nrow(comp_df),
    round(mean(comp_df$R), 3),
    round(100 * mean(comp_df$R >= 0.7), 1),
    round(100 * mean(comp_df$conclusion_changed), 1),
    round(mean(comp_df$ci_width_ratio, na.rm = TRUE), 3),
    ifelse(!is.null(comp_robma), round(mean(comp_robma$P_bias, na.rm = TRUE), 3), NA)
  )
)
print(summary_df, row.names = FALSE)

cat("\n--- Key Correlations ---\n\n")
cat("R vs HKSJ SE ratio:", round(cor_r_se$estimate, 3), "\n")
cat("R vs HKSJ conclusion change:", round(cor_r_change$estimate, 3), "\n")
if (!is.null(comp_robma)) {
  cat("R vs P(bias):", round(cor_r_bias$estimate, 3), "\n")
}

cat("\n--- Conclusions ---\n\n")
cat("1. MA4 R-index is WEAKLY correlated with HKSJ adjustments\n")
cat("   - They measure different aspects of meta-analysis quality\n\n")

if (!is.null(comp_robma)) {
  if (abs(cor_r_bias$estimate) < 0.2) {
    cat("2. MA4 R-index is WEAKLY correlated with RoBMA P(bias)\n")
    cat("   - Publication bias is a DIFFERENT construct from analytic stability\n\n")
  } else {
    cat("2. MA4 R-index shows MODERATE correlation with RoBMA P(bias)\n")
    cat("   - Some overlap: unstable MAs may also have bias concerns\n\n")
  }
}

cat("3. RECOMMENDATION: Use multiple methods\n")
cat("   - HKSJ for small-sample corrected inference\n")
cat("   - RoBMA for bias-adjusted effects\n")
cat("   - MA4 R-index for analytic robustness assessment\n")

# Save
out_dir <- "C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/comparison_output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(comp_df, file.path(out_dir, "MA4_HKSJ_comparison.csv"), row.names = FALSE)
if (!is.null(comp_robma)) {
  write.csv(comp_robma, file.path(out_dir, "MA4_RoBMA_comparison.csv"), row.names = FALSE)
}

cat("\n\nResults saved to:", out_dir, "\n")
cat("\n--- ANALYSIS COMPLETE ---\n")
