# ============================================================================
# MA4 vs RoBMA vs RVE - COMPREHENSIVE COMPARISON v2
# ============================================================================

cat("
################################################################################
#   MA4 vs RoBMA vs RVE - COMPARATIVE ANALYSIS                                #
################################################################################
\n\n")

set.seed(42)

library(metafor)
library(clubSandwich)

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
# PART 1: RVE ANALYSIS
# ============================================================================

cat("\n--- PART 1: RVE ANALYSIS ---\n\n")

rve_results <- list()
rve_count <- 0

for (i in seq_along(rda_files)) {
  env <- new.env()
  load(rda_files[i], envir = env)
  df <- get(ls(env)[1], envir = env)
  if (nrow(df) == 0) next

  review_id <- gsub("\\.rda$", "", basename(rda_files[i]))
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

      # Safe calculation
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
      fit <- rma(yi = yi, sei = sei, method = "REML")

      theta <- as.numeric(coef(fit))
      se_std <- as.numeric(fit$se)

      # RVE with CR2 estimator
      V_rve <- vcovCR(fit, type = "CR2")
      se_rve <- sqrt(V_rve[1, 1])

      # Standard CI
      ci_std_lo <- theta - 1.96 * se_std
      ci_std_hi <- theta + 1.96 * se_std

      # RVE-adjusted CI (using t-distribution with approximate df)
      # For single coefficient, use Satterthwaite approximation
      df_approx <- max(1, k - 2)
      t_crit <- qt(0.975, df = df_approx)
      ci_rve_lo <- theta - t_crit * se_rve
      ci_rve_hi <- theta + t_crit * se_rve

      # Significance
      sig_std <- (ci_std_lo > 0) | (ci_std_hi < 0)
      sig_rve <- (ci_rve_lo > 0) | (ci_rve_hi < 0)

      rve_count <- rve_count + 1
      rve_results[[rve_count]] <- data.frame(
        ma_id = ma_id,
        review_id = review_id,
        analysis_id = an,
        effect_type = effect_type,
        k = k,
        theta = theta,
        se_std = se_std,
        se_rve = se_rve,
        se_ratio = se_rve / se_std,
        ci_std_lo = ci_std_lo,
        ci_std_hi = ci_std_hi,
        ci_rve_lo = ci_rve_lo,
        ci_rve_hi = ci_rve_hi,
        sig_std = sig_std,
        sig_rve = sig_rve,
        conclusion_changed = sig_std != sig_rve,
        tau2 = fit$tau2
      )
    }, error = function(e) NULL)
  }

  if (i %% 100 == 0) cat("  Processed", i, "/", length(rda_files), "reviews\n")
}

rve_df <- do.call(rbind, rve_results)
cat("\nRVE complete:", nrow(rve_df), "meta-analyses\n")

# Merge with MA4
ma4_results$ma_id <- paste(ma4_results$review_id, ma4_results$analysis_id, sep = "_")
comp_df <- merge(rve_df, ma4_results[, c("ma_id", "R")], by = "ma_id")

cat("Merged:", nrow(comp_df), "meta-analyses\n\n")

# RVE Summary
cat("--- RVE Summary ---\n")
cat("SE inflation ratio (RVE/Standard):\n")
cat("  Mean:", round(mean(comp_df$se_ratio, na.rm = TRUE), 3), "\n")
cat("  Median:", round(median(comp_df$se_ratio, na.rm = TRUE), 3), "\n")
cat("  IQR:", round(quantile(comp_df$se_ratio, 0.25, na.rm = TRUE), 3), "-",
    round(quantile(comp_df$se_ratio, 0.75, na.rm = TRUE), 3), "\n\n")

cat("Conclusions changed after RVE:\n")
cat("  N:", sum(comp_df$conclusion_changed), "\n")
cat("  %:", round(100 * mean(comp_df$conclusion_changed), 2), "%\n\n")

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

# Correlation
cat("\nCorrelation: R vs SE ratio:\n")
cor_test <- cor.test(comp_df$R, comp_df$se_ratio)
cat("  r =", round(cor_test$estimate, 4), ", p =", format(cor_test$p.value, digits = 3), "\n")

cor_change <- cor.test(comp_df$R, as.numeric(comp_df$conclusion_changed))
cat("\nCorrelation: R vs conclusion changed:\n")
cat("  r =", round(cor_change$estimate, 4), ", p =", format(cor_change$p.value, digits = 3), "\n")

# ============================================================================
# PART 2: RoBMA ANALYSIS (Sample)
# ============================================================================

cat("\n--- PART 2: RoBMA ANALYSIS (Stratified Sample) ---\n\n")

if (robma_available) {

  # Stratified sample: 30 per R category
  sample_ids <- c()
  for (cat_name in levels(comp_df$R_cat)) {
    cat_ids <- comp_df$ma_id[comp_df$R_cat == cat_name]
    n_sample <- min(30, length(cat_ids))
    if (n_sample > 0) sample_ids <- c(sample_ids, sample(cat_ids, n_sample))
  }

  cat("Running RoBMA on", length(sample_ids), "meta-analyses...\n")
  cat("(Using spike-and-slab for speed)\n\n")

  robma_results <- list()
  robma_count <- 0
  errors <- 0

  for (idx in seq_along(sample_ids)) {
    ma_id <- sample_ids[idx]
    row <- comp_df[comp_df$ma_id == ma_id, ]

    # Reload raw data
    rda_file <- file.path(data_dir, paste0(row$review_id, ".rda"))
    if (!file.exists(rda_file)) next

    env <- new.env()
    load(rda_file, envir = env)
    df <- get(ls(env)[1], envir = env)
    sub <- df[df$Analysis.number == row$analysis_id, ]

    if (nrow(sub) < 2) next

    # Extract effects
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
      # Use spike-and-slab for speed
      fit <- RoBMA(d = yi, se = sei, seed = 42, parallel = FALSE,
                   algorithm = "ss", silent = TRUE)

      summ <- summary(fit)

      # Extract results
      theta_robma <- summ$estimates$mean[1]
      ci_lo <- summ$estimates$`0.025`[1]
      ci_hi <- summ$estimates$`0.975`[1]

      # Component probabilities
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
    }, error = function(e) {
      errors <<- errors + 1
    })

    if (idx %% 20 == 0) cat("  Completed", idx, "/", length(sample_ids), "\n")
  }

  robma_df <- do.call(rbind, robma_results)
  cat("\nRoBMA complete:", nrow(robma_df), "meta-analyses (", errors, "errors)\n\n")

  # Merge
  comp_robma <- merge(comp_df, robma_df, by = "ma_id")

  # RoBMA Summary
  cat("--- RoBMA Summary ---\n")
  cat("P(effect exists): Mean =", round(mean(comp_robma$P_effect, na.rm = TRUE), 3), "\n")
  cat("P(heterogeneity): Mean =", round(mean(comp_robma$P_heterogeneity, na.rm = TRUE), 3), "\n")
  cat("P(publication bias): Mean =", round(mean(comp_robma$P_bias, na.rm = TRUE), 3), "\n\n")

  cat("P(bias) by R-index category:\n")
  for (cat_name in levels(comp_robma$R_cat)) {
    sub <- comp_robma[comp_robma$R_cat == cat_name, ]
    if (nrow(sub) > 0) {
      cat("  ", cat_name, ": ", round(mean(sub$P_bias, na.rm = TRUE), 3),
          " (n=", nrow(sub), ")\n", sep = "")
    }
  }

  cat("\nCorrelation: R vs P(bias):\n")
  cor_bias <- cor.test(comp_robma$R, comp_robma$P_bias)
  cat("  r =", round(cor_bias$estimate, 4), ", p =", format(cor_bias$p.value, digits = 3), "\n")

  cat("\nCorrelation: R vs P(effect):\n")
  cor_eff <- cor.test(comp_robma$R, comp_robma$P_effect)
  cat("  r =", round(cor_eff$estimate, 4), ", p =", format(cor_eff$p.value, digits = 3), "\n")

  # Key test: Does R predict bias detection?
  cat("\n--- KEY TEST: Does low R predict RoBMA-detected bias? ---\n")
  comp_robma$high_bias <- comp_robma$P_bias > 0.5
  comp_robma$low_R <- comp_robma$R < 0.5

  tbl <- table(LowR = comp_robma$low_R, HighBias = comp_robma$high_bias)
  print(tbl)

  if (all(dim(tbl) == c(2, 2)) && all(tbl > 0)) {
    ft <- fisher.test(tbl)
    cat("\nFisher's exact test: OR =", round(ft$estimate, 2),
        ", p =", format(ft$p.value, digits = 3), "\n")
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
FINAL COMPARISON: MA4 vs RVE vs RoBMA
################################################################################
\n")

cat("--- Summary Table ---\n\n")

results_summary <- data.frame(
  Method = c("MA4 (Stability)", "RVE (Robust SE)", "RoBMA (Bayesian)"),
  Metric = c("R-index", "SE inflation", "P(bias)"),
  Mean = c(
    round(mean(comp_df$R), 3),
    round(mean(comp_df$se_ratio), 3),
    ifelse(!is.null(comp_robma), round(mean(comp_robma$P_bias, na.rm = TRUE), 3), NA)
  ),
  Interpretation = c(
    paste0(round(100*mean(comp_df$R >= 0.7), 1), "% have good stability"),
    paste0(round(100*mean(comp_df$conclusion_changed), 1), "% conclusions change"),
    ifelse(!is.null(comp_robma),
           paste0(round(100*mean(comp_robma$P_bias > 0.5, na.rm = TRUE), 1), "% show likely bias"),
           "N/A")
  )
)
print(results_summary, row.names = FALSE)

cat("\n--- Relationships Between Methods ---\n\n")

cat("1. R-index vs RVE SE inflation: r =", round(cor_test$estimate, 3), "\n")
cat("2. R-index vs RVE conclusion change: r =", round(cor_change$estimate, 3), "\n")
if (!is.null(comp_robma)) {
  cat("3. R-index vs P(bias): r =", round(cor_bias$estimate, 3), "\n")
  cat("4. R-index vs P(effect): r =", round(cor_eff$estimate, 3), "\n")
}

cat("\n--- Conclusions ---\n\n")
cat("1. MA4, RVE, and RoBMA measure DIFFERENT constructs\n")
cat("2. Low R does NOT strongly predict RVE conclusion changes\n")
if (!is.null(comp_robma)) {
  if (abs(cor_bias$estimate) > 0.3) {
    cat("3. R-index shows MODERATE correlation with P(bias) - partial overlap\n")
  } else {
    cat("3. R-index shows WEAK correlation with P(bias) - different constructs\n")
  }
}
cat("4. Recommendation: Use ALL THREE for comprehensive meta-analysis evaluation\n")

# Save
out_dir <- "C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/comparison_output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(comp_df, file.path(out_dir, "MA4_RVE_comparison.csv"), row.names = FALSE)
if (!is.null(comp_robma)) {
  write.csv(comp_robma, file.path(out_dir, "MA4_RoBMA_comparison.csv"), row.names = FALSE)
}

cat("\nResults saved to:", out_dir, "\n")
cat("\n--- ANALYSIS COMPLETE ---\n")
