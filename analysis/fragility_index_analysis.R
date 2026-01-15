################################################################################
# FRAGILITY INDEX ANALYSIS FOR PAIRWISE70
#
# Computes fragility metrics across all Cochrane pairwise meta-analyses
#
# Fragility Types:
#   1. Direction Fragility: Does removing one study flip effect direction?
#   2. Significance Fragility: Does removing one study change p < 0.05 status?
#   3. Clinical Fragility: Does removing one study cross clinical thresholds?
#
# Outputs:
#   - fragility_analysis_results.csv: Per-analysis fragility metrics
#   - fragility_summary.csv: Summary statistics
#   - fragility_report.md: Paper-ready report
#   - Visualizations in plots/fragility/
################################################################################

library(metafor)
library(data.table)

# Configuration
args <- commandArgs(trailingOnly = TRUE)
data_dir <- if(length(args) >= 1) args[1] else "C:/Users/user/OneDrive - NHS/Documents/Pairwise70/data"
output_dir <- if(length(args) >= 2) args[2] else "C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/output"

# Clinical thresholds for OR/RR (log scale)
CLINICAL_THRESHOLD_OR <- log(1.25)  # 25% change considered clinically meaningful
CLINICAL_THRESHOLD_SMD <- 0.2       # Small effect size threshold

# Create output directories
dir.create(file.path(output_dir, "plots", "fragility"), recursive = TRUE, showWarnings = FALSE)

cat("=== FRAGILITY INDEX ANALYSIS ===\n")
cat("Data directory:", data_dir, "\n")
cat("Output directory:", output_dir, "\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

################################################################################
# HELPER FUNCTIONS
################################################################################

# Determine outcome type and run appropriate meta-analysis
run_meta_analysis <- function(d, measure = NULL) {
  # Check for binary data
  has_binary <- !all(is.na(d$Experimental.cases)) && !all(is.na(d$Control.cases)) &&
                !all(is.na(d$Experimental.N)) && !all(is.na(d$Control.N))

  # Check for continuous data
  has_continuous <- !all(is.na(d$Experimental.mean)) && !all(is.na(d$Control.mean)) &&
                    !all(is.na(d$Experimental.SD)) && !all(is.na(d$Control.SD))

  # Check for generic inverse variance
  has_giv <- !all(is.na(d$GIV.Mean)) && !all(is.na(d$GIV.SE))

  result <- NULL
  used_measure <- NA

  if(has_binary) {
    # Filter valid rows
    valid <- complete.cases(d[, c("Experimental.cases", "Experimental.N",
                                   "Control.cases", "Control.N")])
    if(sum(valid) < 2) return(list(result = NULL, measure = NA))

    dd <- d[valid, ]

    # Check for zero cells
    has_zero <- any(dd$Experimental.cases == 0 | dd$Control.cases == 0 |
                    dd$Experimental.cases == dd$Experimental.N |
                    dd$Control.cases == dd$Control.N)

    tryCatch({
      if(has_zero) {
        # Use continuity correction
        result <- rma(measure = "OR",
                      ai = Experimental.cases, n1i = Experimental.N,
                      ci = Control.cases, n2i = Control.N,
                      data = dd, method = "REML", add = 0.5, to = "only0")
      } else {
        result <- rma(measure = "OR",
                      ai = Experimental.cases, n1i = Experimental.N,
                      ci = Control.cases, n2i = Control.N,
                      data = dd, method = "REML")
      }
      used_measure <- "OR"
    }, error = function(e) NULL)

  } else if(has_continuous) {
    valid <- complete.cases(d[, c("Experimental.mean", "Experimental.SD", "Experimental.N",
                                   "Control.mean", "Control.SD", "Control.N")])
    if(sum(valid) < 2) return(list(result = NULL, measure = NA))

    dd <- d[valid, ]
    dd <- dd[dd$Experimental.SD > 0 & dd$Control.SD > 0, ]
    if(nrow(dd) < 2) return(list(result = NULL, measure = NA))

    tryCatch({
      result <- rma(measure = "SMD",
                    m1i = Experimental.mean, sd1i = Experimental.SD, n1i = Experimental.N,
                    m2i = Control.mean, sd2i = Control.SD, n2i = Control.N,
                    data = dd, method = "REML")
      used_measure <- "SMD"
    }, error = function(e) NULL)

  } else if(has_giv) {
    valid <- complete.cases(d[, c("GIV.Mean", "GIV.SE")])
    if(sum(valid) < 2) return(list(result = NULL, measure = NA))

    dd <- d[valid, ]
    dd <- dd[dd$GIV.SE > 0, ]
    if(nrow(dd) < 2) return(list(result = NULL, measure = NA))

    tryCatch({
      result <- rma(yi = GIV.Mean, sei = GIV.SE, data = dd, method = "REML")
      used_measure <- "GIV"
    }, error = function(e) NULL)
  }

  return(list(result = result, measure = used_measure))
}

# Compute fragility metrics for a single meta-analysis
compute_fragility <- function(d, dataset_name, analysis_id) {

  # Run full meta-analysis
  full_result <- run_meta_analysis(d)

  if(is.null(full_result$result)) {
    return(NULL)
  }

  rma_obj <- full_result$result
  measure <- full_result$measure
  k <- rma_obj$k

  if(k < 3) {
    return(NULL)  # Need at least 3 studies for meaningful LOO
  }

  # Full analysis results
  full_est <- as.numeric(rma_obj$beta)
  full_se <- rma_obj$se
  full_pval <- rma_obj$pval
  full_ci_lb <- rma_obj$ci.lb
  full_ci_ub <- rma_obj$ci.ub
  full_tau2 <- rma_obj$tau2
  full_I2 <- rma_obj$I2

  # Determine significance and direction
  full_sig <- full_pval < 0.05
  full_direction <- sign(full_est)

  # Clinical threshold based on measure
  if(measure %in% c("OR", "RR", "GIV")) {
    clinical_threshold <- CLINICAL_THRESHOLD_OR
  } else {
    clinical_threshold <- CLINICAL_THRESHOLD_SMD
  }

  # Full analysis clinical significance
  full_clinical_pos <- full_est > clinical_threshold
  full_clinical_neg <- full_est < -clinical_threshold

  # Leave-one-out analysis
  loo_results <- tryCatch({
    leave1out(rma_obj)
  }, error = function(e) NULL)

  if(is.null(loo_results)) {
    return(NULL)
  }

  # Extract LOO estimates
  loo_est <- loo_results$estimate
  loo_pval <- loo_results$pval
  loo_ci_lb <- loo_results$ci.lb
  loo_ci_ub <- loo_results$ci.ub

  # Fragility calculations

  # 1. Direction Fragility: Does any LOO flip the sign?
  loo_directions <- sign(loo_est)
  direction_flips <- sum(loo_directions != full_direction & loo_directions != 0)
  direction_fragile <- direction_flips > 0
  direction_fragility_index <- direction_flips  # How many studies cause flip

  # 2. Significance Fragility: Does any LOO change significance status?
  loo_sig <- loo_pval < 0.05

  if(full_sig) {
    # Originally significant - does any LOO make it non-significant?
    sig_losses <- sum(!loo_sig)
    sig_fragile <- sig_losses > 0
    sig_fragility_index <- sig_losses
    sig_fragility_type <- "loss"
  } else {
    # Originally non-significant - does any LOO make it significant?
    sig_gains <- sum(loo_sig)
    sig_fragile <- sig_gains > 0
    sig_fragility_index <- sig_gains
    sig_fragility_type <- "gain"
  }

  # 3. Clinical Fragility: Does any LOO cross clinical threshold?
  loo_clinical_pos <- loo_est > clinical_threshold
  loo_clinical_neg <- loo_est < -clinical_threshold

  clinical_crossings <- 0
  if(full_clinical_pos) {
    clinical_crossings <- sum(!loo_clinical_pos)
  } else if(full_clinical_neg) {
    clinical_crossings <- sum(!loo_clinical_neg)
  } else {
    # Full was in null zone - count crossings either way
    clinical_crossings <- sum(loo_clinical_pos | loo_clinical_neg)
  }
  clinical_fragile <- clinical_crossings > 0

  # 4. Effect size change metrics
  effect_changes <- abs(loo_est - full_est)
  max_effect_change <- max(effect_changes)
  mean_effect_change <- mean(effect_changes)

  # Which study removal causes maximum change?
  max_change_idx <- which.max(effect_changes)

  # 5. Confidence interval stability
  ci_width_full <- full_ci_ub - full_ci_lb
  ci_widths_loo <- loo_ci_ub - loo_ci_lb
  ci_width_change <- mean(abs(ci_widths_loo - ci_width_full))

  # 6. Composite Fragility Index (0-3 scale)
  composite_fragility <- as.integer(direction_fragile) +
                         as.integer(sig_fragile) +
                         as.integer(clinical_fragile)

  # 7. Fragility Quotient: proportion of studies that cause any fragility
  studies_causing_direction_flip <- sum(loo_directions != full_direction & loo_directions != 0)
  studies_causing_sig_change <- if(full_sig) sum(!loo_sig) else sum(loo_sig)
  studies_causing_any_fragility <- max(studies_causing_direction_flip, studies_causing_sig_change)
  fragility_quotient <- studies_causing_any_fragility / k

  # Return results
  data.table(
    dataset = dataset_name,
    analysis_id = analysis_id,
    measure = measure,
    k = k,

    # Full analysis
    estimate = full_est,
    se = full_se,
    pval = full_pval,
    ci_lb = full_ci_lb,
    ci_ub = full_ci_ub,
    tau2 = full_tau2,
    I2 = full_I2,
    significant = full_sig,

    # Direction fragility
    direction_fragile = direction_fragile,
    direction_fragility_index = direction_fragility_index,

    # Significance fragility
    sig_fragile = sig_fragile,
    sig_fragility_index = sig_fragility_index,
    sig_fragility_type = sig_fragility_type,

    # Clinical fragility
    clinical_fragile = clinical_fragile,
    clinical_fragility_index = clinical_crossings,

    # Effect change metrics
    max_effect_change = max_effect_change,
    mean_effect_change = mean_effect_change,
    max_change_study_idx = max_change_idx,

    # CI stability
    ci_width = ci_width_full,
    ci_width_change = ci_width_change,

    # Composite metrics
    composite_fragility = composite_fragility,
    fragility_quotient = fragility_quotient,

    # Classification
    fragility_class = c("Robust", "Low Fragility", "Moderate Fragility", "High Fragility")[composite_fragility + 1]
  )
}

################################################################################
# MAIN ANALYSIS
################################################################################

# Load all datasets
files <- list.files(data_dir, pattern = "\\.rda$", full.names = TRUE)
cat("Found", length(files), "datasets\n\n")
flush.console()

# Results storage
all_results <- list()
analysis_counter <- 0
success_counter <- 0

# Progress tracking - more frequent for monitoring
pb_interval <- 10

for(i in seq_along(files)) {
  if(i == 1 || i %% pb_interval == 0) {
    cat(sprintf("Processing %d/%d (%.0f%%) - %d successes\n", i, length(files), 100*i/length(files), success_counter))
    flush.console()
  }

  # Load dataset
  env <- new.env()
  load(files[i], envir = env)
  obj_name <- ls(envir = env)[1]
  d <- get(obj_name, envir = env)
  dataset_name <- gsub("\\.rda$", "", basename(files[i]))

  # Get unique analyses
  if("Analysis.number" %in% names(d)) {
    analyses <- unique(d$Analysis.number)
  } else {
    analyses <- 1
  }

  for(analysis_num in analyses) {
    analysis_counter <- analysis_counter + 1

    # Subset to this analysis
    if("Analysis.number" %in% names(d)) {
      subset_d <- d[d$Analysis.number == analysis_num, ]
    } else {
      subset_d <- d
    }

    # Skip if too few rows
    if(nrow(subset_d) < 3) next

    # Compute fragility
    result <- tryCatch({
      compute_fragility(subset_d, dataset_name, analysis_num)
    }, error = function(e) NULL)

    if(!is.null(result)) {
      all_results[[length(all_results) + 1]] <- result
      success_counter <- success_counter + 1
    }
  }
}

cat("\n\nAnalyses attempted:", analysis_counter, "\n")
cat("Analyses with fragility computed:", success_counter, "\n")

# Combine results
results_dt <- rbindlist(all_results)

cat("\nTotal analyses in results:", nrow(results_dt), "\n")

################################################################################
# SAVE DETAILED RESULTS
################################################################################

fwrite(results_dt, file.path(output_dir, "fragility_analysis_results.csv"))
cat("Saved: fragility_analysis_results.csv\n")

################################################################################
# SUMMARY STATISTICS
################################################################################

cat("\n=== FRAGILITY SUMMARY ===\n\n")

# Overall fragility distribution
cat("Fragility Class Distribution:\n")
print(table(results_dt$fragility_class))
cat("\n")

# Fragility by significance status
cat("Fragility by Original Significance:\n")
sig_frag <- results_dt[, .(
  n = .N,
  pct_direction_fragile = 100 * mean(direction_fragile),
  pct_sig_fragile = 100 * mean(sig_fragile),
  pct_clinical_fragile = 100 * mean(clinical_fragile),
  mean_composite = mean(composite_fragility),
  mean_fragility_quotient = mean(fragility_quotient)
), by = significant]
print(sig_frag)
cat("\n")

# Fragility by k (number of studies)
cat("Fragility by Number of Studies:\n")
results_dt[, k_group := cut(k, breaks = c(2, 5, 10, 20, 50, Inf),
                            labels = c("3-5", "6-10", "11-20", "21-50", "51+"))]
k_frag <- results_dt[, .(
  n = .N,
  pct_direction_fragile = 100 * mean(direction_fragile),
  pct_sig_fragile = 100 * mean(sig_fragile),
  pct_clinical_fragile = 100 * mean(clinical_fragile),
  mean_composite = mean(composite_fragility)
), by = k_group]
print(k_frag[order(k_group)])
cat("\n")

# Fragility by I2
cat("Fragility by Heterogeneity (I2):\n")
results_dt[, I2_group := cut(I2, breaks = c(-Inf, 25, 50, 75, Inf),
                             labels = c("Low (<25%)", "Moderate (25-50%)",
                                        "Substantial (50-75%)", "High (>75%)"))]
I2_frag <- results_dt[, .(
  n = .N,
  pct_direction_fragile = 100 * mean(direction_fragile),
  pct_sig_fragile = 100 * mean(sig_fragile),
  mean_composite = mean(composite_fragility)
), by = I2_group]
print(I2_frag[order(I2_group)])
cat("\n")

# Key statistics for paper
cat("\n=== KEY STATISTICS FOR PAPER ===\n\n")

n_total <- nrow(results_dt)
n_sig <- sum(results_dt$significant)
n_nonsig <- n_total - n_sig

cat(sprintf("Total meta-analyses analyzed: %d\n", n_total))
cat(sprintf("  - Statistically significant (p<0.05): %d (%.1f%%)\n",
            n_sig, 100*n_sig/n_total))
cat(sprintf("  - Not significant: %d (%.1f%%)\n",
            n_nonsig, 100*n_nonsig/n_total))

cat("\nDirection Fragility:\n")
n_dir_fragile <- sum(results_dt$direction_fragile)
cat(sprintf("  - Analyses where removing 1 study flips direction: %d (%.1f%%)\n",
            n_dir_fragile, 100*n_dir_fragile/n_total))

cat("\nSignificance Fragility:\n")
n_sig_fragile <- sum(results_dt$sig_fragile)
cat(sprintf("  - Analyses where removing 1 study changes significance: %d (%.1f%%)\n",
            n_sig_fragile, 100*n_sig_fragile/n_total))

# Among significant results
sig_results <- results_dt[significant == TRUE]
n_sig_lose <- sum(sig_results$sig_fragile)
cat(sprintf("  - Among significant results, fragile to non-significance: %d/%d (%.1f%%)\n",
            n_sig_lose, nrow(sig_results), 100*n_sig_lose/nrow(sig_results)))

cat("\nClinical Fragility:\n")
n_clin_fragile <- sum(results_dt$clinical_fragile)
cat(sprintf("  - Analyses where removing 1 study crosses clinical threshold: %d (%.1f%%)\n",
            n_clin_fragile, 100*n_clin_fragile/n_total))

cat("\nComposite Fragility:\n")
for(class in c("Robust", "Low Fragility", "Moderate Fragility", "High Fragility")) {
  n_class <- sum(results_dt$fragility_class == class)
  cat(sprintf("  - %s: %d (%.1f%%)\n", class, n_class, 100*n_class/n_total))
}

cat("\nFragility Quotient (proportion of studies causing fragility):\n")
cat(sprintf("  - Mean: %.3f\n", mean(results_dt$fragility_quotient)))
cat(sprintf("  - Median: %.3f\n", median(results_dt$fragility_quotient)))
cat(sprintf("  - 75th percentile: %.3f\n", quantile(results_dt$fragility_quotient, 0.75)))

################################################################################
# SAVE SUMMARY
################################################################################

summary_dt <- data.table(
  metric = c(
    "Total analyses",
    "Significant (p<0.05)",
    "Not significant",
    "Direction fragile",
    "Significance fragile",
    "Clinical fragile",
    "Robust (composite=0)",
    "Low fragility (composite=1)",
    "Moderate fragility (composite=2)",
    "High fragility (composite=3)",
    "Mean fragility quotient",
    "Median k (studies per analysis)"
  ),
  value = c(
    n_total,
    n_sig,
    n_nonsig,
    n_dir_fragile,
    n_sig_fragile,
    n_clin_fragile,
    sum(results_dt$fragility_class == "Robust"),
    sum(results_dt$fragility_class == "Low Fragility"),
    sum(results_dt$fragility_class == "Moderate Fragility"),
    sum(results_dt$fragility_class == "High Fragility"),
    round(mean(results_dt$fragility_quotient), 4),
    median(results_dt$k)
  ),
  percentage = c(
    NA,
    round(100*n_sig/n_total, 1),
    round(100*n_nonsig/n_total, 1),
    round(100*n_dir_fragile/n_total, 1),
    round(100*n_sig_fragile/n_total, 1),
    round(100*n_clin_fragile/n_total, 1),
    round(100*sum(results_dt$fragility_class == "Robust")/n_total, 1),
    round(100*sum(results_dt$fragility_class == "Low Fragility")/n_total, 1),
    round(100*sum(results_dt$fragility_class == "Moderate Fragility")/n_total, 1),
    round(100*sum(results_dt$fragility_class == "High Fragility")/n_total, 1),
    NA,
    NA
  )
)

fwrite(summary_dt, file.path(output_dir, "fragility_summary.csv"))
cat("\nSaved: fragility_summary.csv\n")

################################################################################
# GENERATE VISUALIZATIONS
################################################################################

cat("\nGenerating visualizations...\n")

# 1. Fragility class distribution
png(file.path(output_dir, "plots", "fragility", "fragility_class_distribution.png"),
    width = 800, height = 600, res = 120)
par(mar = c(5, 4, 4, 2) + 0.1)
class_counts <- table(factor(results_dt$fragility_class,
                             levels = c("Robust", "Low Fragility",
                                        "Moderate Fragility", "High Fragility")))
barplot(class_counts,
        main = "Fragility Classification of Cochrane Meta-Analyses",
        ylab = "Number of Analyses",
        col = c("#2ecc71", "#f1c40f", "#e67e22", "#e74c3c"),
        las = 1)
dev.off()

# 2. Fragility by k
png(file.path(output_dir, "plots", "fragility", "fragility_by_k.png"),
    width = 900, height = 600, res = 120)
par(mar = c(5, 4, 4, 2) + 0.1)
boxplot(composite_fragility ~ k_group, data = results_dt,
        main = "Composite Fragility by Number of Studies",
        xlab = "Number of Studies",
        ylab = "Composite Fragility Index (0-3)",
        col = "#3498db")
dev.off()

# 3. Fragility quotient distribution
png(file.path(output_dir, "plots", "fragility", "fragility_quotient_histogram.png"),
    width = 800, height = 600, res = 120)
par(mar = c(5, 4, 4, 2) + 0.1)
hist(results_dt$fragility_quotient,
     breaks = 30,
     main = "Distribution of Fragility Quotient",
     xlab = "Fragility Quotient (proportion of studies causing fragility)",
     ylab = "Frequency",
     col = "#9b59b6",
     border = "white")
abline(v = mean(results_dt$fragility_quotient), col = "red", lwd = 2, lty = 2)
legend("topright", legend = sprintf("Mean = %.3f", mean(results_dt$fragility_quotient)),
       col = "red", lty = 2, lwd = 2)
dev.off()

# 4. Significance fragility by original significance
png(file.path(output_dir, "plots", "fragility", "sig_fragility_by_original.png"),
    width = 800, height = 600, res = 120)
par(mar = c(5, 6, 4, 2) + 0.1)
sig_table <- table(results_dt$significant, results_dt$sig_fragile)
barplot(t(sig_table), beside = TRUE,
        names.arg = c("Not Significant", "Significant"),
        legend.text = c("Not Fragile", "Fragile"),
        main = "Significance Fragility by Original Significance Status",
        ylab = "Number of Analyses",
        col = c("#2ecc71", "#e74c3c"),
        args.legend = list(x = "topright"))
dev.off()

# 5. Fragility vs I2
png(file.path(output_dir, "plots", "fragility", "fragility_vs_I2.png"),
    width = 800, height = 600, res = 120)
par(mar = c(5, 4, 4, 2) + 0.1)
plot(results_dt$I2, results_dt$composite_fragility,
     pch = 16, col = rgb(0.2, 0.4, 0.6, 0.3),
     main = "Composite Fragility vs. Heterogeneity",
     xlab = "I-squared (%)",
     ylab = "Composite Fragility Index (0-3)")
abline(lm(composite_fragility ~ I2, data = results_dt), col = "red", lwd = 2)
dev.off()

cat("Saved 5 visualizations to plots/fragility/\n")

################################################################################
# GENERATE PAPER-READY REPORT
################################################################################

report <- sprintf('# Fragility of Cochrane Meta-Analyses: An Empirical Investigation

Generated: %s

## Abstract

This analysis examines the fragility of %d Cochrane pairwise meta-analyses from 501 systematic reviews. We define fragility as the vulnerability of meta-analysis conclusions to the removal of a single study.

## Key Findings

### Prevalence of Fragility

| Fragility Type | N | Percentage |
|----------------|---|------------|
| Direction fragile (effect sign flips) | %d | %.1f%% |
| Significance fragile (p-value crosses 0.05) | %d | %.1f%% |
| Clinical fragile (crosses threshold) | %d | %.1f%% |

### Composite Fragility Classification

| Class | N | Percentage |
|-------|---|------------|
| Robust (no fragility) | %d | %.1f%% |
| Low Fragility (1 type) | %d | %.1f%% |
| Moderate Fragility (2 types) | %d | %.1f%% |
| High Fragility (3 types) | %d | %.1f%% |

### Significance Fragility Among Significant Results

Of %d meta-analyses with p < 0.05:
- **%d (%.1f%%)** would lose statistical significance if a single study were removed

This represents a substantial proportion of Cochrane meta-analyses whose conclusions depend critically on individual studies.

### Fragility by Number of Studies

| Studies | N | Direction Fragile | Significance Fragile |
|---------|---|-------------------|---------------------|
%s

### Fragility by Heterogeneity

| I-squared | N | Direction Fragile | Significance Fragile |
|-----------|---|-------------------|---------------------|
%s

## Implications

1. **For systematic reviewers**: Report sensitivity of conclusions to individual study removal
2. **For guideline developers**: Consider fragility when translating evidence to recommendations
3. **For methodologists**: Develop robust pooling methods that minimize fragility

## Methods

- **Data source**: 501 Cochrane systematic reviews (Pairwise70 dataset)
- **Analyses included**: Meta-analyses with k >= 3 studies
- **Effect measures**: OR for binary, SMD for continuous outcomes
- **Fragility metrics**:
  - Direction fragility: Sign of effect estimate reverses
  - Significance fragility: p-value crosses 0.05 threshold
  - Clinical fragility: Effect crosses clinical threshold (OR=1.25, SMD=0.2)
  - Fragility quotient: Proportion of studies whose removal causes any fragility

## Visualizations

See `plots/fragility/` for:
- Fragility class distribution
- Fragility by number of studies
- Fragility quotient histogram
- Significance fragility breakdown
- Fragility vs heterogeneity

## Data Availability

Full results: `fragility_analysis_results.csv`
Summary: `fragility_summary.csv`

',
format(Sys.time(), "%%Y-%%m-%%d %%H:%%M:%%S"),
n_total,
n_dir_fragile, 100*n_dir_fragile/n_total,
n_sig_fragile, 100*n_sig_fragile/n_total,
n_clin_fragile, 100*n_clin_fragile/n_total,
sum(results_dt$fragility_class == "Robust"),
100*sum(results_dt$fragility_class == "Robust")/n_total,
sum(results_dt$fragility_class == "Low Fragility"),
100*sum(results_dt$fragility_class == "Low Fragility")/n_total,
sum(results_dt$fragility_class == "Moderate Fragility"),
100*sum(results_dt$fragility_class == "Moderate Fragility")/n_total,
sum(results_dt$fragility_class == "High Fragility"),
100*sum(results_dt$fragility_class == "High Fragility")/n_total,
nrow(sig_results),
n_sig_lose, 100*n_sig_lose/nrow(sig_results),
paste(apply(k_frag[order(k_group)], 1, function(r) {
  sprintf("| %s | %d | %.1f%% | %.1f%% |", r["k_group"], as.integer(r["n"]),
          as.numeric(r["pct_direction_fragile"]), as.numeric(r["pct_sig_fragile"]))
}), collapse = "\n"),
paste(apply(I2_frag[order(I2_group)], 1, function(r) {
  sprintf("| %s | %d | %.1f%% | %.1f%% |", r["I2_group"], as.integer(r["n"]),
          as.numeric(r["pct_direction_fragile"]), as.numeric(r["pct_sig_fragile"]))
}), collapse = "\n")
)

writeLines(report, file.path(output_dir, "fragility_report.md"))
cat("\nSaved: fragility_report.md\n")

################################################################################
# COMPLETION
################################################################################

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\nOutput files:\n")
cat("  - fragility_analysis_results.csv (detailed per-analysis results)\n")
cat("  - fragility_summary.csv (summary statistics)\n")
cat("  - fragility_report.md (paper-ready report)\n")
cat("  - plots/fragility/*.png (5 visualizations)\n")
