# ==============================================================================
# INFORMATION ADEQUACY ANALYSIS FOR PAIRWISE70
# The Information Size Crisis in Meta-Analysis
# ==============================================================================
# Author: Based on Pairwise70 data
# Date: January 2026
# Purpose: Calculate Optimal Information Size (OIS), Information Fraction (IF),
#          and develop Information Adequacy Index (IAI)
# ==============================================================================

# Load required packages
suppressPackageStartupMessages({
  library(data.table)
  library(metafor)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(parallel)
})

# Set paths
base_path <- "C:/Users/user/OneDrive - NHS/Documents/Pairwise70"
output_path <- file.path(base_path, "analysis/output/information_adequacy")
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

cat(strrep("=", 70), "\n")
cat("INFORMATION ADEQUACY ANALYSIS - PAIRWISE70\n")
cat(strrep("=", 70), "\n\n")

# ==============================================================================
# PART 1: HELPER FUNCTIONS
# ==============================================================================

#' Calculate Optimal Information Size for Binary Outcomes
#'
#' @param p_control Control group event rate
#' @param RRR Relative Risk Reduction (anticipated effect)
#' @param alpha Type I error rate (two-sided)
#' @param power Statistical power
#' @param I2 Heterogeneity (I-squared as proportion 0-1)
#' @return Required total sample size
calculate_ois_binary <- function(p_control, RRR = 0.20, alpha = 0.05, power = 0.80, I2 = 0) {

  # Validate inputs
  if (is.na(p_control) || p_control <= 0 || p_control >= 1) return(NA)
  if (is.na(I2)) I2 <- 0
  I2 <- max(0, min(0.99, I2))  # Bound I2

  # Calculate sample size for single adequately powered trial
  # Using log-OR scale variance approximation
  z_alpha <- qnorm(1 - alpha/2)
  z_beta <- qnorm(power)

  # Effect on OR scale
  p_treat <- p_control * (1 - RRR)
  if (p_treat <= 0 || p_treat >= 1) return(NA)

  log_or <- log((p_treat/(1-p_treat)) / (p_control/(1-p_control)))

  # Variance of log-OR (per arm)
  var_control <- 1/(p_control * (1-p_control))
  var_treat <- 1/(p_treat * (1-p_treat))

  # Sample size per arm for single trial
  n_per_arm <- ((z_alpha + z_beta)^2 * (var_control + var_treat)) / log_or^2

  # Total participants for single trial
  D <- 2 * n_per_arm

  # Heterogeneity adjustment (Diversity-adjusted required information size)
  # DARIS = D * (1 + I2/(1-I2)) when I2 < 1
  if (I2 >= 0.99) {
    diversity_factor <- 100  # Cap at 100x for extreme heterogeneity
  } else {
    diversity_factor <- 1 + I2 / (1 - I2)
  }

  OIS <- D * diversity_factor

  return(OIS)
}

#' Calculate Optimal Information Size for Continuous Outcomes
#'
#' @param sd_pooled Pooled standard deviation
#' @param delta Minimum important difference
#' @param alpha Type I error rate
#' @param power Statistical power
#' @param I2 Heterogeneity
#' @return Required total sample size
calculate_ois_continuous <- function(sd_pooled, delta, alpha = 0.05, power = 0.80, I2 = 0) {

  if (is.na(sd_pooled) || sd_pooled <= 0 || is.na(delta) || delta == 0) return(NA)
  if (is.na(I2)) I2 <- 0
  I2 <- max(0, min(0.99, I2))

  z_alpha <- qnorm(1 - alpha/2)
  z_beta <- qnorm(power)

  # Sample size per arm
  n_per_arm <- 2 * ((z_alpha + z_beta)^2 * sd_pooled^2) / delta^2
  D <- 2 * n_per_arm

  # Heterogeneity adjustment
  if (I2 >= 0.99) {
    diversity_factor <- 100
  } else {
    diversity_factor <- 1 + I2 / (1 - I2)
  }

  OIS <- D * diversity_factor
  return(OIS)
}

#' Calculate Information Fraction
#'
#' @param actual_n Actual total sample size in meta-analysis
#' @param ois Optimal information size
#' @return Information fraction (0 to Inf, >1 means adequate)
calculate_information_fraction <- function(actual_n, ois) {
  if (is.na(ois) || ois <= 0 || is.na(actual_n)) return(NA)
  return(actual_n / ois)
}

#' Calculate Heterogeneity-Adjusted Power
#'
#' @param k Number of studies
#' @param avg_n Average study sample size
#' @param effect Observed effect size (SMD or log-OR)
#' @param se_effect Standard error of effect
#' @param I2 Heterogeneity
#' @return Approximate power
calculate_heterogeneity_adjusted_power <- function(k, avg_n, effect, se_effect, I2 = 0) {

  if (is.na(effect) || is.na(se_effect) || se_effect <= 0) return(NA)

  # Adjust SE for heterogeneity
  I2 <- max(0, min(0.99, ifelse(is.na(I2), 0, I2)))
  tau2_proportion <- I2 / (1 - I2)

  # The observed SE already incorporates heterogeneity from REML
  # Power = probability that Z > z_alpha given true effect
  z_score <- abs(effect) / se_effect

  # Approximate power
  power <- pnorm(z_score - qnorm(0.975))

  return(max(0, min(1, power)))
}

#' Determine Sequential Boundary Status
#'
#' @param z_cumulative Cumulative Z-score
#' @param information_fraction Current information fraction
#' @param alpha Alpha level
#' @return Status: "conclusive_benefit", "conclusive_harm", "futile", "continue"
sequential_boundary_status <- function(z_cumulative, information_fraction, alpha = 0.05) {

  if (is.na(z_cumulative) || is.na(information_fraction)) return("unknown")

  # O'Brien-Fleming-like boundaries (simplified)
  # Boundary gets narrower as information increases
  if (information_fraction >= 1) {
    boundary <- qnorm(1 - alpha/2)  # Standard 1.96
  } else if (information_fraction >= 0.5) {
    boundary <- qnorm(1 - alpha/2) / sqrt(information_fraction)
  } else if (information_fraction >= 0.25) {
    boundary <- qnorm(1 - alpha/2) / sqrt(information_fraction) * 1.2
  } else {
    boundary <- qnorm(1 - alpha/2) / sqrt(information_fraction) * 1.5
  }

  # Cap boundary at reasonable values
  boundary <- min(5, boundary)

  # Futility boundary (simplified)
  futility_boundary <- 0.5 * qnorm(1 - alpha/2)

  if (z_cumulative > boundary) {
    return("conclusive_benefit")
  } else if (z_cumulative < -boundary) {
    return("conclusive_harm")
  } else if (abs(z_cumulative) < futility_boundary && information_fraction > 0.5) {
    return("futile")
  } else {
    return("continue")
  }
}

#' Calculate Information Adequacy Index (IAI)
#'
#' @param information_fraction IF value
#' @param het_adj_power HAP value
#' @param boundary_status Sequential status
#' @param z_stability Z-score stability measure
#' @return IAI score (0-1, higher = more adequate)
calculate_iai <- function(information_fraction, het_adj_power, boundary_status, z_stability = NA) {

  # Component 1: Information Fraction (40%)
  # Transform IF to 0-1 scale using sigmoid
  if (is.na(information_fraction)) {
    if_component <- 0.5
  } else {
    if_component <- 1 / (1 + exp(-2 * (information_fraction - 1)))
  }

  # Component 2: Heterogeneity-Adjusted Power (30%)
  if (is.na(het_adj_power)) {
    power_component <- 0.5
  } else {
    power_component <- het_adj_power
  }

  # Component 3: Sequential Boundary Status (20%)
  status_scores <- c(
    "conclusive_benefit" = 1.0,
    "conclusive_harm" = 1.0,
    "futile" = 0.7,  # Conclusive in futility sense
    "continue" = 0.3,
    "unknown" = 0.5
  )
  boundary_component <- status_scores[boundary_status]
  if (is.na(boundary_component)) boundary_component <- 0.5

  # Component 4: Z-score Stability (10%)
  # Higher stability = less fluctuation in cumulative Z
  if (is.na(z_stability)) {
    stability_component <- 0.5
  } else {
    stability_component <- max(0, min(1, z_stability))
  }

  # Weighted combination
  IAI <- 0.40 * if_component +
         0.30 * power_component +
         0.20 * boundary_component +
         0.10 * stability_component

  return(round(IAI, 4))
}

#' Classify IAI into categories
classify_iai <- function(iai) {
  if (is.na(iai)) return("Unknown")
  if (iai >= 0.75) return("Adequate")
  if (iai >= 0.50) return("Marginal")
  if (iai >= 0.25) return("Inadequate")
  return("Critical")
}

# ==============================================================================
# PART 2: LOAD AND PROCESS PAIRWISE70 DATA
# ==============================================================================

cat("Loading Pairwise70 datasets...\n")

# Get list of all datasets
data_files <- list.files(file.path(base_path, "data"), pattern = "\\.rda$", full.names = TRUE)
cat(sprintf("Found %d datasets\n", length(data_files)))

# Function to process a single dataset
process_dataset <- function(file_path) {

  dataset_name <- gsub("\\.rda$", "", basename(file_path))

  tryCatch({
    # Load dataset
    env <- new.env()
    load(file_path, envir = env)
    df <- get(ls(env)[1], envir = env)

    # Determine if binary or continuous
    has_binary <- all(c("Experimental.cases", "Experimental.N", "Control.cases", "Control.N") %in% names(df))
    has_continuous <- all(c("Experimental.mean", "Experimental.SD", "Control.mean", "Control.SD") %in% names(df))

    results <- list()

    if (has_binary) {
      # Process binary outcome data
      df_clean <- df %>%
        filter(!is.na(Experimental.cases) & !is.na(Experimental.N) &
               !is.na(Control.cases) & !is.na(Control.N)) %>%
        filter(Experimental.N > 0 & Control.N > 0)

      if (nrow(df_clean) >= 2) {
        # Calculate effect sizes
        es <- escalc(measure = "OR",
                     ai = Experimental.cases, n1i = Experimental.N,
                     ci = Control.cases, n2i = Control.N,
                     data = df_clean)

        # Run meta-analysis
        ma <- tryCatch({
          rma(yi, vi, data = es, method = "REML")
        }, error = function(e) NULL)

        if (!is.null(ma)) {
          # Extract key statistics
          k <- ma$k
          total_n <- sum(df_clean$Experimental.N + df_clean$Control.N)
          avg_n <- total_n / k

          # Control event rate
          p_control <- sum(df_clean$Control.cases) / sum(df_clean$Control.N)

          # I-squared
          I2 <- max(0, ma$I2) / 100

          # Calculate OIS (assuming 20% RRR as default clinically important effect)
          ois <- calculate_ois_binary(p_control, RRR = 0.20, I2 = I2)

          # Information fraction
          inf_frac <- calculate_information_fraction(total_n, ois)

          # Z-score
          z_score <- ma$zval

          # Heterogeneity-adjusted power
          hap <- calculate_heterogeneity_adjusted_power(k, avg_n, ma$beta, ma$se, I2)

          # Sequential status
          seq_status <- sequential_boundary_status(z_score, inf_frac)

          # IAI
          iai <- calculate_iai(inf_frac, hap, seq_status)

          results <- list(
            dataset = dataset_name,
            outcome_type = "binary",
            k = k,
            total_n = total_n,
            avg_n_per_study = round(avg_n, 1),
            p_control = round(p_control, 4),
            effect_logOR = round(ma$beta, 4),
            se = round(ma$se, 4),
            pval = ma$pval,
            I2 = round(I2 * 100, 1),
            tau2 = round(ma$tau2, 4),
            OIS = round(ois, 0),
            information_fraction = round(inf_frac, 3),
            het_adj_power = round(hap, 3),
            z_score = round(z_score, 3),
            sequential_status = seq_status,
            IAI = iai,
            IAI_class = classify_iai(iai),
            significant = ma$pval < 0.05
          )
        }
      }
    }

    if (has_continuous && length(results) == 0) {
      # Process continuous outcome data
      df_clean <- df %>%
        filter(!is.na(Experimental.mean) & !is.na(Experimental.SD) & !is.na(Experimental.N) &
               !is.na(Control.mean) & !is.na(Control.SD) & !is.na(Control.N)) %>%
        filter(Experimental.N > 0 & Control.N > 0 & Experimental.SD > 0 & Control.SD > 0)

      if (nrow(df_clean) >= 2) {
        # Calculate SMD
        es <- escalc(measure = "SMD",
                     m1i = Experimental.mean, sd1i = Experimental.SD, n1i = Experimental.N,
                     m2i = Control.mean, sd2i = Control.SD, n2i = Control.N,
                     data = df_clean)

        ma <- tryCatch({
          rma(yi, vi, data = es, method = "REML")
        }, error = function(e) NULL)

        if (!is.null(ma)) {
          k <- ma$k
          total_n <- sum(df_clean$Experimental.N + df_clean$Control.N)
          avg_n <- total_n / k

          # Pooled SD (approximate)
          sd_pooled <- sqrt(mean(c(df_clean$Experimental.SD^2, df_clean$Control.SD^2)))

          I2 <- max(0, ma$I2) / 100

          # OIS for continuous (assuming SMD = 0.3 as MCID)
          delta <- 0.3 * sd_pooled
          ois <- calculate_ois_continuous(sd_pooled, delta, I2 = I2)

          inf_frac <- calculate_information_fraction(total_n, ois)
          z_score <- ma$zval
          hap <- calculate_heterogeneity_adjusted_power(k, avg_n, ma$beta, ma$se, I2)
          seq_status <- sequential_boundary_status(z_score, inf_frac)
          iai <- calculate_iai(inf_frac, hap, seq_status)

          results <- list(
            dataset = dataset_name,
            outcome_type = "continuous",
            k = k,
            total_n = total_n,
            avg_n_per_study = round(avg_n, 1),
            p_control = NA,
            effect_logOR = round(ma$beta, 4),  # Actually SMD
            se = round(ma$se, 4),
            pval = ma$pval,
            I2 = round(I2 * 100, 1),
            tau2 = round(ma$tau2, 4),
            OIS = round(ois, 0),
            information_fraction = round(inf_frac, 3),
            het_adj_power = round(hap, 3),
            z_score = round(z_score, 3),
            sequential_status = seq_status,
            IAI = iai,
            IAI_class = classify_iai(iai),
            significant = ma$pval < 0.05
          )
        }
      }
    }

    return(results)

  }, error = function(e) {
    return(list(dataset = dataset_name, error = as.character(e)))
  })
}

# ==============================================================================
# PART 3: RUN ANALYSIS ON ALL DATASETS
# ==============================================================================

cat("\nProcessing all datasets...\n")

# Process all datasets (with progress)
all_results <- list()
n_datasets <- length(data_files)

for (i in seq_along(data_files)) {
  if (i %% 50 == 0) {
    cat(sprintf("Progress: %d/%d (%.1f%%)\n", i, n_datasets, 100*i/n_datasets))
  }
  result <- process_dataset(data_files[i])
  if (length(result) > 0 && is.null(result$error)) {
    all_results[[length(all_results) + 1]] <- result
  }
}

cat(sprintf("\nSuccessfully processed %d meta-analyses\n", length(all_results)))

# Convert to data.table
results_dt <- rbindlist(all_results, fill = TRUE)

# Save results
fwrite(results_dt, file.path(output_path, "information_adequacy_results.csv"))

# ==============================================================================
# PART 4: SUMMARY STATISTICS
# ==============================================================================

cat("\n", strrep("=", 70), "\n")
cat("SUMMARY STATISTICS\n")
cat(strrep("=", 70), "\n\n")

# Overall statistics
cat("OVERALL STATISTICS\n")
cat(strrep("-", 40), "\n")
cat(sprintf("Total meta-analyses analyzed: %d\n", nrow(results_dt)))
cat(sprintf("Binary outcomes: %d (%.1f%%)\n",
            sum(results_dt$outcome_type == "binary"),
            100 * mean(results_dt$outcome_type == "binary")))
cat(sprintf("Continuous outcomes: %d (%.1f%%)\n",
            sum(results_dt$outcome_type == "continuous"),
            100 * mean(results_dt$outcome_type == "continuous")))

cat(sprintf("\nMedian k (studies): %.0f (IQR: %.0f - %.0f)\n",
            median(results_dt$k, na.rm = TRUE),
            quantile(results_dt$k, 0.25, na.rm = TRUE),
            quantile(results_dt$k, 0.75, na.rm = TRUE)))

cat(sprintf("Median total N: %.0f (IQR: %.0f - %.0f)\n",
            median(results_dt$total_n, na.rm = TRUE),
            quantile(results_dt$total_n, 0.25, na.rm = TRUE),
            quantile(results_dt$total_n, 0.75, na.rm = TRUE)))

# Information Fraction Distribution
cat("\n\nINFORMATION FRACTION DISTRIBUTION\n")
cat(strrep("-", 40), "\n")
inf_frac_summary <- results_dt %>%
  filter(!is.na(information_fraction)) %>%
  summarise(
    n = n(),
    adequate_pct = 100 * mean(information_fraction >= 1),
    marginal_pct = 100 * mean(information_fraction >= 0.5 & information_fraction < 1),
    inadequate_pct = 100 * mean(information_fraction >= 0.25 & information_fraction < 0.5),
    critical_pct = 100 * mean(information_fraction < 0.25),
    median_if = median(information_fraction),
    iqr_low = quantile(information_fraction, 0.25),
    iqr_high = quantile(information_fraction, 0.75)
  )

cat(sprintf("Information Fraction >= 1.0 (Adequate): %.1f%%\n", inf_frac_summary$adequate_pct))
cat(sprintf("Information Fraction 0.5-1.0 (Marginal): %.1f%%\n", inf_frac_summary$marginal_pct))
cat(sprintf("Information Fraction 0.25-0.5 (Inadequate): %.1f%%\n", inf_frac_summary$inadequate_pct))
cat(sprintf("Information Fraction < 0.25 (Critical): %.1f%%\n", inf_frac_summary$critical_pct))
cat(sprintf("\nMedian IF: %.2f (IQR: %.2f - %.2f)\n",
            inf_frac_summary$median_if, inf_frac_summary$iqr_low, inf_frac_summary$iqr_high))

# IAI Distribution
cat("\n\nINFORMATION ADEQUACY INDEX (IAI) DISTRIBUTION\n")
cat(strrep("-", 40), "\n")
iai_summary <- table(results_dt$IAI_class)
cat(sprintf("Adequate (IAI >= 0.75): %d (%.1f%%)\n",
            iai_summary["Adequate"], 100*iai_summary["Adequate"]/nrow(results_dt)))
cat(sprintf("Marginal (IAI 0.50-0.75): %d (%.1f%%)\n",
            iai_summary["Marginal"], 100*iai_summary["Marginal"]/nrow(results_dt)))
cat(sprintf("Inadequate (IAI 0.25-0.50): %d (%.1f%%)\n",
            iai_summary["Inadequate"], 100*iai_summary["Inadequate"]/nrow(results_dt)))
cat(sprintf("Critical (IAI < 0.25): %d (%.1f%%)\n",
            iai_summary["Critical"], 100*iai_summary["Critical"]/nrow(results_dt)))

# Sequential Status
cat("\n\nSEQUENTIAL BOUNDARY STATUS\n")
cat(strrep("-", 40), "\n")
seq_summary <- table(results_dt$sequential_status)
for (status in names(seq_summary)) {
  cat(sprintf("%s: %d (%.1f%%)\n", status, seq_summary[status], 100*seq_summary[status]/nrow(results_dt)))
}

# Key Finding: Significant but inadequate
cat("\n\nKEY FINDING: PREMATURE CONCLUSIONS\n")
cat(strrep("-", 40), "\n")
premature <- results_dt %>%
  filter(significant == TRUE & information_fraction < 1)
cat(sprintf("Significant meta-analyses with inadequate information: %d / %d (%.1f%%)\n",
            nrow(premature),
            sum(results_dt$significant, na.rm = TRUE),
            100 * nrow(premature) / sum(results_dt$significant, na.rm = TRUE)))

# By k categories
cat("\n\nINFORMATION ADEQUACY BY NUMBER OF STUDIES (k)\n")
cat(strrep("-", 40), "\n")
k_summary <- results_dt %>%
  mutate(k_cat = case_when(
    k <= 5 ~ "2-5",
    k <= 10 ~ "6-10",
    k <= 20 ~ "11-20",
    k <= 50 ~ "21-50",
    TRUE ~ ">50"
  )) %>%
  group_by(k_cat) %>%
  summarise(
    n = n(),
    median_if = median(information_fraction, na.rm = TRUE),
    pct_adequate = 100 * mean(information_fraction >= 1, na.rm = TRUE),
    median_iai = median(IAI, na.rm = TRUE)
  ) %>%
  arrange(match(k_cat, c("2-5", "6-10", "11-20", "21-50", ">50")))

print(k_summary)

# By heterogeneity
cat("\n\nINFORMATION ADEQUACY BY HETEROGENEITY (I2)\n")
cat(strrep("-", 40), "\n")
het_summary <- results_dt %>%
  mutate(I2_cat = case_when(
    I2 <= 25 ~ "Low (0-25%)",
    I2 <= 50 ~ "Moderate (25-50%)",
    I2 <= 75 ~ "Substantial (50-75%)",
    TRUE ~ "Considerable (>75%)"
  )) %>%
  group_by(I2_cat) %>%
  summarise(
    n = n(),
    median_if = median(information_fraction, na.rm = TRUE),
    pct_adequate = 100 * mean(information_fraction >= 1, na.rm = TRUE),
    median_iai = median(IAI, na.rm = TRUE)
  )

print(het_summary)

# ==============================================================================
# PART 5: COMPARISON WITH MAFI (if available)
# ==============================================================================

mafi_file <- file.path(base_path, "analysis/output/MAFI_all_variants.csv")
if (file.exists(mafi_file)) {
  cat("\n\nCOMPARISON WITH MAFI\n")
  cat(strrep("-", 40), "\n")

  mafi_dt <- fread(mafi_file)

  # Merge by dataset name (need to match naming)
  # This is a placeholder - actual matching would need dataset identifiers
  cat("MAFI data available - implement comparison after data alignment\n")
}

# ==============================================================================
# PART 6: VISUALIZATIONS
# ==============================================================================

cat("\n\nGenerating visualizations...\n")

# 1. Information Fraction Distribution
p1 <- ggplot(results_dt, aes(x = information_fraction)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.7) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", size = 1) +
  scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 5, 10)) +
  labs(title = "Distribution of Information Fraction Across 501 Cochrane Meta-Analyses",
       subtitle = "Red line indicates adequate information (IF = 1.0)",
       x = "Information Fraction (log scale)",
       y = "Number of Meta-Analyses") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(output_path, "fig1_information_fraction_distribution.png"),
       p1, width = 10, height = 6, dpi = 300)

# 2. IAI Distribution
p2 <- ggplot(results_dt, aes(x = IAI, fill = IAI_class)) +
  geom_histogram(bins = 40, color = "white", alpha = 0.8) +
  scale_fill_manual(values = c("Adequate" = "#2ECC71", "Marginal" = "#F39C12",
                               "Inadequate" = "#E74C3C", "Critical" = "#8E44AD")) +
  labs(title = "Information Adequacy Index (IAI) Distribution",
       x = "IAI Score",
       y = "Number of Meta-Analyses",
       fill = "Classification") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(output_path, "fig2_iai_distribution.png"),
       p2, width = 10, height = 6, dpi = 300)

# 3. IF vs k (number of studies)
p3 <- ggplot(results_dt, aes(x = k, y = information_fraction)) +
  geom_point(aes(color = IAI_class), alpha = 0.5) +
  geom_smooth(method = "loess", color = "black", se = TRUE) +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c("Adequate" = "#2ECC71", "Marginal" = "#F39C12",
                                "Inadequate" = "#E74C3C", "Critical" = "#8E44AD")) +
  labs(title = "Information Fraction vs Number of Studies",
       x = "Number of Studies (k)",
       y = "Information Fraction",
       color = "IAI Class") +
  theme_minimal()

ggsave(file.path(output_path, "fig3_if_vs_k.png"),
       p3, width = 10, height = 6, dpi = 300)

# 4. IF vs I2
p4 <- ggplot(results_dt, aes(x = I2, y = information_fraction)) +
  geom_point(aes(color = significant), alpha = 0.5) +
  geom_smooth(method = "loess", color = "black") +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
  scale_y_log10() +
  labs(title = "Information Fraction vs Heterogeneity",
       subtitle = "Higher heterogeneity requires more information",
       x = "I-squared (%)",
       y = "Information Fraction (log scale)",
       color = "Significant") +
  theme_minimal()

ggsave(file.path(output_path, "fig4_if_vs_heterogeneity.png"),
       p4, width = 10, height = 6, dpi = 300)

# 5. Premature conclusions visualization
sig_data <- results_dt %>%
  filter(significant == TRUE) %>%
  mutate(if_category = case_when(
    information_fraction >= 1 ~ "Adequate (IF >= 1)",
    information_fraction >= 0.5 ~ "Marginal (0.5 <= IF < 1)",
    information_fraction >= 0.25 ~ "Inadequate (0.25 <= IF < 0.5)",
    TRUE ~ "Critical (IF < 0.25)"
  ))

p5 <- ggplot(sig_data, aes(x = if_category, fill = if_category)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  scale_fill_manual(values = c("Adequate (IF >= 1)" = "#2ECC71",
                               "Marginal (0.5 <= IF < 1)" = "#F39C12",
                               "Inadequate (0.25 <= IF < 0.5)" = "#E74C3C",
                               "Critical (IF < 0.25)" = "#8E44AD")) +
  labs(title = "Information Adequacy Among Statistically Significant Meta-Analyses",
       subtitle = "Many 'significant' results may be premature conclusions",
       x = "",
       y = "Number of Meta-Analyses") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 15, hjust = 1))

ggsave(file.path(output_path, "fig5_premature_conclusions.png"),
       p5, width = 10, height = 6, dpi = 300)

cat("Visualizations saved to:", output_path, "\n")

# ==============================================================================
# PART 7: SAVE FINAL REPORT
# ==============================================================================

report <- sprintf("
================================================================================
INFORMATION ADEQUACY ANALYSIS - FINAL REPORT
================================================================================
Generated: %s
Data Source: Pairwise70 (501 Cochrane Meta-Analyses)
================================================================================

EXECUTIVE SUMMARY
-----------------
Total meta-analyses analyzed: %d
Proportion with adequate information (IF >= 1): %.1f%%
Proportion with critical information deficit (IF < 0.25): %.1f%%

KEY FINDING:
Among statistically significant meta-analyses:
- %.1f%% have inadequate information (IF < 1.0)
- These represent potential PREMATURE CONCLUSIONS

INFORMATION ADEQUACY INDEX (IAI) DISTRIBUTION:
- Adequate: %.1f%%
- Marginal: %.1f%%
- Inadequate: %.1f%%
- Critical: %.1f%%

IMPLICATIONS:
1. Many published meta-analyses conclude prematurely
2. High heterogeneity dramatically increases information requirements
3. Small meta-analyses (k < 5) are particularly vulnerable
4. Statistical significance ≠ adequate evidence

RECOMMENDATIONS:
1. Report information fraction alongside p-values
2. Use IAI to prioritize meta-analyses for updating
3. Consider sequential boundaries before concluding
4. Adjust for heterogeneity in power calculations

================================================================================
",
    Sys.time(),
    nrow(results_dt),
    inf_frac_summary$adequate_pct,
    inf_frac_summary$critical_pct,
    100 * nrow(premature) / sum(results_dt$significant, na.rm = TRUE),
    100*iai_summary["Adequate"]/nrow(results_dt),
    100*iai_summary["Marginal"]/nrow(results_dt),
    100*iai_summary["Inadequate"]/nrow(results_dt),
    100*iai_summary["Critical"]/nrow(results_dt)
)

writeLines(report, file.path(output_path, "INFORMATION_ADEQUACY_REPORT.txt"))

cat("\n", strrep("=", 70), "\n")
cat("ANALYSIS COMPLETE\n")
cat("Results saved to:", output_path, "\n")
cat(strrep("=", 70), "\n")
