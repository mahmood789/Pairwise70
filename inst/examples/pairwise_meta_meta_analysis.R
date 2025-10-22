# Meta-Meta-Analysis of Pairwise70 Package
# Validates all 501 Cochrane pairwise meta-analysis datasets
# Extracts heterogeneity patterns and effect size distributions

library(Pairwise70)
library(metafor)

cat("=== PAIRWISE70 META-META-ANALYSIS ===\n\n")
cat("Analyzing all 501 Cochrane pairwise meta-analyses...\n\n")

# Get all datasets
all_datasets <- data(package = "Pairwise70")$results[,3]
total_datasets <- length(all_datasets)

cat(paste0("Total datasets to analyze: ", total_datasets, "\n\n"))

# Initialize results storage
results <- data.frame(
  dataset_name = character(),
  n_studies = integer(),
  outcome_type = character(),
  effect_measure = character(),
  pooled_effect = numeric(),
  pooled_se = numeric(),
  pooled_ci_lb = numeric(),
  pooled_ci_ub = numeric(),
  i2 = numeric(),
  tau2 = numeric(),
  q_stat = numeric(),
  q_pval = numeric(),
  success = logical(),
  error_message = character(),
  stringsAsFactors = FALSE
)

# Track progress
successful <- 0
failed <- 0
binary_count <- 0
continuous_count <- 0
insufficient_data <- 0

cat("Processing datasets...\n")
pb <- txtProgressBar(min = 0, max = total_datasets, style = 3)

for (i in seq_along(all_datasets)) {
  ds_name <- all_datasets[i]

  # Update progress bar
  setTxtProgressBar(pb, i)

  tryCatch({
    # Load dataset
    data(list = ds_name, package = "Pairwise70", envir = environment())
    dataset <- get(ds_name, envir = environment())

    # Determine outcome type (use lowercase column names)
    has_binary_cols <- all(c("Experimental.cases", "Experimental.N", "Control.cases", "Control.N") %in% names(dataset))
    has_continuous_cols <- all(c("Experimental.mean", "Experimental.SD", "Control.mean", "Control.SD") %in% names(dataset))

    # Check which columns have actual data (not just NAs)
    binary_has_data <- FALSE
    continuous_has_data <- FALSE

    if (has_binary_cols) {
      binary_has_data <- any(!is.na(dataset$Experimental.cases)) || any(!is.na(dataset$Control.cases))
    }

    if (has_continuous_cols) {
      continuous_has_data <- any(!is.na(dataset$Experimental.mean)) || any(!is.na(dataset$Control.mean))
    }

    # Prioritize binary over continuous when both exist (most datasets have both columns but only use one)
    if (binary_has_data) {
      # Binary outcome meta-analysis
      outcome_type <- "binary"
      binary_count <- binary_count + 1

      # Remove rows with missing data
      complete_data <- dataset[complete.cases(dataset[, c("Experimental.cases", "Experimental.N", "Control.cases", "Control.N")]), ]

      if (nrow(complete_data) < 2) {
        insufficient_data <- insufficient_data + 1
        results <- rbind(results, data.frame(
          dataset_name = ds_name,
          n_studies = nrow(dataset),
          outcome_type = outcome_type,
          effect_measure = "OR",
          pooled_effect = NA,
          pooled_se = NA,
          pooled_ci_lb = NA,
          pooled_ci_ub = NA,
          i2 = NA,
          tau2 = NA,
          q_stat = NA,
          q_pval = NA,
          success = FALSE,
          error_message = "Insufficient complete data (<2 studies)"
        ))
        next
      }

      # Run meta-analysis with odds ratio
      ma <- rma(
        measure = "OR",
        ai = Experimental.cases,
        n1i = Experimental.N,
        ci = Control.cases,
        n2i = Control.N,
        data = complete_data,
        method = "REML"
      )

      results <- rbind(results, data.frame(
        dataset_name = ds_name,
        n_studies = nrow(complete_data),
        outcome_type = outcome_type,
        effect_measure = "OR",
        pooled_effect = as.numeric(ma$beta),
        pooled_se = ma$se,
        pooled_ci_lb = ma$ci.lb,
        pooled_ci_ub = ma$ci.ub,
        i2 = ma$I2,
        tau2 = ma$tau2,
        q_stat = ma$QE,
        q_pval = ma$QEp,
        success = TRUE,
        error_message = ""
      ))

      successful <- successful + 1

    } else if (continuous_has_data) {
      # Continuous outcome meta-analysis
      outcome_type <- "continuous"
      continuous_count <- continuous_count + 1

      # Remove rows with missing data
      complete_data <- dataset[complete.cases(dataset[, c("Experimental.mean", "Experimental.SD", "Experimental.N",
                                                            "Control.mean", "Control.SD", "Control.N")]), ]

      if (nrow(complete_data) < 2) {
        insufficient_data <- insufficient_data + 1
        results <- rbind(results, data.frame(
          dataset_name = ds_name,
          n_studies = nrow(dataset),
          outcome_type = outcome_type,
          effect_measure = "SMD",
          pooled_effect = NA,
          pooled_se = NA,
          pooled_ci_lb = NA,
          pooled_ci_ub = NA,
          i2 = NA,
          tau2 = NA,
          q_stat = NA,
          q_pval = NA,
          success = FALSE,
          error_message = "Insufficient complete data (<2 studies)"
        ))
        next
      }

      # Run meta-analysis with standardized mean difference
      ma <- rma(
        measure = "SMD",
        m1i = Experimental.mean,
        sd1i = Experimental.SD,
        n1i = Experimental.N,
        m2i = Control.mean,
        sd2i = Control.SD,
        n2i = Control.N,
        data = complete_data,
        method = "REML"
      )

      results <- rbind(results, data.frame(
        dataset_name = ds_name,
        n_studies = nrow(complete_data),
        outcome_type = outcome_type,
        effect_measure = "SMD",
        pooled_effect = as.numeric(ma$beta),
        pooled_se = ma$se,
        pooled_ci_lb = ma$ci.lb,
        pooled_ci_ub = ma$ci.ub,
        i2 = ma$I2,
        tau2 = ma$tau2,
        q_stat = ma$QE,
        q_pval = ma$QEp,
        success = TRUE,
        error_message = ""
      ))

      successful <- successful + 1

    } else {
      # No usable data
      results <- rbind(results, data.frame(
        dataset_name = ds_name,
        n_studies = nrow(dataset),
        outcome_type = "no_data",
        effect_measure = NA,
        pooled_effect = NA,
        pooled_se = NA,
        pooled_ci_lb = NA,
        pooled_ci_ub = NA,
        i2 = NA,
        tau2 = NA,
        q_stat = NA,
        q_pval = NA,
        success = FALSE,
        error_message = "No usable outcome data found"
      ))

      failed <- failed + 1
    }

  }, error = function(e) {
    results <<- rbind(results, data.frame(
      dataset_name = ds_name,
      n_studies = NA,
      outcome_type = NA,
      effect_measure = NA,
      pooled_effect = NA,
      pooled_se = NA,
      pooled_ci_lb = NA,
      pooled_ci_ub = NA,
      i2 = NA,
      tau2 = NA,
      q_stat = NA,
      q_pval = NA,
      success = FALSE,
      error_message = as.character(e$message)
    ))

    failed <<- failed + 1
  })
}

close(pb)

cat("\n\n=== META-META-ANALYSIS COMPLETE ===\n\n")

# Summary statistics
cat("DATASET VALIDATION:\n")
cat(paste0("  Total datasets analyzed: ", total_datasets, "\n"))
cat(paste0("  Successfully analyzed: ", successful, " (", round(100*successful/total_datasets, 1), "%)\n"))
cat(paste0("  Failed: ", failed, "\n"))
cat(paste0("  Insufficient data (<2 studies): ", insufficient_data, "\n\n"))

cat("OUTCOME TYPES:\n")
cat(paste0("  Binary outcomes: ", binary_count, " (", round(100*binary_count/total_datasets, 1), "%)\n"))
cat(paste0("  Continuous outcomes: ", continuous_count, " (", round(100*continuous_count/total_datasets, 1), "%)\n\n"))

# Analyze heterogeneity patterns (only successful meta-analyses)
successful_results <- results[results$success == TRUE, ]

if (nrow(successful_results) > 0) {
  cat("HETEROGENEITY PATTERNS (across ", nrow(successful_results), " successful meta-analyses):\n", sep="")
  cat(paste0("  Median I²: ", round(median(successful_results$i2, na.rm = TRUE), 1), "%\n"))
  cat(paste0("  Mean I²: ", round(mean(successful_results$i2, na.rm = TRUE), 1), "%\n"))
  cat(paste0("  I² range: ", round(min(successful_results$i2, na.rm = TRUE), 1), "% - ",
             round(max(successful_results$i2, na.rm = TRUE), 1), "%\n\n"))

  cat("  I² categories:\n")
  low_het <- sum(successful_results$i2 < 25, na.rm = TRUE)
  moderate_het <- sum(successful_results$i2 >= 25 & successful_results$i2 < 50, na.rm = TRUE)
  substantial_het <- sum(successful_results$i2 >= 50 & successful_results$i2 < 75, na.rm = TRUE)
  considerable_het <- sum(successful_results$i2 >= 75, na.rm = TRUE)

  cat(paste0("    Low (I² < 25%): ", low_het, " (", round(100*low_het/nrow(successful_results), 1), "%)\n"))
  cat(paste0("    Moderate (25% ≤ I² < 50%): ", moderate_het, " (", round(100*moderate_het/nrow(successful_results), 1), "%)\n"))
  cat(paste0("    Substantial (50% ≤ I² < 75%): ", substantial_het, " (", round(100*substantial_het/nrow(successful_results), 1), "%)\n"))
  cat(paste0("    Considerable (I² ≥ 75%): ", considerable_het, " (", round(100*considerable_het/nrow(successful_results), 1), "%)\n\n"))

  cat("STUDY COUNTS:\n")
  cat(paste0("  Median studies per meta-analysis: ", median(successful_results$n_studies, na.rm = TRUE), "\n"))
  cat(paste0("  Mean studies per meta-analysis: ", round(mean(successful_results$n_studies, na.rm = TRUE), 1), "\n"))
  cat(paste0("  Range: ", min(successful_results$n_studies, na.rm = TRUE), " - ",
             max(successful_results$n_studies, na.rm = TRUE), " studies\n\n"))

  # Total studies analyzed
  total_studies <- sum(successful_results$n_studies, na.rm = TRUE)
  cat(paste0("  Total individual studies analyzed: ", format(total_studies, big.mark = ","), "\n\n"))

  # Visualization: I² distribution
  cat("Generating I² distribution histogram...\n")
  png("pairwise_i2_distribution.png", width = 800, height = 600)
  hist(successful_results$i2,
       breaks = 30,
       main = "Distribution of I² Across 501 Cochrane Meta-Analyses",
       xlab = "I² (%)",
       ylab = "Frequency",
       col = "steelblue",
       border = "white")
  abline(v = median(successful_results$i2, na.rm = TRUE), col = "red", lwd = 2, lty = 2)
  legend("topright",
         legend = paste0("Median I² = ", round(median(successful_results$i2, na.rm = TRUE), 1), "%"),
         col = "red", lty = 2, lwd = 2)
  dev.off()
  cat("Saved: pairwise_i2_distribution.png\n\n")

  # Visualization: Study count distribution
  cat("Generating study count distribution histogram...\n")
  png("pairwise_study_count_distribution.png", width = 800, height = 600)
  hist(successful_results$n_studies,
       breaks = 30,
       main = "Distribution of Study Counts Across Meta-Analyses",
       xlab = "Number of Studies",
       ylab = "Frequency",
       col = "forestgreen",
       border = "white")
  abline(v = median(successful_results$n_studies, na.rm = TRUE), col = "red", lwd = 2, lty = 2)
  legend("topright",
         legend = paste0("Median = ", median(successful_results$n_studies, na.rm = TRUE), " studies"),
         col = "red", lty = 2, lwd = 2)
  dev.off()
  cat("Saved: pairwise_study_count_distribution.png\n\n")
}

# Save full results
cat("Saving detailed results...\n")
write.csv(results, "pairwise_meta_meta_results.csv", row.names = FALSE)
cat("Saved: pairwise_meta_meta_results.csv\n\n")

# Show sample of failed datasets for debugging
if (failed > 0) {
  cat("SAMPLE OF FAILED DATASETS:\n")
  failed_results <- results[results$success == FALSE, ]
  print(head(failed_results[, c("dataset_name", "n_studies", "error_message")], 10))
  cat("\n")
}

cat("=== VALIDATION COMPLETE ===\n")
cat("All 501 Pairwise70 datasets have been tested!\n\n")

# Return results invisibly for further analysis
invisible(results)
