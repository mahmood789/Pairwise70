# Merge review-level effects with target-population covariates
root <- normalizePath("C:/Users/user/OneDrive - NHS/Documents/Pairwise70", winslash = "/", mustWork = TRUE)
output_dir <- file.path(root, "analysis", "transportability")

review_path <- file.path(output_dir, "transportability_review_level.csv")
ctgov_path <- file.path(output_dir, "ctgov_target_covariates.csv")

if (!file.exists(ctgov_path)) {
  stop("Missing ctgov_target_covariates.csv. Run 07_ctgov_aggregate.R after fetching registry data.")
}

review_df <- read.csv(review_path, stringsAsFactors = FALSE)
ctgov_df <- read.csv(ctgov_path, stringsAsFactors = FALSE)

merged <- merge(review_df, ctgov_df, by = "review_id", all.x = TRUE)

merged$year_shift <- merged$mean_year - merged$start_year_mean
merged$enrollment_shift <- merged$mean_total_n - merged$enrollment_mean

merge_path <- file.path(output_dir, "transportability_target_merge.csv")
write.csv(merged, merge_path, row.names = FALSE)

summary_path <- file.path(output_dir, "transportability_target_summary.md")
reviews_with_target <- length(unique(merged$review_id[!is.na(merged$trial_count)]))
rows_with_target <- sum(!is.na(merged$trial_count))

summary_lines <- c(
  "# Transportability Target Merge Summary",
  "",
  sprintf("Reviews with target data: %d", reviews_with_target),
  sprintf("Rows with target data: %d", rows_with_target),
  "",
  "## Shift metrics (median)",
  sprintf("- Year shift (review minus registry): %.2f", median(merged$year_shift, na.rm = TRUE)),
  sprintf("- Enrollment shift (review minus registry): %.2f", median(merged$enrollment_shift, na.rm = TRUE))
)
writeLines(summary_lines, summary_path)

cat("Wrote:", merge_path, "\n")
cat("Wrote:", summary_path, "\n")
