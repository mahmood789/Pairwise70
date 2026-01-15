# Transportability domain report
root <- normalizePath("C:/Users/user/OneDrive - NHS/Documents/Pairwise70", winslash = "/", mustWork = TRUE)
output_dir <- file.path(root, "analysis", "transportability")

counts_path <- file.path(output_dir, "transportability_domain_counts.csv")
summary_path <- file.path(output_dir, "transportability_domain_cv_summary.md")
assign_path <- file.path(output_dir, "transportability_domain_assignments.csv")

if (!file.exists(counts_path) || !file.exists(summary_path) || !file.exists(assign_path)) {
  stop("Missing domain outputs. Run 10_transportability_domain_cv.R first.")
}

counts <- read.csv(counts_path, stringsAsFactors = FALSE)
assignments <- read.csv(assign_path, stringsAsFactors = FALSE)

mortality_rate <- if (nrow(assignments) > 0) {
  mean(assignments$mortality_tag, na.rm = TRUE)
} else {
  NA_real_
}

multi_label_rate <- if (nrow(assignments) > 0) {
  mean(assignments$multi_label, na.rm = TRUE)
} else {
  NA_real_
}

secondary_rate <- if (nrow(assignments) > 0) {
  mean(!is.na(assignments$domain_secondary) & assignments$domain_secondary != "", na.rm = TRUE)
} else {
  NA_real_
}

report_path <- file.path(output_dir, "transportability_domain_report.md")
lines <- c(
  "# Transportability Domain Report",
  "",
  "## Domain coverage",
  paste0("- ", counts$domain, ": ", counts$reviews, " reviews (", counts$rows, " rows, ", counts$reviews_in_list, " reviews_in_list)"),
  "",
  sprintf("## Mortality tag prevalence: %.2f", mortality_rate),
  sprintf("## Multi-label prevalence: %.2f", multi_label_rate),
  sprintf("## Secondary domain prevalence: %.2f", secondary_rate),
  "",
  "## Interpretation notes",
  "- Domains are assigned by keyword voting on Analysis.name; the top-scoring domain is primary.",
  "- Multi-label reviews carry a domain_list and optional secondary domain.",
  "- Mortality is treated as an outcome tag, not a domain.",
  "- Low-count domains are reported but excluded from CV if < 30 rows.",
  "",
  "## Recommended target-population covariates (open registries)",
  "- Study timing (start/completion year) to align clinical era.",
  "- Enrollment size distribution to approximate trial size/setting.",
  "- Sex eligibility proportions as a basic inclusion proxy.",
  "- For cardio/cerebrovascular: include hypertension or vascular terms in query refinement.",
  "- For infectious: include pathogen-specific keywords (e.g., HIV, hepatitis).",
  "- For mental health: include disorder-specific terms (e.g., depression, anxiety)."
)

writeLines(lines, report_path)
cat("Wrote:", report_path, "\n")
