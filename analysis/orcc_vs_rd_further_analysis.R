#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[1] else file.path("analysis", "output")
input_path <- file.path(output_dir, "orcc_vs_rd_disagreements.csv")
dataset_out <- file.path(output_dir, "orcc_vs_rd_dataset_mismatch_counts.csv")
report_path <- file.path(output_dir, "orcc_vs_rd_further_report.md")

if (!file.exists(input_path)) {
  stop(paste0("Missing input file: ", input_path))
}

df <- read.csv(input_path, stringsAsFactors = FALSE)

df$dir_mismatch <- df$mismatch_type %in% c("direction_only", "direction_and_significance")
df$sig_mismatch <- df$mismatch_type %in% c("significance_only", "direction_and_significance")
df$any_mismatch <- df$mismatch_type %in% c("direction_only", "significance_only", "direction_and_significance")

k_vals <- as.numeric(df$k)
k_bin <- cut(
  k_vals,
  breaks = c(-Inf, 4, 9, 19, Inf),
  labels = c("2-4", "5-9", "10-19", "20+"),
  right = TRUE
)

dataset_summary <- aggregate(
  cbind(dir_mismatch, sig_mismatch, any_mismatch) ~ dataset_name,
  data = df,
  FUN = sum,
  na.rm = TRUE
)
dataset_counts <- aggregate(
  any_mismatch ~ dataset_name,
  data = df,
  FUN = length
)
names(dataset_counts)[2] <- "paired_analyses"
dataset_summary <- merge(dataset_summary, dataset_counts, by = "dataset_name")
dataset_summary$mismatch_rate <- round(dataset_summary$any_mismatch / dataset_summary$paired_analyses, 3)
dataset_summary <- dataset_summary[order(-dataset_summary$mismatch_rate, -dataset_summary$paired_analyses), ]
write.csv(dataset_summary, dataset_out, row.names = FALSE)

rate <- function(x, n) {
  if (n == 0) return(NA)
  round(x / n, 3)
}

total <- nrow(df)
dir_mismatch_n <- sum(df$dir_mismatch, na.rm = TRUE)
sig_mismatch_n <- sum(df$sig_mismatch, na.rm = TRUE)
any_mismatch_n <- sum(df$any_mismatch, na.rm = TRUE)

sparse_rate <- round(mean(df$any_mismatch[df$sparse_events], na.rm = TRUE), 3)
non_sparse_rate <- round(mean(df$any_mismatch[!df$sparse_events], na.rm = TRUE), 3)
double_zero_rate <- round(mean(df$any_mismatch[df$double_zero], na.rm = TRUE), 3)
non_double_zero_rate <- round(mean(df$any_mismatch[!df$double_zero], na.rm = TRUE), 3)

k_table <- aggregate(any_mismatch ~ k_bin, data = df, FUN = mean, na.rm = TRUE)
k_table$any_mismatch <- round(k_table$any_mismatch, 3)

top_datasets <- dataset_summary[dataset_summary$paired_analyses >= 5, ]
top_datasets <- head(top_datasets, 15)

report_lines <- c(
  "# OR_cc vs RD further analysis",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Overall mismatch rates",
  paste0("- Total paired analyses: ", total),
  paste0("- Any mismatch: ", any_mismatch_n, " (", rate(any_mismatch_n, total), ")"),
  paste0("- Direction mismatch: ", dir_mismatch_n, " (", rate(dir_mismatch_n, total), ")"),
  paste0("- Significance mismatch: ", sig_mismatch_n, " (", rate(sig_mismatch_n, total), ")"),
  "",
  "## Mismatch rates by data sparsity",
  paste0("- Sparse events: ", sparse_rate),
  paste0("- Non-sparse events: ", non_sparse_rate),
  paste0("- Double-zero present: ", double_zero_rate),
  paste0("- No double-zero: ", non_double_zero_rate),
  "",
  "## Mismatch rate by k bin",
  paste0("- ", paste(paste(k_table$k_bin, k_table$any_mismatch, sep = ": "), collapse = "\n- ")),
  "",
  "## Datasets with highest mismatch rate (paired analyses >= 5)",
  paste0(
    "- ",
    paste(
      paste(
        top_datasets$dataset_name,
        "paired=", top_datasets$paired_analyses,
        "any_mismatch=", top_datasets$any_mismatch,
        "rate=", top_datasets$mismatch_rate,
        sep = " "
      ),
      collapse = "\n- "
    )
  ),
  "",
  "## Output files",
  "- orcc_vs_rd_dataset_mismatch_counts.csv",
  "- orcc_vs_rd_further_report.md"
)

writeLines(report_lines, report_path)

cat("Further analysis written to: ", report_path, "\n", sep = "")
