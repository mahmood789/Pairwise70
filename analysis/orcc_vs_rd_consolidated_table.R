#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[1] else file.path("analysis", "output")
rec_path <- file.path(output_dir, "orcc_vs_rd_recommendations.csv")
results_path <- file.path(output_dir, "remediation_analysis_results.csv")
out_path <- file.path(output_dir, "orcc_vs_rd_consolidated_recommendations.csv")
report_path <- file.path(output_dir, "orcc_vs_rd_consolidated_report.md")

if (!file.exists(rec_path)) stop(paste0("Missing file: ", rec_path))
if (!file.exists(results_path)) stop(paste0("Missing file: ", results_path))

rec <- read.csv(rec_path, stringsAsFactors = FALSE)
res <- read.csv(results_path, stringsAsFactors = FALSE)

res <- res[res$method %in% c("OR_cc", "RD"), ]
res$key <- paste(res$dataset_name, res$analysis_key, sep = "|")

orcc <- res[res$method == "OR_cc", ]
rd <- res[res$method == "RD", ]

orcc$key <- paste(orcc$dataset_name, orcc$analysis_key, sep = "|")
rd$key <- paste(rd$dataset_name, rd$analysis_key, sep = "|")

orcc_ok <- orcc[orcc$meta_status == "ok", ]
rd_ok <- rd[rd$meta_status == "ok", ]

orcc_map <- split(orcc_ok, orcc_ok$key)
rd_map <- split(rd_ok, rd_ok$key)

rows <- list()

for (i in seq_len(nrow(rec))) {
  row <- rec[i, ]
  key <- paste(row$dataset_name, row$analysis_key, sep = "|")
  or_row <- orcc_map[[key]]
  rd_row <- rd_map[[key]]

  or_row <- if (!is.null(or_row)) or_row[1, ] else NULL
  rd_row <- if (!is.null(rd_row)) rd_row[1, ] else NULL

  recommended <- row$recommended_method
  chosen <- NULL
  if (recommended == "OR_cc") {
    chosen <- or_row
  } else if (recommended == "RD") {
    chosen <- rd_row
  }

  rows[[length(rows) + 1]] <- data.frame(
    dataset_name = row$dataset_name,
    analysis_key = row$analysis_key,
    analysis_number = row$analysis_number,
    subgroup_number = row$subgroup_number,
    analysis_name = row$analysis_name,
    k = row$k,
    double_zero = row$double_zero,
    sparse_events = row$sparse_events,
    recommended_method = recommended,
    recommendation_reason = row$recommendation_reason,
    chosen_pooled_effect = if (!is.null(chosen)) chosen$pooled_effect else NA,
    chosen_ci_lb = if (!is.null(chosen)) chosen$ci_lb else NA,
    chosen_ci_ub = if (!is.null(chosen)) chosen$ci_ub else NA,
    orcc_pooled_effect = if (!is.null(or_row)) or_row$pooled_effect else NA,
    orcc_ci_lb = if (!is.null(or_row)) or_row$ci_lb else NA,
    orcc_ci_ub = if (!is.null(or_row)) or_row$ci_ub else NA,
    rd_pooled_effect = if (!is.null(rd_row)) rd_row$pooled_effect else NA,
    rd_ci_lb = if (!is.null(rd_row)) rd_row$ci_lb else NA,
    rd_ci_ub = if (!is.null(rd_row)) rd_row$ci_ub else NA,
    stringsAsFactors = FALSE
  )
}

out_df <- do.call(rbind, rows)
write.csv(out_df, out_path, row.names = FALSE)

missing_chosen <- sum(is.na(out_df$chosen_pooled_effect))
report_lines <- c(
  "# OR_cc vs RD consolidated recommendations",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Summary",
  paste0("- Total rows: ", nrow(out_df)),
  paste0("- Missing chosen effect: ", missing_chosen),
  "",
  "## Output files",
  "- orcc_vs_rd_consolidated_recommendations.csv",
  "- orcc_vs_rd_consolidated_report.md"
)
writeLines(report_lines, report_path)

cat("Consolidated recommendations written to: ", report_path, "\n", sep = "")
