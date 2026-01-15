#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[1] else file.path("analysis", "output")
results_path <- file.path(output_dir, "remediation_analysis_results.csv")
out_path <- file.path(output_dir, "orcc_vs_rd_recommendations.csv")
report_path <- file.path(output_dir, "orcc_vs_rd_recommendations_report.md")

if (!file.exists(results_path)) {
  stop(paste0("Missing file: ", results_path))
}

df <- read.csv(results_path, stringsAsFactors = FALSE)
df <- df[df$method %in% c("OR_cc", "RD"), ]

df$key <- paste(df$dataset_name, df$analysis_key, sep = "|")
orcc <- df[df$method == "OR_cc", ]
rd <- df[df$method == "RD", ]

orcc$key <- paste(orcc$dataset_name, orcc$analysis_key, sep = "|")
rd$key <- paste(rd$dataset_name, rd$analysis_key, sep = "|")

all_keys <- unique(c(orcc$key, rd$key))

orcc_map <- split(orcc, orcc$key)
rd_map <- split(rd, rd$key)

rows <- list()

for (key in all_keys) {
  or_row <- orcc_map[[key]]
  rd_row <- rd_map[[key]]

  or_ok <- !is.null(or_row) && any(or_row$meta_status == "ok")
  rd_ok <- !is.null(rd_row) && any(rd_row$meta_status == "ok")

  if (or_ok) {
    or_row <- or_row[or_row$meta_status == "ok", ][1, ]
  } else if (!is.null(or_row)) {
    or_row <- or_row[1, ]
  }

  if (rd_ok) {
    rd_row <- rd_row[rd_row$meta_status == "ok", ][1, ]
  } else if (!is.null(rd_row)) {
    rd_row <- rd_row[1, ]
  }

  dataset_name <- if (!is.null(or_row)) or_row$dataset_name else rd_row$dataset_name
  analysis_key <- if (!is.null(or_row)) or_row$analysis_key else rd_row$analysis_key
  analysis_number <- if (!is.null(or_row)) or_row$analysis_number else rd_row$analysis_number
  subgroup_number <- if (!is.null(or_row)) or_row$subgroup_number else rd_row$subgroup_number
  analysis_name <- if (!is.null(or_row)) or_row$analysis_name else rd_row$analysis_name
  k <- if (!is.null(or_row)) or_row$k else rd_row$k
  double_zero <- if (!is.null(or_row)) or_row$double_zero else rd_row$double_zero
  sparse_events <- if (!is.null(or_row)) or_row$sparse_events else rd_row$sparse_events

  prefer <- "OR_cc"
  reason <- "default_or"
  if (!is.na(double_zero) && double_zero > 0) {
    prefer <- "RD"
    reason <- "double_zero"
  } else if (!is.na(sparse_events) && sparse_events > 0) {
    prefer <- "RD"
    reason <- "sparse_events"
  }

  recommended <- prefer
  if (prefer == "OR_cc" && !or_ok && rd_ok) {
    recommended <- "RD"
    reason <- "fallback_rd"
  } else if (prefer == "RD" && !rd_ok && or_ok) {
    recommended <- "OR_cc"
    reason <- "fallback_or"
  } else if (!or_ok && !rd_ok) {
    recommended <- "none"
    reason <- "no_ok_method"
  }

  rows[[length(rows) + 1]] <- data.frame(
    dataset_name = dataset_name,
    analysis_key = analysis_key,
    analysis_number = analysis_number,
    subgroup_number = subgroup_number,
    analysis_name = analysis_name,
    k = k,
    double_zero = double_zero,
    sparse_events = sparse_events,
    orcc_ok = or_ok,
    rd_ok = rd_ok,
    recommended_method = recommended,
    recommendation_reason = reason,
    orcc_pooled_effect = if (or_ok) or_row$pooled_effect else NA,
    orcc_ci_lb = if (or_ok) or_row$ci_lb else NA,
    orcc_ci_ub = if (or_ok) or_row$ci_ub else NA,
    rd_pooled_effect = if (rd_ok) rd_row$pooled_effect else NA,
    rd_ci_lb = if (rd_ok) rd_row$ci_lb else NA,
    rd_ci_ub = if (rd_ok) rd_row$ci_ub else NA,
    stringsAsFactors = FALSE
  )
}

out_df <- do.call(rbind, rows)
write.csv(out_df, out_path, row.names = FALSE)

method_counts <- table(out_df$recommended_method)
reason_counts <- table(out_df$recommendation_reason)

report_lines <- c(
  "# OR_cc vs RD recommendations",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Recommendation counts",
  paste0("- ", paste(names(method_counts), as.integer(method_counts), sep = ": ", collapse = "\n- ")),
  "",
  "## Recommendation reasons",
  paste0("- ", paste(names(reason_counts), as.integer(reason_counts), sep = ": ", collapse = "\n- ")),
  "",
  "## Output files",
  "- orcc_vs_rd_recommendations.csv",
  "- orcc_vs_rd_recommendations_report.md"
)
writeLines(report_lines, report_path)

cat("Recommendations written to: ", report_path, "\n", sep = "")
