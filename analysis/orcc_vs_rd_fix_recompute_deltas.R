#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[1] else file.path("analysis", "output")

orig_path <- file.path(output_dir, "orcc_vs_rd_consolidated_recommendations.csv")
fix_path <- file.path(output_dir, "orcc_vs_rd_fix_recompute_recommendations.csv")
compare_path <- file.path(output_dir, "orcc_vs_rd_fix_recompute_comparison.csv")
out_path <- file.path(output_dir, "orcc_vs_rd_fix_recompute_deltas.csv")
report_path <- file.path(output_dir, "orcc_vs_rd_fix_recompute_deltas_report.md")

if (!file.exists(orig_path)) stop(paste0("Missing file: ", orig_path))
if (!file.exists(fix_path)) stop(paste0("Missing file: ", fix_path))
if (!file.exists(compare_path)) stop(paste0("Missing file: ", compare_path))

orig <- read.csv(orig_path, stringsAsFactors = FALSE)
fix <- read.csv(fix_path, stringsAsFactors = FALSE)
comp <- read.csv(compare_path, stringsAsFactors = FALSE)

delta <- data.frame(
  dataset_name = comp$dataset_name_orig,
  analysis_key = comp$analysis_key_orig,
  analysis_number = comp$analysis_number_orig,
  subgroup_number = comp$subgroup_number_orig,
  analysis_name = comp$analysis_name_orig,
  method_orig = comp$recommended_method_orig,
  method_fix = comp$recommended_method_fix,
  method_changed = comp$method_changed,
  chosen_effect_orig = comp$orcc_pooled_effect,
  chosen_effect_fix = comp$orcc_pooled_effect,
  chosen_ci_lb_orig = comp$orcc_ci_lb,
  chosen_ci_ub_orig = comp$orcc_ci_ub,
  chosen_ci_lb_fix = comp$orcc_ci_lb,
  chosen_ci_ub_fix = comp$orcc_ci_ub,
  sparse_events = comp$sparse_events_orig,
  double_zero = comp$double_zero_orig,
  stringsAsFactors = FALSE
)

delta$chosen_effect_orig <- ifelse(
  delta$method_orig == "OR_cc",
  comp$orcc_pooled_effect,
  comp$rd_pooled_effect
)
delta$chosen_ci_lb_orig <- ifelse(
  delta$method_orig == "OR_cc",
  comp$orcc_ci_lb,
  comp$rd_ci_lb
)
delta$chosen_ci_ub_orig <- ifelse(
  delta$method_orig == "OR_cc",
  comp$orcc_ci_ub,
  comp$rd_ci_ub
)

delta$chosen_effect_fix <- ifelse(
  delta$method_fix == "OR_cc",
  comp$orcc_pooled_effect,
  comp$rd_pooled_effect
)
delta$chosen_ci_lb_fix <- ifelse(
  delta$method_fix == "OR_cc",
  comp$orcc_ci_lb,
  comp$rd_ci_lb
)
delta$chosen_ci_ub_fix <- ifelse(
  delta$method_fix == "OR_cc",
  comp$orcc_ci_ub,
  comp$rd_ci_ub
)

delta$effect_delta <- suppressWarnings(as.numeric(delta$chosen_effect_fix) - as.numeric(delta$chosen_effect_orig))

write.csv(delta, out_path, row.names = FALSE)

changes_by_dataset <- aggregate(method_changed ~ dataset_name, data = delta, FUN = sum)
names(changes_by_dataset)[2] <- "method_changes"
total_by_dataset <- aggregate(method_changed ~ dataset_name, data = delta, FUN = length)
names(total_by_dataset)[2] <- "total_analyses"
changes_by_dataset <- merge(changes_by_dataset, total_by_dataset, by = "dataset_name")
changes_by_dataset$change_rate <- round(changes_by_dataset$method_changes / changes_by_dataset$total_analyses, 3)
changes_by_dataset <- changes_by_dataset[order(-changes_by_dataset$method_changes, -changes_by_dataset$change_rate), ]

changes_path <- file.path(output_dir, "orcc_vs_rd_fix_recompute_changes_by_dataset.csv")
write.csv(changes_by_dataset, changes_path, row.names = FALSE)

top_changes <- head(changes_by_dataset, 20)
method_counts <- table(delta$method_fix)
changed_n <- sum(delta$method_changed, na.rm = TRUE)

report_lines <- c(
  "# OR_cc vs RD fix-recompute deltas",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Summary",
  paste0("- Total analyses: ", nrow(delta)),
  paste0("- Methods changed: ", changed_n),
  paste0("- Method counts (fix): ", paste(names(method_counts), as.integer(method_counts), sep = ": ", collapse = ", ")),
  "",
  "## Top datasets by method changes",
  paste0(
    "- ",
    paste(
      paste(
        top_changes$dataset_name,
        "changes=", top_changes$method_changes,
        "rate=", top_changes$change_rate,
        sep = " "
      ),
      collapse = "\n- "
    )
  ),
  "",
  "## Output files",
  "- orcc_vs_rd_fix_recompute_deltas.csv",
  "- orcc_vs_rd_fix_recompute_changes_by_dataset.csv",
  "- orcc_vs_rd_fix_recompute_deltas_report.md"
)
writeLines(report_lines, report_path)

cat("Fix-recompute deltas written to: ", report_path, "\n", sep = "")
