#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[1] else file.path("analysis", "output")
data_dir <- if (length(args) >= 2) args[2] else file.path("data")

run_script <- function(script, ...) {
  cmd <- c(script, ...)
  status <- system2("Rscript", cmd, stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(status, "status"))) {
    stop(paste0("Script failed: ", script))
  }
}

scripts <- c(
  "analysis/orcc_vs_rd_summary.R",
  "analysis/orcc_vs_rd_disagreements.R",
  "analysis/orcc_vs_rd_further_analysis.R",
  "analysis/orcc_vs_rd_extract_rows.R",
  "analysis/orcc_vs_rd_dataset_reports.R",
  "analysis/orcc_vs_rd_issue_rankings.R",
  "analysis/orcc_vs_rd_actionable_fixes.R",
  "analysis/orcc_vs_rd_recommendations.R",
  "analysis/orcc_vs_rd_consolidated_table.R",
  "analysis/orcc_vs_rd_fix_recompute.R",
  "analysis/orcc_vs_rd_fix_recompute_deltas.R",
  "analysis/orcc_vs_rd_robma_rve.R"
)

cat("Running OR_cc vs RD pipeline...\n")
cat("Output dir: ", output_dir, "\n", sep = "")
cat("Data dir: ", data_dir, "\n\n", sep = "")

run_script("analysis/orcc_vs_rd_summary.R", output_dir)
run_script("analysis/orcc_vs_rd_disagreements.R", output_dir)
run_script("analysis/orcc_vs_rd_further_analysis.R", output_dir)
run_script("analysis/orcc_vs_rd_extract_rows.R", output_dir, data_dir, "0", "1")
run_script("analysis/orcc_vs_rd_dataset_reports.R", output_dir)
run_script("analysis/orcc_vs_rd_issue_rankings.R", output_dir, "0")
run_script("analysis/orcc_vs_rd_actionable_fixes.R", output_dir)
run_script("analysis/orcc_vs_rd_recommendations.R", output_dir)
run_script("analysis/orcc_vs_rd_consolidated_table.R", output_dir)
run_script("analysis/orcc_vs_rd_fix_recompute.R", data_dir, output_dir)
run_script("analysis/orcc_vs_rd_fix_recompute_deltas.R", output_dir)
run_script("analysis/orcc_vs_rd_robma_rve.R", data_dir, output_dir)

cat("\nPipeline complete.\n")
