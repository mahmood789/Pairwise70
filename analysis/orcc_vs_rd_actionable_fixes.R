#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[1] else file.path("analysis", "output")
issue_path <- file.path(output_dir, "orcc_vs_rd_issue_rankings.csv")
report_path <- file.path(output_dir, "orcc_vs_rd_actionable_fixes_report.md")
summary_path <- file.path(output_dir, "orcc_vs_rd_actionable_fixes.csv")

if (!file.exists(issue_path)) {
  stop(paste0("Missing file: ", issue_path))
}

issues <- read.csv(issue_path, stringsAsFactors = FALSE)

fix_map <- list(
  double_zero = "Prefer RD or GLMM; exclude double-zero from OR; report sensitivity.",
  sparse_events = "Prefer RD/GLMM; consider exact methods; report sensitivity.",
  all_events = "Check event coding; consider RD; use continuity correction if OR/RR.",
  exp_sd_le0 = "Verify SD extraction; impute SD within analysis or from control/median.",
  ctrl_sd_le0 = "Verify SD extraction; impute SD within analysis or from experimental/median.",
  ci_mismatch = "Re-derive SE from CI or variance; check log scale vs raw.",
  ci_inverted = "Swap CI bounds or correct source extraction.",
  variance_negative = "Fix variance calculation; re-derive from SE or CI.",
  giv_se_negative = "Fix SE calculation; re-derive from CI or variance.",
  exp_cases_gt_n = "Correct event counts; cases must be <= N.",
  ctrl_cases_gt_n = "Correct event counts; cases must be <= N.",
  exp_cases_lt0 = "Correct negative event counts.",
  ctrl_cases_lt0 = "Correct negative event counts.",
  exp_n_le0 = "Correct sample size; N must be > 0.",
  ctrl_n_le0 = "Correct sample size; N must be > 0.",
  exp_cases_nonint = "Review extraction/rounding; events should be integers.",
  ctrl_cases_nonint = "Review extraction/rounding; events should be integers.",
  exp_n_nonint = "Review extraction/rounding; sample sizes should be integers.",
  ctrl_n_nonint = "Review extraction/rounding; sample sizes should be integers.",
  duplicate_study = "De-duplicate within analysis (Study + analysis/subgroup).",
  missing_study = "Backfill missing study identifiers before deduplication.",
  weight_negative = "Recompute weights from effect sizes/variances."
)

get_fix <- function(issue) {
  if (issue %in% names(fix_map)) return(fix_map[[issue]])
  "Review data extraction and analysis assumptions."
}

fix_rows <- list()

for (i in seq_len(nrow(issues))) {
  row <- issues[i, ]
  issue_names <- c(row$issue1, row$issue2, row$issue3, row$issue4, row$issue5)
  issue_counts <- c(row$issue1_count, row$issue2_count, row$issue3_count, row$issue4_count, row$issue5_count)

  fixes <- character(0)
  for (j in seq_along(issue_names)) {
    issue <- issue_names[j]
    count <- issue_counts[j]
    if (is.na(issue) || issue == "" || count <= 0) next
    fixes <- c(fixes, paste0(issue, ": ", get_fix(issue)))
  }

  fix_rows[[length(fix_rows) + 1]] <- data.frame(
    dataset_name = row$dataset_name,
    mismatch_rate = row$mismatch_rate,
    paired_analyses = row$paired_analyses,
    fix_1 = ifelse(length(fixes) >= 1, fixes[1], ""),
    fix_2 = ifelse(length(fixes) >= 2, fixes[2], ""),
    fix_3 = ifelse(length(fixes) >= 3, fixes[3], ""),
    fix_4 = ifelse(length(fixes) >= 4, fixes[4], ""),
    fix_5 = ifelse(length(fixes) >= 5, fixes[5], ""),
    stringsAsFactors = FALSE
  )
}

fix_df <- do.call(rbind, fix_rows)
write.csv(fix_df, summary_path, row.names = FALSE)

top <- head(fix_df, 20)
report_lines <- c(
  "# OR_cc vs RD actionable fixes",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Top 20 datasets by mismatch rate",
  paste0(
    "- ",
    paste(
      paste(
        top$dataset_name,
        "mismatch_rate=", top$mismatch_rate,
        "fix_1=", top$fix_1,
        sep = " "
      ),
      collapse = "\n- "
    )
  ),
  "",
  "## Output files",
  "- orcc_vs_rd_actionable_fixes.csv",
  "- orcc_vs_rd_actionable_fixes_report.md"
)
writeLines(report_lines, report_path)

cat("Actionable fixes written to: ", report_path, "\n", sep = "")
