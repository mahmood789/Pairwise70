#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[1] else file.path("analysis", "output")
results_path <- file.path(output_dir, "remediation_analysis_results.csv")
report_path <- file.path(output_dir, "orcc_vs_rd_report.md")
summary_path <- file.path(output_dir, "orcc_vs_rd_summary.csv")

if (!file.exists(results_path)) {
  stop(paste0("Missing results file: ", results_path))
}

df <- read.csv(results_path, stringsAsFactors = FALSE)
df <- df[df$method %in% c("OR_cc", "RD") & df$meta_status == "ok", ]
df$key <- paste(df$dataset_name, df$analysis_key, sep = "|")

orcc <- df[df$method == "OR_cc", ]
rd <- df[df$method == "RD", ]
merged <- merge(orcc, rd, by = "key", suffixes = c("_or", "_rd"))

log_or <- as.numeric(merged$pooled_effect_or)
rd_eff <- as.numeric(merged$pooled_effect_rd)

sign_match <- sign(log_or) == sign(rd_eff)
sign_match[log_or == 0 | rd_eff == 0] <- NA

spearman <- suppressWarnings(cor(log_or, rd_eff, use = "complete.obs", method = "spearman"))

or_sig <- (merged$ci_lb_or > 0) | (merged$ci_ub_or < 0)
rd_sig <- (merged$ci_lb_rd > 0) | (merged$ci_ub_rd < 0)

sparse <- merged$sparse_events_or > 0
double_zero <- merged$double_zero_or > 0

summary_df <- data.frame(
  paired_analyses = nrow(merged),
  sign_match_rate = round(mean(sign_match, na.rm = TRUE), 3),
  spearman_log_or_vs_rd = round(spearman, 3),
  both_sig = sum(or_sig & rd_sig, na.rm = TRUE),
  or_sig_only = sum(or_sig & !rd_sig, na.rm = TRUE),
  rd_sig_only = sum(!or_sig & rd_sig, na.rm = TRUE),
  neither_sig = sum(!or_sig & !rd_sig, na.rm = TRUE),
  sparse_sign_match_rate = round(mean(sign_match[sparse], na.rm = TRUE), 3),
  double_zero_sign_match_rate = round(mean(sign_match[double_zero], na.rm = TRUE), 3),
  stringsAsFactors = FALSE
)

write.csv(summary_df, summary_path, row.names = FALSE)

report_lines <- c(
  "# OR_cc vs RD comparison",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Summary",
  paste0("- Paired analyses: ", summary_df$paired_analyses),
  paste0("- Direction agreement (sign match): ", summary_df$sign_match_rate),
  paste0("- Spearman correlation (log OR vs RD): ", summary_df$spearman_log_or_vs_rd),
  paste0("- Both significant: ", summary_df$both_sig),
  paste0("- OR_cc significant only: ", summary_df$or_sig_only),
  paste0("- RD significant only: ", summary_df$rd_sig_only),
  paste0("- Neither significant: ", summary_df$neither_sig),
  paste0("- Direction agreement (sparse events): ", summary_df$sparse_sign_match_rate),
  paste0("- Direction agreement (double-zero present): ", summary_df$double_zero_sign_match_rate),
  "",
  "## Notes",
  "- OR_cc pooled effects are on the log odds ratio scale.",
  "- RD pooled effects are on the absolute risk difference scale.",
  "- Magnitudes are not directly comparable; use direction and significance checks for contrast.",
  "",
  "## Output files",
  "- orcc_vs_rd_summary.csv",
  "- orcc_vs_rd_report.md"
)

writeLines(report_lines, report_path)

cat("OR_cc vs RD summary written to: ", report_path, "\n", sep = "")
