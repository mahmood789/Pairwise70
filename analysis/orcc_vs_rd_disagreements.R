#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[1] else file.path("analysis", "output")
results_path <- file.path(output_dir, "remediation_analysis_results.csv")
disagree_path <- file.path(output_dir, "orcc_vs_rd_disagreements.csv")
report_path <- file.path(output_dir, "orcc_vs_rd_disagreement_report.md")

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
dir_mismatch <- sign_match == FALSE

or_sig <- !is.na(merged$ci_lb_or) & !is.na(merged$ci_ub_or) &
  ((merged$ci_lb_or > 0) | (merged$ci_ub_or < 0))
rd_sig <- !is.na(merged$ci_lb_rd) & !is.na(merged$ci_ub_rd) &
  ((merged$ci_lb_rd > 0) | (merged$ci_ub_rd < 0))
sig_mismatch <- or_sig != rd_sig
sig_mismatch[is.na(or_sig) | is.na(rd_sig)] <- NA

mismatch_type <- rep("none", nrow(merged))
mismatch_type[dir_mismatch & sig_mismatch] <- "direction_and_significance"
mismatch_type[dir_mismatch & !sig_mismatch] <- "direction_only"
mismatch_type[!dir_mismatch & sig_mismatch] <- "significance_only"
mismatch_type[is.na(dir_mismatch) | is.na(sig_mismatch)] <- "unknown"

sparse <- merged$sparse_events_or > 0
double_zero <- merged$double_zero_or > 0

k_vals <- as.numeric(merged$k_or)
k_bin <- cut(
  k_vals,
  breaks = c(-Inf, 4, 9, 19, Inf),
  labels = c("2-4", "5-9", "10-19", "20+"),
  right = TRUE
)

disagree_df <- data.frame(
  dataset_name = merged$dataset_name_or,
  analysis_key = merged$analysis_key_or,
  analysis_number = merged$analysis_number_or,
  subgroup_number = merged$subgroup_number_or,
  analysis_name = merged$analysis_name_or,
  k = k_vals,
  log_or = log_or,
  log_or_ci_lb = merged$ci_lb_or,
  log_or_ci_ub = merged$ci_ub_or,
  rd = rd_eff,
  rd_ci_lb = merged$ci_lb_rd,
  rd_ci_ub = merged$ci_ub_rd,
  or_sig = or_sig,
  rd_sig = rd_sig,
  direction_match = sign_match,
  mismatch_type = mismatch_type,
  sparse_events = sparse,
  double_zero = double_zero,
  stringsAsFactors = FALSE
)

write.csv(disagree_df, disagree_path, row.names = FALSE)

total <- nrow(disagree_df)
dir_mismatch_n <- sum(disagree_df$direction_match == FALSE, na.rm = TRUE)
sig_mismatch_n <- sum(disagree_df$mismatch_type %in% c("significance_only", "direction_and_significance"), na.rm = TRUE)
both_mismatch_n <- sum(disagree_df$mismatch_type == "direction_and_significance", na.rm = TRUE)
dir_only_n <- sum(disagree_df$mismatch_type == "direction_only", na.rm = TRUE)
sig_only_n <- sum(disagree_df$mismatch_type == "significance_only", na.rm = TRUE)
unknown_n <- sum(disagree_df$mismatch_type == "unknown", na.rm = TRUE)

rate <- function(x) {
  if (is.na(x) || total == 0) return(NA)
  round(x / total, 3)
}

dir_match_rate_sparse <- round(mean(disagree_df$direction_match[disagree_df$sparse_events], na.rm = TRUE), 3)
dir_match_rate_double_zero <- round(mean(disagree_df$direction_match[disagree_df$double_zero], na.rm = TRUE), 3)
sig_mismatch_rate_sparse <- round(mean(
  disagree_df$mismatch_type %in% c("significance_only", "direction_and_significance")
    [disagree_df$sparse_events],
  na.rm = TRUE
), 3)
sig_mismatch_rate_double_zero <- round(mean(
  disagree_df$mismatch_type %in% c("significance_only", "direction_and_significance")
    [disagree_df$double_zero],
  na.rm = TRUE
), 3)

k_table <- as.data.frame(table(k_bin, disagree_df$mismatch_type), stringsAsFactors = FALSE)

dir_disagree <- disagree_df[disagree_df$mismatch_type %in% c("direction_only", "direction_and_significance"), ]
dir_disagree <- dir_disagree[order(-abs(dir_disagree$log_or)), ]
top_dir <- head(dir_disagree, 15)

report_lines <- c(
  "# OR_cc vs RD disagreement analysis",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Summary",
  paste0("- Paired analyses: ", total),
  paste0("- Direction mismatches: ", dir_mismatch_n, " (", rate(dir_mismatch_n), ")"),
  paste0("- Significance mismatches: ", sig_mismatch_n, " (", rate(sig_mismatch_n), ")"),
  paste0("- Direction only: ", dir_only_n),
  paste0("- Significance only: ", sig_only_n),
  paste0("- Both direction and significance: ", both_mismatch_n),
  paste0("- Unknown: ", unknown_n),
  "",
  "## Stratified rates",
  paste0("- Direction match rate (sparse events): ", dir_match_rate_sparse),
  paste0("- Direction match rate (double-zero present): ", dir_match_rate_double_zero),
  paste0("- Significance mismatch rate (sparse events): ", sig_mismatch_rate_sparse),
  paste0("- Significance mismatch rate (double-zero present): ", sig_mismatch_rate_double_zero),
  "",
  "## Mismatch counts by k bin",
  paste0(
    "- ",
    paste(
      paste(k_table$k_bin, k_table$Var2, k_table$Freq, sep = ": "),
      collapse = "\n- "
    )
  ),
  "",
  "## Top direction disagreements (by |log OR|)",
  paste0(
    "- ",
    paste(
      paste(
        top_dir$dataset_name,
        top_dir$analysis_key,
        top_dir$analysis_name,
        "log_or=", round(top_dir$log_or, 3),
        "rd=", round(top_dir$rd, 3),
        "k=", top_dir$k,
        sep = " | "
      ),
      collapse = "\n- "
    )
  ),
  "",
  "## Output files",
  "- orcc_vs_rd_disagreements.csv",
  "- orcc_vs_rd_disagreement_report.md"
)

writeLines(report_lines, report_path)

cat("OR_cc vs RD disagreement outputs written to: ", report_path, "\n", sep = "")
