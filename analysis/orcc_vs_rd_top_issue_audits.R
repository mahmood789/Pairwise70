#!/usr/bin/env Rscript

standardize_names <- function(nms) {
  nms <- gsub("[^A-Za-z0-9]+", ".", nms)
  nms <- gsub("\\.+", ".", nms)
  nms <- gsub("^\\.|\\.$", "", nms)
  nms
}

ascii_clean <- function(x) {
  ifelse(is.na(x), "", iconv(x, to = "ASCII//TRANSLIT", sub = ""))
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

count_nonint <- function(x) {
  sum(!is.na(x) & abs(x - round(x)) > 1e-6)
}

ci_mismatch_count <- function(mean_val, se_val, ci_start, ci_end, tol = 0.01) {
  ok <- !is.na(mean_val) & !is.na(se_val) & !is.na(ci_start) & !is.na(ci_end)
  if (!any(ok)) {
    return(rep(FALSE, length(mean_val)))
  }
  calc_start <- mean_val - 1.96 * se_val
  calc_end <- mean_val + 1.96 * se_val
  ok & (abs(calc_start - ci_start) > tol | abs(calc_end - ci_end) > tol)
}

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[1] else file.path("analysis", "output")
data_dir <- if (length(args) >= 2) args[2] else file.path("data")
top_n <- if (length(args) >= 3) suppressWarnings(as.numeric(args[3])) else 20
min_paired <- if (length(args) >= 4) suppressWarnings(as.numeric(args[4])) else 5
if (is.na(top_n)) top_n <- 20
if (is.na(min_paired) || min_paired < 1) min_paired <- 1

report_summary_path <- file.path(output_dir, "orcc_vs_rd_dataset_reports", "orcc_vs_rd_dataset_report_summary.csv")
if (!file.exists(report_summary_path)) {
  stop(paste0("Missing file: ", report_summary_path))
}

report_summary <- read.csv(report_summary_path, stringsAsFactors = FALSE)
report_summary$analyses <- suppressWarnings(as.numeric(report_summary$analyses))
report_summary$mismatch_rate <- suppressWarnings(as.numeric(report_summary$mismatch_rate))
report_summary <- report_summary[report_summary$analyses >= min_paired, ]
report_summary <- report_summary[order(-report_summary$mismatch_rate, -report_summary$analyses), ]
if (top_n <= 0 || top_n >= nrow(report_summary)) {
  targets <- report_summary$dataset_name
  top_n <- nrow(report_summary)
} else {
  targets <- head(report_summary$dataset_name, top_n)
}

out_dir <- file.path(output_dir, "orcc_vs_rd_top_issue_rows")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

issue_summary_list <- list()

for (ds in targets) {
  rda_path <- file.path(data_dir, paste0(ds, ".rda"))
  if (!file.exists(rda_path)) {
    warning(paste0("Missing dataset: ", rda_path))
    next
  }

  env <- new.env()
  load(rda_path, envir = env)
  obj_names <- ls(env)
  if (length(obj_names) == 0) next
  df <- env[[obj_names[1]]]
  if (!is.data.frame(df)) next

  names(df) <- standardize_names(names(df))

  n_rows <- nrow(df)
  analysis_number <- if ("Analysis.number" %in% names(df)) as.character(df$Analysis.number) else rep("1", n_rows)
  subgroup_number <- if ("Subgroup.number" %in% names(df)) as.character(df$Subgroup.number) else rep(NA, n_rows)
  subgroup_id <- ifelse(is.na(subgroup_number) | subgroup_number == "", "overall", subgroup_number)
  analysis_key <- paste(analysis_number, subgroup_id, sep = "::")

  exp_cases <- safe_num(if ("Experimental.cases" %in% names(df)) df$Experimental.cases else NA)
  exp_n <- safe_num(if ("Experimental.N" %in% names(df)) df$Experimental.N else NA)
  ctrl_cases <- safe_num(if ("Control.cases" %in% names(df)) df$Control.cases else NA)
  ctrl_n <- safe_num(if ("Control.N" %in% names(df)) df$Control.N else NA)

  exp_sd <- safe_num(if ("Experimental.SD" %in% names(df)) df$Experimental.SD else NA)
  ctrl_sd <- safe_num(if ("Control.SD" %in% names(df)) df$Control.SD else NA)

  mean_val <- safe_num(if ("Mean" %in% names(df)) df$Mean else NA)
  giv_mean <- safe_num(if ("GIV.Mean" %in% names(df)) df$GIV.Mean else NA)
  giv_se <- safe_num(if ("GIV.SE" %in% names(df)) df$GIV.SE else NA)
  var_val <- safe_num(if ("Variance" %in% names(df)) df$Variance else NA)
  ci_start <- safe_num(if ("CI.start" %in% names(df)) df$CI.start else NA)
  ci_end <- safe_num(if ("CI.end" %in% names(df)) df$CI.end else NA)

  mean_for_ci <- mean_val
  mean_for_ci[is.na(mean_for_ci)] <- giv_mean[is.na(mean_for_ci)]
  se_from_var <- ifelse(!is.na(var_val) & var_val >= 0, sqrt(var_val), NA_real_)
  se_for_ci <- giv_se
  se_for_ci[is.na(se_for_ci)] <- se_from_var[is.na(se_for_ci)]

  exp_cases_gt_n <- exp_cases > exp_n
  ctrl_cases_gt_n <- ctrl_cases > ctrl_n
  exp_cases_lt0 <- exp_cases < 0
  ctrl_cases_lt0 <- ctrl_cases < 0
  exp_n_le0 <- exp_n <= 0
  ctrl_n_le0 <- ctrl_n <= 0
  exp_cases_nonint <- abs(exp_cases - round(exp_cases)) > 1e-6
  ctrl_cases_nonint <- abs(ctrl_cases - round(ctrl_cases)) > 1e-6
  exp_n_nonint <- abs(exp_n - round(exp_n)) > 1e-6
  ctrl_n_nonint <- abs(ctrl_n - round(ctrl_n)) > 1e-6

  double_zero <- exp_cases == 0 & ctrl_cases == 0
  all_events <- exp_cases == exp_n & ctrl_cases == ctrl_n
  sparse_events <- (exp_cases + ctrl_cases) < 5

  exp_sd_le0 <- exp_sd <= 0
  ctrl_sd_le0 <- ctrl_sd <= 0

  variance_negative <- var_val < 0
  giv_se_negative <- giv_se < 0
  ci_inverted <- ci_start > ci_end
  ci_mismatch <- ci_mismatch_count(mean_for_ci, se_for_ci, ci_start, ci_end)

  any_issue <- exp_cases_gt_n | ctrl_cases_gt_n | exp_cases_lt0 | ctrl_cases_lt0 |
    exp_n_le0 | ctrl_n_le0 | exp_cases_nonint | ctrl_cases_nonint |
    exp_n_nonint | ctrl_n_nonint | double_zero | all_events | sparse_events |
    exp_sd_le0 | ctrl_sd_le0 | variance_negative | giv_se_negative |
    ci_inverted | ci_mismatch

  issue_summary_list[[length(issue_summary_list) + 1]] <- data.frame(
    dataset_name = ds,
    n_rows = n_rows,
    issue_rows = sum(any_issue, na.rm = TRUE),
    exp_cases_gt_n = sum(exp_cases_gt_n, na.rm = TRUE),
    ctrl_cases_gt_n = sum(ctrl_cases_gt_n, na.rm = TRUE),
    exp_cases_lt0 = sum(exp_cases_lt0, na.rm = TRUE),
    ctrl_cases_lt0 = sum(ctrl_cases_lt0, na.rm = TRUE),
    exp_n_le0 = sum(exp_n_le0, na.rm = TRUE),
    ctrl_n_le0 = sum(ctrl_n_le0, na.rm = TRUE),
    exp_cases_nonint = sum(exp_cases_nonint, na.rm = TRUE),
    ctrl_cases_nonint = sum(ctrl_cases_nonint, na.rm = TRUE),
    exp_n_nonint = sum(exp_n_nonint, na.rm = TRUE),
    ctrl_n_nonint = sum(ctrl_n_nonint, na.rm = TRUE),
    double_zero = sum(double_zero, na.rm = TRUE),
    all_events = sum(all_events, na.rm = TRUE),
    sparse_events = sum(sparse_events, na.rm = TRUE),
    exp_sd_le0 = sum(exp_sd_le0, na.rm = TRUE),
    ctrl_sd_le0 = sum(ctrl_sd_le0, na.rm = TRUE),
    variance_negative = sum(variance_negative, na.rm = TRUE),
    giv_se_negative = sum(giv_se_negative, na.rm = TRUE),
    ci_inverted = sum(ci_inverted, na.rm = TRUE),
    ci_mismatch = sum(ci_mismatch, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  if (any(any_issue, na.rm = TRUE)) {
    keep_cols <- c(
      "Study", "Study.year", "Analysis.name", "Analysis.number", "Subgroup.number",
      "Experimental.cases", "Experimental.N", "Control.cases", "Control.N",
      "Experimental.mean", "Experimental.SD", "Control.mean", "Control.SD",
      "GIV.Mean", "GIV.SE", "Mean", "CI.start", "CI.end", "Variance", "Weight"
    )
    keep_cols <- keep_cols[keep_cols %in% names(df)]
    out_df <- df[any_issue, keep_cols, drop = FALSE]
    out_df$analysis_key <- analysis_key[any_issue]
    out_df$dataset_name <- ds
    out_df$exp_cases_gt_n <- exp_cases_gt_n[any_issue]
    out_df$ctrl_cases_gt_n <- ctrl_cases_gt_n[any_issue]
    out_df$exp_cases_lt0 <- exp_cases_lt0[any_issue]
    out_df$ctrl_cases_lt0 <- ctrl_cases_lt0[any_issue]
    out_df$exp_n_le0 <- exp_n_le0[any_issue]
    out_df$ctrl_n_le0 <- ctrl_n_le0[any_issue]
    out_df$exp_cases_nonint <- exp_cases_nonint[any_issue]
    out_df$ctrl_cases_nonint <- ctrl_cases_nonint[any_issue]
    out_df$exp_n_nonint <- exp_n_nonint[any_issue]
    out_df$ctrl_n_nonint <- ctrl_n_nonint[any_issue]
    out_df$double_zero <- double_zero[any_issue]
    out_df$all_events <- all_events[any_issue]
    out_df$sparse_events <- sparse_events[any_issue]
    out_df$exp_sd_le0 <- exp_sd_le0[any_issue]
    out_df$ctrl_sd_le0 <- ctrl_sd_le0[any_issue]
    out_df$variance_negative <- variance_negative[any_issue]
    out_df$giv_se_negative <- giv_se_negative[any_issue]
    out_df$ci_inverted <- ci_inverted[any_issue]
    out_df$ci_mismatch <- ci_mismatch[any_issue]

    char_cols <- vapply(out_df, is.character, logical(1))
    out_df[char_cols] <- lapply(out_df[char_cols], ascii_clean)

    out_path <- file.path(out_dir, paste0(ds, "_issue_rows.csv"))
    write.csv(out_df, out_path, row.names = FALSE)
  }
}

issue_summary <- do.call(rbind, issue_summary_list)
issue_summary_path <- file.path(output_dir, "orcc_vs_rd_top_dataset_issue_audits.csv")
write.csv(issue_summary, issue_summary_path, row.names = FALSE)

report_path <- file.path(output_dir, "orcc_vs_rd_top_dataset_issue_report.md")
report_lines <- c(
  "# OR_cc vs RD top mismatch issue audit",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("- Top N datasets: ", top_n),
  paste0("- Min paired analyses: ", min_paired),
  paste0("- Datasets audited: ", length(targets)),
  "",
  "## Output files",
  "- orcc_vs_rd_top_dataset_issue_audits.csv",
  "- orcc_vs_rd_top_issue_rows/ (per-dataset flagged rows)"
)
writeLines(report_lines, report_path)

cat("Top mismatch issue audit complete.\n")
cat("Summary: ", issue_summary_path, "\n", sep = "")
cat("Rows folder: ", out_dir, "\n", sep = "")
