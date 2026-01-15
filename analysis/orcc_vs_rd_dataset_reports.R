#!/usr/bin/env Rscript

normalize_names <- function(nms) {
  nms <- tolower(nms)
  nms <- gsub("[^a-z0-9]+", ".", nms)
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

get_first <- function(df, col, default = NA) {
  if (!col %in% names(df)) {
    return(default)
  }
  val <- df[[col]]
  if (length(val) == 0) {
    return(default)
  }
  val <- val[!is.na(val)]
  if (length(val) == 0) {
    return(default)
  }
  val[1]
}

rate <- function(x, n) {
  if (is.na(n) || n == 0) return(NA_real_)
  round(x / n, 3)
}

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[1] else file.path("analysis", "output")
mismatch_dir <- if (length(args) >= 2) args[2] else file.path(output_dir, "orcc_vs_rd_mismatch_rows")
report_dir <- file.path(output_dir, "orcc_vs_rd_dataset_reports")
dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

files <- list.files(mismatch_dir, pattern = "_mismatch_rows\\.csv$", full.names = TRUE)
if (length(files) == 0) {
  stop(paste0("No mismatch files found in: ", mismatch_dir))
}

dataset_summary_list <- list()
analysis_summary_list <- list()

for (file in files) {
  df <- read.csv(file, stringsAsFactors = FALSE)
  if (nrow(df) == 0) next

  names(df) <- normalize_names(names(df))

  dataset_name <- unique(df$dataset.name)
  dataset_name <- dataset_name[!is.na(dataset_name)]
  if (length(dataset_name) == 0) {
    dataset_name <- NA_character_
  } else {
    dataset_name <- dataset_name[1]
  }

  if (is.na(dataset_name) || dataset_name == "") {
    dataset_name <- sub("_mismatch_rows\\.csv$", "", basename(file))
  }

  df$analysis.name <- ascii_clean(df$analysis.name)
  df$mismatch.type <- ifelse(is.na(df$mismatch.type), "unknown", df$mismatch.type)

  exp_cases <- safe_num(df$experimental.cases)
  exp_n <- safe_num(df$experimental.n)
  ctrl_cases <- safe_num(df$control.cases)
  ctrl_n <- safe_num(df$control.n)

  valid_binary <- !is.na(exp_cases) & !is.na(exp_n) & !is.na(ctrl_cases) & !is.na(ctrl_n) &
    exp_cases >= 0 & ctrl_cases >= 0 &
    exp_n > 0 & ctrl_n > 0 &
    exp_cases <= exp_n & ctrl_cases <= ctrl_n

  double_zero_rows <- sum(exp_cases == 0 & ctrl_cases == 0, na.rm = TRUE)
  all_events_rows <- sum(exp_cases == exp_n & ctrl_cases == ctrl_n, na.rm = TRUE)
  sparse_rows <- sum((exp_cases + ctrl_cases) < 5, na.rm = TRUE)
  any_zero_rows <- sum(exp_cases == 0 | ctrl_cases == 0, na.rm = TRUE)

  total_events <- sum(exp_cases[valid_binary] + ctrl_cases[valid_binary], na.rm = TRUE)
  total_n <- sum(exp_n[valid_binary] + ctrl_n[valid_binary], na.rm = TRUE)
  event_rate <- if (total_n > 0) round(total_events / total_n, 4) else NA_real_

  analysis_keys <- unique(df$analysis.key)
  per_analysis <- list()
  for (key in analysis_keys) {
    sub <- df[df$analysis.key == key, ]
    if (nrow(sub) == 0) next

    exp_cases_a <- safe_num(sub$experimental.cases)
    exp_n_a <- safe_num(sub$experimental.n)
    ctrl_cases_a <- safe_num(sub$control.cases)
    ctrl_n_a <- safe_num(sub$control.n)

    valid_a <- !is.na(exp_cases_a) & !is.na(exp_n_a) & !is.na(ctrl_cases_a) & !is.na(ctrl_n_a) &
      exp_cases_a >= 0 & ctrl_cases_a >= 0 &
      exp_n_a > 0 & ctrl_n_a > 0 &
      exp_cases_a <= exp_n_a & ctrl_cases_a <= ctrl_n_a

    total_events_a <- sum(exp_cases_a[valid_a] + ctrl_cases_a[valid_a], na.rm = TRUE)
    total_n_a <- sum(exp_n_a[valid_a] + ctrl_n_a[valid_a], na.rm = TRUE)
    event_rate_a <- if (total_n_a > 0) round(total_events_a / total_n_a, 4) else NA_real_

    analysis_summary <- data.frame(
      dataset_name = dataset_name,
      analysis_key = key,
      analysis_number = get_first(sub, "analysis.number", ""),
      subgroup_number = get_first(sub, "subgroup.number", ""),
      analysis_name = ascii_clean(get_first(sub, "analysis.name", "")),
      mismatch_type = get_first(sub, "mismatch.type", "unknown"),
      k_rows = nrow(sub),
      k_valid = sum(valid_a, na.rm = TRUE),
      invalid_rows = sum(!valid_a, na.rm = TRUE),
      double_zero_rows = sum(exp_cases_a == 0 & ctrl_cases_a == 0, na.rm = TRUE),
      all_events_rows = sum(exp_cases_a == exp_n_a & ctrl_cases_a == ctrl_n_a, na.rm = TRUE),
      sparse_rows = sum((exp_cases_a + ctrl_cases_a) < 5, na.rm = TRUE),
      any_zero_rows = sum(exp_cases_a == 0 | ctrl_cases_a == 0, na.rm = TRUE),
      event_rate = event_rate_a,
      log_or = safe_num(get_first(sub, "log.or", NA)),
      rd = safe_num(get_first(sub, "rd", NA)),
      or_sig = get_first(sub, "or.sig", NA),
      rd_sig = get_first(sub, "rd.sig", NA),
      direction_match = get_first(sub, "direction.match", NA),
      stringsAsFactors = FALSE
    )
    per_analysis[[length(per_analysis) + 1]] <- analysis_summary
  }

  analysis_df <- do.call(rbind, per_analysis)
  analysis_summary_list[[length(analysis_summary_list) + 1]] <- analysis_df

  total_analyses <- nrow(analysis_df)
  mismatch_counts <- table(analysis_df$mismatch_type)

  dir_mismatch_n <- sum(analysis_df$mismatch_type %in% c("direction_only", "direction_and_significance"))
  sig_mismatch_n <- sum(analysis_df$mismatch_type %in% c("significance_only", "direction_and_significance"))
  any_mismatch_n <- sum(analysis_df$mismatch_type %in% c("direction_only", "significance_only", "direction_and_significance"))

  dataset_summary_list[[length(dataset_summary_list) + 1]] <- data.frame(
    dataset_name = dataset_name,
    analyses = total_analyses,
    any_mismatch = any_mismatch_n,
    direction_mismatch = dir_mismatch_n,
    significance_mismatch = sig_mismatch_n,
    mismatch_rate = rate(any_mismatch_n, total_analyses),
    double_zero_rows = double_zero_rows,
    sparse_rows = sparse_rows,
    any_zero_rows = any_zero_rows,
    all_events_rows = all_events_rows,
    invalid_binary_rows = sum(!valid_binary, na.rm = TRUE),
    event_rate = event_rate,
    stringsAsFactors = FALSE
  )

  top_dir <- analysis_df[analysis_df$mismatch_type %in% c("direction_only", "direction_and_significance"), ]
  top_dir <- top_dir[order(-abs(top_dir$log_or)), ]
  top_dir <- head(top_dir, 10)

  report_lines <- c(
    paste0("# Dataset report: ", dataset_name),
    "",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "## Mismatch summary",
    paste0("- Analyses with mismatches: ", total_analyses),
    paste0("- Any mismatch: ", any_mismatch_n, " (", rate(any_mismatch_n, total_analyses), ")"),
    paste0("- Direction mismatch: ", dir_mismatch_n, " (", rate(dir_mismatch_n, total_analyses), ")"),
    paste0("- Significance mismatch: ", sig_mismatch_n, " (", rate(sig_mismatch_n, total_analyses), ")"),
    "",
    "## Row-level data signals (mismatch rows only)",
    paste0("- Double-zero rows: ", double_zero_rows),
    paste0("- All-event rows: ", all_events_rows),
    paste0("- Sparse-event rows (<5 total): ", sparse_rows),
    paste0("- Any zero rows: ", any_zero_rows),
    paste0("- Invalid binary rows: ", sum(!valid_binary, na.rm = TRUE)),
    paste0("- Event rate (valid rows): ", event_rate),
    "",
    "## Mismatch types",
    paste0(
      "- ",
      paste(paste(names(mismatch_counts), as.integer(mismatch_counts), sep = ": "), collapse = "\n- ")
    ),
    "",
    "## Top direction mismatches (by |log OR|)",
    paste0(
      "- ",
      paste(
        paste(
          top_dir$analysis_key,
          top_dir$analysis_name,
          paste0("log_or=", round(top_dir$log_or, 3)),
          paste0("rd=", round(top_dir$rd, 3)),
          paste0("k=", top_dir$k_rows),
          sep = " | "
        ),
        collapse = "\n- "
      )
    )
  )

  report_path <- file.path(report_dir, paste0(dataset_name, "_report.md"))
  writeLines(report_lines, report_path)
}

dataset_summary <- do.call(rbind, dataset_summary_list)
analysis_summary <- do.call(rbind, analysis_summary_list)

dataset_summary_path <- file.path(report_dir, "orcc_vs_rd_dataset_report_summary.csv")
analysis_summary_path <- file.path(report_dir, "orcc_vs_rd_analysis_report_summary.csv")
write.csv(dataset_summary, dataset_summary_path, row.names = FALSE)
write.csv(analysis_summary, analysis_summary_path, row.names = FALSE)

index_path <- file.path(report_dir, "orcc_vs_rd_dataset_report_index.md")
index_lines <- c(
  "# OR_cc vs RD dataset reports",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Summary table",
  paste0("- Summary CSV: ", dataset_summary_path),
  paste0("- Analysis CSV: ", analysis_summary_path),
  "",
  "## Reports",
  paste0("- ", paste(paste0(dataset_summary$dataset_name, "_report.md"), collapse = "\n- "))
)
writeLines(index_lines, index_path)

cat("Dataset reports written to: ", report_dir, "\n", sep = "")
