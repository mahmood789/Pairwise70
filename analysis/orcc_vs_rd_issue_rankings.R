#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[1] else file.path("analysis", "output")
top_n <- if (length(args) >= 2) suppressWarnings(as.numeric(args[2])) else 20
if (is.na(top_n) || top_n <= 0) top_n <- 20

issue_path <- file.path(output_dir, "dataset_issue_counts.csv")
inventory_path <- file.path(output_dir, "dataset_inventory.csv")
report_summary_path <- file.path(output_dir, "orcc_vs_rd_dataset_reports", "orcc_vs_rd_dataset_report_summary.csv")

if (!file.exists(issue_path)) stop(paste0("Missing file: ", issue_path))
if (!file.exists(inventory_path)) stop(paste0("Missing file: ", inventory_path))
if (!file.exists(report_summary_path)) stop(paste0("Missing file: ", report_summary_path))

issues <- read.csv(issue_path, stringsAsFactors = FALSE)
inventory <- read.csv(inventory_path, stringsAsFactors = FALSE)
report_summary <- read.csv(report_summary_path, stringsAsFactors = FALSE)

issues <- merge(issues, inventory[, c("dataset_name", "n_rows")], by = "dataset_name", all.x = TRUE)
report_summary <- report_summary[, c("dataset_name", "analyses", "mismatch_rate")]
names(report_summary)[2] <- "paired_analyses"
issues <- merge(issues, report_summary, by = "dataset_name", all.x = TRUE)

issue_cols <- setdiff(names(issues), c("dataset_name", "n_rows", "paired_analyses", "mismatch_rate"))

ranked_list <- list()

for (i in seq_len(nrow(issues))) {
  row <- issues[i, ]
  counts <- as.numeric(row[issue_cols])
  names(counts) <- issue_cols
  counts <- counts[!is.na(counts)]
  if (length(counts) == 0) next

  ord <- order(-counts)
  top_issues <- counts[ord]
  top_issues <- top_issues[top_issues > 0]

  if (length(top_issues) == 0) {
    ranked_list[[length(ranked_list) + 1]] <- data.frame(
      dataset_name = row$dataset_name,
      paired_analyses = row$paired_analyses,
      mismatch_rate = row$mismatch_rate,
      issue1 = "",
      issue1_count = 0,
      issue1_rate = NA,
      issue2 = "",
      issue2_count = 0,
      issue2_rate = NA,
      issue3 = "",
      issue3_count = 0,
      issue3_rate = NA,
      issue4 = "",
      issue4_count = 0,
      issue4_rate = NA,
      issue5 = "",
      issue5_count = 0,
      issue5_rate = NA,
      stringsAsFactors = FALSE
    )
    next
  }

  n_rows <- ifelse(is.na(row$n_rows), NA, row$n_rows)
  top_n_issues <- head(top_issues, 5)
  issue_names <- names(top_n_issues)
  issue_counts <- as.integer(top_n_issues)
  issue_rates <- if (!is.na(n_rows) && n_rows > 0) round(issue_counts / n_rows, 3) else rep(NA, length(issue_counts))

  padded_names <- c(issue_names, rep("", 5 - length(issue_names)))
  padded_counts <- c(issue_counts, rep(0, 5 - length(issue_counts)))
  padded_rates <- c(issue_rates, rep(NA, 5 - length(issue_rates)))

  ranked_list[[length(ranked_list) + 1]] <- data.frame(
    dataset_name = row$dataset_name,
    paired_analyses = row$paired_analyses,
    mismatch_rate = row$mismatch_rate,
    issue1 = padded_names[1],
    issue1_count = padded_counts[1],
    issue1_rate = padded_rates[1],
    issue2 = padded_names[2],
    issue2_count = padded_counts[2],
    issue2_rate = padded_rates[2],
    issue3 = padded_names[3],
    issue3_count = padded_counts[3],
    issue3_rate = padded_rates[3],
    issue4 = padded_names[4],
    issue4_count = padded_counts[4],
    issue4_rate = padded_rates[4],
    issue5 = padded_names[5],
    issue5_count = padded_counts[5],
    issue5_rate = padded_rates[5],
    stringsAsFactors = FALSE
  )
}

ranked_df <- do.call(rbind, ranked_list)
ranked_df <- ranked_df[order(-ranked_df$mismatch_rate, -ranked_df$paired_analyses), ]

ranked_path <- file.path(output_dir, "orcc_vs_rd_issue_rankings.csv")
write.csv(ranked_df, ranked_path, row.names = FALSE)

top_ranked <- head(ranked_df, top_n)

report_path <- file.path(output_dir, "orcc_vs_rd_issue_rankings_report.md")
report_lines <- c(
  "# OR_cc vs RD issue rankings",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  paste0("## Top ", top_n, " datasets by mismatch rate"),
  paste0(
    "- ",
    paste(
      paste(
        top_ranked$dataset_name,
        "mismatch_rate=", top_ranked$mismatch_rate,
        "paired=", top_ranked$paired_analyses,
        "issue1=", top_ranked$issue1,
        "count1=", top_ranked$issue1_count,
        sep = " "
      ),
      collapse = "\n- "
    )
  ),
  "",
  "## Output files",
  "- orcc_vs_rd_issue_rankings.csv",
  "- orcc_vs_rd_issue_rankings_report.md"
)
writeLines(report_lines, report_path)

cat("Issue rankings written to: ", report_path, "\n", sep = "")
