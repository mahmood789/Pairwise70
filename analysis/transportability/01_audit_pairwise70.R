# Audit Pairwise70 datasets for transportability groundwork
root <- normalizePath("C:/Users/user/OneDrive - NHS/Documents/Pairwise70", winslash = "/", mustWork = TRUE)
output_dir <- file.path(root, "analysis", "transportability")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

data_dir <- file.path(root, "data")
files <- list.files(data_dir, pattern = "\\.rda$", full.names = TRUE)

safe_load <- function(path) {
  env <- new.env(parent = emptyenv())
  nm <- load(path, envir = env)
  if (length(nm) == 0) return(NULL)
  obj <- env[[nm[1]]]
  list(name = nm[1], data = obj)
}

pct_missing <- function(x) {
  if (length(x) == 0) return(NA_real_)
  x_chr <- as.character(x)
  mean(is.na(x_chr) | trimws(x_chr) == "")
}

rows <- list()
for (f in files) {
  entry <- safe_load(f)
  if (is.null(entry)) next
  df <- entry$data
  if (!is.data.frame(df)) next

  cols <- names(df)
  has_binary_cols <- all(c("Experimental.cases", "Experimental.N", "Control.cases", "Control.N") %in% cols)
  has_cont_cols <- all(c("Experimental.mean", "Experimental.SD", "Experimental.N",
                         "Control.mean", "Control.SD", "Control.N") %in% cols)

  binary_present <- FALSE
  if (has_binary_cols) {
    binary_vals <- c(df$Experimental.cases, df$Experimental.N, df$Control.cases, df$Control.N)
    binary_present <- any(is.finite(suppressWarnings(as.numeric(binary_vals))))
  }

  cont_present <- FALSE
  if (has_cont_cols) {
    cont_vals <- c(df$Experimental.mean, df$Experimental.SD, df$Control.mean, df$Control.SD)
    cont_present <- any(is.finite(suppressWarnings(as.numeric(cont_vals))))
  }

  n_rows <- nrow(df)
  n_studies <- if ("Study" %in% cols) length(unique(df$Study)) else NA_integer_
  n_comparisons <- if ("Comparison" %in% cols) length(unique(df$Comparison)) else NA_integer_
  n_outcomes <- if ("Outcome" %in% cols) length(unique(df$Outcome)) else NA_integer_
  n_subgroups <- if ("Subgroup" %in% cols) length(unique(df$Subgroup)) else NA_integer_

  study_year <- if ("Study.year" %in% cols) suppressWarnings(as.numeric(df$Study.year)) else numeric(0)

  ctrl_risk <- numeric(0)
  if (binary_present) {
    ctrl_cases <- suppressWarnings(as.numeric(df$Control.cases))
    ctrl_n <- suppressWarnings(as.numeric(df$Control.N))
    ok <- is.finite(ctrl_cases) & is.finite(ctrl_n) & ctrl_n > 0
    ctrl_risk <- ctrl_cases[ok] / ctrl_n[ok]
  }

  ctrl_mean <- numeric(0)
  if (cont_present) {
    ctrl_mean <- suppressWarnings(as.numeric(df$Control.mean))
    ctrl_mean <- ctrl_mean[is.finite(ctrl_mean)]
  }

  outcome_type <- if (binary_present && cont_present) {
    "mixed"
  } else if (binary_present) {
    "binary"
  } else if (cont_present) {
    "continuous"
  } else {
    "unknown"
  }

  has_year <- length(study_year) > 0 && any(is.finite(study_year))

  rows[[length(rows) + 1]] <- data.frame(
    dataset = entry$name,
    file = basename(f),
    n_rows = n_rows,
    n_studies = n_studies,
    n_comparisons = n_comparisons,
    n_outcomes = n_outcomes,
    n_subgroups = n_subgroups,
    has_outcome_col = "Outcome" %in% cols,
    has_comparison_col = "Comparison" %in% cols,
    has_subgroup_col = "Subgroup" %in% cols,
    has_review_doi_col = "review_doi" %in% cols,
    has_review_url_col = "review_url" %in% cols,
    outcome_type = outcome_type,
    study_year_min = if (has_year) min(study_year, na.rm = TRUE) else NA_real_,
    study_year_max = if (has_year) max(study_year, na.rm = TRUE) else NA_real_,
    ctrl_risk_mean = if (length(ctrl_risk) > 0) mean(ctrl_risk) else NA_real_,
    ctrl_risk_median = if (length(ctrl_risk) > 0) median(ctrl_risk) else NA_real_,
    ctrl_mean_mean = if (length(ctrl_mean) > 0) mean(ctrl_mean) else NA_real_,
    missing_study_year = if ("Study.year" %in% cols) pct_missing(df$Study.year) else NA_real_,
    missing_outcome = if ("Outcome" %in% cols) pct_missing(df$Outcome) else NA_real_,
    missing_comparison = if ("Comparison" %in% cols) pct_missing(df$Comparison) else NA_real_,
    missing_review_doi = if ("review_doi" %in% cols) pct_missing(df$review_doi) else NA_real_,
    stringsAsFactors = FALSE
  )
}

audit <- do.call(rbind, rows)
audit_path <- file.path(output_dir, "pairwise70_audit.csv")
write.csv(audit, audit_path, row.names = FALSE)

# Aggregate summary
summary_path <- file.path(output_dir, "pairwise70_audit_summary.md")

outcome_counts <- sort(table(audit$outcome_type), decreasing = TRUE)

summary_lines <- c(
  "# Pairwise70 Audit Summary",
  "",
  sprintf("Datasets: %d", nrow(audit)),
  sprintf("Total rows: %d", sum(audit$n_rows, na.rm = TRUE)),
  "",
  "## Column coverage (datasets with column)",
  sprintf("- Outcome: %d", sum(audit$has_outcome_col, na.rm = TRUE)),
  sprintf("- Comparison: %d", sum(audit$has_comparison_col, na.rm = TRUE)),
  sprintf("- Subgroup: %d", sum(audit$has_subgroup_col, na.rm = TRUE)),
  sprintf("- review_doi: %d", sum(audit$has_review_doi_col, na.rm = TRUE)),
  sprintf("- review_url: %d", sum(audit$has_review_url_col, na.rm = TRUE)),
  "",
  "## Outcome types",
  paste0("- ", names(outcome_counts), ": ", as.integer(outcome_counts)),
  "",
  "## Missingness (median per dataset)",
  sprintf("- Study.year: %.2f", median(audit$missing_study_year, na.rm = TRUE)),
  sprintf("- Outcome: %.2f", median(audit$missing_outcome[audit$has_outcome_col], na.rm = TRUE)),
  sprintf("- Comparison: %.2f", median(audit$missing_comparison[audit$has_comparison_col], na.rm = TRUE)),
  sprintf("- review_doi: %.2f", median(audit$missing_review_doi[audit$has_review_doi_col], na.rm = TRUE)),
  "",
  "## Baseline risk proxies (binary)",
  sprintf("- Mean control risk (median across datasets): %.3f", median(audit$ctrl_risk_mean, na.rm = TRUE)),
  sprintf("- Median control risk (median across datasets): %.3f", median(audit$ctrl_risk_median, na.rm = TRUE)),
  "",
  "## Continuous control mean (where available)",
  sprintf("- Mean control mean (median across datasets): %.3f", median(audit$ctrl_mean_mean, na.rm = TRUE))
)

writeLines(summary_lines, summary_path)

cat("Wrote:", audit_path, "\n")
cat("Wrote:", summary_path, "\n")
