#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop("Package 'metafor' is required. Install with install.packages('metafor').")
  }
})

standardize_names <- function(nms) {
  nms <- gsub("[^A-Za-z0-9]+", ".", nms)
  nms <- gsub("\\.+", ".", nms)
  nms <- gsub("^\\.|\\.$", "", nms)
  nms
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

run_rma <- function(yi, vi, effect_log) {
  keep <- is.finite(yi) & is.finite(vi) & vi > 0
  yi <- yi[keep]
  vi <- vi[keep]
  k_used <- length(yi)
  if (k_used < 2) {
    return(list(
      status = "insufficient_data",
      error = "",
      method = NA,
      k_used = k_used,
      pooled = NA,
      se = NA,
      ci_lb = NA,
      ci_ub = NA,
      pooled_tr = NA,
      ci_lb_tr = NA,
      ci_ub_tr = NA
    ))
  }

  ma <- tryCatch(metafor::rma(yi, vi, method = "REML"), error = function(e) e)
  if (inherits(ma, "error")) {
    ma_fe <- tryCatch(metafor::rma(yi, vi, method = "FE"), error = function(e) e)
    if (inherits(ma_fe, "error")) {
      return(list(
        status = "error",
        error = ma$message,
        method = NA,
        k_used = k_used,
        pooled = NA,
        se = NA,
        ci_lb = NA,
        ci_ub = NA,
        pooled_tr = NA,
        ci_lb_tr = NA,
        ci_ub_tr = NA
      ))
    }
    ma <- ma_fe
    model_method <- "FE"
  } else {
    model_method <- "REML"
  }

  pooled <- as.numeric(ma$beta)
  ci_lb <- ma$ci.lb
  ci_ub <- ma$ci.ub
  if (effect_log) {
    pooled_tr <- exp(pooled)
    ci_lb_tr <- exp(ci_lb)
    ci_ub_tr <- exp(ci_ub)
  } else {
    pooled_tr <- pooled
    ci_lb_tr <- ci_lb
    ci_ub_tr <- ci_ub
  }

  list(
    status = "ok",
    error = "",
    method = model_method,
    k_used = k_used,
    pooled = pooled,
    se = ma$se,
    ci_lb = ci_lb,
    ci_ub = ci_ub,
    pooled_tr = pooled_tr,
    ci_lb_tr = ci_lb_tr,
    ci_ub_tr = ci_ub_tr
  )
}

clean_binary <- function(ai, n1, ci, n2, tol = 1e-6) {
  ai <- safe_num(ai)
  n1 <- safe_num(n1)
  ci <- safe_num(ci)
  n2 <- safe_num(n2)

  nonint_ai <- !is.na(ai) & abs(ai - round(ai)) > tol
  nonint_ci <- !is.na(ci) & abs(ci - round(ci)) > tol
  nonint_n1 <- !is.na(n1) & abs(n1 - round(n1)) > tol
  nonint_n2 <- !is.na(n2) & abs(n2 - round(n2)) > tol

  ai[!nonint_ai] <- round(ai[!nonint_ai])
  ci[!nonint_ci] <- round(ci[!nonint_ci])
  n1[!nonint_n1] <- round(n1[!nonint_n1])
  n2[!nonint_n2] <- round(n2[!nonint_n2])

  invalid <- is.na(ai) | is.na(ci) | is.na(n1) | is.na(n2) |
    ai < 0 | ci < 0 |
    n1 <= 0 | n2 <= 0 |
    ai > n1 | ci > n2 |
    nonint_ai | nonint_ci | nonint_n1 | nonint_n2

  list(
    ai = ai,
    ci = ci,
    n1 = n1,
    n2 = n2,
    invalid = invalid,
    nonint_ai = nonint_ai,
    nonint_ci = nonint_ci,
    nonint_n1 = nonint_n1,
    nonint_n2 = nonint_n2
  )
}

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else file.path("data")
output_dir <- if (length(args) >= 2) args[2] else file.path("analysis", "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

rda_files <- list.files(data_dir, pattern = "\\.rda$", full.names = TRUE)
if (length(rda_files) == 0) {
  stop(paste0("No .rda files found in: ", data_dir))
}

results_list <- list()
drop_summary_list <- list()

cat("Fix-then-recompute OR_cc vs RD\n")
cat("Data dir: ", data_dir, "\n", sep = "")
cat("Output dir: ", output_dir, "\n\n", sep = "")

pb <- txtProgressBar(min = 0, max = length(rda_files), style = 3)

for (i in seq_along(rda_files)) {
  file <- rda_files[i]
  dataset_name <- sub("\\.rda$", "", basename(file))

  env <- new.env()
  load(file, envir = env)
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

  exp_cases <- if ("Experimental.cases" %in% names(df)) df$Experimental.cases else NA
  exp_n <- if ("Experimental.N" %in% names(df)) df$Experimental.N else NA
  ctrl_cases <- if ("Control.cases" %in% names(df)) df$Control.cases else NA
  ctrl_n <- if ("Control.N" %in% names(df)) df$Control.N else NA

  for (key in unique(analysis_key)) {
    idx <- analysis_key == key
    if (!any(idx)) next

    clean <- clean_binary(exp_cases[idx], exp_n[idx], ctrl_cases[idx], ctrl_n[idx])
    ai <- clean$ai
    ci <- clean$ci
    n1 <- clean$n1
    n2 <- clean$n2

    invalid <- clean$invalid
    valid <- !invalid
    k_valid <- sum(valid, na.rm = TRUE)
    if (k_valid < 2) {
      next
    }

    ai <- ai[valid]
    ci <- ci[valid]
    n1 <- n1[valid]
    n2 <- n2[valid]
    bi <- n1 - ai
    di <- n2 - ci

    double_zero <- sum(ai == 0 & ci == 0, na.rm = TRUE)
    sparse_events <- sum((ai + ci) < 5, na.rm = TRUE)

    analysis_name <- if ("Analysis.name" %in% names(df)) {
      nm <- unique(df$Analysis.name[idx])
      nm <- nm[!is.na(nm)]
      if (length(nm) > 0) nm[1] else ""
    } else {
      ""
    }

    base_meta <- list(
      dataset_name = dataset_name,
      analysis_key = key,
      analysis_number = unique(analysis_number[idx])[1],
      subgroup_number = unique(subgroup_id[idx])[1],
      analysis_name = analysis_name,
      k = k_valid,
      double_zero = double_zero,
      sparse_events = sparse_events
    )

    or_ai <- ai
    or_bi <- bi
    or_ci <- ci
    or_di <- di
    if (double_zero > 0) {
      drop00 <- !(or_ai == 0 & or_ci == 0)
      or_ai <- or_ai[drop00]
      or_bi <- or_bi[drop00]
      or_ci <- or_ci[drop00]
      or_di <- or_di[drop00]
    }

    es_or <- tryCatch(
      metafor::escalc(measure = "OR", ai = or_ai, bi = or_bi, ci = or_ci, di = or_di, add = 0.5, to = "only0"),
      error = function(e) NULL
    )
    if (!is.null(es_or)) {
      fit_or <- run_rma(es_or$yi, es_or$vi, TRUE)
      results_list[[length(results_list) + 1]] <- c(
        base_meta,
        list(
          method = "OR_cc",
          model_method = fit_or$method,
          k_used = fit_or$k_used,
          meta_status = fit_or$status,
          meta_error = fit_or$error,
          pooled_effect = fit_or$pooled,
          pooled_se = fit_or$se,
          ci_lb = fit_or$ci_lb,
          ci_ub = fit_or$ci_ub,
          pooled_effect_transformed = fit_or$pooled_tr,
          ci_lb_transformed = fit_or$ci_lb_tr,
          ci_ub_transformed = fit_or$ci_ub_tr
        )
      )
    }

    es_rd <- tryCatch(
      metafor::escalc(measure = "RD", ai = ai, bi = bi, ci = ci, di = di),
      error = function(e) NULL
    )
    if (!is.null(es_rd)) {
      fit_rd <- run_rma(es_rd$yi, es_rd$vi, FALSE)
      results_list[[length(results_list) + 1]] <- c(
        base_meta,
        list(
          method = "RD",
          model_method = fit_rd$method,
          k_used = fit_rd$k_used,
          meta_status = fit_rd$status,
          meta_error = fit_rd$error,
          pooled_effect = fit_rd$pooled,
          pooled_se = fit_rd$se,
          ci_lb = fit_rd$ci_lb,
          ci_ub = fit_rd$ci_ub,
          pooled_effect_transformed = fit_rd$pooled_tr,
          ci_lb_transformed = fit_rd$ci_lb_tr,
          ci_ub_transformed = fit_rd$ci_ub_tr
        )
      )
    }

    drop_summary_list[[length(drop_summary_list) + 1]] <- data.frame(
      dataset_name = dataset_name,
      analysis_key = key,
      k_total = length(clean$ai),
      k_valid = k_valid,
      invalid_rows = sum(invalid, na.rm = TRUE),
      nonint_ai = sum(clean$nonint_ai, na.rm = TRUE),
      nonint_ci = sum(clean$nonint_ci, na.rm = TRUE),
      nonint_n1 = sum(clean$nonint_n1, na.rm = TRUE),
      nonint_n2 = sum(clean$nonint_n2, na.rm = TRUE),
      double_zero = double_zero,
      sparse_events = sparse_events,
      stringsAsFactors = FALSE
    )
  }

  setTxtProgressBar(pb, i)
}

close(pb)

results_df <- do.call(rbind, lapply(results_list, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
drop_df <- do.call(rbind, drop_summary_list)

results_path <- file.path(output_dir, "orcc_vs_rd_fix_recompute_results.csv")
drop_path <- file.path(output_dir, "orcc_vs_rd_fix_recompute_drop_summary.csv")
write.csv(results_df, results_path, row.names = FALSE)
write.csv(drop_df, drop_path, row.names = FALSE)

rec_rows <- list()
orcc <- results_df[results_df$method == "OR_cc", ]
rd <- results_df[results_df$method == "RD", ]
orcc$key <- paste(orcc$dataset_name, orcc$analysis_key, sep = "|")
rd$key <- paste(rd$dataset_name, rd$analysis_key, sep = "|")

orcc_map <- split(orcc, orcc$key)
rd_map <- split(rd, rd$key)
all_keys <- unique(c(orcc$key, rd$key))

for (key in all_keys) {
  or_row <- orcc_map[[key]]
  rd_row <- rd_map[[key]]

  or_ok <- !is.null(or_row) && any(or_row$meta_status == "ok")
  rd_ok <- !is.null(rd_row) && any(rd_row$meta_status == "ok")

  if (or_ok) or_row <- or_row[or_row$meta_status == "ok", ][1, ]
  if (rd_ok) rd_row <- rd_row[rd_row$meta_status == "ok", ][1, ]

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

  rec_rows[[length(rec_rows) + 1]] <- data.frame(
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
    stringsAsFactors = FALSE
  )
}

rec_df <- do.call(rbind, rec_rows)
rec_path <- file.path(output_dir, "orcc_vs_rd_fix_recompute_recommendations.csv")
write.csv(rec_df, rec_path, row.names = FALSE)

original_path <- file.path(output_dir, "orcc_vs_rd_recommendations.csv")
comparison_path <- file.path(output_dir, "orcc_vs_rd_fix_recompute_comparison.csv")
changed <- NA

if (file.exists(original_path)) {
  orig <- read.csv(original_path, stringsAsFactors = FALSE)
  orig$key <- paste(orig$dataset_name, orig$analysis_key, sep = "|")
  rec_df$key <- paste(rec_df$dataset_name, rec_df$analysis_key, sep = "|")
  merged <- merge(orig, rec_df, by = "key", suffixes = c("_orig", "_fix"))
  merged$method_changed <- merged$recommended_method_orig != merged$recommended_method_fix
  changed <- sum(merged$method_changed, na.rm = TRUE)
  write.csv(merged, comparison_path, row.names = FALSE)
}

report_path <- file.path(output_dir, "orcc_vs_rd_fix_recompute_report.md")
method_counts <- table(rec_df$recommended_method)
report_lines <- c(
  "# OR_cc vs RD fix-then-recompute",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Summary",
  paste0("- Analyses with results: ", nrow(results_df)),
  paste0("- Recommendation rows: ", nrow(rec_df)),
  paste0("- Method changed vs original: ", ifelse(is.na(changed), "NA", changed)),
  "",
  "## Recommendation counts",
  paste0("- ", paste(names(method_counts), as.integer(method_counts), sep = ": ", collapse = "\n- ")),
  "",
  "## Output files",
  "- orcc_vs_rd_fix_recompute_results.csv",
  "- orcc_vs_rd_fix_recompute_drop_summary.csv",
  "- orcc_vs_rd_fix_recompute_recommendations.csv",
  "- orcc_vs_rd_fix_recompute_comparison.csv",
  "- orcc_vs_rd_fix_recompute_report.md"
)
writeLines(report_lines, report_path)

cat("Fix-then-recompute complete. Report: ", report_path, "\n", sep = "")
