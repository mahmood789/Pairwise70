#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop("Package 'metafor' is required. Install with install.packages('metafor').")
  }
})

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", cmd_args[grep("^--file=", cmd_args)])
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(file_arg)))
  }
  return(getwd())
}

standardize_names <- function(nms) {
  nms <- gsub("[^A-Za-z0-9]+", ".", nms)
  nms <- gsub("\\.+", ".", nms)
  nms <- gsub("^\\.|\\.$", "", nms)
  nms
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

get_col <- function(df, name) {
  if (name %in% names(df)) {
    return(df[[name]])
  }
  rep(NA, nrow(df))
}

count_nonint <- function(x) {
  sum(!is.na(x) & abs(x - round(x)) > 1e-6)
}

sum_flag <- function(x) {
  sum(!is.na(x) & x, na.rm = TRUE)
}

ci_mismatch_count <- function(mean_val, se_val, ci_start, ci_end, tol = 0.01) {
  ok <- !is.na(mean_val) & !is.na(se_val) & !is.na(ci_start) & !is.na(ci_end)
  if (!any(ok)) {
    return(0)
  }
  calc_start <- mean_val - 1.96 * se_val
  calc_end <- mean_val + 1.96 * se_val
  mismatch <- ok & (abs(calc_start - ci_start) > tol | abs(calc_end - ci_end) > tol)
  sum(mismatch, na.rm = TRUE)
}

variance_mismatch_count <- function(var_val, se_val, rel_tol = 0.05, abs_tol = 1e-6) {
  ok <- !is.na(var_val) & !is.na(se_val)
  if (!any(ok)) {
    return(0)
  }
  diff <- abs(var_val - se_val^2)
  tol <- abs(var_val) * rel_tol + abs_tol
  sum(ok & diff > tol, na.rm = TRUE)
}

issue_count <- function(issue_totals, name) {
  val <- issue_totals$count[issue_totals$issue == name]
  if (length(val) == 0) {
    return(0)
  }
  val[1]
}

script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else file.path(project_root, "data")
output_dir <- if (length(args) >= 2) args[2] else file.path(project_root, "analysis", "output")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

binary_cols <- c("Experimental.cases", "Experimental.N", "Control.cases", "Control.N")
continuous_cols <- c("Experimental.mean", "Experimental.SD", "Experimental.N",
                     "Control.mean", "Control.SD", "Control.N")

ci_tol <- 0.01
var_rel_tol <- 0.05
fe_re_diff_threshold <- 0.2
high_i2_threshold <- 75
sparse_event_threshold <- 5

rda_files <- list.files(data_dir, pattern = "\\.rda$", full.names = TRUE)
if (length(rda_files) == 0) {
  stop(paste0("No .rda files found in: ", data_dir))
}

dataset_summary_list <- list()
issue_summary_list <- list()
analysis_results_list <- list()

column_counts <- integer(0)

cat("Pairwise70 diagnostics\n")
cat("Data dir: ", data_dir, "\n", sep = "")
cat("Output dir: ", output_dir, "\n\n", sep = "")

pb <- txtProgressBar(min = 0, max = length(rda_files), style = 3)

for (i in seq_along(rda_files)) {
  file <- rda_files[i]
  dataset_name <- sub("\\.rda$", "", basename(file))

  env <- new.env()
  load(file, envir = env)
  obj_names <- ls(env)
  if (length(obj_names) == 0) {
    next
  }
  if (length(obj_names) > 1) {
    warning(paste0("Multiple objects in ", dataset_name, "; using first."))
  }
  df <- env[[obj_names[1]]]
  if (!is.data.frame(df)) {
    warning(paste0("Skipping non-data.frame: ", dataset_name))
    next
  }

  names(df) <- standardize_names(names(df))
  cols_unique <- unique(names(df))
  for (nm in cols_unique) {
    if (!nm %in% names(column_counts)) {
      column_counts[nm] <- 0
    }
    column_counts[nm] <- column_counts[nm] + 1
  }

  n_rows <- nrow(df)
  n_cols <- ncol(df)

  analysis_number <- if ("Analysis.number" %in% names(df)) as.character(df$Analysis.number) else rep("1", n_rows)
  subgroup_number <- if ("Subgroup.number" %in% names(df)) as.character(df$Subgroup.number) else rep(NA, n_rows)
  subgroup_id <- ifelse(is.na(subgroup_number) | subgroup_number == "", "overall", subgroup_number)
  analysis_id <- paste(analysis_number, subgroup_id, sep = "::")

  n_analyses <- length(unique(analysis_id))
  n_subgroups <- length(unique(subgroup_id))

  has_binary_cols <- all(binary_cols %in% names(df))
  has_continuous_cols <- all(continuous_cols %in% names(df))

  binary_complete <- if (has_binary_cols) complete.cases(df[, binary_cols]) else rep(FALSE, n_rows)
  continuous_complete <- if (has_continuous_cols) complete.cases(df[, continuous_cols]) else rep(FALSE, n_rows)

  dataset_summary_list[[length(dataset_summary_list) + 1]] <- data.frame(
    dataset_name = dataset_name,
    n_rows = n_rows,
    n_cols = n_cols,
    n_analyses = n_analyses,
    n_subgroups = n_subgroups,
    has_binary_cols = has_binary_cols,
    has_continuous_cols = has_continuous_cols,
    binary_complete_rows = sum(binary_complete),
    continuous_complete_rows = sum(continuous_complete),
    stringsAsFactors = FALSE
  )

  exp_cases <- safe_num(get_col(df, "Experimental.cases"))
  exp_n <- safe_num(get_col(df, "Experimental.N"))
  ctrl_cases <- safe_num(get_col(df, "Control.cases"))
  ctrl_n <- safe_num(get_col(df, "Control.N"))

  exp_mean <- safe_num(get_col(df, "Experimental.mean"))
  exp_sd <- safe_num(get_col(df, "Experimental.SD"))
  ctrl_mean <- safe_num(get_col(df, "Control.mean"))
  ctrl_sd <- safe_num(get_col(df, "Control.SD"))

  giv_mean <- safe_num(get_col(df, "GIV.Mean"))
  giv_se <- safe_num(get_col(df, "GIV.SE"))
  mean_val <- safe_num(get_col(df, "Mean"))
  var_val <- safe_num(get_col(df, "Variance"))
  ci_start <- safe_num(get_col(df, "CI.start"))
  ci_end <- safe_num(get_col(df, "CI.end"))
  weight_val <- safe_num(get_col(df, "Weight"))

  mean_for_ci <- mean_val
  mean_for_ci[is.na(mean_for_ci)] <- giv_mean[is.na(mean_for_ci)]
  se_from_var <- ifelse(!is.na(var_val) & var_val >= 0, sqrt(var_val), NA_real_)
  se_for_ci <- giv_se
  se_for_ci[is.na(se_for_ci)] <- se_from_var[is.na(se_for_ci)]

  study_val <- if ("Study" %in% names(df)) as.character(df$Study) else rep(NA, n_rows)
  missing_study <- sum(is.na(study_val) | study_val == "", na.rm = TRUE)
  dup_key <- paste(analysis_id, study_val, sep = "|")
  duplicate_study <- sum(!is.na(study_val) & duplicated(dup_key))

  issue_summary_list[[length(issue_summary_list) + 1]] <- data.frame(
    dataset_name = dataset_name,
    missing_study = missing_study,
    duplicate_study = duplicate_study,
    exp_cases_gt_n = sum(exp_cases > exp_n, na.rm = TRUE),
    ctrl_cases_gt_n = sum(ctrl_cases > ctrl_n, na.rm = TRUE),
    exp_cases_lt0 = sum(exp_cases < 0, na.rm = TRUE),
    ctrl_cases_lt0 = sum(ctrl_cases < 0, na.rm = TRUE),
    exp_n_le0 = sum(exp_n <= 0, na.rm = TRUE),
    ctrl_n_le0 = sum(ctrl_n <= 0, na.rm = TRUE),
    exp_cases_nonint = count_nonint(exp_cases),
    ctrl_cases_nonint = count_nonint(ctrl_cases),
    exp_n_nonint = count_nonint(exp_n),
    ctrl_n_nonint = count_nonint(ctrl_n),
    double_zero = sum(exp_cases == 0 & ctrl_cases == 0, na.rm = TRUE),
    all_events = sum(exp_cases == exp_n & ctrl_cases == ctrl_n, na.rm = TRUE),
    sparse_events = sum((exp_cases + ctrl_cases) < sparse_event_threshold, na.rm = TRUE),
    exp_sd_le0 = sum(exp_sd <= 0, na.rm = TRUE),
    ctrl_sd_le0 = sum(ctrl_sd <= 0, na.rm = TRUE),
    variance_negative = sum(var_val < 0, na.rm = TRUE),
    giv_se_negative = sum(giv_se < 0, na.rm = TRUE),
    ci_inverted = sum(ci_start > ci_end, na.rm = TRUE),
    ci_mismatch = ci_mismatch_count(mean_for_ci, se_for_ci, ci_start, ci_end, tol = ci_tol),
    variance_se_mismatch = variance_mismatch_count(var_val, giv_se, rel_tol = var_rel_tol),
    weight_negative = sum(weight_val < 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  analysis_keys <- unique(analysis_id)
  for (analysis_key in analysis_keys) {
    idx <- analysis_id == analysis_key
    analysis_df <- df[idx, , drop = FALSE]

    analysis_num <- unique(analysis_number[idx])
    subgroup_num <- unique(subgroup_id[idx])
    analysis_name <- if ("Analysis.name" %in% names(analysis_df)) {
      unique(analysis_df$Analysis.name)
    } else {
      NA
    }
    analysis_name <- analysis_name[!is.na(analysis_name)]
    if (length(analysis_name) == 0) {
      analysis_name <- NA
    } else {
      analysis_name <- analysis_name[1]
    }

    exp_cases_a <- safe_num(get_col(analysis_df, "Experimental.cases"))
    exp_n_a <- safe_num(get_col(analysis_df, "Experimental.N"))
    ctrl_cases_a <- safe_num(get_col(analysis_df, "Control.cases"))
    ctrl_n_a <- safe_num(get_col(analysis_df, "Control.N"))

    exp_mean_a <- safe_num(get_col(analysis_df, "Experimental.mean"))
    exp_sd_a <- safe_num(get_col(analysis_df, "Experimental.SD"))
    ctrl_mean_a <- safe_num(get_col(analysis_df, "Control.mean"))
    ctrl_sd_a <- safe_num(get_col(analysis_df, "Control.SD"))

    giv_mean_a <- safe_num(get_col(analysis_df, "GIV.Mean"))
    giv_se_a <- safe_num(get_col(analysis_df, "GIV.SE"))
    mean_a <- safe_num(get_col(analysis_df, "Mean"))
    var_a <- safe_num(get_col(analysis_df, "Variance"))
    ci_start_a <- safe_num(get_col(analysis_df, "CI.start"))
    ci_end_a <- safe_num(get_col(analysis_df, "CI.end"))

    outcome_type <- "no_data"
    meta_status <- "no_data"
    meta_error <- ""
    k <- 0
    pooled_effect <- NA
    pooled_se <- NA
    ci_lb <- NA
    ci_ub <- NA
    i2 <- NA
    tau2 <- NA
    q_stat <- NA
    q_p <- NA
    fe_re_diff <- NA
    egger_p <- NA
    max_cooks <- NA
    max_dfbeta <- NA
    max_weight_share <- NA
    cc_applied <- FALSE

    double_zero_a <- sum(exp_cases_a == 0 & ctrl_cases_a == 0, na.rm = TRUE)
    all_events_a <- sum(exp_cases_a == exp_n_a & ctrl_cases_a == ctrl_n_a, na.rm = TRUE)
    sparse_events_a <- sum((exp_cases_a + ctrl_cases_a) < sparse_event_threshold, na.rm = TRUE)

    has_binary <- all(binary_cols %in% names(analysis_df))
    has_continuous <- all(continuous_cols %in% names(analysis_df))

    if (has_binary) {
      valid_binary <- !is.na(exp_cases_a) & !is.na(exp_n_a) & !is.na(ctrl_cases_a) & !is.na(ctrl_n_a) &
        exp_cases_a >= 0 & ctrl_cases_a >= 0 &
        exp_n_a > 0 & ctrl_n_a > 0 &
        exp_cases_a <= exp_n_a & ctrl_cases_a <= ctrl_n_a

      k <- sum(valid_binary)
      if (k >= 2) {
        outcome_type <- "binary"
        add_cc <- (double_zero_a > 0 || all_events_a > 0)
        cc_applied <- add_cc
        es <- tryCatch({
          metafor::escalc(
            measure = "OR",
            ai = exp_cases_a[valid_binary],
            bi = exp_n_a[valid_binary] - exp_cases_a[valid_binary],
            ci = ctrl_cases_a[valid_binary],
            di = ctrl_n_a[valid_binary] - ctrl_cases_a[valid_binary],
            add = if (add_cc) 0.5 else 0,
            to = "only0"
          )
        }, error = function(e) NULL)

        if (!is.null(es)) {
          meta_status <- "ok"
          ma_re <- tryCatch(metafor::rma(yi, vi, data = es, method = "REML"), error = function(e) e)
          ma_fe <- tryCatch(metafor::rma(yi, vi, data = es, method = "FE"), error = function(e) e)

          if (inherits(ma_re, "error")) {
            meta_status <- "error"
            meta_error <- ma_re$message
          } else {
            pooled_effect <- as.numeric(ma_re$beta)
            pooled_se <- ma_re$se
            ci_lb <- ma_re$ci.lb
            ci_ub <- ma_re$ci.ub
            i2 <- ma_re$I2
            tau2 <- ma_re$tau2
            q_stat <- ma_re$QE
            q_p <- ma_re$QEp

            if (!inherits(ma_fe, "error")) {
              fe_re_diff <- abs(as.numeric(ma_re$beta) - as.numeric(ma_fe$beta))
            }

            if (k >= 10) {
              reg <- tryCatch(metafor::regtest(ma_re, model = "lm"), error = function(e) NULL)
              if (!is.null(reg)) {
                egger_p <- reg$pval
              }
            }

            if (k >= 3 && k <= 200) {
              inf <- tryCatch(metafor::influence(ma_re), error = function(e) NULL)
              if (!is.null(inf)) {
                max_cooks <- max(inf$cook.d, na.rm = TRUE)
                max_dfbeta <- max(abs(inf$dfbs), na.rm = TRUE)
              }
            }

            weights <- 1 / (es$vi + ma_re$tau2)
            max_weight_share <- max(weights / sum(weights), na.rm = TRUE)
          }
        } else {
          meta_status <- "error"
          meta_error <- "escalc failed"
        }
      } else if (k > 0) {
        outcome_type <- "binary"
        meta_status <- "insufficient_data"
      }
    }

    if (outcome_type == "no_data" && has_continuous) {
      valid_cont <- !is.na(exp_mean_a) & !is.na(exp_sd_a) & !is.na(exp_n_a) &
        !is.na(ctrl_mean_a) & !is.na(ctrl_sd_a) & !is.na(ctrl_n_a) &
        exp_sd_a > 0 & ctrl_sd_a > 0 &
        exp_n_a > 0 & ctrl_n_a > 0

      k <- sum(valid_cont)
      if (k >= 2) {
        outcome_type <- "continuous"
        es <- tryCatch({
          metafor::escalc(
            measure = "SMD",
            m1i = exp_mean_a[valid_cont],
            sd1i = exp_sd_a[valid_cont],
            n1i = exp_n_a[valid_cont],
            m2i = ctrl_mean_a[valid_cont],
            sd2i = ctrl_sd_a[valid_cont],
            n2i = ctrl_n_a[valid_cont]
          )
        }, error = function(e) NULL)

        if (!is.null(es)) {
          meta_status <- "ok"
          ma_re <- tryCatch(metafor::rma(yi, vi, data = es, method = "REML"), error = function(e) e)
          ma_fe <- tryCatch(metafor::rma(yi, vi, data = es, method = "FE"), error = function(e) e)

          if (inherits(ma_re, "error")) {
            meta_status <- "error"
            meta_error <- ma_re$message
          } else {
            pooled_effect <- as.numeric(ma_re$beta)
            pooled_se <- ma_re$se
            ci_lb <- ma_re$ci.lb
            ci_ub <- ma_re$ci.ub
            i2 <- ma_re$I2
            tau2 <- ma_re$tau2
            q_stat <- ma_re$QE
            q_p <- ma_re$QEp

            if (!inherits(ma_fe, "error")) {
              fe_re_diff <- abs(as.numeric(ma_re$beta) - as.numeric(ma_fe$beta))
            }

            if (k >= 10) {
              reg <- tryCatch(metafor::regtest(ma_re, model = "lm"), error = function(e) NULL)
              if (!is.null(reg)) {
                egger_p <- reg$pval
              }
            }

            if (k >= 3 && k <= 200) {
              inf <- tryCatch(metafor::influence(ma_re), error = function(e) NULL)
              if (!is.null(inf)) {
                max_cooks <- max(inf$cook.d, na.rm = TRUE)
                max_dfbeta <- max(abs(inf$dfbs), na.rm = TRUE)
              }
            }

            weights <- 1 / (es$vi + ma_re$tau2)
            max_weight_share <- max(weights / sum(weights), na.rm = TRUE)
          }
        } else {
          meta_status <- "error"
          meta_error <- "escalc failed"
        }
      } else if (k > 0) {
        outcome_type <- "continuous"
        meta_status <- "insufficient_data"
      }
    }

    if (outcome_type == "no_data") {
      yi <- ifelse(!is.na(giv_mean_a), giv_mean_a, mean_a)
      sei <- giv_se_a
      if (all(is.na(sei)) && !all(is.na(var_a))) {
        sei <- sqrt(var_a)
      }
      if (all(is.na(sei)) && !all(is.na(ci_start_a)) && !all(is.na(ci_end_a))) {
        sei <- (ci_end_a - ci_start_a) / (2 * 1.96)
      }

      valid_generic <- !is.na(yi) & !is.na(sei) & sei > 0
      k <- sum(valid_generic)
      if (k >= 2) {
        outcome_type <- "generic"
        meta_status <- "ok"
        ma_re <- tryCatch(metafor::rma(yi = yi[valid_generic], sei = sei[valid_generic], method = "REML"), error = function(e) e)
        ma_fe <- tryCatch(metafor::rma(yi = yi[valid_generic], sei = sei[valid_generic], method = "FE"), error = function(e) e)

        if (inherits(ma_re, "error")) {
          meta_status <- "error"
          meta_error <- ma_re$message
        } else {
          pooled_effect <- as.numeric(ma_re$beta)
          pooled_se <- ma_re$se
          ci_lb <- ma_re$ci.lb
          ci_ub <- ma_re$ci.ub
          i2 <- ma_re$I2
          tau2 <- ma_re$tau2
          q_stat <- ma_re$QE
          q_p <- ma_re$QEp

          if (!inherits(ma_fe, "error")) {
            fe_re_diff <- abs(as.numeric(ma_re$beta) - as.numeric(ma_fe$beta))
          }

          if (k >= 10) {
            reg <- tryCatch(metafor::regtest(ma_re, model = "lm"), error = function(e) NULL)
            if (!is.null(reg)) {
              egger_p <- reg$pval
            }
          }

          if (k >= 3 && k <= 200) {
            inf <- tryCatch(metafor::influence(ma_re), error = function(e) NULL)
            if (!is.null(inf)) {
              max_cooks <- max(inf$cook.d, na.rm = TRUE)
              max_dfbeta <- max(abs(inf$dfbs), na.rm = TRUE)
            }
          }

          weights <- 1 / (ma_re$vi + ma_re$tau2)
          max_weight_share <- max(weights / sum(weights), na.rm = TRUE)
        }
      } else if (k > 0) {
        outcome_type <- "generic"
        meta_status <- "insufficient_data"
      }
    }

    analysis_results_list[[length(analysis_results_list) + 1]] <- data.frame(
      dataset_name = dataset_name,
      analysis_key = analysis_key,
      analysis_number = analysis_num[1],
      subgroup_number = subgroup_num[1],
      analysis_name = ifelse(is.na(analysis_name), "", analysis_name),
      outcome_type = outcome_type,
      k = k,
      n_rows = nrow(analysis_df),
      meta_status = meta_status,
      meta_error = meta_error,
      pooled_effect = pooled_effect,
      pooled_se = pooled_se,
      ci_lb = ci_lb,
      ci_ub = ci_ub,
      i2 = i2,
      tau2 = tau2,
      q_stat = q_stat,
      q_p = q_p,
      fe_re_diff = fe_re_diff,
      egger_p = egger_p,
      max_cooks = max_cooks,
      max_dfbeta = max_dfbeta,
      max_weight_share = max_weight_share,
      cc_applied = cc_applied,
      double_zero = double_zero_a,
      all_events = all_events_a,
      sparse_events = sparse_events_a,
      stringsAsFactors = FALSE
    )
  }

  setTxtProgressBar(pb, i)
}

close(pb)

dataset_summary <- do.call(rbind, dataset_summary_list)
issue_summary <- do.call(rbind, issue_summary_list)
analysis_results <- do.call(rbind, analysis_results_list)

write.csv(dataset_summary, file.path(output_dir, "dataset_inventory.csv"), row.names = FALSE)
write.csv(issue_summary, file.path(output_dir, "dataset_issue_counts.csv"), row.names = FALSE)
write.csv(analysis_results, file.path(output_dir, "analysis_diagnostics_results.csv"), row.names = FALSE)

column_counts_df <- data.frame(column = names(column_counts), count = as.integer(column_counts), stringsAsFactors = FALSE)
column_counts_df <- column_counts_df[order(column_counts_df$count, decreasing = TRUE), ]
write.csv(column_counts_df, file.path(output_dir, "column_frequency.csv"), row.names = FALSE)

issue_totals <- data.frame(
  issue = names(issue_summary)[names(issue_summary) != "dataset_name"],
  count = colSums(issue_summary[, names(issue_summary) != "dataset_name", drop = FALSE], na.rm = TRUE),
  stringsAsFactors = FALSE
)
issue_totals <- issue_totals[order(issue_totals$count, decreasing = TRUE), ]
write.csv(issue_totals, file.path(output_dir, "issue_totals.csv"), row.names = FALSE)

analysis_ok <- analysis_results[analysis_results$meta_status == "ok", ]
analysis_fail <- analysis_results[analysis_results$meta_status != "ok", ]

total_analyses <- nrow(analysis_results)
total_ok <- nrow(analysis_ok)
total_fail <- nrow(analysis_fail)
total_rows <- sum(dataset_summary$n_rows)

high_i2 <- sum(analysis_ok$i2 >= high_i2_threshold, na.rm = TRUE)
egger_sig <- sum(analysis_ok$egger_p < 0.05, na.rm = TRUE)
fe_re_diff_high <- sum(analysis_ok$fe_re_diff > fe_re_diff_threshold, na.rm = TRUE)
high_weight_share <- sum(analysis_ok$max_weight_share > 0.5, na.rm = TRUE)

report_path <- file.path(output_dir, "diagnosis_report.md")
report_lines <- c(
  "# Pairwise70 meta-analysis diagnostics",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Coverage",
  paste0("- Datasets: ", nrow(dataset_summary)),
  paste0("- Total analyses: ", total_analyses),
  paste0("- Total rows: ", total_rows),
  paste0("- Meta-analyses succeeded: ", total_ok),
  paste0("- Meta-analyses failed or insufficient data: ", total_fail),
  "",
  "## Common data issues (row-level counts)",
  paste0("- Double-zero event rows: ", issue_count(issue_totals, "double_zero")),
  paste0("- All-event rows: ", issue_count(issue_totals, "all_events")),
  paste0("- Sparse-event rows (<", sparse_event_threshold, " total events): ", issue_count(issue_totals, "sparse_events")),
  paste0("- Cases > N (experimental): ", issue_count(issue_totals, "exp_cases_gt_n")),
  paste0("- Cases > N (control): ", issue_count(issue_totals, "ctrl_cases_gt_n")),
  paste0("- SD <= 0 (experimental): ", issue_count(issue_totals, "exp_sd_le0")),
  paste0("- SD <= 0 (control): ", issue_count(issue_totals, "ctrl_sd_le0")),
  paste0("- Variance < 0: ", issue_count(issue_totals, "variance_negative")),
  paste0("- SE < 0: ", issue_count(issue_totals, "giv_se_negative")),
  paste0("- CI inverted: ", issue_count(issue_totals, "ci_inverted")),
  paste0("- CI mismatch vs Mean/SE: ", issue_count(issue_totals, "ci_mismatch")),
  "",
  "## Meta-analysis stability signals",
  paste0("- High heterogeneity (I2 >= ", high_i2_threshold, "%): ", high_i2),
  paste0("- Egger test p < 0.05 (k>=10): ", egger_sig),
  paste0("- Large FE vs RE shift (|diff| > ", fe_re_diff_threshold, "): ", fe_re_diff_high),
  paste0("- Max weight share > 0.5: ", high_weight_share),
  "",
  "## Methodological remedies to consider",
  "- Zero-event or all-event studies: continuity corrections (0.5, treatment-arm, or empirical), Peto OR for rare events, or GLMM/beta-binomial models. Double-zero studies may be excluded for OR or analyzed with risk difference.",
  "- Sparse events: use exact methods, GLMM, or switch to risk ratio/risk difference; report sensitivity to different corrections.",
  "- SD missing/zero: verify extraction, impute SD from pooled SD or similar studies, or use mean difference when scales align.",
  "- Cases > N or negative values: fix data integrity before analysis; these rows should be corrected or excluded.",
  "- CI/SE/variance inconsistencies: re-derive SE from CI or variance and verify scale (log vs raw).",
  "- Large FE vs RE shifts or high influence: report sensitivity, consider robust variance or leave-one-out checks.",
  "",
  "## Output files",
  "- dataset_inventory.csv: dataset-level inventory and completeness",
  "- dataset_issue_counts.csv: dataset-level issue counts",
  "- analysis_diagnostics_results.csv: per-analysis meta-results and diagnostics",
  "- column_frequency.csv: column prevalence across datasets",
  "- issue_totals.csv: aggregate issue counts",
  "- diagnosis_report.md: this report"
)

writeLines(report_lines, report_path)

cat("\nDiagnostics complete.\n")
cat("Outputs written to: ", output_dir, "\n", sep = "")
