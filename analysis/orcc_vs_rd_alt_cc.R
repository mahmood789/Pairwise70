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

get_col <- function(df, name) {
  if (name %in% names(df)) {
    return(df[[name]])
  }
  rep(NA, nrow(df))
}

run_rma <- function(yi, vi, effect_log) {
  keep <- is.finite(yi) & is.finite(vi) & vi > 0
  yi <- yi[keep]
  vi <- vi[keep]
  k_used <- length(yi)
  if (k_used < 2) {
    return(list(
      meta_status = "insufficient_data",
      meta_error = "",
      model_method = NA,
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
        meta_status = "error",
        meta_error = ma$message,
        model_method = NA,
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
    meta_status = "ok",
    meta_error = "",
    model_method = model_method,
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

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else file.path("data")
output_dir <- if (length(args) >= 2) args[2] else file.path("analysis", "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

sparse_event_threshold <- 5

rda_files <- list.files(data_dir, pattern = "\\.rda$", full.names = TRUE)
if (length(rda_files) == 0) {
  stop(paste0("No .rda files found in: ", data_dir))
}

cc_methods <- list(
  OR_cc = list(add = 0.5, to = "only0"),
  OR_cc_all = list(add = 0.5, to = "all"),
  OR_cc_0.1 = list(add = 0.1, to = "only0")
)

results_list <- list()

cat("OR_cc vs RD alternative CC analysis\n")
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
  df <- env[[obj_names[1]]]
  if (!is.data.frame(df)) {
    next
  }

  names(df) <- standardize_names(names(df))
  n_rows <- nrow(df)

  analysis_number <- if ("Analysis.number" %in% names(df)) as.character(df$Analysis.number) else rep("1", n_rows)
  subgroup_number <- if ("Subgroup.number" %in% names(df)) as.character(df$Subgroup.number) else rep(NA, n_rows)
  subgroup_id <- ifelse(is.na(subgroup_number) | subgroup_number == "", "overall", subgroup_number)
  analysis_key <- paste(analysis_number, subgroup_id, sep = "::")

  exp_cases <- safe_num(get_col(df, "Experimental.cases"))
  exp_n <- safe_num(get_col(df, "Experimental.N"))
  ctrl_cases <- safe_num(get_col(df, "Control.cases"))
  ctrl_n <- safe_num(get_col(df, "Control.N"))

  for (key in unique(analysis_key)) {
    idx <- analysis_key == key
    if (!any(idx)) next

    ai <- exp_cases[idx]
    n1 <- exp_n[idx]
    ci <- ctrl_cases[idx]
    n2 <- ctrl_n[idx]

    valid <- !is.na(ai) & !is.na(n1) & !is.na(ci) & !is.na(n2) &
      ai >= 0 & ci >= 0 & n1 > 0 & n2 > 0 &
      ai <= n1 & ci <= n2
    if (sum(valid) < 2) {
      next
    }

    ai <- ai[valid]
    bi <- n1[valid] - ai
    ci <- ci[valid]
    di <- n2[valid] - ci

    double_zero <- sum(ai == 0 & ci == 0, na.rm = TRUE)
    sparse_events <- sum((ai + ci) < sparse_event_threshold, na.rm = TRUE)

    analysis_meta <- list(
      dataset_name = dataset_name,
      analysis_key = key,
      analysis_number = unique(analysis_number[idx])[1],
      subgroup_number = unique(subgroup_id[idx])[1],
      analysis_name = if ("Analysis.name" %in% names(df)) {
        unique(df$Analysis.name[idx])[1]
      } else {
        ""
      },
      k = length(ai),
      double_zero = double_zero,
      sparse_events = sparse_events
    )

    es_rd <- tryCatch(
      metafor::escalc(measure = "RD", ai = ai, bi = bi, ci = ci, di = di),
      error = function(e) NULL
    )
    if (!is.null(es_rd)) {
      rd_fit <- run_rma(es_rd$yi, es_rd$vi, FALSE)
      results_list[[length(results_list) + 1]] <- c(
        analysis_meta,
        list(
          method = "RD",
          model_method = rd_fit$model_method,
          meta_status = rd_fit$meta_status,
          meta_error = rd_fit$meta_error,
          pooled_effect = rd_fit$pooled,
          pooled_se = rd_fit$se,
          ci_lb = rd_fit$ci_lb,
          ci_ub = rd_fit$ci_ub,
          pooled_effect_transformed = rd_fit$pooled_tr,
          ci_lb_transformed = rd_fit$ci_lb_tr,
          ci_ub_transformed = rd_fit$ci_ub_tr,
          k_used = rd_fit$k_used
        )
      )
    }

    for (m in names(cc_methods)) {
      cc <- cc_methods[[m]]
      es_or <- tryCatch(
        metafor::escalc(
          measure = "OR",
          ai = ai, bi = bi, ci = ci, di = di,
          add = cc$add, to = cc$to
        ),
        error = function(e) NULL
      )
      if (is.null(es_or)) {
        next
      }
      or_fit <- run_rma(es_or$yi, es_or$vi, TRUE)
      results_list[[length(results_list) + 1]] <- c(
        analysis_meta,
        list(
          method = m,
          model_method = or_fit$model_method,
          meta_status = or_fit$meta_status,
          meta_error = or_fit$meta_error,
          pooled_effect = or_fit$pooled,
          pooled_se = or_fit$se,
          ci_lb = or_fit$ci_lb,
          ci_ub = or_fit$ci_ub,
          pooled_effect_transformed = or_fit$pooled_tr,
          ci_lb_transformed = or_fit$ci_lb_tr,
          ci_ub_transformed = or_fit$ci_ub_tr,
          k_used = or_fit$k_used
        )
      )
    }
  }

  setTxtProgressBar(pb, i)
}

close(pb)

results_df <- do.call(rbind, lapply(results_list, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
results_path <- file.path(output_dir, "orcc_vs_rd_alt_cc_results.csv")
write.csv(results_df, results_path, row.names = FALSE)

summaries <- list()
rd <- results_df[results_df$method == "RD" & results_df$meta_status == "ok", ]
rd$key <- paste(rd$dataset_name, rd$analysis_key, sep = "|")

for (m in names(cc_methods)) {
  or <- results_df[results_df$method == m & results_df$meta_status == "ok", ]
  or$key <- paste(or$dataset_name, or$analysis_key, sep = "|")
  merged <- merge(or, rd, by = "key", suffixes = c("_or", "_rd"))

  log_or <- as.numeric(merged$pooled_effect_or)
  rd_eff <- as.numeric(merged$pooled_effect_rd)
  sign_match <- sign(log_or) == sign(rd_eff)
  sign_match[log_or == 0 | rd_eff == 0] <- NA

  or_sig <- (merged$ci_lb_or > 0) | (merged$ci_ub_or < 0)
  rd_sig <- (merged$ci_lb_rd > 0) | (merged$ci_ub_rd < 0)

  summaries[[length(summaries) + 1]] <- data.frame(
    method = m,
    paired_analyses = nrow(merged),
    sign_match_rate = round(mean(sign_match, na.rm = TRUE), 3),
    spearman = round(cor(log_or, rd_eff, use = "complete.obs", method = "spearman"), 3),
    both_sig = sum(or_sig & rd_sig, na.rm = TRUE),
    or_sig_only = sum(or_sig & !rd_sig, na.rm = TRUE),
    rd_sig_only = sum(!or_sig & rd_sig, na.rm = TRUE),
    neither_sig = sum(!or_sig & !rd_sig, na.rm = TRUE),
    mismatch_rate = round(mean(or_sig != rd_sig | sign_match == FALSE, na.rm = TRUE), 3),
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summaries)
summary_path <- file.path(output_dir, "orcc_vs_rd_alt_cc_summary.csv")
write.csv(summary_df, summary_path, row.names = FALSE)

report_path <- file.path(output_dir, "orcc_vs_rd_alt_cc_report.md")
report_lines <- c(
  "# OR_cc vs RD alternative CC summary",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Summary table",
  paste0(
    "- ",
    paste(
      paste(
        summary_df$method,
        "paired=", summary_df$paired_analyses,
        "sign_match=", summary_df$sign_match_rate,
        "spearman=", summary_df$spearman,
        "mismatch_rate=", summary_df$mismatch_rate,
        sep = " "
      ),
      collapse = "\n- "
    )
  ),
  "",
  "## Output files",
  "- orcc_vs_rd_alt_cc_results.csv",
  "- orcc_vs_rd_alt_cc_summary.csv",
  "- orcc_vs_rd_alt_cc_report.md"
)
writeLines(report_lines, report_path)

cat("Alternative CC summary written to: ", report_path, "\n", sep = "")
