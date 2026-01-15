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

infer_measure <- function(name) {
  if (is.null(name) || is.na(name) || name == "") {
    return(NA_character_)
  }
  x <- tolower(name)
  if (grepl("peto", x)) return("PETO")
  if (grepl("odds ratio", x)) return("OR")
  if (grepl("risk ratio|relative risk", x)) return("RR")
  if (grepl("risk difference", x)) return("RD")
  if (grepl("hazard ratio", x)) return("HR")
  if (grepl("rate ratio", x)) return("RR")
  if (grepl("standardised mean difference|standardized mean difference|\\bsmd\\b", x)) return("SMD")
  if (grepl("mean difference|weighted mean difference|\\bwmd\\b|\\bmd\\b", x)) return("MD")
  NA_character_
}

effect_is_log <- function(measure) {
  measure %in% c("OR", "RR", "HR", "PETO")
}

run_rma <- function(yi, vi, method_label, effect_log, analysis_meta) {
  keep <- is.finite(yi) & is.finite(vi) & vi > 0
  yi <- yi[keep]
  vi <- vi[keep]
  k_used <- length(yi)
  if (k_used < 2) {
    return(c(analysis_meta, list(
      method = method_label,
      model_method = NA,
      k_used = k_used,
      meta_status = "insufficient_data",
      meta_error = "",
      pooled_effect = NA,
      pooled_se = NA,
      ci_lb = NA,
      ci_ub = NA,
      pooled_effect_transformed = NA,
      ci_lb_transformed = NA,
      ci_ub_transformed = NA,
      i2 = NA,
      tau2 = NA,
      q_stat = NA,
      q_p = NA
    )))
  }
  ma <- tryCatch(metafor::rma(yi, vi, method = "REML"), error = function(e) e)
  if (inherits(ma, "error")) {
    ma_fe <- tryCatch(metafor::rma(yi, vi, method = "FE"), error = function(e) e)
    if (inherits(ma_fe, "error")) {
      return(c(analysis_meta, list(
        method = method_label,
        model_method = NA,
        k_used = k_used,
        meta_status = "error",
        meta_error = ma$message,
        pooled_effect = NA,
        pooled_se = NA,
        ci_lb = NA,
        ci_ub = NA,
        pooled_effect_transformed = NA,
        ci_lb_transformed = NA,
        ci_ub_transformed = NA,
        i2 = NA,
        tau2 = NA,
        q_stat = NA,
        q_p = NA
      )))
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

  c(analysis_meta, list(
    method = method_label,
    model_method = model_method,
    k_used = k_used,
    meta_status = "ok",
    meta_error = "",
    pooled_effect = pooled,
    pooled_se = ma$se,
    ci_lb = ci_lb,
    ci_ub = ci_ub,
    pooled_effect_transformed = pooled_tr,
    ci_lb_transformed = ci_lb_tr,
    ci_ub_transformed = ci_ub_tr,
    i2 = ma$I2,
    tau2 = ma$tau2,
    q_stat = ma$QE,
    q_p = ma$QEp
  ))
}

run_glmm_or <- function(ai, bi, ci, di, analysis_meta) {
  keep <- is.finite(ai) & is.finite(bi) & is.finite(ci) & is.finite(di) &
    ai >= 0 & bi >= 0 & ci >= 0 & di >= 0
  ai <- ai[keep]
  bi <- bi[keep]
  ci <- ci[keep]
  di <- di[keep]
  drop00 <- !(ai == 0 & ci == 0)
  ai <- ai[drop00]
  bi <- bi[drop00]
  ci <- ci[drop00]
  di <- di[drop00]
  k_used <- length(ai)
  if (k_used < 2) {
    return(c(analysis_meta, list(
      method = "GLMM_OR",
      model_method = "GLMM",
      k_used = k_used,
      meta_status = "insufficient_data",
      meta_error = "",
      pooled_effect = NA,
      pooled_se = NA,
      ci_lb = NA,
      ci_ub = NA,
      pooled_effect_transformed = NA,
      ci_lb_transformed = NA,
      ci_ub_transformed = NA,
      i2 = NA,
      tau2 = NA,
      q_stat = NA,
      q_p = NA
    )))
  }
  ma <- tryCatch(
    metafor::rma.glmm(
      ai = ai,
      bi = bi,
      ci = ci,
      di = di,
      measure = "OR",
      model = "UM.FS",
      method = "ML",
      add = 0.5,
      to = "only0",
      drop00 = TRUE,
      nAGQ = 1
    ),
    error = function(e) e
  )
  if (inherits(ma, "error")) {
    return(c(analysis_meta, list(
      method = "GLMM_OR",
      model_method = "GLMM",
      k_used = k_used,
      meta_status = "error",
      meta_error = ma$message,
      pooled_effect = NA,
      pooled_se = NA,
      ci_lb = NA,
      ci_ub = NA,
      pooled_effect_transformed = NA,
      ci_lb_transformed = NA,
      ci_ub_transformed = NA,
      i2 = NA,
      tau2 = NA,
      q_stat = NA,
      q_p = NA
    )))
  }

  pooled <- as.numeric(ma$beta)
  ci_lb <- ma$ci.lb
  ci_ub <- ma$ci.ub

  c(analysis_meta, list(
    method = "GLMM_OR",
    model_method = "GLMM",
    k_used = k_used,
    meta_status = "ok",
    meta_error = "",
    pooled_effect = pooled,
    pooled_se = ma$se,
    ci_lb = ci_lb,
    ci_ub = ci_ub,
    pooled_effect_transformed = exp(pooled),
    ci_lb_transformed = exp(ci_lb),
    ci_ub_transformed = exp(ci_ub),
    i2 = NA,
    tau2 = NA,
    q_stat = NA,
    q_p = NA
  ))
}

script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else file.path(project_root, "data")
output_dir <- if (length(args) >= 2) args[2] else file.path(project_root, "analysis", "output")
glmm_max_k_arg <- if (length(args) >= 3) suppressWarnings(as.numeric(args[3])) else NA_real_
clean_dir <- file.path(output_dir, "cleaned_rds")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(clean_dir, showWarnings = FALSE, recursive = TRUE)

sparse_event_threshold <- 5
peto_event_rate_threshold <- 0.1
glmm_max_k <- if (!is.na(glmm_max_k_arg)) glmm_max_k_arg else 30
glmm_enabled <- !is.na(glmm_max_k) && glmm_max_k > 0

binary_cols <- c("Experimental.cases", "Experimental.N", "Control.cases", "Control.N")
continuous_cols <- c("Experimental.mean", "Experimental.SD", "Experimental.N",
                     "Control.mean", "Control.SD", "Control.N")

rda_files <- list.files(data_dir, pattern = "\\.rda$", full.names = TRUE)
if (length(rda_files) == 0) {
  stop(paste0("No .rda files found in: ", data_dir))
}

results_list <- list()
dataset_summary_list <- list()

cat("Pairwise70 remediation pipeline\n")
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
  analysis_id <- paste(analysis_number, subgroup_id, sep = "::")

  df_clean <- df
  exp_sd <- safe_num(get_col(df, "Experimental.SD"))
  ctrl_sd <- safe_num(get_col(df, "Control.SD"))
  exp_sd_clean <- exp_sd
  ctrl_sd_clean <- ctrl_sd
  sd_imputed <- rep(FALSE, n_rows)

  for (analysis_key in unique(analysis_id)) {
    idx <- which(analysis_id == analysis_key)
    exp_sd_a <- exp_sd[idx]
    ctrl_sd_a <- ctrl_sd[idx]
    sd_pool <- c(exp_sd_a, ctrl_sd_a)
    sd_pool <- sd_pool[!is.na(sd_pool) & sd_pool > 0]
    sd_median <- if (length(sd_pool) > 0) median(sd_pool) else NA_real_

    for (j in idx) {
      if (is.na(exp_sd_clean[j]) || exp_sd_clean[j] <= 0) {
        if (!is.na(ctrl_sd_clean[j]) && ctrl_sd_clean[j] > 0) {
          exp_sd_clean[j] <- ctrl_sd_clean[j]
          sd_imputed[j] <- TRUE
        } else if (!is.na(sd_median)) {
          exp_sd_clean[j] <- sd_median
          sd_imputed[j] <- TRUE
        }
      }
      if (is.na(ctrl_sd_clean[j]) || ctrl_sd_clean[j] <= 0) {
        if (!is.na(exp_sd_clean[j]) && exp_sd_clean[j] > 0) {
          ctrl_sd_clean[j] <- exp_sd_clean[j]
          sd_imputed[j] <- TRUE
        } else if (!is.na(sd_median)) {
          ctrl_sd_clean[j] <- sd_median
          sd_imputed[j] <- TRUE
        }
      }
    }
  }

  df_clean$Experimental.SD.cleaned <- exp_sd_clean
  df_clean$Control.SD.cleaned <- ctrl_sd_clean
  df_clean$SD.imputed <- sd_imputed

  saveRDS(df_clean, file.path(clean_dir, paste0(dataset_name, ".rds")))

  dataset_sd_imputed <- sum(sd_imputed, na.rm = TRUE)
  dataset_invalid_binary <- 0
  dataset_invalid_cont <- 0
  analyses_sparse <- 0
  analyses_double_zero <- 0
  analyses_with_cc <- 0
  analyses_with_glmm <- 0

  for (analysis_key in unique(analysis_id)) {
    idx <- analysis_id == analysis_key
    analysis_df <- df_clean[idx, , drop = FALSE]

    analysis_num <- unique(analysis_number[idx])[1]
    subgroup_num <- unique(subgroup_id[idx])[1]
    analysis_name <- if ("Analysis.name" %in% names(analysis_df)) {
      unique(analysis_df$Analysis.name)
    } else {
      NA
    }
    analysis_name <- analysis_name[!is.na(analysis_name)][1]

    exp_cases <- safe_num(get_col(analysis_df, "Experimental.cases"))
    exp_n <- safe_num(get_col(analysis_df, "Experimental.N"))
    ctrl_cases <- safe_num(get_col(analysis_df, "Control.cases"))
    ctrl_n <- safe_num(get_col(analysis_df, "Control.N"))

    exp_mean <- safe_num(get_col(analysis_df, "Experimental.mean"))
    ctrl_mean <- safe_num(get_col(analysis_df, "Control.mean"))
    exp_n_c <- safe_num(get_col(analysis_df, "Experimental.N"))
    ctrl_n_c <- safe_num(get_col(analysis_df, "Control.N"))
    exp_sd_c <- safe_num(get_col(analysis_df, "Experimental.SD.cleaned"))
    ctrl_sd_c <- safe_num(get_col(analysis_df, "Control.SD.cleaned"))

    giv_mean <- safe_num(get_col(analysis_df, "GIV.Mean"))
    giv_se <- safe_num(get_col(analysis_df, "GIV.SE"))
    mean_val <- safe_num(get_col(analysis_df, "Mean"))
    var_val <- safe_num(get_col(analysis_df, "Variance"))
    ci_start <- safe_num(get_col(analysis_df, "CI.start"))
    ci_end <- safe_num(get_col(analysis_df, "CI.end"))

    outcome_type <- "no_data"
    measure_hint <- infer_measure(analysis_name)

    analysis_meta <- list(
      dataset_name = dataset_name,
      analysis_key = analysis_key,
      analysis_number = analysis_num,
      subgroup_number = subgroup_num,
      analysis_name = ifelse(is.na(analysis_name), "", analysis_name),
      k = NA,
      n_rows = nrow(analysis_df),
      outcome_type = NA,
      cc_applied = FALSE,
      sd_imputed_rows = sum(analysis_df$SD.imputed, na.rm = TRUE),
      se_derived_rows = 0,
      double_zero = sum(exp_cases == 0 & ctrl_cases == 0, na.rm = TRUE),
      sparse_events = sum((exp_cases + ctrl_cases) < sparse_event_threshold, na.rm = TRUE)
    )

    has_binary <- all(binary_cols %in% names(analysis_df))
    if (has_binary) {
      valid_binary <- !is.na(exp_cases) & !is.na(exp_n) & !is.na(ctrl_cases) & !is.na(ctrl_n) &
        exp_cases >= 0 & ctrl_cases >= 0 &
        exp_n > 0 & ctrl_n > 0 &
        exp_cases <= exp_n & ctrl_cases <= ctrl_n
      dataset_invalid_binary <- dataset_invalid_binary + sum(!valid_binary, na.rm = TRUE)
      k_bin <- sum(valid_binary)
      if (k_bin >= 2) {
        outcome_type <- "binary"
        analysis_meta$outcome_type <- "binary"
        analysis_meta$k <- k_bin

        ai <- exp_cases[valid_binary]
        bi <- exp_n[valid_binary] - exp_cases[valid_binary]
        ci <- ctrl_cases[valid_binary]
        di <- ctrl_n[valid_binary] - ctrl_cases[valid_binary]

        has_zero_cell <- any(ai == 0 | bi == 0 | ci == 0 | di == 0)
        cc_applied <- has_zero_cell
        analysis_meta$cc_applied <- cc_applied
        if (cc_applied) {
          analyses_with_cc <- analyses_with_cc + 1
        }

        total_events <- sum(ai + ci, na.rm = TRUE)
        total_n <- sum(ai + bi + ci + di, na.rm = TRUE)
        event_rate <- if (total_n > 0) total_events / total_n else NA
        if (analysis_meta$double_zero > 0) {
          analyses_double_zero <- analyses_double_zero + 1
        }
        if (analysis_meta$sparse_events > 0) {
          analyses_sparse <- analyses_sparse + 1
        }

        es_or <- tryCatch(
          metafor::escalc(measure = "OR", ai = ai, bi = bi, ci = ci, di = di, add = 0.5, to = "only0"),
          error = function(e) NULL
        )
        if (!is.null(es_or)) {
          results_list[[length(results_list) + 1]] <- run_rma(es_or$yi, es_or$vi, "OR_cc", TRUE, analysis_meta)
        }

        es_rr <- tryCatch(
          metafor::escalc(measure = "RR", ai = ai, bi = bi, ci = ci, di = di, add = 0.5, to = "only0"),
          error = function(e) NULL
        )
        if (!is.null(es_rr)) {
          results_list[[length(results_list) + 1]] <- run_rma(es_rr$yi, es_rr$vi, "RR_cc", TRUE, analysis_meta)
        }

        es_rd <- tryCatch(
          metafor::escalc(measure = "RD", ai = ai, bi = bi, ci = ci, di = di, add = 0),
          error = function(e) NULL
        )
        if (!is.null(es_rd)) {
          results_list[[length(results_list) + 1]] <- run_rma(es_rd$yi, es_rd$vi, "RD", FALSE, analysis_meta)
        }

        if (!is.na(event_rate) && event_rate <= peto_event_rate_threshold && analysis_meta$double_zero == 0) {
          es_peto <- tryCatch(
            metafor::escalc(measure = "PETO", ai = ai, bi = bi, ci = ci, di = di),
            error = function(e) NULL
          )
          if (!is.null(es_peto)) {
            results_list[[length(results_list) + 1]] <- run_rma(es_peto$yi, es_peto$vi, "PETO", TRUE, analysis_meta)
          }
        }

        if (glmm_enabled && (analysis_meta$sparse_events > 0 || analysis_meta$double_zero > 0) && k_bin <= glmm_max_k) {
          analyses_with_glmm <- analyses_with_glmm + 1
          results_list[[length(results_list) + 1]] <- run_glmm_or(ai, bi, ci, di, analysis_meta)
        }
      }
    }

    if (outcome_type == "no_data") {
      valid_cont <- !is.na(exp_mean) & !is.na(ctrl_mean) &
        !is.na(exp_sd_c) & !is.na(ctrl_sd_c) &
        exp_sd_c > 0 & ctrl_sd_c > 0 &
        !is.na(exp_n_c) & !is.na(ctrl_n_c) &
        exp_n_c > 0 & ctrl_n_c > 0
      dataset_invalid_cont <- dataset_invalid_cont + sum(!valid_cont, na.rm = TRUE)
      k_cont <- sum(valid_cont)
      if (k_cont >= 2) {
        outcome_type <- "continuous"
        analysis_meta$outcome_type <- "continuous"
        analysis_meta$k <- k_cont

        measure_cont <- if (!is.na(measure_hint) && measure_hint %in% c("MD", "SMD")) {
          measure_hint
        } else {
          "SMD"
        }

        es_cont <- tryCatch(
          metafor::escalc(
            measure = measure_cont,
            m1i = exp_mean[valid_cont],
            sd1i = exp_sd_c[valid_cont],
            n1i = exp_n_c[valid_cont],
            m2i = ctrl_mean[valid_cont],
            sd2i = ctrl_sd_c[valid_cont],
            n2i = ctrl_n_c[valid_cont]
          ),
          error = function(e) NULL
        )
        if (!is.null(es_cont)) {
          results_list[[length(results_list) + 1]] <- run_rma(
            es_cont$yi,
            es_cont$vi,
            paste0(measure_cont, "_impSD"),
            FALSE,
            analysis_meta
          )
        }
      }
    }

    if (outcome_type == "no_data") {
      valid_giv <- !is.na(giv_mean) & !is.na(giv_se) & giv_se > 0
      eff_log <- effect_is_log(ifelse(is.na(measure_hint), "GEN", measure_hint))
      if (eff_log) {
        valid_giv <- valid_giv & giv_mean > 0
      }
      k_giv <- sum(valid_giv)
      if (k_giv >= 2) {
        outcome_type <- "generic"
        analysis_meta$outcome_type <- "generic"
        analysis_meta$k <- k_giv

        yi <- if (eff_log) log(giv_mean[valid_giv]) else giv_mean[valid_giv]
        vi <- giv_se[valid_giv]^2
        results_list[[length(results_list) + 1]] <- run_rma(yi, vi, "GEN_GIV", eff_log, analysis_meta)
      } else {
        valid_mean_ci <- !is.na(mean_val) & !is.na(ci_start) & !is.na(ci_end)
        k_ci <- sum(valid_mean_ci)
        if (k_ci >= 2) {
          outcome_type <- "generic"
          analysis_meta$outcome_type <- "generic"
          analysis_meta$k <- k_ci

          eff_log <- effect_is_log(ifelse(is.na(measure_hint), "GEN", measure_hint))
          if (eff_log) {
            valid_mean_ci <- valid_mean_ci & mean_val > 0 & ci_start > 0 & ci_end > 0
            yi <- log(mean_val[valid_mean_ci])
            sei <- (log(ci_end[valid_mean_ci]) - log(ci_start[valid_mean_ci])) / (2 * 1.96)
          } else {
            yi <- mean_val[valid_mean_ci]
            sei <- (ci_end[valid_mean_ci] - ci_start[valid_mean_ci]) / (2 * 1.96)
          }

          analysis_meta$se_derived_rows <- sum(valid_mean_ci)
          vi <- sei^2
          results_list[[length(results_list) + 1]] <- run_rma(yi, vi, "GEN_CI", eff_log, analysis_meta)
        } else {
          valid_mean_var <- !is.na(mean_val) & !is.na(var_val) & var_val > 0
          k_var <- sum(valid_mean_var)
          if (k_var >= 2) {
            outcome_type <- "generic"
            analysis_meta$outcome_type <- "generic"
            analysis_meta$k <- k_var

            eff_log <- effect_is_log(ifelse(is.na(measure_hint), "GEN", measure_hint))
            yi <- if (eff_log) log(mean_val[valid_mean_var]) else mean_val[valid_mean_var]
            vi <- var_val[valid_mean_var]
            analysis_meta$se_derived_rows <- sum(valid_mean_var)
            results_list[[length(results_list) + 1]] <- run_rma(yi, vi, "GEN_VAR", eff_log, analysis_meta)
          }
        }
      }
    }
  }

  dataset_summary_list[[length(dataset_summary_list) + 1]] <- data.frame(
    dataset_name = dataset_name,
    n_rows = n_rows,
    n_analyses = length(unique(analysis_id)),
    sd_imputed_rows = dataset_sd_imputed,
    invalid_binary_rows = dataset_invalid_binary,
    invalid_cont_rows = dataset_invalid_cont,
    analyses_sparse = analyses_sparse,
    analyses_double_zero = analyses_double_zero,
    analyses_with_cc = analyses_with_cc,
    analyses_with_glmm = analyses_with_glmm,
    stringsAsFactors = FALSE
  )

  setTxtProgressBar(pb, i)
}

close(pb)

results_df <- do.call(rbind, lapply(results_list, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
dataset_summary <- do.call(rbind, dataset_summary_list)

write.csv(results_df, file.path(output_dir, "remediation_analysis_results.csv"), row.names = FALSE)
write.csv(dataset_summary, file.path(output_dir, "remediation_dataset_summary.csv"), row.names = FALSE)

method_counts <- table(results_df$method)
status_counts <- table(results_df$meta_status)
outcome_counts <- table(results_df$outcome_type)

report_path <- file.path(output_dir, "remediation_report.md")
report_lines <- c(
  "# Pairwise70 remediation pipeline",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Coverage",
  paste0("- Datasets processed: ", nrow(dataset_summary)),
  paste0("- Total analyses with results: ", nrow(results_df)),
  paste0("- GLMM enabled: ", ifelse(glmm_enabled, "yes", "no")),
  paste0("- GLMM max k: ", ifelse(glmm_enabled, glmm_max_k, "disabled")),
  "",
  "## Method counts",
  paste0("- ", paste(names(method_counts), as.integer(method_counts), sep = ": ", collapse = "\n- ")),
  "",
  "## Meta-analysis status",
  paste0("- ", paste(names(status_counts), as.integer(status_counts), sep = ": ", collapse = "\n- ")),
  "",
  "## Outcome types",
  paste0("- ", paste(names(outcome_counts), as.integer(outcome_counts), sep = ": ", collapse = "\n- ")),
  "",
  "## Outputs",
  "- remediation_analysis_results.csv: per-analysis results across remedial methods",
  "- remediation_dataset_summary.csv: dataset-level remediation counts",
  "- cleaned_rds/: datasets with SD-imputed columns and standardized names"
)
writeLines(report_lines, report_path)

cat("\nRemediation pipeline complete.\n")
cat("Outputs written to: ", output_dir, "\n", sep = "")
