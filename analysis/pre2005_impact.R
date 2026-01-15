#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(metafor)
  library(dplyr)
  library(ggplot2)
  library(lme4)
  library(clubSandwich)
  library(splines)
  library(robumeta)
  library(RoBMA)
})

data_dir <- file.path("data")
out_dir <- file.path("analysis", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
results_path <- file.path(out_dir, "analysis_results.csv")
issues_path <- file.path(out_dir, "issues.csv")
sens_path <- file.path(out_dir, "sensitivity_fixed_influence.csv")
sens_summary_path <- file.path(out_dir, "sensitivity_summary.txt")
sens_summary_csv <- file.path(out_dir, "sensitivity_summary.csv")
extended_summary_path <- file.path(out_dir, "extended_summary_2005.txt")
extended_summary_csv <- file.path(out_dir, "extended_summary_2005.csv")
extended_by_group_csv <- file.path(out_dir, "extended_summary_by_group_2005.csv")
analysis_type_summary_path <- file.path(out_dir, "analysis_type_summary_2005.csv")
analysis_type_summary_all_path <- file.path(out_dir, "analysis_type_summary_2005_all.csv")
sweep_summary_2yr_by_type_path <- file.path(out_dir, "cutoff_sweep_summary_2yr_by_type.csv")
sweep_summary_2yr_by_measure_path <- file.path(out_dir, "cutoff_sweep_summary_2yr_by_measure.csv")
sweep_summary_2yr_by_category_path <- file.path(out_dir, "cutoff_sweep_summary_2yr_by_category.csv")
category_summary_path <- file.path(out_dir, "category_summary_2005.csv")
category_summary_all_path <- file.path(out_dir, "category_summary_2005_all.csv")
category_summary_log_ratio_path <- file.path(out_dir, "category_summary_2005_log_ratio.csv")
category_summary_log_ratio_all_path <- file.path(out_dir, "category_summary_2005_log_ratio_all.csv")
sweep_summary_2yr_by_category_mean_diff_path <- file.path(out_dir, "cutoff_sweep_summary_2yr_by_category_mean_diff.csv")
category_summary_mean_diff_path <- file.path(out_dir, "category_summary_2005_mean_diff.csv")
category_summary_mean_diff_all_path <- file.path(out_dir, "category_summary_2005_mean_diff_all.csv")
sweep_summary_2yr_by_category_log_ratio_path <- file.path(out_dir, "cutoff_sweep_summary_2yr_by_category_log_ratio.csv")
rob_summary_by_category_path <- file.path(out_dir, "rob_low_summary_by_category.csv")
rob_summary_by_category_log_ratio_path <- file.path(out_dir, "rob_low_summary_by_category_log_ratio.csv")
rob_summary_by_category_mean_diff_path <- file.path(out_dir, "rob_low_summary_by_category_mean_diff.csv")
rob_summary_by_type_path <- file.path(out_dir, "rob_low_summary_by_type.csv")
rob_summary_by_measure_path <- file.path(out_dir, "rob_low_summary_by_measure.csv")
sweep_results_path <- file.path(out_dir, "cutoff_sweep_results.csv")
sweep_summary_path <- file.path(out_dir, "cutoff_sweep_summary.csv")
sweep_summary_2yr_path <- file.path(out_dir, "cutoff_sweep_summary_2yr.csv")
rob_results_path <- file.path(out_dir, "rob_low_results.csv")
rob_summary_path <- file.path(out_dir, "rob_low_summary.txt")
rob_summary_csv <- file.path(out_dir, "rob_low_summary.csv")
advanced_results_path <- file.path(out_dir, "advanced_model_results.csv")
advanced_summary_path <- file.path(out_dir, "advanced_model_summary.txt")
advanced_meta_summary_path <- file.path(out_dir, "advanced_model_meta_summary.csv")
advanced_meta_summary_log_ratio_path <- file.path(out_dir, "advanced_model_meta_summary_log_ratio.csv")
advanced_meta_summary_mean_diff_path <- file.path(out_dir, "advanced_model_meta_summary_mean_diff.csv")
advanced_mixed_path <- file.path(out_dir, "advanced_mixed_effects.csv")
advanced_mixed_logit_path <- file.path(out_dir, "advanced_mixed_effects_logit.csv")
rve_summary_path <- file.path(out_dir, "rve_summary.csv")
robma_summary_path <- file.path(out_dir, "robma_summary.txt")
robma_inference_path <- file.path(out_dir, "robma_inference.csv")
robma_models_path <- file.path(out_dir, "robma_models.csv")

thresholds <- c(2005, 2000, 2010)
sweep_cutoffs <- 1950:2024
bias_col <- "Bias.arising.from.the.randomization.process..judgement."
bias_support_col <- "Bias.arising.from.the.randomization.process..support."

is_low_risk <- function(x) {
  x <- tolower(trimws(x))
  x == "low risk"
}

is_non_high <- function(x) {
  x <- tolower(trimws(x))
  x %in% c("low risk", "some concerns")
}

has_support <- function(x) {
  x <- trimws(x)
  !is.na(x) & x != ""
}

categorize_analysis_name <- function(name) {
  if (is.na(name) || is.null(name)) {
    return("other")
  }
  x <- tolower(name)
  if (grepl("mortality|death|survival|fatal", x)) return("mortality")
  if (grepl("adverse|toxicity|side effect|complication|bleed|bleeding|haemorrhage|hemorrhage|infection|sepsis|harm|safety", x)) {
    return("adverse_events")
  }
  if (grepl("pain|analges", x)) return("pain")
  if (grepl("quality of life|qol|well[- ]?being|hrqol", x)) return("quality_of_life")
  if (grepl("length of stay|hospital stay|days in hospital|hospitalization", x)) return("length_of_stay")
  if (grepl("relapse|recurrence|progression", x)) return("relapse")
  if (grepl("response|remission|cure|recovery|resolution", x)) return("response")
  if (grepl("score|index|scale|rating|questionnaire|functional|disability", x)) return("function_score")
  if (grepl("blood pressure|hb|hba1c|cholesterol|lipid|glucose|viral load|titer|titre|concentration|level", x)) {
    return("biomarker")
  }
  if (grepl("cost|economic|resource|budget", x)) return("cost")
  "other"
}

measure_group <- function(measure, effect_is_log) {
  if (is.character(effect_is_log)) {
    effect_is_log <- tolower(effect_is_log) == "true"
  }
  if (is.factor(effect_is_log)) {
    effect_is_log <- tolower(as.character(effect_is_log)) == "true"
  }
  if (!is.na(measure)) {
    if (measure %in% c("RR", "OR", "HR", "PETO")) return("log_ratio")
    if (measure %in% c("MD", "SMD", "RD")) return("mean_diff")
    if (measure %in% c("GEN")) return("generic")
  }
  if (!is.na(effect_is_log) && effect_is_log) return("log_ratio")
  "other"
}

infer_measure <- function(name) {
  if (is.null(name) || is.na(name)) {
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

choose_method <- function(df, analysis_name) {
  measure_hint <- infer_measure(analysis_name)

  has_binary <- sum(
    !is.na(df$Experimental.cases) &
      !is.na(df$Experimental.N) &
      !is.na(df$Control.cases) &
      !is.na(df$Control.N)
  ) >= 2

  if (has_binary) {
    measure <- if (!is.na(measure_hint) && measure_hint %in% c("OR", "RR", "RD", "PETO")) {
      measure_hint
    } else {
      "RR"
    }
    effect_is_log <- measure %in% c("OR", "RR", "PETO")
    return(list(type = "binary", measure = measure, effect_is_log = effect_is_log))
  }

  has_cont <- sum(
    !is.na(df$Experimental.mean) &
      !is.na(df$Experimental.SD) &
      !is.na(df$Experimental.N) &
      !is.na(df$Control.mean) &
      !is.na(df$Control.SD) &
      !is.na(df$Control.N)
  ) >= 2

  if (has_cont) {
    measure <- if (!is.na(measure_hint) && measure_hint %in% c("MD", "SMD")) {
      measure_hint
    } else {
      "MD"
    }
    return(list(type = "continuous", measure = measure, effect_is_log = FALSE))
  }

  has_giv <- sum(!is.na(df$GIV.Mean) & !is.na(df$GIV.SE)) >= 2
  if (has_giv) {
    measure <- if (!is.na(measure_hint)) measure_hint else "GEN"
    effect_is_log <- !(measure %in% c("MD", "SMD", "RD"))
    return(list(type = "giv", measure = measure, effect_is_log = effect_is_log))
  }

  has_mean_ci <- sum(!is.na(df$Mean) & !is.na(df$CI.start) & !is.na(df$CI.end)) >= 2
  if (has_mean_ci) {
    if (!is.na(measure_hint)) {
      measure <- measure_hint
    } else {
      if (any(df$Mean < 0, na.rm = TRUE) || any(df$CI.start < 0, na.rm = TRUE)) {
        measure <- "MD"
      } else {
        measure <- "RR"
      }
    }
    effect_is_log <- measure %in% c("OR", "RR", "HR", "PETO")
    return(list(type = "mean_ci", measure = measure, effect_is_log = effect_is_log))
  }

  has_mean_var <- sum(!is.na(df$Mean) & !is.na(df$Variance)) >= 2
  if (has_mean_var) {
    measure <- if (!is.na(measure_hint)) measure_hint else "GEN"
    effect_is_log <- !(measure %in% c("MD", "SMD", "RD"))
    return(list(type = "mean_var", measure = measure, effect_is_log = effect_is_log))
  }

  NULL
}

extract_es <- function(df, method) {
  df <- df[!is.na(df$Study) & df$Study != "", , drop = FALSE]
  df$Study.year <- suppressWarnings(as.numeric(df$Study.year))
  type <- method$type
  effect_is_log <- method$effect_is_log

  if (type == "binary") {
    ai <- df$Experimental.cases
    n1 <- df$Experimental.N
    ci <- df$Control.cases
    n2 <- df$Control.N
    bi <- n1 - ai
    di <- n2 - ci

    valid <- !is.na(ai) & !is.na(bi) & !is.na(ci) & !is.na(di) &
      ai >= 0 & bi >= 0 & ci >= 0 & di >= 0

    missing_year <- sum(valid & is.na(df$Study.year))
    df_valid <- df[valid & !is.na(df$Study.year), , drop = FALSE]
    if (nrow(df_valid) < 2) {
      return(list(es = data.frame(), missing_year = missing_year, valid_n = sum(valid)))
    }

    es <- escalc(
      measure = method$measure,
      ai = df_valid$Experimental.cases,
      bi = df_valid$Experimental.N - df_valid$Experimental.cases,
      ci = df_valid$Control.cases,
      di = df_valid$Control.N - df_valid$Control.cases,
      add = 0.5,
      to = "only0"
    )

    es_df <- data.frame(
      yi = es$yi,
      vi = es$vi,
      year = df_valid$Study.year
    )
  } else if (type == "continuous") {
    valid <- !is.na(df$Experimental.mean) & !is.na(df$Experimental.SD) &
      !is.na(df$Experimental.N) & !is.na(df$Control.mean) &
      !is.na(df$Control.SD) & !is.na(df$Control.N)

    missing_year <- sum(valid & is.na(df$Study.year))
    df_valid <- df[valid & !is.na(df$Study.year), , drop = FALSE]
    if (nrow(df_valid) < 2) {
      return(list(es = data.frame(), missing_year = missing_year, valid_n = sum(valid)))
    }

    es <- escalc(
      measure = method$measure,
      m1i = df_valid$Experimental.mean,
      sd1i = df_valid$Experimental.SD,
      n1i = df_valid$Experimental.N,
      m2i = df_valid$Control.mean,
      sd2i = df_valid$Control.SD,
      n2i = df_valid$Control.N
    )

    es_df <- data.frame(
      yi = es$yi,
      vi = es$vi,
      year = df_valid$Study.year
    )
  } else if (type == "giv") {
    valid <- !is.na(df$GIV.Mean) & !is.na(df$GIV.SE) & df$GIV.SE > 0
    missing_year <- sum(valid & is.na(df$Study.year))
    df_valid <- df[valid & !is.na(df$Study.year), , drop = FALSE]
    if (nrow(df_valid) < 2) {
      return(list(es = data.frame(), missing_year = missing_year, valid_n = sum(valid)))
    }
    es_df <- data.frame(
      yi = df_valid$GIV.Mean,
      vi = df_valid$GIV.SE^2,
      year = df_valid$Study.year
    )
  } else if (type == "mean_ci") {
    valid <- !is.na(df$Mean) & !is.na(df$CI.start) & !is.na(df$CI.end)
    if (effect_is_log) {
      valid <- valid & df$Mean > 0 & df$CI.start > 0 & df$CI.end > 0
    }
    missing_year <- sum(valid & is.na(df$Study.year))
    df_valid <- df[valid & !is.na(df$Study.year), , drop = FALSE]
    if (nrow(df_valid) < 2) {
      return(list(es = data.frame(), missing_year = missing_year, valid_n = sum(valid)))
    }

    if (effect_is_log) {
      yi <- log(df_valid$Mean)
      sei <- (log(df_valid$CI.end) - log(df_valid$CI.start)) / (2 * 1.96)
    } else {
      yi <- df_valid$Mean
      sei <- (df_valid$CI.end - df_valid$CI.start) / (2 * 1.96)
    }
    vi <- sei^2

    es_df <- data.frame(
      yi = yi,
      vi = vi,
      year = df_valid$Study.year
    )
  } else if (type == "mean_var") {
    valid <- !is.na(df$Mean) & !is.na(df$Variance) & df$Variance > 0
    if (effect_is_log) {
      valid <- valid & df$Mean > 0
    }
    missing_year <- sum(valid & is.na(df$Study.year))
    df_valid <- df[valid & !is.na(df$Study.year), , drop = FALSE]
    if (nrow(df_valid) < 2) {
      return(list(es = data.frame(), missing_year = missing_year, valid_n = sum(valid)))
    }

    yi <- if (effect_is_log) log(df_valid$Mean) else df_valid$Mean
    vi <- df_valid$Variance
    es_df <- data.frame(
      yi = yi,
      vi = vi,
      year = df_valid$Study.year
    )
  } else {
    return(list(es = data.frame(), missing_year = 0, valid_n = 0))
  }

  es_df <- es_df[is.finite(es_df$yi) & is.finite(es_df$vi) & es_df$vi > 0, , drop = FALSE]
  list(es = es_df, missing_year = missing_year, valid_n = sum(valid))
}

fit_meta <- function(es_df, effect_is_log, method = "REML") {
  if (nrow(es_df) < 2) {
    return(list(
      k = nrow(es_df),
      mu = NA_real_,
      se = NA_real_,
      ci.lb = NA_real_,
      ci.ub = NA_real_,
      tau2 = NA_real_,
      i2 = NA_real_,
      q = NA_real_,
      qp = NA_real_,
      pi.lb = NA_real_,
      pi.ub = NA_real_,
      effect = NA_real_,
      effect.lb = NA_real_,
      effect.ub = NA_real_,
      pi.width = NA_real_,
      error = NA_character_
    ))
  }

  res <- tryCatch(
    rma(yi = es_df$yi, vi = es_df$vi, method = method),
    error = function(e) e
  )

  if (inherits(res, "error")) {
    return(list(
      k = nrow(es_df),
      mu = NA_real_,
      se = NA_real_,
      ci.lb = NA_real_,
      ci.ub = NA_real_,
      tau2 = NA_real_,
      i2 = NA_real_,
      q = NA_real_,
      qp = NA_real_,
      pi.lb = NA_real_,
      pi.ub = NA_real_,
      effect = NA_real_,
      effect.lb = NA_real_,
      effect.ub = NA_real_,
      pi.width = NA_real_,
      error = res$message
    ))
  }

  pred <- tryCatch(predict(res), error = function(e) NULL)
  mu <- as.numeric(res$b)
  ci.lb <- as.numeric(res$ci.lb)
  ci.ub <- as.numeric(res$ci.ub)
  pi.lb <- if (!is.null(pred) && !is.null(pred$pi.lb)) pred$pi.lb else NA_real_
  pi.ub <- if (!is.null(pred) && !is.null(pred$pi.ub)) pred$pi.ub else NA_real_
  if (length(pi.lb) == 0) {
    pi.lb <- NA_real_
  }
  if (length(pi.ub) == 0) {
    pi.ub <- NA_real_
  }
  pi.lb <- as.numeric(pi.lb[1])
  pi.ub <- as.numeric(pi.ub[1])

  if (effect_is_log) {
    effect <- exp(mu)
    effect.lb <- exp(ci.lb)
    effect.ub <- exp(ci.ub)
    if (!is.na(pi.lb) && !is.na(pi.ub)) {
      pi.lb <- exp(pi.lb)
      pi.ub <- exp(pi.ub)
    }
  } else {
    effect <- mu
    effect.lb <- ci.lb
    effect.ub <- ci.ub
  }

  pi.width <- if (!is.na(pi.lb) && !is.na(pi.ub)) pi.ub - pi.lb else NA_real_

  list(
    k = nrow(es_df),
    mu = mu,
    se = as.numeric(res$se),
    ci.lb = ci.lb,
    ci.ub = ci.ub,
    tau2 = as.numeric(res$tau2),
    i2 = as.numeric(res$I2),
    q = as.numeric(res$QE),
    qp = as.numeric(res$QEp),
    pi.lb = pi.lb,
    pi.ub = pi.ub,
    effect = effect,
    effect.lb = effect.lb,
    effect.ub = effect.ub,
    pi.width = pi.width,
    error = NA_character_
  )
}

sig_flag <- function(meta) {
  if (is.na(meta$ci.lb) || is.na(meta$ci.ub)) {
    return(NA)
  }
  meta$ci.lb > 0 || meta$ci.ub < 0
}

sign_flag <- function(meta) {
  if (is.na(meta$mu)) {
    return(NA_integer_)
  }
  if (meta$mu > 0) {
    return(1L)
  }
  if (meta$mu < 0) {
    return(-1L)
  }
  0L
}

compute_influence_metrics <- function(es_df, cutoff_year = 2005) {
  if (nrow(es_df) < 3) {
    return(list(
      cook_thresh = NA_real_,
      n_influential = NA_integer_,
      max_cook = NA_real_,
      prop_influ_pre = NA_real_,
      prop_pre_overall = NA_real_,
      mean_cook_pre = NA_real_,
      mean_cook_post = NA_real_
    ))
  }

  res <- tryCatch(
    rma(yi = es_df$yi, vi = es_df$vi, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(res)) {
    return(list(
      cook_thresh = NA_real_,
      n_influential = NA_integer_,
      max_cook = NA_real_,
      prop_influ_pre = NA_real_,
      prop_pre_overall = NA_real_,
      mean_cook_pre = NA_real_,
      mean_cook_post = NA_real_
    ))
  }

  infl <- tryCatch(influence(res), error = function(e) NULL)
  if (is.null(infl)) {
    return(list(
      cook_thresh = NA_real_,
      n_influential = NA_integer_,
      max_cook = NA_real_,
      prop_influ_pre = NA_real_,
      prop_pre_overall = NA_real_,
      mean_cook_pre = NA_real_,
      mean_cook_post = NA_real_
    ))
  }

  cook <- infl$inf$cook.d
  if (is.null(cook)) {
    return(list(
      cook_thresh = NA_real_,
      n_influential = NA_integer_,
      max_cook = NA_real_,
      prop_influ_pre = NA_real_,
      prop_pre_overall = NA_real_,
      mean_cook_pre = NA_real_,
      mean_cook_post = NA_real_
    ))
  }
  k <- length(cook)
  cook_thresh <- 4 / k
  influ_idx <- which(cook > cook_thresh)
  n_influential <- length(influ_idx)

  years <- es_df$year
  prop_pre_overall <- mean(years < cutoff_year, na.rm = TRUE)
  prop_influ_pre <- if (n_influential > 0) {
    mean(years[influ_idx] < cutoff_year, na.rm = TRUE)
  } else {
    NA_real_
  }

  mean_cook_pre <- if (any(years < cutoff_year)) {
    mean(cook[years < cutoff_year], na.rm = TRUE)
  } else {
    NA_real_
  }
  mean_cook_post <- if (any(years >= cutoff_year)) {
    mean(cook[years >= cutoff_year], na.rm = TRUE)
  } else {
    NA_real_
  }

  list(
    cook_thresh = cook_thresh,
    n_influential = n_influential,
    max_cook = max(cook, na.rm = TRUE),
    prop_influ_pre = prop_influ_pre,
    prop_pre_overall = prop_pre_overall,
    mean_cook_pre = mean_cook_pre,
    mean_cook_post = mean_cook_post
  )
}

compute_loo_metrics <- function(es_df, effect_is_log, method = "REML") {
  if (nrow(es_df) < 3) {
    return(list(
      max_abs_delta_effect = NA_real_,
      max_abs_delta_i2 = NA_real_,
      max_abs_delta_tau2 = NA_real_
    ))
  }

  res <- tryCatch(
    rma(yi = es_df$yi, vi = es_df$vi, method = method),
    error = function(e) NULL
  )
  if (is.null(res)) {
    return(list(
      max_abs_delta_effect = NA_real_,
      max_abs_delta_i2 = NA_real_,
      max_abs_delta_tau2 = NA_real_
    ))
  }

  loo <- tryCatch(leave1out(res), error = function(e) NULL)
  if (is.null(loo)) {
    return(list(
      max_abs_delta_effect = NA_real_,
      max_abs_delta_i2 = NA_real_,
      max_abs_delta_tau2 = NA_real_
    ))
  }

  full_effect <- if (effect_is_log) exp(as.numeric(res$b)) else as.numeric(res$b)
  loo_effect <- if (effect_is_log) exp(loo$estimate) else loo$estimate

  list(
    max_abs_delta_effect = max(abs(loo_effect - full_effect), na.rm = TRUE),
    max_abs_delta_i2 = max(abs(loo$I2 - as.numeric(res$I2)), na.rm = TRUE),
    max_abs_delta_tau2 = max(abs(loo$tau2 - as.numeric(res$tau2)), na.rm = TRUE)
  )
}

read_dataset <- function(path) {
  env <- new.env(parent = emptyenv())
  obj <- load(path, envir = env)
  df <- env[[obj[1]]]
  names(df) <- gsub(" ", ".", names(df))
  df
}

reuse_results <- file.exists(results_path) && Sys.getenv("PAIRWISE70_RECOMPUTE") != "1"
if (reuse_results) {
  results_df <- read.csv(results_path, stringsAsFactors = FALSE)
  issues_df <- if (file.exists(issues_path)) {
    read.csv(issues_path, stringsAsFactors = FALSE)
  } else {
    data.frame()
  }
} else {
  rda_files <- list.files(data_dir, pattern = "\\.rda$", full.names = TRUE)
  results <- list()
  issues <- list()

  for (i in seq_along(rda_files)) {
    path <- rda_files[i]
    dataset_id <- tools::file_path_sans_ext(basename(path))
    df <- read_dataset(path)

    if (!all(c("Analysis.group", "Analysis.number", "Analysis.name") %in% names(df))) {
      issues[[length(issues) + 1]] <- data.frame(
        dataset_id = dataset_id,
        analysis_id = NA_character_,
        issue = "Missing analysis identifiers"
      )
      next
    }

    review_doi <- if ("review_doi" %in% names(df)) {
      unique(na.omit(df$review_doi))[1]
    } else {
      NA_character_
    }

    analysis_keys <- unique(df[, c("Analysis.group", "Analysis.number", "Analysis.name")])

    for (j in seq_len(nrow(analysis_keys))) {
      key <- analysis_keys[j, ]
      df_sub <- df[
        df$Analysis.group == key$Analysis.group &
          df$Analysis.number == key$Analysis.number &
          df$Analysis.name == key$Analysis.name,
        ,
        drop = FALSE
      ]

      method <- choose_method(df_sub, key$Analysis.name)
      if (is.null(method)) {
        issues[[length(issues) + 1]] <- data.frame(
          dataset_id = dataset_id,
          analysis_id = paste(key$Analysis.group, key$Analysis.number, key$Analysis.name, sep = "::"),
          issue = "No usable data for effect sizes"
        )
        next
      }

      es_info <- extract_es(df_sub, method)
      es_df <- es_info$es
      if (nrow(es_df) < 2) {
        next
      }

      for (threshold in thresholds) {
        es_all <- es_df
        es_post <- es_all[es_all$year >= threshold, , drop = FALSE]
        es_pre <- es_all[es_all$year < threshold, , drop = FALSE]

        meta_all <- fit_meta(es_all, method$effect_is_log)
        meta_post <- fit_meta(es_post, method$effect_is_log)

        results[[length(results) + 1]] <- data.frame(
          dataset_id = dataset_id,
          review_doi = review_doi,
          analysis_group = key$Analysis.group,
          analysis_number = key$Analysis.number,
          analysis_name = key$Analysis.name,
          data_type = method$type,
          measure = method$measure,
          effect_is_log = method$effect_is_log,
          threshold = threshold,
          k_all = meta_all$k,
          k_post = meta_post$k,
          k_pre = nrow(es_pre),
          k_missing_year = es_info$missing_year,
          frac_pre = if (meta_all$k > 0) nrow(es_pre) / meta_all$k else NA_real_,
          i2_all = meta_all$i2,
          i2_post = meta_post$i2,
          tau2_all = meta_all$tau2,
          tau2_post = meta_post$tau2,
          q_all = meta_all$q,
          q_post = meta_post$q,
          pi_width_all = meta_all$pi.width,
          pi_width_post = meta_post$pi.width,
          effect_all = meta_all$effect,
          effect_post = meta_post$effect,
          delta_i2 = if (!is.na(meta_all$i2) && !is.na(meta_post$i2)) meta_post$i2 - meta_all$i2 else NA_real_,
          delta_tau2 = if (!is.na(meta_all$tau2) && !is.na(meta_post$tau2)) meta_post$tau2 - meta_all$tau2 else NA_real_,
          delta_pi_width = if (!is.na(meta_all$pi.width) && !is.na(meta_post$pi.width)) meta_post$pi.width - meta_all$pi.width else NA_real_,
          delta_effect = if (!is.na(meta_all$effect) && !is.na(meta_post$effect)) meta_post$effect - meta_all$effect else NA_real_,
          error_all = meta_all$error,
          error_post = meta_post$error
        )
      }
    }

    if (i %% 25 == 0) {
      message("Processed ", i, " / ", length(rda_files), " datasets")
    }
  }

  results_df <- bind_rows(results)
  issues_df <- bind_rows(issues)

  write.csv(results_df, results_path, row.names = FALSE)
  write.csv(issues_df, issues_path, row.names = FALSE)
}

summary_df <- results_df %>%
  group_by(threshold) %>%
  summarize(
    n_analyses = sum(!is.na(delta_i2) & !is.na(delta_tau2)),
    prop_i2_increase = ifelse(n_analyses > 0, mean(delta_i2 > 0, na.rm = TRUE), NA_real_),
    median_delta_i2 = ifelse(n_analyses > 0, median(delta_i2, na.rm = TRUE), NA_real_),
    mean_delta_i2 = ifelse(n_analyses > 0, mean(delta_i2, na.rm = TRUE), NA_real_),
    prop_tau2_increase = ifelse(n_analyses > 0, mean(delta_tau2 > 0, na.rm = TRUE), NA_real_),
    median_delta_tau2 = ifelse(n_analyses > 0, median(delta_tau2, na.rm = TRUE), NA_real_),
    mean_delta_tau2 = ifelse(n_analyses > 0, mean(delta_tau2, na.rm = TRUE), NA_real_),
    median_delta_pi_width = ifelse(n_analyses > 0, median(delta_pi_width, na.rm = TRUE), NA_real_),
    median_delta_effect = ifelse(n_analyses > 0, median(delta_effect, na.rm = TRUE), NA_real_),
    .groups = "drop"
  )

write.csv(summary_df, file.path(out_dir, "threshold_summary.csv"), row.names = FALSE)

main <- results_df %>%
  filter(threshold == 2005, !is.na(delta_i2), !is.na(delta_tau2))

if (nrow(main) > 0) {
  i2_test <- suppressWarnings(wilcox.test(main$delta_i2, mu = 0))
  tau2_test <- suppressWarnings(wilcox.test(main$delta_tau2, mu = 0))

  main <- main %>%
    mutate(i2_increase = delta_i2 > 0)

  glm_fit <- suppressWarnings(glm(i2_increase ~ frac_pre + log1p(k_all), data = main, family = binomial()))
  lm_fit <- suppressWarnings(lm(delta_i2 ~ frac_pre + log1p(k_all), data = main))

  summary_path <- file.path(out_dir, "summary_2005.txt")
  sink(summary_path)
  cat("Pre-2005 impact analysis (threshold 2005)\n")
  cat("Analyses with valid paired metrics:", nrow(main), "\n")
  cat("Proportion with higher I2 after excluding pre-2005:", mean(main$i2_increase), "\n")
  cat("Median delta I2:", median(main$delta_i2), "\n")
  cat("Median delta tau2:", median(main$delta_tau2), "\n")
  cat("Wilcoxon test delta I2 p-value:", i2_test$p.value, "\n")
  cat("Wilcoxon test delta tau2 p-value:", tau2_test$p.value, "\n\n")
  cat("Logistic model (I2 increase ~ frac_pre + log1p(k_all)):\n")
  print(summary(glm_fit))
  cat("\nLinear model (delta I2 ~ frac_pre + log1p(k_all)):\n")
  print(summary(lm_fit))
  sink()

  p1 <- ggplot(main, aes(x = delta_i2)) +
    geom_histogram(bins = 30, fill = "#377eb8", color = "white") +
    theme_minimal() +
    labs(title = "Delta I2 (post-2005 minus all)", x = "Delta I2", y = "Count")
  ggsave(file.path(out_dir, "plots", "delta_i2_2005.png"), p1, width = 7, height = 4)

  p2 <- ggplot(main, aes(x = i2_all, y = i2_post)) +
    geom_point(alpha = 0.4) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
    theme_minimal() +
    labs(title = "I2 all vs post-2005", x = "I2 (all)", y = "I2 (post-2005)")
  ggsave(file.path(out_dir, "plots", "i2_scatter_2005.png"), p2, width = 5, height = 5)

  p3 <- ggplot(main, aes(x = frac_pre, y = delta_i2)) +
    geom_point(alpha = 0.4) +
    geom_smooth(method = "lm", se = TRUE, color = "#e41a1c") +
    theme_minimal() +
    labs(title = "Delta I2 vs fraction pre-2005", x = "Fraction pre-2005", y = "Delta I2")
  ggsave(file.path(out_dir, "plots", "delta_i2_vs_frac_pre_2005.png"), p3, width = 6, height = 4)
}

reuse_sens <- file.exists(sens_path) && Sys.getenv("PAIRWISE70_RECOMPUTE_SENS") != "1"
if (reuse_sens) {
  sens_df <- read.csv(sens_path, stringsAsFactors = FALSE)
} else {
  key_df <- results_df %>%
    filter(threshold == 2005, k_pre > 0, k_post >= 2, !is.na(delta_i2)) %>%
    distinct(dataset_id, analysis_group, analysis_number, analysis_name)

  rda_files_all <- list.files(data_dir, pattern = "\\.rda$", full.names = TRUE)
  rda_map <- setNames(rda_files_all, tools::file_path_sans_ext(basename(rda_files_all)))
  keys_by_dataset <- split(key_df, key_df$dataset_id)

  sens_list <- list()
  dataset_ids <- names(keys_by_dataset)

  for (i in seq_along(dataset_ids)) {
    dataset_id <- dataset_ids[i]
    path <- rda_map[[dataset_id]]
    if (is.null(path)) {
      next
    }

    df <- read_dataset(path)
    key_set <- keys_by_dataset[[dataset_id]]

    for (j in seq_len(nrow(key_set))) {
      key <- key_set[j, ]
      df_sub <- df[
        df$Analysis.group == key$analysis_group &
          df$Analysis.number == key$analysis_number &
          df$Analysis.name == key$analysis_name,
        ,
        drop = FALSE
      ]

      method <- choose_method(df_sub, key$analysis_name)
      if (is.null(method)) {
        next
      }

      es_info <- extract_es(df_sub, method)
      es_df <- es_info$es
      if (nrow(es_df) < 2) {
        next
      }

      es_post <- es_df[es_df$year >= 2005, , drop = FALSE]
      es_pre <- es_df[es_df$year < 2005, , drop = FALSE]
      if (nrow(es_post) < 2) {
        next
      }

      meta_re_all <- fit_meta(es_df, method$effect_is_log, method = "REML")
      meta_re_post <- fit_meta(es_post, method$effect_is_log, method = "REML")
      meta_re_pre <- fit_meta(es_pre, method$effect_is_log, method = "REML")
      meta_fe_all <- fit_meta(es_df, method$effect_is_log, method = "FE")
      meta_fe_post <- fit_meta(es_post, method$effect_is_log, method = "FE")
      meta_fe_pre <- fit_meta(es_pre, method$effect_is_log, method = "FE")

      loo_all <- compute_loo_metrics(es_df, method$effect_is_log, method = "REML")
      loo_post <- compute_loo_metrics(es_post, method$effect_is_log, method = "REML")
      infl_all <- compute_influence_metrics(es_df, cutoff_year = 2005)

      sig_all <- sig_flag(meta_re_all)
      sig_post <- sig_flag(meta_re_post)
      sig_pre <- sig_flag(meta_re_pre)
      sign_all <- sign_flag(meta_re_all)
      sign_post <- sign_flag(meta_re_post)
      sign_pre <- sign_flag(meta_re_pre)

      sens_list[[length(sens_list) + 1]] <- data.frame(
        dataset_id = dataset_id,
        analysis_group = key$analysis_group,
        analysis_number = key$analysis_number,
        analysis_name = key$analysis_name,
        data_type = method$type,
        measure = method$measure,
        effect_is_log = method$effect_is_log,
        k_all = meta_re_all$k,
        k_post = meta_re_post$k,
        k_pre = nrow(es_pre),
        frac_pre = if (meta_re_all$k > 0) nrow(es_pre) / meta_re_all$k else NA_real_,
        i2_re_all = meta_re_all$i2,
        i2_re_post = meta_re_post$i2,
        i2_re_pre = meta_re_pre$i2,
        delta_i2_re = if (!is.na(meta_re_all$i2) && !is.na(meta_re_post$i2)) meta_re_post$i2 - meta_re_all$i2 else NA_real_,
        delta_i2_post_pre = if (!is.na(meta_re_pre$i2) && !is.na(meta_re_post$i2)) meta_re_post$i2 - meta_re_pre$i2 else NA_real_,
        tau2_re_all = meta_re_all$tau2,
        tau2_re_post = meta_re_post$tau2,
        tau2_re_pre = meta_re_pre$tau2,
        delta_tau2_re = if (!is.na(meta_re_all$tau2) && !is.na(meta_re_post$tau2)) meta_re_post$tau2 - meta_re_all$tau2 else NA_real_,
        delta_tau2_post_pre = if (!is.na(meta_re_pre$tau2) && !is.na(meta_re_post$tau2)) meta_re_post$tau2 - meta_re_pre$tau2 else NA_real_,
        effect_re_all = meta_re_all$effect,
        effect_re_post = meta_re_post$effect,
        effect_re_pre = meta_re_pre$effect,
        delta_effect_re = if (!is.na(meta_re_all$effect) && !is.na(meta_re_post$effect)) meta_re_post$effect - meta_re_all$effect else NA_real_,
        delta_effect_post_pre = if (!is.na(meta_re_pre$effect) && !is.na(meta_re_post$effect)) meta_re_post$effect - meta_re_pre$effect else NA_real_,
        effect_fe_all = meta_fe_all$effect,
        effect_fe_post = meta_fe_post$effect,
        effect_fe_pre = meta_fe_pre$effect,
        delta_effect_fe = if (!is.na(meta_fe_all$effect) && !is.na(meta_fe_post$effect)) meta_fe_post$effect - meta_fe_all$effect else NA_real_,
        delta_effect_fe_post_pre = if (!is.na(meta_fe_pre$effect) && !is.na(meta_fe_post$effect)) meta_fe_post$effect - meta_fe_pre$effect else NA_real_,
        sig_all = sig_all,
        sig_post = sig_post,
        sig_pre = sig_pre,
        sign_all = sign_all,
        sign_post = sign_post,
        sign_pre = sign_pre,
        max_abs_delta_effect_loo_all = loo_all$max_abs_delta_effect,
        max_abs_delta_effect_loo_post = loo_post$max_abs_delta_effect,
        max_abs_delta_i2_loo_all = loo_all$max_abs_delta_i2,
        max_abs_delta_i2_loo_post = loo_post$max_abs_delta_i2,
        max_abs_delta_tau2_loo_all = loo_all$max_abs_delta_tau2,
        max_abs_delta_tau2_loo_post = loo_post$max_abs_delta_tau2,
        cook_thresh = infl_all$cook_thresh,
        n_influential = infl_all$n_influential,
        max_cook = infl_all$max_cook,
        prop_influ_pre = infl_all$prop_influ_pre,
        prop_pre_overall = infl_all$prop_pre_overall,
        mean_cook_pre = infl_all$mean_cook_pre,
        mean_cook_post = infl_all$mean_cook_post,
        stringsAsFactors = FALSE
      )
    }

    if (i %% 25 == 0) {
      message("Sensitivity processed ", i, " / ", length(dataset_ids), " datasets")
    }
  }

  sens_df <- bind_rows(sens_list)
  write.csv(sens_df, sens_path, row.names = FALSE)
}

if (exists("sens_df") && nrow(sens_df) > 0) {
  summarize_sens <- function(df, min_k_post) {
    df_sub <- df[df$k_post >= min_k_post, , drop = FALSE]
    n <- nrow(df_sub)
    df_influ <- df_sub[!is.na(df_sub$prop_influ_pre), , drop = FALSE]
    n_influ <- nrow(df_influ)
    data.frame(
      min_k_post = min_k_post,
      n_analyses = n,
      n_with_influential = n_influ,
      prop_delta_i2_gt0 = if (n > 0) mean(df_sub$delta_i2_re > 0, na.rm = TRUE) else NA_real_,
      median_delta_i2 = if (n > 0) median(df_sub$delta_i2_re, na.rm = TRUE) else NA_real_,
      mean_delta_i2 = if (n > 0) mean(df_sub$delta_i2_re, na.rm = TRUE) else NA_real_,
      prop_delta_tau2_gt0 = if (n > 0) mean(df_sub$delta_tau2_re > 0, na.rm = TRUE) else NA_real_,
      median_delta_tau2 = if (n > 0) median(df_sub$delta_tau2_re, na.rm = TRUE) else NA_real_,
      mean_delta_tau2 = if (n > 0) mean(df_sub$delta_tau2_re, na.rm = TRUE) else NA_real_,
      median_delta_effect_re = if (n > 0) median(df_sub$delta_effect_re, na.rm = TRUE) else NA_real_,
      median_delta_effect_fe = if (n > 0) median(df_sub$delta_effect_fe, na.rm = TRUE) else NA_real_,
      median_max_abs_delta_effect_loo_all = if (n > 0) median(df_sub$max_abs_delta_effect_loo_all, na.rm = TRUE) else NA_real_,
      median_max_abs_delta_effect_loo_post = if (n > 0) median(df_sub$max_abs_delta_effect_loo_post, na.rm = TRUE) else NA_real_,
      median_max_abs_delta_i2_loo_all = if (n > 0) median(df_sub$max_abs_delta_i2_loo_all, na.rm = TRUE) else NA_real_,
      median_max_abs_delta_i2_loo_post = if (n > 0) median(df_sub$max_abs_delta_i2_loo_post, na.rm = TRUE) else NA_real_,
      mean_mean_cook_pre = if (n > 0) mean(df_sub$mean_cook_pre, na.rm = TRUE) else NA_real_,
      mean_mean_cook_post = if (n > 0) mean(df_sub$mean_cook_post, na.rm = TRUE) else NA_real_,
      mean_mean_cook_delta = if (n > 0) mean(df_sub$mean_cook_pre - df_sub$mean_cook_post, na.rm = TRUE) else NA_real_,
      mean_prop_influ_pre = if (n_influ > 0) mean(df_influ$prop_influ_pre, na.rm = TRUE) else NA_real_,
      mean_prop_pre_overall = if (n_influ > 0) mean(df_influ$prop_pre_overall, na.rm = TRUE) else NA_real_,
      mean_prop_influ_minus_overall = if (n_influ > 0) mean(df_influ$prop_influ_pre - df_influ$prop_pre_overall, na.rm = TRUE) else NA_real_,
      prop_influ_pre_gt_overall = if (n_influ > 0) mean((df_influ$prop_influ_pre - df_influ$prop_pre_overall) > 0, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  sens_summary_df <- bind_rows(
    summarize_sens(sens_df, 2),
    summarize_sens(sens_df, 5),
    summarize_sens(sens_df, 10)
  )

  write.csv(sens_summary_df, sens_summary_csv, row.names = FALSE)

  sink(sens_summary_path)
  cat("Sensitivity analysis (fixed vs random, influence, leave-one-out)\n")
  cat("Analyses limited to those with pre-2005 data and k_post >= 2\n\n")
  print(sens_summary_df)
  sink()
}

if (exists("sens_df") && nrow(sens_df) > 0) {
  summarize_extended <- function(df, min_k_post) {
    df_sub <- df[df$k_post >= min_k_post, , drop = FALSE]
    n <- nrow(df_sub)
    df_sig <- df_sub[!is.na(df_sub$sig_all) & !is.na(df_sub$sig_post), , drop = FALSE]
    df_sign <- df_sub[!is.na(df_sub$sign_all) & !is.na(df_sub$sign_post), , drop = FALSE]
    df_pre_ok <- df_sub[df_sub$k_pre >= 2, , drop = FALSE]

    data.frame(
      min_k_post = min_k_post,
      n_analyses = n,
      prop_delta_i2_gt0 = if (n > 0) mean(df_sub$delta_i2_re > 0, na.rm = TRUE) else NA_real_,
      prop_abs_delta_i2_le_5 = if (n > 0) mean(abs(df_sub$delta_i2_re) <= 5, na.rm = TRUE) else NA_real_,
      prop_abs_delta_i2_le_10 = if (n > 0) mean(abs(df_sub$delta_i2_re) <= 10, na.rm = TRUE) else NA_real_,
      prop_delta_tau2_gt0 = if (n > 0) mean(df_sub$delta_tau2_re > 0, na.rm = TRUE) else NA_real_,
      median_delta_i2 = if (n > 0) median(df_sub$delta_i2_re, na.rm = TRUE) else NA_real_,
      mean_delta_i2 = if (n > 0) mean(df_sub$delta_i2_re, na.rm = TRUE) else NA_real_,
      median_delta_tau2 = if (n > 0) median(df_sub$delta_tau2_re, na.rm = TRUE) else NA_real_,
      mean_delta_tau2 = if (n > 0) mean(df_sub$delta_tau2_re, na.rm = TRUE) else NA_real_,
      median_delta_effect_re = if (n > 0) median(df_sub$delta_effect_re, na.rm = TRUE) else NA_real_,
      median_delta_effect_fe = if (n > 0) median(df_sub$delta_effect_fe, na.rm = TRUE) else NA_real_,
      prop_sig_change_all_post = if (nrow(df_sig) > 0) mean(df_sig$sig_all != df_sig$sig_post, na.rm = TRUE) else NA_real_,
      prop_sig_loss_all_post = if (nrow(df_sig) > 0) mean(df_sig$sig_all & !df_sig$sig_post, na.rm = TRUE) else NA_real_,
      prop_sig_gain_all_post = if (nrow(df_sig) > 0) mean(!df_sig$sig_all & df_sig$sig_post, na.rm = TRUE) else NA_real_,
      prop_sign_flip_all_post = if (nrow(df_sign) > 0) mean(df_sign$sign_all * df_sign$sign_post < 0, na.rm = TRUE) else NA_real_,
      n_pre_ok = nrow(df_pre_ok),
      median_delta_i2_post_pre = if (nrow(df_pre_ok) > 0) median(df_pre_ok$delta_i2_post_pre, na.rm = TRUE) else NA_real_,
      mean_delta_i2_post_pre = if (nrow(df_pre_ok) > 0) mean(df_pre_ok$delta_i2_post_pre, na.rm = TRUE) else NA_real_,
      median_delta_tau2_post_pre = if (nrow(df_pre_ok) > 0) median(df_pre_ok$delta_tau2_post_pre, na.rm = TRUE) else NA_real_,
      mean_delta_tau2_post_pre = if (nrow(df_pre_ok) > 0) mean(df_pre_ok$delta_tau2_post_pre, na.rm = TRUE) else NA_real_,
      median_delta_effect_post_pre = if (nrow(df_pre_ok) > 0) median(df_pre_ok$delta_effect_post_pre, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  extended_summary_df <- bind_rows(
    summarize_extended(sens_df, 2),
    summarize_extended(sens_df, 5),
    summarize_extended(sens_df, 10)
  )

  write.csv(extended_summary_df, extended_summary_csv, row.names = FALSE)

  sens_grouped <- sens_df %>%
    mutate(
      k_all_bin = cut(k_all, breaks = c(-Inf, 5, 10, 20, 50, Inf),
                      labels = c("<=5", "6-10", "11-20", "21-50", ">50")),
      frac_pre_bin = cut(frac_pre, breaks = c(0, 0.25, 0.5, 0.75, 1),
                         include.lowest = TRUE)
    )

  summarize_group <- function(df, min_k_post, group_type, group_vars) {
    df_sub <- df[df$k_post >= min_k_post, , drop = FALSE]
    if (nrow(df_sub) == 0) {
      return(data.frame())
    }
    df_sub %>%
      group_by(across(all_of(group_vars))) %>%
      summarize(
        min_k_post = min_k_post,
        n_analyses = n(),
        prop_delta_i2_gt0 = mean(delta_i2_re > 0, na.rm = TRUE),
        median_delta_i2 = median(delta_i2_re, na.rm = TRUE),
        prop_sig_change_all_post = mean(sig_all != sig_post, na.rm = TRUE),
        prop_sign_flip_all_post = mean(sign_all * sign_post < 0, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        group_type = group_type,
        group_value = do.call(paste, c(across(all_of(group_vars)), sep = " | "))
      ) %>%
      select(group_type, group_value, min_k_post, n_analyses, prop_delta_i2_gt0,
             median_delta_i2, prop_sig_change_all_post, prop_sign_flip_all_post)
  }

  group_summary <- bind_rows(
    summarize_group(sens_grouped, 2, "data_type", c("data_type")),
    summarize_group(sens_grouped, 2, "measure", c("measure")),
    summarize_group(sens_grouped, 2, "data_type+measure", c("data_type", "measure")),
    summarize_group(sens_grouped, 2, "k_all_bin", c("k_all_bin")),
    summarize_group(sens_grouped, 2, "frac_pre_bin", c("frac_pre_bin"))
  ) %>%
    filter(n_analyses >= 30)

  write.csv(group_summary, extended_by_group_csv, row.names = FALSE)

  sink(extended_summary_path)
  cat("Extended 2005 summary (significance changes, direction flips, pre vs post)\n\n")
  print(extended_summary_df)
  cat("\nGroup summaries (filtered to n >= 30):\n")
  print(group_summary)
  sink()

  summarize_type <- function(df, group_vars, group_type_label) {
    df_sub <- df[df$k_post >= 2, , drop = FALSE]
    if (nrow(df_sub) == 0) {
      return(data.frame())
    }
    df_sub %>%
      group_by(across(all_of(group_vars))) %>%
      summarize(
        n_analyses = n(),
        prop_delta_i2_gt0 = mean(delta_i2_re > 0, na.rm = TRUE),
        prop_abs_delta_i2_le_5 = mean(abs(delta_i2_re) <= 5, na.rm = TRUE),
        prop_abs_delta_i2_le_10 = mean(abs(delta_i2_re) <= 10, na.rm = TRUE),
        median_delta_i2 = median(delta_i2_re, na.rm = TRUE),
        mean_delta_i2 = mean(delta_i2_re, na.rm = TRUE),
        prop_delta_tau2_gt0 = mean(delta_tau2_re > 0, na.rm = TRUE),
        median_delta_tau2 = median(delta_tau2_re, na.rm = TRUE),
        mean_delta_tau2 = mean(delta_tau2_re, na.rm = TRUE),
        median_delta_effect_re = median(delta_effect_re, na.rm = TRUE),
        median_delta_effect_fe = median(delta_effect_fe, na.rm = TRUE),
        prop_sig_change = mean(sig_all != sig_post, na.rm = TRUE),
        prop_sig_loss = mean(sig_all & !sig_post, na.rm = TRUE),
        prop_sig_gain = mean(!sig_all & sig_post, na.rm = TRUE),
        prop_sign_flip = mean(sign_all * sign_post < 0, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        group_type = group_type_label,
        group_value = do.call(paste, c(across(all_of(group_vars)), sep = " | "))
      ) %>%
      select(group_type, group_value, everything())
  }

  type_summary_all <- bind_rows(
    summarize_type(sens_df, c("data_type"), "data_type"),
    summarize_type(sens_df, c("measure"), "measure"),
    summarize_type(sens_df, c("data_type", "measure"), "data_type+measure")
  )

  write.csv(type_summary_all, analysis_type_summary_all_path, row.names = FALSE)

  type_summary_df <- type_summary_all %>%
    filter(n_analyses >= 30)

  write.csv(type_summary_df, analysis_type_summary_path, row.names = FALSE)

  category_df <- sens_df %>%
    mutate(outcome_category = vapply(analysis_name, categorize_analysis_name, character(1)))

  category_summary_all <- category_df %>%
    filter(k_post >= 2) %>%
    group_by(outcome_category) %>%
    summarize(
      n_analyses = n(),
      prop_delta_i2_gt0 = mean(delta_i2_re > 0, na.rm = TRUE),
      prop_abs_delta_i2_le_5 = mean(abs(delta_i2_re) <= 5, na.rm = TRUE),
      prop_abs_delta_i2_le_10 = mean(abs(delta_i2_re) <= 10, na.rm = TRUE),
      median_delta_i2 = median(delta_i2_re, na.rm = TRUE),
      mean_delta_i2 = mean(delta_i2_re, na.rm = TRUE),
      prop_sig_change = mean(sig_all != sig_post, na.rm = TRUE),
      prop_sign_flip = mean(sign_all * sign_post < 0, na.rm = TRUE),
      .groups = "drop"
    )

  write.csv(category_summary_all, category_summary_all_path, row.names = FALSE)

  category_summary <- category_summary_all %>%
    filter(n_analyses >= 30)
  write.csv(category_summary, category_summary_path, row.names = FALSE)

  category_df_lr <- category_df %>%
    mutate(measure_group = mapply(measure_group, measure, effect_is_log)) %>%
    filter(measure_group == "log_ratio")

  category_summary_lr_all <- category_df_lr %>%
    filter(k_post >= 2) %>%
    group_by(outcome_category) %>%
    summarize(
      n_analyses = n(),
      prop_delta_i2_gt0 = mean(delta_i2_re > 0, na.rm = TRUE),
      prop_abs_delta_i2_le_5 = mean(abs(delta_i2_re) <= 5, na.rm = TRUE),
      prop_abs_delta_i2_le_10 = mean(abs(delta_i2_re) <= 10, na.rm = TRUE),
      median_delta_i2 = median(delta_i2_re, na.rm = TRUE),
      mean_delta_i2 = mean(delta_i2_re, na.rm = TRUE),
      prop_sig_change = mean(sig_all != sig_post, na.rm = TRUE),
      prop_sign_flip = mean(sign_all * sign_post < 0, na.rm = TRUE),
      .groups = "drop"
    )

  write.csv(category_summary_lr_all, category_summary_log_ratio_all_path, row.names = FALSE)

  category_summary_lr <- category_summary_lr_all %>%
    filter(n_analyses >= 30)
  write.csv(category_summary_lr, category_summary_log_ratio_path, row.names = FALSE)

  category_df_md <- category_df %>%
    mutate(measure_group = mapply(measure_group, measure, effect_is_log)) %>%
    filter(measure_group == "mean_diff")

  category_summary_md_all <- category_df_md %>%
    filter(k_post >= 2) %>%
    group_by(outcome_category) %>%
    summarize(
      n_analyses = n(),
      prop_delta_i2_gt0 = mean(delta_i2_re > 0, na.rm = TRUE),
      prop_abs_delta_i2_le_5 = mean(abs(delta_i2_re) <= 5, na.rm = TRUE),
      prop_abs_delta_i2_le_10 = mean(abs(delta_i2_re) <= 10, na.rm = TRUE),
      median_delta_i2 = median(delta_i2_re, na.rm = TRUE),
      mean_delta_i2 = mean(delta_i2_re, na.rm = TRUE),
      prop_sig_change = mean(sig_all != sig_post, na.rm = TRUE),
      prop_sign_flip = mean(sign_all * sign_post < 0, na.rm = TRUE),
      .groups = "drop"
    )

  write.csv(category_summary_md_all, category_summary_mean_diff_all_path, row.names = FALSE)

  category_summary_md <- category_summary_md_all %>%
    filter(n_analyses >= 30)
  write.csv(category_summary_md, category_summary_mean_diff_path, row.names = FALSE)
}

reuse_sweep <- file.exists(sweep_results_path) && Sys.getenv("PAIRWISE70_RECOMPUTE_SWEEP") != "1"
if (reuse_sweep) {
  sweep_df <- read.csv(sweep_results_path, stringsAsFactors = FALSE)
} else {
  rda_files <- list.files(data_dir, pattern = "\\.rda$", full.names = TRUE)
  existing_df <- if (file.exists(sweep_results_path)) {
    read.csv(sweep_results_path, stringsAsFactors = FALSE)
  } else {
    data.frame()
  }
  existing_cutoffs <- if (nrow(existing_df) > 0) unique(existing_df$cutoff) else integer()
  cutoffs_to_run <- setdiff(sweep_cutoffs, existing_cutoffs)

  if (length(cutoffs_to_run) == 0) {
    sweep_df <- existing_df
  } else {
    existing_keys <- if (nrow(existing_df) > 0) {
      subset_df <- existing_df[existing_df$cutoff %in% cutoffs_to_run, , drop = FALSE]
      if (nrow(subset_df) > 0) {
        paste(subset_df$dataset_id, subset_df$analysis_group, subset_df$analysis_number,
              subset_df$analysis_name, subset_df$cutoff, sep = "|")
      } else {
        character()
      }
    } else {
      character()
    }

    sweep_has_header <- file.exists(sweep_results_path)

    for (i in seq_along(rda_files)) {
      path <- rda_files[i]
      dataset_id <- tools::file_path_sans_ext(basename(path))
      df <- read_dataset(path)

      if (!all(c("Analysis.group", "Analysis.number", "Analysis.name") %in% names(df))) {
        next
      }

      review_doi <- if ("review_doi" %in% names(df)) {
        unique(na.omit(df$review_doi))[1]
      } else {
        NA_character_
      }

      analysis_keys <- unique(df[, c("Analysis.group", "Analysis.number", "Analysis.name")])

      sweep_list <- list()
      for (j in seq_len(nrow(analysis_keys))) {
        key <- analysis_keys[j, ]
        df_sub <- df[
          df$Analysis.group == key$Analysis.group &
            df$Analysis.number == key$Analysis.number &
            df$Analysis.name == key$Analysis.name,
          ,
          drop = FALSE
        ]

        method <- choose_method(df_sub, key$Analysis.name)
        if (is.null(method)) {
          next
        }

        es_info <- extract_es(df_sub, method)
        es_df <- es_info$es
        if (nrow(es_df) < 2) {
          next
        }

        meta_all <- fit_meta(es_df, method$effect_is_log, method = "REML")
        sig_all <- sig_flag(meta_all)
        sign_all <- sign_flag(meta_all)

        for (cutoff in cutoffs_to_run) {
          es_post <- es_df[es_df$year >= cutoff, , drop = FALSE]
          es_pre <- es_df[es_df$year < cutoff, , drop = FALSE]
          if (nrow(es_post) < 2 || nrow(es_pre) < 1) {
            next
          }

          meta_post <- fit_meta(es_post, method$effect_is_log, method = "REML")
          sig_post <- sig_flag(meta_post)
          sign_post <- sign_flag(meta_post)

          sweep_list[[length(sweep_list) + 1]] <- data.frame(
            dataset_id = dataset_id,
            review_doi = review_doi,
            analysis_group = key$Analysis.group,
            analysis_number = key$Analysis.number,
            analysis_name = key$Analysis.name,
            data_type = method$type,
            measure = method$measure,
            effect_is_log = method$effect_is_log,
            cutoff = cutoff,
            k_all = meta_all$k,
            k_post = meta_post$k,
            k_pre = nrow(es_pre),
            frac_pre = if (meta_all$k > 0) nrow(es_pre) / meta_all$k else NA_real_,
            i2_all = meta_all$i2,
            i2_post = meta_post$i2,
            tau2_all = meta_all$tau2,
            tau2_post = meta_post$tau2,
            effect_all = meta_all$effect,
            effect_post = meta_post$effect,
            delta_i2 = if (!is.na(meta_all$i2) && !is.na(meta_post$i2)) meta_post$i2 - meta_all$i2 else NA_real_,
            delta_tau2 = if (!is.na(meta_all$tau2) && !is.na(meta_post$tau2)) meta_post$tau2 - meta_all$tau2 else NA_real_,
            delta_effect = if (!is.na(meta_all$effect) && !is.na(meta_post$effect)) meta_post$effect - meta_all$effect else NA_real_,
            sig_all = sig_all,
            sig_post = sig_post,
            sign_all = sign_all,
            sign_post = sign_post,
            stringsAsFactors = FALSE
          )
        }
      }

      if (length(sweep_list) > 0) {
        sweep_df_chunk <- bind_rows(sweep_list)
        if (length(existing_keys) > 0) {
          chunk_keys <- paste(
            sweep_df_chunk$dataset_id,
            sweep_df_chunk$analysis_group,
            sweep_df_chunk$analysis_number,
            sweep_df_chunk$analysis_name,
            sweep_df_chunk$cutoff,
            sep = "|"
          )
          sweep_df_chunk <- sweep_df_chunk[!chunk_keys %in% existing_keys, , drop = FALSE]
        }
        if (nrow(sweep_df_chunk) > 0) {
          write.table(
            sweep_df_chunk,
            sweep_results_path,
            sep = ",",
            row.names = FALSE,
            col.names = !sweep_has_header,
            append = sweep_has_header
          )
          sweep_has_header <- TRUE
        }
      }

      if (i %% 25 == 0) {
        message("Cutoff sweep processed ", i, " / ", length(rda_files), " datasets")
      }
    }

    sweep_df <- if (file.exists(sweep_results_path)) {
      read.csv(sweep_results_path, stringsAsFactors = FALSE)
    } else {
      data.frame()
    }
  }
}

if (exists("sweep_df") && nrow(sweep_df) > 0) {
  sweep_summary_df <- sweep_df %>%
    group_by(cutoff) %>%
    summarize(
      n_analyses = n(),
      prop_delta_i2_gt0 = mean(delta_i2 > 0, na.rm = TRUE),
      prop_abs_delta_i2_le_5 = mean(abs(delta_i2) <= 5, na.rm = TRUE),
      prop_abs_delta_i2_le_10 = mean(abs(delta_i2) <= 10, na.rm = TRUE),
      median_delta_i2 = median(delta_i2, na.rm = TRUE),
      mean_delta_i2 = mean(delta_i2, na.rm = TRUE),
      prop_delta_tau2_gt0 = mean(delta_tau2 > 0, na.rm = TRUE),
      median_delta_tau2 = median(delta_tau2, na.rm = TRUE),
      mean_delta_tau2 = mean(delta_tau2, na.rm = TRUE),
      prop_sig_change = mean(sig_all != sig_post, na.rm = TRUE),
      prop_sig_loss = mean(sig_all & !sig_post, na.rm = TRUE),
      prop_sig_gain = mean(!sig_all & sig_post, na.rm = TRUE),
      prop_sign_flip = mean(sign_all * sign_post < 0, na.rm = TRUE),
      .groups = "drop"
    )

  write.csv(sweep_summary_df, sweep_summary_path, row.names = FALSE)

  sweep_summary_2yr <- sweep_summary_df %>%
    mutate(cutoff = as.integer(cutoff)) %>%
    filter(cutoff %% 2 == 0)
  write.csv(sweep_summary_2yr, sweep_summary_2yr_path, row.names = FALSE)

  sweep_2yr <- sweep_df %>%
    mutate(cutoff = as.integer(cutoff)) %>%
    filter(cutoff %% 2 == 0)

  if (nrow(sweep_2yr) > 0) {
    sweep_2yr_by_type <- sweep_2yr %>%
      group_by(cutoff, data_type) %>%
      summarize(
        n_analyses = n(),
        prop_delta_i2_gt0 = mean(delta_i2 > 0, na.rm = TRUE),
        prop_abs_delta_i2_le_5 = mean(abs(delta_i2) <= 5, na.rm = TRUE),
        prop_abs_delta_i2_le_10 = mean(abs(delta_i2) <= 10, na.rm = TRUE),
        median_delta_i2 = median(delta_i2, na.rm = TRUE),
        mean_delta_i2 = mean(delta_i2, na.rm = TRUE),
        prop_sig_change = mean(sig_all != sig_post, na.rm = TRUE),
        prop_sign_flip = mean(sign_all * sign_post < 0, na.rm = TRUE),
        .groups = "drop"
      )
    write.csv(sweep_2yr_by_type, sweep_summary_2yr_by_type_path, row.names = FALSE)

    sweep_2yr_by_measure <- sweep_2yr %>%
      group_by(cutoff, measure) %>%
      summarize(
        n_analyses = n(),
        prop_delta_i2_gt0 = mean(delta_i2 > 0, na.rm = TRUE),
        prop_abs_delta_i2_le_5 = mean(abs(delta_i2) <= 5, na.rm = TRUE),
        prop_abs_delta_i2_le_10 = mean(abs(delta_i2) <= 10, na.rm = TRUE),
        median_delta_i2 = median(delta_i2, na.rm = TRUE),
        mean_delta_i2 = mean(delta_i2, na.rm = TRUE),
        prop_sig_change = mean(sig_all != sig_post, na.rm = TRUE),
        prop_sign_flip = mean(sign_all * sign_post < 0, na.rm = TRUE),
        .groups = "drop"
      )
    write.csv(sweep_2yr_by_measure, sweep_summary_2yr_by_measure_path, row.names = FALSE)

    sweep_2yr_cat <- sweep_2yr %>%
      mutate(outcome_category = vapply(analysis_name, categorize_analysis_name, character(1)))
    sweep_2yr_by_category <- sweep_2yr_cat %>%
      group_by(cutoff, outcome_category) %>%
      summarize(
        n_analyses = n(),
        prop_delta_i2_gt0 = mean(delta_i2 > 0, na.rm = TRUE),
        prop_abs_delta_i2_le_5 = mean(abs(delta_i2) <= 5, na.rm = TRUE),
        prop_abs_delta_i2_le_10 = mean(abs(delta_i2) <= 10, na.rm = TRUE),
        median_delta_i2 = median(delta_i2, na.rm = TRUE),
        mean_delta_i2 = mean(delta_i2, na.rm = TRUE),
        prop_sig_change = mean(sig_all != sig_post, na.rm = TRUE),
        prop_sign_flip = mean(sign_all * sign_post < 0, na.rm = TRUE),
        .groups = "drop"
      )
    write.csv(sweep_2yr_by_category, sweep_summary_2yr_by_category_path, row.names = FALSE)

    sweep_2yr_cat_lr <- sweep_2yr_cat %>%
      mutate(measure_group = mapply(measure_group, measure, effect_is_log)) %>%
      filter(measure_group == "log_ratio")
    sweep_2yr_by_category_lr <- sweep_2yr_cat_lr %>%
      group_by(cutoff, outcome_category) %>%
      summarize(
        n_analyses = n(),
        prop_delta_i2_gt0 = mean(delta_i2 > 0, na.rm = TRUE),
        prop_abs_delta_i2_le_5 = mean(abs(delta_i2) <= 5, na.rm = TRUE),
        prop_abs_delta_i2_le_10 = mean(abs(delta_i2) <= 10, na.rm = TRUE),
        median_delta_i2 = median(delta_i2, na.rm = TRUE),
        mean_delta_i2 = mean(delta_i2, na.rm = TRUE),
        prop_sig_change = mean(sig_all != sig_post, na.rm = TRUE),
        prop_sign_flip = mean(sign_all * sign_post < 0, na.rm = TRUE),
        .groups = "drop"
      )
    write.csv(sweep_2yr_by_category_lr, sweep_summary_2yr_by_category_log_ratio_path, row.names = FALSE)

    sweep_2yr_cat_md <- sweep_2yr_cat %>%
      mutate(measure_group = mapply(measure_group, measure, effect_is_log)) %>%
      filter(measure_group == "mean_diff")
    sweep_2yr_by_category_md <- sweep_2yr_cat_md %>%
      group_by(cutoff, outcome_category) %>%
      summarize(
        n_analyses = n(),
        prop_delta_i2_gt0 = mean(delta_i2 > 0, na.rm = TRUE),
        prop_abs_delta_i2_le_5 = mean(abs(delta_i2) <= 5, na.rm = TRUE),
        prop_abs_delta_i2_le_10 = mean(abs(delta_i2) <= 10, na.rm = TRUE),
        median_delta_i2 = median(delta_i2, na.rm = TRUE),
        mean_delta_i2 = mean(delta_i2, na.rm = TRUE),
        prop_sig_change = mean(sig_all != sig_post, na.rm = TRUE),
        prop_sign_flip = mean(sign_all * sign_post < 0, na.rm = TRUE),
        .groups = "drop"
      )
    write.csv(sweep_2yr_by_category_md, sweep_summary_2yr_by_category_mean_diff_path, row.names = FALSE)
  }
}

rob_results_path_use <- rob_results_path
recompute_rob <- Sys.getenv("PAIRWISE70_RECOMPUTE_ROB") == "1"
if (recompute_rob && file.exists(rob_results_path_use)) {
  removed <- tryCatch(file.remove(rob_results_path_use), warning = function(w) FALSE, error = function(e) FALSE)
  if (!isTRUE(removed)) {
    rob_results_path_use <- file.path(out_dir, "rob_low_results_recompute.csv")
  }
}

reuse_rob <- file.exists(rob_results_path_use) && !recompute_rob
if (reuse_rob) {
  rob_df <- read.csv(rob_results_path_use, stringsAsFactors = FALSE)
} else {
  rda_files <- list.files(data_dir, pattern = "\\.rda$", full.names = TRUE)
  resume_rob <- file.exists(rob_results_path_use) && Sys.getenv("PAIRWISE70_ROB_RESUME") == "1"
  if (!resume_rob && file.exists(rob_results_path_use)) {
    file.remove(rob_results_path_use)
  }
  processed_ids <- character()
  if (resume_rob) {
    processed_ids <- unique(read.csv(rob_results_path_use, stringsAsFactors = FALSE)$dataset_id)
  }
  rob_has_header <- file.exists(rob_results_path_use)
  rob_filters <- list(
    low_only = list(
      cols = c(bias_col),
      fun = function(df) is_low_risk(df[[bias_col]])
    ),
    non_high = list(
      cols = c(bias_col),
      fun = function(df) is_non_high(df[[bias_col]])
    ),
    low_only_with_support = list(
      cols = c(bias_col, bias_support_col),
      fun = function(df) is_low_risk(df[[bias_col]]) & has_support(df[[bias_support_col]])
    )
  )

  for (i in seq_along(rda_files)) {
    path <- rda_files[i]
    dataset_id <- tools::file_path_sans_ext(basename(path))
    if (dataset_id %in% processed_ids) {
      next
    }
    df <- read_dataset(path)

    if (!all(c("Analysis.group", "Analysis.number", "Analysis.name") %in% names(df))) {
      next
    }
    if (!(bias_col %in% names(df))) {
      next
    }

    review_doi <- if ("review_doi" %in% names(df)) {
      unique(na.omit(df$review_doi))[1]
    } else {
      NA_character_
    }

    analysis_keys <- unique(df[, c("Analysis.group", "Analysis.number", "Analysis.name")])

    rob_list <- list()
    for (j in seq_len(nrow(analysis_keys))) {
      key <- analysis_keys[j, ]
      df_sub <- df[
        df$Analysis.group == key$Analysis.group &
          df$Analysis.number == key$Analysis.number &
          df$Analysis.name == key$Analysis.name,
        ,
        drop = FALSE
      ]

      method <- choose_method(df_sub, key$Analysis.name)
      if (is.null(method)) {
        next
      }

      for (rob_name in names(rob_filters)) {
        rob_def <- rob_filters[[rob_name]]
        if (!all(rob_def$cols %in% names(df_sub))) {
          next
        }
        df_low <- df_sub[rob_def$fun(df_sub), , drop = FALSE]
        if (nrow(df_low) < 2) {
          next
        }

        es_info <- extract_es(df_low, method)
        es_df <- es_info$es
        if (nrow(es_df) < 2) {
          next
        }

        es_post <- es_df[es_df$year >= 2005, , drop = FALSE]
        es_pre <- es_df[es_df$year < 2005, , drop = FALSE]
        if (nrow(es_post) < 2 || nrow(es_pre) < 1) {
          next
        }

        meta_all <- fit_meta(es_df, method$effect_is_log, method = "REML")
        meta_post <- fit_meta(es_post, method$effect_is_log, method = "REML")
        meta_pre <- fit_meta(es_pre, method$effect_is_log, method = "REML")

        rob_list[[length(rob_list) + 1]] <- data.frame(
          dataset_id = dataset_id,
          review_doi = review_doi,
          analysis_group = key$Analysis.group,
          analysis_number = key$Analysis.number,
          analysis_name = key$Analysis.name,
          data_type = method$type,
          measure = method$measure,
          effect_is_log = method$effect_is_log,
          rob_type = rob_name,
          k_all = meta_all$k,
          k_post = meta_post$k,
          k_pre = nrow(es_pre),
          frac_pre = if (meta_all$k > 0) nrow(es_pre) / meta_all$k else NA_real_,
          i2_all = meta_all$i2,
          i2_post = meta_post$i2,
          i2_pre = meta_pre$i2,
          delta_i2 = if (!is.na(meta_all$i2) && !is.na(meta_post$i2)) meta_post$i2 - meta_all$i2 else NA_real_,
          delta_i2_post_pre = if (!is.na(meta_pre$i2) && !is.na(meta_post$i2)) meta_post$i2 - meta_pre$i2 else NA_real_,
          tau2_all = meta_all$tau2,
          tau2_post = meta_post$tau2,
          tau2_pre = meta_pre$tau2,
          delta_tau2 = if (!is.na(meta_all$tau2) && !is.na(meta_post$tau2)) meta_post$tau2 - meta_all$tau2 else NA_real_,
          delta_tau2_post_pre = if (!is.na(meta_pre$tau2) && !is.na(meta_post$tau2)) meta_post$tau2 - meta_pre$tau2 else NA_real_,
          effect_all = meta_all$effect,
          effect_post = meta_post$effect,
          effect_pre = meta_pre$effect,
          delta_effect = if (!is.na(meta_all$effect) && !is.na(meta_post$effect)) meta_post$effect - meta_all$effect else NA_real_,
          sig_all = sig_flag(meta_all),
          sig_post = sig_flag(meta_post),
          sign_all = sign_flag(meta_all),
          sign_post = sign_flag(meta_post),
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(rob_list) > 0) {
      rob_df_chunk <- bind_rows(rob_list)
      write.table(
        rob_df_chunk,
        rob_results_path_use,
        sep = ",",
        row.names = FALSE,
        col.names = !rob_has_header,
        append = rob_has_header
      )
      rob_has_header <- TRUE
    }

    if (i %% 25 == 0) {
      message("ROB sweep processed ", i, " / ", length(rda_files), " datasets")
    }
  }

  rob_df <- if (file.exists(rob_results_path_use)) {
    read.csv(rob_results_path_use, stringsAsFactors = FALSE)
  } else {
    data.frame()
  }
}

if (exists("rob_df") && nrow(rob_df) > 0) {
  summarize_rob <- function(df, min_k_post) {
    df_sub <- df[df$k_post >= min_k_post, , drop = FALSE]
    if (nrow(df_sub) == 0) {
      return(data.frame())
    }
    df_sub %>%
      group_by(rob_type) %>%
      summarize(
        min_k_post = min_k_post,
        n_analyses = n(),
        prop_delta_i2_gt0 = mean(delta_i2 > 0, na.rm = TRUE),
        prop_abs_delta_i2_le_5 = mean(abs(delta_i2) <= 5, na.rm = TRUE),
        prop_abs_delta_i2_le_10 = mean(abs(delta_i2) <= 10, na.rm = TRUE),
        median_delta_i2 = median(delta_i2, na.rm = TRUE),
        mean_delta_i2 = mean(delta_i2, na.rm = TRUE),
        prop_sig_change = mean(sig_all != sig_post, na.rm = TRUE),
        prop_sig_loss = mean(sig_all & !sig_post, na.rm = TRUE),
        prop_sig_gain = mean(!sig_all & sig_post, na.rm = TRUE),
        prop_sign_flip = mean(sign_all * sign_post < 0, na.rm = TRUE),
        .groups = "drop"
      )
  }

  rob_summary_df <- bind_rows(
    summarize_rob(rob_df, 2),
    summarize_rob(rob_df, 5),
    summarize_rob(rob_df, 10)
  )

  write.csv(rob_summary_df, rob_summary_csv, row.names = FALSE)

  sink(rob_summary_path)
  cat("Randomization bias filter summary (pre-2005 exclusion)\n\n")
  print(rob_summary_df)
  sink()

  rob_by_type <- rob_df %>%
    group_by(rob_type, data_type) %>%
    summarize(
      n_analyses = n(),
      prop_delta_i2_gt0 = mean(delta_i2 > 0, na.rm = TRUE),
      prop_abs_delta_i2_le_5 = mean(abs(delta_i2) <= 5, na.rm = TRUE),
      prop_abs_delta_i2_le_10 = mean(abs(delta_i2) <= 10, na.rm = TRUE),
      median_delta_i2 = median(delta_i2, na.rm = TRUE),
      mean_delta_i2 = mean(delta_i2, na.rm = TRUE),
      prop_sig_change = mean(sig_all != sig_post, na.rm = TRUE),
      prop_sign_flip = mean(sign_all * sign_post < 0, na.rm = TRUE),
      .groups = "drop"
    )
  write.csv(rob_by_type, rob_summary_by_type_path, row.names = FALSE)

  rob_by_measure <- rob_df %>%
    group_by(rob_type, measure) %>%
    summarize(
      n_analyses = n(),
      prop_delta_i2_gt0 = mean(delta_i2 > 0, na.rm = TRUE),
      prop_abs_delta_i2_le_5 = mean(abs(delta_i2) <= 5, na.rm = TRUE),
      prop_abs_delta_i2_le_10 = mean(abs(delta_i2) <= 10, na.rm = TRUE),
      median_delta_i2 = median(delta_i2, na.rm = TRUE),
      mean_delta_i2 = mean(delta_i2, na.rm = TRUE),
      prop_sig_change = mean(sig_all != sig_post, na.rm = TRUE),
      prop_sign_flip = mean(sign_all * sign_post < 0, na.rm = TRUE),
      .groups = "drop"
    )
  write.csv(rob_by_measure, rob_summary_by_measure_path, row.names = FALSE)

  rob_cat <- rob_df %>%
    mutate(outcome_category = vapply(analysis_name, categorize_analysis_name, character(1)))
  rob_by_category <- rob_cat %>%
    group_by(rob_type, outcome_category) %>%
    summarize(
      n_analyses = n(),
      prop_delta_i2_gt0 = mean(delta_i2 > 0, na.rm = TRUE),
      prop_abs_delta_i2_le_5 = mean(abs(delta_i2) <= 5, na.rm = TRUE),
      prop_abs_delta_i2_le_10 = mean(abs(delta_i2) <= 10, na.rm = TRUE),
      median_delta_i2 = median(delta_i2, na.rm = TRUE),
      mean_delta_i2 = mean(delta_i2, na.rm = TRUE),
      prop_sig_change = mean(sig_all != sig_post, na.rm = TRUE),
      prop_sign_flip = mean(sign_all * sign_post < 0, na.rm = TRUE),
      .groups = "drop"
    )
  write.csv(rob_by_category, rob_summary_by_category_path, row.names = FALSE)

  rob_cat_lr <- rob_cat %>%
    mutate(measure_group = mapply(measure_group, measure, effect_is_log)) %>%
    filter(measure_group == "log_ratio")
  rob_by_category_lr <- rob_cat_lr %>%
    group_by(rob_type, outcome_category) %>%
    summarize(
      n_analyses = n(),
      prop_delta_i2_gt0 = mean(delta_i2 > 0, na.rm = TRUE),
      prop_abs_delta_i2_le_5 = mean(abs(delta_i2) <= 5, na.rm = TRUE),
      prop_abs_delta_i2_le_10 = mean(abs(delta_i2) <= 10, na.rm = TRUE),
      median_delta_i2 = median(delta_i2, na.rm = TRUE),
      mean_delta_i2 = mean(delta_i2, na.rm = TRUE),
      prop_sig_change = mean(sig_all != sig_post, na.rm = TRUE),
      prop_sign_flip = mean(sign_all * sign_post < 0, na.rm = TRUE),
      .groups = "drop"
    )
  write.csv(rob_by_category_lr, rob_summary_by_category_log_ratio_path, row.names = FALSE)

  rob_cat_md <- rob_cat %>%
    mutate(measure_group = mapply(measure_group, measure, effect_is_log)) %>%
    filter(measure_group == "mean_diff")
  rob_by_category_md <- rob_cat_md %>%
    group_by(rob_type, outcome_category) %>%
    summarize(
      n_analyses = n(),
      prop_delta_i2_gt0 = mean(delta_i2 > 0, na.rm = TRUE),
      prop_abs_delta_i2_le_5 = mean(abs(delta_i2) <= 5, na.rm = TRUE),
      prop_abs_delta_i2_le_10 = mean(abs(delta_i2) <= 10, na.rm = TRUE),
      median_delta_i2 = median(delta_i2, na.rm = TRUE),
      mean_delta_i2 = mean(delta_i2, na.rm = TRUE),
      prop_sig_change = mean(sig_all != sig_post, na.rm = TRUE),
      prop_sign_flip = mean(sign_all * sign_post < 0, na.rm = TRUE),
      .groups = "drop"
    )
  write.csv(rob_by_category_md, rob_summary_by_category_mean_diff_path, row.names = FALSE)
}

advanced_results_path_use <- advanced_results_path
recompute_advanced <- Sys.getenv("PAIRWISE70_RECOMPUTE_ADVANCED") == "1"
advanced_results_alt <- file.path(out_dir, "advanced_model_results_recompute.csv")
if (!recompute_advanced && file.exists(advanced_results_alt)) {
  use_alt <- !file.exists(advanced_results_path_use)
  if (!use_alt) {
    use_alt <- file.info(advanced_results_alt)$size > file.info(advanced_results_path_use)$size
  }
  if (use_alt) {
    advanced_results_path_use <- advanced_results_alt
  }
}
if (recompute_advanced && file.exists(advanced_results_path_use)) {
  removed <- tryCatch(file.remove(advanced_results_path_use), warning = function(w) FALSE, error = function(e) FALSE)
  if (!isTRUE(removed)) {
    advanced_results_path_use <- file.path(out_dir, "advanced_model_results_recompute.csv")
  }
}

reuse_advanced <- file.exists(advanced_results_path_use) && !recompute_advanced
if (reuse_advanced) {
  advanced_df <- read.csv(advanced_results_path_use, stringsAsFactors = FALSE)
} else {
  rda_files <- list.files(data_dir, pattern = "\\.rda$", full.names = TRUE)
  advanced_has_header <- file.exists(advanced_results_path_use)

  for (i in seq_along(rda_files)) {
    path <- rda_files[i]
    dataset_id <- tools::file_path_sans_ext(basename(path))
    df <- read_dataset(path)

    if (!all(c("Analysis.group", "Analysis.number", "Analysis.name") %in% names(df))) {
      next
    }

    review_doi <- if ("review_doi" %in% names(df)) {
      unique(na.omit(df$review_doi))[1]
    } else {
      NA_character_
    }

    analysis_keys <- unique(df[, c("Analysis.group", "Analysis.number", "Analysis.name")])
    chunk <- list()

    for (j in seq_len(nrow(analysis_keys))) {
      key <- analysis_keys[j, ]
      df_sub <- df[
        df$Analysis.group == key$Analysis.group &
          df$Analysis.number == key$Analysis.number &
          df$Analysis.name == key$Analysis.name,
        ,
        drop = FALSE
      ]

      method <- choose_method(df_sub, key$Analysis.name)
      if (is.null(method)) {
        next
      }

      es_info <- extract_es(df_sub, method)
      es_df <- es_info$es
      if (nrow(es_df) < 2) {
        next
      }

      years <- es_df$year
      k <- nrow(es_df)
      year_unique <- length(unique(years))
      year_min <- min(years, na.rm = TRUE)
      year_max <- max(years, na.rm = TRUE)
      pre <- as.integer(years < 2005)
      k_pre <- sum(pre == 1, na.rm = TRUE)
      k_post <- sum(pre == 0, na.rm = TRUE)

      res_all <- tryCatch(rma(yi = es_df$yi, vi = es_df$vi, method = "REML"), error = function(e) NULL)

      slope_year <- slope_year_se <- slope_year_p <- slope_year_ci_lb <- slope_year_ci_ub <- NA_real_
      aic_linear <- aic_spline <- delta_aic <- NA_real_

      if (k >= 5 && year_unique >= 3) {
        year_c <- years - mean(years, na.rm = TRUE)
        res_year <- tryCatch(rma(yi = es_df$yi, vi = es_df$vi, mods = ~ year_c, method = "REML"), error = function(e) NULL)
        if (!is.null(res_year)) {
          slope_year <- as.numeric(res_year$beta[2])
          slope_year_se <- as.numeric(res_year$se[2])
          slope_year_p <- as.numeric(res_year$pval[2])
          slope_year_ci_lb <- as.numeric(res_year$ci.lb[2])
          slope_year_ci_ub <- as.numeric(res_year$ci.ub[2])
        }

        if (k >= 10 && year_unique >= 5 && !is.null(res_year)) {
          res_spline <- tryCatch(
            rma(yi = es_df$yi, vi = es_df$vi, mods = ~ ns(year_c, df = 3), method = "REML"),
            error = function(e) NULL
          )
          if (!is.null(res_spline)) {
            aic_linear <- AIC(res_year)
            aic_spline <- AIC(res_spline)
            delta_aic <- aic_spline - aic_linear
          }
        }
      }

      slope_pre <- slope_pre_se <- slope_pre_p <- slope_pre_ci_lb <- slope_pre_ci_ub <- NA_real_
      if (k_pre >= 2 && k_post >= 2) {
        res_pre <- tryCatch(rma(yi = es_df$yi, vi = es_df$vi, mods = ~ pre, method = "REML"), error = function(e) NULL)
        if (!is.null(res_pre)) {
          slope_pre <- as.numeric(res_pre$beta[2])
          slope_pre_se <- as.numeric(res_pre$se[2])
          slope_pre_p <- as.numeric(res_pre$pval[2])
          slope_pre_ci_lb <- as.numeric(res_pre$ci.lb[2])
          slope_pre_ci_ub <- as.numeric(res_pre$ci.ub[2])
        }
      }

      egger_p_all <- egger_p_post <- NA_real_
      if (!is.null(res_all) && k >= 10) {
        egger_all <- tryCatch(regtest(res_all), error = function(e) NULL)
        if (!is.null(egger_all)) {
          egger_p_all <- as.numeric(egger_all$pval)
        }
      }
      if (k_post >= 10) {
        es_post <- es_df[pre == 0, , drop = FALSE]
        res_post <- tryCatch(rma(yi = es_post$yi, vi = es_post$vi, method = "REML"), error = function(e) NULL)
        if (!is.null(res_post)) {
          egger_post <- tryCatch(regtest(res_post), error = function(e) NULL)
          if (!is.null(egger_post)) {
            egger_p_post <- as.numeric(egger_post$pval)
          }
        }
      }

      chunk[[length(chunk) + 1]] <- data.frame(
        dataset_id = dataset_id,
        review_doi = review_doi,
        analysis_group = key$Analysis.group,
        analysis_number = key$Analysis.number,
        analysis_name = key$Analysis.name,
        data_type = method$type,
        measure = method$measure,
        effect_is_log = method$effect_is_log,
        k = k,
        k_pre = k_pre,
        k_post = k_post,
        year_min = year_min,
        year_max = year_max,
        year_unique = year_unique,
        slope_year = slope_year,
        slope_year_se = slope_year_se,
        slope_year_p = slope_year_p,
        slope_year_ci_lb = slope_year_ci_lb,
        slope_year_ci_ub = slope_year_ci_ub,
        slope_pre = slope_pre,
        slope_pre_se = slope_pre_se,
        slope_pre_p = slope_pre_p,
        slope_pre_ci_lb = slope_pre_ci_lb,
        slope_pre_ci_ub = slope_pre_ci_ub,
        aic_linear = aic_linear,
        aic_spline = aic_spline,
        delta_aic = delta_aic,
        egger_p_all = egger_p_all,
        egger_p_post = egger_p_post,
        stringsAsFactors = FALSE
      )
    }

    if (length(chunk) > 0) {
      chunk_df <- bind_rows(chunk)
      write.table(
        chunk_df,
        advanced_results_path_use,
        sep = ",",
        row.names = FALSE,
        col.names = !advanced_has_header,
        append = advanced_has_header
      )
      advanced_has_header <- TRUE
    }

    if (i %% 25 == 0) {
      message("Advanced models processed ", i, " / ", length(rda_files), " datasets")
    }
  }

  advanced_df <- if (file.exists(advanced_results_path_use)) {
    read.csv(advanced_results_path_use, stringsAsFactors = FALSE)
  } else {
    data.frame()
  }
}

if (exists("advanced_df") && nrow(advanced_df) > 0) {
  advanced_df <- advanced_df %>%
    mutate(
      outcome_category = vapply(analysis_name, categorize_analysis_name, character(1)),
      measure_group = mapply(measure_group, measure, effect_is_log)
    )

  meta_summarize <- function(df, slope_col, se_col, group_cols, model_label, transform) {
    df_use <- df[!is.na(df[[slope_col]]) & !is.na(df[[se_col]]) & df[[se_col]] > 0, , drop = FALSE]
    if (nrow(df_use) == 0) {
      return(data.frame())
    }
    groups <- df_use %>%
      group_by(across(all_of(group_cols))) %>%
      group_split()

    out <- list()
    for (g in groups) {
      group_vals <- g[1, group_cols, drop = FALSE]
      yi <- g[[slope_col]]
      vi <- g[[se_col]]^2
      if (length(yi) < 5) {
        next
      }
      res <- tryCatch(rma(yi = yi, vi = vi, method = "REML"), error = function(e) NULL)
      if (is.null(res)) {
        next
      }
      est <- as.numeric(res$b)
      ci.lb <- as.numeric(res$ci.lb)
      ci.ub <- as.numeric(res$ci.ub)
      extra <- data.frame(
        model = model_label,
        k = length(yi),
        estimate = est,
        se = as.numeric(res$se),
        ci.lb = ci.lb,
        ci.ub = ci.ub,
        p = as.numeric(res$pval),
        tau2 = as.numeric(res$tau2),
        i2 = as.numeric(res$I2),
        stringsAsFactors = FALSE
      )
      if (transform == "log_ratio_10y") {
        extra$ratio_10y <- exp(10 * est)
        extra$ratio_10y_ci_lb <- exp(10 * ci.lb)
        extra$ratio_10y_ci_ub <- exp(10 * ci.ub)
      } else if (transform == "log_ratio") {
        extra$ratio <- exp(est)
        extra$ratio_ci_lb <- exp(ci.lb)
        extra$ratio_ci_ub <- exp(ci.ub)
      }
      out[[length(out) + 1]] <- cbind(group_vals, extra)
    }
    bind_rows(out)
  }

  meta_year <- bind_rows(
    meta_summarize(advanced_df, "slope_year", "slope_year_se", c("measure_group"), "year_slope_by_measure_group", "log_ratio_10y"),
    meta_summarize(advanced_df, "slope_year", "slope_year_se", c("data_type"), "year_slope_by_type", "none"),
    meta_summarize(advanced_df, "slope_year", "slope_year_se", c("outcome_category"), "year_slope_by_category", "none")
  )

  meta_pre <- bind_rows(
    meta_summarize(advanced_df, "slope_pre", "slope_pre_se", c("measure_group"), "pre_indicator_by_measure_group", "log_ratio"),
    meta_summarize(advanced_df, "slope_pre", "slope_pre_se", c("data_type"), "pre_indicator_by_type", "none"),
    meta_summarize(advanced_df, "slope_pre", "slope_pre_se", c("outcome_category"), "pre_indicator_by_category", "none")
  )

  meta_summary_df <- bind_rows(meta_year, meta_pre)
  write.csv(meta_summary_df, advanced_meta_summary_path, row.names = FALSE)

  advanced_lr <- advanced_df %>%
    filter(measure_group == "log_ratio")
  meta_year_lr <- meta_summarize(advanced_lr, "slope_year", "slope_year_se",
                                 c("outcome_category"), "year_slope_by_category_log_ratio", "log_ratio_10y")
  meta_pre_lr <- meta_summarize(advanced_lr, "slope_pre", "slope_pre_se",
                                c("outcome_category"), "pre_indicator_by_category_log_ratio", "log_ratio")
  meta_summary_lr <- bind_rows(meta_year_lr, meta_pre_lr)
  write.csv(meta_summary_lr, advanced_meta_summary_log_ratio_path, row.names = FALSE)

  advanced_md <- advanced_df %>%
    filter(measure_group == "mean_diff")
  meta_year_md <- meta_summarize(advanced_md, "slope_year", "slope_year_se",
                                 c("outcome_category"), "year_slope_by_category_mean_diff", "none")
  meta_pre_md <- meta_summarize(advanced_md, "slope_pre", "slope_pre_se",
                                c("outcome_category"), "pre_indicator_by_category_mean_diff", "none")
  meta_summary_md <- bind_rows(meta_year_md, meta_pre_md)
  write.csv(meta_summary_md, advanced_meta_summary_mean_diff_path, row.names = FALSE)

  rve_recompute <- Sys.getenv("PAIRWISE70_RECOMPUTE_RVE") == "1"
  if (rve_recompute || !file.exists(rve_summary_path)) {
    rve_extract <- function(reg_row, name) {
      if (!is.null(reg_row) && name %in% names(reg_row)) {
        return(as.numeric(reg_row[[name]]))
      }
      NA_real_
    }

    run_rve <- function(df, label, slope_col, se_col) {
      df_use <- df[!is.na(df[[slope_col]]) & !is.na(df[[se_col]]) & df[[se_col]] > 0, , drop = FALSE]
      if (nrow(df_use) < 5 || length(unique(df_use$dataset_id)) < 2) {
        return(data.frame())
      }

      fit <- tryCatch(
        robu(df_use[[slope_col]] ~ 1,
             data = df_use,
             studynum = df_use$dataset_id,
             var.eff = df_use[[se_col]]^2,
             modelweights = "HIER",
             small = TRUE),
        error = function(e) NULL
      )
      if (is.null(fit) || is.null(fit$reg_table)) {
        return(data.frame())
      }

      reg <- as.data.frame(fit$reg_table)
      reg_row <- reg[1, , drop = FALSE]
      estimate <- rve_extract(reg_row, "beta")
      se <- rve_extract(reg_row, "SE")
      tval <- rve_extract(reg_row, "t")
      df_rve <- rve_extract(reg_row, "df")
      pval <- rve_extract(reg_row, "p")
      ci.lb <- rve_extract(reg_row, "CI.L")
      ci.ub <- rve_extract(reg_row, "CI.U")

      data.frame(
        measure_group = label,
        slope = slope_col,
        k = nrow(df_use),
        clusters = length(unique(df_use$dataset_id)),
        estimate = estimate,
        se = se,
        t = tval,
        df = df_rve,
        p = pval,
        ci.lb = ci.lb,
        ci.ub = ci.ub,
        ratio = if (label == "log_ratio" && slope_col == "slope_pre") exp(estimate) else NA_real_,
        ratio_ci_lb = if (label == "log_ratio" && slope_col == "slope_pre") exp(ci.lb) else NA_real_,
        ratio_ci_ub = if (label == "log_ratio" && slope_col == "slope_pre") exp(ci.ub) else NA_real_,
        ratio_10y = if (label == "log_ratio" && slope_col == "slope_year") exp(10 * estimate) else NA_real_,
        ratio_10y_ci_lb = if (label == "log_ratio" && slope_col == "slope_year") exp(10 * ci.lb) else NA_real_,
        ratio_10y_ci_ub = if (label == "log_ratio" && slope_col == "slope_year") exp(10 * ci.ub) else NA_real_,
        stringsAsFactors = FALSE
      )
    }

    rve_targets <- list(
      list(label = "log_ratio", df = advanced_df %>% filter(measure_group == "log_ratio")),
      list(label = "mean_diff", df = advanced_df %>% filter(measure_group == "mean_diff"))
    )
    rve_specs <- list(
      list(slope_col = "slope_pre", se_col = "slope_pre_se"),
      list(slope_col = "slope_year", se_col = "slope_year_se")
    )

    rve_out <- list()
    for (target in rve_targets) {
      for (spec in rve_specs) {
        rve_out[[length(rve_out) + 1]] <- run_rve(
          target$df,
          target$label,
          spec$slope_col,
          spec$se_col
        )
      }
    }
    rve_df <- bind_rows(rve_out)
    if (nrow(rve_df) > 0) {
      write.csv(rve_df, rve_summary_path, row.names = FALSE)
    }
  }

  robma_recompute <- Sys.getenv("PAIRWISE70_RECOMPUTE_ROBMA") == "1"
  if (robma_recompute || !file.exists(robma_summary_path)) {
    set.seed(1)
    robma_targets <- list(
      list(label = "log_ratio_pre", df = advanced_df %>% filter(measure_group == "log_ratio"),
           slope_col = "slope_pre", se_col = "slope_pre_se"),
      list(label = "log_ratio_year", df = advanced_df %>% filter(measure_group == "log_ratio"),
           slope_col = "slope_year", se_col = "slope_year_se"),
      list(label = "mean_diff_pre", df = advanced_df %>% filter(measure_group == "mean_diff"),
           slope_col = "slope_pre", se_col = "slope_pre_se"),
      list(label = "mean_diff_year", df = advanced_df %>% filter(measure_group == "mean_diff"),
           slope_col = "slope_year", se_col = "slope_year_se")
    )

    robma_inference <- list()
    robma_models <- list()

    sink(robma_summary_path)
    cat("RoBMA summary for pre-2005 analyses\n\n")
    for (target in robma_targets) {
      df_use <- target$df
      df_use <- df_use[!is.na(df_use[[target$slope_col]]) &
                         !is.na(df_use[[target$se_col]]) &
                         df_use[[target$se_col]] > 0, , drop = FALSE]
      cat("\n== ", target$label, " ==\n", sep = "")
      cat("Analyses used: ", nrow(df_use), "\n", sep = "")
      if (nrow(df_use) < 10) {
        cat("Skipped (insufficient data).\n")
        next
      }

      fit <- tryCatch(
        RoBMA(y = df_use[[target$slope_col]], se = df_use[[target$se_col]]),
        error = function(e) e
      )

      if (inherits(fit, "error")) {
        cat("RoBMA failed: ", fit$message, "\n", sep = "")
        next
      }

      print(fit)
      fit_summary <- tryCatch(summary(fit), error = function(e) NULL)
      if (!is.null(fit_summary)) {
        print(fit_summary)
        if (is.list(fit_summary) && "Inference" %in% names(fit_summary)) {
          inf <- tryCatch(as.data.frame(fit_summary$Inference), error = function(e) NULL)
          if (!is.null(inf)) {
            inf$label <- target$label
            robma_inference[[length(robma_inference) + 1]] <- inf
          }
        }
        if (is.list(fit_summary) && "Models" %in% names(fit_summary)) {
          mods <- tryCatch(as.data.frame(fit_summary$Models), error = function(e) NULL)
          if (!is.null(mods)) {
            mods$label <- target$label
            robma_models[[length(robma_models) + 1]] <- mods
          }
        }
      }
    }
    sink()

    if (length(robma_inference) > 0) {
      write.csv(bind_rows(robma_inference), robma_inference_path, row.names = FALSE)
    }
    if (length(robma_models) > 0) {
      write.csv(bind_rows(robma_models), robma_models_path, row.names = FALSE)
    }
  }

  summary_path <- advanced_summary_path
  sink(summary_path)
  cat("Advanced model summary\n\n")
  cat("Analyses with year slope:", sum(!is.na(advanced_df$slope_year)), "\n")
  cat("Analyses with pre2005 slope:", sum(!is.na(advanced_df$slope_pre)), "\n")
  cat("Analyses with spline comparison:", sum(!is.na(advanced_df$delta_aic)), "\n")
  cat("Spline better (delta AIC <= -2):", mean(advanced_df$delta_aic <= -2, na.rm = TRUE), "\n")
  cat("Spline strongly better (delta AIC <= -10):", mean(advanced_df$delta_aic <= -10, na.rm = TRUE), "\n")
  cat("Egger p<0.10 all:", mean(advanced_df$egger_p_all < 0.10, na.rm = TRUE), "\n")
  cat("Egger p<0.10 post-2005:", mean(advanced_df$egger_p_post < 0.10, na.rm = TRUE), "\n")
  sink()
}

if (exists("sens_df") && nrow(sens_df) > 0) {
  analysis_df <- sens_df %>%
    filter(k_post >= 2, k_pre > 0, !is.na(delta_i2_re)) %>%
    mutate(
      outcome_category = vapply(analysis_name, categorize_analysis_name, character(1)),
      measure_group = mapply(measure_group, measure, effect_is_log),
      log_k_all = log1p(k_all)
    )

  category_counts <- table(analysis_df$outcome_category)
  small_cats <- names(category_counts[category_counts < 30])
  analysis_df$outcome_category <- ifelse(analysis_df$outcome_category %in% small_cats, "other_small", analysis_df$outcome_category)

  fit_lmer <- tryCatch(
    lmer(delta_i2_re ~ frac_pre + log_k_all + data_type + measure_group + outcome_category + (1 | dataset_id),
         data = analysis_df, REML = FALSE),
    error = function(e) NULL
  )

  if (!is.null(fit_lmer)) {
    vc <- tryCatch(vcovCR(fit_lmer, cluster = analysis_df$dataset_id, type = "CR2"), error = function(e) NULL)
    if (!is.null(vc)) {
      coef_tab <- coef_test(fit_lmer, vcov = vc, test = "Satterthwaite")
      write.csv(as.data.frame(coef_tab), advanced_mixed_path, row.names = FALSE)
    }
  }

  analysis_df$delta_i2_gt0 <- analysis_df$delta_i2_re > 0
  fit_glmer <- tryCatch(
    glmer(delta_i2_gt0 ~ frac_pre + log_k_all + data_type + measure_group + outcome_category + (1 | dataset_id),
          data = analysis_df, family = binomial()),
    error = function(e) NULL
  )

  if (!is.null(fit_glmer)) {
    vc_g <- tryCatch(vcovCR(fit_glmer, cluster = analysis_df$dataset_id, type = "CR2"), error = function(e) NULL)
    if (!is.null(vc_g)) {
      coef_tab_g <- coef_test(fit_glmer, vcov = vc_g, test = "Satterthwaite")
      write.csv(as.data.frame(coef_tab_g), advanced_mixed_logit_path, row.names = FALSE)
    }
  } else {
    fit_glm <- tryCatch(
      glm(delta_i2_gt0 ~ frac_pre + log_k_all + data_type + measure_group + outcome_category,
          data = analysis_df, family = binomial()),
      error = function(e) NULL
    )
    if (!is.null(fit_glm)) {
      vc_g <- tryCatch(vcovCR(fit_glm, cluster = analysis_df$dataset_id, type = "CR2"), error = function(e) NULL)
      if (!is.null(vc_g)) {
        coef_tab_g <- coef_test(fit_glm, vcov = vc_g, test = "Satterthwaite")
        write.csv(as.data.frame(coef_tab_g), advanced_mixed_logit_path, row.names = FALSE)
      }
    }
  }
}
