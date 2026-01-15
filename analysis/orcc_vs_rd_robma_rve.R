#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop("Package 'metafor' is required. Install with install.packages('metafor').")
  }
  if (!requireNamespace("clubSandwich", quietly = TRUE)) {
    stop("Package 'clubSandwich' is required. Install with install.packages('clubSandwich').")
  }
})

robma_available <- requireNamespace("RoBMA", quietly = TRUE)

standardize_names <- function(nms) {
  nms <- gsub("[^A-Za-z0-9]+", ".", nms)
  nms <- gsub("\\.+", ".", nms)
  nms <- gsub("^\\.|\\.$", "", nms)
  nms
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
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
    invalid = invalid
  )
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
    ci_ub_tr = ci_ub_tr,
    model = ma
  )
}

run_rve <- function(model, cluster) {
  if (is.null(model)) {
    return(list(status = "no_model"))
  }
  if (length(unique(cluster)) < 2) {
    return(list(status = "insufficient_clusters"))
  }
  out <- tryCatch({
    ct <- clubSandwich::coef_test(model, vcov = "CR2", cluster = cluster)
    est <- ct$beta[1]
    se <- ct$SE[1]
    df <- ct$df[1]
    tcrit <- qt(0.975, df)
    ci_lb <- est - tcrit * se
    ci_ub <- est + tcrit * se
    list(status = "ok", est = est, se = se, df = df, ci_lb = ci_lb, ci_ub = ci_ub)
  }, error = function(e) {
    list(status = "error", error = e$message)
  })
  out
}

run_robma <- function(yi, sei, seed = 1, sample = 1000, burnin = 500, chains = 1) {
  if (!robma_available) {
    return(list(status = "not_installed"))
  }
  out <- tryCatch({
    fit <- RoBMA::RoBMA(
      y = yi,
      se = sei,
      transformation = "none",
      algorithm = "ss",
      sample = sample,
      burnin = burnin,
      chains = chains,
      silent = TRUE,
      seed = seed
    )
    post <- suppressWarnings(as.numeric(RoBMA::extract_posterior(fit, parameter = "mu")))
    post <- post[is.finite(post)]
    if (length(post) < 10) {
      return(list(status = "insufficient_posterior"))
    }
    list(
      status = "ok",
      mu = mean(post),
      ci_lb = quantile(post, 0.025),
      ci_ub = quantile(post, 0.975)
    )
  }, error = function(e) {
    list(status = "error", error = e$message)
  })
  out
}

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else file.path("data")
output_dir <- if (length(args) >= 2) args[2] else file.path("analysis", "output")
max_k_robma <- if (length(args) >= 3) suppressWarnings(as.numeric(args[3])) else 0
if (is.na(max_k_robma)) max_k_robma <- 0

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

rda_files <- list.files(data_dir, pattern = "\\.rda$", full.names = TRUE)
if (length(rda_files) == 0) {
  stop(paste0("No .rda files found in: ", data_dir))
}

results_list <- list()

cat("ROBMA and RVE comparison\n")
cat("Data dir: ", data_dir, "\n", sep = "")
cat("Output dir: ", output_dir, "\n", sep = "")
cat("RoBMA available: ", robma_available, "\n", sep = "")
cat("Max k for RoBMA (0 = no limit): ", max_k_robma, "\n\n", sep = "")

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
  study <- if ("Study" %in% names(df)) as.character(df$Study) else rep(NA, n_rows)

  for (key in unique(analysis_key)) {
    idx <- analysis_key == key
    if (!any(idx)) next

    clean <- clean_binary(exp_cases[idx], exp_n[idx], ctrl_cases[idx], ctrl_n[idx])
    valid <- !clean$invalid
    k_valid <- sum(valid, na.rm = TRUE)
    if (k_valid < 2) {
      next
    }

    ai <- clean$ai[valid]
    ci <- clean$ci[valid]
    n1 <- clean$n1[valid]
    n2 <- clean$n2[valid]
    bi <- n1 - ai
    di <- n2 - ci
    cluster <- study[idx][valid]
    cluster[is.na(cluster) | cluster == ""] <- as.character(seq_along(cluster))[is.na(cluster) | cluster == ""]

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
      cluster <- cluster[drop00]
    }

    es_or <- tryCatch(
      metafor::escalc(measure = "OR", ai = or_ai, bi = or_bi, ci = or_ci, di = or_di, add = 0.5, to = "only0"),
      error = function(e) NULL
    )
    or_fit <- list(status = "error")
    if (!is.null(es_or)) {
      or_fit <- run_rma(es_or$yi, es_or$vi, TRUE)
    }

    rve_fit <- list(status = "no_model")
    if (or_fit$status == "ok") {
      rve_fit <- run_rve(or_fit$model, cluster)
    }

    robma_fit <- list(status = "not_run")
    if (robma_available && or_fit$status == "ok") {
      if (max_k_robma <= 0 || k_valid <= max_k_robma) {
        sei <- sqrt(es_or$vi)
        robma_fit <- run_robma(es_or$yi, sei, seed = 1, sample = 1000, burnin = 500, chains = 1)
      } else {
        robma_fit <- list(status = "skipped_max_k")
      }
    } else if (!robma_available) {
      robma_fit <- list(status = "not_installed")
    }

    es_rd <- tryCatch(
      metafor::escalc(measure = "RD", ai = ai, bi = bi, ci = ci, di = di),
      error = function(e) NULL
    )
    rd_fit <- list(status = "error")
    if (!is.null(es_rd)) {
      rd_fit <- run_rma(es_rd$yi, es_rd$vi, FALSE)
    }

    results_list[[length(results_list) + 1]] <- c(
      base_meta,
      list(
        orcc_status = or_fit$status,
        orcc_pooled = or_fit$pooled,
        orcc_ci_lb = or_fit$ci_lb,
        orcc_ci_ub = or_fit$ci_ub,
        rd_status = rd_fit$status,
        rd_pooled = rd_fit$pooled,
        rd_ci_lb = rd_fit$ci_lb,
        rd_ci_ub = rd_fit$ci_ub,
        rve_status = rve_fit$status,
        rve_pooled = rve_fit$est,
        rve_se = rve_fit$se,
        rve_df = rve_fit$df,
        rve_ci_lb = rve_fit$ci_lb,
        rve_ci_ub = rve_fit$ci_ub,
        robma_status = robma_fit$status,
        robma_mu = robma_fit$mu,
        robma_ci_lb = robma_fit$ci_lb,
        robma_ci_ub = robma_fit$ci_ub
      )
    )
  }

  setTxtProgressBar(pb, i)
}

close(pb)

results_df <- do.call(rbind, lapply(results_list, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
out_path <- file.path(output_dir, "orcc_vs_rd_robma_rve_results.csv")
write.csv(results_df, out_path, row.names = FALSE)

robma_ok <- results_df[results_df$robma_status == "ok", ]
rve_ok <- results_df[results_df$rve_status == "ok", ]

robma_cor <- if (nrow(robma_ok) > 5) {
  suppressWarnings(cor(as.numeric(robma_ok$orcc_pooled), as.numeric(robma_ok$robma_mu), use = "complete.obs"))
} else {
  NA
}

rve_cor <- if (nrow(rve_ok) > 5) {
  suppressWarnings(cor(as.numeric(rve_ok$orcc_pooled), as.numeric(rve_ok$rve_pooled), use = "complete.obs"))
} else {
  NA
}

or_sig <- (results_df$orcc_ci_lb > 0) | (results_df$orcc_ci_ub < 0)
rve_sig <- (results_df$rve_ci_lb > 0) | (results_df$rve_ci_ub < 0)
robma_sig <- (results_df$robma_ci_lb > 0) | (results_df$robma_ci_ub < 0)

sig_compare <- data.frame(
  or_vs_rve_both_sig = sum(or_sig & rve_sig, na.rm = TRUE),
  or_sig_only = sum(or_sig & !rve_sig, na.rm = TRUE),
  rve_sig_only = sum(!or_sig & rve_sig, na.rm = TRUE),
  or_vs_robma_both_sig = sum(or_sig & robma_sig, na.rm = TRUE),
  robma_sig_only = sum(!or_sig & robma_sig, na.rm = TRUE),
  or_sig_only_robma = sum(or_sig & !robma_sig, na.rm = TRUE),
  stringsAsFactors = FALSE
)

report_path <- file.path(output_dir, "orcc_vs_rd_robma_rve_report.md")
report_lines <- c(
  "# OR_cc vs ROBMA vs RVE comparison",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Coverage",
  paste0("- Total analyses: ", nrow(results_df)),
  paste0("- RVE ok: ", sum(results_df$rve_status == "ok")),
  paste0("- RoBMA ok: ", sum(results_df$robma_status == "ok")),
  paste0("- RoBMA skipped (max_k): ", sum(results_df$robma_status == "skipped_max_k")),
  "",
  "## Correlations",
  paste0("- OR_cc vs RVE (log OR): ", round(rve_cor, 3)),
  paste0("- OR_cc vs RoBMA (log OR): ", round(robma_cor, 3)),
  "",
  "## Significance agreement (log OR scale)",
  paste0("- OR vs RVE both significant: ", sig_compare$or_vs_rve_both_sig),
  paste0("- OR significant only: ", sig_compare$or_sig_only),
  paste0("- RVE significant only: ", sig_compare$rve_sig_only),
  paste0("- OR vs RoBMA both significant: ", sig_compare$or_vs_robma_both_sig),
  paste0("- RoBMA significant only: ", sig_compare$robma_sig_only),
  paste0("- OR significant only (RoBMA): ", sig_compare$or_sig_only_robma),
  "",
  "## Output files",
  "- orcc_vs_rd_robma_rve_results.csv",
  "- orcc_vs_rd_robma_rve_report.md"
)
writeLines(report_lines, report_path)

cat("ROBMA/RVE comparison written to: ", report_path, "\n", sep = "")
