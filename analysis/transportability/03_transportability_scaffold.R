# Transportability modeling scaffold for Pairwise70
suppressWarnings({
  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop("Package 'metafor' is required. Install with install.packages('metafor').")
  }
})

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

as_num <- function(x) suppressWarnings(as.numeric(x))

align_cols <- function(d, cols) {
  missing <- setdiff(cols, names(d))
  if (length(missing) > 0) d[missing] <- NA
  d[, cols, drop = FALSE]
}

make_analysis_key <- function(df) {
  group <- if ("Analysis.group" %in% names(df)) df[["Analysis.group"]] else NA
  number <- if ("Analysis.number" %in% names(df)) df[["Analysis.number"]] else NA
  name <- if ("Analysis.name" %in% names(df)) df[["Analysis.name"]] else NA
  subgroup <- if ("Subgroup" %in% names(df)) df[["Subgroup"]] else NA
  paste(group, number, name, subgroup, sep = "|")
}

calc_effects <- function(df) {
  df$Experimental.cases <- as_num(df$Experimental.cases)
  df$Experimental.N <- as_num(df$Experimental.N)
  df$Control.cases <- as_num(df$Control.cases)
  df$Control.N <- as_num(df$Control.N)
  df$Experimental.mean <- as_num(df$Experimental.mean)
  df$Experimental.SD <- as_num(df$Experimental.SD)
  df$Control.mean <- as_num(df$Control.mean)
  df$Control.SD <- as_num(df$Control.SD)

  # Prefer reported effects (GIV/Mean/CI/Variance) to preserve original scale
  has_giv <- all(c("GIV.Mean", "GIV.SE") %in% names(df)) &&
    any(is.finite(as_num(df[["GIV.Mean"]]))) && any(is.finite(as_num(df[["GIV.SE"]])))
  has_mean <- "Mean" %in% names(df) && any(is.finite(as_num(df[["Mean"]])))

  if (has_giv || has_mean) {
    yi <- if (has_giv) as_num(df[["GIV.Mean"]]) else as_num(df[["Mean"]])
    se <- if (has_giv) as_num(df[["GIV.SE"]]) else {
      if ("Variance" %in% names(df) && any(is.finite(as_num(df[["Variance"]])))) {
        sqrt(as_num(df[["Variance"]]))
      } else if (all(c("CI.start", "CI.end") %in% names(df))) {
        (as_num(df[["CI.end"]]) - as_num(df[["CI.start"]])) / (2 * 1.96)
      } else {
        rep(NA_real_, nrow(df))
      }
    }
    use <- is.finite(yi) & is.finite(se) & se > 0
    d <- df[use, , drop = FALSE]
    if (nrow(d) > 0) {
      d$yi <- yi[use]
      d$vi <- se[use]^2
      d$effect_class <- "reported"
      d$outcome_type <- "reported"
      return(d)
    }
  }

  # Binary rows (fallback)
  bin_ok <- is.finite(df$Experimental.cases) & is.finite(df$Experimental.N) &
    is.finite(df$Control.cases) & is.finite(df$Control.N) &
    df$Experimental.N > 0 & df$Control.N > 0 &
    df$Experimental.cases >= 0 & df$Control.cases >= 0 &
    df$Experimental.cases <= df$Experimental.N &
    df$Control.cases <= df$Control.N

  bin <- df[bin_ok, , drop = FALSE]
  if (nrow(bin) > 0) {
    bi <- bin$Experimental.N - bin$Experimental.cases
    di <- bin$Control.N - bin$Control.cases
    bin_es <- metafor::escalc(measure = "OR",
                              ai = bin$Experimental.cases,
                              bi = bi,
                              ci = bin$Control.cases,
                              di = di,
                              add = 0.5, to = "all")
    bin$yi <- bin_es$yi
    bin$vi <- bin_es$vi
    bin$effect_class <- "or"
    bin$outcome_type <- "binary"
    return(bin)
  }

  # Continuous rows (fallback)
  cont_ok <- is.finite(df$Experimental.mean) & is.finite(df$Experimental.SD) &
    is.finite(df$Experimental.N) & is.finite(df$Control.mean) &
    is.finite(df$Control.SD) & is.finite(df$Control.N) &
    df$Experimental.N > 1 & df$Control.N > 1 &
    df$Experimental.SD > 0 & df$Control.SD > 0

  cont <- df[cont_ok, , drop = FALSE]
  if (nrow(cont) > 0) {
    cont_es <- metafor::escalc(measure = "SMD",
                               m1i = cont$Experimental.mean,
                               sd1i = cont$Experimental.SD,
                               n1i = cont$Experimental.N,
                               m2i = cont$Control.mean,
                               sd2i = cont$Control.SD,
                               n2i = cont$Control.N,
                               vtype = "UB")
    cont$yi <- cont_es$yi
    cont$vi <- cont_es$vi
    cont$effect_class <- "smd"
    cont$outcome_type <- "continuous"
    return(cont)
  }

  NULL
}

review_path <- file.path(output_dir, "transportability_review_level.csv")
rebuild_review <- Sys.getenv("REBUILD_REVIEW", "0") == "1"
review_df <- NULL

if (file.exists(review_path) && !rebuild_review) {
  review_df <- read.csv(review_path)
} else {
  rows <- list()
  for (f in files) {
    entry <- safe_load(f)
    if (is.null(entry)) next
    df <- entry$data
    if (!is.data.frame(df)) next

    df$dataset <- entry$name
    df$review_id <- sub("_pub.*$", "", entry$name)
    df$study_year <- if ("Study.year" %in% names(df)) as_num(df$Study.year) else NA_real_
    df$analysis_key <- make_analysis_key(df)

    groups <- split(df, df$analysis_key)
    for (g in groups) {
      effects <- calc_effects(g)
      if (is.null(effects)) next
      effects$analysis_key <- g$analysis_key[1]
      effects$analysis_name <- if ("Analysis.name" %in% names(g)) g$Analysis.name[1] else NA_character_
      effects$analysis_group <- if ("Analysis.group" %in% names(g)) g$Analysis.group[1] else NA_character_
      effects$analysis_number <- if ("Analysis.number" %in% names(g)) g$Analysis.number[1] else NA_character_
      effects$subgroup <- if ("Subgroup" %in% names(g)) g$Subgroup[1] else NA_character_
      rows[[length(rows) + 1]] <- effects
    }
  }

  all_cols <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, align_cols, cols = all_cols)
  all_effects <- do.call(rbind, rows)

  # Review-level meta-analysis per analysis_key and effect class
  review_summaries <- list()
  keys <- unique(all_effects[, c("review_id", "analysis_key", "effect_class")])
  for (i in seq_len(nrow(keys))) {
    rid <- keys$review_id[i]
    key <- keys$analysis_key[i]
    cls <- keys$effect_class[i]
    df <- all_effects[all_effects$review_id == rid &
        all_effects$analysis_key == key &
        all_effects$effect_class == cls, , drop = FALSE]
    if (nrow(df) < 2) next

    fit <- try(metafor::rma.uni(yi = df$yi, vi = df$vi, method = "REML"), silent = TRUE)
    if (inherits(fit, "try-error")) next

    ctrl_risk <- NA_real_
    ok <- is.finite(df$Control.cases) & is.finite(df$Control.N) & df$Control.N > 0
    if (any(ok)) ctrl_risk <- mean(df$Control.cases[ok] / df$Control.N[ok])

    total_n <- NA_real_
    if ("Experimental.N" %in% names(df) && "Control.N" %in% names(df)) {
      total_n <- mean(df$Experimental.N + df$Control.N, na.rm = TRUE)
    }

    review_summaries[[length(review_summaries) + 1]] <- data.frame(
      review_id = rid,
      dataset = df$dataset[1],
      analysis_key = key,
      analysis_name = df$analysis_name[1],
      analysis_group = df$analysis_group[1],
      analysis_number = df$analysis_number[1],
      subgroup = df$subgroup[1],
      effect_class = cls,
      outcome_type = df$outcome_type[1],
      k = nrow(df),
      yi = fit$b[1],
      se = fit$se[1],
      tau2 = fit$tau2,
      mean_year = mean(df$study_year, na.rm = TRUE),
      mean_ctrl_risk = ctrl_risk,
      mean_total_n = total_n,
      stringsAsFactors = FALSE
    )
  }

  review_df <- do.call(rbind, review_summaries)
  write.csv(review_df, review_path, row.names = FALSE)
}

# Align predictions when model drops columns
predict_from_fit <- function(fit, mods, data) {
  X <- model.matrix(mods, data = data)
  colnames(X)[colnames(X) == "(Intercept)"] <- "intrcpt"
  beta <- as.numeric(fit$beta)
  names(beta) <- rownames(fit$beta)
  X_use <- matrix(0, nrow = nrow(X), ncol = length(beta))
  colnames(X_use) <- names(beta)
  shared <- intersect(colnames(X), colnames(X_use))
  X_use[, shared] <- X[, shared, drop = FALSE]
  as.numeric(X_use %*% beta)
}

# Transportability CV (leave-one-review-out) per effect type
cv_rows <- list()
k_folds <- 10
set.seed(1)
for (cls in unique(review_df$effect_class)) {
  df <- review_df[review_df$effect_class == cls, , drop = FALSE]

  if (cls == "or") {
    mods <- ~ mean_year + mean_ctrl_risk + mean_total_n
    covars <- c("mean_year", "mean_ctrl_risk", "mean_total_n")
  } else {
    mods <- ~ mean_year + mean_total_n
    covars <- c("mean_year", "mean_total_n")
  }

  keep <- complete.cases(df[, c("yi", "se", covars), drop = FALSE])
  df <- df[keep, , drop = FALSE]
  if (nrow(df) < 5) next

  review_ids <- unique(df$review_id)
  if (length(review_ids) < 3) next
  folds <- sample(rep(seq_len(min(k_folds, length(review_ids))), length.out = length(review_ids)))
  for (k in seq_len(max(folds))) {
    test_ids <- review_ids[folds == k]
    train <- df[!(df$review_id %in% test_ids), , drop = FALSE]
    test <- df[df$review_id %in% test_ids, , drop = FALSE]
    if (nrow(test) == 0 || nrow(train) < 5) next

    fit <- try(metafor::rma.mv(yi = train$yi, V = train$se^2, mods = mods,
                               random = ~1 | review_id, data = train, method = "REML"), silent = TRUE)
    if (inherits(fit, "try-error")) next

    pred <- predict_from_fit(fit, mods, test)
    cv_rows[[length(cv_rows) + 1]] <- data.frame(
      review_id = test$review_id,
      effect_class = cls,
      observed = test$yi,
      predicted = pred,
      error = pred - test$yi,
      stringsAsFactors = FALSE
    )
  }
}

cv_df <- if (length(cv_rows) > 0) do.call(rbind, cv_rows) else data.frame()
cv_path <- file.path(output_dir, "transportability_cv_results.csv")
write.csv(cv_df, cv_path, row.names = FALSE)

# CV summary
summary_path <- file.path(output_dir, "transportability_cv_summary.md")
summary_lines <- c(
  "# Transportability CV Summary",
  "",
  sprintf("Reviews modeled: %d", length(unique(review_df$review_id))),
  "",
  "## CV metrics (by effect class)"
)

if (!is.null(cv_df) && nrow(cv_df) > 0) {
  for (cls in unique(cv_df$effect_class)) {
    sub <- cv_df[cv_df$effect_class == cls, , drop = FALSE]
    rmse <- sqrt(mean(sub$error^2, na.rm = TRUE))
    mae <- mean(abs(sub$error), na.rm = TRUE)
    summary_lines <- c(
      summary_lines,
      sprintf("- %s: RMSE = %.3f, MAE = %.3f (n = %d)", cls, rmse, mae, nrow(sub))
    )
  }
} else {
  summary_lines <- c(summary_lines, "- No CV results (insufficient reviews)")
}

writeLines(summary_lines, summary_path)

cat("Wrote:", review_path, "\n")
cat("Wrote:", cv_path, "\n")
cat("Wrote:", summary_path, "\n")
