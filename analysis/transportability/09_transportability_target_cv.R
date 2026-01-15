# Transportability model with registry covariates (10-fold CV)
suppressWarnings({
  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop("Package 'metafor' is required. Install with install.packages('metafor').")
  }
})

root <- normalizePath("C:/Users/user/OneDrive - NHS/Documents/Pairwise70", winslash = "/", mustWork = TRUE)
output_dir <- file.path(root, "analysis", "transportability")

merge_path <- file.path(output_dir, "transportability_target_merge.csv")
if (!file.exists(merge_path)) {
  stop("Missing transportability_target_merge.csv. Run 08_transportability_merge.R first.")
}

merged <- read.csv(merge_path, stringsAsFactors = FALSE)

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

cv_rows <- list()
set.seed(1)

for (cls in unique(merged$effect_class)) {
  df <- merged[merged$effect_class == cls, , drop = FALSE]

  if (cls == "or") {
    covars <- c("mean_year", "mean_total_n", "mean_ctrl_risk",
                "start_year_mean", "enrollment_mean", "sex_all_pct")
  } else {
    covars <- c("mean_year", "mean_total_n", "start_year_mean", "enrollment_mean", "sex_all_pct")
  }

  keep <- complete.cases(df[, c("yi", "se", covars), drop = FALSE])
  df <- df[keep, , drop = FALSE]
  if (nrow(df) < 10) next

  review_ids <- unique(df$review_id)
  if (length(review_ids) < 3) next

  k_folds <- min(10, length(review_ids))
  folds <- sample(rep(seq_len(k_folds), length.out = length(review_ids)))
  mods <- as.formula(paste("~", paste(covars, collapse = " + ")))

  for (k in seq_len(k_folds)) {
    test_ids <- review_ids[folds == k]
    train <- df[!(df$review_id %in% test_ids), , drop = FALSE]
    test <- df[df$review_id %in% test_ids, , drop = FALSE]
    if (nrow(test) == 0 || nrow(train) < 10) next

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
cv_path <- file.path(output_dir, "transportability_target_cv_results.csv")
write.csv(cv_df, cv_path, row.names = FALSE)

summary_path <- file.path(output_dir, "transportability_target_cv_summary.md")
summary_lines <- c(
  "# Transportability Target CV Summary",
  "",
  sprintf("Rows modeled: %d", nrow(cv_df)),
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
  summary_lines <- c(summary_lines, "- No CV results (insufficient data)")
}

writeLines(summary_lines, summary_path)

cat("Wrote:", cv_path, "\n")
cat("Wrote:", summary_path, "\n")
