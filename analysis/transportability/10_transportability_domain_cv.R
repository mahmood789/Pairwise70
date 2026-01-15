# Domain-stratified transportability CV using registry covariates
suppressWarnings({
  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop("Package 'metafor' is required. Install with install.packages('metafor').")
  }
})

root <- normalizePath("C:/Users/user/OneDrive - NHS/Documents/Pairwise70", winslash = "/", mustWork = TRUE)
output_dir <- file.path(root, "analysis", "transportability")

merge_path <- file.path(output_dir, "transportability_target_merge.csv")
terms_path <- file.path(output_dir, "pairwise70_analysis_terms.csv")
if (!file.exists(merge_path)) stop("Missing transportability_target_merge.csv")
if (!file.exists(terms_path)) stop("Missing pairwise70_analysis_terms.csv")

merged <- read.csv(merge_path, stringsAsFactors = FALSE)
terms <- read.csv(terms_path, stringsAsFactors = FALSE)
terms <- terms[terms$source %in% c("Analysis.name", "Subgroup"), , drop = FALSE]
terms$review_id <- sub("_pub.*$", "", terms$dataset)

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

normalize_text <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9 ]", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

build_pattern <- function(keys) {
  paste0("\\b(", paste(keys, collapse = "|"), ")\\b")
}

categories <- list(
  cardiovascular = c("cardiac", "cardio", "heart", "coronary", "angina", "myocard", "ischemi",
                     "arrhythm", "hypertens", "blood pressure", "cardiomyopathy", "heart failure",
                     "atrial fibrillation", "af", "valvular", "anticoag", "thrombo", "embol"),
  cerebrovascular = c("stroke", "cerebro", "transient ischemic", "tia", "intracran", "carotid"),
  respiratory = c("asthma", "copd", "respir", "lung", "pulmonary", "pneumonia", "bronch"),
  infectious = c("infection", "infect", "sepsis", "hiv", "hepatitis", "tb", "tuberculosis", "influenza",
                 "covid", "malaria", "antibiotic", "antiviral"),
  oncology = c("cancer", "tumor", "tumour", "carcinoma", "oncology", "chemotherapy", "radiation"),
  gastrointestinal = c("gastro", "hepatic", "liver", "cirrhosis", "ulcer", "colitis", "crohn",
                       "pancreat", "bowel", "diarr"),
  musculoskeletal = c("arthritis", "oste", "bone", "fracture", "joint", "rheumat", "spine", "back pain"),
  mental_health = c("depress", "anxiety", "schizo", "bipolar", "psych", "mental", "ptsd"),
  maternal_neonatal = c("pregnan", "maternal", "neonat", "infant", "neonate", "birth", "labor", "caesarean"),
  renal = c("renal", "kidney", "dialysis", "nephro"),
  endocrine_metabolic = c("diabetes", "glyc", "insulin", "thyroid", "obesity", "metabolic"),
  neurology = c("neuro", "epilepsy", "seizure", "parkinson", "alzheimer", "dementia", "headache",
                "migraine", "multiple sclerosis"),
  dermatology = c("skin", "dermat", "psoriasis", "eczema", "rash"),
  ophthalmology = c("eye", "ocular", "glaucoma", "retina", "macular", "vision"),
  dentistry = c("dental", "tooth", "teeth", "oral", "periodont"),
  urology = c("urinary", "uro", "bladder", "prostat", "kidney stone"),
  pain = c("pain", "analges", "opioid")
)

patterns <- lapply(categories, build_pattern)

tag_patterns <- list(
  mortality = build_pattern(c("mortality", "death", "survival"))
)

score_review <- function(texts) {
  counts <- setNames(rep(0, length(categories)), names(categories))
  tags <- setNames(rep(FALSE, length(tag_patterns)), names(tag_patterns))
  for (txt in texts) {
    txt <- normalize_text(txt)
    if (txt == "") next
    for (cat in names(patterns)) {
      if (grepl(patterns[[cat]], txt)) counts[cat] <- counts[cat] + 1
    }
    for (tag in names(tag_patterns)) {
      if (grepl(tag_patterns[[tag]], txt)) tags[tag] <- TRUE
    }
  }
  list(counts = counts, tags = tags)
}

review_ids <- unique(terms$review_id)
review_domains <- data.frame(
  review_id = review_ids,
  domain_primary = "other",
  domain_secondary = NA_character_,
  domain_list = NA_character_,
  multi_label = FALSE,
  mortality_tag = FALSE,
  top_score = 0,
  second_score = 0,
  stringsAsFactors = FALSE
)

for (i in seq_along(review_ids)) {
  rid <- review_ids[i]
  texts <- terms$term[terms$review_id == rid]
  scored <- score_review(texts)
  scores <- scored$counts
  top_val <- max(scores)
  review_domains$mortality_tag[i] <- isTRUE(scored$tags["mortality"])
  if (top_val == 0) {
    review_domains$domain_primary[i] <- "other"
    review_domains$domain_list[i] <- "other"
    next
  }
  ordered <- sort(scores, decreasing = TRUE)
  top_cat <- names(ordered)[1]
  second_val <- if (length(ordered) > 1) ordered[2] else 0

  review_domains$top_score[i] <- top_val
  review_domains$second_score[i] <- second_val

  threshold <- max(2, ceiling(top_val * 0.6))
  candidates <- names(scores)[scores >= threshold]
  if (!(top_cat %in% candidates)) candidates <- c(top_cat, candidates)
  candidates <- unique(candidates)

  review_domains$domain_primary[i] <- top_cat
  review_domains$domain_list[i] <- paste(candidates, collapse = ";")
  review_domains$multi_label[i] <- length(candidates) > 1
  if (length(candidates) > 1) {
    review_domains$domain_secondary[i] <- candidates[2]
  }
}

merged <- merge(merged, review_domains, by = "review_id", all.x = TRUE)
merged$domain_primary[is.na(merged$domain_primary)] <- "other"
merged$domain_list[is.na(merged$domain_list)] <- "other"
merged$multi_label[is.na(merged$multi_label)] <- FALSE

count_path <- file.path(output_dir, "transportability_domain_counts.csv")
count_df <- as.data.frame(table(merged$domain_primary), stringsAsFactors = FALSE)
colnames(count_df) <- c("domain", "rows")
review_counts <- aggregate(review_id ~ domain_primary, data = review_domains, FUN = function(x) length(unique(x)))
colnames(review_counts) <- c("domain", "reviews")
count_df <- merge(count_df, review_counts, by = "domain", all.x = TRUE)

multi_counts <- as.data.frame(table(unlist(strsplit(review_domains$domain_list, ";"))), stringsAsFactors = FALSE)
colnames(multi_counts) <- c("domain", "reviews_in_list")
count_df <- merge(count_df, multi_counts, by = "domain", all.x = TRUE)
write.csv(count_df, count_path, row.names = FALSE)

assign_path <- file.path(output_dir, "transportability_domain_assignments.csv")
write.csv(review_domains, assign_path, row.names = FALSE)

cv_rows <- list()
set.seed(1)

for (domain in unique(merged$domain_primary)) {
  for (cls in unique(merged$effect_class)) {
    df <- merged[merged$domain_primary == domain & merged$effect_class == cls, , drop = FALSE]

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

    k_folds <- min(5, length(review_ids))
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
        domain = domain,
        effect_class = cls,
        review_id = test$review_id,
        observed = test$yi,
        predicted = pred,
        error = pred - test$yi,
        stringsAsFactors = FALSE
      )
    }
  }
}

cv_df <- if (length(cv_rows) > 0) do.call(rbind, cv_rows) else data.frame()
cv_path <- file.path(output_dir, "transportability_domain_cv_results.csv")
write.csv(cv_df, cv_path, row.names = FALSE)

summary_path <- file.path(output_dir, "transportability_domain_cv_summary.md")
summary_lines <- c(
  "# Transportability Domain CV Summary",
  "",
  sprintf("Rows modeled: %d", nrow(cv_df)),
  "",
  "## Domain counts",
  paste0("- ", count_df$domain, ": ", count_df$rows, " rows, ", count_df$reviews, " reviews, ", count_df$reviews_in_list, " reviews_in_list"),
  "",
  "## CV metrics (by domain and effect class)"
)

if (!is.null(cv_df) && nrow(cv_df) > 0) {
  domains <- unique(cv_df$domain)
  for (d in domains) {
    sub_d <- cv_df[cv_df$domain == d, , drop = FALSE]
    for (cls in unique(sub_d$effect_class)) {
      sub <- sub_d[sub_d$effect_class == cls, , drop = FALSE]
      rmse <- sqrt(mean(sub$error^2, na.rm = TRUE))
      mae <- mean(abs(sub$error), na.rm = TRUE)
      summary_lines <- c(
        summary_lines,
        sprintf("- %s | %s: RMSE = %.3f, MAE = %.3f (n = %d)", d, cls, rmse, mae, nrow(sub))
      )
    }
  }
} else {
  summary_lines <- c(summary_lines, "- No CV results (insufficient data)")
}

writeLines(summary_lines, summary_path)

cat("Wrote:", count_path, "\n")
cat("Wrote:", assign_path, "\n")
cat("Wrote:", cv_path, "\n")
cat("Wrote:", summary_path, "\n")
