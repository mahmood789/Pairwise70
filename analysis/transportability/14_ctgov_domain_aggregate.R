# Aggregate ClinicalTrials.gov domain data into target covariates
suppressWarnings({
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required. Install with install.packages('jsonlite').")
  }
})

root <- normalizePath("C:/Users/user/OneDrive - NHS/Documents/Pairwise70", winslash = "/", mustWork = TRUE)
output_dir <- file.path(root, "analysis", "transportability")
external_dir <- file.path(output_dir, "external", "ctgov", "domain_raw")

json_files <- list.files(external_dir, pattern = "\\.json$", full.names = TRUE)

as_num <- function(x) suppressWarnings(as.numeric(x))

weighted_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(mean(x, na.rm = TRUE))
  sum(x[ok] * w[ok]) / sum(w[ok])
}

parse_year <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  x <- as.character(x[1])
  if (is.na(x) || nchar(x) < 4) return(NA_real_)
  as_num(substr(x, 1, 4))
}

flatten_studies <- function(studies) {
  if (is.data.frame(studies)) {
    return(split(studies, seq_len(nrow(studies))))
  }
  if (is.list(studies)) {
    if (length(studies) == 0) return(list())
    if (all(vapply(studies, is.data.frame, logical(1)))) {
      all_cols <- unique(unlist(lapply(studies, names)))
      aligned <- lapply(studies, function(d) {
        missing <- setdiff(all_cols, names(d))
        if (length(missing) > 0) d[missing] <- NA
        d[, all_cols, drop = FALSE]
      })
      combined <- do.call(rbind, aligned)
      return(split(combined, seq_len(nrow(combined))))
    }
    return(studies)
  }
  list()
}

pick_field <- function(study, keys) {
  for (k in keys) {
    if (k %in% names(study)) return(study[[k]])
  }
  NULL
}

normalize_scalar <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_character_)
  if (is.list(x)) {
    if (length(x) == 0) return(NA_character_)
    x <- x[[1]]
  }
  if (length(x) == 0) return(NA_character_)
  as.character(x[1])
}

rows <- list()
for (f in json_files) {
  payload <- jsonlite::fromJSON(f, flatten = TRUE)
  domain <- tools::file_path_sans_ext(basename(f))
  studies <- payload$studies
  if (is.null(studies) || length(studies) == 0) next

  study_list <- flatten_studies(studies)
  if (length(study_list) == 0) next

  start_dates <- character(0)
  completion_dates <- character(0)
  enrollment <- numeric(0)
  sex <- character(0)
  min_age <- character(0)
  max_age <- character(0)

  for (study in study_list) {
    start_date <- normalize_scalar(pick_field(
      study,
      c("statusModule.startDateStruct.date", "protocolSection.statusModule.startDateStruct.date",
        "statusModule.startDateStruct", "protocolSection.statusModule.startDateStruct")
    ))
    completion_date <- normalize_scalar(pick_field(
      study,
      c("statusModule.completionDateStruct.date", "protocolSection.statusModule.completionDateStruct.date",
        "statusModule.completionDateStruct", "protocolSection.statusModule.completionDateStruct")
    ))
    enroll_val <- normalize_scalar(pick_field(
      study,
      c("designModule.enrollmentInfo.count", "protocolSection.designModule.enrollmentInfo.count")
    ))
    sex_val <- normalize_scalar(pick_field(
      study,
      c("eligibilityModule.sex", "protocolSection.eligibilityModule.sex")
    ))
    min_val <- normalize_scalar(pick_field(
      study,
      c("eligibilityModule.minimumAge", "protocolSection.eligibilityModule.minimumAge")
    ))
    max_val <- normalize_scalar(pick_field(
      study,
      c("eligibilityModule.maximumAge", "protocolSection.eligibilityModule.maximumAge")
    ))

    start_dates <- c(start_dates, start_date)
    completion_dates <- c(completion_dates, completion_date)
    enrollment <- c(enrollment, as_num(enroll_val))
    sex <- c(sex, sex_val)
    min_age <- c(min_age, min_val)
    max_age <- c(max_age, max_val)
  }

  start_year <- vapply(start_dates, parse_year, numeric(1))
  completion_year <- vapply(completion_dates, parse_year, numeric(1))
  weights <- ifelse(is.finite(enrollment) & enrollment > 0, enrollment, NA_real_)

  rows[[length(rows) + 1]] <- data.frame(
    domain = domain,
    trial_count = length(start_dates),
    start_year_mean = weighted_mean(start_year, weights),
    completion_year_mean = weighted_mean(completion_year, weights),
    enrollment_median = median(enrollment, na.rm = TRUE),
    enrollment_mean = mean(enrollment, na.rm = TRUE),
    sex_all_pct = weighted_mean(as.numeric(sex == "ALL"), weights),
    sex_female_pct = weighted_mean(as.numeric(sex == "FEMALE"), weights),
    sex_male_pct = weighted_mean(as.numeric(sex == "MALE"), weights),
    min_age_missing = mean(is.na(min_age) | min_age == "", na.rm = TRUE),
    max_age_missing = mean(is.na(max_age) | max_age == "", na.rm = TRUE),
    start_year_mean_unweighted = mean(start_year, na.rm = TRUE),
    completion_year_mean_unweighted = mean(completion_year, na.rm = TRUE),
    sex_all_pct_unweighted = mean(sex == "ALL", na.rm = TRUE),
    sex_female_pct_unweighted = mean(sex == "FEMALE", na.rm = TRUE),
    sex_male_pct_unweighted = mean(sex == "MALE", na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

agg <- do.call(rbind, rows)
agg_path <- file.path(output_dir, "ctgov_domain_covariates.csv")
write.csv(agg, agg_path, row.names = FALSE)

summary_path <- file.path(output_dir, "ctgov_domain_summary.md")
summary_lines <- c(
  "# ClinicalTrials.gov Domain Summary",
  "",
  paste0("- ", agg$domain, ": ", agg$trial_count, " trials")
)
writeLines(summary_lines, summary_path)

cat("Wrote:", agg_path, "\n")
cat("Wrote:", summary_path, "\n")
