# Fetch ClinicalTrials.gov data for each review query
suppressWarnings({
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required. Install with install.packages('jsonlite').")
  }
})

root <- normalizePath("C:/Users/user/OneDrive - NHS/Documents/Pairwise70", winslash = "/", mustWork = TRUE)
output_dir <- file.path(root, "analysis", "transportability")
query_path <- file.path(output_dir, "ctgov_query_terms.csv")

external_dir <- file.path(output_dir, "external", "ctgov")
raw_dir <- file.path(external_dir, "raw")
if (!dir.exists(raw_dir)) dir.create(raw_dir, recursive = TRUE)

queries <- read.csv(query_path, stringsAsFactors = FALSE)

start_index <- as.integer(Sys.getenv("CTGOV_START_INDEX", 1))
if (is.na(start_index) || start_index < 1) start_index <- 1

max_reviews_raw <- Sys.getenv("CTGOV_MAX_REVIEWS", "")
max_reviews <- if (nzchar(max_reviews_raw)) as.integer(max_reviews_raw) else nrow(queries)
if (is.na(max_reviews) || max_reviews < 1) max_reviews <- nrow(queries)

max_pages <- as.integer(Sys.getenv("CTGOV_MAX_PAGES", 3))
page_size <- as.integer(Sys.getenv("CTGOV_PAGE_SIZE", 100))
sleep_ms <- as.integer(Sys.getenv("CTGOV_SLEEP_MS", 0))
if (is.na(sleep_ms) || sleep_ms < 0) sleep_ms <- 0

fields <- c(
  "protocolSection.identificationModule.nctId",
  "protocolSection.identificationModule.briefTitle",
  "protocolSection.conditionsModule.conditions",
  "protocolSection.statusModule.startDateStruct",
  "protocolSection.statusModule.completionDateStruct",
  "protocolSection.designModule.studyType",
  "protocolSection.designModule.enrollmentInfo",
  "protocolSection.eligibilityModule.sex",
  "protocolSection.eligibilityModule.minimumAge",
  "protocolSection.eligibilityModule.maximumAge",
  "protocolSection.eligibilityModule.healthyVolunteers"
)

fetch_query <- function(query_term, max_pages = 3L, page_size = 100L) {
  base <- "https://clinicaltrials.gov/api/v2/studies"
  all_studies <- list()
  page_token <- NULL

  for (i in seq_len(max_pages)) {
    url <- paste0(
      base,
      "?query.term=", utils::URLencode(query_term, reserved = TRUE),
      "&pageSize=", page_size,
      "&fields=", utils::URLencode(paste(fields, collapse = ","), reserved = TRUE)
    )
    if (!is.null(page_token)) {
      url <- paste0(url, "&pageToken=", utils::URLencode(page_token, reserved = TRUE))
    }

    resp <- tryCatch(jsonlite::fromJSON(url), error = function(e) NULL)
    if (is.null(resp)) {
      url_fallback <- paste0(
        base,
        "?query.term=", utils::URLencode(query_term, reserved = TRUE),
        "&pageSize=", page_size
      )
      if (!is.null(page_token)) {
        url_fallback <- paste0(url_fallback, "&pageToken=", utils::URLencode(page_token, reserved = TRUE))
      }
      resp <- tryCatch(jsonlite::fromJSON(url_fallback), error = function(e) NULL)
    }
    if (is.null(resp)) {
      message("Failed to fetch query: ", query_term)
      break
    }
    if (!is.null(resp$studies)) {
      all_studies <- c(all_studies, resp$studies)
    }

    if (is.null(resp$nextPageToken)) break
    page_token <- resp$nextPageToken
  }

  list(query = query_term, studies = all_studies)
}

if (start_index > nrow(queries)) {
  message("CTGOV_START_INDEX exceeds number of reviews. Exiting.")
  quit(status = 0)
}

end_index <- min(nrow(queries), start_index + max_reviews - 1L)
for (i in seq(from = start_index, to = end_index)) {
  rid <- queries$review_id[i]
  query <- queries$query_term[i]
  if (is.na(query) || trimws(query) == "") next

  out_path <- file.path(raw_dir, paste0(rid, ".json"))
  if (file.exists(out_path)) next

  payload <- fetch_query(query, max_pages = max_pages, page_size = page_size)
  jsonlite::write_json(payload, out_path, auto_unbox = TRUE, pretty = TRUE)
  message("Saved: ", out_path)
  if (sleep_ms > 0) {
    Sys.sleep(sleep_ms / 1000)
  }
}
