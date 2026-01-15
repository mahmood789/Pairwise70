# Fetch ClinicalTrials.gov data for domain queries
suppressWarnings({
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required. Install with install.packages('jsonlite').")
  }
})

root <- normalizePath("C:/Users/user/OneDrive - NHS/Documents/Pairwise70", winslash = "/", mustWork = TRUE)
output_dir <- file.path(root, "analysis", "transportability")
query_path <- file.path(output_dir, "ctgov_domain_query_terms.csv")

external_dir <- file.path(output_dir, "external", "ctgov", "domain_raw")
if (!dir.exists(external_dir)) dir.create(external_dir, recursive = TRUE)

queries <- read.csv(query_path, stringsAsFactors = FALSE)

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

max_pages <- as.integer(Sys.getenv("CTGOV_DOMAIN_MAX_PAGES", 5))
page_size <- as.integer(Sys.getenv("CTGOV_DOMAIN_PAGE_SIZE", 200))
sleep_ms <- as.integer(Sys.getenv("CTGOV_DOMAIN_SLEEP_MS", 200))
use_short <- Sys.getenv("CTGOV_DOMAIN_USE_SHORT", "1") == "1"
force_refresh <- Sys.getenv("CTGOV_DOMAIN_FORCE", "0") == "1"

filter_domains <- Sys.getenv("CTGOV_DOMAIN_LIST", "")
if (nzchar(filter_domains)) {
  keep <- trimws(unlist(strsplit(filter_domains, ",")))
  queries <- queries[queries$domain %in% keep, , drop = FALSE]
}

fetch_query <- function(query_term, max_pages, page_size) {
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

for (i in seq_len(nrow(queries))) {
  domain <- queries$domain[i]
  query <- if (use_short && "query_term_short" %in% names(queries)) {
    queries$query_term_short[i]
  } else {
    queries$query_term[i]
  }
  if (is.na(query) || trimws(query) == "") next

  out_path <- file.path(external_dir, paste0(domain, ".json"))
  if (file.exists(out_path) && !force_refresh) next

  payload <- fetch_query(query, max_pages = max_pages, page_size = page_size)
  jsonlite::write_json(payload, out_path, auto_unbox = TRUE, pretty = TRUE)
  message("Saved: ", out_path)
  if (sleep_ms > 0) {
    Sys.sleep(sleep_ms / 1000)
  }
}
