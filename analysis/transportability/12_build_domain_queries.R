# Build ClinicalTrials.gov query terms by domain
root <- normalizePath("C:/Users/user/OneDrive - NHS/Documents/Pairwise70", winslash = "/", mustWork = TRUE)
output_dir <- file.path(root, "analysis", "transportability")
terms_path <- file.path(output_dir, "pairwise70_analysis_terms.csv")
assign_path <- file.path(output_dir, "transportability_domain_assignments.csv")

terms <- read.csv(terms_path, stringsAsFactors = FALSE)
terms <- terms[terms$source %in% c("Analysis.name", "Subgroup"), , drop = FALSE]
terms$review_id <- sub("_pub.*$", "", terms$dataset)

assign <- read.csv(assign_path, stringsAsFactors = FALSE)

stopwords <- c("with", "without", "trial", "trials", "study", "studies", "randomized", "randomised",
               "patients", "patient", "disease", "treatment", "therapy", "versus", "vs", "using",
               "effect", "effects", "outcome", "outcomes", "rate", "rates", "risk", "group", "groups",
               "change", "changes", "prevention", "management", "drug", "drugs", "acute", "chronic",
               "compared", "comparison", "dose", "doses", "daily", "placebo", "control", "controlled",
               "mortality", "survival", "overall", "adverse")

normalize_terms <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z ]", " ", x)
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  if (x == "") return(character(0))
  tokens <- unlist(strsplit(x, " ", fixed = TRUE))
  tokens <- tokens[nchar(tokens) >= 4]
  tokens <- tokens[!tokens %in% stopwords]
  tokens
}

build_query <- function(tokens, n = 4) {
  if (length(tokens) == 0) return(NA_character_)
  tokens <- tokens[seq_len(min(n, length(tokens)))]
  paste(tokens, collapse = " ")
}

assign$review_id <- as.character(assign$review_id)
terms$review_id <- as.character(terms$review_id)

queries <- list()
for (domain in unique(assign$domain_primary)) {
  reviews <- assign$review_id[assign$domain_primary == domain]
  domain_terms <- terms$term[terms$review_id %in% reviews]
  tokens <- unlist(lapply(domain_terms, normalize_terms))
  if (length(tokens) == 0) next

  freq <- sort(table(tokens), decreasing = TRUE)
  top_tokens <- names(freq)

  queries[[length(queries) + 1]] <- data.frame(
    domain = domain,
    review_count = length(unique(reviews)),
    top_tokens = paste(top_tokens[seq_len(min(15, length(top_tokens)))], collapse = ", "),
    query_term = build_query(top_tokens, n = 4),
    query_term_short = build_query(top_tokens, n = 2),
    stringsAsFactors = FALSE
  )
}

query_df <- do.call(rbind, queries)
query_path <- file.path(output_dir, "ctgov_domain_query_terms.csv")
write.csv(query_df, query_path, row.names = FALSE)

cat("Wrote:", query_path, "\n")
