# Build ClinicalTrials.gov query terms from Pairwise70 Analysis.name
root <- normalizePath("C:/Users/user/OneDrive - NHS/Documents/Pairwise70", winslash = "/", mustWork = TRUE)
output_dir <- file.path(root, "analysis", "transportability")
terms_path <- file.path(output_dir, "pairwise70_analysis_terms.csv")

terms <- read.csv(terms_path, stringsAsFactors = FALSE)
terms <- terms[terms$source == "Analysis.name", , drop = FALSE]

stopwords <- c("with", "without", "trial", "trials", "study", "studies", "randomized", "randomised",
               "patients", "patient", "disease", "treatment", "therapy", "versus", "vs", "using",
               "effect", "effects", "outcome", "outcomes", "rate", "rates", "risk", "group", "groups",
               "change", "changes", "prevention", "management", "drug", "drugs", "acute", "chronic",
               "compared", "comparison", "dose", "doses", "daily", "placebo", "control", "controlled",
               "mortality", "survival", "efficacy", "safety", "events", "event", "overall")

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

review_id <- sub("_pub.*$", "", terms$dataset)
terms$review_id <- review_id

query_rows <- list()
for (rid in unique(terms$review_id)) {
  vals <- terms$term[terms$review_id == rid]
  tokens <- unlist(lapply(vals, normalize_terms))
  if (length(tokens) == 0) next
  freq <- sort(table(tokens), decreasing = TRUE)
  top_tokens <- names(freq)
  query_rows[[length(query_rows) + 1]] <- data.frame(
    review_id = rid,
    top_tokens = paste(top_tokens[seq_len(min(12, length(top_tokens)))], collapse = ", "),
    query_term = build_query(top_tokens, n = 4),
    query_term_short = build_query(top_tokens, n = 2),
    stringsAsFactors = FALSE
  )
}

query_df <- do.call(rbind, query_rows)
query_path <- file.path(output_dir, "ctgov_query_terms.csv")
write.csv(query_df, query_path, row.names = FALSE)

# Global token list
all_tokens <- unlist(lapply(terms$term, normalize_terms))
all_freq <- sort(table(all_tokens), decreasing = TRUE)
all_df <- data.frame(token = names(all_freq), count = as.integer(all_freq), stringsAsFactors = FALSE)
write.csv(all_df, file.path(output_dir, "ctgov_global_tokens.csv"), row.names = FALSE)

cat("Wrote:", query_path, "\n")
cat("Wrote:", file.path(output_dir, "ctgov_global_tokens.csv"), "\n")
