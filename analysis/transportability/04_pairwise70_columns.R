# Column coverage and analysis terms for Pairwise70
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

`%||%` <- function(a, b) if (is.null(a)) b else a

col_counts <- list()
term_rows <- list()

for (f in files) {
  entry <- safe_load(f)
  if (is.null(entry)) next
  df <- entry$data
  if (!is.data.frame(df)) next

  cols <- names(df)
  for (c in cols) {
    col_counts[[c]] <- (col_counts[[c]] %||% 0) + 1
  }

  analysis_group <- if ("Analysis.group" %in% cols) unique(na.omit(df$Analysis.group)) else character(0)
  analysis_name <- if ("Analysis.name" %in% cols) unique(na.omit(df$Analysis.name)) else character(0)
  subgroup <- if ("Subgroup" %in% cols) unique(na.omit(df$Subgroup)) else character(0)

  add_terms <- function(vals, source) {
    if (length(vals) == 0) return()
    vals <- vals[vals != ""]
    if (length(vals) == 0) return()
    term_rows[[length(term_rows) + 1]] <<- data.frame(
      dataset = entry$name,
      source = source,
      term = as.character(vals),
      stringsAsFactors = FALSE
    )
  }

  add_terms(analysis_group, "Analysis.group")
  add_terms(analysis_name, "Analysis.name")
  add_terms(subgroup, "Subgroup")
}

col_df <- data.frame(
  column = names(col_counts),
  datasets = as.integer(unlist(col_counts)),
  stringsAsFactors = FALSE
)
col_df <- col_df[order(-col_df$datasets, col_df$column), , drop = FALSE]
write.csv(col_df, file.path(output_dir, "pairwise70_column_coverage.csv"), row.names = FALSE)

if (length(term_rows) > 0) {
  terms <- do.call(rbind, term_rows)
  write.csv(terms, file.path(output_dir, "pairwise70_analysis_terms.csv"), row.names = FALSE)
}

cat("Wrote column coverage and term extracts\n")
