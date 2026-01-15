# Convert Pairwise70 RDA files to CSV format for Python usage
# =============================================================

library(data.table)

# Set working directory to data folder
setwd("C:/Users/user/OneDrive - NHS/Documents/Pairwise70/data")

# Get all RDA files
rda_files <- list.files(pattern = "\\.rda$", full.names = FALSE)

cat(sprintf("Found %d RDA files to convert\n\n", length(rda_files)))

# Create index data
index_data <- data.table(
  dataset_name = character(),
  n_studies = integer(),
  outcome_type = character(),
  review_doi = character(),
  stringsAsFactors = FALSE
)

# Convert each file
for (i in seq_along(rda_files)) {
  rda_file <- rda_files[i]

  tryCatch({
    # Load the RDA file
    load(rda_file)

    # Get the dataset name
    dataset_name <- gsub("\\.rda$", "", rda_file)

    # Get the data object (the variable with the same name as the file)
    data_obj <- get(dataset_name)

    # Ensure it's a data.table
    if (!is.data.table(data_obj)) {
      data_obj <- as.data.table(data_obj)
    }

    # Get dimensions
    n_studies <- nrow(data_obj)
    n_cols <- ncol(data_obj)

    # Determine outcome type
    outcome_type <- "unknown"
    if ("Experimental.cases" %in% names(data_obj) &&
        "Control.cases" %in% names(data_obj)) {
      outcome_type <- "binary"
    } else if ("Experimental.mean" %in% names(data_obj) &&
               "Control.mean" %in% names(data_obj)) {
      outcome_type <- "continuous"
    } else if ("Mean" %in% names(data_obj)) {
      outcome_type <- "effect_size"
    }

    # Get DOI
    review_doi <- "NA"
    if ("review_doi" %in% names(data_obj)) {
      review_doi <- as.character(data_obj$review_doi[1])
    }

    # Write to CSV
    csv_file <- gsub("\\.rda$", ".csv", rda_file)
    fwrite(data_obj, csv_file)

    # Add to index
    index_data <- rbind(index_data, data.table(
      dataset_name = dataset_name,
      n_studies = n_studies,
      outcome_type = outcome_type,
      review_doi = review_doi
    ))

    cat(sprintf("[%d/%d] Converted: %s (%d studies, %s)\n",
                i, length(rda_files), dataset_name, n_studies, outcome_type))

  }, error = function(e) {
    cat(sprintf("Error converting %s: %s\n", rda_file, e$message))
  })
}

# Save index
setwd("C:/Users/user/OneDrive - NHS/Documents/Pairwise70")
fwrite(index_data, "data_index.csv")

cat(sprintf("\n=== Conversion Complete ===\n"))
cat(sprintf("Total datasets: %d\n", nrow(index_data)))
cat(sprintf("Binary outcomes: %d\n", sum(index_data$outcome_type == "binary")))
cat(sprintf("Continuous outcomes: %d\n", sum(index_data$outcome_type == "continuous")))
cat(sprintf("Effect size data: %d\n", sum(index_data$outcome_type == "effect_size")))
cat(sprintf("\nIndex saved to: C:/Users/user/OneDrive - NHS/Documents/Pairwise70/data_index.csv\n"))
