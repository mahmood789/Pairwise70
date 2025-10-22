# Automated conversion of all 521 Cochrane pairwise reviews to RDA format
# This replaces the manual 23-dataset script with full automation

library(tools)

setwd('C:/Users/user/OneDrive - NHS/Documents/Pairwise70')

# Discover all pairwise CSV files
pairwise_dir <- 'C:/Users/user/OneDrive - NHS/Documents/CochraneDataExtractor/data/pairwise'
all_csv_files <- list.files(pairwise_dir, pattern = "data-rows\\.csv$", full.names = TRUE)

cat(paste0("=== COCHRANE PAIRWISE DATA CONVERSION ===\n"))
cat(paste0("Found ", length(all_csv_files), " CSV files to convert\n\n"))

# Track statistics
converted <- 0
skipped <- 0
errors <- 0

for (csv_file in all_csv_files) {
  # Extract dataset name from filename
  # Example: 10_1002_14651858_CD002042_pub6_CD002042-data-rows.csv -> CD002042_pub6_data
  filename <- basename(csv_file)

  # Extract CD number and publication version
  parts <- strsplit(filename, "_")[[1]]

  # Find CD number (starts with CD)
  cd_idx <- which(grepl("^CD", parts))

  if (length(cd_idx) == 0) {
    cat(paste0("SKIP: ", filename, " (no CD number found)\n"))
    skipped <- skipped + 1
    next
  }

  cd_num <- parts[cd_idx[1]]

  # Check if there's a pub version
  pub_idx <- cd_idx[1] + 1
  if (pub_idx <= length(parts) && grepl("^pub[0-9]", parts[pub_idx])) {
    dataset_name <- paste0(cd_num, "_", parts[pub_idx], "_data")
  } else {
    dataset_name <- paste0(cd_num, "_data")
  }

  # Check if already exists
  rda_path <- file.path('data', paste0(dataset_name, '.rda'))
  if (file.exists(rda_path)) {
    cat(paste0("EXISTS: ", dataset_name, "\n"))
    skipped <- skipped + 1
    next
  }

  tryCatch({
    # Read CSV
    data <- read.csv(csv_file, stringsAsFactors = FALSE)

    # Skip if empty
    if (nrow(data) == 0) {
      cat(paste0("SKIP: ", dataset_name, " (empty dataset)\n"))
      skipped <- skipped + 1
      next
    }

    # Assign to variable name
    assign(dataset_name, data)

    # Save as RDA
    save(list = dataset_name, file = rda_path)

    cat(paste0("CONVERTED: ", dataset_name, " (", nrow(data), " rows)\n"))
    converted <- converted + 1

  }, error = function(e) {
    cat(paste0("ERROR: ", dataset_name, " - ", e$message, "\n"))
    errors <- errors + 1
  })
}

cat("\n=== CONVERSION SUMMARY ===\n")
cat(paste0("Total CSV files found: ", length(all_csv_files), "\n"))
cat(paste0("Successfully converted: ", converted, "\n"))
cat(paste0("Already existed: ", skipped, "\n"))
cat(paste0("Errors: ", errors, "\n"))
cat(paste0("Total RDA files: ", length(list.files('data', pattern = '\\.rda$')), "\n"))
cat("\nConversion complete!\n")
