# Auto-generate .Rd documentation files for all 501 Pairwise70 datasets
# Each documentation file will include:
# - Dataset description with Cochrane review DOI and title
# - Column descriptions (ALL columns introspected from actual data)
# - Format information
# - Source and provenance

library(tools)

cat("=== GENERATING DOCUMENTATION FOR ALL 501 DATASETS ===\n\n")

# Get all datasets from data/ directory
data_files <- list.files("data", pattern = "\\.rda$", full.names = FALSE)
all_datasets <- sub("\\.rda$", "", data_files)
cat(paste0("Total datasets to document: ", length(all_datasets), "\n\n"))

# Create man directory if it doesn't exist
man_dir <- "man"
if (!dir.exists(man_dir)) {
  dir.create(man_dir)
}

# Helper function to escape special Rd characters
escape_rd <- function(text) {
  # Escape % which is the comment character in Rd syntax
  text <- gsub("%", "\\\\%", text, fixed = TRUE)
  # Escape backslash (must come first to avoid double-escaping)
  # text <- gsub("\\", "\\\\", text, fixed = TRUE)  # Usually not needed in descriptions
  return(text)
}

# Function to get column description based on column name
get_column_description <- function(col_name) {
  desc_map <- list(
    "Study" = "Study identifier (typically author and year)",
    "Study.year" = "Publication year of the study",
    "Experimental.cases" = "Number of events in experimental/treatment group (binary outcomes)",
    "Experimental.N" = "Total participants in experimental/treatment group",
    "Control.cases" = "Number of events in control group (binary outcomes)",
    "Control.N" = "Total participants in control group",
    "Experimental.mean" = "Mean value in experimental group (continuous outcomes)",
    "Experimental.SD" = "Standard deviation in experimental group (continuous outcomes)",
    "Control.mean" = "Mean value in control group (continuous outcomes)",
    "Control.SD" = "Standard deviation in control group (continuous outcomes)",
    "Subgroup" = "Subgroup classification (if applicable)",
    "Subgroup.number" = "Numeric identifier for subgroup classification",
    "review_doi" = "Digital Object Identifier (DOI) for the Cochrane systematic review",
    "review_url" = "URL to the Cochrane review on the Cochrane Library",
    "Analysis.name" = "Name/description of the meta-analysis comparison",
    "Analysis.group" = "Group/category of the meta-analysis",
    "Analysis.number" = "Numeric identifier for the analysis within the review",
    "Comparison" = "Description of the treatment comparison",
    "Outcome" = "Description of the outcome measure",
    "CI.start" = "Lower bound of 95\\% confidence interval for effect estimate",
    "CI.end" = "Upper bound of 95\\% confidence interval for effect estimate",
    "GIV.Mean" = "Generic inverse variance pooled mean effect estimate",
    "GIV.SE" = "Standard error of the generic inverse variance effect estimate",
    "Mean" = "Pooled mean effect estimate",
    "O.E" = "Observed minus expected events (O-E statistic)",
    "Variance" = "Variance of the effect estimate",
    "Weight" = "Weight assigned to the study in the meta-analysis",
    "Footnotes" = "Additional notes or footnotes from the original Cochrane review",
    "Applicability" = "Assessment of applicability of the study to the review question",
    "Overall.bias..judgement." = "Overall risk of bias judgment (Cochrane Risk of Bias tool)",
    "Overall.bias..support." = "Supporting rationale for overall risk of bias judgment",
    "Bias.arising.from.the.randomization.process..judgement." = "Risk of bias judgment for randomization domain",
    "Bias.arising.from.the.randomization.process..support." = "Supporting rationale for randomization bias judgment",
    "Bias.due.to.deviations.from.intended.interventions..judgement." = "Risk of bias judgment for deviations domain",
    "Bias.due.to.deviations.from.intended.interventions..support." = "Supporting rationale for deviations bias judgment",
    "Bias.due.to.missing.outcome.data..judgement." = "Risk of bias judgment for missing data domain",
    "Bias.due.to.missing.outcome.data..support." = "Supporting rationale for missing data bias judgment",
    "Bias.in.measurement.of.the.outcome..judgement." = "Risk of bias judgment for outcome measurement domain",
    "Bias.in.measurement.of.the.outcome..support." = "Supporting rationale for outcome measurement bias judgment",
    "Bias.in.selection.of.the.reported.result..judgement." = "Risk of bias judgment for selective reporting domain",
    "Bias.in.selection.of.the.reported.result..support." = "Supporting rationale for selective reporting bias judgment"
  )

  # Return description from map, or generic description if not found
  # Escape special Rd characters in the final description
  if (col_name %in% names(desc_map)) {
    return(escape_rd(desc_map[[col_name]]))
  } else {
    return(escape_rd(paste("Additional metadata column from Cochrane review:", col_name)))
  }
}

# Track progress
documented <- 0
failed <- 0

cat("Generating .Rd files...\n")
pb <- txtProgressBar(min = 0, max = length(all_datasets), style = 3)

for (i in seq_along(all_datasets)) {
  ds_name <- all_datasets[i]

  setTxtProgressBar(pb, i)

  tryCatch({
    # Load dataset to extract metadata
    load(file.path("data", paste0(ds_name, ".rda")), envir = environment())
    dataset <- get(ds_name, envir = environment())

    # Get ALL column names from the dataset
    all_columns <- names(dataset)

    # Extract review DOI and URL from first row (all rows have same metadata)
    review_doi <- ""
    review_url <- ""

    if ("review_doi" %in% names(dataset) && nrow(dataset) > 0) {
      review_doi <- as.character(dataset$review_doi[1])
      if (is.na(review_doi)) review_doi <- ""
    }

    if ("review_url" %in% names(dataset) && nrow(dataset) > 0) {
      review_url <- as.character(dataset$review_url[1])
      if (is.na(review_url)) review_url <- ""
    }

    # Extract Cochrane review ID from dataset name
    cochrane_id <- sub("_.*", "", ds_name)

    # Create title
    title <- paste0("Cochrane Review ", cochrane_id, " - Pairwise Meta-Analysis Data")

    # Create description
    description <- paste0(
      "Pairwise meta-analysis dataset extracted from Cochrane Systematic Review ", cochrane_id, ". ",
      "Contains study-level data with treatment comparisons, outcome measures, and effect sizes. ",
      "Data includes both binary outcomes (events/totals) and continuous outcomes (mean/SD/N) where applicable."
    )

    # Add DOI to description if available
    if (nchar(review_doi) > 0) {
      description <- paste0(
        description, "\n\n",
        "\\strong{Cochrane Review DOI:} \\url{https://doi.org/", review_doi, "}"
      )
    }

    # Count studies and determine outcome type
    n_studies <- nrow(dataset)

    has_binary <- any(!is.na(dataset$Experimental.cases)) || any(!is.na(dataset$Control.cases))
    has_continuous <- any(!is.na(dataset$Experimental.mean)) || any(!is.na(dataset$Control.mean))

    outcome_type <- if (has_binary) "binary" else if (has_continuous) "continuous" else "mixed"

    # Add dataset statistics
    description <- paste0(
      description, "\n\n",
      "\\strong{Dataset Statistics:}\n",
      "\\itemize{\n",
      "  \\item Number of studies: ", n_studies, "\n",
      "  \\item Outcome type: ", outcome_type, "\n",
      "  \\item Columns: ", ncol(dataset), "\n",
      "}"
    )

    # Create format description with ALL columns explicitly documented
    format_desc <- paste0(
      "A data frame with ", n_studies, " rows (studies) and ", ncol(dataset), " columns:\n",
      "\\describe{\n"
    )

    # Add explicit \item{} for EVERY column
    for (col_name in all_columns) {
      col_desc <- get_column_description(col_name)
      format_desc <- paste0(format_desc, "  \\item{", col_name, "}{", col_desc, "}\n")
    }

    format_desc <- paste0(format_desc, "}")

    # Create source description
    source_desc <- paste0(
      "Systematically extracted from the Cochrane Library (\\url{https://www.cochranelibrary.com/}) ",
      "using official Cochrane data export functionality. ",
      "Part of the Pairwise70 package containing 501 Cochrane pairwise meta-analysis datasets."
    )

    # Create example code
    example_code <- paste0(
      "# Load the dataset\n",
      "data(", ds_name, ")\n\n",
      "# View dataset structure\n",
      "str(", ds_name, ")\n",
      "head(", ds_name, ")\n\n",
      "# View review metadata\n",
      "cat(\"Cochrane Review DOI:\", ", ds_name, "$review_doi[1], \"\\n\")\n\n"
    )

    # Add example meta-analysis code based on outcome type
    if (has_binary) {
      example_code <- paste0(
        example_code,
        "# Example: Run meta-analysis with metafor (binary outcomes)\n",
        "# library(metafor)\n",
        "# result <- rma(measure = \"OR\",\n",
        "#               ai = Experimental.cases,\n",
        "#               n1i = Experimental.N,\n",
        "#               ci = Control.cases,\n",
        "#               n2i = Control.N,\n",
        "#               data = ", ds_name, ",\n",
        "#               method = \"REML\")\n",
        "# summary(result)\n"
      )
    } else if (has_continuous) {
      example_code <- paste0(
        example_code,
        "# Example: Run meta-analysis with metafor (continuous outcomes)\n",
        "# library(metafor)\n",
        "# result <- rma(measure = \"SMD\",\n",
        "#               m1i = Experimental.mean,\n",
        "#               sd1i = Experimental.SD,\n",
        "#               n1i = Experimental.N,\n",
        "#               m2i = Control.mean,\n",
        "#               sd2i = Control.SD,\n",
        "#               n2i = Control.N,\n",
        "#               data = ", ds_name, ",\n",
        "#               method = \"REML\")\n",
        "# summary(result)\n"
      )
    }

    # Generate .Rd file content
    rd_content <- paste0(
      "\\name{", ds_name, "}\n",
      "\\alias{", ds_name, "}\n",
      "\\docType{data}\n",
      "\\title{", title, "}\n",
      "\\description{\n",
      description, "\n",
      "}\n",
      "\\usage{data(\"", ds_name, "\")}\n",
      "\\format{\n",
      format_desc, "\n",
      "}\n",
      "\\source{\n",
      source_desc, "\n",
      "}\n",
      "\\examples{\n",
      example_code,
      "}\n",
      "\\keyword{datasets}\n"
    )

    # Write .Rd file
    rd_file <- file.path(man_dir, paste0(ds_name, ".Rd"))
    writeLines(rd_content, rd_file)

    documented <- documented + 1

  }, error = function(e) {
    cat("\nError documenting", ds_name, ":", e$message, "\n")
    failed <- failed + 1
  })
}

close(pb)

cat("\n\n=== DOCUMENTATION GENERATION COMPLETE ===\n\n")
cat(paste0("Successfully documented: ", documented, " datasets\n"))
cat(paste0("Failed: ", failed, " datasets\n"))
cat(paste0("\nTotal .Rd files created: ", length(list.files(man_dir, pattern = "\\.Rd$")), "\n"))
cat("\nAll dataset documentation has been generated in the man/ directory!\n")
cat("\nNote: All columns are now explicitly documented to satisfy R CMD check requirements.\n")
