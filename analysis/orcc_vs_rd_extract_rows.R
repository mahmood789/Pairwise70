#!/usr/bin/env Rscript

standardize_names <- function(nms) {
  nms <- gsub("[^A-Za-z0-9]+", ".", nms)
  nms <- gsub("\\.+", ".", nms)
  nms <- gsub("^\\.|\\.$", "", nms)
  nms
}

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[1] else file.path("analysis", "output")
data_dir <- if (length(args) >= 2) args[2] else file.path("data")
top_n <- if (length(args) >= 3) suppressWarnings(as.numeric(args[3])) else NA_real_
min_paired <- if (length(args) >= 4) suppressWarnings(as.numeric(args[4])) else 1
if (is.na(min_paired) || min_paired < 1) {
  min_paired <- 1
}

disagree_path <- file.path(output_dir, "orcc_vs_rd_disagreements.csv")
dataset_path <- file.path(output_dir, "orcc_vs_rd_dataset_mismatch_counts.csv")
out_dir <- file.path(output_dir, "orcc_vs_rd_mismatch_rows")
manifest_path <- file.path(output_dir, "orcc_vs_rd_mismatch_manifest.csv")

if (!file.exists(disagree_path)) {
  stop(paste0("Missing file: ", disagree_path))
}
if (!file.exists(dataset_path)) {
  stop(paste0("Missing file: ", dataset_path))
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

disagree <- read.csv(disagree_path, stringsAsFactors = FALSE)
dataset_summary <- read.csv(dataset_path, stringsAsFactors = FALSE)

dataset_summary <- dataset_summary[dataset_summary$paired_analyses >= min_paired, ]
dataset_summary <- dataset_summary[order(-dataset_summary$mismatch_rate, -dataset_summary$paired_analyses), ]
if (is.na(top_n) || top_n <= 0 || top_n >= nrow(dataset_summary)) {
  target_datasets <- dataset_summary$dataset_name
} else {
  target_datasets <- head(dataset_summary$dataset_name, top_n)
}

mismatch_types <- c("direction_only", "significance_only", "direction_and_significance")
disagree <- disagree[disagree$dataset_name %in% target_datasets &
  disagree$mismatch_type %in% mismatch_types, ]

if (nrow(disagree) == 0) {
  stop("No mismatches found for the selected datasets.")
}

manifest <- disagree
write.csv(manifest, manifest_path, row.names = FALSE)

for (ds in unique(disagree$dataset_name)) {
  ds_rows <- disagree[disagree$dataset_name == ds, ]
  ds_rows <- ds_rows[!duplicated(ds_rows$analysis_key), ]

  rda_path <- file.path(data_dir, paste0(ds, ".rda"))
  if (!file.exists(rda_path)) {
    warning(paste0("Missing dataset file: ", rda_path))
    next
  }

  env <- new.env()
  load(rda_path, envir = env)
  obj_names <- ls(env)
  if (length(obj_names) == 0) {
    next
  }
  df <- env[[obj_names[1]]]
  if (!is.data.frame(df)) {
    next
  }

  names(df) <- standardize_names(names(df))
  n_rows <- nrow(df)

  analysis_number <- if ("Analysis.number" %in% names(df)) as.character(df$Analysis.number) else rep("1", n_rows)
  subgroup_number <- if ("Subgroup.number" %in% names(df)) as.character(df$Subgroup.number) else rep(NA, n_rows)
  subgroup_id <- ifelse(is.na(subgroup_number) | subgroup_number == "", "overall", subgroup_number)
  analysis_key <- paste(analysis_number, subgroup_id, sep = "::")
  df$analysis_key <- analysis_key

  idx <- match(df$analysis_key, ds_rows$analysis_key)
  keep <- !is.na(idx)
  df_out <- df[keep, , drop = FALSE]
  map <- ds_rows[idx[keep], ]

  df_out$dataset_name <- ds
  df_out$mismatch_type <- map$mismatch_type
  df_out$direction_match <- map$direction_match
  df_out$or_sig <- map$or_sig
  df_out$rd_sig <- map$rd_sig
  df_out$log_or <- map$log_or
  df_out$rd <- map$rd
  df_out$log_or_ci_lb <- map$log_or_ci_lb
  df_out$log_or_ci_ub <- map$log_or_ci_ub
  df_out$rd_ci_lb <- map$rd_ci_lb
  df_out$rd_ci_ub <- map$rd_ci_ub
  df_out$sparse_events <- map$sparse_events
  df_out$double_zero <- map$double_zero

  out_path <- file.path(out_dir, paste0(ds, "_mismatch_rows.csv"))
  write.csv(df_out, out_path, row.names = FALSE)
}

cat("Mismatch row extraction complete.\n")
cat("Manifest: ", manifest_path, "\n", sep = "")
cat("Output folder: ", out_dir, "\n", sep = "")
