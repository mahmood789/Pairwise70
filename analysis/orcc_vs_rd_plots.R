#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[1] else file.path("analysis", "output")
plots_dir <- file.path(output_dir, "plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

disagree_path <- file.path(output_dir, "orcc_vs_rd_disagreements.csv")
if (!file.exists(disagree_path)) {
  stop(paste0("Missing file: ", disagree_path))
}

df <- read.csv(disagree_path, stringsAsFactors = FALSE)
df <- df[!is.na(df$log_or) & !is.na(df$rd), ]

log_or <- as.numeric(df$log_or)
rd <- as.numeric(df$rd)

png(file.path(plots_dir, "orcc_vs_rd_scatter.png"), width = 900, height = 700)
plot(
  log_or, rd,
  pch = 16,
  col = rgb(0, 0, 0, 0.2),
  xlab = "Log OR (OR_cc)",
  ylab = "Risk difference (RD)",
  main = "OR_cc vs RD (paired analyses)"
)
abline(h = 0, v = 0, col = "gray70", lwd = 1)
dev.off()

mismatch <- df$mismatch_type %in% c("direction_only", "significance_only", "direction_and_significance")
png(file.path(plots_dir, "orcc_vs_rd_scatter_mismatch.png"), width = 900, height = 700)
plot(
  log_or, rd,
  pch = 16,
  col = rgb(0.2, 0.2, 0.2, 0.2),
  xlab = "Log OR (OR_cc)",
  ylab = "Risk difference (RD)",
  main = "OR_cc vs RD (mismatches highlighted)"
)
points(log_or[mismatch], rd[mismatch], pch = 16, col = rgb(0.9, 0.2, 0.2, 0.6))
abline(h = 0, v = 0, col = "gray70", lwd = 1)
legend("topright", legend = c("All", "Mismatch"), pch = 16,
       col = c(rgb(0.2, 0.2, 0.2, 0.3), rgb(0.9, 0.2, 0.2, 0.6)))
dev.off()

sparse <- df$sparse_events
double_zero <- df$double_zero
group <- ifelse(sparse & double_zero, "both",
  ifelse(sparse, "sparse_only",
    ifelse(double_zero, "double_zero_only", "none")
  )
)
cols <- c(
  none = rgb(0.2, 0.2, 0.2, 0.2),
  sparse_only = rgb(0.2, 0.6, 0.9, 0.5),
  double_zero_only = rgb(0.6, 0.4, 0.9, 0.5),
  both = rgb(0.9, 0.4, 0.2, 0.6)
)

png(file.path(plots_dir, "orcc_vs_rd_scatter_sparse_doublezero.png"), width = 900, height = 700)
plot(
  log_or, rd,
  pch = 16,
  col = cols[group],
  xlab = "Log OR (OR_cc)",
  ylab = "Risk difference (RD)",
  main = "OR_cc vs RD (sparse/double-zero)"
)
abline(h = 0, v = 0, col = "gray70", lwd = 1)
legend(
  "topright",
  legend = c("none", "sparse_only", "double_zero_only", "both"),
  pch = 16,
  col = cols[c("none", "sparse_only", "double_zero_only", "both")]
)
dev.off()

k_vals <- as.numeric(df$k)
k_bin <- cut(
  k_vals,
  breaks = c(-Inf, 4, 9, 19, Inf),
  labels = c("2-4", "5-9", "10-19", "20+"),
  right = TRUE
)
mismatch_rate <- tapply(mismatch, k_bin, function(x) mean(x, na.rm = TRUE))

png(file.path(plots_dir, "orcc_vs_rd_mismatch_rate_k.png"), width = 700, height = 500)
barplot(
  mismatch_rate,
  ylim = c(0, max(mismatch_rate, na.rm = TRUE) * 1.2),
  col = "steelblue",
  ylab = "Mismatch rate",
  xlab = "k bin",
  main = "Mismatch rate by k bin"
)
dev.off()

mismatch_rate_sparse <- c(
  sparse = mean(mismatch[sparse], na.rm = TRUE),
  non_sparse = mean(mismatch[!sparse], na.rm = TRUE),
  double_zero = mean(mismatch[double_zero], na.rm = TRUE),
  no_double_zero = mean(mismatch[!double_zero], na.rm = TRUE)
)

png(file.path(plots_dir, "orcc_vs_rd_mismatch_rate_sparsity.png"), width = 800, height = 500)
barplot(
  mismatch_rate_sparse,
  col = "darkorange",
  ylab = "Mismatch rate",
  main = "Mismatch rate by sparsity/double-zero"
)
dev.off()

cat("Plots written to: ", plots_dir, "\n", sep = "")
