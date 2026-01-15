# Pairwise70 meta-analysis diagnostics

This script inventories and diagnoses all 501 Cochrane pairwise datasets and produces
dataset-level and per-analysis diagnostics for common meta-analysis issues.

## Requirements
- R (>= 4.0)
- Package: metafor

## Run
From the repo root:

```r
Rscript analysis/meta_diagnostics.R
```

Optional args:
- First arg: data directory (defaults to `./data`)
- Second arg: output directory (defaults to `./analysis/output`)

Example:

```r
Rscript analysis/meta_diagnostics.R "C:/path/to/data" "C:/path/to/output"
```

## Outputs (written to the output directory)
- `dataset_inventory.csv`: dataset-level inventory and completeness
- `dataset_issue_counts.csv`: dataset-level issue counts
- `analysis_diagnostics_results.csv`: per-analysis meta-results and diagnostics
- `column_frequency.csv`: column prevalence across datasets
- `issue_totals.csv`: aggregate issue counts
- `diagnosis_report.md`: human-readable report with suggested remedies
