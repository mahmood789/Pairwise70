# Pairwise70 remediation pipeline

This script applies remediation rules across all datasets and runs alternative
meta-analysis methods to handle rare events, sparse data, and missing SDs.

## Requirements
- R (>= 4.0)
- Package: metafor

## Run
From the repo root:

```r
Rscript analysis/meta_remediation.R
```

Optional args:
- First arg: data directory (defaults to `./data`)
- Second arg: output directory (defaults to `./analysis/output`)
- Third arg: GLMM max k (defaults to 30; set to 0 to disable GLMM)

Example:

```r
Rscript analysis/meta_remediation.R "C:/path/to/data" "C:/path/to/output" 30
```

## Remediation rules (summary)
- Binary outcomes: OR and RR with continuity correction (0.5, only zero cells),
  RD without correction, Peto OR for low event rate analyses, and GLMM OR
  for sparse/double-zero analyses (k <= 50).
- Continuous outcomes: SD imputation (within-analysis median or opposite arm),
  then MD or SMD based on analysis name.
- Generic outcomes: use GIV.Mean/SE when available; otherwise derive SE from
  CI or variance based on inferred measure.

## Outputs
- `remediation_analysis_results.csv`: per-analysis results across remedial methods
- `remediation_dataset_summary.csv`: dataset-level remediation counts
- `cleaned_rds/`: datasets with SD-imputed columns and standardized names
- `remediation_report.md`: top-level summary
