# Pairwise70 meta-analysis diagnostics

Generated: 2026-01-01 13:49:10

## Coverage
- Datasets: 501
- Total analyses: 8087
- Total rows: 86492
- Meta-analyses succeeded: 6240
- Meta-analyses failed or insufficient data: 1847

## Common data issues (row-level counts)
- Double-zero event rows: 37727
- All-event rows: 944
- Sparse-event rows (<5 total events): 48563
- Cases > N (experimental): 1
- Cases > N (control): 1
- SD <= 0 (experimental): 52120
- SD <= 0 (control): 52077
- Variance < 0: 0
- SE < 0: 0
- CI inverted: 0
- CI mismatch vs Mean/SE: 64628

## Meta-analysis stability signals
- High heterogeneity (I2 >= 75%): 495
- Egger test p < 0.05 (k>=10): 423
- Large FE vs RE shift (|diff| > 0.2): 463
- Max weight share > 0.5: 1515

## Methodological remedies to consider
- Zero-event or all-event studies: continuity corrections (0.5, treatment-arm, or empirical), Peto OR for rare events, or GLMM/beta-binomial models. Double-zero studies may be excluded for OR or analyzed with risk difference.
- Sparse events: use exact methods, GLMM, or switch to risk ratio/risk difference; report sensitivity to different corrections.
- SD missing/zero: verify extraction, impute SD from pooled SD or similar studies, or use mean difference when scales align.
- Cases > N or negative values: fix data integrity before analysis; these rows should be corrected or excluded.
- CI/SE/variance inconsistencies: re-derive SE from CI or variance and verify scale (log vs raw).
- Large FE vs RE shifts or high influence: report sensitivity, consider robust variance or leave-one-out checks.

## Output files
- dataset_inventory.csv: dataset-level inventory and completeness
- dataset_issue_counts.csv: dataset-level issue counts
- analysis_diagnostics_results.csv: per-analysis meta-results and diagnostics
- column_frequency.csv: column prevalence across datasets
- issue_totals.csv: aggregate issue counts
- diagnosis_report.md: this report
