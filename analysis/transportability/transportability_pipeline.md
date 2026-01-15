# Transportability Pipeline (Open Registries)

This pipeline uses Pairwise70 study data and open trial registries to approximate target-population covariate distributions.

## Outputs
- Pairwise audit: `analysis/transportability/pairwise70_audit.csv`
- Review-level summaries: `analysis/transportability/transportability_review_level.csv`
- Query terms for registries: `analysis/transportability/ctgov_query_terms.csv`
- Registry covariates: `analysis/transportability/ctgov_target_covariates.csv`
- Final merged dataset: `analysis/transportability/transportability_target_merge.csv`

## Steps
1) Audit and review-level summaries
   - `01_audit_pairwise70.R`
   - `03_transportability_scaffold.R`
   - Optional rebuild: set `REBUILD_REVIEW=1` to regenerate review-level summaries.

2) Build registry query terms
   - `05_build_ctgov_queries.R`
   - Edit `ctgov_query_terms.csv` if you want to refine terms or add domain-specific keywords.

3) Fetch ClinicalTrials.gov registry data (open)
   - `06_ctgov_fetch.R`
   - Writes raw JSON to `analysis/transportability/external/ctgov/raw/`
   - Optional limits (environment variables):
     - `CTGOV_START_INDEX` (default: 1)
     - `CTGOV_MAX_REVIEWS` (default: all)
     - `CTGOV_MAX_PAGES` (default: 3)
     - `CTGOV_PAGE_SIZE` (default: 100)
     - `CTGOV_SLEEP_MS` (default: 0)

4) Aggregate registry covariates
   - `07_ctgov_aggregate.R`

5) Merge with review-level effects
   - `08_transportability_merge.R`

6) Transportability CV using registry covariates
   - `09_transportability_target_cv.R`
   - CV uses leave-one-review-out folds and `rma.mv` with a review-level random intercept.

7) Domain-stratified transportability CV
   - `10_transportability_domain_cv.R`

8) Domain report
   - `11_transportability_domain_report.R`

9) Domain-specific registry queries
   - `12_build_domain_queries.R`

10) Domain-specific registry fetch and aggregate
   - `13_ctgov_domain_fetch.R`
   - `14_ctgov_domain_aggregate.R`
   - Optional limits (environment variables):
     - `CTGOV_DOMAIN_USE_SHORT` (default: 1)
     - `CTGOV_DOMAIN_FORCE` (default: 0)
     - `CTGOV_DOMAIN_MAX_PAGES` (default: 5)
     - `CTGOV_DOMAIN_PAGE_SIZE` (default: 200)
     - `CTGOV_DOMAIN_SLEEP_MS` (default: 200)
     - `CTGOV_DOMAIN_LIST` (comma-separated domains)

## Notes
- ClinicalTrials.gov API rate limits apply; `06_ctgov_fetch.R` uses pagination and caches results per review.
- If you want to target a specific clinical domain, update `ctgov_query_terms.csv` to keep queries precise.
- This pipeline does not compute baseline risk from registry data; use external outcome registries or EHR datasets if needed.
