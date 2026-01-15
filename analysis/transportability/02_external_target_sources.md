# External Target-Population Data Sources (Transportability)

This list prioritizes sources with either open access or a clear request mechanism. Access terms and exact fields vary by source.

| Source | Type | Access | Potential use | Link |
| --- | --- | --- | --- | --- |
| WHO ICTRP "WHO data set" | Trial registry metadata | Open download | Trial population descriptors, eligibility, setting/region proxies | https://www.who.int/tools/clinical-trials-registry-platform/network/who-data-set |
| ClinicalTrials.gov Data API | Trial registry metadata/results | Open API | Population descriptors, intervention taxonomy, baseline risk proxies | https://clinicaltrials.gov/data-api |
| OpenTrials | Linked trial registry/meta-docs | Open | Cross-link trial identifiers, dates, conditions, settings | https://opentrials.net/about/index.html |
| Vivli | IPD repository | Request/gated | IPD for target-population covariates and baseline risk | https://vivli.org/ |
| YODA Project | IPD repository | Request/gated | IPD in selected domains | https://yoda.yale.edu/ |
| CSDR | IPD repository | Request/gated | IPD in selected domains | https://clinicalstudydatarequest.com/ |
| MIMIC-IV | EHR (ICU) | Credentialed (PhysioNet) | External population baseline risk, covariate distributions | https://physionet.org/content/mimiciv/2.2/ |
| NHANES | Population survey | Open download | Population covariates and baseline risk proxies | https://wwwn.cdc.gov/nchs/nhanes/Default.aspx |
| Synthea | Synthetic EHR | Open | Sensitivity testing for transportability assumptions | https://synthea.mitre.org/ |
| EpistasisLab ClinicalDataSources | Catalog | Open | Find additional open-access target datasets | https://github.com/EpistasisLab/ClinicalDataSources |

## Selection logic
- Prioritize open registries and surveys for baseline population covariates.
- Add IPD repositories where domain-specific covariates are essential.
- Use synthetic populations to test robustness when real target data are unavailable.

## Next actions
- Choose 1-2 external sources aligned with the clinical domain.
- Map baseline risk and covariates to Pairwise70 fields (Study.year, Subgroup, Analysis.name, review_doi, etc.).
- Record access steps and data provenance for reproducibility.
