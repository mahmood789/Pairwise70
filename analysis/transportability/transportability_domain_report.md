# Transportability Domain Report

## Domain coverage
- cardiovascular: 41 reviews (1074 rows, 52 reviews_in_list)
- cerebrovascular: 16 reviews (608 rows, 16 reviews_in_list)
- dentistry: 10 reviews (340 rows, 11 reviews_in_list)
- dermatology: 10 reviews (524 rows, 11 reviews_in_list)
- endocrine_metabolic: 8 reviews (293 rows, 10 reviews_in_list)
- gastrointestinal: 14 reviews (402 rows, 18 reviews_in_list)
- infectious: 41 reviews (806 rows, 48 reviews_in_list)
- maternal_neonatal: 39 reviews (1802 rows, 40 reviews_in_list)
- mental_health: 37 reviews (1106 rows, 44 reviews_in_list)
- musculoskeletal: 17 reviews (805 rows, 18 reviews_in_list)
- neurology: 9 reviews (210 rows, 9 reviews_in_list)
- oncology: 22 reviews (470 rows, 23 reviews_in_list)
- ophthalmology: 8 reviews (69 rows, 9 reviews_in_list)
- other: 128 reviews (1638 rows, 128 reviews_in_list)
- pain: 57 reviews (1186 rows, 70 reviews_in_list)
- renal: 18 reviews (395 rows, 28 reviews_in_list)
- respiratory: 18 reviews (543 rows, 20 reviews_in_list)
- urology: 8 reviews (394 rows, 10 reviews_in_list)

## Mortality tag prevalence: 0.40
## Multi-label prevalence: 0.10
## Secondary domain prevalence: 0.10

## Interpretation notes
- Domains are assigned by keyword voting on Analysis.name; the top-scoring domain is primary.
- Multi-label reviews carry a domain_list and optional secondary domain.
- Mortality is treated as an outcome tag, not a domain.
- Low-count domains are reported but excluded from CV if < 30 rows.

## Recommended target-population covariates (open registries)
- Study timing (start/completion year) to align clinical era.
- Enrollment size distribution to approximate trial size/setting.
- Sex eligibility proportions as a basic inclusion proxy.
- For cardio/cerebrovascular: include hypertension or vascular terms in query refinement.
- For infectious: include pathogen-specific keywords (e.g., HIV, hepatitis).
- For mental health: include disorder-specific terms (e.g., depression, anxiety).
