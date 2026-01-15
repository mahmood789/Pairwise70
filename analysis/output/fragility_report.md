# Fragility of Cochrane Meta-Analyses: An Empirical Investigation

Generated: %Y-%m-%d %H:%M:%S

## Abstract

This analysis examines the fragility of 4424 Cochrane pairwise meta-analyses from 501 systematic reviews. We define fragility as the vulnerability of meta-analysis conclusions to the removal of a single study.

## Key Findings

### Sample Characteristics

- **Total meta-analyses**: 4424
- **Median studies per analysis**: 9 (range: 3 - 494)
- **Statistically significant (p<0.05)**: 1072 (24.2%)

### Prevalence of Fragility

| Fragility Type | N | Percentage |
|----------------|---|------------|
| Direction fragile (effect sign flips) | 1154 | 26.1% |
| Significance fragile (p-value crosses 0.05) | 648 | 14.6% |
| Clinical fragile (crosses threshold) | 1275 | 28.8% |

### Composite Fragility Classification

| Class | N | Percentage |
|-------|---|------------|
| Robust (no fragility) | 1947 | 44.0% |
| Low Fragility (1 type) | 1701 | 38.4% |
| Moderate Fragility (2 types) | 628 | 14.2% |
| High Fragility (3 types) | 40 | 0.9% |

### Significance Fragility Among Significant Results

Of 1072 meta-analyses with p < 0.05:
- **321 (29.9%)** would lose statistical significance if a single study were removed

This represents a substantial proportion of Cochrane meta-analyses whose conclusions depend critically on individual studies.

### Fragility by Number of Studies

| Studies | N | Direction Fragile | Significance Fragile | Mean Composite |
|---------|---|-------------------|---------------------|----------------|
| 3-5 | 1410 | 34.4% | 17.2% | 0.92 |
| 6-10 | 1152 | 28.9% | 15.9% | 0.77 |
| 11-20 |  938 | 25.5% | 13.9% | 0.64 |
| 21-50 |  637 | 16.2% | 13.7% | 0.47 |
| 51+ |  287 | 6.4% | 6.7% | 0.21 |

### Fragility by Heterogeneity

| I-squared | N | Direction Fragile | Significance Fragile | Mean Composite |
|-----------|---|-------------------|---------------------|----------------|
| Low (<25%) | 3261 | 28.2% | 9.2% | 0.65 |
| Moderate (25-50%) |  436 | 22.5% | 31.4% | 0.91 |
| Substantial (50-75%) |  401 | 25.3% | 30.4% | 0.90 |
| High (>75%) |  326 | 20.1% | 33.0% | 0.84 |

## Implications

1. **For systematic reviewers**: Report sensitivity of conclusions to individual study removal
2. **For guideline developers**: Consider fragility when translating evidence to recommendations
3. **For methodologists**: Develop robust pooling methods that minimize fragility
4. **For journal editors**: Consider requiring fragility assessment as standard reporting

## Methods

- **Data source**: 501 Cochrane systematic reviews (Pairwise70 dataset)
- **Analyses included**: Meta-analyses with k >= 3 studies
- **Effect measures**: OR for binary, SMD for continuous outcomes
- **Fragility metrics**:
  - Direction fragility: Sign of effect estimate reverses
  - Significance fragility: p-value crosses 0.05 threshold
  - Clinical fragility: Effect crosses clinical threshold (OR=1.25, SMD=0.2)
  - Fragility quotient: Proportion of studies whose removal causes any fragility
  - Composite fragility: Sum of fragility types (0-3 scale)

## Visualizations

See `plots/fragility/` for:
1. Fragility class distribution
2. Fragility by number of studies (boxplot)
3. Fragility quotient histogram
4. Significance fragility breakdown
5. Fragility vs heterogeneity
6. Fragility rates by study count (bar chart)

## Data Availability

Full results: `fragility_analysis_results.csv` (4424 rows)
Summary: `fragility_summary.csv`


