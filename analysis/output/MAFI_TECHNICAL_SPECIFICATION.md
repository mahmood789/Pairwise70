# MAFI: Meta-Analysis Fragility Index
## Technical Specification and Validation Report

**Version:** 1.0
**Date:** January 2026
**Authors:** [Your name]

---

## Executive Summary

MAFI (Meta-Analysis Fragility Index) is a novel, comprehensive measure of meta-analysis robustness that addresses key limitations of existing fragility indices. Validated on 4,424 Cochrane meta-analyses, MAFI provides a single interpretable score (0-1) combining direction, significance, and clinical fragility dimensions.

---

## 1. Background and Rationale

### 1.1 Limitations of Existing Fragility Indices

The traditional Fragility Index (Atal et al., 2019) has several limitations:

| Limitation | Traditional FI | MAFI Solution |
|------------|----------------|---------------|
| Only binary outcomes | Yes | Works for any outcome |
| Single dimension (significance) | Yes | 5 dimensions combined |
| Ignores heterogeneity | Yes | Heterogeneity penalty |
| No risk prediction | Yes | Predictive risk scoring |
| No clinical threshold | Often | Customizable thresholds |
| Absolute count only | Yes | Rates and composite score |

### 1.2 Evidence Base

MAFI is derived from empirical analysis of:
- **4,424 pairwise meta-analyses**
- **501 Cochrane systematic reviews**
- **86,492 study-level observations**
- **38.3 million participants**

---

## 2. MAFI Components

### 2.1 Core Components (5 Dimensions)

| Component | Symbol | Definition | Weight |
|-----------|--------|------------|--------|
| Direction Fragility Index | DFI | Studies that flip effect direction | 30% |
| Significance Fragility Index | SFI | Studies that change p < 0.05 status | 25% |
| Clinical Fragility Index | CFI | Studies that cross clinical threshold | 20% |
| Effect Stability Index | ESI | Maximum proportional effect change | 15% |
| CI Stability Index | CISI | Studies that flip CI inclusion of null | 10% |

### 2.2 Penalty Terms

| Penalty | Formula | Maximum |
|---------|---------|---------|
| Heterogeneity | (I² / 100) × 0.2 | 20% |
| Sample Size | max(0, (1 - k/20) × 0.3) | 30% |

### 2.3 Final MAFI Calculation

```
MAFI_core = 0.30×DFI_rate + 0.25×SFI_rate + 0.20×CFI_rate +
            0.15×min(ESI,1) + 0.10×CISI_rate

MAFI = min(1, MAFI_core + het_penalty + k_penalty)
```

---

## 3. MAFI Classification

| MAFI Score | Class | Interpretation |
|------------|-------|----------------|
| 0.00 - 0.15 | Robust | Conclusions stable; high confidence |
| 0.15 - 0.30 | Low Fragility | Generally reliable; note sensitivity |
| 0.30 - 0.50 | Moderate Fragility | Interpret with caution |
| 0.50 - 1.00 | High Fragility | Strong caveats required |

---

## 4. Validation Results

### 4.1 Distribution on Cochrane Data (n=4,424)

| Class | N | Percentage |
|-------|---|------------|
| Robust | 1,186 | 26.8% |
| Low Fragility | 1,591 | 36.0% |
| Moderate Fragility | 1,109 | 25.1% |
| High Fragility | 430 | 9.7% |

**Summary Statistics:**
- Mean MAFI: 0.260 (SD: 0.175)
- Median MAFI: 0.242
- Range: 0.000 - 1.000

### 4.2 Correlation with Existing Measures

| Measure | Correlation with MAFI |
|---------|----------------------|
| Composite Fragility (0-3) | r = 0.610 |
| Fragility Quotient | r = 0.673 |
| Direction Fragility Index | r = 0.122 |
| Significance Fragility Index | r = 0.119 |

### 4.3 Predictive Model Performance

| Model | AUC | Pseudo R² |
|-------|-----|-----------|
| Direction Fragility | 0.837 | 0.267 |
| Significance Fragility | 0.773 | 0.141 |

---

## 5. Minimum k Recommendations

Based on empirical thresholds from Cochrane analysis:

| Target Fragility Level | Direction < | Significance < | Minimum k |
|------------------------|-------------|----------------|-----------|
| Very Low | 10% | 10% | 35 |
| Low | 20% | 15% | 20 |
| Moderate | 30% | 20% | 10 |

**Key Finding:** Meta-analyses with k ≥ 20 studies achieve 55.7% robustness rate, compared to 34.5% for k = 3-4.

---

## 6. Risk Prediction Equations

### 6.1 Direction Fragility Risk
```
logit(p) = -0.86 - 0.64×log(k) + 0.02×I² + 0.21×τ² - 5.62×|effect|
```

### 6.2 Significance Fragility Risk
```
logit(p) = -0.87 - 0.56×log(k) + 0.02×I² - 0.01×τ² + 1.37×significant
```

### 6.3 Combined Risk Score
```
Risk = 0.5 × p_direction + 0.5 × p_significance
```

---

## 7. Clinical Thresholds by Outcome Type

| Measure | Default Threshold | Clinical Interpretation |
|---------|-------------------|------------------------|
| Odds Ratio (OR) | 1.25 | Minimum clinically important |
| Risk Ratio (RR) | 1.25 | Minimum clinically important |
| Standardized Mean Difference (SMD) | 0.2 | Small effect (Cohen) |
| Mean Difference (MD) | User-defined | Domain-specific |

---

## 8. Comparison with Atal et al. (2019)

| Aspect | Atal et al. (2019) | MAFI |
|--------|-------------------|------|
| Sample | 906 Cochrane MAs | 4,424 Cochrane MAs |
| Scope | Binary outcomes only | Any outcome type |
| Focus | Event modifications | Leave-one-out analysis |
| Output | Single count | Composite 0-1 score |
| Classification | FI ≤ 5 = fragile | 4-level classification |
| Risk Prediction | No | Yes (AUC 0.77-0.84) |
| Median (significant) | 12 | Varies by component |

---

## 9. Implementation

### 9.1 R Function

```r
# Calculate MAFI for a meta-analysis
mafi_result <- calculate_MAFI(
  yi = effect_sizes,
  vi = variances,
  measure = "SMD",
  clinical_threshold = 0.2
)

# Access results
mafi_result$MAFI           # 0-1 score
mafi_result$MAFI_class     # Classification
mafi_result$interpretation # Human-readable summary
mafi_result$risk_score     # Predictive risk
```

### 9.2 Required Inputs

- `yi`: Vector of effect sizes
- `vi`: Vector of variances (or `sei` for standard errors)
- `measure`: "OR", "RR", "SMD", or "MD"
- `clinical_threshold`: Optional custom threshold

### 9.3 Outputs

- MAFI score (0-1)
- MAFI class (4 levels)
- Component indices (DFI, SFI, CFI, ESI, CISI)
- Predictive risk scores
- Human-readable interpretation

---

## 10. Recommendations for Reporting

### 10.1 Minimum Reporting

1. Report MAFI score and classification
2. Report number of fragility-inducing studies
3. State if k < 10 (recommend caveat)

### 10.2 Extended Reporting

1. All component indices (DFI, SFI, CFI, ESI, CISI)
2. Risk prediction scores
3. Heterogeneity penalty applied
4. Leave-one-out forest plot

### 10.3 GRADE Integration

| MAFI Class | GRADE Downgrade Suggestion |
|------------|----------------------------|
| Robust | None for fragility |
| Low Fragility | Consider -1 if other concerns |
| Moderate Fragility | -1 for imprecision/inconsistency |
| High Fragility | -2 for serious concerns |

---

## 11. Limitations

1. Leave-one-out analysis computationally intensive for large k
2. Risk prediction equations derived from Cochrane data (may differ for non-Cochrane)
3. Clinical thresholds require domain expertise for MD outcomes
4. Does not replace publication bias assessment

---

## 12. Future Directions

1. Web-based MAFI calculator
2. Network meta-analysis extension
3. Bayesian MAFI implementation
4. Living review integration
5. Machine learning enhancement of risk prediction

---

## 13. References

1. Atal I, Porcher R, Boutron I, Ravaud P. The statistical significance of meta-analyses is frequently fragile: definition of a Fragility Index for meta-analyses. J Clin Epidemiol. 2019;111:32-40.

2. Walsh M, Srinathan SK, McAuley DF, et al. The statistical significance of randomized controlled trial results is frequently fragile: a case for a Fragility Index. J Clin Epidemiol. 2014;67(6):622-628.

3. Ahmed W, Fowler RA, McCredie VA. Does sample size matter when interpreting the Fragility Index? Crit Care Med. 2016;44(11):e1142-e1143.

4. Caldwell DM, Ades AE, Higgins JPT. Simultaneous comparison of multiple treatments: combining direct and indirect evidence. BMJ. 2005;331(7521):897-900.

5. Cochrane Handbook for Systematic Reviews of Interventions. Version 6.4, 2023.

---

## Appendix A: MAFI Distribution Visualization

```
MAFI Score Distribution (n=4,424 Cochrane Meta-Analyses)

Robust (0-0.15):        ████████████████ 26.8%
Low Fragility (0.15-0.30):  ██████████████████████ 36.0%
Moderate (0.30-0.50):   ████████████████ 25.1%
High Fragility (0.50+): ██████ 9.7%
```

## Appendix B: Files Generated

1. `MAFI_development.R` - Core implementation
2. `MAFI_validated_results.csv` - Validated scores for 4,424 MAs
3. `MAFI_TECHNICAL_SPECIFICATION.md` - This document

---

*MAFI was developed using the Pairwise70 dataset of 501 Cochrane systematic reviews.*
