# MA4 v1.0.1 Comprehensive Statistical Analysis Report
## Pairwise70 Cochrane Meta-Analysis Collection

**Date:** January 2026
**Dataset:** 5,088 meta-analyses from 486 Cochrane systematic reviews

---

## Executive Summary

This report presents a comprehensive statistical analysis of meta-analytic stability using the MA4 framework (theta, sigma, tau, R) applied to 5,088 pairwise meta-analyses extracted from the Pairwise70 Cochrane collection.

### Key Findings

1. **Overall Stability**: Mean R = 0.674, Median R = 0.673
2. **High Stability Rate**: 31.4% of meta-analyses achieve R >= 0.8
3. **Sign Instability**: 19.9% have R < 0.5 (potential sign flip risk)
4. **Validation**: MA4 achieves **perfect correlation (r = 1.000)** with metafor REML estimates
5. **Best Predictive Model**: Explains 20.5% of variance in R (R-squared = 0.205)

---

## 1. Data Overview

| Metric | Value |
|--------|-------|
| Total meta-analyses | 5,088 |
| Unique Cochrane reviews | 486 |
| Effect types | logRR (binary), GIV (pre-computed), MD (continuous) |
| Number of studies (k) range | 2 - 167 |
| Median k | 5 |
| Mean k | 8.4 |

### Effect Type Distribution

| Effect Type | Count | Percentage |
|-------------|-------|------------|
| logRR (binary outcomes) | 3,556 | 69.9% |
| GIV (generic inverse variance) | 1,350 | 26.5% |
| MD (mean difference) | 182 | 3.6% |

---

## 2. R Stability Distribution

### Summary Statistics

| Statistic | Value |
|-----------|-------|
| Minimum | 0.022 |
| Q1 (25th percentile) | 0.500 |
| Median | 0.673 |
| Mean | 0.674 |
| Q3 (75th percentile) | 0.880 |
| Maximum | 1.000 |
| Standard Deviation | 0.224 |

### R Category Distribution

| Category | Count | Percentage |
|----------|-------|------------|
| Very Low (R < 0.3) | 140 | 2.8% |
| Low (0.3 - 0.5) | 1,803 | 35.4% |
| Moderate (0.5 - 0.7) | 819 | 16.1% |
| Good (0.7 - 0.8) | 730 | 14.3% |
| High (0.8 - 0.9) | 390 | 7.7% |
| Excellent (0.9 - 1.0) | 1,206 | 23.7% |

---

## 3. Predictors of Stability

### 3.1 Correlation Analysis

| Variable | Pearson r | Spearman rho | Interpretation |
|----------|-----------|--------------|----------------|
| log(k) | -0.162 | -0.214 | More studies = slightly lower R |
| I2 proxy | -0.169 | -0.286 | More heterogeneity = lower R |
| log(tau) | -0.166 | 0.000 | Complex relationship |
| abs(theta) | 0.039 | 0.296 | Larger effects = higher R |
| log(sigma) | 0.246 | 0.253 | Larger SE = higher R |

### 3.2 Best Regression Model

**Model:** R ~ log(k) + effect_type + log(tau) + |theta| + theta_near_zero + log(k):effect_type

| Predictor | Coefficient | Std. Error | t-value | p-value |
|-----------|------------|------------|---------|---------|
| (Intercept) | 0.508 | 0.014 | 36.97 | < 0.001 |
| log(k) | 0.014 | 0.007 | 2.10 | 0.036 |
| effect_type = logRR | 0.251 | 0.014 | 17.59 | < 0.001 |
| effect_type = MD | 0.225 | 0.035 | 6.43 | < 0.001 |
| log(tau) | -0.011 | 0.001 | -16.89 | < 0.001 |
| abs(theta) | 0.002 | 0.001 | 3.13 | 0.002 |
| theta_near_zero | -0.091 | 0.007 | -13.97 | < 0.001 |
| log(k):logRR | -0.068 | 0.008 | -8.97 | < 0.001 |
| log(k):MD | -0.097 | 0.020 | -4.84 | < 0.001 |

**Model Performance:**
- R-squared: 0.205
- Adjusted R-squared: 0.203
- Residual SE: 0.200
- AIC: -1917.4 (best among all models tested)

---

## 4. Subgroup Analyses

### 4.1 R by Number of Studies (k)

| k Category | n | Mean R | Median R | SD |
|------------|---|--------|----------|-----|
| k = 2 | 1,124 | 0.823 | 1.000 | 0.229 |
| k = 3 | 725 | 0.602 | 0.500 | 0.217 |
| k = 4-5 | 943 | 0.617 | 0.624 | 0.207 |
| k = 6-10 | 1,095 | 0.636 | 0.660 | 0.200 |
| k = 11-20 | 824 | 0.647 | 0.682 | 0.193 |
| k = 21-50 | 317 | 0.684 | 0.727 | 0.193 |
| k > 50 | 60 | 0.634 | 0.618 | 0.186 |

**Key Insight**: k=2 shows highest R (often R=1) because there's no perturbation diversity; k=3-5 shows lowest R due to high sensitivity to single study removal.

### 4.2 R by Effect Type

| Effect Type | n | Mean R | Median R | SD |
|-------------|---|--------|----------|-----|
| logRR | 3,556 | 0.718 | 0.740 | 0.234 |
| MD | 182 | 0.648 | 0.682 | 0.249 |
| GIV | 1,350 | 0.562 | 0.500 | 0.140 |

**Key Insight**: Binary outcome meta-analyses (logRR) show highest stability; GIV analyses show lowest due to theta often being exactly 0.

### 4.3 R by Heterogeneity Level

| Heterogeneity | n | Mean R | Median R |
|---------------|---|--------|----------|
| Low (tau < 0.1) | 3,690 | 0.696 | 0.699 |
| Moderate (0.1-0.3) | 518 | 0.625 | 0.661 |
| High (0.3-0.5) | 335 | 0.615 | 0.657 |
| Very High (tau > 0.5) | 545 | 0.605 | 0.655 |

### 4.4 R by Effect Size Magnitude

| Effect Magnitude | n | Mean R | Median R |
|------------------|---|--------|----------|
| Negligible (|theta| < 0.1) | 2,410 | 0.625 | 0.500 |
| Small (0.1-0.3) | 1,191 | 0.698 | 0.725 |
| Medium (0.3-0.5) | 575 | 0.727 | 0.739 |
| Large (0.5-1.0) | 541 | 0.733 | 0.722 |
| Very Large (> 1.0) | 371 | 0.742 | 0.741 |

**Key Insight**: Effects near zero have substantially lower R due to sign flip vulnerability.

---

## 5. Sign Flip Analysis

### 5.1 Prevalence
- **R <= 0.5 (potential sign instability)**: 1,943 meta-analyses (38.2%)

### 5.2 Sign Flip by Effect Type

| Effect Type | n with R <= 0.5 | Percentage |
|-------------|-----------------|------------|
| GIV | 976 | 72.3% |
| logRR | 916 | 25.8% |
| MD | 51 | 28.0% |

### 5.3 Sign Flip by Theta Near Zero

| Category | n with R <= 0.5 | Percentage |
|----------|-----------------|------------|
| Near zero (|theta| < 0.1) | 1,486 | 61.7% |
| Away from zero | 457 | 17.1% |

**Chi-square test**: X² = 1066.8, df = 1, p < 2.2e-16

---

## 6. Validation Against metafor

MA4 v1.0.1 was validated against the R metafor package (version 4.8-0) using REML estimation.

### Results (n = 100 random meta-analyses)

| Metric | Correlation | Mean Abs. Diff | Max Abs. Diff |
|--------|-------------|----------------|---------------|
| Theta (effect size) | 1.000000 | 0.000000 | 0.000006 |
| Sigma (SE) | 1.000000 | 0.000001 | 0.000031 |
| Tau (heterogeneity) | 1.000000 | 0.000053 | 0.002278 |

**Exact matches (diff < 0.01)**: 100/100 (100%)

---

## 7. Advanced Modeling

### 7.1 Quantile Regression

Coefficients differ across quantiles, revealing heterogeneous effects:

| Predictor | 25th Percentile | Median | 75th Percentile | OLS |
|-----------|-----------------|--------|-----------------|-----|
| log(k) | 0.000 | 0.000 | -0.059 | -0.037 |
| logRR effect | -0.016 | 0.173 | 0.235 | 0.138 |
| MD effect | -0.089 | 0.181 | 0.198 | 0.084 |
| log(tau) | -0.013 | -0.015 | -0.010 | -0.011 |
| theta_near_zero | -0.222 | -0.174 | 0.000 | -0.096 |

**Key Insight**: The theta_near_zero penalty is strongest at lower quantiles of R, confirming that near-zero effects predominantly drive low stability.

### 7.2 Hierarchical Model

Mixed effects model with random intercepts for Cochrane reviews:

- **Intraclass Correlation (ICC)**: 0.132
- **Interpretation**: 13.2% of variance in R is explained by between-review differences

### 7.3 Classification Model

Multinomial logistic regression for R categories (Low/Medium/High):

- **Training Accuracy**: 63.6%
- **Best predicted category**: Low R (sensitivity = 64.3%)

---

## 8. Practical Recommendations

Based on this comprehensive analysis:

1. **Low-k Meta-Analyses (k < 5)**: Interpret stability with caution; R values are highly sensitive to individual study influence

2. **Near-Zero Effects**: Effects with |theta| < 0.1 have 3.6x higher risk of sign instability; consider this when interpreting borderline significant results

3. **GIV Analyses**: Pre-computed effect sizes show systematically lower stability; prefer raw data when available

4. **High Heterogeneity**: tau > 0.5 reduces mean R by ~0.09 compared to homogeneous analyses

5. **Binary vs Continuous**: logRR analyses (binary outcomes) are most stable; MD analyses show moderate stability

6. **Optimal Stability Zone**: Meta-analyses with k >= 10, |theta| > 0.3, and tau < 0.3 have mean R > 0.75

---

## 9. Files Generated

| File | Description |
|------|-------------|
| `ma4_results_pairwise70.csv` | Full MA4 results for 5,088 meta-analyses |
| `ma4_metafor_validation.csv` | Validation comparison with metafor |
| `ma4_modeling_summary.txt` | Summary statistics text file |
| `ma4_analysis.R` | Core MA4 analysis script |
| `ma4_full_modeling.R` | Comprehensive regression modeling |
| `ma4_advanced_modeling.R` | Advanced modeling (quantile, mixed, classification) |

---

## 10. Technical Notes

### MA4 v1.0.1 Specification
- **Tau² Estimator**: REML via golden-section search (fallback: PM → DL → 0)
- **Perturbation Set Π v1**: FE, RE-DL, RE-PM, LOO-worst, drop-noisiest-20%
- **Sign Flip Penalty**: F = 0.5 if baseline sign = 0 or any perturbation flips sign
- **R Formula**: R = clip(R₀ × F), where R₀ = 1/(1 + Δ), Δ = max|θ_π - θ|/(1.96σ)

### Software Environment
- R version 4.5.2
- metafor package 4.8-0
- quantreg package (for quantile regression)
- lme4 package (for mixed effects)

---

*Report generated by MA4 v1.0.1 Analysis Pipeline*
