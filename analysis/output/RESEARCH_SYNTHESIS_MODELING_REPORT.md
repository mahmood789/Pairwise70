# Research Synthesis Modeling Framework
## Comprehensive Analysis Report

**Date:** January 2026
**Dataset:** 4,424 Cochrane Meta-Analyses from Pairwise70

---

## Executive Summary

This report presents a comprehensive modeling framework for research synthesis, developed using 4,424 Cochrane meta-analyses. The framework includes frequentist and Bayesian approaches to predict meta-analysis fragility and provides actionable recommendations for evidence synthesis.

---

## 1. Meta-Regression Models

### 1.1 Univariate Analysis (Effect Size Predictors)

| Predictor | Coefficient | SE | P-value | R² |
|-----------|-------------|-----|---------|-----|
| k | -0.0001 | 0.0005 | 0.753 | 0.0002 |
| log_k | -0.0048 | 0.0217 | 0.826 | 0.0001 |
| I² | -0.0010 | 0.0008 | 0.208 | 0.0034 |
| τ² | -0.0857 | 0.0525 | 0.103 | 0.0056 |

**Finding:** Individual predictors explain minimal variance in effect sizes alone.

### 1.2 Multivariate Fragility Models

#### Direction Fragility Model

| Variable | Estimate | SE | z-value | P-value |
|----------|----------|-----|---------|---------|
| Intercept | 1.407 | 0.124 | 11.33 | <0.001 |
| log(k) | -0.632 | 0.050 | -12.67 | <0.001 |
| I² | 0.019 | 0.002 | 10.19 | <0.001 |
| τ² | 0.204 | 0.020 | 10.46 | <0.001 |
| |Effect| | -5.528 | 0.308 | -17.97 | <0.001 |
| Significant | -3.320 | 0.510 | -6.51 | <0.001 |

#### Significance Fragility Model

| Variable | Estimate | SE | z-value | P-value |
|----------|----------|-----|---------|---------|
| Intercept | -1.407 | 0.131 | -10.74 | <0.001 |
| log(k) | -0.565 | 0.059 | -9.63 | <0.001 |
| I² | 0.020 | 0.001 | 13.27 | <0.001 |
| τ² | -0.005 | 0.015 | -0.36 | 0.722 |
| |Effect| | 0.006 | 0.114 | 0.05 | 0.958 |
| Significant | 1.371 | 0.123 | 11.11 | <0.001 |

---

## 2. Heterogeneity Partitioning

### 2.1 Heterogeneity Distribution by Meta-Analysis Size

| k Group | N MAs | Mean I² | Median I² | Mean τ² | % High Het |
|---------|-------|---------|-----------|---------|------------|
| 3-5 | 876 | 14.2% | 0% | 0.368 | 7.5% |
| 6-10 | 1,152 | 15.7% | 0% | 0.180 | 6.9% |
| 11-20 | 938 | 17.2% | 0% | 0.235 | 6.1% |
| 21-50 | 637 | 21.8% | 0% | 0.184 | 9.3% |
| >50 | 287 | 25.9% | 7.7% | 0.257 | 11.8% |

**Key Finding:** Larger meta-analyses tend to have higher heterogeneity but more stable conclusions.

### 2.2 Fragility by Heterogeneity Level

| I² Category | N | Direction Fragile | Significance Fragile | Mean Effect |
|-------------|---|-------------------|---------------------|-------------|
| Low (<25%) | 3,261 | 28.2% | 9.2% | 0.267 |
| Moderate (25-50%) | 436 | 22.5% | 31.4% | 0.441 |
| Substantial (50-75%) | 401 | 25.3% | 30.4% | 0.534 |
| Considerable (>75%) | 326 | 20.1% | 33.0% | 0.787 |

---

## 3. Predictive Modeling Performance

### 3.1 Cross-Validated AUC Comparison

| Model | Mean CV AUC | SD |
|-------|-------------|-----|
| Logistic Regression | 0.688 | 0.023 |
| LASSO | 0.688 | 0.025 |
| Random Forest | 0.685 | 0.028 |

### 3.2 Information Criteria

| Model | AIC | BIC | Deviance |
|-------|-----|-----|----------|
| Null | 4076.0 | 4082.0 | 4074.0 |
| k only | 3926.7 | 3938.8 | 3922.7 |
| k + I² | 3847.4 | 3865.4 | 3841.4 |
| Full | 3735.1 | 3771.2 | 3723.1 |

**Best Model:** Full model (log_k + I² + τ² + |effect| + significant)

---

## 4. Bayesian Analysis

### 4.1 Posterior Estimates (95% Credible Intervals)

| Parameter | MAP | 95% CI | P(β > 0) |
|-----------|-----|--------|----------|
| Intercept | 1.134 | [0.946, 1.322] | 1.000 |
| log(k) | -0.610 | [-0.686, -0.533] | 0.000 |
| I²/100 | 1.510 | [1.210, 1.810] | 1.000 |
| τ² | 0.385 | [0.219, 0.550] | 1.000 |
| |Effect| | -1.431 | [-1.642, -1.220] | 0.000 |

### 4.2 Bayesian Model Averaging

| Model | Predictors | BIC | Posterior Prob |
|-------|------------|-----|----------------|
| M8 | log_k + I² + τ² + |effect| | 5359.9 | 0.999 |
| M6 | log_k + I² + |effect| | 5373.7 | 0.001 |
| Others | - | - | <0.001 |

**BMA AUC:** 0.693

### 4.3 Prior Sensitivity Analysis

| Prior Type | log_k | I² | τ² | |Effect| |
|------------|-------|----|----|---------|
| Weakly Informative | -0.610 | 1.510 | 0.385 | -1.431 |
| Informative | -0.596 | 1.405 | 0.388 | -1.364 |
| Diffuse | -0.614 | 1.548 | 0.383 | -1.454 |

**Conclusion:** Results are robust to prior specification (CV < 0.05 for all coefficients).

### 4.4 Model Calibration

**Posterior Predictive p-value:** 0.507 (excellent fit)

---

## 5. Risk Stratification

### 5.1 Bayesian Risk Categories

| Risk Category | N | % | Actual Fragility | Predicted Prob | Uncertainty |
|---------------|---|---|-----------------|----------------|-------------|
| Very Low | 453 | 10.5% | 16.8% | 0.141 | 0.051 |
| Low | 1,638 | 38.0% | 30.2% | 0.309 | 0.050 |
| Moderate | 1,746 | 40.5% | 47.4% | 0.492 | 0.052 |
| High | 420 | 9.7% | 73.8% | 0.672 | 0.073 |
| Very High | 59 | 1.4% | 88.1% | 0.860 | 0.082 |

---

## 6. Key Findings and Recommendations

### 6.1 Strongest Predictors of Fragility

1. **Number of Studies (log_k):** Most protective factor
   - Each doubling of k reduces fragility odds by 45%
   - Odds Ratio: 0.54 [0.50, 0.59]

2. **Effect Magnitude (|estimate|):** Larger effects more stable
   - Each unit increase reduces direction fragility odds by 77%
   - Odds Ratio: 0.24 [0.19, 0.30]

3. **Heterogeneity (I²):** Increases fragility risk
   - Each 10% increase in I² increases fragility odds by 16%
   - Odds Ratio: 1.16 [1.13, 1.20]

4. **Between-Study Variance (τ²):** Increases fragility
   - Odds Ratio: 1.47 [1.25, 1.73]

### 6.2 Practical Recommendations

#### For Systematic Reviewers:
- **Include k ≥ 10 studies** for reliable conclusions
- **Report fragility indices** alongside pooled estimates
- **Account for heterogeneity** when interpreting fragility

#### For Guideline Developers:
- **Use risk stratification** to weight evidence quality
- **Apply uncertainty-aware** recommendations for moderate-risk MAs
- **Require additional evidence** for high-fragility conclusions

#### For Methodologists:
- **Preferred model:** Logistic regression (simple, interpretable, competitive)
- **Bayesian approach** provides uncertainty quantification
- **BMA unnecessary:** Full model dominates (99.9% posterior probability)

---

## 7. Files Generated

| File | Description |
|------|-------------|
| `RESEARCH_SYNTHESIS_MODELING.R` | Frequentist modeling framework |
| `BAYESIAN_SYNTHESIS_MODELING.R` | Bayesian modeling framework |
| `model_comparison.csv` | Model AUC comparison |
| `cv_results.csv` | Cross-validation results |
| `final_model_coefficients.csv` | Final logistic model coefficients |
| `bayesian_posterior_summary.csv` | Bayesian posterior estimates |
| `bayesian_model_comparison.csv` | BMA model comparison |
| `bayesian_prior_sensitivity.csv` | Prior sensitivity results |

---

## 8. Technical Notes

- **Software:** R 4.5.2 with metafor, data.table, glmnet, ranger
- **Dataset:** Pairwise70 (501 Cochrane systematic reviews, 4,424 meta-analyses)
- **Validation:** 10-fold cross-validation for frequentist; posterior predictive checks for Bayesian
- **Priors:** Weakly informative Normal(0, 1) for coefficients, Normal(0, 2.5) for intercept

---

*Report generated as part of comprehensive research synthesis methodology development.*
