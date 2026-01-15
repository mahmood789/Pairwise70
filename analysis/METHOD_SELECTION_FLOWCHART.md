# Method Selection Flowchart for Advanced Pooling Methods

## Visual Flowchart (ASCII)

```
                            ┌─────────────────────┐
                            │   START: Have MA    │
                            │   with k studies    │
                            └──────────┬──────────┘
                                       │
                                       ▼
                            ┌─────────────────────┐
                            │     Is k ≥ 5?       │
                            └──────────┬──────────┘
                                       │
                    ┌──────────────────┴──────────────────┐
                    │ NO                                   │ YES
                    ▼                                      ▼
         ┌─────────────────────┐              ┌─────────────────────┐
         │   Use HKSJ only     │              │  Suspect outliers?  │
         │ (methods unstable   │              │  (check residuals)  │
         │  for very small k)  │              └──────────┬──────────┘
         └─────────────────────┘                         │
                                          ┌──────────────┴──────────────┐
                                          │ YES                          │ NO
                                          ▼                              ▼
                               ┌─────────────────────┐      ┌─────────────────────┐
                               │      Use SIT_v3     │      │ Suspect publication │
                               │  - Auto-trims       │      │      bias?          │
                               │  - Bootstrap SE     │      └──────────┬──────────┘
                               │  - Report both full │                 │
                               │    and trimmed      │      ┌──────────┴──────────┐
                               └─────────────────────┘      │ YES                  │ NO
                                                            ▼                      ▼
                                               ┌─────────────────────┐ ┌─────────────────────┐
                                               │      Use SIT_v3     │ │   High I² (>75%)?   │
                                               │  + Trim-and-fill    │ └──────────┬──────────┘
                                               │  + Selection models │            │
                                               └─────────────────────┘ ┌──────────┴──────────┐
                                                                       │ YES                  │ NO
                                                                       ▼                      ▼
                                                          ┌─────────────────────┐ ┌─────────────────────┐
                                                          │     Use ARP_v3      │ │     Use MWM_v3      │
                                                          │  - Multiple         │ │  (default choice)   │
                                                          │    estimators       │ │  - Best overall     │
                                                          │  - Robust to        │ │    RMSE             │
                                                          │    heterogeneity    │ │  - Stability        │
                                                          └─────────────────────┘ │    weighting        │
                                                                                  └─────────────────────┘
```

## Decision Table

| Condition | Recommended Method | Rationale |
|-----------|-------------------|-----------|
| k < 5 | HKSJ | LOO unstable; methods fall back automatically |
| k ≥ 5, suspected outliers | SIT_v3 | Identifies and trims outliers; bootstrap SE |
| k ≥ 5, suspected pub bias | SIT_v3 + trim-and-fill | SIT mitigates some bias; combine with other methods |
| k ≥ 5, high heterogeneity (I² > 75%) | ARP_v3 | Combines multiple estimators; robust |
| k ≥ 5, standard conditions | MWM_v3 | Best overall RMSE; stability weighting |
| Any k, confirmatory analysis | Consider RoBMA | Principled Bayesian inference |
| Many meta-analyses | MWM_v3 or ARP_v3 | Computational efficiency |

## Quick Reference

### When to Use Each Method

**MWM_v3 (MAFI-Weighted Meta-Analysis)**
- Default choice for most meta-analyses
- Best when: k ≥ 5, no obvious outliers, moderate heterogeneity
- Advantage: Best overall RMSE, stability weighting

**SIT_v3 (Sequential Influence Trimming)**
- Use when: Suspected outliers OR publication bias
- Advantage: Auto-identifies influential studies
- Note: Report both full and trimmed analyses

**ARP_v3 (Adaptive Robust Pooling)**
- Use when: High heterogeneity, want robustness
- Advantage: Combines REML, DL, PM estimators
- Note: Good when unsure which estimator is best

**HKSJ (Hartung-Knapp-Sidik-Jonkman)**
- Use when: k < 5 (automatic fallback)
- Advantage: Widely validated, conservative
- Note: Our methods use HKSJ for very small k

**RVE (Robust Variance Estimation)**
- Use when: Correlated effects, multi-arm trials
- Advantage: Handles dependence
- Note: Requires cluster structure

## Diagnostic Checks Before Choosing

### 1. Check for Outliers
```r
fit <- rma(yi, vi)
rstudent(fit)  # |z| > 2.5 suggests outlier
```

### 2. Check for Publication Bias
```r
regtest(fit)  # Egger's test
funnel(fit)   # Visual inspection
```

### 3. Check Heterogeneity
```r
fit$I2   # I² statistic
fit$QE   # Cochran's Q
```

### 4. Compare Methods
```r
source("R/advanced_pooling_v3.R")
compare_methods_v3(yi, vi)  # Side-by-side comparison
```

## Implementation Flowchart (Code)

```r
select_method <- function(yi, vi) {
  k <- length(yi)

  # Step 1: Check k
  if (k < 5) {
    return("HKSJ")
  }

  # Step 2: Fit baseline model
  fit <- rma(yi, vi, method = "REML")

  # Step 3: Check for outliers
  z_scores <- abs(rstudent(fit)$z)
  has_outliers <- any(z_scores > 2.5)

  # Step 4: Check for publication bias
  egger <- regtest(fit)
  has_pubbias <- egger$pval < 0.10

  # Step 5: Check heterogeneity
  high_het <- fit$I2 > 75

  # Step 6: Select method
  if (has_outliers || has_pubbias) {
    return("SIT_v3")
  } else if (high_het) {
    return("ARP_v3")
  } else {
    return("MWM_v3")
  }
}
```

## Reporting Recommendations

1. **Always report** standard REML/HKSJ alongside novel methods
2. **For SIT**: Report both full (k=original) and trimmed (k=final) results
3. **Sensitivity**: If methods disagree substantially, discuss why
4. **Transparency**: State which method was pre-specified vs. exploratory

---

## References

- Hartung, J., & Knapp, G. (2001). On tests of the overall treatment effect in meta-analysis. *Statistics in Medicine*, 20(12), 1771-1782.
- IntHout, J., et al. (2014). The Hartung-Knapp-Sidik-Jonkman method. *BMC Medical Research Methodology*, 14(1), 25.
- Viechtbauer, W., & Cheung, M. W. L. (2010). Outlier and influence diagnostics for meta-analysis. *Research Synthesis Methods*, 1(2), 112-125.
