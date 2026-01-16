# V4 Advanced Pooling Methods - Quick Guide

**Package:** Pairwise70 v2.0.0
**Date:** January 15, 2026

---

## Overview

This package now includes **13 advanced pooling methods** for meta-analysis:
- **V3 Methods (3):** MWM, SIT, ARP
- **V4 Methods (10):** WRD, CBM, RBM, SWA, TAS, EVE, PVM, AEM, SPE, SMS

---

## V4 Methods by Category

### Category A: Robustness Methods

| Method | Full Name | Best For | Min k |
|--------|-----------|----------|-------|
| **WRD** | Winsorized Robust Downsampling | Extreme outliers | 4 |
| **CBM** | Clustering-Based Meta-Analysis | Multimodal distributions | 6 |
| **RBM** | Residual-Based Weighting | Automatic downweighting | 4 |

### Category B: Bias Correction Methods

| Method | Full Name | Best For | Min k |
|--------|-----------|----------|-------|
| **SWA** | Selection-Weight Adjustment | Publication bias | 10 |
| **TAS** | Trim-and-Shrink | Publication bias + stability | 10 |

### Category C: Advanced Variance Methods

| Method | Full Name | Best For | Min k |
|--------|-----------|----------|-------|
| **EVE** | Empirical Bayes Variance Estimation | Better τ² estimation | 3 |
| **PVM** | Profile Variance Minimization | Optimal CI coverage | 5 |

### Category D: Ensemble & Adaptive Methods

| Method | Full Name | Best For | Min k |
|--------|-----------|----------|-------|
| **AEM** | Adaptive Ensemble Meta-Analysis | Smart method combination | 5 |
| **SPE** | Stochastic Point Estimator | Full uncertainty quantification | 3 |

### Category E: Specialized Methods

| Method | Full Name | Best For | Min k |
|--------|-----------|----------|-------|
| **SMS** | Small-Meta Shrinkage | Small meta-analyses (k < 10) | 3 |

---

## Quick Start

```r
# Install and load
devtools::install_github("mahmood789/Pairwise70")
library(Pairwise70)

# Example data
yi <- c(0.25, 0.30, 0.28, 0.35, 0.22, 0.31, 0.29, 0.33, 0.27, 0.26)
vi <- c(0.02, 0.01, 0.02, 0.01, 0.03, 0.02, 0.01, 0.02, 0.02, 0.01)

# Run a single method
result <- wrd_meta(yi, vi)
print(result$estimate)  # Point estimate
print(result$ci_lb)     # Lower CI
print(result$ci_ub)     # Upper CI

# Compare all methods
comparison <- compare_methods_v4(yi, vi)
print(comparison)
```

---

## Method Recommendations

| Scenario | Recommended Method | Reason |
|----------|-------------------|--------|
| Standard meta-analysis | MWM or AEM | Best overall performance |
| Suspected outliers | WRD or SIT | Robust to extreme values |
| Multimodal effects | CBM | Detects subgroups |
| Publication bias concern | SWA or TAS | Direct bias correction |
| Small k (k < 10) | SMS or MWM | Optimized for small samples |
| High heterogeneity | RBM or EVE | Handles τ² uncertainty |
| Confirmatory analysis | HKSJ | Well-established |
| Exploratory analysis | AEM | Adaptive, comprehensive |

---

## Output Structure

All methods return a list with:

```r
list(
  estimate = numeric,     # Point estimate
  se = numeric,          # Standard error
  ci_lb = numeric,       # Lower 95% CI
  ci_ub = numeric,       # Upper 95% CI
  pvalue = numeric,      # P-value for H₀: θ = 0
  tau2 = numeric,        # Heterogeneity variance
  method = character,    # Method name
  k = integer,           # Number of studies
  ...                    # Method-specific diagnostics
)
```

---

## Method-Specific Diagnostics

### WRD
- `n_winsorized`: Number of winsorized studies
- `huber_weights`: Final Huber weights
- `se_method`: "bootstrap" or "sandwich"

### CBM
- `n_clusters`: Number of clusters detected
- `cluster_sizes`: Size of each cluster
- `cluster_assignments`: Which cluster each study belongs to
- `cluster_estimates`: Estimate for each cluster

### RBM
- `weights`: Final study weights
- `converged`: Whether iteration converged
- `iterations`: Number of iterations

### SWA
- `selection_function`: "step" or "continuous"
- `selection_probabilities`: Probability for each study
- `n_boot`: Number of bootstrap iterations

### TAS
- `n_imputed`: Number of imputed studies
- `tf_estimate`: Trim-and-fill estimate
- `shrunk_estimate`: Shrunken estimate
- `shrinkage_factor`: Shrinkage strength

### EVE
- `tau2_estimates`: All 5 τ² estimates (REML, DL, PM, ML, SJ)
- `tau2_eb`: Empirical Bayes combined τ²

### PVM
- `tau2_optimal`: Selected τ² from profiling
- `achieved_coverage`: Actual coverage achieved
- `coverage_by_tau`: Coverage for each τ² value

### AEM
- `features`: Data features used for weighting
- `base_estimates`: Estimates from base methods
- `ensemble_weights`: Final ensemble weights

### SPE
- `tau2_samples`: MCMC samples of τ²
- `theta_samples`: MCMC samples of θ
- `accept_rate`: MCMC acceptance rate
- `n_samples`: Number of MCMC samples

### SMS
- `shrinkage_factor`: Amount of shrinkage applied
- `standard_estimate`: Unshrunken estimate
- `shrunk_estimate`: Final shrunken estimate

---

## Performance Summary (Based on Simulation Design)

| Method | Expected Strengths |
|--------|-------------------|
| WRD | Handles extreme outliers without removing studies |
| CBM | Detects and models subgroups |
| RBM | Automatic downweighting of influential studies |
| SWA | Direct publication bias modeling |
| TAS | Conservative bias correction |
| EVE | Stable τ² estimation |
| PVM | Optimizes for target coverage |
| AEM | Adaptive, good overall performance |
| SPE | Full uncertainty quantification |
| SMS | Optimized for small k |

---

## Common Pitfalls

1. **Too few studies**: Most methods need k ≥ 5. Use HKSJ for k < 5.
2. **Publication bias**: SWA and TAS need k ≥ 10.
3. **Computation time**: PVM and SPE are slower. Use reduced iterations for quick analysis.
4. **Convergence**: RBM may not converge with all identical effect sizes.

---

## Citation

If you use V4 methods, please cite:

```
Arai M. (2026). Pairwise70: Comprehensive Cochrane Pairwise Meta-Analysis
Dataset Collection with 13 Advanced Pooling Methods. R package version 2.0.0.
https://github.com/mahmood789/Pairwise70
```

---

## References

- Huber, P. J. (1964). Robust estimation of a location parameter.
- Vevea, J. L., & Hedges, L. V. (1995). A general linear model for estimating effect size heterogeneity.
- Efron, B., & Morris, C. (1975). Data analysis using Stein's estimator.
- Gelman, A. (2006). Prior distributions for variance parameters in hierarchical models.
- Wolpert, D. H. (1992). Stacked generalization.

---

**Last Updated:** January 15, 2026
