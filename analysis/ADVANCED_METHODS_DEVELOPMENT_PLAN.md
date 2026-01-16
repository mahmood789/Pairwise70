# Advanced Pooling Methods Development Plan v2.0
## New Methods & Comprehensive Testing Framework

**Date Created:** January 15, 2026
**Status:** Planning Phase
**Package:** Pairwise70 v1.1.0+

---

## Executive Summary

This plan outlines the development of **10 novel advanced pooling methods** for meta-analysis, building on the existing MWM, SIT, and ARP methods. Each new method addresses specific limitations identified in the current implementation and targets unexplored areas in meta-analytic methodology.

---

## Part 1: New Methods to Develop

### Category A: Robustness Methods (3 new methods)

#### A1. **WRD** - Winsorized Robust Downsampling
**Purpose:** Address extreme outliers without removing studies entirely

**Algorithm:**
1. Identify studies with |studentized residual| > 3.0
2. Instead of trimming, winsorize to 3.0 (replace extreme values)
3. Apply robust downsampling: overweight small studies by sqrt(N)
4. Use Huber M-estimator for final pooling
5. Bootstrap SE correction (n=499)

**Parameters:**
- `winsor_threshold = 3.0`
- `downweight_power = 0.5` (sqrt)
- `huber_k = 1.345`

**Expected Benefits:**
- Retains all studies (addresses SIT's selection problem)
- Robust to extreme outliers
- Better small-study representation

---

#### A2. **CBM** - Clustering-Based Meta-Analysis
**Purpose:** Handle multimodal effect distributions and subgroups

**Algorithm:**
1. Apply Gaussian Mixture Model to effect sizes
2. Determine optimal number of clusters (BIC)
3. Estimate separate τ² for each cluster
4. Pool within clusters, then combine using cluster weights
5. Cross-validate cluster stability (bootstrap)

**Parameters:**
- `max_clusters = 5`
- `min_cluster_size = 3`
- `cluster_method = "mixture"` or `"hierarchical"`

**Expected Benefits:**
- Detects hidden subgroups
- Handles heterogeneous populations
- Provides cluster-level diagnostics

---

#### A3. **RBM** - Residual-Based Weighting Method
**Purpose:** Iteratively reweight based on residual patterns

**Algorithm:**
1. Initial REML fit
2. Calculate squared studentized residuals: r²i
3. Update weights: wi = 1 / (vi × max(1, r²i/median(r²)))
4. Refit with new weights (iteratively until convergence)
5. Use sandwich variance estimator for SE
6. T-distribution CIs with adjusted df

**Parameters:**
- `max_iterations = 25`
- `convergence_tol = 1e-6`
- `weight_cap = 10` (max weight multiplier)

**Expected Benefits:**
- Automatically downweights poorly fitting studies
- No hard trimming (avoids selection problem)
- Theoretically grounded (iterative reweighted least squares)

---

### Category B: Bias Correction Methods (2 new methods)

#### B1. **SWA** - Selection-Weight Adjustment
**Purpose:** Model publication bias directly without regression

**Algorithm:**
1. Fit weight-function model (Vevea & Hedges, 1995)
2. Model selection probability as function of p-value
3. Adjust each study's contribution by inverse selection probability
4. Estimate corrected effect via weighted estimator
5. Delta method for variance

**Parameters:**
- `selection_function = "step"` or `"continuous"`
- `p_cutoff = 0.05` (for step function)
- `n_boot = 1000` (for SE estimation)

**Expected Benefits:**
- Directly models publication process
- Better than PET-PEESE for strong publication bias
- Provides selection probability estimates

---

#### B2. **TAS** - Trim-and-Shrink
**Purpose:** Combination of trimming and shrinkage estimation

**Algorithm:**
1. Run standard trim-and-fill
2. Instead of just adding imputed studies:
   - Shrink imputed studies toward pooled mean
   - Apply empirical Bayes shrinkage factor
3. Final estimate: weighted combination of original + shrunk
4. Adjust SE for imputation uncertainty

**Parameters:**
- `shrinkage_prior = "empirical_bayes"`
- `trim_estimator = "R0"` or `"L0"` or `"Q0"`
- `shrinkage_weight = 0.5`

**Expected Benefits:**
- More conservative than standard trim-and-fill
- Accounts for imputation uncertainty
- Combines frequency + Bayesian approaches

---

### Category C: Advanced Variance Methods (2 new methods)

#### C1. **EVE** - Empirical Bayes Variance Estimation
**Purpose:** Improved τ² estimation with shrinkage

**Algorithm:**
1. Compute multiple τ² estimators (REML, DL, PM, ML, SJ)
2. Use empirical Bayes to combine them
3. Prior: inverse-gamma on τ²
4. Posterior mean of τ² as final estimate
5. Full Bayesian propagation of uncertainty

**Parameters:**
- `prior_shape = 0.001`
- `prior_rate = 0.001` (weakly informative)
- `n_mcmc = 5000` (for uncertainty quantification)

**Expected Benefits:**
- Better τ² estimation for small k
- Quantifies uncertainty in τ²
- Stable even with extreme heterogeneity

---

#### C2. **PVM** - Profile Variance Minimization
**Purpose:** Choose τ² that directly optimizes CI properties

**Algorithm:**
1. Define objective: maximize coverage while minimizing CI width
2. Profile over τ² grid
3. Choose τ²* that achieves 95% empirical coverage
4. Use bootstrap to estimate coverage function
5. Interpolate for continuous τ² estimate

**Parameters:**
- `tau_grid_min = 0`
- `tau_grid_max = 1`
- `n_grid_points = 100`
- `coverage_target = 0.95`

**Expected Benefits:**
- Directly optimizes for target coverage
- Data-driven τ² selection
- Theoretically motivated CI construction

---

### Category D: Ensemble & Adaptive Methods (2 new methods)

#### D1. **AEM** - Adaptive Ensemble Meta-Analysis
**Purpose:** Smartly combine methods based on data characteristics

**Algorithm:**
1. Compute data features: k, mean(SE), skew(yi), kurtosis(yi), I²
2. Train meta-model on simulation to predict best method
3. Apply adaptive weights: larger weight to predicted best method
4. Use stacking (super-learner) approach
5. Cross-validation to prevent overfitting

**Parameters:**
- `base_methods = c("REML", "HKSJ", "MWM", "SIT", "ARP")`
- `meta_model = "random_forest"` or `"elastic_net"`
- `n_folds = 5`

**Expected Benefits:**
- Automatically selects appropriate method
- Better average performance across scenarios
- Machine learning grounded

---

#### D2. **SPE** - Stochastic Point Estimator
**Purpose:** Quantify estimation uncertainty via stochastic sampling

**Algorithm:**
1. Treat τ² as random variable with posterior
2. Sample τ²(s) from posterior (MCMC)
3. For each τ²(s), compute pooled estimate θ̂(s)
4. Final estimate: median of θ̂(s)
5. SE: std dev of θ̂(s) (includes τ² uncertainty)
6. CI: quantile method

**Parameters:**
- `n_samples = 10000`
- `burnin = 1000`
- `tau_prior = "half_cauchy"` or `"inv_gamma"`

**Expected Benefits:**
- Fully propagates τ² uncertainty to pooled estimate
- Proper uncertainty quantification
- Natural Bayesian interpretation

---

### Category E: Specialized Methods (1 new method)

#### E1. **SMS** - Small-Meta Shrinkage
**Purpose:** Specialized for k < 10 meta-analyses

**Algorithm:**
1. Use hierarchical shrinkage prior on θ
2. Prior: θ ~ N(μ₀, τ²global) where τ²global estimated from all meta-analyses
3. Empirical Bayes: borrow strength across meta-analyses
4. Shrinkage toward global mean based on k and precision
5. Posterior predictive for new study

**Parameters:**
- `global_prior = "hierarchical"`
- `shrinkage_intensity = "adaptive"`
- `min_k_for_shrinkage = 10`

**Expected Benefits:**
- Optimized for small k (common in practice)
- Borrows strength from historical data
- Reduces overfitting in small samples

---

## Part 2: Comprehensive Testing Framework

### 2.1 Simulation Infrastructure

**File:** `analysis/simulation/Comprehensive_Testing_Framework.R`

```r
# Core simulation function
run_comprehensive_simulation <- function(
    methods = list_all_methods(),
    scenarios = all_scenarios(),
    n_sim = 1000,
    n_cores = parallel::detectCores() - 1,
    seed = 20260115
) {
    # Parallel simulation across all scenarios
    # Return: tidy data frame with all results
}
```

---

### 2.2 Expanded Scenarios (25 total)

#### Baseline Scenarios (5)
| ID | k | τ² | Effect | Description |
|----|---|----|----|-------------|
| B1 | 5 | 0.05 | 0.0 | Null effect, small k |
| B2 | 10 | 0.05 | 0.3 | Standard, medium k |
| B3 | 20 | 0.05 | 0.5 | Standard, large k |
| B4 | 50 | 0.05 | 0.3 | Large k |
| B5 | 100 | 0.05 | 0.3 | Very large k |

#### Heterogeneity Scenarios (5)
| ID | k | τ² | Effect | Description |
|----|---|----|----|-------------|
| H1 | 10 | 0.0 | 0.3 | No heterogeneity |
| H2 | 10 | 0.10 | 0.3 | Moderate heterogeneity |
| H3 | 10 | 0.20 | 0.3 | High heterogeneity |
| H4 | 10 | 0.50 | 0.3 | Very high heterogeneity |
| H5 | 10 | 1.00 | 0.3 | Extreme heterogeneity |

#### Outlier Scenarios (5)
| ID | k | τ² | Outlier Type | Description |
|----|---|----|----|-------------|
| O1 | 10 | 0.05 | 1 study at +3SD | Single mild outlier |
| O2 | 10 | 0.05 | 1 study at +5SD | Single extreme outlier |
| O3 | 10 | 0.05 | 2 studies at ±3SD | Two opposing outliers |
| O4 | 15 | 0.05 | 3 studies at +3SD | Multiple outliers |
| O5 | 10 | 0.05 | 10% at +4SD | Clustered outliers |

#### Publication Bias Scenarios (5)
| ID | k | Bias Mechanism | Severity |
|----|---|----|----|
| PB1 | 20 | p < 0.10 excluded | Mild |
| PB2 | 20 | p < 0.05 excluded | Moderate |
| PB3 | 20 | p < 0.01 excluded | Strong |
| PB4 | 20 | Step function at p=0.05 | Classic |
| PB5 | 20 | Continuous weight function | Realistic |

#### Small Study Scenarios (3)
| ID | k | Mean N | Description |
|----|---|----|----|
| S1 | 10 | N=20 per arm | Very small studies |
| S2 | 10 | N=50 per arm | Small studies |
| S3 | 10 | N=20-500 (bimodal) | Mixed sizes |

#### Distribution Scenarios (2)
| ID | True Effect Dist | Description |
|----|----|----|
| D1 | t(3) | Heavy-tailed effects |
| D2 | Skewed (lognormal) | Skewed effects |

---

### 2.3 Performance Metrics (8 metrics)

For each method-scenario combination, compute:

1. **Bias**: |mean(θ̂) - θ_true|
2. **RMSE**: sqrt(mean((θ̂ - θ_true)²))
3. **Coverage**: Proportion of CIs containing θ_true
4. **CI Width**: Mean(CI_ub - CI_lb)
5. **Coverage × Width**: Product (efficiency-adjusted coverage)
6. **Type I Error**: When θ_true = 0, proportion p < 0.05
7. **Power**: When θ_true ≠ 0, proportion p < 0.05
8. **Stability**: Variance of θ̂ across bootstrap samples

---

### 2.4 Testing Phases

**Phase 1: Unit Tests (Week 1)**
- Test each function with toy data
- Check for NA handling, edge cases
- Validate output structure
- File: `tests/testthat/test-new-methods.R`

**Phase 2: Small-Scale Validation (Week 2)**
- 100 iterations × 25 scenarios
- Quick sanity check on all methods
- Identify catastrophic failures early
- File: `analysis/simulation/Quick_Validation_V4.R`

**Phase 3: Medium-Scale Testing (Week 3-4)**
- 500 iterations × 25 scenarios
- Compare new methods to existing (MWM, SIT, ARP)
- Initial performance ranking
- File: `analysis/simulation/Medium_Scale_V4.R`

**Phase 4: Full Simulation Study (Week 5-6)**
- 2000 iterations × 25 scenarios = 50,000 runs per method
- Definitive performance comparison
- Generate publication-ready tables/figures
- File: `analysis/simulation/Full_Simulation_V4.R`

**Phase 5: Real Data Validation (Week 7)**
- Apply to all 501 Pairwise70 datasets
- Compare to REML/HKSJ on real data
- Identify where methods disagree
- Case studies of interesting discrepancies
- File: `analysis/simulation/Real_Data_Validation.R`

---

### 2.5 Benchmarking & Diagnostics

**Computational Benchmark**
```r
benchmark_methods <- function(yi, vi, n_reps = 100) {
    # Measure runtime for each method
    # Report: mean time, std dev, memory usage
}
```

**Convergence Monitoring**
```r
check_convergence <- function(method_result) {
    # Return: convergence_status, warnings, diagnostics
    # Auto-flag problematic fits
}
```

**Robustness Dashboard**
```r
generate_diagnostics <- function(method_result) {
    # Create diagnostic plots:
    # - Weight distribution
    # - Influence plot
    # - Q-Q plot of residuals
    # - Coverage by scenario
}
```

---

## Part 3: Implementation Timeline

### Week 1: Setup & Unit Tests
- [ ] Set up testing infrastructure
- [ ] Write unit tests for all 10 new methods
- [ ] Set up parallel simulation framework
- [ ] Create result storage templates

### Week 2: Implement Robustness Methods (A1-A3)
- [ ] Implement WRD (Winsorized Robust Downsampling)
- [ ] Implement CBM (Clustering-Based Meta-Analysis)
- [ ] Implement RBM (Residual-Based Weighting)
- [ ] Run quick validation (100 iterations)

### Week 3: Implement Bias Correction (B1-B2)
- [ ] Implement SWA (Selection-Weight Adjustment)
- [ ] Implement TAS (Trim-and-Shrink)
- [ ] Run quick validation

### Week 4: Implement Variance Methods (C1-C2)
- [ ] Implement EVE (Empirical Bayes Variance Estimation)
- [ ] Implement PVM (Profile Variance Minimization)
- [ ] Run quick validation

### Week 5: Implement Ensemble Methods (D1-D2)
- [ ] Implement AEM (Adaptive Ensemble)
- [ ] Implement SPE (Stochastic Point Estimator)
- [ ] Run quick validation

### Week 6: Implement Specialized Method & Full Testing
- [ ] Implement SMS (Small-Meta Shrinkage)
- [ ] Begin full simulation (2000 iterations)

### Weeks 7-8: Complete Full Simulation
- [ ] Complete full simulation study
- [ ] Analyze results
- [ ] Generate tables and figures

### Week 9: Real Data Validation
- [ ] Apply to all 501 Pairwise70 datasets
- [ ] Case study analysis
- [ ] Comparison to ground truth

### Week 10: Documentation & Package Integration
- [ ] Write comprehensive documentation
- [ ] Update man pages
- [ ] Create vignettes
- [ ] Prepare package v2.0

---

## Part 4: Success Criteria

A method is considered successful if it meets **at least 3 of 5** criteria:

1. **Coverage**: 93-97% average coverage across scenarios
2. **RMSE**: ≤10% improvement over REML in ≥30% of scenarios
3. **Robustness**: ≤5% increase in RMSE in outlier scenarios
4. **Efficiency**: CI width ≤110% of HKSJ in standard scenarios
5. **Stability**: Convergence rate ≥99% across all scenarios

**Overall project success:**
- At least 7 of 10 new methods meet success criteria
- At least 2 new methods outperform existing MWM/SIT/ARP
- Complete simulation study with publication-ready results

---

## Part 5: File Structure

```
Pairwise70/
├── R/
│   ├── advanced_pooling_v4.R          # All 10 new methods
│   └── testing_framework.R             # Testing utilities
├── analysis/
│   ├── simulation/
│   │   ├── Comprehensive_Testing_Framework.R
│   │   ├── Quick_Validation_V4.R
│   │   ├── Medium_Scale_V4.R
│   │   ├── Full_Simulation_V4.R
│   │   └── Real_Data_Validation.R
│   └── results/
│       ├── v4_unit_tests.rds
│       ├── v4_quick_validation.rds
│       ├── v4_medium_scale.rds
│       ├── v4_full_simulation.rds
│       └── v4_real_data.rds
├── tests/
│   └── testthat/
│       └── test-new-methods.R
└── man/
    ├── wrd_meta.Rd
    ├── cbm_meta.Rd
    ├── rbw_meta.Rd
    ├── swa_meta.Rd
    ├── tas_meta.Rd
    ├── eve_meta.Rd
    ├── pvm_meta.Rd
    ├── aem_meta.Rd
    ├── spe_meta.Rd
    └── sms_meta.Rd
```

---

## Part 6: References for Implementation

1. **Vevea & Hedges (1995)** - Weight-function model for publication bias
2. **Huber (1964)** - Robust M-estimation
3. **Efron (2010)** - Empirical Bayes methods
4. **Van der Pas et al. (2014)** - Shrinkage priors for variance
5. **Hedges (1992)** - Modeling publication selection
6. **Rücker et al. (2011)** - Network meta-analysis clustering
7. **White (2011)** - Multivariate meta-analysis
8. **Jackson et al. (2010)** - Multivariate meta-analysis for heterogeneity
9. **Borenstein et al. (2010)** - Small sample corrections
10. **Liu et al. (2017)** - Clustering in meta-analysis

---

## Part 7: Risk Assessment & Mitigation

| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| Method convergence failure | Medium | High | Fallback to HKSJ |
| Simulation too slow | High | Medium | Parallel processing |
| New methods underperform | Medium | High | Document anyway (negative results) |
| Computational resources | Low | Medium | Use cloud if needed |
| Overfitting to simulation | Medium | High | Separate training/test scenarios |

---

**Document Version:** 1.0
**Status:** Ready for Implementation
**Next Step:** Begin Week 1 tasks - Setup & Unit Tests

---

*Generated with Claude Code*
