# Revision Plan: Addressing Editorial Review

**Status:** Action items prioritized by importance
**Estimated Total Time:** 2-3 weeks

---

## PRIORITY 1: CRITICAL (Must Do)

### 1.1 Extend Simulation to 1000+ Iterations

**Issue:** 200 iterations insufficient for reliable coverage estimation
**Monte Carlo SE at n=200:** ±1.54% → Cannot distinguish 94% from 96%
**Target:** 1000 iterations minimum (MC SE: ±0.69%)

**Action:**
```r
# Update Simulation_V2_200.R → Simulation_V2_1000.R
N_SIM <- 1000  # Was 200

# Add MC standard errors to output
mc_se_coverage <- sqrt(coverage * (1-coverage) / N_SIM)
```

**Time:** 1-2 days (mostly computation)

---

### 1.2 Add Comparison Methods

**Missing:** RVE, RoBMA (Bayesian model averaging)

**Action:**
```r
# Add to simulation:

# Robust Variance Estimation (clubSandwich package)
rve_result <- clubSandwich::coef_test(rma_fit, vcov = "CR2")

# RoBMA comparison (for subset - computationally expensive)
robma_fit <- RoBMA::RoBMA(yi, sei, seed = 42)
```

**Implementation notes:**
- RoBMA is slow (~30 sec per fit) - run on subset or separate
- RVE requires cluster structure - may not apply to all scenarios

**Time:** 3-5 days

---

### 1.3 Add Very Small k Scenarios

**Issue:** k=3,4 common in practice but not tested

**Action:**
```r
# Add scenarios to simulation:
list(name = "S_k3", k = 3, tau2 = 0.05, ...),
list(name = "S_k4", k = 4, tau2 = 0.05, ...),
```

**Expected issues:**
- LOO unstable for k=3 (33% removed)
- T-distribution very wide (t(1,0.975) = 12.7!)
- May need method modifications for very small k

**Time:** 1-2 days

---

### 1.4 Sensitivity Analysis for Parameters

**Issue:** Parameters appear ad hoc

**Parameters to test:**

| Parameter | Current | Test Range |
|-----------|---------|------------|
| stability_weight | 0.3 | 0.1, 0.2, 0.3, 0.4, 0.5 |
| SIT threshold | 2.5 | 2.0, 2.5, 3.0, 3.5 |
| SIT max_trim | 0.2 | 0.1, 0.2, 0.3 |
| UBSF blend | 60/40 | 70/30, 60/40, 50/50, 40/60 |

**Action:**
```r
# Grid search over parameter combinations
# Evaluate: RMSE, Coverage, Bias
# Report: Which parameters stable, which sensitive
```

**Time:** 2-3 days

---

### 1.5 Address Post-Selection Inference in SIT

**Issue:** Trimming based on residuals then using same data inflates Type I error

**Options:**

**Option A: Data Splitting**
- Split data: half for selection, half for inference
- Problem: Loses power with small k

**Option B: Selective Inference Theory**
- Use conditional inference (Lee et al., 2016)
- Complex to implement

**Option C: Bootstrap Correction**
```r
# Bootstrap the entire SIT procedure
for (b in 1:B) {
  boot_idx <- sample(1:k, replace = TRUE)
  sit_boot <- sit(yi[boot_idx], vi[boot_idx])
  estimates[b] <- sit_boot$estimate
}
se_corrected <- sd(estimates)
```

**Option D: Honest Acknowledgment**
- State limitation clearly
- Recommend SIT for exploration, not confirmatory inference

**Recommended:** Option C (bootstrap) + Option D (acknowledge limitation)

**Time:** 2-3 days

---

## PRIORITY 2: IMPORTANT (Should Do)

### 2.1 Add Large-k Scenarios

**Action:** Add k=50, k=100 scenarios
**Purpose:** Verify methods don't break at scale
**Expected:** Methods should converge to REML

**Time:** 0.5 days

---

### 2.2 Real Data Applications

**Current:** Only BCG (k=13, clean data)

**Needed:**
1. **Outlier example:** Find Cochrane MA with obvious outlier
2. **Publication bias example:** Asymmetric funnel plot
3. **Small k example:** k=4-5 studies
4. **Heterogeneous example:** High I²

**Sources:**
- Pairwise70 dataset (501 MAs available)
- metadat package

**Time:** 2-3 days

---

### 2.3 Theory/Limitations Section

**Write section acknowledging:**
1. Methods are heuristic, not theoretically optimal
2. No formal proof of coverage probability
3. Post-selection inference problem in SIT
4. Rubin's rules assume independence (violated in ARP)
5. Parameters tuned empirically

**Time:** 1 day

---

### 2.4 Method Selection Flowchart

```
                    Start
                      │
                      ▼
              Is k ≥ 5? ──No──► Use HKSJ only
                      │           (methods unstable)
                     Yes
                      │
                      ▼
         Suspect outliers? ──Yes──► SIT
                      │
                     No
                      │
                      ▼
         Suspect pub bias? ──Yes──► SIT
                      │
                     No
                      │
                      ▼
            Want robustness? ──Yes──► ARP
                      │
                     No
                      │
                      ▼
                    MWM
              (default choice)
```

**Time:** 0.5 days

---

## PRIORITY 3: DESIRABLE (Nice to Have)

### 3.1 Computational Benchmarking

```r
microbenchmark::microbenchmark(
  REML = rma(yi, vi),
  MWM = mafi_weighted_ma(yi, vi),
  SIT = sequential_influence_trimming(yi, vi),
  ARP = adaptive_robust_pooling(yi, vi),
  times = 100
)
```

### 3.2 Unit Tests

```r
# tests/testthat/test-advanced-pooling.R
test_that("MWM returns valid output", {
  result <- mafi_weighted_ma(c(0.5, 0.6, 0.4), c(0.1, 0.1, 0.1))
  expect_true(result$ci_lb < result$estimate)
  expect_true(result$estimate < result$ci_ub)
})
```

### 3.3 Include UBSF and EMA in Package

**Decision:** Include with warnings, or exclude with explanation
**Recommendation:** Exclude - not validated sufficiently

---

## SPECIFIC TECHNICAL FIXES

### Fix A: MWM - Handle scale() NA

```r
# Current (can produce NA):
combined_influence <- 0.5 * scale(influence)[, 1] + 0.5 * scale(abs(student_resid))[, 1]

# Fixed:
safe_scale <- function(x) {
  if (sd(x) == 0) return(rep(0, length(x)))
  (x - mean(x)) / sd(x)
}
combined_influence <- 0.5 * safe_scale(influence) + 0.5 * safe_scale(abs(student_resid))
```

### Fix B: SIT - Document Adaptive Threshold

```r
# Add to documentation:
# Adaptive threshold = base_threshold * (1 + 0.5 * (10 - min(k, 10)) / 10)
# For k=5: 2.5 * 1.25 = 3.125
# For k=10+: 2.5 * 1.0 = 2.5
# Rationale: More conservative for very small k
```

### Fix C: ARP - Handle Missing Estimators

```r
# Current: Falls through silently
# Fixed: Adjust variance for number of estimators
if (length(estimators) < 3) {
  warning(sprintf("Only %d estimators converged", length(estimators)))
  # Increase SE to account for reduced combining
  se <- se * sqrt(3 / length(estimators))
}
```

---

## REVISED SIMULATION PROTOCOL

```r
# Simulation_V2_1000_REVISED.R

# Parameters
N_SIM <- 1000
TRUE_EFFECT <- 0.3

# Extended scenarios
scenarios <- list(
  # Very small k (NEW)
  list(name = "k3_standard", k = 3, tau2 = 0.05),
  list(name = "k4_standard", k = 4, tau2 = 0.05),

  # Original scenarios
  list(name = "k5_standard", k = 5, tau2 = 0.05),
  list(name = "k10_standard", k = 10, tau2 = 0.05),
  list(name = "k20_standard", k = 20, tau2 = 0.05),

  # Large k (NEW)
  list(name = "k50_standard", k = 50, tau2 = 0.05),
  list(name = "k100_standard", k = 100, tau2 = 0.05),

  # Heterogeneity
  list(name = "k10_no_het", k = 10, tau2 = 0),
  list(name = "k10_high_het", k = 10, tau2 = 0.20),

  # Outliers
  list(name = "k10_outlier", k = 10, tau2 = 0.05, outlier = TRUE),
  list(name = "k15_outlier", k = 15, tau2 = 0.05, outlier = TRUE),

  # Publication bias
  list(name = "k15_pubbias_mild", k = 15, tau2 = 0.05, bias = 0.3),
  list(name = "k15_pubbias_mod", k = 15, tau2 = 0.05, bias = 0.5),
  list(name = "k15_pubbias_severe", k = 15, tau2 = 0.05, bias = 0.8)
)

# Methods to compare
methods <- c("REML", "HKSJ", "MWM", "SIT", "ARP", "RVE")  # Add RVE

# Output with MC standard errors
results <- results[, .(
  Bias = mean(est - true),
  Bias_MCSE = sd(est - true) / sqrt(.N),
  RMSE = sqrt(mean((est - true)^2)),
  Coverage = mean(ci_lb <= true & ci_ub >= true),
  Coverage_MCSE = sqrt(Coverage * (1 - Coverage) / .N)
), by = .(scenario, method)]
```

---

## TIMELINE

| Week | Tasks |
|------|-------|
| Week 1 | 1.1 (simulation), 1.3 (small k), 2.1 (large k) |
| Week 2 | 1.2 (comparisons), 1.4 (sensitivity), 1.5 (SIT fix) |
| Week 3 | 2.2 (real data), 2.3 (limitations), manuscript revision |

---

## SUCCESS CRITERIA

After revisions, we should be able to claim:

1. ✓ Coverage validated with 1000+ iterations (MC SE < 0.7%)
2. ✓ Methods compared to RVE (and ideally RoBMA)
3. ✓ Performance verified for k=3 to k=100
4. ✓ Sensitivity analysis shows robust parameters
5. ✓ Post-selection inference addressed (bootstrap + acknowledgment)
6. ✓ Real data examples demonstrate practical utility
7. ✓ Limitations clearly stated

---

**Next Step:** Start with Priority 1.1 (extend simulation to 1000 iterations)
