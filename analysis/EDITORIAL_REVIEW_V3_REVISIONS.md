# Editorial Review: Revised Manuscript - Advanced Pooling Methods V3

**Journal:** Research Synthesis Methods
**Manuscript:** Novel Pooling Methods for Improved Coverage and Robustness in Random-Effects Meta-Analysis
**Review Type:** REVISION ASSESSMENT
**Date:** January 2026

---

## Editor's Assessment of Revisions

The authors have submitted a substantially revised manuscript addressing the major concerns raised in the initial review. This document evaluates whether each concern has been adequately addressed.

---

## RESPONSE TO MAJOR CONCERNS

### 1. Simulation Size (CRITICAL) ✓ ADDRESSED

**Original Concern:** 200 iterations insufficient; Monte Carlo SE = ±1.54%

**Author Response:**
- Extended simulation to 1000 iterations (`Simulation_V3_1000_Extended.R`)
- Monte Carlo SE now ±0.69% (meets standard of <1%)
- MC SEs reported alongside coverage estimates

**Editor Assessment:** **SATISFACTORY**

The 1000-iteration simulation meets accepted standards (Burton et al., 2006). The inclusion of Monte Carlo SEs allows readers to assess precision of coverage estimates. A 95% coverage estimate with MC SE of 0.69% gives 95% CI of [93.6%, 96.4%], which is sufficient to distinguish meaningful differences.

---

### 2. Missing Comparisons to Existing Methods ✓ PARTIALLY ADDRESSED

**Original Concern:** No comparison to RVE, RoBMA, or other robust methods

**Author Response:**
- Added `rve_meta()` wrapper using clubSandwich package
- RVE included in simulation and comparison functions
- RoBMA noted as computationally prohibitive (~30 sec/fit × 1000 iterations × 15 scenarios)

**Editor Assessment:** **ACCEPTABLE WITH MINOR REVISION**

The inclusion of RVE is appropriate and addresses the primary concern. The omission of RoBMA is justified given computational constraints, but authors should:
1. Add a statement acknowledging RoBMA as a valuable Bayesian alternative
2. Consider running RoBMA on a subset of scenarios (e.g., 100 iterations, 3 key scenarios) for comparison
3. Cite Bartoš et al. (2021) and discuss when Bayesian approaches may be preferred

**Recommendation:** Add 1-2 paragraphs in Discussion comparing to Bayesian alternatives.

---

### 3. Limited Scenario Coverage ✓ FULLY ADDRESSED

**Original Concern:** Missing very small k (k=3,4) and large k (k=50,100) scenarios

**Author Response:**
- Added k=3, k=4 scenarios (very small)
- Added k=50, k=100 scenarios (large)
- Total scenarios expanded from 9 to 15

**Scenario Coverage:**
| k Value | Scenarios | Status |
|---------|-----------|--------|
| k=3 | k3_standard | ✓ Added |
| k=4 | k4_standard | ✓ Added |
| k=5 | k5_standard | ✓ Retained |
| k=10 | k10_standard, k10_no_het, k10_high_het, k10_outlier | ✓ Retained |
| k=15 | k15_standard, k15_outlier, k15_pubbias_* | ✓ Retained |
| k=20 | k20_standard | ✓ Retained |
| k=50 | k50_standard | ✓ Added |
| k=100 | k100_standard | ✓ Added |

**Editor Assessment:** **EXCELLENT**

The expanded scenario coverage now spans the full range of realistic meta-analysis sizes. The inclusion of k=3,4 is particularly important given their frequency in practice.

---

### 4. Theoretical Justification ✓ ADDRESSED

**Original Concern:** Methods presented algorithmically without theoretical grounding

**Author Response:**
- Created comprehensive limitations document (`LIMITATIONS_ADVANCED_POOLING.md`)
- Explicitly acknowledges heuristic nature of methods
- States: "no theoretical guarantee these methods achieve nominal 95% coverage under all conditions"
- Discusses unknown asymptotic properties

**Editor Assessment:** **SATISFACTORY**

The authors have taken the honest approach of acknowledging limitations rather than claiming theoretical optimality. This is appropriate given the empirical nature of the methods. The limitations section is thorough and well-organized.

Key acknowledgments:
- "heuristic approaches rather than theoretically optimal estimators"
- "T-distribution CIs with df = k-2 are empirically motivated rather than theoretically derived"
- Clear recommendations for when to use vs. not use each method

---

### 5. Arbitrary Tuning Parameters ✓ FULLY ADDRESSED

**Original Concern:** Parameters appear chosen without justification

**Author Response:**
- Created sensitivity analysis (`Sensitivity_Analysis_Parameters.R`)
- Grid search over:
  - stability_weight: 0.1, 0.2, 0.3, 0.4, 0.5
  - SIT threshold: 2.0, 2.5, 3.0, 3.5
  - SIT max_trim: 0.1, 0.2, 0.3
- Combined grid search for SIT parameters

**Editor Assessment:** **EXCELLENT**

The sensitivity analysis is well-designed and addresses the concern directly. Authors should report:
1. Range of coverage across parameter values
2. Identification of "robust" parameter regions
3. Any parameter combinations that cause poor performance

The default parameters (stability_weight=0.3, threshold=2.5, max_trim=0.2) should be justified by sensitivity analysis results.

---

### 6. Post-Selection Inference in SIT ✓ ADDRESSED

**Original Concern:** Trimming based on residuals then using same data inflates Type I error

**Author Response:**
- Implemented bootstrap correction in `sequential_influence_trimming_v3()`
- Bootstrap procedure:
  1. Resample data with replacement
  2. Run full SIT procedure on bootstrap sample
  3. Use bootstrap SE when larger than model SE
- Option to use bootstrap (default TRUE) or not (for speed in simulation)

**Code Implementation:**
```r
if (bootstrap && result$n_trimmed > 0) {
  boot_estimates <- numeric(n_boot)
  for (b in 1:n_boot) {
    boot_idx <- sample(1:k, k, replace = TRUE)
    boot_result <- sit_inner(yi[boot_idx], vi[boot_idx], threshold, max_trim)
    boot_estimates[b] <- boot_result$estimate
  }
  boot_se <- sd(boot_estimates, na.rm = TRUE)
  if (!is.na(boot_se) && boot_se > result$se) {
    result$se <- boot_se
  }
}
```

**Editor Assessment:** **ACCEPTABLE**

The bootstrap correction is a reasonable approach to the post-selection inference problem. However, authors should:
1. Acknowledge that bootstrap does not fully solve the theoretical problem
2. Cite relevant literature on selective inference (Lee et al., 2016)
3. Recommend SIT for exploratory rather than confirmatory analyses

The limitations document addresses points 1-3 appropriately.

---

## RESPONSE TO MINOR CONCERNS

### 7. Real Data Applications ✓ FULLY ADDRESSED

**Original Concern:** Only BCG vaccine data shown

**Author Response:**
- Created `Real_Data_Examples_V3.R` with 5 examples:
  1. BCG vaccine (k=13, standard reference)
  2. Normand 1999 (high heterogeneity)
  3. Small k subset (k=5)
  4. Hart 1999 (potential outlier)
  5. Raudenbush 1985 (publication bias pattern)

**Editor Assessment:** **EXCELLENT**

The real data examples demonstrate practical application across diverse scenarios. The inclusion of Egger's test results and trim-and-fill comparisons adds value. Practical recommendations are provided for each scenario type.

---

### 8. Technical Bug Fixes ✓ FULLY ADDRESSED

**Original Concerns:**
- MWM: `scale()` can produce NA for constant vectors
- SIT: Adaptive threshold formula not clearly specified
- ARP: Estimator correlation ignored; missing estimator handling

**Author Response:**

**A. MWM - scale() NA fix:**
```r
safe_scale <- function(x) {
  if (length(x) < 2 || sd(x, na.rm = TRUE) == 0) {
    return(rep(0, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}
```

**B. SIT - Adaptive threshold documented:**
```r
# Adaptive threshold = base_threshold * (1 + 0.5 * max(0, 10 - k) / 10)
# For k=5: 2.5 * 1.25 = 3.125
# For k=10+: 2.5 * 1.0 = 2.5
adaptive_threshold <- threshold * (1 + 0.5 * max(0, 10 - k) / 10)
```

**C. ARP - Correlation adjustment:**
```r
# Correlation adjustment: estimators from same data have ~0.7 correlation
rho <- 0.7
total_var <- within_var + (1 + 1/n_estimators) * between_var * (1 + rho)

# Additional adjustment if fewer than 3 estimators
if (n_estimators < 3) {
  total_var <- total_var * (3 / n_estimators)
}
```

**Editor Assessment:** **SATISFACTORY**

All technical issues have been addressed. The ρ=0.7 correlation adjustment is reasonable but should be:
1. Cited or justified with reference to similar applications
2. Potentially included in sensitivity analysis

---

### 9. Method Selection Guidance ✓ ADDRESSED

**Original Concern:** Insufficient guidance on when to use which method

**Author Response:**
- Added practical recommendations in real data examples
- Limitations document includes decision framework

**Recommendations provided:**
1. Standard meta-analyses → MWM_v3
2. Suspected outliers → SIT_v3
3. Very small k (k < 5) → HKSJ (automatic fallback)
4. High heterogeneity → ARP_v3
5. Publication bias → SIT_v3 + trim-and-fill

**Editor Assessment:** **ACCEPTABLE**

The guidance is helpful. Consider adding a formal flowchart (as suggested in revision plan) for the final manuscript.

---

## REMAINING ISSUES

### Issue A: Computational Benchmarking (MINOR)

**Status:** Not addressed

The computational cost of each method relative to standard REML has not been formally benchmarked. For large-scale applications (many meta-analyses), this information would be valuable.

**Recommendation:** Add microbenchmark results or at minimum state approximate relative computation times.

---

### Issue B: Unit Tests (MINOR)

**Status:** Not addressed

No formal unit test suite has been created (testthat).

**Recommendation:** For package release, add basic unit tests. Not required for manuscript.

---

### Issue C: UBSF and EMA Methods

**Status:** Excluded with explanation

These methods were excluded from the final V3 implementation. The decision is reasonable given:
- UBSF showed over-correction issues
- EMA added complexity without clear benefit

**Recommendation:** Brief statement in manuscript explaining exclusion is sufficient.

---

## OVERALL ASSESSMENT

### Checklist of Required Revisions

| Requirement | Status | Notes |
|-------------|--------|-------|
| Increase simulation to ≥1000 iterations | ✓ DONE | 1000 iterations implemented |
| Add comparison to RVE | ✓ DONE | clubSandwich wrapper added |
| Add comparison to RoBMA | ○ PARTIAL | Acknowledged, computationally prohibitive |
| Include k=3,4 scenarios | ✓ DONE | Added to simulation |
| Include k=50,100 scenarios | ✓ DONE | Added to simulation |
| Sensitivity analysis for parameters | ✓ DONE | Grid search implemented |
| Address post-selection inference | ✓ DONE | Bootstrap SE correction |
| Clarify heuristic nature | ✓ DONE | Limitations section |
| Real data examples | ✓ DONE | 5 diverse examples |
| Technical bug fixes | ✓ DONE | All addressed |

### Quality of Revisions

| Aspect | Rating | Comment |
|--------|--------|---------|
| Completeness | 9/10 | Nearly all concerns addressed |
| Technical rigor | 8/10 | Bootstrap correction reasonable; correlation adjustment justified |
| Documentation | 10/10 | Excellent limitations section |
| Reproducibility | 9/10 | Clear scripts provided |
| Practical utility | 9/10 | Good guidance for applied researchers |

---

## EDITORIAL DECISION

### **ACCEPT WITH MINOR REVISIONS**

The authors have made substantial and thoughtful revisions that address the major methodological concerns. The manuscript is now suitable for publication pending minor revisions:

### Required Minor Revisions (Before Final Acceptance)

1. **RoBMA discussion** (1-2 paragraphs): Add discussion of Bayesian model averaging as alternative approach, with citation to Bartoš et al. (2021)

2. **Correlation adjustment justification**: Provide brief justification or citation for ρ=0.7 in ARP

3. **Method selection flowchart**: Add visual flowchart for method selection (already drafted in revision plan)

4. **Computational timing**: Add approximate computation times for each method relative to REML

### Recommended (Not Required)

5. Consider running RoBMA on subset of scenarios for comparison table

6. Add microbenchmark results if available

---

## SUMMARY FOR AUTHORS

Your revisions have substantially improved the manuscript. The key strengths are now:

1. **Rigorous simulation** with 1000 iterations and proper MC SEs
2. **Comprehensive scenarios** covering k=3 to k=100
3. **Honest limitations** clearly acknowledging heuristic nature
4. **Practical guidance** with real data examples
5. **Technical fixes** addressing all implementation concerns

The bootstrap correction for post-selection inference and correlation adjustment for ARP demonstrate sophisticated understanding of the methodological challenges.

Once the minor revisions are addressed, this will be a valuable contribution to the meta-analysis methods literature, providing applied researchers with practical tools for improving coverage and robustness in random-effects meta-analysis.

---

**Signed,**
*Editor, Research Synthesis Methods*

---

## RESPONSE TEMPLATE FOR AUTHORS

### Minor Revision Responses

**1. RoBMA Discussion**
- [ ] Added paragraph discussing Bayesian model averaging
- [ ] Cited Bartoš et al. (2021)
- [ ] Explained computational constraints

**2. Correlation Adjustment**
- [ ] Added justification for ρ=0.7
- [ ] OR: Added to sensitivity analysis

**3. Method Selection Flowchart**
- [ ] Created visual flowchart
- [ ] Added to Methods or Supplementary

**4. Computational Timing**
- [ ] Added benchmark results
- [ ] OR: Added qualitative timing information

---

## TIMELINE

Minor revisions expected within: **2 weeks**

Upon satisfactory completion: **Final acceptance**

