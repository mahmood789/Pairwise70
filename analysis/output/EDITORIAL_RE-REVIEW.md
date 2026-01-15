# Editorial Re-Review: MAFI and Research Synthesis Modeling Framework

**Journal:** Research Synthesis Methods
**Manuscript:** "MAFI: A Novel Multi-Dimensional Fragility Index for Meta-Analysis"
**Review Type:** Re-Review Following Major Revisions
**Date:** January 2026

---

## OVERALL ASSESSMENT

**Previous Recommendation:** Major Revisions Required
**Current Recommendation:** Minor Revisions Required

The authors have made substantial efforts to address the editorial concerns. The weight sensitivity analysis, review-level cross-validation, and mixed-effects modeling are particularly commendable. However, several issues emerged from the revisions that require attention before acceptance.

**Revised Scores:**
| Criterion | Original | Revised | Change |
|-----------|----------|---------|--------|
| Significance | 7/10 | 8/10 | +1 |
| Methodological Rigor | 6/10 | 8/10 | +2 |
| Presentation Quality | 5/10 | 7/10 | +2 |
| **Overall** | **6/10** | **7.7/10** | **+1.7** |

---

## RESPONSE TO MAJOR CONCERNS

### 1. Weight Justification - SATISFACTORILY ADDRESSED

**Original Concern:** Arbitrary weight selection (30/25/20/15/10%)

**Response Quality:** Excellent

The sensitivity analysis is thorough and convincing:
- All ±10% perturbations yield correlations >0.996 with base MAFI
- Maximum classification change: 4.8% (Effect Stability -10%)
- Mean classification change: 3.1%

This demonstrates MAFI is robust to reasonable weight variations. The conclusion that weights are not driving the results is well-supported.

**Remaining Issue:** The empirical weights derived from logistic regression reveal a concern:

| Component | Original Weight | Empirical Weight |
|-----------|-----------------|------------------|
| Direction (DFI) | 0.30 | 0.171 |
| Significance (SFI) | 0.25 | 0.303 |
| Clinical (CFI) | 0.20 | **0.001** |
| Effect Stability | 0.15 | **0.001** |
| CI Stability | 0.10 | **0.523** |

**Critical Finding:** Clinical (CFI) and Effect Stability contribute essentially zero predictive value empirically. This raises questions:
1. Should these components be retained in MAFI?
2. Is the 0.89 correlation between original and empirical MAFI acceptable given these discrepancies?

**Required Action:** Discuss whether a 3-component MAFI (DFI + SFI + CI Stability) might be more parsimonious and equally valid.

---

### 2. Penalty Parameters - PARTIALLY ADDRESSED

**Original Concern:** Arbitrary penalty thresholds

**Response Quality:** Good

The empirical derivation is appreciated:
- Heterogeneity penalty: (I²/100) × 0.17 (derived from logistic coefficient)
- k threshold: Fragility stabilizes around k=20-30 (empirically shown)

The fragility-by-k analysis is particularly compelling:

| k Range | Fragility Rate |
|---------|---------------|
| 2-5 | 49.2% |
| 10-15 | 41.9% |
| 30-50 | 25.4% |
| >100 | 5.9% |

**Remaining Issue:** The empirical k threshold reported as "k=3" appears to be a computational artifact (first point where derivative < 0.005). The data clearly show fragility continues declining until k≈50-100. Recommend revising to "k=30-50" as the practical threshold.

---

### 3. Circular Validation - SATISFACTORILY ADDRESSED

**Original Concern:** MAFI validated on same data used to derive it

**Response Quality:** Excellent

Review-level cross-validation properly addresses the clustering concern:
- 10-fold CV at systematic review level (473 reviews)
- AUC: 0.687 (SD: 0.034)
- Consistent with MA-level CV (0.688)

This demonstrates the model generalizes across independent reviews, not just independent MAs. The similar AUC values suggest minimal optimism bias from within-review correlation.

**No Further Action Required**

---

### 4. Leave-One-Out Assumption - APPROPRIATELY ACKNOWLEDGED

**Original Concern:** Equal study importance assumption

**Response Quality:** Good

The limitation is honestly acknowledged:
- Study-level sample sizes unavailable in aggregated data
- Risk of bias not consistently available
- Weighted LOO would require individual participant data

The framing as "study-count-based" vs "precision-weighted" fragility is appropriate.

**Minor Issue:** The SE quartile analysis shows a non-monotonic relationship (highest fragility in Q2, not Q4). This is unexpected and deserves brief discussion.

---

### 5. Clustering - SATISFACTORILY ADDRESSED

**Original Concern:** MAs nested within reviews not properly handled

**Response Quality:** Excellent

The mixed-effects analysis is comprehensive:
- **ICC = 0.161** (16.1% of variance between reviews)
- Random intercept model converged successfully
- Fixed effects consistent between mixed and standard GLM

| Variable | Mixed Effect (SE) | Fixed Effect (SE) |
|----------|------------------|-------------------|
| log_k | -0.681 (0.048) | -0.623 (0.039) |
| I² | 1.977 (0.153) | 1.887 (0.135) |
| τ² | 0.050 (0.014) | 0.043 (0.013) |
| |Effect| | -1.436 (0.116) | -1.348 (0.105) |

The coefficient inflation in mixed model (5-10%) is expected and appropriate. Standard errors are correctly larger in the mixed model.

**No Further Action Required**

---

### 6. Missing Data - SATISFACTORILY DOCUMENTED

**Original Concern:** 108 exclusions unexplained

**Response Quality:** Good

Documentation is clear:
- 108 MAs (2.4%) missing fragility calculations
- Missingness not associated with k (p = 0.095)
- Missing cases have slightly larger k (21.7 vs 18.1) and effect sizes (0.46 vs 0.34)

**Minor Issue:** The p = 0.095 for k is borderline. The effect size difference (0.46 vs 0.34) was not tested. Recommend complete-case sensitivity analysis or brief discussion of potential bias direction.

---

### 7. Atal Comparison - APPROPRIATELY CLARIFIED

**Original Concern:** Incomplete comparison with Atal et al. (2019)

**Response Quality:** Good

The conceptual distinction is well-articulated:
- **Atal:** Manipulation-based (event reassignment within studies)
- **MAFI:** Exclusion-based (study removal from meta-analysis)

These measure complementary aspects of fragility:
- Atal: Robustness to data errors
- MAFI: Robustness to study selection

The 31.3% significance-fragile (MAFI) vs 46% FI≤5 (Atal) difference is expected given different constructs.

**Minor Issue:** A formal correlation between SFI and traditional FI on the same dataset would strengthen this section, but is not required.

---

### 8. Simplification - PROVIDED AS REQUESTED

**Original Concern:** Consider simpler alternative

**Response Quality:** Good

MAFI-Simple (DFI + SFI only, 50/50 weights) is offered:
- Correlation with full MAFI: 0.889
- Classification agreement: 66.5%

**Issue:** The 66.5% agreement is lower than expected. One-third of MAs would be classified differently. This needs discussion - is MAFI-Simple an acceptable substitute or a fundamentally different measure?

---

## NEW CONCERNS FROM REVISIONS

### A. Empirical Weight Discrepancy (MODERATE)

The empirical weights suggest a radically different structure:
- CI Stability dominates (52% vs 10% original)
- Clinical and Effect Stability contribute nothing (0.1% each)

This wasn't adequately addressed. Options:
1. Acknowledge and justify retaining original weights on theoretical grounds
2. Propose empirical MAFI as the preferred version
3. Present both and let users choose

### B. Convergence Warnings (MINOR)

The revision script produced:
```
Warning: glm.fit: algorithm did not converge
Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
```

These indicate separation issues in the logistic models used for empirical weight derivation. The empirical weights may be unstable. Consider penalized regression (Firth) or note this limitation.

### C. MAFI-Simple Agreement (MODERATE)

66.5% classification agreement between full and simple MAFI is concerning. If users are offered both:
- Which should be primary?
- When would disagreement matter?
- Could disagreement itself be informative?

---

## REVISED CHECKLIST

| Item | Status | Notes |
|------|--------|-------|
| Weight sensitivity | ✓ Complete | Robust to ±10% |
| Empirical weights | ⚠ Discuss | CFI/Effect near zero |
| Penalty derivation | ✓ Complete | Empirical basis provided |
| k threshold | ⚠ Revise | "k=3" appears erroneous |
| Review-level CV | ✓ Complete | AUC 0.687 stable |
| LOO limitation | ✓ Acknowledged | Appropriately framed |
| ICC reported | ✓ Complete | 16.1% between-review |
| Mixed model | ✓ Complete | Coefficients consistent |
| Missing data | ✓ Documented | 2.4%, appears random |
| Atal comparison | ✓ Complete | Complementary measures |
| MAFI-Simple | ⚠ Discuss | 66.5% agreement low |

---

## SPECIFIC REVISIONS REQUIRED

### Must Address (Before Acceptance)

1. **Discuss empirical weight discrepancy** (Section 1.2):
   - Why do CFI and Effect Stability contribute nothing empirically?
   - Should original weights be retained on theoretical grounds?
   - Consider removing or downweighting these components

2. **Correct k threshold statement**:
   - Change "k=3" to empirically supported "k≈30-50"
   - Or clarify the algorithm used

3. **Address MAFI-Simple agreement**:
   - Discuss implications of 34% classification disagreement
   - Provide guidance on when to use each version

### Should Address (Recommended)

4. **Note convergence issues** in empirical weight derivation
5. **Test effect size difference** in missing data analysis
6. **Discuss non-monotonic SE-fragility relationship**

### Could Address (Optional)

7. Formal correlation with Atal FI on overlapping subset
8. Calibration plots for predictive models
9. Case studies demonstrating MAFI application

---

## QUESTIONS FOR AUTHORS

1. Given that CFI and Effect Stability contribute <1% empirically, what is the justification for their 35% combined weight in MAFI?

2. The k=3 threshold for fragility stabilization contradicts the empirical data showing fragility declining through k=100. Can you clarify this discrepancy?

3. Would you recommend MAFI-Full or MAFI-Simple for routine use? Under what circumstances would each be preferred?

4. The mixed-effects model shows 16% between-review variance. Does this suggest review-level factors (topic, methodology, team) systematically affect fragility?

5. Have you considered a weighted average of DFI, SFI, and CI Stability only, given the empirical weight analysis?

---

## RECOMMENDATION SUMMARY

| Aspect | Original | After Revision |
|--------|----------|----------------|
| Weight justification | Poor | Good (with caveat) |
| Penalty parameters | Poor | Good |
| Validation approach | Moderate | Excellent |
| Clustering handling | Poor | Excellent |
| Missing data | Poor | Good |
| Comparisons | Moderate | Good |
| Transparency | Moderate | Very Good |

**Decision:** Minor Revisions

The authors have made substantial improvements. The remaining issues are clarifications rather than fundamental methodological problems. Once the discrepancy between original and empirical weights is addressed (either by justification or revision), the manuscript will be suitable for publication.

---

## SUGGESTED TIMELINE

**Required revisions:** 2-4 weeks
- Address empirical weight discrepancy
- Correct k threshold
- Discuss MAFI-Simple agreement

**Review turnaround:** 1 week after resubmission

---

*Reviewed by: [Editor, Research Synthesis Methods]*
*Conflict of Interest: None declared*

---

## APPENDIX: SUMMARY STATISTICS FROM REVISIONS

### Weight Sensitivity Analysis
- 10 perturbation scenarios tested
- All correlations with base MAFI: 0.996-0.999
- Classification changes: 2.3%-4.8%
- **Conclusion:** Robust

### Review-Level Cross-Validation
- 10 folds, 47-48 reviews per fold
- AUC range: 0.634-0.731
- Mean AUC: 0.687 (SD: 0.034)
- **Conclusion:** Stable, no overfitting

### Mixed-Effects Model
- ICC: 0.161
- All fixed effects significant (p < 0.001)
- Coefficient changes vs GLM: 5-10% (expected)
- **Conclusion:** Clustering appropriately handled

### Fragility by Study Count
| k | N | Fragility Rate |
|---|---|---------------|
| 2-5 | 1,409 | 49.2% |
| 6-10 | 1,124 | 44.4% |
| 11-15 | 580 | 41.9% |
| 16-20 | 315 | 34.3% |
| 21-30 | 354 | 33.1% |
| 31-50 | 252 | 25.4% |
| 51-100 | 181 | 17.1% |
| >100 | 101 | 5.9% |

**Key Finding:** Fragility decreases monotonically with k; stabilizes around k=50-100.
