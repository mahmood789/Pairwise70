# Editorial Review: Advanced Pooling Methods for Meta-Analysis

**Journal:** Research Synthesis Methods
**Manuscript:** Novel Pooling Methods for Improved Coverage and Robustness in Random-Effects Meta-Analysis
**Decision:** MAJOR REVISIONS REQUIRED

---

## Editor's Summary

This manuscript introduces five novel pooling methods (MWM, SIT, ARP, UBSF, EMA) for random-effects meta-analysis, addressing the well-known under-coverage problem of standard REML estimation. The work demonstrates improved coverage (95-96% vs 92.6%) and robustness to outliers through simulation. While the contribution is potentially valuable, several methodological and presentation issues require attention before publication.

---

## STRENGTHS

### 1. Addresses Important Problem
- REML under-coverage is well-documented (IntHout et al., 2014; Hartung & Knapp, 2001)
- The gap between nominal 95% and actual ~92% coverage is clinically meaningful
- Outlier handling in meta-analysis remains an unsolved problem

### 2. Sound Core Ideas
- T-distribution CIs for small k is theoretically justified
- Studentized residuals superior to raw residuals for influence detection
- Combining multiple estimators (ARP) has precedent (Sidik & Jonkman, 2007)

### 3. Simulation Design
- Multiple scenarios covering realistic conditions
- Inclusion of outlier and publication bias scenarios
- Performance metrics appropriate (RMSE, coverage, bias)

### 4. Practical Implementation
- Methods integrated into R package
- Clear function interfaces
- Comparison function for applied researchers

---

## MAJOR CONCERNS

### 1. Insufficient Simulation Size (CRITICAL)

**Issue:** 200 iterations per scenario is inadequate for coverage estimation.

**Details:**
- Monte Carlo SE for 95% coverage at n=200: sqrt(0.95×0.05/200) = 1.54%
- Cannot distinguish 94% from 96% coverage reliably
- Standard practice: ≥1000 iterations (Burton et al., 2006; Morris et al., 2019)

**Required:** Increase to minimum 1000 iterations, preferably 5000 for publication-quality results.

### 2. Missing Comparisons to Existing Methods

**Issue:** No comparison to established robust methods.

**Missing comparisons:**
- HKSJ-modified estimators (Hartung-Knapp-Sidik-Jonkman)
- Robust variance estimation (RVE; Hedges et al., 2010)
- Bayesian model averaging (RoBMA; Bartoš et al., 2021)
- Permutation-based inference (Follmann & Proschan, 1999)
- Profile likelihood CIs

**Required:** Include at least HKSJ (already done), RVE, and one Bayesian comparator.

### 3. Limited Scenario Coverage

**Issue:** Simulation scenarios don't cover full range of realistic conditions.

**Missing scenarios:**
- Very small k (k=3, k=4) - common in practice
- Large k (k=50, k=100) - to check for unexpected behavior
- Varying n per study (not just 30-150)
- Non-normal true effects (skewed, heavy-tailed)
- Correlated effects (multi-arm trials)
- Missing data / selective reporting beyond simple publication bias

**Required:** Add at minimum k=3,4 scenarios and large-k validation.

### 4. Theoretical Justification Lacking

**Issue:** Methods are presented algorithmically without theoretical grounding.

**Questions unanswered:**
- What is the theoretical coverage probability of MWM?
- Under what conditions does SIT maintain valid inference after trimming?
- Is the 60/40 blend in UBSF principled or ad hoc?
- What are the asymptotic properties of these estimators?

**Required:** Either derive theoretical properties OR clearly state these are heuristic methods requiring empirical validation.

### 5. Arbitrary Tuning Parameters

**Issue:** Several key parameters appear chosen without justification.

| Parameter | Value | Justification Provided |
|-----------|-------|----------------------|
| stability_weight | 0.3 | "0.5 was too aggressive" |
| SIT threshold | 2.5 | None |
| SIT max_trim | 0.2 | None |
| UBSF blend | 60/40 | None |
| UBSF max bias_weight | 0.5 | None |

**Required:**
- Sensitivity analysis across parameter ranges
- Or principled selection criteria (e.g., cross-validation)
- Or acknowledgment as limitation

---

## MINOR CONCERNS

### 6. Presentation Issues

- Methods section mixes algorithm with rationale
- No flowchart for complex methods (UBSF, EMA)
- Inconsistent notation (τ² vs tau2)

### 7. Reproducibility

- Random seed provided (good)
- But no pre-registration of simulation protocol
- Results files saved but analysis code for figures not shown

### 8. Real Data Application

- Only BCG vaccine data shown
- No demonstration on challenging real-world examples
- No case study showing when methods differ substantively

### 9. Statistical Reporting

- P-values reported to 4 decimal places (over-precision)
- No confidence intervals on simulation performance metrics
- No formal hypothesis tests comparing methods

### 10. Software Considerations

- UBSF and EMA not included in final package (only MWM, SIT, ARP)
- No unit tests shown
- No benchmarking of computational time

---

## SPECIFIC TECHNICAL ISSUES

### A. MWM Method

1. **Leave-one-out instability:** For k=5, removing one study is 20% of data - LOO estimates may be unstable
2. **Weight normalization:** Combined influence score uses `scale()` which can produce NA for constant vectors
3. **Effective n:** Formula `1/sum(w²)` assumes independence - not appropriate for correlated weights

### B. SIT Method

1. **Post-selection inference:** Trimming based on residuals then using same data for inference inflates Type I error
2. **SE inflation:** The 0.1×n_trimmed/k adjustment is ad hoc - should derive from theory
3. **Adaptive threshold:** Formula not clearly specified - what is the functional form?

### C. ARP Method

1. **Estimator correlation:** REML, DL, PM are correlated (same data) - Rubin's rules assume independence
2. **Missing estimator handling:** What if PM fails? Falls back to fewer estimators without adjustment
3. **Equal prior weighting:** No justification for treating all estimators equally a priori

### D. UBSF Method (V2)

1. **PET-PEESE selection:** The p≥0.05 rule is known to be problematic (Stanley, 2017)
2. **Blend weights:** 60/40 split chosen to fix one failure case - may break others
3. **Conflict resolution:** "Use more conservative" may be biased toward null

---

## REQUIRED REVISIONS

### Essential (Must Address)

1. **Increase simulation to ≥1000 iterations**
2. **Add comparison to RVE and/or RoBMA**
3. **Include very small k scenarios (k=3,4)**
4. **Sensitivity analysis for tuning parameters**
5. **Address post-selection inference in SIT**
6. **Clarify that methods are heuristic, not theoretically optimal**

### Strongly Recommended

7. Add large-k scenarios (k=50+) to check behavior
8. Include real data examples beyond BCG
9. Add computational benchmarking
10. Provide confidence intervals on simulation metrics
11. Create method selection flowchart

### Desirable

12. Pre-register extended simulation protocol
13. Add unit tests to package
14. Include UBSF and EMA in package (or explain exclusion)
15. Formal theoretical analysis (if possible)

---

## REVIEWER RECOMMENDATIONS

### Reviewer 1 (Statistician)
"The core idea of combining stability weighting with robust estimation is sound, but the execution lacks rigor. The simulation is underpowered, and the post-selection inference problem in SIT is not addressed. Major revision required with theoretical clarification."

### Reviewer 2 (Applied Meta-Analyst)
"Practically useful contribution that addresses real pain points in meta-analysis. However, the guidance on when to use which method is insufficient. Need real data examples showing method selection in practice. The 200-iteration simulation is not convincing - I've run simulations and this is too few."

### Reviewer 3 (Methodologist)
"Interesting methods but the comparison landscape is incomplete. How do these compare to existing robust alternatives? The field has moved toward Bayesian approaches (RoBMA) - authors should position their frequentist contributions relative to this. Also concerned about parameter tuning - these look like they were tuned to the simulation scenarios."

---

## SUGGESTED REVISION TIMELINE

| Task | Estimated Effort |
|------|-----------------|
| Extend simulation to 1000 iterations | 1-2 days (computation) |
| Add RVE/RoBMA comparison | 3-5 days |
| Sensitivity analysis | 2-3 days |
| Small-k and large-k scenarios | 1-2 days |
| Real data applications | 2-3 days |
| Theory/limitations section | 2-3 days |
| Manuscript revision | 3-5 days |
| **Total** | **2-3 weeks** |

---

## EDITOR'S DECISION

**MAJOR REVISIONS**

The contribution addresses an important problem and shows promise, but requires substantial additional work before publication. The simulation evidence is insufficient, comparisons incomplete, and theoretical grounding absent.

If authors can address the essential revisions, particularly expanding the simulation and adding comparisons to established methods, this could be a valuable contribution to the meta-analysis methods literature.

The authors should pay particular attention to:
1. Honest assessment of when these methods fail
2. Clear guidance for applied researchers
3. Positioning relative to Bayesian alternatives

---

**Signed,**
*Editor, Research Synthesis Methods*

---

## AUTHOR RESPONSE TEMPLATE

### Response to Major Concerns

**1. Simulation Size**
- [ ] Extended to 1000 iterations
- [ ] Added Monte Carlo SEs to coverage estimates
- [ ] Results: [summarize any changes]

**2. Missing Comparisons**
- [ ] Added RVE comparison
- [ ] Added RoBMA comparison
- [ ] Results: [summarize findings]

**3. Limited Scenarios**
- [ ] Added k=3, k=4 scenarios
- [ ] Added k=50, k=100 scenarios
- [ ] Results: [summarize any edge cases]

**4. Theoretical Justification**
- [ ] Added limitations section acknowledging heuristic nature
- [ ] Or: Added theoretical derivations for [which methods]

**5. Tuning Parameters**
- [ ] Added sensitivity analysis (Table X)
- [ ] Results: Methods robust across [range] OR sensitive to [parameter]

### Response to Technical Issues

**A. MWM:** [response]
**B. SIT:** [response]
**C. ARP:** [response]
**D. UBSF:** [response]
