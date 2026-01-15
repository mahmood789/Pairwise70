# Editorial Review: Research Synthesis Methods

## Manuscript: "The Information Size Crisis in Meta-Analysis"

**Editor:** Anonymous (Simulated Review)
**Date:** January 2026
**Decision:** MAJOR REVISIONS REQUIRED

---

## Overall Assessment

This manuscript addresses an important and genuinely overlooked topic in meta-analysis methodology. The application of Trial Sequential Analysis concepts at scale to 492 Cochrane meta-analyses represents a substantial contribution. However, several methodological concerns must be addressed before publication can be recommended.

**Recommendation:** Major Revision

---

## MAJOR ISSUES (Must Address)

### 1. Arbitrary Choice of Effect Size for OIS Calculation

**Problem:** The OIS calculation assumes a 20% Relative Risk Reduction (RRR) for all binary outcomes and SMD = 0.3 for continuous outcomes. This is a critical assumption that drives the entire analysis.

**Concern:**
- Different clinical contexts have vastly different minimal clinically important differences (MCIDs)
- A 20% RRR for mortality is very different from 20% RRR for symptom improvement
- This assumption is not justified empirically or clinically

**Required Action:**
1. Conduct sensitivity analysis with RRR = 10%, 15%, 25%, 30%
2. Consider using the observed effect size as the target (what information would be needed to detect *this* effect?)
3. Alternatively, use clinical thresholds from GRADE or similar frameworks
4. Report how results change across different assumed effect sizes

### 2. IAI Component Weights Lack Justification

**Problem:** The IAI uses fixed weights (40% IF, 30% HAP, 20% SBCS, 10% CZS) without empirical or theoretical justification.

**Concern:**
- This mirrors the criticism your MAFI work addressed through sensitivity analysis
- Why is Information Fraction weighted twice as much as Sequential Status?
- No cross-validation or data-driven weight optimization presented

**Required Action:**
1. Provide theoretical justification for weight choices
2. Conduct weight sensitivity analysis (±10%, ±20% perturbations)
3. Consider empirical weight derivation using predictive modeling
4. Report correlation between components (are they redundant?)

### 3. Sequential Boundary Implementation Oversimplified

**Problem:** The O'Brien-Fleming boundary implementation is a simplified approximation, not the proper alpha-spending function.

**From Code:**
```r
if (information_fraction >= 1) {
  boundary <- qnorm(1 - alpha/2)
} else if (information_fraction >= 0.5) {
  boundary <- qnorm(1 - alpha/2) / sqrt(information_fraction)
}
```

**Concern:**
- True O'Brien-Fleming uses: `Z = c * Φ⁻¹(1 - α/2)` where c depends on number of looks
- The implementation doesn't properly account for discrete number of interim analyses
- Futility boundary of 0.5 × z_α/2 is arbitrary

**Required Action:**
1. Use proper alpha-spending function (Lan-DeMets)
2. Or explicitly state this is an approximation and validate against TSA software
3. Provide comparison with results from Copenhagen Trial Unit's TSA software on subset

### 4. Heterogeneity-Diversity Confusion

**Problem:** The OIS adjustment formula uses I² directly, but the correct adjustment uses D² (diversity), not I².

**Technical Issue:**
```
OIS = D × (1 + I²/(1-I²))  [Current implementation]
```

Should consider:
```
DARIS = D × Diversity = D × (1 + D²/(1-D²))
```

Where D² ≠ I² when studies have unequal sizes.

**Required Action:**
1. Clarify whether I² or D² is appropriate for your context
2. Cite Wetterslev et al. (2009) on DARIS calculation
3. Consider calculating proper diversity measure

### 5. Unit of Analysis Problem

**Problem:** Analysis is at the dataset level (n=492), but datasets contain multiple comparisons/outcomes. The proposal mentions 4,424 meta-analyses but only 492 were analyzed.

**Concern:**
- Why the discrepancy?
- Are you analyzing only primary comparisons?
- How were datasets with multiple outcomes handled?
- This needs clarification

**Required Action:**
1. Explain selection of 492 datasets from 501 available
2. Clarify handling of multiple comparisons within reviews
3. Consider expanding to all 4,424 comparison-outcome pairs

---

## MODERATE ISSUES

### 6. Missing Comparison with Existing TSA Literature

**Concern:** The manuscript doesn't compare findings with previous TSA meta-research:
- Imberger et al. (2016) found 14% false positives in Cochrane reviews
- Castellini et al. (2018) systematic review of TSA applications
- How do your 24.3% "premature conclusions" compare?

**Required Action:**
Add comparison with published TSA findings; discuss concordance/discordance.

### 7. Confidence Intervals Missing

**Key Finding:** "24.3% of significant meta-analyses have inadequate information"

**Problem:** No confidence interval or uncertainty quantification provided.

**Required Action:**
- Bootstrap confidence intervals for key proportions
- Consider Bayesian credible intervals

### 8. Counterintuitive Heterogeneity Finding

**Result:** Low I² (0-25%) shows LOWER adequacy (44.4%) than high I² (>75%, 62.0%)

**This is opposite to expectation.** High heterogeneity should require MORE information, making adequacy harder to achieve.

**Possible Explanation:**
- Selection bias? High-I² reviews may have more studies/larger samples
- Check if k or N differs systematically by I² category

**Required Action:**
Investigate and explain this finding; it may indicate a flaw in the methodology.

### 9. "Critical" Category Definition

**Problem:** IAI < 0.25 labeled "Critical" but no validation that this threshold is meaningful.

**Required Action:**
- What happens to conclusions in "Critical" vs "Adequate" meta-analyses?
- External validation needed (e.g., do Critical MAs reverse on update more often?)

---

## MINOR ISSUES

### 10. Continuous Outcomes Underrepresented
Only 24/492 (4.9%) are continuous. Consider whether findings generalize.

### 11. Code Errors in Output
Warning messages visible in output (`..count..` deprecated). Clean for final version.

### 12. Missing MAFI Comparison
Stated objective but "implement comparison after data alignment" appears in output. This should be completed.

### 13. Visualization Improvements Needed
- Fig 5 uses deprecated ggplot2 syntax
- Consider adding forest plot of most/least adequate MAs

### 14. Reproducibility
- Random seed not set
- Package versions not documented
- Consider R Markdown for reproducibility

---

## STRENGTHS (Acknowledge in Revision)

1. **Novel Contribution:** First large-scale empirical assessment of information adequacy
2. **Important Finding:** 24.3% premature conclusion rate is clinically concerning
3. **Clear Presentation:** IAI framework is intuitive and potentially useful
4. **Strong Data:** 492 real Cochrane meta-analyses provide credibility
5. **k-Crisis Pattern:** Clear demonstration that k < 10 is problematic

---

## SPECIFIC QUESTIONS FOR AUTHORS

1. Why was 20% RRR chosen? Is this evidence-based or arbitrary?
2. How does IAI perform when the true effect is null?
3. Have you validated against the Copenhagen TSA software?
4. What is the correlation between IAI and MAFI?
5. How should practitioners use IAI in combination with other metrics?

---

## REQUIRED REVISIONS SUMMARY

| Priority | Issue | Action Required |
|----------|-------|-----------------|
| Critical | RRR assumption | Sensitivity analysis across effect sizes |
| Critical | IAI weights | Justify or optimize empirically |
| Critical | Sequential boundaries | Validate against proper TSA |
| Major | I² vs D² | Clarify diversity adjustment |
| Major | 492 vs 4,424 | Explain unit of analysis |
| Major | Heterogeneity paradox | Investigate and explain |
| Moderate | CI for proportions | Add uncertainty quantification |
| Moderate | TSA literature comparison | Compare with Imberger et al. |
| Minor | MAFI comparison | Complete as stated |
| Minor | Code cleanup | Remove warnings |

---

## DECISION

**Major Revision Required**

The core contribution is valuable and addresses a real gap in the literature. However, the methodological concerns—particularly around the assumed effect size and weight justification—must be addressed before this work can be recommended for publication in Research Synthesis Methods.

The counterintuitive finding regarding heterogeneity requires investigation, as it may indicate a fundamental issue with the methodology.

I look forward to reviewing a revised version that addresses these concerns.

---

**Editor**
Research Synthesis Methods (Simulated)

*This review is intended to strengthen the manuscript before submission.*
