# Manuscript Outline
## Predicting Meta-Analysis Fragility: A Meta-Epidemiological Study

**Target Journal:** Research Synthesis Methods
**Word Limit:** ~5,000 words (check journal guidelines)

---

## TITLE PAGE

**Title:** Predicting Meta-Analysis Fragility: A Meta-Epidemiological Study of 3,556 Cochrane Intervention Reviews

**Running head:** Predicting Meta-Analysis Fragility

**Authors:** [Names and affiliations]

**Corresponding author:** [Contact details]

**Word count:** [Abstract: 248; Main text: ~4,500]

**Tables:** 3 (main) + 3 (supplementary)
**Figures:** 4 (main) + 8 (supplementary)

---

## ABSTRACT (248 words)
See SUBMISSION_ELEMENTS.md

---

## 1. INTRODUCTION (~600 words)

### 1.1 The Problem of Evidence Fragility
- Meta-analyses are cornerstone of evidence-based medicine
- Conclusions can be sensitive to removal of individual studies
- "Fragility" = vulnerability to single-study exclusion
- Current practice: fragility assessed post-hoc, if at all

### 1.2 Existing Approaches
- Fragility Index (Walsh et al., 2014) - limited to binary outcomes
- Influence diagnostics (Viechtbauer & Cheung, 2010)
- Leave-one-out analyses - descriptive, not predictive
- Gap: no prospective prediction tools

### 1.3 Study Objectives
1. Characterize fragility prevalence across Cochrane meta-analyses
2. Identify predictors of fragility
3. Develop and validate a predictive model
4. Assess publication bias and its relationship to fragility
5. Compare fragility across clinical outcome domains

---

## 2. METHODS (~1,200 words)

### 2.1 Data Source and Extraction
- Cochrane Library systematic reviews (n = 501)
- Meta-analyses extracted: 5,088 total
- Study-level data: effect sizes, variances, sample sizes
- Restriction to logRR for primary analysis (n = 3,556)
- Rationale for effect type restriction

### 2.2 Fragility Metric (R Stability Score)
- Definition: R = 1 - max(individual study influence)
- Influence = change in pooled estimate when study removed
- Interpretation: R < 0.5 = highly fragile
- Advantages over binary Fragility Index

### 2.3 Outcome Domain Classification
- 31 outcome categories based on keywords
- Hierarchical classification algorithm
- Validation: reduction of "Other" category to 23.2%

### 2.4 Predictive Modeling
- Random forest classifier
- Features: k, effect size, SE, tau, I-squared, domain
- 10-fold stratified cross-validation
- Class weighting for imbalance (30% fragile)
- Threshold optimization: Youden's J index

### 2.5 Publication Bias Assessment
- Subsample: meta-analyses with k >= 10 (n = 86)
- Egger's regression test
- Trim-and-fill method
- Sensitivity analysis excluding outliers

### 2.6 Statistical Analysis
- ROC and PR curves (pROC package)
- Calibration assessment by decile
- Domain comparisons with 95% CIs
- FDR correction (Benjamini-Hochberg)
- Effect sizes: Cohen's d
- Software: R 4.5.2, metafor, randomForest, caret

---

## 3. RESULTS (~1,500 words)

### 3.1 Dataset Characteristics
- 501 systematic reviews → 5,088 meta-analyses
- Primary analysis: 3,556 logRR meta-analyses
- Median k = 5 studies per MA (IQR: 3-9)
- 31 outcome domains; largest: Other Clinical (23.2%), Adverse Events (15.2%)

**Table 1: Dataset characteristics**

### 3.2 Fragility Prevalence
- Overall: 25.8% highly fragile (R < 0.5)
- Distribution: right-skewed, median R = 0.75
- Highly stable (R > 0.9): 18.3%

**Figure 2: Stability distribution** (or incorporate into Figure 1)

### 3.3 Predictive Model Performance
- AUC-ROC: 0.788 (95% CI: 0.77-0.81)
- Optimal threshold: 0.35 (Youden's J = 0.423)
- At optimal: Sensitivity 65.7%, Specificity 76.6%
- Balanced accuracy: 71.1%
- Alternative thresholds for different applications

**Figure 2: ROC curve with optimal threshold**
**Table 2: Threshold selection guide**

### 3.4 Predictors of Fragility
- Top 3: Effect magnitude, SE, tau
- Variable importance (Gini impurity)
- Interpretation: larger effects with higher uncertainty = more fragile

**Figure 3 or Supplementary: Feature importance**

### 3.5 Domain-Specific Patterns
- Significantly MORE stable (FDR < 0.05):
  - Quality of Life: R = 0.81 (d = 0.39)
  - Clinical Scores: R = 0.82 (d = 0.45)
  - Mental Health: R = 0.79 (d = 0.29)
  - Musculoskeletal: R = 0.81 (d = 0.45)
- Significantly LESS stable:
  - Mortality: R = 0.67 (d = -0.22)

**Figure 3: Domain forest plot**
**Table 3: Domain-specific results**

### 3.6 Publication Bias
- 86 meta-analyses tested (k >= 10)
- Egger positive: 34.9% (p < 0.1)
- Trim-and-fill: 66.2% had imputed studies
- Median imputed: 5.5 studies
- Effect direction change: 38.2%
- Outlier handling documented

### 3.7 Model Calibration
- Mean absolute calibration error: 0.033
- Good calibration at extremes
- Slight overconfidence at 0.6-0.8 range

**Figure 4: Calibration plot**

---

## 4. DISCUSSION (~1,200 words)

### 4.1 Principal Findings
- One-quarter of Cochrane MAs are highly fragile
- Fragility predictable with moderate accuracy (AUC 0.79)
- Mortality less stable; patient-reported outcomes more stable
- Publication bias common but effect changes modest

### 4.2 Comparison with Previous Literature
- Higher fragility than anticipated based on prior small studies
- Consistent with concerns raised by Ioannidis, Fanelli
- Domain differences align with outcome measurement reliability

### 4.3 Clinical and Methodological Implications
- Guideline developers: check fragility before strong recommendations
- Systematic reviewers: report sensitivity analyses
- Journal editors: require influence diagnostics
- Threshold selection depends on application:
  - Screening: 0.25 (high sensitivity)
  - Confirmation: 0.50 (high specificity)
  - Default: 0.35 (balanced)

### 4.4 Strengths
- Largest meta-epidemiological study of fragility
- Rigorous cross-validation
- Threshold optimization for practical use
- Comprehensive domain analysis
- Publication bias with corrected estimates

### 4.5 Limitations
1. Restricted to Cochrane logRR meta-analyses
2. Bias subsample only 8.2% of eligible MAs
3. Model probabilities require interpretation as relative scores
4. Calibration imperfect at mid-range
5. Domain classification based on keywords, not clinical review
6. Cannot distinguish true fragility from sampling variability

### 4.6 Future Directions
- Validation in non-Cochrane reviews
- Extension to diagnostic accuracy meta-analyses
- Web-based calculator implementation
- Prospective application in living reviews

---

## 5. CONCLUSIONS (~150 words)

- 25.8% of Cochrane meta-analyses are highly fragile
- Predictive model enables prospective identification (AUC 0.79)
- Threshold 0.35 recommended for balanced performance
- Mortality outcomes warrant particular caution
- Routine fragility assessment should become standard practice
- Tool available for systematic reviewers and guideline developers

---

## ACKNOWLEDGEMENTS
[If applicable]

---

## REFERENCES (~40-50 references)

Key citations to include:
1. Cochrane Handbook
2. Walsh et al. (2014) - Fragility Index
3. Viechtbauer & Cheung (2010) - Influence diagnostics
4. Egger et al. (1997) - Publication bias test
5. Duval & Tweedie (2000) - Trim and fill
6. Ioannidis (2005) - Why most findings are false
7. IntHout et al. (2016) - Small-study effects
8. Sterne et al. (2011) - ROBINS-I
9. Higgins et al. (2003) - I-squared
10. DerSimonian & Laird (1986) - Random effects

---

## TABLES

### Table 1: Dataset Characteristics
| Characteristic | Value |
|----------------|-------|
| Systematic reviews | 501 |
| Meta-analyses (total) | 5,088 |
| Meta-analyses (logRR) | 3,556 |
| Outcome domains | 31 |
| Median k per MA | 5 (IQR: 3-9) |
| Fragile (R < 0.5) | 25.8% |

### Table 2: Predictive Model Performance
[From SUBMISSION_ELEMENTS.md - Threshold guide]

### Table 3: Domain-Specific Fragility
[From SUBMISSION_ELEMENTS.md - Supplementary Table S2]

---

## FIGURES

### Main Figures (4)
1. Study flow / PRISMA diagram
2. ROC curve with optimal threshold
3. Domain forest plot
4. Threshold optimization OR Calibration plot

### Supplementary Figures (select from 12 available)
- PR curve
- Calibration plot
- Feature importance
- Stability distribution
- Others as needed

---

## SUPPLEMENTARY MATERIALS

1. Supplementary Table S1: Full threshold optimization results
2. Supplementary Table S2: All domain results (31 domains)
3. Supplementary Table S3: Publication bias detailed results
4. Supplementary Figures S1-S8
5. R code availability statement
6. Sensitivity analyses (if any)

---

## WRITING CHECKLIST

- [ ] Introduction: States gap and objectives clearly
- [ ] Methods: Reproducible detail
- [ ] Results: Matches methods order
- [ ] Discussion: Addresses limitations honestly
- [ ] Abstract: Accurate summary of main findings
- [ ] Tables: Self-explanatory with footnotes
- [ ] Figures: Publication quality, clear legends
- [ ] References: Complete and formatted
- [ ] Supplementary: Organized and referenced in text
- [ ] Word count: Within limit
- [ ] Author contributions: CRediT format
- [ ] Data/code availability: Statements included
- [ ] COI/Funding: Declared

---

*Outline generated: January 4, 2026*
*For submission to: Research Synthesis Methods*
