# Manuscript Submission Elements
## Research Synthesis Methods

**Title:** Predicting Meta-Analysis Fragility: A Meta-Epidemiological Study of 3,556 Cochrane Intervention Reviews

---

## ABSTRACT (250 words max)

**Background:** Meta-analyses vary in their robustness to the removal of individual studies, yet no validated tools exist to predict this fragility prospectively.

**Objectives:** To characterize fragility patterns across Cochrane intervention meta-analyses, develop a predictive model for identifying fragile evidence, and assess publication bias prevalence.

**Methods:** We analyzed 3,556 meta-analyses (log risk ratio) from 501 Cochrane systematic reviews. Fragility was quantified using the R stability metric (1 minus maximum study influence), with R < 0.5 indicating high fragility. A random forest classifier with 10-fold stratified cross-validation was developed to predict fragility. Publication bias was assessed using Egger's regression and trim-and-fill in meta-analyses with k >= 10 studies.

**Results:** Overall, 25.8% of meta-analyses were highly fragile (R < 0.5). The predictive model achieved AUC 0.79 (95% CI: 0.77-0.81). At the optimal threshold (0.35), sensitivity was 65.7% and specificity 76.6% (balanced accuracy 71.1%). Effect magnitude, standard error, and heterogeneity were the strongest predictors. Mortality outcomes were significantly less stable (R = 0.67, d = -0.22, p_adj = 0.002), while quality of life (R = 0.81, d = 0.39) and mental health outcomes (R = 0.79, d = 0.29) were more stable. Among 86 meta-analyses tested, 34.9% showed significant funnel plot asymmetry (Egger p < 0.1), with trim-and-fill imputing a median of 5.5 studies.

**Conclusions:** One-quarter of Cochrane meta-analyses are highly fragile. Our predictive model enables prospective identification of vulnerable evidence, supporting more cautious interpretation of fragile findings.

**Word count:** 248

---

## KEYWORDS (5-7)

1. Meta-analysis
2. Fragility index
3. Evidence synthesis
4. Publication bias
5. Cochrane reviews
6. Predictive modeling
7. Meta-epidemiology

---

## DATA AVAILABILITY STATEMENT

The meta-analysis data used in this study were derived from publicly available Cochrane systematic reviews accessed via the Cochrane Library (www.cochranelibrary.com). The processed dataset containing 5,088 meta-analyses with computed stability metrics is available from the corresponding author upon reasonable request. Individual study-level data remain subject to the original Cochrane review licensing terms.

---

## CODE AVAILABILITY STATEMENT

All R code for data processing, statistical analysis, and figure generation is available at [GitHub repository URL] under an MIT license. The analysis was conducted using R version 4.5.2 with the following key packages: metafor (v4.x) for meta-analysis, randomForest (v4.7) for predictive modeling, pROC (v1.18) for ROC analysis, and ggplot2 (v3.5) for visualization. A reproducible analysis script (REVISION_4_EDITORIAL_FINAL.R) and final corrections script (FINAL_CORRECTIONS_BMJ.R) are provided.

---

## CONFLICT OF INTEREST STATEMENT

The authors declare no conflicts of interest. No funding sources had any role in the design, analysis, interpretation, or writing of this manuscript.

---

## FUNDING STATEMENT

[Option 1 - If unfunded:]
This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

[Option 2 - If funded, template:]
This work was supported by [Funding Agency] under Grant [Number]. The funder had no role in study design, data collection and analysis, decision to publish, or preparation of the manuscript.

---

## AUTHOR CONTRIBUTIONS (CRediT Format)

[Author 1]: Conceptualization, Methodology, Software, Formal Analysis, Writing - Original Draft, Visualization

[Author 2]: Validation, Writing - Review & Editing, Supervision

[Additional authors as applicable:]
- Data Curation: [Name]
- Project Administration: [Name]
- Funding Acquisition: [Name]

---

## PRISMA-STYLE FLOW DIAGRAM

```
                    Cochrane Library Database
                              |
                              v
            +----------------------------------+
            | Systematic reviews identified    |
            | n = 501                          |
            +----------------------------------+
                              |
                              v
            +----------------------------------+
            | Meta-analyses extracted          |
            | n = 5,088                        |
            +----------------------------------+
                              |
          +-------------------+-------------------+
          |                   |                   |
          v                   v                   v
    +-----------+       +-----------+       +-----------+
    | Log RR    |       | MD        |       | GIV       |
    | n = 3,556 |       | n = 1,127 |       | n = 405   |
    | (69.9%)   |       | (22.1%)   |       | (8.0%)    |
    +-----------+       +-----------+       +-----------+
          |
          v (Primary Analysis)
    +----------------------------------+
    | Predictive modeling sample       |
    | n = 2,851 (complete cases)       |
    | - Stable: 1,989 (69.8%)          |
    | - Fragile: 862 (30.2%)           |
    +----------------------------------+
          |
          v
    +----------------------------------+
    | Publication bias subsample       |
    | (k >= 10 studies)                |
    | n = 86 meta-analyses             |
    +----------------------------------+
```

---

## SUPPLEMENTARY MATERIALS

### Supplementary Table S1: Threshold Selection Guide

| Use Case | Threshold | Sensitivity | Specificity | PPV | NPV | Balanced Accuracy |
|----------|-----------|-------------|-------------|-----|-----|-------------------|
| Maximum sensitivity | 0.20 | 83.2% | 56.7% | 45.4% | 88.6% | 69.9% |
| Screening/triage | 0.25 | 77.8% | 64.4% | 48.6% | 87.0% | 71.1% |
| Optimal (Youden) | 0.35 | 65.7% | 76.6% | 54.9% | 83.7% | 71.1% |
| Balanced-specific | 0.40 | 59.4% | 81.3% | 57.9% | 82.2% | 70.3% |
| Confirmation | 0.50 | 47.2% | 88.4% | 63.9% | 79.4% | 67.8% |

**Recommendations:**
- Systematic review triage: Use 0.25 (catches 78% of fragile meta-analyses)
- Guideline development: Use 0.35 (balanced, Youden-optimal)
- Individual study assessment: Use 0.50 (higher confidence required)
- Sensitivity analysis: Report results at multiple thresholds

### Supplementary Table S2: Domain-Specific Fragility

| Outcome Domain | n | Mean R | 95% CI | Cohen's d | p (FDR-adjusted) | Direction |
|----------------|---|--------|--------|-----------|------------------|-----------|
| Quality of Life | 189 | 0.809 | 0.776-0.842 | 0.39 | <0.001 | More stable |
| Clinical Scores | 82 | 0.820 | 0.771-0.870 | 0.45 | 0.002 | More stable |
| Mental Health | 139 | 0.791 | 0.749-0.833 | 0.29 | 0.006 | More stable |
| Musculoskeletal | 55 | 0.807 | 0.755-0.860 | 0.45 | 0.008 | More stable |
| Mortality | 281 | 0.670 | 0.644-0.695 | -0.22 | 0.002 | Less stable |
| Overall | 3,556 | 0.718 | 0.710-0.725 | - | - | Reference |

### Supplementary Table S3: Publication Bias Results (Corrected)

| Metric | Value |
|--------|-------|
| Meta-analyses tested (k >= 10) | 86 |
| Proportion of eligible MAs tested | 8.2% |
| Egger's test positive (p < 0.1) | 34.9% |
| Egger's test positive (p < 0.05) | 31.4% |
| Trim-and-fill: any imputation | 66.2% |
| Median studies imputed | 5.5 |
| Mean studies imputed (when >0) | 15.3 |
| Effect direction change after T&F | 38.2% |

*Note: Six meta-analyses with >100 imputed studies were excluded as outliers from mean calculations.*

### Supplementary Figure S1: ROC Curve
File: `figures/ROC_curve_optimized.png`

### Supplementary Figure S2: Precision-Recall Curve
File: `figures/PR_curve.png`

### Supplementary Figure S3: Threshold Optimization
File: `figures/threshold_optimization.png`

### Supplementary Figure S4: Domain Forest Plot
File: `figures/domain_forest_plot.png`

### Supplementary Figure S5: Calibration Plot
File: `figures/calibration_plot.png`

---

## REPORTING CHECKLIST

Research Synthesis Methods recommends following reporting guidelines. Consider:

- [ ] PRISMA 2020 (for systematic review elements)
- [ ] TRIPOD (for predictive model reporting)
- [ ] STROBE (for observational/epidemiological elements)

---

## SUGGESTED REVIEWERS (Optional)

1. [Expert in meta-epidemiology]
2. [Expert in fragility indices]
3. [Expert in Cochrane methodology]
4. [Statistician with publication bias expertise]

---

## COVER LETTER TEMPLATE

Dear Editor,

We are pleased to submit our manuscript entitled "Predicting Meta-Analysis Fragility: A Meta-Epidemiological Study of 3,556 Cochrane Intervention Reviews" for consideration in Research Synthesis Methods.

This study addresses a critical gap in evidence synthesis methodology: the prospective identification of fragile meta-analyses. Using the largest dataset of its kind (3,556 Cochrane meta-analyses), we demonstrate that one-quarter of meta-analyses are highly vulnerable to the removal of single studies, and we provide a validated predictive model (AUC 0.79) to identify such fragility before conclusions are drawn.

Key contributions include:
1. Quantification of fragility prevalence across clinical domains
2. A practical threshold selection guide for different applications
3. Evidence that mortality outcomes are less stable while quality-of-life outcomes are more stable
4. Documentation of publication bias patterns with corrected estimates

This work has direct implications for systematic reviewers, guideline developers, and clinicians interpreting meta-analytic evidence. We believe it aligns well with the scope of Research Synthesis Methods and will be of interest to your readership.

The manuscript has not been published previously and is not under consideration elsewhere. All authors have approved the submitted version.

We suggest the following potential reviewers: [Names]

Thank you for considering our submission.

Sincerely,
[Corresponding Author]

---

*Document generated: January 4, 2026*
*For: Pairwise70 Meta-Epidemiological Project*
