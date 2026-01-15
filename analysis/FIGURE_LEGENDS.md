# Figure Legends

## Main Manuscript Figures

### Figure 1: Study Flow Diagram
**File:** Create from PRISMA template in SUBMISSION_ELEMENTS.md

Flow diagram showing the selection of meta-analyses from Cochrane systematic reviews. Starting with 501 systematic reviews containing 5,088 meta-analyses, we restricted the primary analysis to 3,556 meta-analyses reporting log risk ratios (logRR). After excluding cases with missing data, 2,851 meta-analyses were included in predictive modeling. Publication bias testing was conducted on 86 meta-analyses with k >= 10 studies.

---

### Figure 2: ROC Curve with Optimal Threshold
**File:** `figures/ROC_curve_optimized.png`

Receiver Operating Characteristic (ROC) curve for the random forest classifier predicting meta-analysis fragility (R < 0.5). The area under the curve (AUC) was 0.788 (95% CI: 0.77-0.81). The optimal threshold of 0.35, determined by maximizing Youden's J index, is indicated by the red point, yielding 65.7% sensitivity and 76.6% specificity. The diagonal dashed line represents chance performance (AUC = 0.5).

---

### Figure 3: Domain-Specific Fragility Forest Plot
**File:** `figures/domain_forest_plot.png`

Forest plot showing mean stability (R) scores by outcome domain with 95% confidence intervals. Domains are sorted by effect size (Cohen's d) relative to the overall mean (R = 0.718, vertical dashed line). Domains with significantly different stability after FDR correction (p < 0.05) are highlighted. Quality of life, clinical scores, mental health, and musculoskeletal outcomes showed significantly greater stability, while mortality outcomes showed significantly lower stability.

---

### Figure 4: Threshold Optimization
**File:** `figures/threshold_optimization.png`

Sensitivity and specificity as functions of the classification threshold for predicting meta-analysis fragility. The optimal threshold (0.35) maximizing Youden's J index is indicated by the vertical dashed line. At this threshold, sensitivity is 65.7% and specificity is 76.6%, yielding a balanced accuracy of 71.1%. Alternative thresholds for different applications are shown: 0.25 for screening (78% sensitivity) and 0.50 for confirmation (88% specificity).

---

## Supplementary Figures

### Supplementary Figure S1: Precision-Recall Curve
**File:** `figures/PR_curve.png`

Precision-Recall curve for the fragility prediction model. Given the class imbalance (30.2% fragile), the PR curve provides additional insight into model performance. The area under the PR curve (AUPRC) is shown. The baseline (random classifier) is indicated by the horizontal dashed line at the prevalence level (0.302).

---

### Supplementary Figure S2: Model Calibration Plot
**File:** `figures/calibration_plot.png`

Calibration plot comparing predicted probabilities of fragility against observed fragility rates across deciles of predicted risk. Perfect calibration is indicated by the diagonal line. The model shows good calibration at the extremes (lowest and highest risk deciles) with slight overconfidence at mid-range probabilities (0.6-0.8), where predicted probabilities exceed observed rates by approximately 7%.

---

### Supplementary Figure S3: Fragility Atlas
**File:** `figures/01_fragility_atlas.png`

Heatmap showing the distribution of meta-analysis stability (R) across outcome domains. Each cell represents the proportion of meta-analyses in each stability category (Highly Fragile: R < 0.5, Moderate: 0.5-0.7, Stable: 0.7-0.9, Highly Stable: R > 0.9) within each domain. Darker colors indicate higher proportions.

---

### Supplementary Figure S4: Stability Distribution
**File:** `figures/02_stability_distribution.png`

Histogram showing the overall distribution of the R stability metric across all 3,556 logRR meta-analyses. The vertical dashed line indicates the fragility threshold (R = 0.5). Approximately 25.8% of meta-analyses fall below this threshold, classified as highly fragile.

---

### Supplementary Figure S5: Heterogeneity Taxonomy
**File:** `figures/03_heterogeneity_taxonomy.png`

Classification of meta-analyses by heterogeneity pattern, showing the relationship between statistical heterogeneity (I-squared, tau) and fragility. Categories include: Low heterogeneity/Stable, Low heterogeneity/Fragile, High heterogeneity/Stable, and High heterogeneity/Fragile.

---

### Supplementary Figure S6: Feature Importance
**File:** `figures/04_feature_importance.png`

Variable importance plot from the random forest classifier, ranked by mean decrease in Gini impurity. The top predictors of fragility are: (1) effect size magnitude (Gini = 279), (2) standard error (Gini = 232), and (3) heterogeneity tau (Gini = 229). These three variables account for the majority of predictive power.

---

### Supplementary Figure S7: Effect Size Calibration
**File:** `figures/05_effect_calibration.png`

Relationship between effect size magnitude and stability across outcome domains, illustrating domain-specific patterns in the effect-fragility relationship.

---

### Supplementary Figure S8: Fragility vs Effect Size
**File:** `figures/06_fragility_vs_effect.png`

Scatter plot showing the relationship between absolute effect size (|theta|) and stability (R) for all meta-analyses. A LOESS smoothing curve illustrates the nonlinear relationship, with larger effect sizes generally associated with lower stability.

---

### Supplementary Figure S7: Summary Dashboard
**File:** `figures/07_dashboard.png`

Multi-panel summary dashboard combining key findings: (A) stability distribution, (B) domain comparison, (C) predictor importance, and (D) ROC curve. This figure provides an overview of the main results for quick reference.

---

## Figure Specifications for Journal Submission

| Requirement | Specification |
|-------------|---------------|
| Format | PNG (provided) or TIFF/EPS (convert if required) |
| Resolution | 300 DPI minimum |
| Color mode | RGB |
| Maximum width | 180 mm (full page) or 85 mm (single column) |
| Font | Arial or Helvetica, minimum 8 pt |
| Line weight | Minimum 0.5 pt |

## File Checklist

- [x] ROC_curve_optimized.png
- [x] PR_curve.png
- [x] threshold_optimization.png
- [x] domain_forest_plot.png
- [x] calibration_plot.png
- [x] 01_fragility_atlas.png
- [x] 02_stability_distribution.png
- [x] 03_heterogeneity_taxonomy.png
- [x] 04_feature_importance.png
- [x] 05_effect_calibration.png
- [x] 06_fragility_vs_effect.png
- [x] 07_dashboard.png

**Total: 12 figures available**

---

*Generated: January 4, 2026*
