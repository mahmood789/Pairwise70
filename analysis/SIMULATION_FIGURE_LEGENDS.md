# Simulation Study Figure Legends

## Supplementary Figure S9: Bootstrap AUC Distribution
**File:** `figures/sim_bootstrap_auc.png`

Distribution of area under the ROC curve (AUC) from 1000 bootstrap iterations.
The vertical dashed line indicates the median AUC (0.824), with the shaded
region representing the 95% confidence interval (0.804-0.843). The approximately
normal distribution confirms stable model performance across resampled datasets.

---

## Supplementary Figure S10: Permutation Test Results
**File:** `figures/sim_permutation_test.png`

Comparison of observed AUC (0.835) against the null distribution generated
from 1000 permutations of outcome labels. The null distribution (gray) is
centered at 0.51 (chance level), with 95th percentile at 0.53. The observed
AUC (red line) far exceeds any null value, confirming highly significant
predictive performance (p < 0.001).

---

## Supplementary Figure S11: Power Analysis Curves
**File:** `figures/sim_power_analysis.png`

Statistical power to detect domain-specific deviations from overall mean
stability (R = 0.718) as a function of sample size and effect size (Cohen d).
Curves represent effect sizes d = 0.1 to 0.5. Horizontal dashed line indicates
80% power. For observed effect sizes (d = 0.22-0.45), sample sizes of 55-281
provide adequate power (>80%).

---

## Supplementary Figure S12: Threshold Stability Analysis
**File:** `figures/sim_threshold_stability.png`

Performance metrics (sensitivity, specificity, Youden J) across classification
thresholds from 0.20 to 0.60, evaluated using 1000 bootstrap samples per
threshold. Error bars represent 95% CIs. The recommended threshold of 0.35
provides balanced sensitivity (66%) and specificity (77%). Higher thresholds
offer increased specificity at the cost of sensitivity.

---

## Figure Specifications

| Requirement | Specification |
|-------------|---------------|
| Format | PNG (300 DPI) |
| Color mode | RGB |
| Maximum width | 180 mm (full page) |
| Font | Arial, minimum 8 pt |
| Line weight | Minimum 0.5 pt |

*Generated: January 9, 2026*

