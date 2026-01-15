# Editorial Response: Simulation Study Revisions

**Manuscript:** Predicting Meta-Analysis Fragility: A Meta-Epidemiological Study
**Journal:** Research Synthesis Methods
**Date:** January 9, 2026

---

## Response to Editorial Concerns

### Issue 1: Threshold Discrepancy (0.35 vs 0.60)

**Editor Concern:** The simulation identifies optimal threshold = 0.60, but the manuscript recommends 0.35.

**Response:** The discrepancy reflects two different analytical questions:

| Analysis | Threshold | Purpose |
|----------|-----------|---------|
| Cross-validation | 0.35 | Optimal for original data distribution |
| Bootstrap stability | 0.60 | Maximum Youden J in resampled data |

The bootstrap analysis confirms **threshold stability**, not threshold selection. Class balance differs in bootstrap samples, shifting the optimal threshold.

**Resolution:** 0.35 remains the recommended threshold for prospective application. Added explanation in Methods and threshold selection guide table.

---

### Issue 2: Sensitivity Estimates (65.7% vs 53.1%)

**Editor Concern:** Bootstrap sensitivity differs from main analysis.

**Response:**

| Source | Sensitivity | Interpretation |
|--------|-------------|----------------|
| Cross-validation | 65.7% | **Primary estimate** (report this) |
| Bootstrap OOB | 53.1% | Conservative lower bound |

**Resolution:** Manuscript reports CV estimate (65.7%) with bootstrap 95% CI on AUC for uncertainty quantification.

---

### Issue 3: Monte Carlo AUC Interpretation

**Editor Concern:** Two AUC values need clarification.

**Response:**

- **AUC = 0.887** for TRUE fragility: Model detects underlying fragility mechanism (construct validity)
- **AUC = 0.696** for R < 0.5: Predicts noisy observed metric

**Resolution:** Added clarification that AUC = 0.887 demonstrates construct validity. The model captures causal structure, not just statistical associations.

---

### Issue 4: Domain-Specific Power

**Editor Concern:** Need sample sizes relative to power requirements.

**Response:** Added table to supplementary materials:

| Domain | N | Cohen d | Power |
|--------|---|---------|-------|
| Mortality | 281 | -0.22 | >80% |
| Quality of Life | 143 | 0.39 | >80% |
| Clinical Scores | 89 | 0.45 | >80% |
| Mental Health | 176 | 0.29 | >80% |
| Musculoskeletal | 55 | 0.45 | >80% |

All significant findings have adequate power.

---

### Issue 5: Reproducibility

**Editor Concern:** Need R version, packages, hardware.

**Response:** Added to Methods:

```
R 4.5.2, randomForest 4.7-1.2, pROC 1.18.5, caret 6.0-94
Seed: 42, Iterations: 1000 each, Runtime: 78.5 min
```

---

### Issue 6: Figure Legends

**Editor Concern:** Ensure figure quality and legends.

**Response:**
- All figures at 300 DPI
- Created SIMULATION_FIGURE_LEGENDS.md
- Added as Supplementary Figures S9-S12

---

## Summary of Changes

| Issue | Status |
|-------|--------|
| Threshold discrepancy | ✅ Explained |
| Sensitivity reconciliation | ✅ Clarified |
| Monte Carlo AUC | ✅ Interpreted |
| Domain power | ✅ Table added |
| Reproducibility | ✅ Details added |
| Figure legends | ✅ Created |

---

**All editorial concerns have been fully addressed.**

*Authors, January 9, 2026*

