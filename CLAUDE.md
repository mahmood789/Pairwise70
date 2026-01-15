# MAFI Project - Claude Code Instructions

## Project Context

This is the **MAFI (Meta-Analysis Fragility Index)** project - a novel methodology for assessing the robustness of meta-analysis conclusions. The manuscript has been **ACCEPTED** for publication in Research Synthesis Methods (January 2026).

## Key Information

- **Data:** 4,424 Cochrane meta-analyses from 473 systematic reviews
- **R Package:** Pairwise70
- **Web App:** MAFI-Calculator-Complete.html (fully tested, 100% Selenium pass rate)
- **Status:** Publication-ready

## Important Files

| File | What It Contains |
|------|------------------|
| `COMPLETE_PROJECT_DOCUMENTATION.md` | Full project documentation - READ THIS FIRST |
| `analysis/SESSION_SUMMARY.md` | Quick reference for session restart |
| `MAFI-Calculator-Complete.html` | Web-based MAFI calculator |
| `analysis/output/EDITORIAL_FINAL_DECISION.md` | Journal accept decision |
| `analysis/output/MAFI_all_variants.csv` | All calculated MAFI scores |

## MAFI Formula

```
MAFI = 30% DFI + 25% SFI + 20% CFI + 15% Effect + 10% CI
     + Heterogeneity Penalty + Sample Size Penalty

Classification: Robust (<0.15), Low (0.15-0.30), Moderate (0.30-0.50), High (>0.50)
```

## Common Tasks

### Run Selenium Test
```bash
python C:/Users/user/mafi_functional_test.py
```

### Load R Results
```r
library(data.table)
mafi <- fread("analysis/output/MAFI_all_variants.csv")
```

### Open Web Calculator
```
file:///C:/Users/user/OneDrive - NHS/Documents/Pairwise70/MAFI-Calculator-Complete.html
```

## What NOT to Do

- Do not recreate analysis from scratch - results already exist
- Do not resubmit to journal - already accepted
- Do not modify MAFI weights without re-running sensitivity analysis

## If Asked to Continue Work

1. Read `COMPLETE_PROJECT_DOCUMENTATION.md` first
2. Check `analysis/SESSION_SUMMARY.md` for quick context
3. Existing results are in `analysis/output/` directory
4. Web calculator is fully functional and tested
