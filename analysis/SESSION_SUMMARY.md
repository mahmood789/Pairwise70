# Pairwise70 Analysis Session Summary
## Quick Reference for Conversation Restart

**Last Updated:** January 2026 (Post-Selenium Testing)
**Full Documentation:** See `COMPLETE_PROJECT_DOCUMENTATION.md` in parent directory
**Status:** ACCEPTED for publication in Research Synthesis Methods

---

## Project Status: COMPLETE

### What Was Accomplished

1. **MAFI Development** - Novel 5-component fragility index
2. **Validation** - 4,424 Cochrane meta-analyses
3. **Editorial Review** - All 8 major + 3 minor concerns addressed
4. **Web Calculator** - Full HTML/JS implementation with case studies
5. **Selenium Testing** - 100% pass rate (16/16 tests)

---

## Quick Start

### To View Results
```r
library(data.table)
mafi <- fread("C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/output/MAFI_all_variants.csv")
summary(mafi$MAFI_5comp)
table(mafi$MAFI_class)
```

### To Use Web Calculator
Open in browser:
```
C:/Users/user/OneDrive - NHS/Documents/Pairwise70/MAFI-Calculator-Complete.html
```

### To Run Selenium Tests
```bash
python C:/Users/user/mafi_functional_test.py
```

---

## Key Files

| File | Purpose |
|------|---------|
| `COMPLETE_PROJECT_DOCUMENTATION.md` | Full project documentation |
| `MAFI-Calculator-Complete.html` | Web calculator with case studies |
| `analysis/output/MAFI_all_variants.csv` | All MAFI scores |
| `analysis/output/EDITORIAL_FINAL_DECISION.md` | Accept decision |

---

## MAFI Formula (Quick Reference)

```
MAFI = 30% DFI + 25% SFI + 20% CFI + 15% Effect + 10% CI
     + (I²/100 × 0.20) + max(0, (1-k/20) × 0.30)

Classification:
  0.00-0.15: Robust
  0.15-0.30: Low Fragility
  0.30-0.50: Moderate Fragility
  0.50-1.00: High Fragility
```

---

## Key Results

- **n = 4,424** meta-analyses from 473 reviews
- **26.8% Robust**, 36.0% Low, 25.1% Moderate, 9.7% High
- **AUC = 0.687** (review-level CV)
- **ICC = 16.1%** (between-review variance)

---

## Test Scripts Available

| Script | Purpose |
|--------|---------|
| `C:/Users/user/mafi_selenium_test.py` | Basic Selenium test |
| `C:/Users/user/mafi_detailed_test.py` | Detailed test (older) |
| `C:/Users/user/mafi_functional_test.py` | Comprehensive test (100% pass) |

---

## Bug Fix Applied (January 2026)

Fixed `switchMainTab()` in MAFI-Calculator-Complete.html:
- Issue: `event.target` undefined when called programmatically
- Fix: Find tab by matching onclick attribute content

---

*For full details, see COMPLETE_PROJECT_DOCUMENTATION.md*
