# Quick Reference: Advanced Pooling Methods
## Pairwise70 Package

---

## TL;DR - Which Method to Use

| Situation | Method | Command |
|-----------|--------|---------|
| **Default choice** | MWM | `mafi_weighted_ma(yi, vi)` |
| **Suspected outliers** | SIT | `sequential_influence_trimming(yi, vi)` |
| **Publication bias** | SIT | `sequential_influence_trimming(yi, vi)` |
| **Want robustness** | ARP | `adaptive_robust_pooling(yi, vi)` |
| **Compare all** | - | `compare_pooling_methods(yi, vi)` |

---

## Performance Summary (200-iteration simulation)

| Method | RMSE | Coverage | Best For |
|--------|------|----------|----------|
| **MWM** | 0.0980 | 96.2% | Best overall RMSE |
| REML | 0.0986 | 92.6% | Baseline (under-covers) |
| HKSJ | 0.0986 | 94.8% | Standard adjustment |
| ARP | 0.0986 | 95.8% | General robustness |
| **SIT** | 0.1064 | 95.3% | Outliers (11% better) |

---

## Quick Usage

```r
library(Pairwise70)
library(metafor)

# Your data
dat <- escalc(measure = "RR", ai = tpos, bi = tneg,
              ci = cpos, di = cneg, data = your_data)

# Run all methods
compare_pooling_methods(dat$yi, dat$vi)

# Or individual methods
mwm <- mafi_weighted_ma(dat$yi, dat$vi)
sit <- sequential_influence_trimming(dat$yi, dat$vi)
arp <- adaptive_robust_pooling(dat$yi, dat$vi)
```

---

## Key Files

| File | Purpose |
|------|---------|
| `R/advanced_pooling.R` | Package functions (use this) |
| `analysis/Advanced_Pooling_Methods_V2.R` | Full V2 implementation |
| `analysis/DOCUMENTATION_ADVANCED_POOLING_METHODS.md` | Complete documentation |

---

## V2 Improvements Over V1

1. **T-distribution CIs** (df = k-2) - Fixes coverage
2. **Studentized residuals** - Better outlier detection
3. **Reduced stability weight** (0.5 → 0.3) - Less over-adjustment
4. **Blended bias correction** (60% TF + 40% PET-PEESE) - Prevents over-correction

---

## Method Details

### MWM (MAFI-Weighted Meta-Analysis)
- Downweights influential studies
- Uses leave-one-out + studentized residuals
- Best overall RMSE

### SIT (Sequential Influence Trimming)
- Auto-removes outliers (threshold: 2.5 SD)
- Max 20% studies trimmed
- Best for contaminated data

### ARP (Adaptive Robust Pooling)
- Combines REML + DL + PM
- Rubin's rules for variance
- Good general choice

---

**Created:** January 2026
