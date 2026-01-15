# Context for Future Chats: Pairwise70 Advanced Pooling Methods

**Last Updated:** January 2026
**Status:** COMPLETE - All development tasks finished

---

## What Was Done

We developed and validated **novel pooling methods for meta-analysis** that outperform standard REML:

### Methods Created (in order of recommendation)

1. **MWM (MAFI-Weighted MA)** - Best overall RMSE (0.0980), 96.2% coverage
   - Uses leave-one-out stability + studentized residuals
   - Downweights influential studies

2. **SIT (Sequential Influence Trimming)** - Best for outliers, near-zero bias
   - Auto-removes outliers using studentized residuals
   - 11% better RMSE than REML for contaminated data

3. **ARP (Adaptive Robust Pooling)** - Good general robustness
   - Combines REML + DL + PM estimators
   - Uses Rubin's rules for variance

### Key Improvements Over Standard Methods

| Issue | Standard REML | Our Methods |
|-------|---------------|-------------|
| Coverage | 92.6% | 95-96% |
| Outlier handling | None | SIT trims automatically |
| Pub bias | None | SIT achieves near-zero bias |

---

## File Locations

### Use These Files:

| Purpose | File |
|---------|------|
| **Package functions** | `Pairwise70/R/advanced_pooling.R` |
| **Full V2 code** | `analysis/Advanced_Pooling_Methods_V2.R` |
| **Documentation** | `analysis/DOCUMENTATION_ADVANCED_POOLING_METHODS.md` |
| **Quick reference** | `analysis/QUICK_REFERENCE_POOLING_METHODS.md` |
| **Changelog** | `analysis/CHANGELOG_ADVANCED_POOLING.md` |

### Simulation Results:

| File | Contents |
|------|----------|
| `analysis/results/simulation_v2_200_raw.csv` | All 9,000 results |
| `analysis/results/simulation_v2_200_overall.csv` | Summary table |

---

## Quick Usage

```r
library(Pairwise70)
library(metafor)

# Compare all methods
compare_pooling_methods(yi, vi)

# Individual methods
mwm <- mafi_weighted_ma(yi, vi)      # Best RMSE
sit <- sequential_influence_trimming(yi, vi)  # Best for outliers
arp <- adaptive_robust_pooling(yi, vi)  # Robust
```

---

## V1 → V2 Evolution

**V1 Issues (identified via 1000-iteration simulation):**
- MWM: RMSE worst, stability weight too high
- SIT: Failed outlier detection (used raw residuals)
- UBSF: Over-corrected bias
- All: Coverage below 95% (used normal CIs)

**V2 Fixes:**
1. T-distribution CIs (df = k-2)
2. Studentized residuals (not raw)
3. Stability weight 0.5 → 0.3
4. UBSF: Blend PET-PEESE (40%) + trim-and-fill (60%)

---

## Simulation Design

- 200 iterations × 9 scenarios = 1,800 per method
- TRUE_EFFECT = 0.3
- Scenarios: Standard, outlier, high heterogeneity, pub bias

---

## What's Left To Do (Future Work)

1. [ ] 1000+ iteration simulation for publication
2. [ ] Compare formally to RoBMA
3. [ ] Add vignette to package
4. [ ] Test on all 501 Pairwise70 datasets
5. [ ] Parallel processing for speed

---

## Key Technical Details

### Why T-distribution?
- Normal CIs assume infinite df, under-cover for small k
- T with df=k-2 fixes coverage

### Why studentized residuals?
- Normalize by expected variance sqrt(vi + tau²)
- Detects true outliers vs just imprecise studies

### Why blend bias corrections?
- PET-PEESE alone over-corrects
- Trim-and-fill more conservative
- 60/40 blend balances both

---

**To continue this work, read:**
1. `DOCUMENTATION_ADVANCED_POOLING_METHODS.md` (complete details)
2. `Advanced_Pooling_Methods_V2.R` (full implementation)
3. `R/advanced_pooling.R` (package version)
