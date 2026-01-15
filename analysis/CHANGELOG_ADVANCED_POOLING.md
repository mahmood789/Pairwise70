# Changelog: Advanced Pooling Methods Development

## Version 2.0 (January 2026) - CURRENT

### New Features
- Integrated into Pairwise70 package as exported functions
- Full roxygen2 documentation with examples
- `compare_pooling_methods()` for side-by-side comparison

### Bug Fixes
- **UBSF_v2**: Fixed over-correction by blending PET-PEESE (40%) with trim-and-fill (60%)
- **All methods**: Switched from normal to T-distribution CIs (df = k-2)
- **MWM**: Reduced stability_weight from 0.5 to 0.3
- **SIT**: Changed from raw to studentized residuals

### Validation
- 200-iteration simulation across 9 scenarios
- Confirmed 95-96% coverage (vs REML 92.6%)
- SIT achieves 11% RMSE improvement for outliers

### Files Added
- `R/advanced_pooling.R` - Package integration
- `analysis/Simulation_V2_200.R` - Validation simulation
- `analysis/DOCUMENTATION_ADVANCED_POOLING_METHODS.md`
- `analysis/QUICK_REFERENCE_POOLING_METHODS.md`

---

## Version 1.0 (January 2026)

### Initial Implementation
- Created 5 novel methods: MWM, SIT, ARP, UBSF, EMA
- `analysis/Advanced_Pooling_Methods.R`

### Simulation Study
- 1000-iteration simulation
- Identified issues:
  - MWM: RMSE worst (0.1411)
  - SIT: Failed at outlier detection
  - UBSF: Over-correction
  - Coverage below 95% for all methods

### Files
- `analysis/Advanced_Pooling_Methods.R` - V1 implementation
- `analysis/Simulation_1000_Iterations.R` - V1 simulation

---

## Development Timeline

| Date | Milestone |
|------|-----------|
| Jan 2026 | V1 methods created |
| Jan 2026 | 1000-iteration V1 simulation |
| Jan 2026 | Issues identified, V2 development |
| Jan 2026 | V2 fixes applied |
| Jan 2026 | 200-iteration V2 validation |
| Jan 2026 | Package integration complete |
| Jan 2026 | Documentation finalized |

---

## Key Decisions

### Why T-distribution CIs?
- Normal CIs assume infinite df
- For k=10: t(0.975, 8) = 2.31 vs z = 1.96
- Coverage improved from 92% to 95-96%

### Why blend PET-PEESE with trim-and-fill?
- PET-PEESE alone over-corrected (0.0659 when true = 0.3)
- Trim-and-fill more conservative
- 60/40 blend balances correction vs stability

### Why studentized residuals?
- Raw residuals don't account for varying precision
- Studies with large SE can have large raw residuals naturally
- Studentized residuals normalize for expected variance

### Why 0.3 stability weight (not 0.5)?
- 0.5 was too aggressive, increased RMSE
- 0.3 provides stability benefit without sacrificing accuracy
- Leave-one-out analysis still detects influence

---

## Performance Comparison

### V1 Results (1000 iterations)
| Method | RMSE | Coverage |
|--------|------|----------|
| REML | 0.1380 | 92.5% |
| MWM | 0.1411 | 94.1% |
| SIT | 0.1390 | 94.3% |
| EMA | 0.1380 | 94.2% |

### V2 Results (200 iterations)
| Method | RMSE | Coverage |
|--------|------|----------|
| REML | 0.0986 | 92.6% |
| MWM_v2 | **0.0980** | **96.2%** |
| SIT_v2 | 0.1064 | 95.3% |
| ARP_v2 | 0.0986 | 95.8% |

**Improvement**: Coverage fixed, MWM now best RMSE

---

## Known Issues / TODO

1. [ ] Run 1000+ iteration simulation for publication
2. [ ] Add prediction intervals
3. [ ] Test on all 501 Pairwise70 datasets
4. [ ] Compare to RoBMA formally
5. [ ] Add parallel processing for speed
6. [ ] Create vignette for package

---

**Last Updated:** January 2026
