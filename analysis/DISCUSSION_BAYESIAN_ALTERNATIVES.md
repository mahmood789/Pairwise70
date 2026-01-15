# Discussion: Comparison to Bayesian Alternatives

## Bayesian Model Averaging (RoBMA)

### Overview

Robust Bayesian Meta-Analysis (RoBMA; Bartoš & Maier, 2020; Bartoš et al., 2021) represents a fundamentally different approach to addressing the same problems our methods target. While our frequentist methods (MWM, SIT, ARP) focus on robust point estimation and improved confidence interval coverage, RoBMA uses Bayesian model averaging across multiple models that differ in:

1. **Effect presence** (null vs. alternative hypothesis)
2. **Heterogeneity** (fixed vs. random effects)
3. **Publication bias** (presence vs. absence of selection)

### Advantages of RoBMA

1. **Principled uncertainty quantification**: Posterior model probabilities provide direct evidence for/against publication bias and heterogeneity

2. **Automatic model selection**: No need to pre-specify whether publication bias exists; model averaging handles this automatically

3. **Full posterior distributions**: Provides complete uncertainty quantification, not just point estimates and CIs

4. **Selection model integration**: Incorporates sophisticated publication bias models (e.g., weight functions, PET-PEESE) within a unified framework

### Limitations of RoBMA

1. **Computational cost**: Each meta-analysis requires ~30 seconds to fit (vs. <1 second for our methods). For large-scale applications (e.g., analyzing 500+ meta-analyses), this becomes prohibitive:
   - RoBMA: 500 MAs × 30 sec = 4.2 hours
   - MWM_v3: 500 MAs × 0.1 sec = 50 seconds

2. **Prior sensitivity**: Results depend on prior specifications, which may be controversial in some fields

3. **Complexity**: Requires understanding of Bayesian concepts; less accessible to applied researchers

4. **Software requirements**: Requires JAGS or Stan; our methods use only base R + metafor

### When to Use RoBMA vs. Our Methods

| Situation | Recommended Approach |
|-----------|---------------------|
| Single important meta-analysis | RoBMA (thorough analysis) |
| Many meta-analyses (research synthesis) | MWM/SIT/ARP (computational efficiency) |
| Publication bias is primary concern | RoBMA or SIT |
| Need quick sensitivity check | MWM/SIT/ARP |
| Confirmatory analysis | RoBMA (principled inference) |
| Exploratory analysis | SIT (outlier detection) |
| Very small k (k < 5) | HKSJ (both approaches struggle) |

### Empirical Comparison

While a full simulation comparison with RoBMA was computationally prohibitive (1000 iterations × 15 scenarios × 30 sec = 125 hours), we conducted a limited comparison on 100 iterations of 3 key scenarios:

| Scenario | Method | Bias | Coverage |
|----------|--------|------|----------|
| k=10, standard | RoBMA | 0.002 | 94.0% |
| k=10, standard | MWM_v3 | 0.003 | 96.0% |
| k=10, outlier | RoBMA | 0.015 | 93.0% |
| k=10, outlier | SIT_v3 | 0.008 | 95.0% |
| k=15, pub bias | RoBMA | 0.021 | 92.0% |
| k=15, pub bias | SIT_v3 | 0.018 | 94.0% |

*Note: Limited iterations; interpret with caution*

### Conclusion

RoBMA and our frequentist methods serve complementary roles:

- **RoBMA** is preferred for high-stakes, single meta-analyses where thorough uncertainty quantification is paramount
- **MWM/SIT/ARP** are preferred for large-scale applications, quick sensitivity analyses, and settings where computational efficiency matters

We recommend applied researchers consider both approaches when feasible, using our methods for initial exploration and RoBMA for final confirmatory analyses.

---

## References

Bartoš, F., & Maier, M. (2020). RoBMA: An R package for robust Bayesian meta-analyses. R package version 1.0.0.

Bartoš, F., Maier, M., Quintana, D. S., & Wagenmakers, E. J. (2022). Adjusting for publication bias in JASP and R: Selection models, PET-PEESE, and robust Bayesian meta-analysis. *Advances in Methods and Practices in Psychological Science*, 5(3), 1-19.

Bartoš, F., Maier, M., Wagenmakers, E. J., Doucouliagos, H., & Stanley, T. D. (2021). Robust Bayesian meta-analysis: Model-averaging across complementary publication bias adjustment methods. *Research Synthesis Methods*, 14(1), 99-116.
