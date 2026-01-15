#' Real Data Examples for Advanced Pooling Methods V3
#'
#' Demonstrates methods on challenging real-world meta-analyses from Pairwise70.
#' Per editorial request: outliers, publication bias, small k, heterogeneous.
#'
#' @author Pairwise70 Team
#' @date January 2026

suppressPackageStartupMessages({
  library(metafor)
  library(data.table)
})

source("C:/Users/user/OneDrive - NHS/Documents/Pairwise70/R/advanced_pooling_v3.R")

cat("=" %>% rep(70) %>% paste(collapse = ""), "\n")
cat("REAL DATA EXAMPLES: ADVANCED POOLING METHODS V3\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n\n")

# === EXAMPLE 1: BCG VACCINE (k=13, Standard Reference) ===

cat("EXAMPLE 1: BCG Vaccine Meta-Analysis\n")
cat("-" %>% rep(50) %>% paste(collapse = ""), "\n")
cat("Context: Classic meta-analysis of BCG vaccine effectiveness\n")
cat("Studies: k=13, moderate heterogeneity (I2~92%)\n\n")

data(dat.bcg)
dat_bcg <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg, data = dat.bcg)

bcg_comparison <- compare_methods_v3(dat_bcg$yi, dat_bcg$vi)
print(bcg_comparison)

# Check for outliers using SIT
sit_bcg <- sequential_influence_trimming_v3(dat_bcg$yi, dat_bcg$vi, bootstrap = TRUE)
cat(sprintf("\nSIT Diagnosis: %d studies trimmed as potential outliers\n", sit_bcg$n_trimmed))
if (sit_bcg$n_trimmed > 0) {
  cat(sprintf("Trimmed study indices: %s\n", paste(sit_bcg$trimmed_studies, collapse = ", ")))
}

# === EXAMPLE 2: ANTI-DEPRESSANTS (High Heterogeneity) ===

cat("\n\nEXAMPLE 2: High Heterogeneity Dataset\n")
cat("-" %>% rep(50) %>% paste(collapse = ""), "\n")
cat("Context: Simulated high I2 scenario based on antidepressant meta-analyses\n\n")

# Use Normand 1999 data (high heterogeneity)
data(dat.normand1999)
dat_normand <- dat.normand1999[, c("m1i", "sd1i", "n1i", "m2i", "sd2i", "n2i")]
dat_high_het <- escalc(measure = "SMD",
                       m1i = m1i, sd1i = sd1i, n1i = n1i,
                       m2i = m2i, sd2i = sd2i, n2i = n2i,
                       data = dat_normand)

fit_reml <- rma(yi, vi, data = dat_high_het)
cat(sprintf("Dataset: k=%d, I2=%.1f%%, tau2=%.4f\n",
            fit_reml$k, fit_reml$I2, fit_reml$tau2))

high_het_comparison <- compare_methods_v3(dat_high_het$yi, dat_high_het$vi)
print(high_het_comparison)

# === EXAMPLE 3: SMALL k (k=5) ===

cat("\n\nEXAMPLE 3: Small Sample Meta-Analysis (k=5)\n")
cat("-" %>% rep(50) %>% paste(collapse = ""), "\n")
cat("Context: First 5 studies from BCG - demonstrates small k behavior\n\n")

dat_small <- dat_bcg[1:5, ]
cat(sprintf("Dataset: k=%d\n", nrow(dat_small)))

small_comparison <- compare_methods_v3(dat_small$yi, dat_small$vi)
print(small_comparison)

cat("\nNote: For k=5, T-distribution CIs are appropriately wider.\n")
cat("df = k-2 = 3, t(0.975, 3) = 3.18 vs normal z = 1.96\n")

# === EXAMPLE 4: POTENTIAL OUTLIER ===

cat("\n\nEXAMPLE 4: Meta-Analysis with Potential Outlier\n")
cat("-" %>% rep(50) %>% paste(collapse = ""), "\n")
cat("Context: Hart 1999 dataset - one study may be outlying\n\n")

data(dat.hart1999)
dat_hart <- escalc(measure = "OR", ai = x1i, n1i = n1i, ci = x2i, n2i = n2i, data = dat.hart1999)

fit_hart <- rma(yi, vi, data = dat_hart)
cat(sprintf("Dataset: k=%d, I2=%.1f%%\n", fit_hart$k, fit_hart$I2))

# Check studentized residuals
rstudent_hart <- rstudent(fit_hart)
cat("\nStudentized residuals (|z| > 2.5 suggests outlier):\n")
for (i in 1:nrow(dat_hart)) {
  if (abs(rstudent_hart$z[i]) > 2.0) {
    cat(sprintf("  Study %d: z = %.2f %s\n", i, rstudent_hart$z[i],
                ifelse(abs(rstudent_hart$z[i]) > 2.5, "***OUTLIER***", "(borderline)")))
  }
}

hart_comparison <- compare_methods_v3(dat_hart$yi, dat_hart$vi)
print(hart_comparison)

sit_hart <- sequential_influence_trimming_v3(dat_hart$yi, dat_hart$vi, bootstrap = TRUE)
cat(sprintf("\nSIT trimmed %d studies: %s\n", sit_hart$n_trimmed,
            ifelse(sit_hart$n_trimmed > 0, paste(sit_hart$trimmed_studies, collapse = ", "), "none")))
cat(sprintf("SIT estimate: %.4f (SE: %.4f, bootstrap: %s)\n",
            sit_hart$estimate, sit_hart$se, sit_hart$se_method))

# === EXAMPLE 5: PUBLICATION BIAS PATTERN ===

cat("\n\nEXAMPLE 5: Potential Publication Bias\n")
cat("-" %>% rep(50) %>% paste(collapse = ""), "\n")
cat("Context: Raudenbush 1985 - teacher expectancy effects\n\n")

data(dat.raudenbush1985)
dat_rauden <- dat.raudenbush1985

fit_rauden <- rma(yi, vi, data = dat_rauden)
cat(sprintf("Dataset: k=%d, I2=%.1f%%\n", fit_rauden$k, fit_rauden$I2))

# Egger's test for asymmetry
egger <- regtest(fit_rauden)
cat(sprintf("Egger's test for funnel asymmetry: z=%.2f, p=%.4f\n",
            egger$zval, egger$pval))

if (egger$pval < 0.1) {
  cat(">>> Evidence of funnel plot asymmetry (potential publication bias)\n")
}

rauden_comparison <- compare_methods_v3(dat_rauden$yi, dat_rauden$vi)
print(rauden_comparison)

# Trim and fill
tf <- trimfill(fit_rauden)
cat(sprintf("\nTrim-and-fill: Added %d studies, adjusted estimate = %.4f\n",
            tf$k0, as.numeric(coef(tf))))

# === SUMMARY TABLE ===

cat("\n\n" %>% rep(70) %>% paste(collapse = ""), "\n")
cat("SUMMARY: METHOD PERFORMANCE ACROSS EXAMPLES\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n\n")

examples <- list(
  BCG = list(yi = dat_bcg$yi, vi = dat_bcg$vi, context = "Standard"),
  High_Het = list(yi = dat_high_het$yi, vi = dat_high_het$vi, context = "High I2"),
  Small_k = list(yi = dat_small$yi, vi = dat_small$vi, context = "k=5"),
  Outlier = list(yi = dat_hart$yi, vi = dat_hart$vi, context = "Potential outlier"),
  Pub_Bias = list(yi = dat_rauden$yi, vi = dat_rauden$vi, context = "Funnel asymmetry")
)

summary_results <- list()

for (ex_name in names(examples)) {
  ex <- examples[[ex_name]]

  reml <- rma(ex$yi, ex$vi, method = "REML")
  mwm <- mafi_weighted_ma_v3(ex$yi, ex$vi)
  sit <- sequential_influence_trimming_v3(ex$yi, ex$vi, bootstrap = FALSE)

  summary_results[[ex_name]] <- data.table(
    Example = ex_name,
    Context = ex$context,
    k = length(ex$yi),
    REML = sprintf("%.3f [%.3f, %.3f]", as.numeric(coef(reml)), reml$ci.lb, reml$ci.ub),
    MWM_v3 = sprintf("%.3f [%.3f, %.3f]", mwm$estimate, mwm$ci_lb, mwm$ci_ub),
    SIT_v3 = sprintf("%.3f [%.3f, %.3f]", sit$estimate, sit$ci_lb, sit$ci_ub),
    SIT_trimmed = sit$n_trimmed
  )
}

summary_table <- rbindlist(summary_results)
print(summary_table)

# === PRACTICAL RECOMMENDATIONS ===

cat("\n\nPRACTICAL RECOMMENDATIONS\n")
cat("-" %>% rep(50) %>% paste(collapse = ""), "\n")

cat("
1. For STANDARD meta-analyses (moderate k, low-moderate I2):
   -> Use MWM_v3 for best overall performance

2. For meta-analyses with SUSPECTED OUTLIERS:
   -> Use SIT_v3 to identify and handle outliers
   -> Report both full and trimmed analyses

3. For VERY SMALL k (k < 5):
   -> Use HKSJ (MWM_v3 and SIT_v3 fall back automatically)
   -> Be cautious interpreting any results

4. For HIGH HETEROGENEITY (I2 > 75%):
   -> ARP_v3 provides robust estimation
   -> Consider MWM_v3 for stability weighting

5. For PUBLICATION BIAS concerns:
   -> SIT_v3 can mitigate some bias effects
   -> Also use trim-and-fill and selection models

")

# === SAVE RESULTS ===

results_dir <- "C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/results"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

fwrite(summary_table, file.path(results_dir, "real_data_examples_summary.csv"))
cat("\nResults saved to: results/real_data_examples_summary.csv\n")
