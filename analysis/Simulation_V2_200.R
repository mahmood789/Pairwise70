################################################################################
# 200-ITERATION SIMULATION: V2 METHODS (Extended)
################################################################################

library(metafor)
library(data.table)

source("C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/Advanced_Pooling_Methods_V2.R")

set.seed(42)

cat(strrep("=", 70), "\n")
cat("200-ITERATION V2 SIMULATION\n")
cat(strrep("=", 70), "\n\n")

N_SIM <- 200
TRUE_EFFECT <- 0.3
results <- data.table()

run_methods <- function(yi, vi) {
  out <- list()

  out$REML <- tryCatch({
    fit <- rma(yi, vi, method = "REML")
    c(coef(fit), fit$se, fit$ci.lb, fit$ci.ub)
  }, error = function(e) c(NA, NA, NA, NA))

  out$HKSJ <- tryCatch({
    fit <- rma(yi, vi, method = "REML", test = "knha")
    c(coef(fit), fit$se, fit$ci.lb, fit$ci.ub)
  }, error = function(e) c(NA, NA, NA, NA))

  out$MWM_v2 <- tryCatch({
    r <- mafi_weighted_ma_v2(yi, vi)
    c(r$estimate, r$se, r$ci_lb, r$ci_ub)
  }, error = function(e) c(NA, NA, NA, NA))

  out$SIT_v2 <- tryCatch({
    r <- sequential_influence_trimming_v2(yi, vi)
    c(r$estimate, r$se, r$ci_lb, r$ci_ub)
  }, error = function(e) c(NA, NA, NA, NA))

  out$ARP_v2 <- tryCatch({
    r <- adaptive_robust_pooling_v2(yi, vi)
    c(r$estimate, r$se, r$ci_lb, r$ci_ub)
  }, error = function(e) c(NA, NA, NA, NA))

  out
}

scenarios <- list(
  # Standard
  list(name = "S1_Standard_k10", k = 10, tau2 = 0.05, outlier = FALSE, bias = 0),
  list(name = "S2_Standard_k20", k = 20, tau2 = 0.05, outlier = FALSE, bias = 0),
  list(name = "S3_Small_k5", k = 5, tau2 = 0.05, outlier = FALSE, bias = 0),
  # Heterogeneity
  list(name = "S4_HighHet", k = 10, tau2 = 0.20, outlier = FALSE, bias = 0),
  list(name = "S5_NoHet", k = 10, tau2 = 0, outlier = FALSE, bias = 0),
  # Outliers
  list(name = "S6_Outlier_k10", k = 10, tau2 = 0.05, outlier = TRUE, bias = 0),
  list(name = "S7_Outlier_k15", k = 15, tau2 = 0.05, outlier = TRUE, bias = 0),
  # Publication bias
  list(name = "S8_PubBias_Mild", k = 15, tau2 = 0.05, outlier = FALSE, bias = 0.3),
  list(name = "S9_PubBias_Mod", k = 15, tau2 = 0.05, outlier = FALSE, bias = 0.5)
)

for (s in seq_along(scenarios)) {
  scen <- scenarios[[s]]
  cat(sprintf("%s (%d iters) ... ", scen$name, N_SIM))
  flush.console()

  for (i in 1:N_SIM) {
    n <- sample(30:150, scen$k, replace = TRUE)
    vi <- 2 / n
    sei <- sqrt(vi)
    theta <- rnorm(scen$k, TRUE_EFFECT, sqrt(scen$tau2))
    yi <- rnorm(scen$k, theta, sei)

    if (scen$outlier) {
      idx <- sample(scen$k, 1)
      yi[idx] <- TRUE_EFFECT + 3 * sqrt(scen$tau2 + vi[idx])
    }

    if (scen$bias > 0) {
      yi <- yi + scen$bias * (sei - mean(sei)) / sd(sei) * 0.1
    }

    res <- run_methods(yi, vi)

    for (m in names(res)) {
      results <- rbind(results, data.table(
        scenario = scen$name, iter = i, method = m,
        est = res[[m]][1], se = res[[m]][2],
        ci_lb = res[[m]][3], ci_ub = res[[m]][4],
        true = TRUE_EFFECT
      ))
    }
  }
  cat("done\n")
}

################################################################################
# RESULTS
################################################################################

cat("\n")
cat(strrep("=", 70), "\n")
cat("PERFORMANCE SUMMARY\n")
cat(strrep("=", 70), "\n\n")

# Overall
overall <- results[, .(
  Bias = round(mean(est - true, na.rm = TRUE), 4),
  RMSE = round(sqrt(mean((est - true)^2, na.rm = TRUE)), 4),
  Coverage = round(mean(ci_lb <= true & ci_ub >= true, na.rm = TRUE) * 100, 1),
  Width = round(mean(ci_ub - ci_lb, na.rm = TRUE), 3)
), by = method]

cat("OVERALL (1800 sims per method):\n")
print(overall[order(RMSE)])

# By category
cat("\n\nSTANDARD SCENARIOS (S1-S3, S5):\n")
std <- results[grepl("Standard|Small|NoHet", scenario), .(
  Bias = round(mean(est - true, na.rm = TRUE), 4),
  RMSE = round(sqrt(mean((est - true)^2, na.rm = TRUE)), 4),
  Coverage = round(mean(ci_lb <= true & ci_ub >= true, na.rm = TRUE) * 100, 1)
), by = method]
print(std[order(RMSE)])

cat("\n\nOUTLIER SCENARIOS (S6-S7):\n")
out <- results[grepl("Outlier", scenario), .(
  Bias = round(mean(est - true, na.rm = TRUE), 4),
  RMSE = round(sqrt(mean((est - true)^2, na.rm = TRUE)), 4),
  Coverage = round(mean(ci_lb <= true & ci_ub >= true, na.rm = TRUE) * 100, 1)
), by = method]
print(out[order(RMSE)])

cat("\n\nHIGH HETEROGENEITY (S4):\n")
het <- results[grepl("HighHet", scenario), .(
  Bias = round(mean(est - true, na.rm = TRUE), 4),
  RMSE = round(sqrt(mean((est - true)^2, na.rm = TRUE)), 4),
  Coverage = round(mean(ci_lb <= true & ci_ub >= true, na.rm = TRUE) * 100, 1)
), by = method]
print(het[order(RMSE)])

cat("\n\nPUBLICATION BIAS (S8-S9):\n")
pb <- results[grepl("PubBias", scenario), .(
  Bias = round(mean(est - true, na.rm = TRUE), 4),
  RMSE = round(sqrt(mean((est - true)^2, na.rm = TRUE)), 4),
  Coverage = round(mean(ci_lb <= true & ci_ub >= true, na.rm = TRUE) * 100, 1)
), by = method]
print(pb[order(abs(Bias))])

################################################################################
# RANKINGS
################################################################################

cat("\n")
cat(strrep("=", 70), "\n")
cat("METHOD RANKINGS\n")
cat(strrep("=", 70), "\n\n")

cat("Best by scenario category:\n")
cat(sprintf("  Standard:  %s (RMSE=%.4f, Cov=%.1f%%)\n",
            std[which.min(RMSE), method], min(std$RMSE), std[which.min(RMSE), Coverage]))
cat(sprintf("  Outliers:  %s (RMSE=%.4f, Cov=%.1f%%)\n",
            out[which.min(RMSE), method], min(out$RMSE), out[which.min(RMSE), Coverage]))
cat(sprintf("  High Het:  %s (RMSE=%.4f, Cov=%.1f%%)\n",
            het[which.min(RMSE), method], min(het$RMSE), het[which.min(RMSE), Coverage]))
cat(sprintf("  Pub Bias:  %s (Bias=%.4f, Cov=%.1f%%)\n",
            pb[which.min(abs(Bias)), method], pb[which.min(abs(Bias)), Bias],
            pb[which.min(abs(Bias)), Coverage]))

cat("\nOverall winner by RMSE: ", overall[which.min(RMSE), method], "\n")
cat("Best coverage (closest to 95%): ", overall[which.min(abs(Coverage - 95)), method], "\n")

################################################################################
# V2 vs STANDARD COMPARISON
################################################################################

cat("\n")
cat(strrep("=", 70), "\n")
cat("V2 IMPROVEMENTS\n")
cat(strrep("=", 70), "\n\n")

reml_rmse <- overall[method == "REML", RMSE]
reml_cov <- overall[method == "REML", Coverage]

cat("RMSE improvement over REML:\n")
for (m in c("HKSJ", "MWM_v2", "SIT_v2", "ARP_v2")) {
  mr <- overall[method == m, RMSE]
  pct <- (reml_rmse - mr) / reml_rmse * 100
  cat(sprintf("  %s: %+.2f%% (%s)\n", m, pct, ifelse(pct > 0, "better", "worse")))
}

cat("\nCoverage improvement over REML:\n")
for (m in c("HKSJ", "MWM_v2", "SIT_v2", "ARP_v2")) {
  mc <- overall[method == m, Coverage]
  cat(sprintf("  %s: %.1f%% (REML: %.1f%%)\n", m, mc, reml_cov))
}

# Save results
output_dir <- "C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/results"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
fwrite(results, file.path(output_dir, "simulation_v2_200_raw.csv"))
fwrite(overall, file.path(output_dir, "simulation_v2_200_overall.csv"))

cat("\n")
cat(strrep("=", 70), "\n")
cat("SIMULATION COMPLETE - Results saved\n")
cat(strrep("=", 70), "\n")
