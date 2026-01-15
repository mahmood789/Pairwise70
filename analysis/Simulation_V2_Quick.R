################################################################################
# QUICK 100-ITERATION SIMULATION: V2 METHODS
# Fast validation for immediate feedback
################################################################################

library(metafor)
library(data.table)

source("C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/Advanced_Pooling_Methods_V2.R")

set.seed(42)

cat(strrep("=", 70), "\n")
cat("QUICK 100-ITERATION V2 SIMULATION\n")
cat(strrep("=", 70), "\n\n")

N_SIM <- 100
TRUE_EFFECT <- 0.3

results <- data.table()

################################################################################
# SIMPLIFIED METHOD RUNNER
################################################################################

run_methods <- function(yi, vi) {
  k <- length(yi)
  out <- list()

  # REML
  out$REML <- tryCatch({
    fit <- rma(yi, vi, method = "REML")
    c(coef(fit), fit$se, fit$ci.lb, fit$ci.ub)
  }, error = function(e) c(NA, NA, NA, NA))

  # HKSJ
  out$HKSJ <- tryCatch({
    fit <- rma(yi, vi, method = "REML", test = "knha")
    c(coef(fit), fit$se, fit$ci.lb, fit$ci.ub)
  }, error = function(e) c(NA, NA, NA, NA))

  # MWM_v2
  out$MWM_v2 <- tryCatch({
    r <- mafi_weighted_ma_v2(yi, vi)
    c(r$estimate, r$se, r$ci_lb, r$ci_ub)
  }, error = function(e) c(NA, NA, NA, NA))

  # SIT_v2
  out$SIT_v2 <- tryCatch({
    r <- sequential_influence_trimming_v2(yi, vi)
    c(r$estimate, r$se, r$ci_lb, r$ci_ub)
  }, error = function(e) c(NA, NA, NA, NA))

  # EMA_v2
  out$EMA_v2 <- tryCatch({
    r <- ensemble_meta_analysis_v2(yi, vi)
    c(r$estimate, r$se, r$ci_lb, r$ci_ub)
  }, error = function(e) c(NA, NA, NA, NA))

  out
}

################################################################################
# SCENARIOS
################################################################################

scenarios <- list(
  list(name = "Standard k=10", k = 10, tau2 = 0.05, outlier = FALSE, bias = 0),
  list(name = "High het k=10", k = 10, tau2 = 0.15, outlier = FALSE, bias = 0),
  list(name = "Small k=5", k = 5, tau2 = 0.05, outlier = FALSE, bias = 0),
  list(name = "Outlier k=10", k = 10, tau2 = 0.05, outlier = TRUE, bias = 0),
  list(name = "Pub bias k=15", k = 15, tau2 = 0.05, outlier = FALSE, bias = 0.5)
)

for (s in seq_along(scenarios)) {
  scen <- scenarios[[s]]
  cat(sprintf("Scenario %d: %s ... ", s, scen$name))
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
        scenario = s, name = scen$name, iter = i, method = m,
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
cat("RESULTS\n")
cat(strrep("=", 70), "\n\n")

# Overall
overall <- results[, .(
  Bias = round(mean(est - true, na.rm = TRUE), 4),
  RMSE = round(sqrt(mean((est - true)^2, na.rm = TRUE)), 4),
  Coverage = round(mean(ci_lb <= true & ci_ub >= true, na.rm = TRUE) * 100, 1),
  Width = round(mean(ci_ub - ci_lb, na.rm = TRUE), 3)
), by = method]

cat("OVERALL (500 sims per method):\n")
print(overall[order(RMSE)])

# By scenario
cat("\n\nBY SCENARIO:\n")
for (s in 1:5) {
  cat(sprintf("\n%s:\n", scenarios[[s]]$name))
  scen_res <- results[scenario == s, .(
    Bias = round(mean(est - true, na.rm = TRUE), 4),
    RMSE = round(sqrt(mean((est - true)^2, na.rm = TRUE)), 4),
    Coverage = round(mean(ci_lb <= true & ci_ub >= true, na.rm = TRUE) * 100, 1)
  ), by = method]
  print(scen_res[order(RMSE)])
}

################################################################################
# KEY FINDINGS
################################################################################

cat("\n")
cat(strrep("=", 70), "\n")
cat("KEY FINDINGS\n")
cat(strrep("=", 70), "\n\n")

# Best RMSE
best_rmse <- overall[which.min(RMSE)]
cat(sprintf("Best overall RMSE: %s (%.4f)\n", best_rmse$method, best_rmse$RMSE))

# Best coverage
best_cov <- overall[which.min(abs(Coverage - 95))]
cat(sprintf("Best coverage: %s (%.1f%%)\n", best_cov$method, best_cov$Coverage))

# RMSE improvement over REML
reml_rmse <- overall[method == "REML", RMSE]
for (m in c("HKSJ", "MWM_v2", "SIT_v2", "EMA_v2")) {
  mr <- overall[method == m, RMSE]
  pct <- (reml_rmse - mr) / reml_rmse * 100
  cat(sprintf("%s vs REML: %+.1f%% RMSE\n", m, pct))
}

# Outlier scenario winner
out_res <- results[scenario == 4, .(RMSE = sqrt(mean((est - true)^2, na.rm = TRUE))), by = method]
best_out <- out_res[which.min(RMSE)]
cat(sprintf("\nBest for outliers: %s (RMSE=%.4f)\n", best_out$method, best_out$RMSE))

# Pub bias scenario winner
pb_res <- results[scenario == 5, .(Bias = mean(est - true, na.rm = TRUE)), by = method]
best_pb <- pb_res[which.min(abs(Bias))]
cat(sprintf("Best for pub bias: %s (Bias=%.4f)\n", best_pb$method, best_pb$Bias))

cat("\n")
cat(strrep("=", 70), "\n")
cat("SIMULATION COMPLETE\n")
cat(strrep("=", 70), "\n")
