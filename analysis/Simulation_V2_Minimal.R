################################################################################
# MINIMAL 50-ITERATION SIMULATION: V2 METHODS
################################################################################

library(metafor)
library(data.table)

source("C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis/Advanced_Pooling_Methods_V2.R")

set.seed(42)

cat("MINIMAL V2 SIMULATION (50 iterations)\n")
cat(strrep("=", 50), "\n\n")

N_SIM <- 50
TRUE_EFFECT <- 0.3
results <- data.table()

# Only test core methods (no EMA which is slow)
run_core_methods <- function(yi, vi) {
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

  out
}

# 3 core scenarios
scenarios <- list(
  list(name = "Standard", k = 10, tau2 = 0.05, outlier = FALSE),
  list(name = "Outlier", k = 10, tau2 = 0.05, outlier = TRUE),
  list(name = "High_het", k = 10, tau2 = 0.20, outlier = FALSE)
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

    res <- run_core_methods(yi, vi)

    for (m in names(res)) {
      results <- rbind(results, data.table(
        scenario = scen$name, iter = i, method = m,
        est = res[[m]][1], ci_lb = res[[m]][3], ci_ub = res[[m]][4]
      ))
    }
  }
  cat("done\n")
}

cat("\n")
cat(strrep("=", 50), "\n")
cat("RESULTS\n")
cat(strrep("=", 50), "\n\n")

# Overall
overall <- results[, .(
  Bias = round(mean(est - TRUE_EFFECT, na.rm = TRUE), 4),
  RMSE = round(sqrt(mean((est - TRUE_EFFECT)^2, na.rm = TRUE)), 4),
  Coverage = round(mean(ci_lb <= TRUE_EFFECT & ci_ub >= TRUE_EFFECT, na.rm = TRUE) * 100, 1)
), by = method]

cat("OVERALL:\n")
print(overall[order(RMSE)])

# By scenario
for (sc in unique(results$scenario)) {
  cat(sprintf("\n%s:\n", sc))
  scen_res <- results[scenario == sc, .(
    Bias = round(mean(est - TRUE_EFFECT, na.rm = TRUE), 4),
    RMSE = round(sqrt(mean((est - TRUE_EFFECT)^2, na.rm = TRUE)), 4),
    Coverage = round(mean(ci_lb <= TRUE_EFFECT & ci_ub >= TRUE_EFFECT, na.rm = TRUE) * 100, 1)
  ), by = method]
  print(scen_res[order(RMSE)])
}

cat("\n")
cat(strrep("=", 50), "\n")
cat("COMPLETE\n")
cat(strrep("=", 50), "\n")
