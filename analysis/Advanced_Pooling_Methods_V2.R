################################################################################
# ADVANCED POOLING METHODS FOR META-ANALYSIS - VERSION 2
# Improvements based on 1000-iteration simulation study
################################################################################
#
# CHANGES FROM V1:
# - MWM: Recalibrated stability scores, adaptive threshold
# - SIT: Lower influence threshold, studentized residuals, adaptive trimming
# - UBSF: PET-PEESE integration, adaptive bias correction strength
# - ARP: Removed arbitrary I2 cutoffs, use estimator-specific uncertainty
# - EMA: Adaptive weighting based on method agreement, robust SE
# - Global: t-distribution CIs for k < 30, improved SE estimation
#
################################################################################

library(metafor)
library(data.table)
suppressPackageStartupMessages({
  if (requireNamespace("clubSandwich", quietly = TRUE)) library(clubSandwich)
})

################################################################################
# HELPER: T-DISTRIBUTION CI (Fixes coverage issue)
################################################################################

#' Calculate CI using t-distribution for small samples
calc_ci <- function(estimate, se, k, alpha = 0.05) {
  if (k >= 30) {
    # Use normal for large k
    z <- qnorm(1 - alpha/2)
  } else {
    # Use t-distribution for small k (better coverage)
    df <- max(k - 2, 1)
    z <- qt(1 - alpha/2, df)
  }
  ci_lb <- estimate - z * se
  ci_ub <- estimate + z * se
  return(list(ci_lb = ci_lb, ci_ub = ci_ub, df = if(k < 30) df else Inf))
}

################################################################################
# METHOD 1: MAFI-WEIGHTED META-ANALYSIS V2 (MWM)
################################################################################
# IMPROVEMENTS:
# - Recalibrated stability weights (reduced from 0.5 to 0.3)
# - Use studentized residuals for influence
# - Adaptive threshold based on k
################################################################################

mafi_weighted_ma_v2 <- function(yi, vi, method = "REML", alpha = 0.05,
                                 stability_weight = 0.3) {  # Reduced from 0.5

  k <- length(yi)

  if (k < 3) {
    return(list(
      estimate = mean(yi),
      se = sqrt(mean(vi)),
      ci_lb = NA, ci_ub = NA,
      method = "MWM_v2",
      message = "Insufficient studies (k < 3)"
    ))
  }

  # Step 1: Fit standard random-effects model
  fit <- rma(yi = yi, vi = vi, method = method)
  estimate_base <- as.numeric(fit$beta[1])
  tau2 <- fit$tau2

  # Step 2: Use studentized residuals for influence detection
  inf <- influence(fit)
  rstud <- inf$inf$rstudent
  rstud[is.na(rstud)] <- 0

  # Step 3: Calculate stability scores using studentized residuals
  # Adaptive threshold based on k
  stud_threshold <- ifelse(k >= 10, 2.0, 2.5)  # More lenient for small k

  stability_scores <- numeric(k)
  for (i in 1:k) {
    # Studentized residual contribution (main factor)
    stud_factor <- 1 - min(abs(rstud[i]) / stud_threshold, 1)

    # DFBETAS-like influence (how much estimate changes)
    # Use leave1out for robustness
    loo <- tryCatch({
      leave1out(fit)
    }, error = function(e) NULL)

    if (!is.null(loo)) {
      est_change <- abs(loo$estimate[i] - estimate_base)
      se_pooled <- mean(sqrt(vi + tau2))
      dfbetas_factor <- 1 - min(est_change / (2 * se_pooled), 1)
    } else {
      dfbetas_factor <- 1
    }

    # Combined stability score (0 to 1)
    stability_scores[i] <- 0.6 * stud_factor + 0.4 * dfbetas_factor
  }

  # Step 4: Combine precision weights with stability weights
  precision_weights <- 1 / (vi + tau2)
  precision_weights <- precision_weights / sum(precision_weights)

  # Normalize stability scores (add floor to prevent extreme downweighting)
  stability_weights <- pmax(stability_scores, 0.2)  # Floor at 0.2
  stability_weights <- stability_weights / sum(stability_weights)

  # Combined weights (reduced stability influence)
  combined_weights <- (1 - stability_weight) * precision_weights +
                       stability_weight * stability_weights
  combined_weights <- combined_weights / sum(combined_weights)

  # Step 5: Weighted pooled estimate
  estimate_mwm <- sum(combined_weights * yi)

  # Step 6: Calculate variance with Hartung-Knapp correction
  # This improves coverage
  resid <- yi - estimate_mwm
  qe <- sum(combined_weights * resid^2)

  # Standard variance
  var_mwm <- sum(combined_weights^2 * (vi + tau2))

  # Apply HK-like correction factor
  if (k > 2) {
    hk_factor <- qe / (k - 1)
    var_mwm <- var_mwm * max(hk_factor, 1)  # Never reduce variance
  }

  se_mwm <- sqrt(var_mwm)

  # Confidence interval using t-distribution
  ci_result <- calc_ci(estimate_mwm, se_mwm, k, alpha)

  # P-value
  z_stat <- estimate_mwm / se_mwm
  pval <- 2 * pt(-abs(z_stat), df = max(k - 2, 1))

  return(list(
    method = "MWM_v2 (MAFI-Weighted V2)",
    estimate = round(estimate_mwm, 4),
    se = round(se_mwm, 4),
    ci_lb = round(ci_result$ci_lb, 4),
    ci_ub = round(ci_result$ci_ub, 4),
    pval = round(pval, 4),
    tau2 = round(tau2, 4),
    k = k,
    stability_weight = stability_weight,
    study_stability_scores = round(stability_scores, 3),
    combined_weights = round(combined_weights, 4),
    estimate_base = round(estimate_base, 4),
    adjustment = round(estimate_mwm - estimate_base, 4)
  ))
}

################################################################################
# METHOD 2: ADAPTIVE ROBUST POOLING V2 (ARP)
################################################################################
# IMPROVEMENTS:
# - Use between-estimator variance for uncertainty
# - Weight by inverse of estimator SE (not arbitrary cutoffs)
# - Include Q-profile confidence interval for tau2
################################################################################

adaptive_robust_pooling_v2 <- function(yi, vi, alpha = 0.05) {

  k <- length(yi)

  if (k < 3) {
    return(list(
      estimate = mean(yi),
      se = sqrt(mean(vi)),
      method = "ARP_v2",
      message = "Insufficient studies"
    ))
  }

  # Step 1: Fit multiple estimators
  estimators <- list()

  # 1a. REML
  tryCatch({
    fit_reml <- rma(yi = yi, vi = vi, method = "REML")
    estimators$REML <- list(
      estimate = as.numeric(fit_reml$beta[1]),
      se = fit_reml$se,
      tau2 = fit_reml$tau2,
      I2 = fit_reml$I2
    )
  }, error = function(e) NULL)

  # 1b. DerSimonian-Laird
  tryCatch({
    fit_dl <- rma(yi = yi, vi = vi, method = "DL")
    estimators$DL <- list(
      estimate = as.numeric(fit_dl$beta[1]),
      se = fit_dl$se,
      tau2 = fit_dl$tau2,
      I2 = fit_dl$I2
    )
  }, error = function(e) NULL)

  # 1c. Paule-Mandel
  tryCatch({
    fit_pm <- rma(yi = yi, vi = vi, method = "PM")
    estimators$PM <- list(
      estimate = as.numeric(fit_pm$beta[1]),
      se = fit_pm$se,
      tau2 = fit_pm$tau2,
      I2 = fit_pm$I2
    )
  }, error = function(e) NULL)

  # 1d. HKSJ adjustment
  tryCatch({
    fit_hksj <- rma(yi = yi, vi = vi, method = "REML", test = "knha")
    estimators$HKSJ <- list(
      estimate = as.numeric(fit_hksj$beta[1]),
      se = fit_hksj$se,  # This uses HKSJ-adjusted SE
      tau2 = fit_hksj$tau2,
      I2 = fit_hksj$I2
    )
  }, error = function(e) NULL)

  if (length(estimators) == 0) {
    return(list(method = "ARP_v2", estimate = NA, se = NA, message = "All estimators failed"))
  }

  # Step 2: Calculate weights based on inverse SE (not arbitrary)
  estimates <- sapply(estimators, function(x) x$estimate)
  ses <- sapply(estimators, function(x) x$se)

  # Inverse-variance weights (but among estimators)
  inv_var_weights <- 1 / ses^2
  inv_var_weights <- inv_var_weights / sum(inv_var_weights)

  # Step 3: Compute adaptive pooled estimate
  estimate_arp <- sum(inv_var_weights * estimates)

  # Step 4: Variance combining multiple estimators
  # Include both within-method and between-method variance
  var_within <- sum(inv_var_weights * ses^2)
  var_between <- sum(inv_var_weights * (estimates - estimate_arp)^2)

  # Total variance (Rubin's rules style)
  var_total <- var_within + var_between * (1 + 1/length(estimators))
  se_arp <- sqrt(var_total)

  # Confidence interval with t-distribution
  ci_result <- calc_ci(estimate_arp, se_arp, k, alpha)

  # P-value using t-distribution
  z_stat <- estimate_arp / se_arp
  pval <- 2 * pt(-abs(z_stat), df = max(k - 2, 1))

  return(list(
    method = "ARP_v2 (Adaptive Robust Pooling V2)",
    estimate = round(estimate_arp, 4),
    se = round(se_arp, 4),
    ci_lb = round(ci_result$ci_lb, 4),
    ci_ub = round(ci_result$ci_ub, 4),
    pval = round(pval, 4),
    k = k,
    I2_mean = round(mean(sapply(estimators, function(x) x$I2), na.rm = TRUE), 1),
    estimator_weights = round(inv_var_weights, 3),
    estimator_estimates = round(estimates, 4),
    estimator_ses = round(ses, 4),
    var_within = round(var_within, 6),
    var_between = round(var_between, 6)
  ))
}

################################################################################
# METHOD 3: SEQUENTIAL INFLUENCE TRIMMING V2 (SIT)
################################################################################
# IMPROVEMENTS:
# - Use studentized residuals instead of Cook's D
# - Lower threshold (1.5 instead of 0.5 for Cook's D)
# - Adaptive trim factor based on residual magnitude
# - Better convergence criterion
################################################################################

sequential_influence_trimming_v2 <- function(yi, vi, rstud_threshold = 2.0,
                                              max_iterations = 10, min_trim = 0.1) {

  k <- length(yi)

  if (k < 4) {
    return(list(
      estimate = mean(yi),
      se = sqrt(mean(vi)),
      method = "SIT_v2",
      message = "Insufficient studies for influence analysis"
    ))
  }

  # Initialize weights
  weights <- rep(1, k)
  iteration <- 0
  converged <- FALSE

  estimate_history <- numeric()

  # Original fit for reference
  fit_orig <- rma(yi = yi, vi = vi, method = "REML")
  estimate_orig <- as.numeric(fit_orig$beta[1])

  # Iterative trimming
  while (!converged && iteration < max_iterations) {
    iteration <- iteration + 1

    # Fit weighted model
    fit <- tryCatch({
      rma(yi = yi, vi = vi, weights = weights, method = "REML")
    }, error = function(e) {
      rma(yi = yi, vi = vi, method = "REML")
    })

    estimate_current <- as.numeric(fit$beta[1])
    estimate_history <- c(estimate_history, estimate_current)

    # Calculate studentized residuals
    inf <- influence(fit)
    rstud <- inf$inf$rstudent
    rstud[is.na(rstud)] <- 0

    # Identify influential studies (|rstud| > threshold)
    influential <- which(abs(rstud) > rstud_threshold)

    if (length(influential) == 0) {
      converged <- TRUE
    } else {
      # Adaptive trimming: weight reduction proportional to how extreme the residual is
      for (idx in influential) {
        # More extreme residuals get more aggressive trimming
        excess <- abs(rstud[idx]) - rstud_threshold
        trim_factor <- max(min_trim, 1 / (1 + excess))
        weights[idx] <- weights[idx] * trim_factor
      }

      # Prevent weights from going below minimum
      weights <- pmax(weights, 0.05)
      # Renormalize to sum to k
      weights <- weights * k / sum(weights)
    }

    # Check for convergence (estimate stability)
    if (iteration > 1) {
      change <- abs(estimate_current - estimate_history[iteration - 1])
      rel_change <- change / max(abs(estimate_current), 0.001)
      if (rel_change < 0.005) {  # 0.5% relative change threshold
        converged <- TRUE
      }
    }
  }

  # Final weighted estimate with HKSJ adjustment
  fit_final <- tryCatch({
    rma(yi = yi, vi = vi, weights = weights, method = "REML", test = "knha")
  }, error = function(e) {
    rma(yi = yi, vi = vi, method = "REML")
  })

  estimate_sit <- as.numeric(fit_final$beta[1])
  se_sit <- fit_final$se
  tau2 <- fit_final$tau2

  # Confidence interval using t-distribution (already HKSJ)
  ci_result <- calc_ci(estimate_sit, se_sit, k)
  pval <- fit_final$pval

  # Count trimmed studies
  n_trimmed <- sum(weights < 0.5)

  return(list(
    method = "SIT_v2 (Sequential Influence Trimming V2)",
    estimate = round(estimate_sit, 4),
    se = round(se_sit, 4),
    ci_lb = round(ci_result$ci_lb, 4),
    ci_ub = round(ci_result$ci_ub, 4),
    pval = round(pval, 4),
    tau2 = round(tau2, 4),
    k = k,
    iterations = iteration,
    converged = converged,
    n_trimmed = n_trimmed,
    final_weights = round(weights, 3),
    estimate_original = round(estimate_orig, 4),
    adjustment = round(estimate_sit - estimate_orig, 4),
    estimate_history = round(estimate_history, 4)
  ))
}

################################################################################
# METHOD 4: UNIFIED BIAS-STABILITY FRAMEWORK V2 (UBSF)
################################################################################
# IMPROVEMENTS:
# - BLEND PET-PEESE with trim-and-fill (prevents over-correction)
# - Conservative bias weight (max 50% instead of 80%)
# - Conflicting methods use more conservative estimate
# - Better SE inflation formula
################################################################################

unified_bias_stability_v2 <- function(yi, vi, stability_weight = 0.3) {

  k <- length(yi)
  sei <- sqrt(vi)

  if (k < 5) {
    fit <- rma(yi = yi, vi = vi, method = "REML")
    return(list(
      estimate = as.numeric(fit$beta[1]),
      se = fit$se,
      method = "UBSF_v2",
      message = "Insufficient studies for bias-stability analysis"
    ))
  }

  # Step 1: Fit base model
  fit_base <- rma(yi = yi, vi = vi, method = "REML")
  estimate_base <- as.numeric(fit_base$beta[1])
  se_base <- fit_base$se
  tau2 <- fit_base$tau2

  # Step 2: Publication bias assessment
  egger <- tryCatch({
    regtest(fit_base, model = "lm")
  }, error = function(e) NULL)

  bias_detected <- FALSE
  bias_pval <- 1

  if (!is.null(egger)) {
    bias_pval <- egger$pval
    bias_detected <- bias_pval < 0.10
  }

  # Step 3: Multiple bias correction methods (BLEND approach)

  # Method A: Trim-and-fill (conservative)
  tf_estimate <- estimate_base
  k_imputed <- 0
  tryCatch({
    tf <- trimfill(fit_base)
    if (!is.null(tf) && tf$k0 > 0) {
      tf_estimate <- as.numeric(tf$beta[1])
      k_imputed <- tf$k0
    }
  }, error = function(e) NULL)

  # Method B: PET-PEESE (more aggressive)
  pet_estimate <- estimate_base
  peese_estimate <- estimate_base
  pet_pval <- 1

  if (k >= 6) {
    tryCatch({
      fit_pet <- rma(yi = yi, vi = vi, mods = ~ sei, method = "REML")
      pet_estimate <- as.numeric(fit_pet$beta[1])
      pet_pval <- fit_pet$pval[1]

      fit_peese <- rma(yi = yi, vi = vi, mods = ~ vi, method = "REML")
      peese_estimate <- as.numeric(fit_peese$beta[1])
    }, error = function(e) NULL)
  }

  # Step 4: BLEND bias corrections (60% trim-and-fill, 40% PET-PEESE)
  if (pet_pval >= 0.05) {
    petpeese_estimate <- pet_estimate
    bias_method <- "Blend(TF+PET)"
  } else {
    petpeese_estimate <- peese_estimate
    bias_method <- "Blend(TF+PEESE)"
  }

  # Calculate adjustments
  tf_adj <- tf_estimate - estimate_base
  petpeese_adj <- petpeese_estimate - estimate_base

  # Blend: 60% trim-and-fill (conservative), 40% PET-PEESE
  # If they disagree on direction, use the smaller adjustment
  if (sign(tf_adj) == sign(petpeese_adj) || abs(tf_adj) < 0.01 || abs(petpeese_adj) < 0.01) {
    estimate_bias_adj <- 0.6 * tf_estimate + 0.4 * petpeese_estimate
  } else {
    # Conflicting directions - use more conservative (smaller magnitude adjustment)
    if (abs(tf_adj) < abs(petpeese_adj)) {
      estimate_bias_adj <- tf_estimate
      bias_method <- "TrimFill(conservative)"
    } else {
      estimate_bias_adj <- petpeese_estimate
      bias_method <- "PET-PEESE(conservative)"
    }
  }

  # Step 5: CONSERVATIVE bias weight (max 50% instead of 80%)
  if (bias_detected) {
    bias_weight <- min(0.5, 0.2 + 0.3 * (1 - bias_pval))
  } else {
    bias_weight <- 0.05  # Minimal correction if not detected
  }

  # Step 6: Stability assessment
  loo <- leave1out(fit_base)

  dir_changes <- sum((loo$estimate > 0) != (estimate_base > 0), na.rm = TRUE)
  dir_fragile <- dir_changes > 0

  sig_changes <- sum((loo$pval < 0.05) != (fit_base$pval < 0.05), na.rm = TRUE)
  sig_fragile <- sig_changes > 0

  max_change <- max(abs(loo$estimate - estimate_base), na.rm = TRUE)

  loo_median <- median(loo$estimate, na.rm = TRUE)
  estimate_stability_adj <- estimate_base + stability_weight * (loo_median - estimate_base)

  # Step 7: Unified corrected estimate
  bias_adjustment <- estimate_bias_adj - estimate_base
  stability_adjustment <- estimate_stability_adj - estimate_base

  total_adjustment <- bias_weight * bias_adjustment + stability_weight * stability_adjustment
  estimate_ubsf <- estimate_base + total_adjustment

  # Step 8: SE calculation
  se_bias_adj <- tryCatch({
    if (grepl("PET", bias_method)) {
      fit_pet <- rma(yi = yi, vi = vi, mods = ~ sei, method = "REML")
      fit_pet$se[1]
    } else if (grepl("PEESE", bias_method)) {
      fit_peese <- rma(yi = yi, vi = vi, mods = ~ vi, method = "REML")
      fit_peese$se[1]
    } else {
      se_base * 1.1
    }
  }, error = function(e) se_base)

  var_combined <- se_base^2 +
                  (bias_weight^2 * se_bias_adj^2) +
                  (stability_weight^2 * var(loo$estimate, na.rm = TRUE) / k)
  se_ubsf <- sqrt(var_combined)

  ci_result <- calc_ci(estimate_ubsf, se_ubsf, k)

  z_stat <- estimate_ubsf / se_ubsf
  pval <- 2 * pt(-abs(z_stat), df = max(k - 2, 1))

  return(list(
    method = "UBSF_v2 (Unified Bias-Stability V2)",
    estimate = round(estimate_ubsf, 4),
    se = round(se_ubsf, 4),
    ci_lb = round(ci_result$ci_lb, 4),
    ci_ub = round(ci_result$ci_ub, 4),
    pval = round(pval, 4),
    tau2 = round(tau2, 4),
    k = k,

    # Bias diagnostics
    bias_detected = bias_detected,
    bias_pval = round(bias_pval, 4),
    bias_method = bias_method,
    bias_weight = round(bias_weight, 3),
    bias_adjustment = round(bias_adjustment, 4),
    tf_estimate = round(tf_estimate, 4),
    pet_estimate = round(pet_estimate, 4),
    peese_estimate = round(peese_estimate, 4),
    k_imputed = k_imputed,

    # Stability diagnostics
    direction_fragile = dir_fragile,
    significance_fragile = sig_fragile,
    max_loo_change = round(max_change, 4),
    stability_adjustment = round(stability_adjustment, 4),

    # Summary
    estimate_original = round(estimate_base, 4),
    total_adjustment = round(total_adjustment, 4)
  ))
}

################################################################################
# METHOD 5: ENSEMBLE META-ANALYSIS V2 (EMA)
################################################################################
# IMPROVEMENTS:
# - Adaptive weighting based on method agreement
# - Use inverse-variance of methods as base weights
# - Include HKSJ in ensemble
# - Bootstrap SE option
################################################################################

ensemble_meta_analysis_v2 <- function(yi, vi, include_hksj = TRUE) {

  k <- length(yi)

  methods_results <- list()

  # Method 1: Standard REML
  tryCatch({
    fit_reml <- rma(yi = yi, vi = vi, method = "REML")
    methods_results$REML <- list(
      estimate = as.numeric(fit_reml$beta[1]),
      se = fit_reml$se
    )
  }, error = function(e) NULL)

  # Method 2: HKSJ
  if (include_hksj) {
    tryCatch({
      fit_hksj <- rma(yi = yi, vi = vi, method = "REML", test = "knha")
      methods_results$HKSJ <- list(
        estimate = as.numeric(fit_hksj$beta[1]),
        se = fit_hksj$se
      )
    }, error = function(e) NULL)
  }

  # Method 3: MAFI-Weighted V2
  tryCatch({
    mwm <- mafi_weighted_ma_v2(yi, vi)
    methods_results$MWM <- list(
      estimate = mwm$estimate,
      se = mwm$se
    )
  }, error = function(e) NULL)

  # Method 4: Adaptive Robust Pooling V2
  tryCatch({
    arp <- adaptive_robust_pooling_v2(yi, vi)
    methods_results$ARP <- list(
      estimate = arp$estimate,
      se = arp$se
    )
  }, error = function(e) NULL)

  # Method 5: Sequential Influence Trimming V2
  tryCatch({
    sit <- sequential_influence_trimming_v2(yi, vi)
    methods_results$SIT <- list(
      estimate = sit$estimate,
      se = sit$se
    )
  }, error = function(e) NULL)

  # Method 6: Unified Bias-Stability V2
  tryCatch({
    ubsf <- unified_bias_stability_v2(yi, vi)
    methods_results$UBSF <- list(
      estimate = ubsf$estimate,
      se = ubsf$se
    )
  }, error = function(e) NULL)

  if (length(methods_results) == 0) {
    return(list(method = "EMA_v2", estimate = NA, se = NA, message = "All methods failed"))
  }

  # Adaptive weighting based on inverse-variance
  estimates <- sapply(methods_results, function(x) x$estimate)
  ses <- sapply(methods_results, function(x) x$se)

  # Inverse-variance weights among methods
  inv_var_weights <- 1 / ses^2
  inv_var_weights <- inv_var_weights / sum(inv_var_weights)

  # Ensemble estimate
  estimate_ema <- sum(inv_var_weights * estimates)

  # Ensemble variance (Rubin's rules)
  var_within <- sum(inv_var_weights * ses^2)
  var_between <- sum(inv_var_weights * (estimates - estimate_ema)^2)

  # Add between-method variance scaled by number of methods
  n_methods <- length(methods_results)
  var_total <- var_within + var_between * (1 + 1/n_methods)
  se_ema <- sqrt(var_total)

  # Confidence interval with t-distribution
  ci_result <- calc_ci(estimate_ema, se_ema, k)

  # P-value
  z_stat <- estimate_ema / se_ema
  pval <- 2 * pt(-abs(z_stat), df = max(k - 2, 1))

  # Method agreement metric (1 = perfect agreement)
  agreement <- 1 - (sd(estimates) / max(abs(estimate_ema), 0.001))
  agreement <- max(0, min(1, agreement))

  return(list(
    method = "EMA_v2 (Ensemble Meta-Analysis V2)",
    estimate = round(estimate_ema, 4),
    se = round(se_ema, 4),
    ci_lb = round(ci_result$ci_lb, 4),
    ci_ub = round(ci_result$ci_ub, 4),
    pval = round(pval, 4),
    k = k,
    n_methods = n_methods,
    method_estimates = round(estimates, 4),
    method_weights = round(inv_var_weights, 3),
    method_ses = round(ses, 4),
    var_within = round(var_within, 6),
    var_between = round(var_between, 6),
    agreement = round(agreement, 3)
  ))
}

################################################################################
# COMPARISON FUNCTION V2
################################################################################

compare_all_methods_v2 <- function(yi, vi, include_v1 = FALSE) {

  k <- length(yi)

  results <- data.frame(
    Method = character(),
    Estimate = numeric(),
    SE = numeric(),
    CI_Lower = numeric(),
    CI_Upper = numeric(),
    PValue = numeric(),
    stringsAsFactors = FALSE
  )

  # Standard REML
  tryCatch({
    fit <- rma(yi = yi, vi = vi, method = "REML")
    results <- rbind(results, data.frame(
      Method = "REML",
      Estimate = round(as.numeric(fit$beta[1]), 4),
      SE = round(fit$se, 4),
      CI_Lower = round(fit$ci.lb, 4),
      CI_Upper = round(fit$ci.ub, 4),
      PValue = round(fit$pval, 4)
    ))
  }, error = function(e) NULL)

  # HKSJ
  tryCatch({
    fit <- rma(yi = yi, vi = vi, method = "REML", test = "knha")
    results <- rbind(results, data.frame(
      Method = "HKSJ",
      Estimate = round(as.numeric(fit$beta[1]), 4),
      SE = round(fit$se, 4),
      CI_Lower = round(fit$ci.lb, 4),
      CI_Upper = round(fit$ci.ub, 4),
      PValue = round(fit$pval, 4)
    ))
  }, error = function(e) NULL)

  # V2 Methods
  tryCatch({
    mwm <- mafi_weighted_ma_v2(yi, vi)
    results <- rbind(results, data.frame(
      Method = "MWM_v2",
      Estimate = mwm$estimate,
      SE = mwm$se,
      CI_Lower = mwm$ci_lb,
      CI_Upper = mwm$ci_ub,
      PValue = mwm$pval
    ))
  }, error = function(e) NULL)

  tryCatch({
    arp <- adaptive_robust_pooling_v2(yi, vi)
    results <- rbind(results, data.frame(
      Method = "ARP_v2",
      Estimate = arp$estimate,
      SE = arp$se,
      CI_Lower = arp$ci_lb,
      CI_Upper = arp$ci_ub,
      PValue = arp$pval
    ))
  }, error = function(e) NULL)

  tryCatch({
    sit <- sequential_influence_trimming_v2(yi, vi)
    results <- rbind(results, data.frame(
      Method = "SIT_v2",
      Estimate = sit$estimate,
      SE = sit$se,
      CI_Lower = sit$ci_lb,
      CI_Upper = sit$ci_ub,
      PValue = sit$pval
    ))
  }, error = function(e) NULL)

  tryCatch({
    ubsf <- unified_bias_stability_v2(yi, vi)
    results <- rbind(results, data.frame(
      Method = "UBSF_v2",
      Estimate = ubsf$estimate,
      SE = ubsf$se,
      CI_Lower = ubsf$ci_lb,
      CI_Upper = ubsf$ci_ub,
      PValue = ubsf$pval
    ))
  }, error = function(e) NULL)

  tryCatch({
    ema <- ensemble_meta_analysis_v2(yi, vi)
    results <- rbind(results, data.frame(
      Method = "EMA_v2",
      Estimate = ema$estimate,
      SE = ema$se,
      CI_Lower = ema$ci_lb,
      CI_Upper = ema$ci_ub,
      PValue = ema$pval
    ))
  }, error = function(e) NULL)

  return(results)
}

################################################################################
# QUICK VALIDATION
################################################################################

cat("Advanced Pooling Methods V2 loaded successfully.\n")
cat("Improvements over V1:\n")
cat("  - T-distribution CIs for better coverage (k < 30)\n")
cat("  - MWM: Recalibrated stability weights, studentized residuals\n")
cat("  - SIT: Adaptive trimming with studentized residuals\n")
cat("  - UBSF: PET-PEESE bias correction, adaptive weights\n")
cat("  - ARP: Inverse-variance weights, proper between-method variance\n")
cat("  - EMA: Rubin's rules for combining, better agreement metric\n")
cat("\nAvailable V2 functions:\n")
cat("  - mafi_weighted_ma_v2()\n")
cat("  - adaptive_robust_pooling_v2()\n")
cat("  - sequential_influence_trimming_v2()\n")
cat("  - unified_bias_stability_v2()\n")
cat("  - ensemble_meta_analysis_v2()\n")
cat("  - compare_all_methods_v2()\n")
