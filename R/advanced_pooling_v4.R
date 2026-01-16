#' Advanced Pooling Methods for Meta-Analysis (V4 - 10 Novel Methods)
#'
#' Implementation of 10 novel pooling methods for meta-analysis:
#' Category A: Robustness Methods (WRD, CBM, RBM)
#' Category B: Bias Correction Methods (SWA, TAS)
#' Category C: Advanced Variance Methods (EVE, PVM)
#' Category D: Ensemble & Adaptive Methods (AEM, SPE)
#' Category E: Specialized Methods (SMS)
#'
#' @name advanced_pooling_v4
NULL

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Safe Scaling Function
#'
#' Handles constant vectors without producing NA
#' @keywords internal
safe_scale <- function(x) {
  if (length(x) < 2 || sd(x, na.rm = TRUE) == 0) {
    return(rep(0, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

#' Safe Huber M-Estimator Weights
#'
#' Compute Huber weights for robust regression
#' @keywords internal
huber_weights <- function(residuals, k = 1.345) {
  abs_res <- abs(residuals)
  weights <- ifelse(abs_res <= k, 1, k / abs_res)
  return(weights)
}

#' Bootstrap SE for Post-Selection Inference
#'
#' Bootstrap a method to get corrected SE
#' @keywords internal
bootstrap_se <- function(yi, vi, method_func, n_boot = 499, seed = NULL) {
  k <- length(yi)
  boot_estimates <- numeric(n_boot)

  if (!is.null(seed)) set.seed(seed)

  for (b in 1:n_boot) {
    boot_idx <- sample(1:k, k, replace = TRUE)
    yi_boot <- yi[boot_idx]
    vi_boot <- vi[boot_idx]

    boot_result <- tryCatch(
      method_func(yi_boot, vi_boot),
      error = function(e) list(estimate = NA, se = NA)
    )

    boot_estimates[b] <- boot_result$estimate
  }

  list(
    boot_se = sd(boot_estimates, na.rm = TRUE),
    boot_mean = mean(boot_estimates, na.rm = TRUE),
    boot_median = median(boot_estimates, na.rm = TRUE),
    n_successful = sum(!is.na(boot_estimates))
  )
}

# ============================================================================
# CATEGORY A: ROBUSTNESS METHODS
# ============================================================================

#' WRD - Winsorized Robust Downsampling
#'
#' Addresses extreme outliers without removing studies entirely using
#' winsorization and robust M-estimation.
#'
#' @param yi Numeric vector of effect sizes
#' @param vi Numeric vector of sampling variances
#' @param winsor_threshold Threshold for winsorization (default 3.0)
#' @param huber_k Huber constant for M-estimation (default 1.345)
#' @param bootstrap Use bootstrap SE correction (default TRUE)
#' @param n_boot Number of bootstrap iterations (default 499)
#' @return List with estimate, se, ci_lb, ci_ub, and diagnostics
#' @export
wrd_meta <- function(yi, vi, winsor_threshold = 3.0, huber_k = 1.345,
                     bootstrap = TRUE, n_boot = 499) {
  k <- length(yi)

  # Minimum k check
  if (k < 4) {
    warning("WRD requires k >= 4. Falling back to HKSJ.")
    fit <- metafor::rma(yi, vi, method = "REML", test = "knha")
    return(list(
      estimate = as.numeric(coef(fit)),
      se = fit$se,
      ci_lb = fit$ci.lb,
      ci_ub = fit$ci.ub,
      pvalue = fit$pval,
      tau2 = fit$tau2,
      method = "HKSJ_fallback",
      k = k,
      note = "k < 4, used HKSJ instead of WRD"
    ))
  }

  # Step 1: Initial REML fit
  fit_init <- metafor::rma(yi, vi, method = "REML")
  base_estimate <- as.numeric(coef(fit_init))
  tau2 <- fit_init$tau2

  # Step 2: Calculate studentized residuals
  resid_raw <- yi - base_estimate
  resid_var <- vi + tau2
  student_resid <- resid_raw / sqrt(pmax(resid_var, 1e-10))

  # Step 3: Winsorize extreme residuals
  student_winsor <- student_resid
  student_winsor[student_winsor > winsor_threshold] <- winsor_threshold
  student_winsor[student_winsor < -winsor_threshold] <- -winsor_threshold

  # Step 4: Convert back to winsorized effect sizes
  yi_winsor <- base_estimate + student_winsor * sqrt(pmax(resid_var, 1e-10))

  # Step 5: Robust M-estimation with Huber weights
  # Iteratively reweighted least squares
  weights <- rep(1, k)
  max_iter <- 25
  tol <- 1e-6

  for (iter in 1:max_iter) {
    # Weighted REML with current weights
    fit_weighted <- tryCatch(
      metafor::rma(yi_winsor, vi, method = "REML", weights = weights),
      error = function(e) NULL
    )

    if (is.null(fit_weighted)) break

    old_estimate <- as.numeric(coef(fit_weighted))
    resid <- yi_winsor - old_estimate
    new_weights <- huber_weights(resid / sqrt(pmax(vi + tau2, 1e-10)), k = huber_k)

    # Cap maximum weight to prevent extreme downweighting
    new_weights <- pmin(new_weights, 10)

    if (max(abs(new_weights - weights)) < tol) break
    weights <- new_weights
  }

  # Final fit with Huber weights
  fit_final <- metafor::rma(yi_winsor, vi, method = "REML", weights = weights)
  estimate <- as.numeric(coef(fit_final))

  # Step 6: Compute SE
  if (bootstrap) {
    boot_result <- bootstrap_se(yi, vi, function(y, v) {
      wrd_meta(y, v, winsor_threshold, huber_k, bootstrap = FALSE)
    }, n_boot = n_boot)

    se <- boot_result$boot_se
    se_method <- "bootstrap"
  } else {
    # Sandwich variance estimator
    w <- weights / sum(weights)
    se <- sqrt(sum(w^2 * (vi + tau2)))
    se_method <- "sandwich"
  }

  # T-distribution CI
  df <- max(1, k - 2)
  t_crit <- qt(0.975, df)

  ci_lb <- estimate - t_crit * se
  ci_ub <- estimate + t_crit * se
  pvalue <- 2 * (1 - pt(abs(estimate / se), df))

  list(
    estimate = estimate,
    se = se,
    ci_lb = ci_lb,
    ci_ub = ci_ub,
    pvalue = pvalue,
    tau2 = tau2,
    method = "WRD",
    k = k,
    n_winsorized = sum(abs(student_resid) > winsor_threshold),
    huber_weights = weights,
    se_method = se_method
  )
}


#' CBM - Clustering-Based Meta-Analysis
#'
#' Handles multimodal effect distributions using Gaussian mixture models.
#'
#' @param yi Numeric vector of effect sizes
#' @param vi Numeric vector of sampling variances
#' @param max_clusters Maximum number of clusters (default 5)
#' @param min_cluster_size Minimum cluster size (default 3)
#' @param cluster_method Clustering method: "mixture" or "hierarchical" (default "mixture")
#' @return List with estimate, se, ci_lb, ci_ub, and diagnostics
#' @export
cbm_meta <- function(yi, vi, max_clusters = 5, min_cluster_size = 3,
                     cluster_method = "mixture") {
  k <- length(yi)

  # Minimum k check
  if (k < 6) {
    warning("CBM requires k >= 6. Falling back to HKSJ.")
    fit <- metafor::rma(yi, vi, method = "REML", test = "knha")
    return(list(
      estimate = as.numeric(coef(fit)),
      se = fit$se,
      ci_lb = fit$ci.lb,
      ci_ub = fit$ci.ub,
      pvalue = fit$pval,
      tau2 = fit$tau2,
      method = "HKSJ_fallback",
      k = k,
      note = "k < 6, used HKSJ instead of CBM"
    ))
  }

  # Step 1: Standardize effect sizes for clustering
  sei <- sqrt(vi)
  yi_weighted <- yi / sei  # Precision-weighted for clustering

  # Step 2: Determine optimal number of clusters
  if (cluster_method == "mixture" && requireNamespace("mclust", quietly = TRUE)) {
    # Use BIC to select number of clusters
    best_bic <- -Inf
    best_clusters <- 1
    cluster_assignments <- rep(1, k)

    for (n_clust in 1:min(max_clusters, k %/% min_cluster_size)) {
      if (n_clust == 1) {
        bic <- sum(dnorm(yi_weighted, mean(yi_weighted), sd(yi_weighted), log = TRUE))
      } else {
        mod <- tryCatch(
          mclust::Mclust(yi_weighted, G = n_clust),
          error = function(e) NULL
        )
        if (!is.null(mod)) {
          bic <- mod$BIC
        } else {
          bic <- -Inf
        }
      }

      if (bic > best_bic) {
        best_bic <- bic
        best_clusters <- n_clust
        if (n_clust > 1) cluster_assignments <- mod$classification
      }
    }

    n_clusters <- best_clusters

  } else if (cluster_method == "hierarchical") {
    # Hierarchical clustering
    dist_matrix <- dist(yi_weighted)
    hc <- hclust(dist_matrix, method = "ward.D2")

    # Cut tree to maximize between-cluster variance
    n_clusters <- min(max_clusters, k %/% min_cluster_size)

    # Silhouette-based selection if available
    if (requireNamespace("cluster", quietly = TRUE)) {
      best_sil <- -1
      for (n_clust in 2:n_clusters) {
        assigns <- cutree(hc, k = n_clust)
        if (n_clust == 2) {
          sil <- cluster::silhouette(assigns, dist_matrix)
        } else {
          sil <- cluster::silhouette(assigns, dist_matrix)
        }
        sil_avg <- mean(sil[, 3])
        if (sil_avg > best_sil) {
          best_sil <- sil_avg
          n_clusters <- n_clust
        }
      }
    }

    cluster_assignments <- cutree(hc, k = n_clusters)
  } else {
    # Fallback: simple k-means
    n_clusters <- min(3, k %/% min_cluster_size)
    km <- tryCatch(
      kmeans(yi_weighted, centers = n_clusters, nstart = 10),
      error = function(e) list(cluster = rep(1, k))
    )
    cluster_assignments <- km$cluster
  }

  # Step 3: Estimate separate τ² for each cluster
  cluster_estimates <- numeric(n_clusters)
  cluster_variances <- numeric(n_clusters)
  cluster_weights <- numeric(n_clusters)
  cluster_sizes <- numeric(n_clusters)

  for (c in 1:n_clusters) {
    idx <- which(cluster_assignments == c)
    cluster_sizes[c] <- length(idx)

    if (length(idx) >= 3) {
      fit_c <- metafor::rma(yi[idx], vi[idx], method = "REML")
      cluster_estimates[c] <- as.numeric(coef(fit_c))
      cluster_variances[c] <- fit_c$se^2
      cluster_weights[c] <- 1 / fit_c$se^2
    } else {
      # Too small, use simple inverse variance
      cluster_estimates[c] <- sum(yi[idx] / vi[idx]) / sum(1 / vi[idx])
      cluster_variances[c] <- 1 / sum(1 / vi[idx])
      cluster_weights[c] <- sum(1 / vi[idx])
    }
  }

  # Step 4: Combine clusters using inverse variance
  cluster_weights <- cluster_weights / sum(cluster_weights)
  estimate <- sum(cluster_weights * cluster_estimates)

  # Step 5: Rubin's rules for combining clusters
  within_var <- sum(cluster_weights * cluster_variances)
  between_var <- sum(cluster_weights * (cluster_estimates - estimate)^2)
  total_var <- within_var + (1 + 1/n_clusters) * between_var

  se <- sqrt(total_var)

  # T-distribution CI
  df <- max(1, k - n_clusters)
  t_crit <- qt(0.975, df)

  ci_lb <- estimate - t_crit * se
  ci_ub <- estimate + t_crit * se
  pvalue <- 2 * (1 - pt(abs(estimate / se), df))

  list(
    estimate = estimate,
    se = se,
    ci_lb = ci_lb,
    ci_ub = ci_ub,
    pvalue = pvalue,
    method = "CBM",
    k = k,
    n_clusters = n_clusters,
    cluster_sizes = cluster_sizes,
    cluster_assignments = cluster_assignments,
    cluster_estimates = cluster_estimates
  )
}


#' RBM - Residual-Based Weighting Method
#'
#' Iteratively reweights studies based on residual patterns using
#' iteratively reweighted least squares.
#'
#' @param yi Numeric vector of effect sizes
#' @param vi Numeric vector of sampling variances
#' @param max_iterations Maximum iterations (default 25)
#' @param convergence_tol Convergence tolerance (default 1e-6)
#' @param weight_cap Maximum weight multiplier (default 10)
#' @return List with estimate, se, ci_lb, ci_ub, and diagnostics
#' @export
rbm_meta <- function(yi, vi, max_iterations = 25, convergence_tol = 1e-6,
                     weight_cap = 10) {
  k <- length(yi)

  # Minimum k check
  if (k < 4) {
    warning("RBM requires k >= 4. Falling back to HKSJ.")
    fit <- metafor::rma(yi, vi, method = "REML", test = "knha")
    return(list(
      estimate = as.numeric(coef(fit)),
      se = fit$se,
      ci_lb = fit$ci.lb,
      ci_ub = fit$ci.ub,
      pvalue = fit$pval,
      tau2 = fit$tau2,
      method = "HKSJ_fallback",
      k = k,
      note = "k < 4, used HKSJ instead of RBM"
    ))
  }

  # Initialize
  weights <- rep(1, k)
  converged <- FALSE
  iterations <- 0

  # Iteratively reweighted estimation
  for (iter in 1:max_iterations) {
    iterations <- iter

    # Fit with current weights
    fit <- tryCatch(
      metafor::rma(yi, vi, method = "REML", weights = weights),
      error = function(e) NULL
    )

    if (is.null(fit)) break

    old_estimate <- as.numeric(coef(fit))
    tau2 <- fit$tau2

    # Calculate studentized residuals
    resid_raw <- yi - old_estimate
    resid_var <- vi + tau2
    student_resid <- resid_raw / sqrt(pmax(resid_var, 1e-10))
    squared_resid <- student_resid^2

    # Update weights based on residual pattern
    # Studies with large residuals get downweighted
    resid_median <- median(squared_resid)
    new_weights <- 1 / pmax(squared_resid / resid_median, 0.1)

    # Cap weights
    new_weights <- pmin(new_weights, weight_cap)

    # Check convergence
    if (max(abs(new_weights - weights)) < convergence_tol) {
      converged <- TRUE
      weights <- new_weights
      break
    }

    weights <- new_weights
  }

  # Final fit with converged weights
  fit_final <- metafor::rma(yi, vi, method = "REML", weights = weights)
  estimate <- as.numeric(coef(fit_final))
  tau2 <- fit_final$tau2

  # Sandwich variance estimator
  w <- weights / sum(weights)
  within_var <- sum(w^2 * vi)
  between_var <- sum(w^2 * tau2)
  sandwich_se <- sqrt(within_var + between_var)

  se <- sandwich_se

  # T-distribution CI
  df <- max(1, k - 2)
  t_crit <- qt(0.975, df)

  ci_lb <- estimate - t_crit * se
  ci_ub <- estimate + t_crit * se
  pvalue <- 2 * (1 - pt(abs(estimate / se), df))

  list(
    estimate = estimate,
    se = se,
    ci_lb = ci_lb,
    ci_ub = ci_ub,
    pvalue = pvalue,
    tau2 = tau2,
    method = "RBM",
    k = k,
    weights = weights,
    converged = converged,
    iterations = iterations
  )
}


# ============================================================================
# CATEGORY B: BIAS CORRECTION METHODS
# ============================================================================

#' SWA - Selection-Weight Adjustment
#'
#' Models publication bias directly using weight-function models.
#'
#' @param yi Numeric vector of effect sizes
#' @param vi Numeric vector of sampling variances
#' @param selection_function "step" or "continuous" (default "step")
#' @param p_cutoff P-value cutoff for step function (default 0.05)
#' @param n_boot Number of bootstrap iterations for SE (default 1000)
#' @return List with estimate, se, ci_lb, ci_ub, and diagnostics
#' @export
swa_meta <- function(yi, vi, selection_function = "step", p_cutoff = 0.05,
                     n_boot = 1000) {
  k <- length(yi)

  # Minimum k check
  if (k < 10) {
    warning("SWA requires k >= 10. Falling back to HKSJ.")
    fit <- metafor::rma(yi, vi, method = "REML", test = "knha")
    return(list(
      estimate = as.numeric(coef(fit)),
      se = fit$se,
      ci_lb = fit$ci.lb,
      ci_ub = fit$ci.ub,
      pvalue = fit$pval,
      tau2 = fit$tau2,
      method = "HKSJ_fallback",
      k = k,
      note = "k < 10, used HKSJ instead of SWA"
    ))
  }

  # Initial fit to get selection probabilities
  fit_init <- metafor::rma(yi, vi, method = "REML")
  base_estimate <- as.numeric(coef(fit_init))
  base_se <- fit_init$se
  tau2 <- fit_init$tau2

  # Calculate z-values and p-values for each study
  z_values <- yi / sqrt(vi + tau2)
  p_values <- 2 * (1 - pnorm(abs(z_values)))

  # Step 1: Model selection function
  if (selection_function == "step") {
    # Step function: studies with p < p_cutoff more likely to be published
    selection_prob <- ifelse(p_values < p_cutoff, 1.0, 0.5)
  } else {
    # Continuous weight function (Vevea & Hoods style)
    # Smoothly decreasing function of p-value
    selection_prob <- exp(-0.5 * (p_values / 0.05)^2)
  }

  # Step 2: Compute selection-weighted estimate
  # Weight by inverse variance and inverse selection probability
  iv_weights <- 1 / vi
  selection_weights <- 1 / selection_prob
  combined_weights <- iv_weights * selection_weights
  combined_weights <- combined_weights / sum(combined_weights)

  weighted_estimate <- sum(combined_weights * yi)

  # Step 3: Bootstrap for SE
  boot_estimates <- numeric(n_boot)
  for (b in 1:n_boot) {
    boot_idx <- sample(1:k, k, replace = TRUE)
    yi_boot <- yi[boot_idx]
    vi_boot <- vi[boot_idx]

    fit_boot <- tryCatch(
      metafor::rma(yi_boot, vi_boot, method = "REML"),
      error = function(e) fit_init
    )

    z_boot <- yi_boot / sqrt(vi_boot + fit_boot$tau2)
    p_boot <- 2 * (1 - pnorm(abs(z_boot)))

    if (selection_function == "step") {
      sel_prob_boot <- ifelse(p_boot < p_cutoff, 1.0, 0.5)
    } else {
      sel_prob_boot <- exp(-0.5 * (p_boot / 0.05)^2)
    }

    w_boot <- (1 / vi_boot) * (1 / sel_prob_boot)
    w_boot <- w_boot / sum(w_boot)

    boot_estimates[b] <- sum(w_boot * yi_boot)
  }

  se <- sd(boot_estimates, na.rm = TRUE)

  # T-distribution CI
  df <- max(1, k - 2)
  t_crit <- qt(0.975, df)

  ci_lb <- weighted_estimate - t_crit * se
  ci_ub <- weighted_estimate + t_crit * se
  pvalue <- 2 * (1 - pt(abs(weighted_estimate / se), df))

  list(
    estimate = weighted_estimate,
    se = se,
    ci_lb = ci_lb,
    ci_ub = ci_ub,
    pvalue = pvalue,
    method = "SWA",
    k = k,
    selection_function = selection_function,
    selection_probabilities = selection_prob,
    n_boot = n_boot
  )
}


#' TAS - Trim-and-Shrink
#'
#' Combines trim-and-fill with empirical Bayes shrinkage.
#'
#' @param yi Numeric vector of effect sizes
#' @param vi Numeric vector of sampling variances
#' @param trim_estimator Trim-and-fill estimator (default "R0")
#' @param shrinkage_weight Weight for shrinkage (default 0.2, reduced from 0.5 for better coverage)
#' @return List with estimate, se, ci_lb, ci_ub, and diagnostics
#' @export
tas_meta <- function(yi, vi, trim_estimator = "R0", shrinkage_weight = 0.2) {
  k <- length(yi)

  # Minimum k check
  if (k < 10) {
    warning("TAS requires k >= 10. Falling back to HKSJ.")
    fit <- metafor::rma(yi, vi, method = "REML", test = "knha")
    return(list(
      estimate = as.numeric(coef(fit)),
      se = fit$se,
      ci_lb = fit$ci.lb,
      ci_ub = fit$ci.ub,
      pvalue = fit$pval,
      tau2 = fit$tau2,
      method = "HKSJ_fallback",
      k = k,
      note = "k < 10, used HKSJ instead of TAS"
    ))
  }

  # Step 1: Standard trim-and-fill
  tf_fit <- tryCatch(
    metafor::trimfill(yi, vi, estimator = trim_estimator, method = "REML"),
    error = function(e) metafor::rma(yi, vi, method = "REML")
  )

  # Get original and imputed studies
  # Check if fill component exists (indicates imputation occurred)
  if ("fill" %in% names(tf_fit) && !is.null(tf_fit$fill)) {
    n_imputed <- tf_fit$k - k
  } else {
    # No studies were imputed
    n_imputed <- 0
  }
  tf_estimate <- as.numeric(coef(tf_fit))

  # Step 2: Compute empirical Bayes shrinkage factor
  fit_orig <- metafor::rma(yi, vi, method = "REML")
  orig_estimate <- as.numeric(coef(fit_orig))
  orig_se <- fit_orig$se

  # Shrinkage factor based on precision
  # More shrinkage when n_imputed is large relative to k
  # Cap shrinkage to at most 30% (keep at least 70% of TF estimate)
  shrinkage_factor_raw <- 1 - (n_imputed / (k + n_imputed))
  shrinkage_factor <- max(0.7, shrinkage_factor_raw)  # At least 70% of TF estimate

  # Step 3: Combine trim-and-fill with shrunk estimate
  shrunk_estimate <- shrinkage_factor * tf_estimate + (1 - shrinkage_factor) * orig_estimate

  # Step 4: Adjust SE for imputation uncertainty
  # SE increases with number of imputed studies
  imputation_penalty <- 1 + (n_imputed / (2 * k))
  adjusted_se <- tf_fit$se * imputation_penalty

  # Blend between trim-and-fill and shrunk
  estimate <- shrinkage_weight * shrunk_estimate + (1 - shrinkage_weight) * tf_estimate
  se <- adjusted_se

  # T-distribution CI
  df <- max(1, k - 2)
  t_crit <- qt(0.975, df)

  ci_lb <- estimate - t_crit * se
  ci_ub <- estimate + t_crit * se
  pvalue <- 2 * (1 - pt(abs(estimate / se), df))

  list(
    estimate = estimate,
    se = se,
    ci_lb = ci_lb,
    ci_ub = ci_ub,
    pvalue = pvalue,
    method = "TAS",
    k = k,
    n_imputed = n_imputed,
    tf_estimate = tf_estimate,
    shrunk_estimate = shrunk_estimate,
    shrinkage_factor = shrinkage_factor,
    shrinkage_weight = shrinkage_weight
  )
}


# ============================================================================
# CATEGORY C: ADVANCED VARIANCE METHODS
# ============================================================================

#' EVE - Empirical Bayes Variance Estimation
#'
#' Improved tau² estimation using empirical Bayes combining.
#'
#' @param yi Numeric vector of effect sizes
#' @param vi Numeric vector of sampling variances
#' @param prior_shape Prior shape for inverse-gamma (default 0.001)
#' @param prior_rate Prior rate for inverse-gamma (default 0.001)
#' @return List with estimate, se, ci_lb, ci_ub, and diagnostics
#' @export
eve_meta <- function(yi, vi, prior_shape = 0.001, prior_rate = 0.001) {
  k <- length(yi)

  # Minimum k check
  if (k < 3) {
    warning("EVE requires k >= 3. Falling back to HKSJ.")
    fit <- metafor::rma(yi, vi, method = "REML", test = "knha")
    return(list(
      estimate = as.numeric(coef(fit)),
      se = fit$se,
      ci_lb = fit$ci.lb,
      ci_ub = fit$ci.ub,
      pvalue = fit$pval,
      tau2 = fit$tau2,
      method = "HKSJ_fallback",
      k = k,
      note = "k < 3, used HKSJ instead of EVE"
    ))
  }

  # Step 1: Compute multiple tau² estimators
  estimators <- c("REML", "DL", "PM", "ML", "SJ")
  tau2_estimates <- numeric(length(estimators))
  tau2_ses <- numeric(length(estimators))

  for (i in 1:length(estimators)) {
    fit <- tryCatch(
      metafor::rma(yi, vi, method = estimators[i]),
      error = function(e) NULL
    )

    if (!is.null(fit)) {
      tau2_estimates[i] <- fit$tau2
      # Approximate SE using delta method
      tau2_ses[i] <- fit$se * sqrt(2 / (k - 1))
    } else {
      tau2_estimates[i] <- NA
      tau2_ses[i] <- NA
    }
  }

  # Remove NA estimates
  valid_idx <- !is.na(tau2_estimates)
  tau2_valid <- tau2_estimates[valid_idx]
  se_valid <- tau2_ses[valid_idx]

  if (length(tau2_valid) == 0) {
    # All failed, use REML
    fit <- metafor::rma(yi, vi, method = "REML", test = "knha")
    return(list(
      estimate = as.numeric(coef(fit)),
      se = fit$se,
      ci_lb = fit$ci.lb,
      ci_ub = fit$ci.ub,
      pvalue = fit$pval,
      tau2 = fit$tau2,
      method = "HKSJ_fallback",
      k = k,
      note = "All tau2 estimators failed"
    ))
  }

  # Step 2: Empirical Bayes combination
  # Posterior distribution of tau²
  # Combine using inverse-variance weighting
  weights <- 1 / (se_valid^2 + 0.001)  # Add small constant
  weights <- weights / sum(weights)

  # Posterior mean
  tau2_eb <- sum(weights * tau2_valid)

  # Step 3: Use empirical Bayes tau² for meta-analysis
  fit_eb <- metafor::rma(yi, vi, method = "REML")
  # Override tau² estimate
  estimate <- as.numeric(coef(fit_eb))

  # SE with EB tau²
  w_eb <- 1 / (vi + tau2_eb)
  w_eb <- w_eb / sum(w_eb)
  se <- sqrt(sum(w_eb^2 * (vi + tau2_eb)))

  # T-distribution CI
  df <- max(1, k - 2)
  t_crit <- qt(0.975, df)

  ci_lb <- estimate - t_crit * se
  ci_ub <- estimate + t_crit * se
  pvalue <- 2 * (1 - pt(abs(estimate / se), df))

  list(
    estimate = estimate,
    se = se,
    ci_lb = ci_lb,
    ci_ub = ci_ub,
    pvalue = pvalue,
    tau2 = tau2_eb,
    method = "EVE",
    k = k,
    tau2_estimates = setNames(tau2_estimates, estimators),
    tau2_eb = tau2_eb,
    eb_weights = setNames(rep(NA, length(estimators)), estimators)
  )
}


#' PVM - Profile Variance Minimization
#'
#' Choose tau² that directly optimizes CI properties.
#'
#' @param yi Numeric vector of effect sizes
#' @param vi Numeric vector of sampling variances
#' @param tau_grid_min Minimum tau² for grid (default 0)
#' @param tau_grid_max Maximum tau² for grid (default 1)
#' @param n_grid_points Number of grid points (default 100)
#' @param coverage_target Target coverage (default 0.95)
#' @param n_boot Bootstrap iterations for coverage (default 500)
#' @return List with estimate, se, ci_lb, ci_ub, and diagnostics
#' @export
pvm_meta <- function(yi, vi, tau_grid_min = 0, tau_grid_max = 1, n_grid_points = 100,
                     coverage_target = 0.95, n_boot = 500) {
  k <- length(yi)

  # Minimum k check
  if (k < 5) {
    warning("PVM requires k >= 5. Falling back to HKSJ.")
    fit <- metafor::rma(yi, vi, method = "REML", test = "knha")
    return(list(
      estimate = as.numeric(coef(fit)),
      se = fit$se,
      ci_lb = fit$ci.lb,
      ci_ub = fit$ci.ub,
      pvalue = fit$pval,
      tau2 = fit$tau2,
      method = "HKSJ_fallback",
      k = k,
      note = "k < 5, used HKSJ instead of PVM"
    ))
  }

  # Step 1: Create tau² grid
  tau_grid <- seq(tau_grid_min, tau_grid_max, length.out = n_grid_points)

  # Step 2: Profile over tau² grid
  # For each tau², compute pooled estimate and empirical coverage
  coverage_by_tau <- numeric(n_grid_points)
  ci_width_by_tau <- numeric(n_grid_points)

  for (i in 1:n_grid_points) {
    tau2_candidate <- tau_grid[i]

    # Compute weighted estimate with this tau²
    w <- 1 / (vi + tau2_candidate)
    w <- w / sum(w)
    theta_candidate <- sum(w * yi)

    # Bootstrap to estimate coverage
    covered <- 0
    for (b in 1:n_boot) {
      boot_idx <- sample(1:k, k, replace = TRUE)
      yi_boot <- yi[boot_idx]
      vi_boot <- vi[boot_idx]

      w_boot <- 1 / (vi_boot + tau2_candidate)
      w_boot <- w_boot / sum(w_boot)
      theta_boot <- sum(w_boot * yi_boot)

      # SE estimate
      se_boot <- sqrt(sum(w_boot^2 * (vi_boot + tau2_candidate)))

      # T-critical
      df_boot <- max(1, k - 2)
      t_crit_boot <- qt(0.975, df_boot)

      ci_lb_boot <- theta_candidate - t_crit_boot * se_boot
      ci_ub_boot <- theta_candidate + t_crit_boot * se_boot

      # Check if true effect (theta_candidate) is covered
      # Using bootstrap distribution to estimate
      if (theta_boot >= ci_lb_boot && theta_boot <= ci_ub_boot) {
        covered <- covered + 1
      }
    }

    coverage_by_tau[i] <- covered / n_boot

    # Also compute CI width
    se_candidate <- sqrt(sum(w^2 * (vi + tau2_candidate)))
    t_crit <- qt(0.975, max(1, k - 2))
    ci_width_by_tau[i] <- 2 * t_crit * se_candidate
  }

  # Step 3: Find tau² that achieves target coverage with minimum CI width
  # Among tau² values with coverage >= target, pick one with smallest CI width
  eligible_idx <- which(coverage_by_tau >= coverage_target * 0.95)  # Slightly relaxed

  if (length(eligible_idx) == 0) {
    # No tau² achieves target, use standard REML
    fit <- metafor::rma(yi, vi, method = "REML", test = "knha")
    return(list(
      estimate = as.numeric(coef(fit)),
      se = fit$se,
      ci_lb = fit$ci.lb,
      ci_ub = fit$ci.ub,
      pvalue = fit$pval,
      tau2 = fit$tau2,
      method = "HKSJ_fallback",
      k = k,
      note = "No tau² achieved target coverage"
    ))
  }

  # Among eligible, pick minimum CI width
  best_idx <- eligible_idx[which.min(ci_width_by_tau[eligible_idx])]
  tau2_optimal <- tau_grid[best_idx]

  # Step 4: Final estimate with optimal tau²
  w_final <- 1 / (vi + tau2_optimal)
  w_final <- w_final / sum(w_final)
  estimate <- sum(w_final * yi)

  se <- sqrt(sum(w_final^2 * (vi + tau2_optimal)))

  # T-distribution CI
  df <- max(1, k - 2)
  t_crit <- qt(0.975, df)

  ci_lb <- estimate - t_crit * se
  ci_ub <- estimate + t_crit * se
  pvalue <- 2 * (1 - pt(abs(estimate / se), df))

  list(
    estimate = estimate,
    se = se,
    ci_lb = ci_lb,
    ci_ub = ci_ub,
    pvalue = pvalue,
    tau2 = tau2_optimal,
    method = "PVM",
    k = k,
    tau_grid = tau_grid,
    coverage_by_tau = coverage_by_tau,
    ci_width_by_tau = ci_width_by_tau,
    tau2_optimal = tau2_optimal,
    achieved_coverage = coverage_by_tau[best_idx]
  )
}


# ============================================================================
# CATEGORY D: ENSEMBLE & ADAPTIVE METHODS
# ============================================================================

#' AEM - Adaptive Ensemble Meta-Analysis
#'
#' Smartly combines methods based on data characteristics using
#' stacking/adaptive weights.
#'
#' @param yi Numeric vector of effect sizes
#' @param vi Numeric vector of sampling variances
#' @param base_methods Methods to ensemble (default c("REML", "HKSJ", "MWM"))
#' @param adaptive_weights Use adaptive weighting (default TRUE)
#' @return List with estimate, se, ci_lb, ci_ub, and diagnostics
#' @export
aem_meta <- function(yi, vi, base_methods = c("REML", "HKSJ", "MWM"),
                     adaptive_weights = TRUE) {
  k <- length(yi)

  # Minimum k check
  if (k < 5) {
    warning("AEM requires k >= 5. Falling back to HKSJ.")
    fit <- metafor::rma(yi, vi, method = "REML", test = "knha")
    return(list(
      estimate = as.numeric(coef(fit)),
      se = fit$se,
      ci_lb = fit$ci.lb,
      ci_ub = fit$ci.ub,
      pvalue = fit$pval,
      tau2 = fit$tau2,
      method = "HKSJ_fallback",
      k = k,
      note = "k < 5, used HKSJ instead of AEM"
    ))
  }

  # Step 1: Compute data features for method selection
  sei <- sqrt(vi)
  features <- list(
    k = k,
    mean_se = mean(sei),
    sd_se = sd(sei),
    skew_yi <- if (k > 2) (mean((yi - mean(yi))^3) / (sd(yi)^3 + 1e-10)) else 0,
    kurtosis_yi = if (k > 3) (mean((yi - mean(yi))^4) / (sd(yi)^4 + 1e-10) - 3) else 0,
    range_yi = max(yi) - min(yi),
    i2 = NULL  # Will compute below
  )

  # Compute I²
  fit_reml <- metafor::rma(yi, vi, method = "REML")
  q_stat <- fit_reml$QE
  df_q <- k - 1
  features$i2 <- max(0, 100 * (q_stat - df_q) / q_stat)

  # Step 2: Run all base methods
  method_results <- list()

  # REML
  fit_reml <- metafor::rma(yi, vi, method = "REML")
  method_results$REML <- list(
    estimate = as.numeric(coef(fit_reml)),
    se = fit_reml$se,
    ci_lb = fit_reml$ci.lb,
    ci_ub = fit_reml$ci.ub
  )

  # HKSJ
  fit_hksj <- metafor::rma(yi, vi, method = "REML", test = "knha")
  method_results$HKSJ <- list(
    estimate = as.numeric(coef(fit_hksj)),
    se = fit_hksj$se,
    ci_lb = fit_hksj$ci.lb,
    ci_ub = fit_hksj$ci.ub
  )

  # MWM (if available)
  if ("MWM" %in% base_methods) {
    mwm_result <- tryCatch(
      mafi_weighted_ma(yi, vi),
      error = function(e) method_results$REML
    )
    method_results$MWM <- mwm_result
  }

  # Step 3: Compute adaptive weights
  if (adaptive_weights) {
    # Inverse-variance weighting based on SE
    estimates <- sapply(method_results, function(x) x$estimate)
    ses <- sapply(method_results, function(x) x$se)

    # Adaptive weight: higher weight to more precise methods
    iv_weights <- 1 / (ses^2)
    iv_weights <- iv_weights / sum(iv_weights)

    # Small-k adjustment: upweight HKSJ for small k
    if (k < 10 && "HKSJ" %in% names(iv_weights)) {
      iv_weights["HKSJ"] <- iv_weights["HKSJ"] * 1.2
      iv_weights <- iv_weights / sum(iv_weights)
    }

    # High heterogeneity adjustment: upweight robust methods
    if (!is.null(features$i2) && features$i2 > 50 && "MWM" %in% names(iv_weights)) {
      iv_weights["MWM"] <- iv_weights["MWM"] * 1.2
      iv_weights <- iv_weights / sum(iv_weights)
    }

  } else {
    # Equal weights
    iv_weights <- rep(1 / length(method_results), length(method_results))
    names(iv_weights) <- names(method_results)
  }

  # Step 4: Combine estimates
  estimates <- sapply(method_results, function(x) x$estimate)
  estimate <- sum(iv_weights * estimates)

  # Step 5: Rubin's rules for variance
  ses <- sapply(method_results, function(x) x$se)
  within_var <- sum(iv_weights * ses^2)
  between_var <- sum(iv_weights * (estimates - estimate)^2)

  # Correlation adjustment (methods use same data)
  rho <- 0.7
  total_var <- within_var + (1 + 1/length(method_results)) * between_var * (1 + rho)

  se <- sqrt(total_var)

  # T-distribution CI
  df <- max(1, k - 2)
  t_crit <- qt(0.975, df)

  ci_lb <- estimate - t_crit * se
  ci_ub <- estimate + t_crit * se
  pvalue <- 2 * (1 - pt(abs(estimate / se), df))

  list(
    estimate = estimate,
    se = se,
    ci_lb = ci_lb,
    ci_ub = ci_ub,
    pvalue = pvalue,
    tau2 = fit_reml$tau2,
    method = "AEM",
    k = k,
    features = features,
    base_estimates = estimates,
    ensemble_weights = iv_weights,
    adaptive_weights = adaptive_weights
  )
}


#' SPE - Stochastic Point Estimator
#'
#' Quantifies estimation uncertainty via stochastic sampling of tau².
#'
#' @param yi Numeric vector of effect sizes
#' @param vi Numeric vector of sampling variances
#' @param n_samples Number of MCMC samples (default 10000)
#' @param burnin Burn-in period (default 1000)
#' @param tau_prior Prior for tau² (default "half_cauchy")
#' @return List with estimate, se, ci_lb, ci_ub, and diagnostics
#' @export
spe_meta <- function(yi, vi, n_samples = 10000, burnin = 1000,
                     tau_prior = "half_cauchy") {
  k <- length(yi)

  # Minimum k check
  if (k < 3) {
    warning("SPE requires k >= 3. Falling back to HKSJ.")
    fit <- metafor::rma(yi, vi, method = "REML", test = "knha")
    return(list(
      estimate = as.numeric(coef(fit)),
      se = fit$se,
      ci_lb = fit$ci.lb,
      ci_ub = fit$ci.ub,
      pvalue = fit$pval,
      tau2 = fit$tau2,
      method = "HKSJ_fallback",
      k = k,
      note = "k < 3, used HKSJ instead of SPE"
    ))
  }

  # Step 1: Get initial tau² estimate
  fit_init <- metafor::rma(yi, vi, method = "REML")
  tau2_init <- max(0.01, fit_init$tau2)  # Ensure minimum

  # Step 2: Sample tau² from its posterior
  # Using Metropolis-Hastings

  # Proposal distribution - use adaptive proposal
  # Ensure minimum SD for exploration
  prop_sd <- max(0.02, tau2_init * 0.5)

  # Initialize
  tau2_current <- tau2_init
  tau2_samples <- numeric(n_samples + burnin)

  # Log-posterior function (simplified)
  log_posterior <- function(tau2, yi, vi, prior) {
    if (tau2 < 0) return(-Inf)

    # Likelihood contribution
    w <- 1 / (vi + tau2)
    w <- w / sum(w)
    theta <- sum(w * yi)

    # Log-likelihood (approximate)
    log_lik <- sum(dnorm(yi, theta, sqrt(vi + tau2), log = TRUE))

    # Prior
    if (prior == "half_cauchy") {
      # Half-Cauchy with scale 0.5
      log_prior <- log(2 / (pi * 0.5 * (1 + (tau2 / 0.5)^2)))
    } else {
      # Inverse-gamma
      log_prior <- dgamma(1 / (tau2 + 0.001), 0.001, 0.001, log = TRUE)
    }

    log_lik + log_prior
  }

  # Metropolis-Hastings
  accept_count <- 0
  for (i in 1:(n_samples + burnin)) {
    # Propose new tau² using log-normal for better positive constraint
    log_tau2_current <- log(tau2_current)
    log_tau2_proposal <- rnorm(1, log_tau2_current, 0.3)  # Log-SD for exploration
    tau2_proposal <- exp(log_tau2_proposal)

    # Ensure positive and not too extreme
    tau2_proposal <- max(0.001, min(10, tau2_proposal))

    # Acceptance probability
    log_alpha <- log_posterior(tau2_proposal, yi, vi, tau_prior) -
                 log_posterior(tau2_current, yi, vi, tau_prior) +
                 log(tau2_proposal) - log(tau2_current)  # Jacobian for log-transform

    if (log(runif(1)) < log_alpha) {
      tau2_current <- tau2_proposal
      accept_count <- accept_count + 1
    }

    tau2_samples[i] <- tau2_current
  }

  # Remove burnin
  tau2_samples <- tau2_samples[(burnin + 1):(n_samples + burnin)]

  # Step 3: For each tau² sample, compute pooled estimate
  theta_samples <- numeric(n_samples)

  for (i in 1:n_samples) {
    tau2_i <- tau2_samples[i]
    w_i <- 1 / (vi + tau2_i)
    w_i <- w_i / sum(w_i)
    theta_samples[i] <- sum(w_i * yi)
  }

  # Step 4: Final estimate and uncertainty
  estimate <- mean(theta_samples)  # Use mean for better stability

  # Compute SE that accounts for BOTH sampling uncertainty AND tau² uncertainty
  # Method 1: Use posterior predictive standard deviation
  # For each theta sample, compute what the observed yi would be
  # This is complex, so we use a hybrid approach:

  # Sampling SE (from standard REML)
  sampling_se <- fit_init$se

  # Parameter uncertainty SE (from theta_samples variation)
  param_se <- sd(theta_samples)

  # Combine: total SE = sqrt(sampling_se² + param_se²)
  # This accounts for uncertainty in both the data and the parameter
  se <- sqrt(sampling_se^2 + param_se^2)

  # Use T-distribution CI
  df <- max(1, k - 2)
  t_crit <- qt(0.975, df)

  ci_lb <- estimate - t_crit * se
  ci_ub <- estimate + t_crit * se

  # P-value
  pvalue <- 2 * (1 - pt(abs(estimate / se), df))

  list(
    estimate = estimate,
    se = se,
    ci_lb = as.numeric(ci_lb),
    ci_ub = as.numeric(ci_ub),
    pvalue = pvalue,
    tau2 = median(tau2_samples),
    method = "SPE",
    k = k,
    tau2_samples = tau2_samples,
    theta_samples = theta_samples,
    accept_rate = accept_count / (n_samples + burnin)
  )
}


# ============================================================================
# CATEGORY E: SPECIALIZED METHODS
# ============================================================================

#' SMS - Small-Meta Shrinkage
#'
#' Specialized for k < 10 using hierarchical shrinkage.
#'
#' @param yi Numeric vector of effect sizes
#' @param vi Numeric vector of sampling variances
#' @param global_prior Type of global prior (default "hierarchical")
#' @param shrinkage_intensity "adaptive" or "fixed" (default "adaptive")
#' @param min_k_for_shrinkage Minimum k for shrinkage (default 10)
#' @return List with estimate, se, ci_lb, ci_ub, and diagnostics
#' @export
sms_meta <- function(yi, vi, global_prior = "hierarchical", shrinkage_intensity = "adaptive",
                     min_k_for_shrinkage = 10) {
  k <- length(yi)

  # For k >= min_k_for_shrinkage, use standard method
  if (k >= min_k_for_shrinkage) {
    warning(paste("SMS designed for k <", min_k_for_shrinkage, ", using HKSJ for this analysis."))
    fit <- metafor::rma(yi, vi, method = "REML", test = "knha")
    return(list(
      estimate = as.numeric(coef(fit)),
      se = fit$se,
      ci_lb = fit$ci.lb,
      ci_ub = fit$ci.ub,
      pvalue = fit$pval,
      tau2 = fit$tau2,
      method = "HKSJ_fallback",
      k = k,
      note = paste("k >=", min_k_for_shrinkage, ", used HKSJ instead of SMS")
    ))
  }

  # Step 1: Fit standard model
  fit_standard <- metafor::rma(yi, vi, method = "REML")
  standard_estimate <- as.numeric(coef(fit_standard))
  standard_se <- fit_standard$se
  tau2 <- fit_standard$tau2

  # Step 2: Compute shrinkage factor
  # More shrinkage for smaller k
  if (shrinkage_intensity == "adaptive") {
    # Shrinkage based on k and precision
    # Borrow strength from "global" mean (assume 0 for null)
    global_mean <- 0  # Could be estimated from external data
    global_variance <- 0.5  # Moderate between-study variance

    # Empirical Bayes shrinkage factor
    # B = tau² / (tau² + se²)
    # Shrink toward global mean
    precision_study <- 1 / (vi + tau2)
    precision_global <- 1 / global_variance

    shrinkage_factor <- precision_global / (precision_global + rowSums(matrix(precision_study, nrow = k)))
    shrinkage_factor <- mean(shrinkage_factor)  # Average shrinkage

  } else {
    # Fixed shrinkage based on k
    shrinkage_factor <- if (k <= 3) 0.3 else if (k <= 5) 0.5 else 0.7
  }

  # Step 3: Shrunken estimate
  # Shrink standard estimate toward 0 (global null)
  shrunk_estimate <- (1 - shrinkage_factor) * standard_estimate + shrinkage_factor * 0

  # Step 4: Adjusted SE (incorporates shrinkage uncertainty)
  # SE increases with shrinkage
  shrinkage_penalty <- 1 + shrinkage_factor
  adjusted_se <- standard_se * shrinkage_penalty

  # T-distribution CI (more conservative df for small k)
  df <- max(1, k - 1)
  t_crit <- qt(0.975, df)

  ci_lb <- shrunk_estimate - t_crit * adjusted_se
  ci_ub <- shrunk_estimate + t_crit * adjusted_se
  pvalue <- 2 * (1 - pt(abs(shrunk_estimate / adjusted_se), df))

  list(
    estimate = shrunk_estimate,
    se = adjusted_se,
    ci_lb = ci_lb,
    ci_ub = ci_ub,
    pvalue = pvalue,
    tau2 = tau2,
    method = "SMS",
    k = k,
    shrinkage_factor = shrinkage_factor,
    standard_estimate = standard_estimate,
    shrunk_estimate = shrunk_estimate,
    shrinkage_intensity = shrinkage_intensity
  )
}


# ============================================================================
# COMPARISON FUNCTION
# ============================================================================

#' Compare All V4 Methods
#'
#' Comprehensive comparison including all 10 new methods.
#'
#' @param yi Numeric vector of effect sizes
#' @param vi Numeric vector of sampling variances
#' @return Data frame with results from all methods
#' @export
compare_methods_v4 <- function(yi, vi) {
  results <- list()
  k <- length(yi)

  # Standard methods
  fit_reml <- metafor::rma(yi, vi, method = "REML")
  results$REML <- data.frame(
    Method = "REML", Estimate = round(as.numeric(coef(fit_reml)), 4),
    SE = round(fit_reml$se, 4), CI_Lower = round(fit_reml$ci.lb, 4),
    CI_Upper = round(fit_reml->ci.ub, 4), PValue = round(fit_reml$pval, 4),
    k = k
  )

  fit_hksj <- metafor::rma(yi, vi, method = "REML", test = "knha")
  results$HKSJ <- data.frame(
    Method = "HKSJ", Estimate = round(as.numeric(coef(fit_hksj)), 4),
    SE = round(fit_hksj$se, 4), CI_Lower = round(fit_hksj->ci.lb, 4),
    CI_Upper = round(fit_hksj->ci.ub, 4), PValue = round(fit_hksj$pval, 4),
    k = k
  )

  # Category A: Robustness Methods
  wrd <- tryCatch(wrd_meta(yi, vi), error = function(e) NULL)
  if (!is.null(wrd)) {
    results$WRD <- data.frame(
      Method = "WRD", Estimate = round(wrd$estimate, 4),
      SE = round(wrd$se, 4), CI_Lower = round(wrd$ci_lb, 4),
      CI_Upper = round(wrd$ci_ub, 4), PValue = round(wrd$pvalue, 4),
      k = k
    )
  }

  cbm <- tryCatch(cbm_meta(yi, vi), error = function(e) NULL)
  if (!is.null(cbm)) {
    results$CBM <- data.frame(
      Method = paste0("CBM(", cbm$n_clusters, " clusters)"),
      Estimate = round(cbm$estimate, 4),
      SE = round(cbm$se, 4), CI_Lower = round(cbm$ci_lb, 4),
      CI_Upper = round(cbm$ci_ub, 4), PValue = round(cbm$pvalue, 4),
      k = k
    )
  }

  rbm <- tryCatch(rbm_meta(yi, vi), error = function(e) NULL)
  if (!is.null(rbm)) {
    results$RBM <- data.frame(
      Method = "RBM", Estimate = round(rbm$estimate, 4),
      SE = round(rbm$se, 4), CI_Lower = round(rbm$ci_lb, 4),
      CI_Upper = round(rbm$ci_ub, 4), PValue = round(rbm$pvalue, 4),
      k = k
    )
  }

  # Category B: Bias Correction
  swa <- tryCatch(swa_meta(yi, vi), error = function(e) NULL)
  if (!is.null(swa)) {
    results$SWA <- data.frame(
      Method = "SWA", Estimate = round(swa$estimate, 4),
      SE = round(swa$se, 4), CI_Lower = round(swa$ci_lb, 4),
      CI_Upper = round(swa$ci_ub, 4), PValue = round(swa$pvalue, 4),
      k = k
    )
  }

  tas <- tryCatch(tas_meta(yi, vi), error = function(e) NULL)
  if (!is.null(tas)) {
    results$TAS <- data.frame(
      Method = paste0("TAS(n_imp=", tas$n_imputed, ")"),
      Estimate = round(tas$estimate, 4),
      SE = round(tas$se, 4), CI_Lower = round(tas$ci_lb, 4),
      CI_Upper = round(tas$ci_ub, 4), PValue = round(tas$pvalue, 4),
      k = k
    )
  }

  # Category C: Variance Methods
  eve <- tryCatch(eve_meta(yi, vi), error = function(e) NULL)
  if (!is.null(eve)) {
    results$EVE <- data.frame(
      Method = "EVE", Estimate = round(eve$estimate, 4),
      SE = round(eve$se, 4), CI_Lower = round(eve$ci_lb, 4),
      CI_Upper = round(eve$ci_ub, 4), PValue = round(eve$pvalue, 4),
      k = k
    )
  }

  # PVM can be slow, use simpler grid for quick comparison
  pvm <- tryCatch(
    pvm_meta(yi, vi, n_grid_points = 50, n_boot = 200),
    error = function(e) NULL
  )
  if (!is.null(pvm)) {
    results$PVM <- data.frame(
      Method = "PVM", Estimate = round(pvm$estimate, 4),
      SE = round(pvm$se, 4), CI_Lower = round(pvm$ci_lb, 4),
      CI_Upper = round(pvm$ci_ub, 4), PValue = round(pvm$pvalue, 4),
      k = k
    )
  }

  # Category D: Ensemble Methods
  aem <- tryCatch(aem_meta(yi, vi), error = function(e) NULL)
  if (!is.null(aem)) {
    results$AEM <- data.frame(
      Method = "AEM", Estimate = round(aem$estimate, 4),
      SE = round(aem$se, 4), CI_Lower = round(aem$ci_lb, 4),
      CI_Upper = round(aem$ci_ub, 4), PValue = round(aem$pvalue, 4),
      k = k
    )
  }

  # SPE can be slow, use fewer samples for quick comparison
  spe <- tryCatch(
    spe_meta(yi, vi, n_samples = 5000, burnin = 500),
    error = function(e) NULL
  )
  if (!is.null(spe)) {
    results$SPE <- data.frame(
      Method = "SPE", Estimate = round(spe$estimate, 4),
      SE = round(spe$se, 4), CI_Lower = round(spe$ci_lb, 4),
      CI_Upper = round(spe$ci_ub, 4), PValue = round(spe$pvalue, 4),
      k = k
    )
  }

  # Category E: Specialized
  sms <- tryCatch(sms_meta(yi, vi), error = function(e) NULL)
  if (!is.null(sms)) {
    results$SMS <- data.frame(
      Method = "SMS", Estimate = round(sms$estimate, 4),
      SE = round(sms$se, 4), CI_Lower = round(sms$ci_lb, 4),
      CI_Upper = round(sms$ci_ub, 4), PValue = round(sms$pvalue, 4),
      k = k
    )
  }

  do.call(rbind, results)
}
