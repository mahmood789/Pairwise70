# ============================================================================
# Comprehensive Testing Framework for Advanced Pooling Methods V4
# ============================================================================
# Testing 10 novel methods across 25 scenarios with 8 performance metrics
# ============================================================================

# Required packages
required_packages <- c("metafor", "parallel", "data.table")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing package:", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# Optional packages (for specific methods)
optional_packages <- c("mclust", "cluster", "boot")
for (pkg in optional_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Optional package not installed:", pkg, "- some methods may use fallback"))
  }
}

library(metafor)
library(parallel)
library(data.table)

# Source the new methods
source("../../R/advanced_pooling_v4.R", chdir = TRUE)
# For MWM reference - source existing methods if available
tryCatch(source("../../R/advanced_pooling.R", chdir = TRUE), error = function(e) {
  message("Note: V3 methods not available - running V4 only")
})

# ============================================================================
# SIMULATION SCENARIOS (25 total)
# ============================================================================

#' Get all simulation scenarios
#'
#' @return List of scenario definitions
get_all_scenarios <- function() {
  list(
    # Baseline Scenarios (5)
    list(id = "B1", name = "Null effect, small k", k = 5, tau2 = 0.05, true_effect = 0.0,
         type = "baseline", outlier = NULL, pub_bias = NULL),
    list(id = "B2", name = "Standard, medium k", k = 10, tau2 = 0.05, true_effect = 0.3,
         type = "baseline", outlier = NULL, pub_bias = NULL),
    list(id = "B3", name = "Standard, large k", k = 20, tau2 = 0.05, true_effect = 0.5,
         type = "baseline", outlier = NULL, pub_bias = NULL),
    list(id = "B4", name = "Large k", k = 50, tau2 = 0.05, true_effect = 0.3,
         type = "baseline", outlier = NULL, pub_bias = NULL),
    list(id = "B5", name = "Very large k", k = 100, tau2 = 0.05, true_effect = 0.3,
         type = "baseline", outlier = NULL, pub_bias = NULL),

    # Heterogeneity Scenarios (5)
    list(id = "H1", name = "No heterogeneity", k = 10, tau2 = 0.0, true_effect = 0.3,
         type = "heterogeneity", outlier = NULL, pub_bias = NULL),
    list(id = "H2", name = "Moderate heterogeneity", k = 10, tau2 = 0.10, true_effect = 0.3,
         type = "heterogeneity", outlier = NULL, pub_bias = NULL),
    list(id = "H3", name = "High heterogeneity", k = 10, tau2 = 0.20, true_effect = 0.3,
         type = "heterogeneity", outlier = NULL, pub_bias = NULL),
    list(id = "H4", name = "Very high heterogeneity", k = 10, tau2 = 0.50, true_effect = 0.3,
         type = "heterogeneity", outlier = NULL, pub_bias = NULL),
    list(id = "H5", name = "Extreme heterogeneity", k = 10, tau2 = 1.00, true_effect = 0.3,
         type = "heterogeneity", outlier = NULL, pub_bias = NULL),

    # Outlier Scenarios (5)
    list(id = "O1", name = "Single mild outlier", k = 10, tau2 = 0.05, true_effect = 0.3,
         type = "outlier", outlier = list(n = 1, shift = 3 * sqrt(0.05)), pub_bias = NULL),
    list(id = "O2", name = "Single extreme outlier", k = 10, tau2 = 0.05, true_effect = 0.3,
         type = "outlier", outlier = list(n = 1, shift = 5 * sqrt(0.05)), pub_bias = NULL),
    list(id = "O3", name = "Two opposing outliers", k = 10, tau2 = 0.05, true_effect = 0.3,
         type = "outlier", outlier = list(n = 2, shift = c(3, -3) * sqrt(0.05)), pub_bias = NULL),
    list(id = "O4", name = "Multiple outliers", k = 15, tau2 = 0.05, true_effect = 0.3,
         type = "outlier", outlier = list(n = 3, shift = rep(3, 3) * sqrt(0.05)), pub_bias = NULL),
    list(id = "O5", name = "Clustered outliers", k = 10, tau2 = 0.05, true_effect = 0.3,
         type = "outlier", outlier = list(n = 1, shift = 4 * sqrt(0.05)), pub_bias = NULL),

    # Publication Bias Scenarios (5)
    list(id = "PB1", name = "Mild pub bias (p<0.10)", k = 20, tau2 = 0.05, true_effect = 0.3,
         type = "pub_bias", outlier = NULL, pub_bias = list(type = "step", cutoff = 0.10)),
    list(id = "PB2", name = "Moderate pub bias (p<0.05)", k = 20, tau2 = 0.05, true_effect = 0.3,
         type = "pub_bias", outlier = NULL, pub_bias = list(type = "step", cutoff = 0.05)),
    list(id = "PB3", name = "Strong pub bias (p<0.01)", k = 20, tau2 = 0.05, true_effect = 0.3,
         type = "pub_bias", outlier = NULL, pub_bias = list(type = "step", cutoff = 0.01)),
    list(id = "PB4", name = "Classic step function", k = 20, tau2 = 0.05, true_effect = 0.3,
         type = "pub_bias", outlier = NULL, pub_bias = list(type = "continuous")),
    list(id = "PB5", name = "Continuous weight function", k = 20, tau2 = 0.05, true_effect = 0.3,
         type = "pub_bias", outlier = NULL, pub_bias = list(type = "weighted", weight = "severity")),

    # Small Study Scenarios (3)
    list(id = "S1", name = "Very small studies", k = 10, tau2 = 0.05, true_effect = 0.3,
         type = "small_study", mean_n = 20, outlier = NULL, pub_bias = NULL),
    list(id = "S2", name = "Small studies", k = 10, tau2 = 0.05, true_effect = 0.3,
         type = "small_study", mean_n = 50, outlier = NULL, pub_bias = NULL),
    list(id = "S3", name = "Mixed study sizes", k = 10, tau2 = 0.05, true_effect = 0.3,
         type = "small_study", mean_n = c(20, 500), outlier = NULL, pub_bias = NULL),

    # Distribution Scenarios (2)
    list(id = "D1", name = "Heavy-tailed effects", k = 10, tau2 = 0.05, true_effect = 0.3,
         type = "distribution", effect_dist = "t", effect_df = 3, outlier = NULL, pub_bias = NULL),
    list(id = "D2", name = "Skewed effects", k = 10, tau2 = 0.05, true_effect = 0.3,
         type = "distribution", effect_dist = "lognormal", outlier = NULL, pub_bias = NULL)
  )
}

#' Generate simulated data for a scenario
#'
#' @param scenario Scenario definition from get_all_scenarios()
#' @param seed Random seed for reproducibility
#' @return List with yi (effect sizes) and vi (variances)
generate_scenario_data <- function(scenario, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  k <- scenario$k
  tau2 <- scenario$tau2
  true_effect <- scenario$true_effect

  # Generate true effect sizes
  if (!is.null(scenario$effect_dist)) {
    if (scenario$effect_dist == "t") {
      df <- scenario$effect_df
      true_effects <- true_effect + rt(k, df) * sqrt(tau2)
    } else if (scenario$effect_dist == "lognormal") {
      # Skewed distribution
      true_effects <- true_effect + rlnorm(k, 0, sqrt(tau2)) - exp(tau2/2)
    }
  } else {
    # Normal distribution
    true_effects <- rnorm(k, true_effect, sqrt(tau2))
  }

  # Generate sample sizes and variances
  if (!is.null(scenario$mean_n)) {
    if (length(scenario$mean_n) == 2) {
      # Bimodal distribution
      n_per_arm <- ifelse(runif(k) < 0.5, scenario$mean_n[1], scenario$mean_n[2])
    } else {
      n_per_arm <- rpois(k, scenario$mean_n) + 20  # Ensure minimum
    }
    # Approximate variance for log OR or SMD
    vi <- 4 / n_per_arm + tau2 * 0.1
  } else {
    # Standard variance
    vi <- rexp(k, rate = 10) + 0.01
    vi <- pmax(vi, 0.01)  # Minimum variance
  }

  # Add outliers if specified
  if (!is.null(scenario$outlier)) {
    n_outliers <- scenario$outlier$n
    shift <- scenario$outlier$shift
    outlier_idx <- sample(1:k, min(n_outliers, k))
    true_effects[outlier_idx] <- true_effects[outlier_idx] + shift[1:min(n_outliers, length(shift))]
  }

  # Observed effect sizes
  yi <- rnorm(k, true_effects, sqrt(vi))

  list(yi = yi, vi = vi, true_effects = true_effects, true_effect = true_effect)
}

#' Apply publication bias to data
#'
#' @param data Data list from generate_scenario_data()
#' @param pub_bias Publication bias specification
#' @return Filtered/inflated data
apply_publication_bias <- function(data, pub_bias) {
  if (is.null(pub_bias)) return(data)

  yi <- data$yi
  vi <- data$vi

  # Calculate z-values
  z <- yi / sqrt(vi)
  p <- 2 * (1 - pnorm(abs(z)))

  if (pub_bias$type == "step") {
    # Exclude studies with p > cutoff
    keep <- p <= pub_bias$cutoff
    # Keep at least 50% of studies
    if (sum(keep) < length(yi) / 2) {
      keep <- rank(p) <= floor(length(yi) / 2)
    }
    data$yi <- yi[keep]
    data$vi <- vi[keep]
  } else if (pub_bias$type == "continuous") {
    # Weight function: probability of publication
    pub_prob <- 1 / (1 + exp(-2 * (1 - p/0.05)))  # Sigmoid at p=0.05
    keep <- runif(length(yi)) < pub_prob
    data$yi <- yi[keep]
    data$vi <- vi[keep]
  }

  data
}

# ============================================================================
# METHODS TO TEST
# ============================================================================

#' Get all methods to test
#'
#' @return Named list of method functions
get_all_methods <- function() {
  list(
    # Standard methods
    REML = function(yi, vi) {
      fit <- rma(yi, vi, method = "REML")
      list(estimate = coef(fit), se = fit$se, ci_lb = fit$ci.lb,
           ci_ub = fit$ci.ub, pvalue = fit$pval, tau2 = fit$tau2)
    },

    HKSJ = function(yi, vi) {
      fit <- rma(yi, vi, method = "REML", test = "knha")
      list(estimate = coef(fit), se = fit$se, ci_lb = fit$ci.lb,
           ci_ub = fit$ci.ub, pvalue = fit$pval, tau2 = fit$tau2)
    },

    # V3 methods
    MWM = function(yi, vi) {
      mafi_weighted_ma(yi, vi)
    },

    # V4 methods
    WRD = function(yi, vi) {
      wrd_meta(yi, vi, bootstrap = FALSE)  # Fast for simulation
    },

    CBM = function(yi, vi) {
      cbm_meta(yi, vi)
    },

    RBM = function(yi, vi) {
      rbm_meta(yi, vi)
    },

    SWA = function(yi, vi) {
      swa_meta(yi, vi, n_boot = 200)  # Fast for simulation
    },

    TAS = function(yi, vi) {
      tas_meta(yi, vi)
    },

    EVE = function(yi, vi) {
      eve_meta(yi, vi)
    },

    PVM = function(yi, vi) {
      pvm_meta(yi, vi, n_grid_points = 30, n_boot = 100)  # Fast
    },

    AEM = function(yi, vi) {
      aem_meta(yi, vi)
    },

    SPE = function(yi, vi) {
      spe_meta(yi, vi, n_samples = 2000, burnin = 200)  # Fast
    },

    SMS = function(yi, vi) {
      sms_meta(yi, vi)
    }
  )
}

# ============================================================================
# SIMULATION ENGINE
# ============================================================================

#' Run single simulation iteration
#'
#' @param scenario Scenario definition
#' @param methods Named list of method functions
#' @param seed Random seed
#' @return Data frame with results
run_single_simulation <- function(scenario, methods, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Generate data
  data <- generate_scenario_data(scenario)

  # Apply publication bias if specified
  if (!is.null(scenario$pub_bias)) {
    data <- apply_publication_bias(data, scenario$pub_bias)
  }

  yi <- data$yi
  vi <- data$vi
  true_effect <- data$true_effect

  # Check minimum k for each method
  if (length(yi) < 3) {
    return(NULL)  # Skip if too few studies
  }

  results <- list()

  for (method_name in names(methods)) {
    result <- tryCatch({
      res <- methods[[method_name]](yi, vi)
      # Ensure estimate is a single numeric value
      res$estimate <- as.numeric(res$estimate)[1]
      res$se <- as.numeric(res$se)[1]
      res$ci_lb <- as.numeric(res$ci_lb)[1]
      res$ci_ub <- as.numeric(res$ci_ub)[1]
      res$pvalue <- as.numeric(res$pvalue)[1]
      res$tau2 <- as.numeric(res$tau2)[1]
      res
    }, error = function(e) {
      list(estimate = NA_real_, se = NA_real_, ci_lb = NA_real_, ci_ub = NA_real_,
           pvalue = NA_real_, tau2 = NA_real_, error = as.character(e))
    })

    results[[method_name]] <- data.frame(
      scenario = scenario$id,
      method = method_name,
      estimate = result$estimate,
      se = result$se,
      ci_lb = result$ci_lb,
      ci_ub = result$ci_ub,
      pvalue = result$pvalue,
      tau2 = result$tau2,
      true_effect = true_effect[1],
      covered = if (!is.na(result$ci_lb) && !is.na(result$ci_ub)) {
        true_effect[1] >= result$ci_lb & true_effect[1] <= result$ci_ub
      } else NA,
      ci_width = if (!is.na(result$ci_lb) && !is.na(result$ci_ub)) {
        result$ci_ub - result$ci_lb
      } else NA_real_,
      error = if (!is.null(result$error)) result$error else "",
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results)
}

#' Run parallel simulation study
#'
#' @param n_sim Number of iterations per scenario
#' @param scenarios List of scenarios (default: all)
#' @param n_cores Number of cores for parallel processing
#' @param seed Random seed
#' @return Data frame with all results
run_comprehensive_simulation <- function(
    n_sim = 100,
    scenarios = get_all_scenarios(),
    n_cores = parallel::detectCores() - 1,
    seed = 20260115
) {
  set.seed(seed)

  # Create parameter grid
  params <- expand.grid(
    scenario_idx = 1:length(scenarios),
    sim = 1:n_sim,
    stringsAsFactors = FALSE
  )

  # Create seeds for each iteration
  iteration_seeds <- sample(1:1e6, nrow(params))

  # Define simulation function for parallel processing
  sim_function <- function(i, scenarios, methods, seeds) {
    scenario <- scenarios[[params$scenario_idx[i]]]
    run_single_simulation(scenario, methods, seeds[i])
  }

  # Run in parallel
  message(paste("Running", nrow(params), "simulations across", n_cores, "cores..."))

  # Windows doesn't support mclapply, use parLapply or regular lapply
  if (.Platform$OS.type == "windows") {
    if (n_cores > 1) {
      cl <- makeCluster(n_cores)
      clusterExport(cl, c("params", "scenarios", "get_all_methods", "iteration_seeds",
                          "run_single_simulation"))
      results_list <- parLapply(cl, 1:nrow(params), function(i) {
        scenario <- scenarios[[params$scenario_idx[i]]]
        run_single_simulation(scenario, get_all_methods(), iteration_seeds[i])
      })
      stopCluster(cl)
    } else {
      results_list <- lapply(1:nrow(params), function(i) {
        scenario <- scenarios[[params$scenario_idx[i]]]
        run_single_simulation(scenario, get_all_methods(), iteration_seeds[i])
      })
    }
  } else {
    results_list <- mclapply(
      1:nrow(params),
      function(i) sim_function(i, scenarios, get_all_methods(), iteration_seeds),
      mc.cores = n_cores
    )
  }

  # Combine results
  results_dt <- rbindlist(results_list, fill = TRUE)
  results_dt <- as.data.frame(results_dt)

  message(paste("Completed:", nrow(results_dt), "results from", nrow(params), "simulations"))

  results_dt
}

# ============================================================================
# PERFORMANCE METRICS
# ============================================================================

#' Compute performance metrics for simulation results
#'
#' @param results Data frame from run_comprehensive_simulation()
#' @return Data frame with aggregated metrics
compute_performance_metrics <- function(results) {
  # Convert to data.table if not already
  if (!is.data.table(results)) {
    results <- as.data.table(results)
  }

  metrics <- results[, .(
    bias = abs(mean(estimate, na.rm = TRUE) - mean(true_effect[1], na.rm = TRUE)),
    rmse = sqrt(mean((estimate - true_effect)^2, na.rm = TRUE)),
    coverage = mean(covered == TRUE, na.rm = TRUE),
    ci_width = mean(ci_width, na.rm = TRUE),
    coverage_x_width = mean(covered == TRUE, na.rm = TRUE) * mean(ci_width, na.rm = TRUE),
    type_i_error = if (mean(unique(true_effect)) == 0) mean(pvalue < 0.05, na.rm = TRUE) else NA_real_,
    power = if (mean(unique(true_effect)) != 0) mean(pvalue < 0.05, na.rm = TRUE) else NA_real_,
    convergence_rate = mean(error == "", na.rm = TRUE),
    n_sim = .N
  ), by = .(scenario, method)]

  # Add scenario type and name
  scenarios_info <- setNames(
    lapply(get_all_scenarios(), function(s) data.frame(
      type = s$type,
      name = s$name,
      k = s$k,
      tau2 = s$tau2
    )),
    sapply(get_all_scenarios(), function(s) s$id)
  )

  for (sc in unique(metrics$scenario)) {
    metrics$type[metrics$scenario == sc] <- scenarios_info[[sc]]$type
    metrics$name[metrics$scenario == sc] <- scenarios_info[[sc]]$name
    metrics$k[metrics$scenario == sc] <- scenarios_info[[sc]]$k
    metrics$tau2[metrics$scenario == sc] <- scenarios_info[[sc]]$tau2
  }

  metrics
}

# ============================================================================
# QUICK VALIDATION (100 iterations)
# ============================================================================

#' Run quick validation of all methods
#'
#' @param n_sim Number of iterations (default 100)
#' @param n_cores Number of cores (default 4)
#' @return List with results and metrics
run_quick_validation <- function(n_sim = 100, n_cores = 4) {
  message("=== QUICK VALIDATION ===")
  message(paste("Iterations:", n_sim, "x", length(get_all_scenarios()), "scenarios"))

  results <- run_comprehensive_simulation(
    n_sim = n_sim,
    scenarios = get_all_scenarios(),
    n_cores = n_cores
  )

  metrics <- compute_performance_metrics(results)

  # Summary by method
  method_summary <- metrics[, .(
    mean_bias = mean(bias, na.rm = TRUE),
    mean_rmse = mean(rmse, na.rm = TRUE),
    mean_coverage = mean(coverage, na.rm = TRUE),
    mean_ci_width = mean(ci_width, na.rm = TRUE),
    mean_convergence = mean(convergence_rate, na.rm = TRUE)
  ), by = method]

  method_summary <- method_summary[order(mean_rmse)]

  print(method_summary)

  list(results = results, metrics = metrics, method_summary = method_summary)
}

# ============================================================================
# REAL DATA VALIDATION
# ============================================================================

#' Run methods on all Pairwise70 datasets
#'
#' @param max_datasets Maximum number of datasets to analyze (default 50)
#' @return Data frame with results
run_real_data_validation <- function(max_datasets = 50) {
  message("=== REAL DATA VALIDATION ===")

  # Get available datasets
  available_data <- data(package = "Pairwise70")$results[, 3]

  # Sample datasets if needed
  if (length(available_data) > max_datasets) {
    datasets_to_use <- sample(available_data, max_datasets)
  } else {
    datasets_to_use <- available_data
  }

  message(paste("Analyzing", length(datasets_to_use), "datasets..."))

  results_list <- list()

  for (i in seq_along(datasets_to_use)) {
    ds_name <- datasets_to_use[i]

    tryCatch({
      # Load dataset
      data(list = ds_name, envir = environment())
      dat <- get(ds_name, envir = environment())

      # Check if it has required columns
      if (!all(c("Experimental.cases", "Experimental.N", "Control.cases", "Control.N") %in% names(dat))) {
        next
      }

      # Calculate effect sizes (log odds ratio)
      dat <- escalc(measure = "OR",
                    ai = Experimental.cases, n1i = Experimental.N,
                    ci = Control.cases, n2i = Control.N,
                    data = dat)

      yi <- dat$yi
      vi <- dat$vi

      if (length(yi) < 3) next

      # Run all methods
      comparison <- compare_methods_v4(yi, vi)

      results_list[[ds_name]] <- cbind(dataset = ds_name, k = length(yi), comparison)

    }, error = function(e) {
      message(paste("Error in", ds_name, ":", e$message))
    })
  }

  results_dt <- rbindlist(results_list, fill = TRUE)
  results_dt <- as.data.frame(results_dt)

  message(paste("Completed analysis of", length(unique(results_dt$dataset)), "datasets"))

  results_dt
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

if (FALSE) {  # Set to TRUE to run
  # Quick validation
  quick_results <- run_quick_validation(n_sim = 100, n_cores = 4)

  # Save results
  saveRDS(quick_results$results, file = "../results/v4_quick_validation.rds")
  saveRDS(quick_results$metrics, file = "../results/v4_quick_validation_metrics.rds")

  # Real data validation
  real_data_results <- run_real_data_validation(max_datasets = 50)
  saveRDS(real_data_results, file = "../results/v4_real_data.rds")
}
