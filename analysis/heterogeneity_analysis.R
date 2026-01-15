#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(metafor)
  library(MASS)
})

set.seed(123)

out_dir <- file.path("analysis", "heterogeneity_output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE)

rds_dir <- file.path("analysis", "output", "cleaned_rds")
files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)

if (length(files) == 0) {
  stop("No cleaned RDS files found in analysis/output/cleaned_rds.")
}

map_bias_score <- function(x) {
  x <- trimws(tolower(as.character(x)))
  ifelse(x == "low risk", 0,
         ifelse(x == "some concerns", 0.5,
                ifelse(x == "high risk", 1, NA_real_)))
}

build_domain_text <- function(df) {
  fields <- c("Analysis.group", "Analysis.name", "Subgroup", "Applicability", "review.url", "review.doi")
  parts <- c()
  for (field in fields) {
    if (field %in% names(df)) {
      vals <- unique(na.omit(as.character(df[[field]])))
      if (length(vals) > 0) {
        parts <- c(parts, vals[seq_len(min(length(vals), 5))])
      }
    }
  }
  tolower(paste(parts, collapse = " | "))
}

assign_domain <- function(text) {
  if (grepl("pregnan|postpartum|perinatal|neonat|obstet|caesarean|labor|labour|delivery|breastfeed", text)) return("pregnancy")
  if (grepl("cardiac|heart|myocard|angina|stroke|hypertens|vascular|thrombo|atrial|statin|lipid", text)) return("cardio")
  if (grepl("cancer|carcinoma|tumou?r|neoplasm|malignan|chemotherap|radiotherap", text)) return("oncology")
  if (grepl("depress|anxiet|schizophren|bipolar|psych|ptsd|mental", text)) return("mental_health")
  if (grepl("parkinson|alzheimer|epilep|migraine|neuro|dementia|multiple sclerosis", text)) return("neuro")
  if (grepl("infect|antibiot|antiviral|vaccine|hiv|malaria|tuberc|pneumonia|sepsis", text)) return("infectious")
  if (grepl("asthma|copd|pulmon|bronch|respirat", text)) return("respiratory")
  if (grepl("diabet|glyc|insulin|metabolic|obes", text)) return("metabolic")
  if (grepl("arthrit|osteo|fracture|bone|joint|osteopor", text)) return("musculoskeletal")
  if (grepl("renal|kidney|dialysis|nephro", text)) return("renal")
  if (grepl("liver|hepat|biliary|gastro|colon|bowel|pancrea", text)) return("gastro_hepatic")
  if (grepl("skin|psoriasis|eczema|dermat", text)) return("dermatology")
  "other"
}

assign_intervention_class <- function(text) {
  if (grepl("surg|surgery|operative|lapar|resection|transplant", text)) return("surgical")
  if (grepl("device|implant|stent|catheter|prosthe|ventilat|pacemaker", text)) return("device")
  if (grepl("diagnos|screen|imaging|ultrasound|ct scan|mri|test accuracy", text)) return("diagnostic")
  if (grepl("rehab|exercise|education|counsel|psychotherap|behavior|behaviour|lifestyle", text)) return("behavioral")
  if (grepl("drug|pharm|medication|antibiot|antiviral|vaccine|dose|tablet|pill|injection", text)) return("pharmacologic")
  "other"
}

skewness <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  m <- mean(x)
  s <- sd(x)
  if (s == 0) return(NA_real_)
  mean((x - m)^3) / (s^3)
}

kurtosis_excess <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4) return(NA_real_)
  m <- mean(x)
  s <- sd(x)
  if (s == 0) return(NA_real_)
  mean((x - m)^4) / (s^4) - 3
}

detect_outcome_type <- function(df) {
  has_bin_cols <- all(c("Experimental.cases", "Experimental.N", "Control.cases", "Control.N") %in% names(df))
  has_cont_cols <- all(c("Experimental.mean", "Experimental.SD", "Experimental.N",
                         "Control.mean", "Control.SD", "Control.N") %in% names(df))
  has_giv_cols <- all(c("GIV.Mean", "GIV.SE") %in% names(df))
  has_mean_var <- all(c("Mean", "Variance") %in% names(df))

  bin_has_data <- has_bin_cols &&
    (any(!is.na(df[["Experimental.cases"]])) || any(!is.na(df[["Control.cases"]])))
  cont_has_data <- has_cont_cols &&
    (any(!is.na(df[["Experimental.mean"]])) || any(!is.na(df[["Control.mean"]])))
  giv_has_data <- has_giv_cols && any(!is.na(df[["GIV.Mean"]])) && any(!is.na(df[["GIV.SE"]]))
  meanvar_has_data <- has_mean_var && any(!is.na(df[["Mean"]]))

  if (bin_has_data) return("binary")
  if (cont_has_data) return("continuous")
  if (giv_has_data) return("giv")
  if (meanvar_has_data) return("mean_var")
  "none"
}

extract_effects <- function(df) {
  outcome_type <- detect_outcome_type(df)
  effect_measure <- NA_character_
  es <- NULL

  if (outcome_type == "binary") {
    use <- complete.cases(df[, c("Experimental.cases", "Experimental.N", "Control.cases", "Control.N")])
    d <- df[use, , drop = FALSE]
    if (nrow(d) >= 2) {
      valid <- d[["Experimental.cases"]] >= 0 &
        d[["Control.cases"]] >= 0 &
        d[["Experimental.N"]] > 0 &
        d[["Control.N"]] > 0 &
        d[["Experimental.cases"]] <= d[["Experimental.N"]] &
        d[["Control.cases"]] <= d[["Control.N"]]
      d <- d[valid, , drop = FALSE]
    }
    if (nrow(d) >= 2) {
      es_tmp <- escalc(measure = "OR",
                       ai = d[["Experimental.cases"]],
                       n1i = d[["Experimental.N"]],
                       ci = d[["Control.cases"]],
                       n2i = d[["Control.N"]],
                       add = 0.5, to = "all")
      es <- cbind(d, es_tmp[, c("yi", "vi")])
      effect_measure <- "OR"
    }
  } else if (outcome_type == "continuous") {
    sd1 <- if ("Experimental.SD.cleaned" %in% names(df)) df[["Experimental.SD.cleaned"]] else df[["Experimental.SD"]]
    sd2 <- if ("Control.SD.cleaned" %in% names(df)) df[["Control.SD.cleaned"]] else df[["Control.SD"]]
    use <- complete.cases(df[, c("Experimental.mean", "Experimental.N", "Control.mean", "Control.N")]) &
      !is.na(sd1) & !is.na(sd2)
    d <- df[use, , drop = FALSE]
    sd1 <- sd1[use]
    sd2 <- sd2[use]
    if (nrow(d) >= 2) {
      es_tmp <- escalc(measure = "SMD",
                       m1i = d[["Experimental.mean"]],
                       sd1i = sd1,
                       n1i = d[["Experimental.N"]],
                       m2i = d[["Control.mean"]],
                       sd2i = sd2,
                       n2i = d[["Control.N"]])
      es <- cbind(d, es_tmp[, c("yi", "vi")])
      effect_measure <- "SMD"
    }
  } else if (outcome_type %in% c("giv", "mean_var")) {
    yi <- if ("GIV.Mean" %in% names(df)) df[["GIV.Mean"]] else rep(NA_real_, nrow(df))
    if (all(is.na(yi)) && "Mean" %in% names(df)) yi <- df[["Mean"]]

    se <- if ("GIV.SE" %in% names(df)) df[["GIV.SE"]] else rep(NA_real_, nrow(df))
    if (all(is.na(se)) && "Variance" %in% names(df)) se <- sqrt(df[["Variance"]])

    if (all(is.na(se)) && all(c("CI.start", "CI.end") %in% names(df))) {
      se <- (df[["CI.end"]] - df[["CI.start"]]) / (2 * 1.96)
    }

    use <- !is.na(yi) & !is.na(se) & is.finite(se) & se > 0
    d <- df[use, , drop = FALSE]
    if (nrow(d) >= 2) {
      es <- d
      es[["yi"]] <- yi[use]
      es[["vi"]] <- se[use]^2
      effect_measure <- "GIV"
    }
  }

  if (is.null(es)) {
    return(list(data = NULL, outcome_type = outcome_type, effect_measure = effect_measure))
  }

  list(data = es, outcome_type = outcome_type, effect_measure = effect_measure)
}

fit_t_distribution <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 10) return(NULL)
  nll <- function(par) {
    df <- exp(par[1]) + 1
    scale <- exp(par[2])
    -sum(dt(x / scale, df = df, log = TRUE) - log(scale))
  }
  fit <- optim(c(log(10), log(sd(x))), nll, method = "BFGS")
  if (fit$convergence != 0) return(NULL)
  df <- exp(fit$par[1]) + 1
  scale <- exp(fit$par[2])
  ll <- -fit$value
  list(df = df, scale = scale, loglik = ll, k = 2)
}

fit_skew_normal <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 10) return(NULL)
  nll <- function(par) {
    xi <- par[1]
    omega <- exp(par[2])
    alpha <- par[3]
    z <- (x - xi) / omega
    dens <- 2 / omega * dnorm(z) * pnorm(alpha * z)
    if (any(dens <= 0)) return(Inf)
    -sum(log(dens))
  }
  fit <- optim(c(mean(x), log(sd(x)), 0), nll, method = "BFGS")
  if (fit$convergence != 0) return(NULL)
  xi <- fit$par[1]
  omega <- exp(fit$par[2])
  alpha <- fit$par[3]
  ll <- -fit$value
  list(xi = xi, omega = omega, alpha = alpha, loglik = ll, k = 3)
}

fit_half_normal <- function(x) {
  x <- x[is.finite(x) & x >= 0]
  if (length(x) < 10) return(NULL)
  nll <- function(par) {
    sigma <- exp(par[1])
    -sum(log(sqrt(2 / pi)) - log(sigma) - (x^2) / (2 * sigma^2))
  }
  fit <- optim(log(sd(x)), nll, method = "BFGS")
  if (fit$convergence != 0) return(NULL)
  sigma <- exp(fit$par[1])
  ll <- -fit$value
  list(sigma = sigma, loglik = ll, k = 1)
}

fit_half_t <- function(x) {
  x <- x[is.finite(x) & x >= 0]
  if (length(x) < 10) return(NULL)
  nll <- function(par) {
    df <- exp(par[1]) + 1
    scale <- exp(par[2])
    -sum(log(2) + dt(x / scale, df = df, log = TRUE) - log(scale))
  }
  fit <- optim(c(log(10), log(sd(x))), nll, method = "BFGS")
  if (fit$convergence != 0) return(NULL)
  df <- exp(fit$par[1]) + 1
  scale <- exp(fit$par[2])
  ll <- -fit$value
  list(df = df, scale = scale, loglik = ll, k = 2)
}

fit_gmm2 <- function(x, max_iter = 200, tol = 1e-6) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 10) return(NULL)
  init <- quantile(x, c(0.25, 0.75))
  mu1 <- init[1]
  mu2 <- init[2]
  sd1 <- sd(x)
  sd2 <- sd(x)
  if (!is.finite(sd1) || sd1 == 0 || !is.finite(sd2) || sd2 == 0) return(NULL)
  pi1 <- 0.5
  ll_prev <- -Inf
  for (i in seq_len(max_iter)) {
    dens1 <- pi1 * dnorm(x, mean = mu1, sd = sd1)
    dens2 <- (1 - pi1) * dnorm(x, mean = mu2, sd = sd2)
    denom <- dens1 + dens2
    if (any(!is.finite(denom)) || any(denom <= 0)) return(NULL)
    w1 <- dens1 / denom
    w2 <- 1 - w1
    pi1 <- mean(w1)
    mu1 <- sum(w1 * x) / sum(w1)
    mu2 <- sum(w2 * x) / sum(w2)
    sd1 <- sqrt(sum(w1 * (x - mu1)^2) / sum(w1))
    sd2 <- sqrt(sum(w2 * (x - mu2)^2) / sum(w2))
    if (!is.finite(sd1) || sd1 == 0 || !is.finite(sd2) || sd2 == 0) return(NULL)
    ll <- sum(log(denom))
    if (!is.finite(ll)) return(NULL)
    if (abs(ll - ll_prev) < tol) break
    ll_prev <- ll
  }
  if (!is.finite(ll_prev)) return(NULL)
  list(pi1 = pi1, mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2, loglik = ll_prev, k = 5)
}

ma_summary <- list()
residuals_all <- c()
residuals_by_type <- list()
within_ma_resid <- list()
mixture_summary <- list()
meta_reg_results <- list()
pred_interval_results <- list()

for (file in files) {
  ds_id <- sub("\\.rds$", "", basename(file))
  df <- readRDS(file)
  domain_text <- build_domain_text(df)
  domain_primary <- assign_domain(domain_text)
  intervention_class <- assign_intervention_class(domain_text)

  res <- extract_effects(df)
  if (is.null(res$data)) {
    ma_summary[[ds_id]] <- data.frame(
      ma_id = ds_id,
      k = 0,
      outcome_type = res$outcome_type,
      effect_measure = res$effect_measure,
      domain_primary = domain_primary,
      intervention_class = intervention_class,
      tau2_reml = NA_real_,
      tau2_dl = NA_real_,
      tau2_pm = NA_real_,
      i2 = NA_real_,
      median_se = NA_real_,
      mean_year = NA_real_,
      mean_total_n = NA_real_,
      prop_high_risk = NA_real_,
      stringsAsFactors = FALSE
    )
    next
  }

  dat <- res$data
  yi <- dat[["yi"]]
  vi <- dat[["vi"]]
  k <- length(yi)
  if (k < 2) {
    ma_summary[[ds_id]] <- data.frame(
      ma_id = ds_id,
      k = k,
      outcome_type = res$outcome_type,
      effect_measure = res$effect_measure,
      domain_primary = domain_primary,
      intervention_class = intervention_class,
      tau2_reml = NA_real_,
      tau2_dl = NA_real_,
      tau2_pm = NA_real_,
      i2 = NA_real_,
      median_se = NA_real_,
      mean_year = NA_real_,
      mean_total_n = NA_real_,
      prop_high_risk = NA_real_,
      stringsAsFactors = FALSE
    )
    next
  }

  base_reml <- tryCatch(rma.uni(yi = yi, vi = vi, method = "REML"), error = function(e) NULL)
  base_dl <- tryCatch(rma.uni(yi = yi, vi = vi, method = "DL"), error = function(e) NULL)
  base_pm <- tryCatch(rma.uni(yi = yi, vi = vi, method = "PM"), error = function(e) NULL)

  if (is.null(base_reml)) {
    ma_summary[[ds_id]] <- data.frame(
      ma_id = ds_id,
      k = k,
      outcome_type = res$outcome_type,
      effect_measure = res$effect_measure,
      domain_primary = domain_primary,
      intervention_class = intervention_class,
      tau2_reml = NA_real_,
      tau2_dl = NA_real_,
      tau2_pm = NA_real_,
      i2 = NA_real_,
      median_se = median(sqrt(vi), na.rm = TRUE),
      mean_year = NA_real_,
      mean_total_n = NA_real_,
      prop_high_risk = NA_real_,
      stringsAsFactors = FALSE
    )
    next
  }

  tau2_reml <- base_reml$tau2
  tau2_dl <- if (!is.null(base_dl)) base_dl$tau2 else NA_real_
  tau2_pm <- if (!is.null(base_pm)) base_pm$tau2 else NA_real_

  study_year <- if ("Study.year" %in% names(dat)) suppressWarnings(as.numeric(dat[["Study.year"]])) else NA_real_
  total_n <- rep(NA_real_, k)
  if (all(c("Experimental.N", "Control.N") %in% names(dat))) {
    total_n <- suppressWarnings(as.numeric(dat[["Experimental.N"]]) + as.numeric(dat[["Control.N"]]))
  }
  mean_total_n <- if (all(is.na(total_n))) NA_real_ else mean(total_n, na.rm = TRUE)

  prop_high_risk <- NA_real_
  if ("Overall.bias.judgement" %in% names(dat)) {
    bias <- map_bias_score(dat[["Overall.bias.judgement"]])
    prop_high_risk <- if (all(is.na(bias))) NA_real_ else mean(bias == 1, na.rm = TRUE)
  }

  ma_summary[[ds_id]] <- data.frame(
    ma_id = ds_id,
    k = k,
    outcome_type = res$outcome_type,
    effect_measure = res$effect_measure,
    domain_primary = domain_primary,
    intervention_class = intervention_class,
    tau2_reml = tau2_reml,
    tau2_dl = tau2_dl,
    tau2_pm = tau2_pm,
    i2 = base_reml$I2,
    median_se = median(sqrt(vi), na.rm = TRUE),
    mean_year = if (all(is.na(study_year))) NA_real_ else mean(study_year, na.rm = TRUE),
    mean_total_n = mean_total_n,
    prop_high_risk = prop_high_risk,
    stringsAsFactors = FALSE
  )

  z <- (yi - as.numeric(base_reml$beta)) / sqrt(vi + tau2_reml)
  residuals_all <- c(residuals_all, z)
  residuals_by_type[[res$outcome_type]] <- c(residuals_by_type[[res$outcome_type]], z)

  if (k >= 5) {
    z_sample <- if (length(z) > 5000) sample(z, 5000) else z
    shapiro_p <- tryCatch(shapiro.test(z_sample)$p.value, error = function(e) NA_real_)
    within_ma_resid[[ds_id]] <- data.frame(
      ma_id = ds_id,
      k = k,
      shapiro_p = shapiro_p,
      skewness = skewness(z),
      kurtosis_excess = kurtosis_excess(z),
      prop_abs_gt2 = mean(abs(z) > 2, na.rm = TRUE),
      prop_abs_gt3 = mean(abs(z) > 3, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  if (k >= 8) {
    gmm1_ll <- sum(dnorm(yi, mean = mean(yi), sd = sd(yi), log = TRUE))
    gmm2_fit <- fit_gmm2(yi)
    if (!is.null(gmm2_fit)) {
      bic1 <- -2 * gmm1_ll + 2 * log(k)
      bic2 <- -2 * gmm2_fit$loglik + gmm2_fit$k * log(k)
      mixture_summary[[ds_id]] <- data.frame(
        ma_id = ds_id,
        k = k,
        bic_1comp = bic1,
        bic_2comp = bic2,
        prefer_2comp = bic2 < bic1,
        stringsAsFactors = FALSE
      )
    }
  }

  cov_df <- data.frame(
    year = study_year,
    log_n = if (all(is.na(total_n))) NA_real_ else log(total_n),
    bias_score = if ("Overall.bias.judgement" %in% names(dat)) map_bias_score(dat[["Overall.bias.judgement"]]) else NA_real_
  )

  cov_use <- data.frame(yi = yi, vi = vi)
  mod_terms <- c()
  if (sum(!is.na(cov_df$year)) >= 5 && length(unique(na.omit(cov_df$year))) > 1) {
    cov_use$year_c <- cov_df$year - mean(cov_df$year, na.rm = TRUE)
    mod_terms <- c(mod_terms, "year_c")
  }
  if (sum(!is.na(cov_df$log_n)) >= 5 && length(unique(na.omit(cov_df$log_n))) > 1) {
    cov_use$log_n_c <- cov_df$log_n - mean(cov_df$log_n, na.rm = TRUE)
    mod_terms <- c(mod_terms, "log_n_c")
  }
  if (sum(!is.na(cov_df$bias_score)) >= 5 && length(unique(na.omit(cov_df$bias_score))) > 1) {
    cov_use$bias_score <- cov_df$bias_score
    mod_terms <- c(mod_terms, "bias_score")
  }

  if (length(mod_terms) > 0) {
    cov_use <- cov_use[complete.cases(cov_use[, c("yi", "vi", mod_terms), drop = FALSE]), , drop = FALSE]
    if (nrow(cov_use) >= 5) {
      mod_formula <- as.formula(paste("~", paste(mod_terms, collapse = " + ")))
      mod_fit <- tryCatch(rma.uni(yi = cov_use$yi, vi = cov_use$vi, mods = mod_formula, data = cov_use, method = "REML"),
                          error = function(e) NULL)
      if (!is.null(mod_fit)) {
        tau2_mod <- mod_fit$tau2
        r2_tau <- if (is.finite(tau2_reml) && tau2_reml > 0) max(0, 1 - tau2_mod / tau2_reml) else NA_real_

        term_names <- colnames(mod_fit$X)
        if (is.null(term_names) || length(term_names) != length(mod_fit$beta)) {
          term_names <- paste0("beta", seq_along(mod_fit$beta))
        }
        coef_df <- data.frame(
          ma_id = ds_id,
          term = term_names,
          estimate = as.numeric(mod_fit$beta),
          se = mod_fit$se,
          tau2_base = tau2_reml,
          tau2_mod = tau2_mod,
          r2_tau = r2_tau,
          stringsAsFactors = FALSE
        )
        meta_reg_results[[ds_id]] <- coef_df

        zcrit <- qnorm(0.975)
        median_vi <- median(cov_use$vi, na.rm = TRUE)
        pi_width_standard <- 2 * zcrit * sqrt(tau2_reml + median_vi)
        pi_width_conditional <- 2 * zcrit * sqrt(tau2_mod + median_vi)

        cover_std <- rep(NA, nrow(cov_use))
        cover_cond <- rep(NA, nrow(cov_use))
        for (idx in seq_len(nrow(cov_use))) {
          train <- cov_use[-idx, , drop = FALSE]
          test <- cov_use[idx, , drop = FALSE]
          base_cv <- tryCatch(rma.uni(yi = train$yi, vi = train$vi, method = "REML"), error = function(e) NULL)
          mod_cv <- tryCatch(rma.uni(yi = train$yi, vi = train$vi, mods = mod_formula, data = train, method = "REML"),
                             error = function(e) NULL)
          if (!is.null(base_cv)) {
            mu_hat <- as.numeric(base_cv$beta)
            se_pi <- sqrt(test$vi + base_cv$tau2)
            pi_low <- mu_hat - zcrit * se_pi
            pi_high <- mu_hat + zcrit * se_pi
            cover_std[idx] <- test$yi >= pi_low & test$yi <= pi_high
          }
          if (!is.null(mod_cv)) {
            newmods <- as.matrix(test[, mod_terms, drop = FALSE])
            pred <- tryCatch(predict(mod_cv, newmods = newmods)$pred, error = function(e) NA_real_)
            se_pi <- sqrt(test$vi + mod_cv$tau2)
            pi_low <- pred - zcrit * se_pi
            pi_high <- pred + zcrit * se_pi
            cover_cond[idx] <- test$yi >= pi_low & test$yi <= pi_high
          }
        }

        pred_interval_results[[ds_id]] <- data.frame(
          ma_id = ds_id,
          k = nrow(cov_use),
          pi_width_standard = pi_width_standard,
          pi_width_conditional = pi_width_conditional,
          width_ratio = pi_width_conditional / pi_width_standard,
          cv_coverage_standard = if (all(is.na(cover_std))) NA_real_ else mean(cover_std, na.rm = TRUE),
          cv_coverage_conditional = if (all(is.na(cover_cond))) NA_real_ else mean(cover_cond, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

ma_summary_df <- do.call(rbind, ma_summary)
write.csv(ma_summary_df, file.path(out_dir, "tables", "ma_summary.csv"), row.names = FALSE)

residuals_all <- residuals_all[is.finite(residuals_all)]
resid_sample <- if (length(residuals_all) > 20000) sample(residuals_all, 20000) else residuals_all

resid_summary <- data.frame(
  n = length(residuals_all),
  mean = mean(residuals_all),
  sd = sd(residuals_all),
  skewness = skewness(residuals_all),
  kurtosis_excess = kurtosis_excess(residuals_all),
  prop_abs_gt2 = mean(abs(residuals_all) > 2),
  prop_abs_gt3 = mean(abs(residuals_all) > 3),
  shapiro_p = if (length(resid_sample) >= 3) tryCatch(shapiro.test(resid_sample)$p.value, error = function(e) NA_real_) else NA_real_,
  stringsAsFactors = FALSE
)
write.csv(resid_summary, file.path(out_dir, "tables", "residual_summary.csv"), row.names = FALSE)

if (length(residuals_all) > 5) {
  png(file.path(out_dir, "figures", "residual_hist.png"), width = 900, height = 600)
  hist(residuals_all, breaks = 50, main = "Pooled standardized residuals", xlab = "z", col = "steelblue", border = "white")
  abline(v = 0, col = "red", lwd = 2, lty = 2)
  dev.off()

  png(file.path(out_dir, "figures", "residual_qq.png"), width = 900, height = 600)
  qqnorm(residuals_all, main = "QQ plot of pooled standardized residuals")
  qqline(residuals_all, col = "red", lwd = 2)
  dev.off()
}

within_ma_resid_df <- if (length(within_ma_resid) > 0) do.call(rbind, within_ma_resid) else data.frame()
write.csv(within_ma_resid_df, file.path(out_dir, "tables", "within_ma_residuals.csv"), row.names = FALSE)

mixture_df <- if (length(mixture_summary) > 0) do.call(rbind, mixture_summary) else data.frame()
write.csv(mixture_df, file.path(out_dir, "tables", "mixture_model_summary.csv"), row.names = FALSE)

tau_vals <- sqrt(ma_summary_df$tau2_reml)
tau_vals <- tau_vals[is.finite(tau_vals) & tau_vals >= 0]
tau_vals_pos <- tau_vals[tau_vals > 0]
tau_summary <- data.frame(
  n = length(tau_vals),
  median = median(tau_vals),
  iqr = IQR(tau_vals),
  p90 = quantile(tau_vals, 0.9),
  stringsAsFactors = FALSE
)

if (length(tau_vals) > 10) {
  ln_fit <- tryCatch(MASS::fitdistr(tau_vals_pos, "lognormal"), error = function(e) NULL)
  gamma_fit <- tryCatch(MASS::fitdistr(tau_vals_pos, "gamma"), error = function(e) NULL)
  hn_fit <- fit_half_normal(tau_vals)
  ht_fit <- fit_half_t(tau_vals)

  dist_fits <- data.frame(
    dist = character(),
    loglik = numeric(),
    k = integer(),
    aic = numeric(),
    stringsAsFactors = FALSE
  )
  if (!is.null(ln_fit)) {
    dist_fits <- rbind(dist_fits, data.frame(dist = "lognormal", loglik = ln_fit$loglik, k = 2,
                                             aic = -2 * ln_fit$loglik + 2 * 2))
  }
  if (!is.null(gamma_fit)) {
    dist_fits <- rbind(dist_fits, data.frame(dist = "gamma", loglik = gamma_fit$loglik, k = 2,
                                             aic = -2 * gamma_fit$loglik + 2 * 2))
  }
  if (!is.null(hn_fit)) {
    dist_fits <- rbind(dist_fits, data.frame(dist = "half_normal", loglik = hn_fit$loglik, k = hn_fit$k,
                                             aic = -2 * hn_fit$loglik + 2 * hn_fit$k))
  }
  if (!is.null(ht_fit)) {
    dist_fits <- rbind(dist_fits, data.frame(dist = "half_t", loglik = ht_fit$loglik, k = ht_fit$k,
                                             aic = -2 * ht_fit$loglik + 2 * ht_fit$k))
  }

  write.csv(dist_fits, file.path(out_dir, "tables", "tau_distribution_fits.csv"), row.names = FALSE)
}

if (length(tau_vals) > 5) {
  png(file.path(out_dir, "figures", "tau_distribution.png"), width = 900, height = 600)
  hist(tau_vals, breaks = 40, main = "Distribution of tau (REML)", xlab = "tau", col = "forestgreen", border = "white")
  dev.off()
}

tau_by_type <- aggregate(tau2_reml ~ outcome_type, data = ma_summary_df, FUN = function(x) {
  x <- sqrt(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
})
write.csv(tau_by_type, file.path(out_dir, "tables", "tau_by_outcome_type.csv"), row.names = FALSE)

if (nrow(ma_summary_df) > 0) {
  base_df <- ma_summary_df
  base_df <- base_df[is.finite(base_df$tau2_reml) & is.finite(base_df$k) &
                       is.finite(base_df$median_se) & base_df$k > 0 & base_df$median_se > 0, , drop = FALSE]
  base_model <- lm(log(pmax(sqrt(base_df$tau2_reml), 1e-6)) ~ outcome_type + log(k + 1) +
                     log(median_se + 1e-6) + domain_primary + intervention_class,
                   data = base_df)
  base_coef <- summary(base_model)$coefficients
  base_coef_df <- data.frame(
    term = rownames(base_coef),
    estimate = base_coef[, 1],
    se = base_coef[, 2],
    t = base_coef[, 3],
    p = base_coef[, 4],
    model = "base",
    stringsAsFactors = FALSE
  )
  write.csv(base_coef_df, file.path(out_dir, "tables", "tau_regression_base.csv"), row.names = FALSE)
}

write.csv(
  unique(ma_summary_df[, c("ma_id", "domain_primary", "intervention_class")]),
  file.path(out_dir, "tables", "ma_domain_mapping.csv"),
  row.names = FALSE
)

for (otype in names(residuals_by_type)) {
  zvals <- residuals_by_type[[otype]]
  zvals <- zvals[is.finite(zvals)]
  if (length(zvals) < 10) next
  z_sample <- if (length(zvals) > 20000) sample(zvals, 20000) else zvals

  resid_sum <- data.frame(
    outcome_type = otype,
    n = length(zvals),
    mean = mean(zvals),
    sd = sd(zvals),
    skewness = skewness(zvals),
    kurtosis_excess = kurtosis_excess(zvals),
    prop_abs_gt2 = mean(abs(zvals) > 2),
    prop_abs_gt3 = mean(abs(zvals) > 3),
    shapiro_p = if (length(z_sample) >= 3) tryCatch(shapiro.test(z_sample)$p.value, error = function(e) NA_real_) else NA_real_,
    stringsAsFactors = FALSE
  )
  write.csv(resid_sum,
            file.path(out_dir, "tables", paste0("residual_summary_", otype, ".csv")),
            row.names = FALSE)

  png(file.path(out_dir, "figures", paste0("residual_hist_", otype, ".png")), width = 900, height = 600)
  hist(zvals, breaks = 50, main = paste("Pooled standardized residuals:", otype),
       xlab = "z", col = "steelblue", border = "white")
  abline(v = 0, col = "red", lwd = 2, lty = 2)
  dev.off()

  png(file.path(out_dir, "figures", paste0("residual_qq_", otype, ".png")), width = 900, height = 600)
  qqnorm(zvals, main = paste("QQ plot of pooled residuals:", otype))
  qqline(zvals, col = "red", lwd = 2)
  dev.off()

  dist_res_type <- list()
  norm_ll <- sum(dnorm(z_sample, mean = mean(z_sample), sd = sd(z_sample), log = TRUE))
  dist_res_type[["normal"]] <- list(loglik = norm_ll, k = 2)
  t_fit <- fit_t_distribution(z_sample)
  if (!is.null(t_fit)) dist_res_type[["t"]] <- t_fit
  sn_fit <- fit_skew_normal(z_sample)
  if (!is.null(sn_fit)) dist_res_type[["skew_normal"]] <- sn_fit
  mix_fit <- fit_gmm2(z_sample)
  if (!is.null(mix_fit)) dist_res_type[["mixture_2norm"]] <- mix_fit

  dist_table <- do.call(rbind, lapply(names(dist_res_type), function(nm) {
    ll <- dist_res_type[[nm]]$loglik
    k <- dist_res_type[[nm]]$k
    data.frame(
      outcome_type = otype,
      dist = nm,
      loglik = ll,
      k = k,
      aic = -2 * ll + 2 * k,
      stringsAsFactors = FALSE
    )
  }))
  write.csv(dist_table,
            file.path(out_dir, "tables", paste0("residual_distribution_fits_", otype, ".csv")),
            row.names = FALSE)
}

for (otype in unique(ma_summary_df$outcome_type)) {
  sub <- ma_summary_df[ma_summary_df$outcome_type == otype, , drop = FALSE]
  tau_vals_type <- sqrt(sub$tau2_reml)
  tau_vals_type <- tau_vals_type[is.finite(tau_vals_type) & tau_vals_type >= 0]
  if (length(tau_vals_type) < 10) next
  tau_pos_type <- tau_vals_type[tau_vals_type > 0]

  tau_sum <- data.frame(
    outcome_type = otype,
    n = length(tau_vals_type),
    median = median(tau_vals_type),
    iqr = IQR(tau_vals_type),
    p90 = quantile(tau_vals_type, 0.9),
    stringsAsFactors = FALSE
  )
  write.csv(tau_sum,
            file.path(out_dir, "tables", paste0("tau_summary_", otype, ".csv")),
            row.names = FALSE)

  dist_fits_type <- data.frame(dist = character(), loglik = numeric(), k = integer(), aic = numeric(),
                               stringsAsFactors = FALSE)
  ln_fit <- tryCatch(MASS::fitdistr(tau_pos_type, "lognormal"), error = function(e) NULL)
  gamma_fit <- tryCatch(MASS::fitdistr(tau_pos_type, "gamma"), error = function(e) NULL)
  hn_fit <- fit_half_normal(tau_vals_type)
  ht_fit <- fit_half_t(tau_vals_type)

  if (!is.null(ln_fit)) {
    dist_fits_type <- rbind(dist_fits_type, data.frame(dist = "lognormal", loglik = ln_fit$loglik, k = 2,
                                                       aic = -2 * ln_fit$loglik + 4))
  }
  if (!is.null(gamma_fit)) {
    dist_fits_type <- rbind(dist_fits_type, data.frame(dist = "gamma", loglik = gamma_fit$loglik, k = 2,
                                                       aic = -2 * gamma_fit$loglik + 4))
  }
  if (!is.null(hn_fit)) {
    dist_fits_type <- rbind(dist_fits_type, data.frame(dist = "half_normal", loglik = hn_fit$loglik, k = hn_fit$k,
                                                       aic = -2 * hn_fit$loglik + 2 * hn_fit$k))
  }
  if (!is.null(ht_fit)) {
    dist_fits_type <- rbind(dist_fits_type, data.frame(dist = "half_t", loglik = ht_fit$loglik, k = ht_fit$k,
                                                       aic = -2 * ht_fit$loglik + 2 * ht_fit$k))
  }
  write.csv(dist_fits_type,
            file.path(out_dir, "tables", paste0("tau_distribution_fits_", otype, ".csv")),
            row.names = FALSE)

  png(file.path(out_dir, "figures", paste0("tau_distribution_", otype, ".png")), width = 900, height = 600)
  hist(tau_vals_type, breaks = 40, main = paste("Distribution of tau (REML):", otype),
       xlab = "tau", col = "forestgreen", border = "white")
  dev.off()
}

meta_reg_df <- if (length(meta_reg_results) > 0) {
  do.call(rbind, meta_reg_results)
} else {
  data.frame(
    ma_id = character(),
    term = character(),
    estimate = numeric(),
    se = numeric(),
    tau2_base = numeric(),
    tau2_mod = numeric(),
    r2_tau = numeric(),
    stringsAsFactors = FALSE
  )
}
write.csv(meta_reg_df, file.path(out_dir, "tables", "meta_regression_coefficients.csv"), row.names = FALSE)

if (nrow(meta_reg_df) > 0) {
  pooled_list <- lapply(split(meta_reg_df, meta_reg_df$term), function(d) {
    if (nrow(d) < 3) return(NULL)
    ma <- tryCatch(rma.uni(yi = d$estimate, sei = d$se, method = "REML"), error = function(e) NULL)
    if (is.null(ma)) return(NULL)
    data.frame(
      term = unique(d$term),
      pooled_estimate = as.numeric(ma$beta),
      pooled_se = ma$se,
      tau2 = ma$tau2,
      k = ma$k,
      stringsAsFactors = FALSE
    )
  })
  pooled_df <- do.call(rbind, pooled_list)
  write.csv(pooled_df, file.path(out_dir, "tables", "meta_regression_pooled.csv"), row.names = FALSE)
}

pred_df <- if (length(pred_interval_results) > 0) {
  do.call(rbind, pred_interval_results)
} else {
  data.frame(
    ma_id = character(),
    k = integer(),
    pi_width_standard = numeric(),
    pi_width_conditional = numeric(),
    width_ratio = numeric(),
    cv_coverage_standard = numeric(),
    cv_coverage_conditional = numeric(),
    stringsAsFactors = FALSE
  )
}
write.csv(pred_df, file.path(out_dir, "tables", "prediction_interval_summary.csv"), row.names = FALSE)

dist_res <- list()
if (length(resid_sample) >= 10) {
  norm_ll <- sum(dnorm(resid_sample, mean = mean(resid_sample), sd = sd(resid_sample), log = TRUE))
  dist_res[["normal"]] <- list(loglik = norm_ll, k = 2)
  t_fit <- fit_t_distribution(resid_sample)
  if (!is.null(t_fit)) dist_res[["t"]] <- t_fit
  sn_fit <- fit_skew_normal(resid_sample)
  if (!is.null(sn_fit)) dist_res[["skew_normal"]] <- sn_fit
  mix_fit <- fit_gmm2(resid_sample)
  if (!is.null(mix_fit)) dist_res[["mixture_2norm"]] <- mix_fit
}

if (length(dist_res) > 0) {
  dist_table <- do.call(rbind, lapply(names(dist_res), function(nm) {
    ll <- dist_res[[nm]]$loglik
    k <- dist_res[[nm]]$k
    data.frame(
      dist = nm,
      loglik = ll,
      k = k,
      aic = -2 * ll + 2 * k,
      stringsAsFactors = FALSE
    )
  }))
  write.csv(dist_table, file.path(out_dir, "tables", "residual_distribution_fits.csv"), row.names = FALSE)
}

write.csv(tau_summary, file.path(out_dir, "tables", "tau_summary.csv"), row.names = FALSE)
