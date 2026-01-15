# MA4 v1.0.1 Analysis for Pairwise70 Dataset Collection
# Runs MA4 on all 501 Cochrane review datasets

# ============================================================================
# MA4 v1.0.1 Core Implementation
# ============================================================================

MA4_EPS <- 1e-12

sgn <- function(x, eps = MA4_EPS) {
  if (is.na(x)) return(NA_integer_)
  if (abs(x) <= eps) return(0L)
  if (x > 0) return(1L)
  -1L
}

var_ddof1 <- function(x) {
  n <- length(x)
  if (n < 2) return(0)
  m <- mean(x)
  sum((x - m)^2) / (n - 1)
}

theta_sigma_q_nll <- function(y, se, tau2) {
  v <- se^2 + tau2
  w <- 1 / v
  sum_w <- sum(w)
  theta <- sum(w * y) / sum_w
  r <- y - theta
  Q <- sum(w * r * r)
  sigma <- sqrt(1 / sum_w)
  nll <- sum(log(v)) + log(sum_w) + Q
  list(theta = theta, sigma = sigma, Q = Q, nll = nll)
}

tau2_reml_golden <- function(y, se, tol = 1e-10, max_iter = 200) {
  k <- length(y)
  if (k < 2) return(0)

  tmax <- max(1e-6, 10 * var_ddof1(y))
  a <- 0
  b <- tmax

  gr <- (sqrt(5) - 1) / 2
  c <- b - gr * (b - a)
  d <- a + gr * (b - a)

  fc <- theta_sigma_q_nll(y, se, c)$nll
  fd <- theta_sigma_q_nll(y, se, d)$nll

  for (it in seq_len(max_iter)) {
    if (abs(b - a) <= tol * (1 + abs(a) + abs(b))) break
    if (fc < fd) {
      b <- d; d <- c; fd <- fc
      c <- b - gr * (b - a)
      fc <- theta_sigma_q_nll(y, se, c)$nll
    } else {
      a <- c; c <- d; fc <- fd
      d <- a + gr * (b - a)
      fd <- theta_sigma_q_nll(y, se, d)$nll
    }
  }
  max(0, 0.5 * (a + b))
}

tau2_dl <- function(y, se) {
  k <- length(y)
  if (k < 2) return(0)

  w <- 1 / (se^2)
  sum_w <- sum(w)
  theta_fe <- sum(w * y) / sum_w
  Q <- sum(w * (y - theta_fe)^2)

  C <- sum_w - (sum(w^2) / sum_w)
  if (C <= 0) return(0)
  max(0, (Q - (k - 1)) / C)
}

tau2_pm <- function(y, se, tol = 1e-10, max_iter = 200) {
  k <- length(y)
  if (k < 2) return(0)
  target <- k - 1

  Qof <- function(t) theta_sigma_q_nll(y, se, t)$Q

  a <- 0
  Qa <- Qof(a)
  if (Qa <= target) return(0)

  b <- max(1e-6, 10 * var_ddof1(y))
  Qb <- Qof(b)

  while (Qb > target && b < 1e6) {
    b <- b * 2
    Qb <- Qof(b)
  }
  if (Qb > target) return(b)

  for (it in seq_len(max_iter)) {
    if (abs(b - a) <= tol * (1 + abs(a) + abs(b))) break
    m <- 0.5 * (a + b)
    Qm <- Qof(m)
    if (Qm > target) a <- m else b <- m
  }
  0.5 * (a + b)
}

tau2_baseline <- function(y, se) {
  t <- tau2_reml_golden(y, se)
  used <- "REML_golden_v1"
  if (!is.finite(t)) {
    t <- tau2_pm(y, se)
    used <- "PM_bisection_v1"
  }
  if (!is.finite(t)) {
    t <- tau2_dl(y, se)
    used <- "DL_closed_v1"
  }
  if (!is.finite(t)) {
    t <- 0
    used <- "zero"
  }
  list(tau2 = t, used = used)
}

compute_ma4 <- function(y, se) {
  stopifnot(length(y) == length(se), length(y) >= 1)

  k <- length(y)
  if (k == 1) {
    return(list(
      theta = y[1], sigma = se[1], tau = 0, R = 0,
      tau_estimator_used = "zero",
      R_status = "undefined_k1"
    ))
  }

  base <- tau2_baseline(y, se)
  tau2 <- base$tau2
  used <- base$used

  base_stats <- theta_sigma_q_nll(y, se, tau2)
  theta <- base_stats$theta
  sigma <- base_stats$sigma
  tau <- sqrt(max(0, tau2))

  # Perturbation set Π v1
  thetas <- c(
    theta_sigma_q_nll(y, se, 0)$theta,
    theta_sigma_q_nll(y, se, tau2_dl(y, se))$theta,
    theta_sigma_q_nll(y, se, tau2_pm(y, se))$theta
  )

  # LOO-worst
  if (k >= 3) {
    worst_dev <- -Inf
    worst_th <- NA_real_
    for (j in seq_len(k)) {
      idx <- setdiff(seq_len(k), j)
      y2 <- y[idx]; se2 <- se[idx]
      t2 <- tau2_baseline(y2, se2)$tau2
      th2 <- theta_sigma_q_nll(y2, se2, t2)$theta
      dev <- abs(th2 - theta)
      if (dev > worst_dev) { worst_dev <- dev; worst_th <- th2 }
    }
    thetas <- c(thetas, worst_th)
  }

  # Drop-noisiest-20%
  if (k >= 3) {
    idx <- order(se, seq_len(k))
    m <- max(1, ceiling(0.2 * k))
    keep <- idx[seq_len(k - m)]
    if (length(keep) >= 2) {
      y3 <- y[keep]; se3 <- se[keep]
      t3 <- tau2_baseline(y3, se3)$tau2
      thetas <- c(thetas, theta_sigma_q_nll(y3, se3, t3)$theta)
    }
  }

  if (!(sigma > 0)) {
    return(list(
      theta = theta, sigma = sigma, tau = tau, R = 0,
      tau_estimator_used = used,
      R_status = "sigma_nonpositive"
    ))
  }

  denom <- 1.96 * sigma
  delta <- max(abs(thetas - theta) / denom)
  R0 <- 1 / (1 + delta)

  base_sign <- sgn(theta)
  flip <- (base_sign == 0)
  if (!flip) {
    for (th in thetas) {
      st <- sgn(th)
      if (st != 0 && st != base_sign) { flip <- TRUE; break }
    }
  }
  F <- if (flip) 0.5 else 1.0
  R <- max(0, min(1, R0 * F))

  list(
    theta = theta, sigma = sigma, tau = tau, R = R,
    tau_estimator_used = used,
    R_status = "ok"
  )
}

# ============================================================================
# Cochrane-Specific Data Extraction
# ============================================================================

# Detect outcome type and extract y, se
extract_y_se <- function(df) {
  # Check if GIV.Mean and GIV.SE are available and valid
  if ("GIV.Mean" %in% names(df) && "GIV.SE" %in% names(df)) {
    y <- as.numeric(df$GIV.Mean)
    se <- as.numeric(df$GIV.SE)
    valid <- !is.na(y) & !is.na(se) & is.finite(y) & is.finite(se) & se > 0
    if (sum(valid) >= 2) {
      return(list(y = y[valid], se = se[valid], type = "GIV", valid_idx = which(valid)))
    }
  }

  # Check for binary outcome (events/totals)
  has_binary <- all(c("Experimental.cases", "Experimental.N",
                       "Control.cases", "Control.N") %in% names(df))
  if (has_binary) {
    et <- as.numeric(df$Experimental.cases)
    nt <- as.numeric(df$Experimental.N)
    ec <- as.numeric(df$Control.cases)
    nc <- as.numeric(df$Control.N)

    valid <- !is.na(et) & !is.na(nt) & !is.na(ec) & !is.na(nc) &
             nt > 0 & nc > 0 & et >= 0 & ec >= 0 & et <= nt & ec <= nc

    if (sum(valid) >= 2) {
      et <- et[valid]; nt <- nt[valid]; ec <- ec[valid]; nc <- nc[valid]

      # Compute log(RR) with continuity correction
      a <- et; b <- nt - et; c <- ec; d <- nc - ec
      need_cc <- (a == 0) | (b == 0) | (c == 0) | (d == 0)
      cc <- ifelse(need_cc, 0.5, 0.0)
      a <- a + cc; b <- b + cc; c <- c + cc; d <- d + cc

      rr <- (a / (a + b)) / (c / (c + d))
      y <- log(rr)
      se <- sqrt(1/a - 1/(a + b) + 1/c - 1/(c + d))

      finite <- is.finite(y) & is.finite(se) & se > 0
      if (sum(finite) >= 2) {
        return(list(y = y[finite], se = se[finite], type = "logRR", valid_idx = which(valid)[finite]))
      }
    }
  }

  # Check for continuous outcome (mean/SD/N)
  has_cont <- all(c("Experimental.mean", "Experimental.SD", "Experimental.N",
                    "Control.mean", "Control.SD", "Control.N") %in% names(df))
  if (has_cont) {
    mt <- as.numeric(df$Experimental.mean)
    sdt <- as.numeric(df$Experimental.SD)
    nt <- as.numeric(df$Experimental.N)
    mc <- as.numeric(df$Control.mean)
    sdc <- as.numeric(df$Control.SD)
    nc <- as.numeric(df$Control.N)

    valid <- !is.na(mt) & !is.na(sdt) & !is.na(nt) &
             !is.na(mc) & !is.na(sdc) & !is.na(nc) &
             sdt >= 0 & sdc >= 0 & nt > 0 & nc > 0

    if (sum(valid) >= 2) {
      mt <- mt[valid]; sdt <- sdt[valid]; nt <- nt[valid]
      mc <- mc[valid]; sdc <- sdc[valid]; nc <- nc[valid]

      y <- mt - mc
      se <- sqrt((sdt^2)/nt + (sdc^2)/nc)

      finite <- is.finite(y) & is.finite(se) & se > 0
      if (sum(finite) >= 2) {
        return(list(y = y[finite], se = se[finite], type = "MD", valid_idx = which(valid)[finite]))
      }
    }
  }

  NULL
}

# ============================================================================
# Main Analysis Runner
# ============================================================================

run_ma4_on_pairwise70 <- function(data_dir, output_file = "ma4_results.csv") {

  rda_files <- list.files(data_dir, pattern = "\\.rda$", full.names = TRUE)
  cat("Found", length(rda_files), "Cochrane review datasets\n")

  results <- list()
  n_meta_analyses <- 0
  n_errors <- 0

  for (i in seq_along(rda_files)) {
    file <- rda_files[i]
    review_id <- gsub("_data\\.rda$", "", basename(file))

    # Load dataset
    e <- new.env()
    tryCatch({
      load(file, envir = e)
    }, error = function(err) {
      cat("Error loading", file, ":", err$message, "\n")
      n_errors <<- n_errors + 1
    })

    objs <- ls(envir = e)
    if (length(objs) == 0) next

    df <- get(objs[1], envir = e)
    if (!is.data.frame(df)) next

    # Get DOI if available
    doi <- if ("review_doi" %in% names(df)) df$review_doi[1] else NA

    # Split by Analysis.number (each represents a distinct meta-analysis)
    if ("Analysis.number" %in% names(df)) {
      analyses <- split(df, df$Analysis.number)
    } else {
      analyses <- list(df)
    }

    for (analysis_num in names(analyses)) {
      analysis_df <- analyses[[analysis_num]]

      # Further split by Applicability to get only OVERALL or single-group analyses
      # Focus on studies that are part of the main analysis
      if ("Applicability" %in% names(analysis_df)) {
        main_df <- analysis_df[grepl("OVERALL|^$", analysis_df$Applicability, ignore.case = TRUE), ]
        if (nrow(main_df) < 2) main_df <- analysis_df
      } else {
        main_df <- analysis_df
      }

      # Remove duplicate studies (keep unique Study entries)
      if ("Study" %in% names(main_df)) {
        main_df <- main_df[!duplicated(main_df$Study), ]
      }

      # Extract y and se
      extraction <- extract_y_se(main_df)
      if (is.null(extraction)) next
      if (length(extraction$y) < 2) next

      n_meta_analyses <- n_meta_analyses + 1

      # Run MA4
      ma4_result <- tryCatch({
        compute_ma4(extraction$y, extraction$se)
      }, error = function(err) {
        list(theta = NA, sigma = NA, tau = NA, R = NA,
             tau_estimator_used = "error", R_status = err$message)
      })

      # Get analysis name
      analysis_name <- if ("Analysis.name" %in% names(main_df)) main_df$Analysis.name[1] else NA

      results[[length(results) + 1]] <- data.frame(
        review_id = review_id,
        analysis_number = analysis_num,
        analysis_name = analysis_name,
        doi = doi,
        k = length(extraction$y),
        effect_type = extraction$type,
        theta = ma4_result$theta,
        sigma = ma4_result$sigma,
        tau = ma4_result$tau,
        R = ma4_result$R,
        tau_estimator = ma4_result$tau_estimator_used,
        R_status = ma4_result$R_status,
        stringsAsFactors = FALSE
      )
    }

    rm(e)
    if (i %% 50 == 0) {
      cat("Processed", i, "/", length(rda_files), "files,",
          n_meta_analyses, "meta-analyses so far\n")
    }
  }

  # Combine results
  results_df <- do.call(rbind, results)

  cat("\n=== MA4 Analysis Complete ===\n")
  cat("Reviews processed:", length(rda_files), "\n")
  cat("Meta-analyses found:", n_meta_analyses, "\n")
  cat("Errors:", n_errors, "\n")

  # Save results
  write.csv(results_df, output_file, row.names = FALSE)
  cat("Results saved to:", output_file, "\n")

  # Summary statistics
  cat("\n=== R Stability Distribution ===\n")
  R_values <- results_df$R[!is.na(results_df$R)]
  cat("Min R:", min(R_values), "\n")
  cat("Median R:", median(R_values), "\n")
  cat("Mean R:", mean(R_values), "\n")
  cat("Max R:", max(R_values), "\n")
  cat("R >= 0.8:", sum(R_values >= 0.8), "(", round(100*mean(R_values >= 0.8), 1), "%)\n")
  cat("R >= 0.5:", sum(R_values >= 0.5), "(", round(100*mean(R_values >= 0.5), 1), "%)\n")
  cat("R < 0.5:", sum(R_values < 0.5), "(", round(100*mean(R_values < 0.5), 1), "%)\n")

  invisible(results_df)
}

# ============================================================================
# RUN ANALYSIS
# ============================================================================

data_dir <- "C:/Users/user/OneDrive - NHS/Documents/Pairwise70/data"
output_dir <- "C:/Users/user/OneDrive - NHS/Documents/Pairwise70/analysis"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_file <- file.path(output_dir, "ma4_results_pairwise70.csv")

cat("Starting MA4 v1.0.1 Analysis on Pairwise70 Dataset Collection\n")
cat("============================================================\n\n")

results <- run_ma4_on_pairwise70(data_dir, output_file)

# Additional summary by effect type
cat("\n=== Results by Effect Type ===\n")
for (et in unique(results$effect_type)) {
  subset_r <- results$R[results$effect_type == et & !is.na(results$R)]
  cat(et, ": n =", length(subset_r),
      ", mean R =", round(mean(subset_r), 3),
      ", median R =", round(median(subset_r), 3), "\n")
}
