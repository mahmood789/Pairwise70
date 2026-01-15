library(MASS)
fit_half_normal <- function(x) {
  x <- x[is.finite(x) & x >= 0]
  if (length(x) < 10) return(NULL)
  nll <- function(par) {
    sigma <- exp(par[1])
    -sum(log(sqrt(2 / pi)) - log(sigma) - (x^2) / (2 * sigma^2))
  }
  fit <- optim(log(sd(x)), nll, method = "BFGS")
  if (fit$convergence != 0) return(NULL)
  ll <- -fit$value
  list(loglik = ll, k = 1)
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
  ll <- -fit$value
  list(loglik = ll, k = 2)
}

d <- read.csv('analysis/heterogeneity_output/tables/ma_summary.csv')
tau <- sqrt(d[['tau2_reml']])
tau <- tau[is.finite(tau) & tau >= 0]
tau_pos <- tau[tau > 0]

dist_fits <- data.frame(dist = character(), loglik = numeric(), k = integer(), aic = numeric(), stringsAsFactors = FALSE)
ln_fit <- tryCatch(fitdistr(tau_pos, 'lognormal'), error = function(e) NULL)
gamma_fit <- tryCatch(fitdistr(tau_pos, 'gamma'), error = function(e) NULL)
hn_fit <- fit_half_normal(tau)
ht_fit <- fit_half_t(tau)

if (!is.null(ln_fit)) {
  dist_fits <- rbind(dist_fits, data.frame(dist = 'lognormal', loglik = ln_fit$loglik, k = 2, aic = -2 * ln_fit$loglik + 4))
}
if (!is.null(gamma_fit)) {
  dist_fits <- rbind(dist_fits, data.frame(dist = 'gamma', loglik = gamma_fit$loglik, k = 2, aic = -2 * gamma_fit$loglik + 4))
}
if (!is.null(hn_fit)) {
  dist_fits <- rbind(dist_fits, data.frame(dist = 'half_normal', loglik = hn_fit$loglik, k = hn_fit$k, aic = -2 * hn_fit$loglik + 2 * hn_fit$k))
}
if (!is.null(ht_fit)) {
  dist_fits <- rbind(dist_fits, data.frame(dist = 'half_t', loglik = ht_fit$loglik, k = ht_fit$k, aic = -2 * ht_fit$loglik + 2 * ht_fit$k))
}

write.csv(dist_fits, 'analysis/heterogeneity_output/tables/tau_distribution_fits.csv', row.names = FALSE)
