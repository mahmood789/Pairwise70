library(metafor)
files <- list.files('analysis/output/cleaned_rds', pattern='\.rds$', full.names=TRUE)
count_mod <- 0
count_attempt <- 0
for (f in files) {
  df <- readRDS(f)
  # quick effect detection: use binary if possible
  has_bin <- all(c('Experimental.cases','Experimental.N','Control.cases','Control.N') %in% names(df))
  bin_has_data <- has_bin && (any(!is.na(df[['Experimental.cases']])) || any(!is.na(df[['Control.cases']])) )
  if (!bin_has_data) next
  use <- complete.cases(df[, c('Experimental.cases','Experimental.N','Control.cases','Control.N')])
  d <- df[use, , drop=FALSE]
  if (nrow(d) < 5) next
  valid <- d[['Experimental.cases']] >= 0 & d[['Control.cases']] >= 0 & d[['Experimental.N']] > 0 & d[['Control.N']] > 0 &
    d[['Experimental.cases']] <= d[['Experimental.N']] & d[['Control.cases']] <= d[['Control.N']]
  d <- d[valid, , drop=FALSE]
  if (nrow(d) < 5) next
  es <- escalc(measure='OR', ai=d[['Experimental.cases']], n1i=d[['Experimental.N']], ci=d[['Control.cases']], n2i=d[['Control.N']], add=0.5, to='all')
  yi <- es[['yi']]; vi <- es[['vi']]
  year <- suppressWarnings(as.numeric(es[['Study.year']]))
  cov_use <- data.frame(yi=yi, vi=vi)
  mod_terms <- c()
  if (sum(!is.na(year)) >= 5 && length(unique(na.omit(year))) > 1) {
    cov_use$year_c <- year - mean(year, na.rm=TRUE)
    mod_terms <- c(mod_terms, 'year_c')
  }
  if (length(mod_terms) == 0) next
  cov_use <- cov_use[complete.cases(cov_use[, c('yi','vi', mod_terms), drop=FALSE]), , drop=FALSE]
  if (nrow(cov_use) < 5) next
  count_attempt <- count_attempt + 1
  mod_formula <- as.formula(paste('~', paste(mod_terms, collapse=' + ')))
  mod_fit <- tryCatch(rma.uni(yi=cov_use$yi, vi=cov_use$vi, mods=mod_formula, data=cov_use, method='REML'), error=function(e) NULL)
  if (!is.null(mod_fit)) count_mod <- count_mod + 1
}
cat('attempted:', count_attempt, '\n')
cat('succeeded:', count_mod, '\n')
