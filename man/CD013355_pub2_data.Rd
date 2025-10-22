\name{CD013355_pub2_data}
\alias{CD013355_pub2_data}
\docType{data}
\title{Cochrane Review CD013355 - Pairwise Meta-Analysis Data}
\description{
Pairwise meta-analysis dataset extracted from Cochrane Systematic Review CD013355. Contains study-level data with treatment comparisons, outcome measures, and effect sizes. Data includes both binary outcomes (events/totals) and continuous outcomes (mean/SD/N) where applicable.

\strong{Cochrane Review DOI:} \url{https://doi.org/10.1002/14651858.CD013355.pub2}

\strong{Dataset Statistics:}
\itemize{
  \item Number of studies: 21
  \item Outcome type: binary
  \item Columns: 27
}
}
\usage{data("CD013355_pub2_data")}
\format{
A data frame with 21 rows (studies) and 27 columns:
\describe{
  \item{Analysis.group}{Group/category of the meta-analysis}
  \item{Analysis.number}{Numeric identifier for the analysis within the review}
  \item{Analysis.name}{Name/description of the meta-analysis comparison}
  \item{Subgroup}{Subgroup classification (if applicable)}
  \item{Subgroup.number}{Numeric identifier for subgroup classification}
  \item{Applicability}{Assessment of applicability of the study to the review question}
  \item{Study}{Study identifier (typically author and year)}
  \item{Study.year}{Publication year of the study}
  \item{GIV.Mean}{Generic inverse variance pooled mean effect estimate}
  \item{GIV.SE}{Standard error of the generic inverse variance effect estimate}
  \item{Experimental.mean}{Mean value in experimental group (continuous outcomes)}
  \item{Experimental.SD}{Standard deviation in experimental group (continuous outcomes)}
  \item{Experimental.cases}{Number of events in experimental/treatment group (binary outcomes)}
  \item{Experimental.N}{Total participants in experimental/treatment group}
  \item{Control.mean}{Mean value in control group (continuous outcomes)}
  \item{Control.SD}{Standard deviation in control group (continuous outcomes)}
  \item{Control.cases}{Number of events in control group (binary outcomes)}
  \item{Control.N}{Total participants in control group}
  \item{O.E}{Observed minus expected events (O-E statistic)}
  \item{Variance}{Variance of the effect estimate}
  \item{Weight}{Weight assigned to the study in the meta-analysis}
  \item{Mean}{Pooled mean effect estimate}
  \item{CI.start}{Lower bound of 95\\\% confidence interval for effect estimate}
  \item{CI.end}{Upper bound of 95\\\% confidence interval for effect estimate}
  \item{Footnotes}{Additional notes or footnotes from the original Cochrane review}
  \item{review_url}{URL to the Cochrane review on the Cochrane Library}
  \item{review_doi}{Digital Object Identifier (DOI) for the Cochrane systematic review}
}
}
\source{
Systematically extracted from the Cochrane Library (\url{https://www.cochranelibrary.com/}) using official Cochrane data export functionality. Part of the Pairwise70 package containing 501 Cochrane pairwise meta-analysis datasets.
}
\examples{
# Load the dataset
data(CD013355_pub2_data)

# View dataset structure
str(CD013355_pub2_data)
head(CD013355_pub2_data)

# View review metadata
cat("Cochrane Review DOI:", CD013355_pub2_data$review_doi[1], "\n")

# Example: Run meta-analysis with metafor (binary outcomes)
# library(metafor)
# result <- rma(measure = "OR",
#               ai = Experimental.cases,
#               n1i = Experimental.N,
#               ci = Control.cases,
#               n2i = Control.N,
#               data = CD013355_pub2_data,
#               method = "REML")
# summary(result)
}
\keyword{datasets}

