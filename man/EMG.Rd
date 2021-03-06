% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataDoc.R
\docType{data}
\name{EMG}
\alias{EMG}
\title{Single-subject time series of facial electromyography data}
\format{A data frame with 695 rows and 4 variables}
\usage{
data(EMG)
}
\description{
A dataset obtained and analyzed in Yang and Chow (2010).
}
\details{
Reference:
 Yang, M-S. & Chow, S-M. (2010). Using state-space models with regime switching to
 represent the dynamics of facial electromyography (EMG) data. Psychometrika, 74(4), 744-771

The variables are as follows:

\itemize{
  \item id. ID of the participant (= 1 in this case, over 695 time points)
  \item time Time in seconds
  \item iEMG. Observed integrated facial electromyograhy data
  \item SelfReport. Covariate - the individual's concurrent self-reports
}
}
\keyword{datasets}

