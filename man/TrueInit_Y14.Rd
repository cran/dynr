% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataDoc.R
\docType{data}
\name{TrueInit_Y14}
\alias{TrueInit_Y14}
\title{Simulated multilevel multi-subject time series of a Van der Pol Oscillator}
\format{A data frame with 60,000 rows and 10 variables}
\usage{
data(TrueInit_Y14)
}
\description{
A dataset simulated using methods described in the reference below.
 
 Reference:
 Chow, S., Lu, Z., Sherwood, A., and Zhu, H. (2016). Fitting Nonlinear Ordinary
 Differential Equation Models with Random Effects and Unknown Initial Conditions
 Using the Stochastic Approximation Expectation-Maximization (SAEM) Algorithm.
 Psychometrika, 81(1), 102-134.
}
\details{
The variables are as follows:

\itemize{
  \item batch. Batch number from simulation
  \item kk. Unclear
  \item trueInit. True initial condition
  \item id. Person ID
  \item time. Continuous time of measurement
  \item y1. Observed score 1
  \item y2. Observed score 2
  \item y3. Observed score 3
  \item co1. Covariate 1
  \item co2. Covariate 2
}
}
\keyword{datasets}

