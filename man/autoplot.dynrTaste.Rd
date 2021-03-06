% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynrPlot.R
\name{autoplot.dynrTaste}
\alias{autoplot.dynrTaste}
\title{The ggplot of the outliers estimates.}
\usage{
\method{autoplot}{dynrTaste}(object, numSubjDemo = 2, idtoPlot = NULL,
  names.state = NULL, names.observed = NULL, ...)
}
\arguments{
\item{object}{A dynrTaste object.}

\item{numSubjDemo}{The number of subjects, who have 
largest joint chi-square statistic, to be selected  for plotting.}

\item{idtoPlot}{Values of the ID variable to plot.}

\item{names.state}{(optional) The names of the states to be plotted, which should be a subset of the state.names slot of the measurement slot of dynrModel. If NULL, the t statistic plots for all state variables will be included.}

\item{names.observed}{(optional) The names of the observed variables to be plotted, which should be a subset of the obs.names slot of the measurement slot of dynrModel. If NULL, the t statistic plots for all observed variables will be included.}

\item{...}{Place holder for other arguments. Please do not use.}
}
\value{
a list of ggplot objects for each ID. 
The plots of chi-square statistics (joint and independent),
and the plots of t statistic for \code{names.state} and \code{names.observed} will be included.
Users can modify the ggplot objects using ggplot grammar.
If a \code{filename} is provided, a pdf of plots will be saved additionally.
}
\description{
The ggplot of the outliers estimates.
}

