% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynrVersion.R
\name{dynr.version}
\alias{dynr.version}
\title{Current Version String}
\usage{
dynr.version(verbose = TRUE)
}
\arguments{
\item{verbose}{If TRUE, print detailed information to the console (default)

This function returns a string with the current version number of dynr. Optionally (with verbose = TRUE (the default)), it prints a message containing the version of R and the platform.  The primary purpose of the function is for bug reporting.}
}
\value{
A (length-one) object of class 'package_version'
}
\description{
Current Version String
}
\examples{
dynr.version()
dynr.version(verbose=FALSE)
packageVersion("dynr")
}

