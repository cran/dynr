% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynrRecipe.R
\name{dynr.ldl}
\alias{dynr.ldl}
\title{LDL Decomposition for Matrices}
\usage{
dynr.ldl(x)
}
\arguments{
\item{x}{a numeric matrix

This is a wrapper function around the \code{\link{chol}} function.
The goal is to factor a square, symmetric, positive (semi-)definite matrix into the product of a lower triangular matrix, a diagonal matrix, and the transpose of the lower triangular matrix.
The value returned is a lower triangular matrix with the elements of D on the diagonal.}
}
\value{
A matrix
}
\description{
LDL Decomposition for Matrices
}

