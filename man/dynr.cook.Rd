% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynrCook.R
\name{dynr.cook}
\alias{dynr.cook}
\title{Cook a dynr model to estimate its free parameters}
\usage{
dynr.cook(dynrModel, conf.level = 0.95, infile, optimization_flag = TRUE,
  hessian_flag = TRUE, verbose = TRUE, weight_flag = FALSE,
  debug_flag = FALSE, perturb_flag = FALSE)
}
\arguments{
\item{dynrModel}{a dynr model compiled using dynr.model, consisting of recipes for submodels, 
starting values, parameter names, and C code for each submodel}

\item{conf.level}{a cumulative proportion indicating the level of desired confidence intervals for
the final parameter estimates (default is .95)}

\item{infile}{(not required for models specified through the recipe functions) the name of a file 
that has the C codes for all dynr submodels for those interested in specifying a model directly in C}

\item{optimization_flag}{a flag (TRUE/FALSE) indicating whether optimization is to be done.}

\item{hessian_flag}{a flag (TRUE/FALSE) indicating whether the Hessian matrix is to be calculated.}

\item{verbose}{a flag (TRUE/FALSE) indicating whether more detailed intermediate output during the 
estimation process should be printed}

\item{weight_flag}{a flag (TRUE/FALSE) indicating whether the negative log likelihood function should 
be weighted by the length of the time series for each individual}

\item{debug_flag}{a flag (TRUE/FALSE) indicating whether users want additional dynr output that can 
be used for diagnostic purposes}

\item{perturb_flag}{a flag (TRUE/FLASE) indicating whether to perturb the latent states during estimation. Only useful for ensemble forecasting.}
}
\value{
Object of class dynrCook.
}
\description{
Cook a dynr model to estimate its free parameters
}
\details{
Free parameter estimation uses the SLSQP routine from NLOPT.

The typical items returned in the cooked model are the filtered and smoothed latent variable estimates. 
\code{eta_smooth_final}, \code{error_cov_smooth_final} and \code{pr_t_given_T} are respectively 
time-varying smoothed latent variable mean estimates, smoothed error covariance estimates, 
and smoothed regime probability. 
\code{eta_filtered}, \code{error_cov_filtered} and \code{pr_t_given_t} are respectively 
time-varying filtered latent variable mean estimates, filtered error covariance matrix estimates, 
and filtered regime probability.
Note that if \code{theta.formula} is provided in \code{dynrModel@dynamics}, this assumes  that random effects are present in the dynamic equation. This would call an internal function to insert the random effect components as additional state variables. In this case, the last set of elements (rows) in \code{eta_smooth_final} would contain the estimated random effect components.

When \code{debug_flag} is TRUE, then additional information is passed into the cooked model. 
\code{eta_predicted}, \code{error_cov_predicted}, \code{innov_vec}, and \code{residual_cov} are respectively 
time-varying predicted latent variable mean estimates, predicted error covariance matrix estimates, the error/residual estimates (innovation vector),
and the error/residual covariance matrix estimates.

The exit flag given after optimization has finished is from the SLSQP optimizer.  Generally, error codes have negative values and successful codes have positive values.  However, codes 5 and 6 do not indicate the model converged, but rather simply ran out of iterations or time, respectively.  A more full description of each code is available at \url{https://nlopt.readthedocs.io/en/latest/NLopt_Reference/#return-values} and is also listed in the table below.

\tabular{lcl}{
NLOPT Term \tab Numeric Code \tab Description \cr
SUCCESS \tab 1 \tab Generic success return value. \cr
STOPVAL_REACHED \tab 2 \tab Optimization stopped because stopval (above) was reached. \cr
FTOL_REACHED \tab 3 \tab Optimization stopped because ftol_rel or ftol_abs (above) was reached. \cr
XTOL_REACHED \tab 4 \tab Optimization stopped because xtol_rel or xtol_abs (above) was reached. \cr
MAXEVAL_REACHED \tab 5 \tab Optimization stopped because maxeval (above) was reached. \cr
MAXTIME_REACHED \tab 6 \tab Optimization stopped because maxtime (above) was reached. \cr
FAILURE \tab -1 \tab Generic failure code. \cr
INVALID_ARGS \tab -2 \tab Invalid arguments (e.g. lower bounds are bigger than upper bounds, an unknown algorithm was specified, etcetera). \cr
OUT_OF_MEMORY \tab -3 \tab Ran out of memory. \cr
ROUNDOFF_LIMITED \tab -4 \tab Halted because roundoff errors limited progress. (In this case, the optimization still typically returns a useful result.) \cr
FORCED_STOP \tab -5 \tab Halted because of a forced termination: the user called nlopt_force_stop(opt) on the optimization's nlopt_opt object opt from the user's objective function or constraints. \cr
NONFINITE_FIT \tab -6 \tab Fit function is not finite (i.e., is NA, NaN, Inf or -Inf). \cr
}
The last row of this table corresponding to an exit code of -6, is not from NLOPT, but rather is specific to the \code{dynr} package.
}
\examples{
# Minimal model
require(dynr)

meas <- prep.measurement(
	values.load=matrix(c(1, 0), 1, 2),
	params.load=matrix(c('fixed', 'fixed'), 1, 2),
	state.names=c("Position","Velocity"),
	obs.names=c("y1"))

ecov <- prep.noise(
	values.latent=diag(c(0, 1), 2),
	params.latent=diag(c('fixed', 'dnoise'), 2),
	values.observed=diag(1.5, 1),
	params.observed=diag('mnoise', 1))

initial <- prep.initial(
	values.inistate=c(0, 1),
	params.inistate=c('inipos', 'fixed'),
	values.inicov=diag(1, 2),
	params.inicov=diag('fixed', 2))

dynamics <- prep.matrixDynamics(
	values.dyn=matrix(c(0, -0.1, 1, -0.2), 2, 2),
	params.dyn=matrix(c('fixed', 'spring', 'fixed', 'friction'), 2, 2),
	isContinuousTime=TRUE)

data(Oscillator)
data <- dynr.data(Oscillator, id="id", time="times", observed="y1")

model <- dynr.model(dynamics=dynamics, measurement=meas,
	noise=ecov, initial=initial, data=data)

\dontrun{
# Now cook the model!
cook <- dynr.cook(model,
	verbose=FALSE, optimization_flag=FALSE, hessian_flag=FALSE)
}
}
\seealso{
\code{\link{autoplot}}, \code{\link{coef}}, \code{\link{confint}},
\code{\link{deviance}}, \code{\link{initialize}}, \code{\link{logLik}},
\code{\link{names}}, \code{\link{nobs}}, \code{\link{plot}}, \code{\link{print}},
\code{\link{show}}, \code{\link{summary}}, \code{\link{vcov}}.
}

