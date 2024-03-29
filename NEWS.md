# dynr 0.1.16-XX
* 2021-04-12
* MAJOR BUG FIX with smoothed latent state covariance 

# dynr 0.1.15-1
* 2019-10-04
* Multiple imputation for missing data with dynr.mi() function
* New demo for multiple imputation called MILinearDiscrete
* New demo for time-varying parameters called SDETVP
* New wrapper functions for computing smoothed derivative estimates using penalized B-splines implemented with the fda package
* New plotting functions, including functions to generate diagnostic plots from smoothed derivative estimates, and plot phase portraits
* New demo for computing and visualizing smoothed derivative estimates in GetDerivs
* New functionality in dynr.cook() to estimate continuous-time dynamic models with mixed effects through use of theta.formula
* MAJOR BUG FIXES with missing data

# dynr 0.1.14-9
* 2019-04-01
* Outlier detection with dynr.taste() function
* Oulier removal and re-fitting with dynr.tast2() function
* New demo for outliers called OutlierDetection
* We now allow 1-regime recipe parts to co-exist with n-regime parts
* Lots of error checking was added around matching the number of regimes
* Many cases of the doDykstra error are now safely caught
* Generally cleaned up the error handling on models that failed to converge
* Shorten several demos to run faster

# dynr 0.1.13-4
* 2018-09-21
* You don't need to install R on Windows to C:/R anymore! The default (C:/Program Files/R) now works.
* New demo for time-varying parameter (TVP) models
* Several new vignettes covering a range of topics
* The deviation form of regimes now displays properly in plotFormula
* No longer require 'outfile' specification in dynr.model()
* Fix a pointer addressing issue that could have caused crashes

# dynr 0.1.12-5
* 2018-02-08
* New 'verbose' argument to dynr.cook turns off printing of optimization history
* New demo for Process Factor Analysis (PFA)
* Regime-switching printing in plotFormula() with new 'printRS' argument
* Greatly improved convergence rates for all models
* Allow full initial covariance estimation
* Fixed major bug in regime-switching matrix dynamics that formerly crashed R

# dynr 0.1.11-8
* 2017-08-21
* Noise printing by plotFormula() function
* Fixed innovation vector computation for larger than 1-dimensional observations

# dynr 0.1.11-2
* 2017-06-16
* New demo for a linear oscillator with time-varying parameters
* Fixed printex output for covariates and deviation form of the initial conditions
* Fixed memory leak for intercepts in measurement models

# dynr 0.1.10
* 2017-05-19
* Use of individual-level covariates in the initial conditions. See ?prep.initial for details.
* Deviation form of regime-switching models.  Seee ?prep.regimes for details.
* Access to the predicted, filtered, and smoothed latent variable estimates, and other by-products from the regime-switching extended Kalman filter in the 'cooked' model.
* We now allow calculation of the negative log-likelihood value, the hessian matrix, and the predicted, filtered, and smoothed latent variable estimates at fixed parameter values without parameter optimization.
* Beta version of a multiple imputation procedure.  See ?dynr.mi for details.
* Fixed a rounding bug that improves free parameter optimization, especially for models with many observed variables.
* Improved documentation throughout
* Added more examples in the help pages

# dynr 0.1.9
* 2017-02-21
* A new demo example is added to replicate the results from Yang & Chow (2010) paper.
* Some standard S3 methods are added for the dynrCook class object.
* autoplot() is added as an alias for dynr.ggplot().
* dynr.data() now automatically handles ts class objects and equally spaced data with missingness.
* Changes are made to accommodate the new release of ggplot2. 

# dynr 0.1.8
* 2016-08-12
* In single-regime models, free parameters for intercepts and covariate effects
    in the measurement model can now be properly estimated.
* Standard errors are more frequently returned
* Flags indicate problematic standard errors.
* Warning messages are more helpful regarding standard errors.
* A weight flag allows easier convergence of multi-subject models.
* Several new plotting features.

# dynr 0.1.7
* 2016-06-07
* Initial release to CRAN!