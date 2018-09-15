### R code from vignette source 'NonlinearContinuousTimeModels.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: setup
###################################################
# library(knitr)
# render_sweave()
# set.seed(22)
# knit_hooks$set(crop = hook_pdfcrop)
# opts_chunk$set(fig.path = 'figures/plots-', warning = FALSE, fig.align = 'center', width.cutoff = 80, fig.show = 'hold', eval = TRUE, echo = TRUE, message = FALSE, background = "white", prompt = FALSE, highlight = FALSE, comment = NA, tidy = FALSE, out.truncate = 80)
# options(replace.assign = TRUE, width = 80, prompt = "R> ", scipen = 12, digits = 3,crop=TRUE)

#setwd("~/Desktop/Repos/dynr/vignettes/") #set this working directory!
Sys.setenv(TEXINPUTS = getwd(),
  BIBINPUTS = getwd(),
  BSTINPUTS = getwd())


###################################################
### code chunk number 2: NonlinearContinuousTimeModels.Rnw:114-115 (eval = FALSE)
###################################################
## file.edit(system.file("demo", "NonlinearODE.R", package = "dynr"))


###################################################
### code chunk number 3: NonlinearContinuousTimeModels.Rnw:119-126 (eval = FALSE)
###################################################
## #------------------------------------------------------------------------------
## # Example 1: Nonlinear Continuous-time Models
## #------------------------------------------------------------------------------
## require(dynr)
## # ---- Read in the data ----
## data(PPsim)
## PPdata <- dynr.data(PPsim, id = "id", time = "time", observed = c("x", "y"))


###################################################
### code chunk number 4: NonlinearContinuousTimeModels.Rnw:146-153 (eval = FALSE)
###################################################
## # dynamics
## preyFormula <- prey ~ a * prey - b * prey * predator
## predFormula <- predator ~ - c * predator + d * prey * predator
## ppFormula <- list(preyFormula, predFormula)
## ppDynamics <- prep.formulaDynamics(formula = ppFormula,
##   startval = c(a = 2.1, c = 0.8, b = 1.9, d = 1.1),
##   isContinuousTime = TRUE)


###################################################
### code chunk number 5: NonlinearContinuousTimeModels.Rnw:160-161 (eval = FALSE)
###################################################
## file.edit(system.file("demo", "RSNonlinearDiscrete.R", package = "dynr"))


###################################################
### code chunk number 6: NonlinearContinuousTimeModels.Rnw:174-186 (eval = FALSE)
###################################################
## # Measurement (factor loadings)
## meas <- prep.measurement(
##   values.load = diag(1, 2),
##   obs.names = c('x', 'y'),
##   state.names = c('prey', 'predator'))
## 
## # alternatively, use prep.loadings
## meas <- prep.loadings(
##   map = list(
##     prey = "x",
##     predator = "y"),
##   params = NULL)


###################################################
### code chunk number 7: NonlinearContinuousTimeModels.Rnw:191-198 (eval = FALSE)
###################################################
## #measurement and dynamics covariances
## mdcov <- prep.noise(
##   values.latent = diag(0, 2),
##   params.latent = diag(c("fixed", "fixed"), 2),
##   values.observed = diag(rep(0.3, 2)),
##   params.observed = diag(c("var_1", "var_2"), 2)
## )


###################################################
### code chunk number 8: NonlinearContinuousTimeModels.Rnw:207-214 (eval = FALSE)
###################################################
## # Initial conditions on the latent state and covariance
## initial <- prep.initial(
##   values.inistate = c(3, 1),
##   params.inistate = c("fixed", "fixed"),
##   values.inicov = diag(c(0.01, 0.01)), 
##   params.inicov = diag("fixed", 2)
## )


###################################################
### code chunk number 9: NonlinearContinuousTimeModels.Rnw:219-224 (eval = FALSE)
###################################################
## #constraints
## trans <- prep.tfun(formula.trans = list(a ~ exp(a), b ~ exp(b), 
##                                         c ~ exp(c), d ~ exp(d)),
##   formula.inv = list(a ~ log(a), b ~ log(b), 
##                      c ~ log(c), d ~ log(d)))


###################################################
### code chunk number 10: NonlinearContinuousTimeModels.Rnw:234-255 (eval = FALSE)
###################################################
## #------------------------------------------------------------------------------
## # Cooking materials
## 
## # Put all the recipes together in a Model Specification
## model2.1 <- dynr.model(dynamics = ppDynamics, 
##   measurement = meas, noise = mdcov, 
##   initial = initial, transform = trans, 
##   data = PPdata,
##   outfile = "NonlinearODE.c")
## 
## # Check the model specification
## printex(model2.1, 
##   ParameterAs = model2.1$param.names,
##   show = FALSE, printInit = TRUE,
##   outFile = "NonlinearODE.tex")
## #tools::texi2pdf("NonlinearODE.tex")
## #system(paste(getOption("pdfviewer"), "NonlinearODE.pdf"))
## 
## # Estimate free parameters
## res2.1 <- dynr.cook(dynrModel = model2.1)
## 


###################################################
### code chunk number 11: NonlinearContinuousTimeModels.Rnw:262-276 (eval = FALSE)
###################################################
## # Examine results
## # True parameter values a = 2, b = 2, c = 1, d = 1
## summary(res2.1)
## #       names parameters       s.e.  t-value  ci.lower  ci.upper 
## # a         a  1.9637320 0.06946322 28.27010 1.8275866 2.0998774 
## # c         c  1.0023304 0.03062620 32.72788 0.9423042 1.0623567 
## # b         b  1.9327832 0.06216237 31.09250 1.8109472 2.0546192 
## # d         d  0.9608279 0.02628627 36.55246 0.9093078 1.0123481 
## # var_1 var_1  0.2399578 0.01089095 22.03277 0.2186119 0.2613036 
## # var_2 var_2  0.2380317 0.01072899 22.18585 0.2170033 0.2590601 
## # 
## # -2 log-likelihood value at convergence = 2843.19
## # AIC = 2855.19
## # BIC = 2884.64


###################################################
### code chunk number 12: NonlinearContinuousTimeModels.Rnw:282-289 (eval = FALSE)
###################################################
## #------------------------------------------------------------------------------
## # Example 2: Regime-Switching Nonlinear Continuous-time Model
## #------------------------------------------------------------------------------
## # ---- Read in the data ----
## data("RSPPsim")
## data <- dynr.data(RSPPsim, id = "id", time = "time",
##   observed = c("x", "y"), covariate = "cond")


###################################################
### code chunk number 13: NonlinearContinuousTimeModels.Rnw:293-296 (eval = FALSE)
###################################################
## cPreyFormula <- prey ~ a * prey - e * prey ^ 2 - b * prey * predator
## cPredFormula <- predator ~ f * predator - c * predator ^ 2 + d * prey * predator
## cpFormula <- list(cPreyFormula, cPredFormula)


###################################################
### code chunk number 14: NonlinearContinuousTimeModels.Rnw:300-304 (eval = FALSE)
###################################################
## rsFormula <- list(ppFormula, cpFormula)
## dynm <- prep.formulaDynamics(formula = rsFormula,
##    startval = c(a = 2.1, c = 3, b = 1.2, d = 1.2, e = 1, f = 2),
##    isContinuousTime = TRUE)


###################################################
### code chunk number 15: NonlinearContinuousTimeModels.Rnw:341-359 (eval = FALSE)
###################################################
## # Regime-switching function
## # The RS model assumes that each element of the transition probability 
## # matrix (TPM) can be expressed as a linear predictor (lp).
## # LPM = 
## # lp(p11) ~ 1 + x1 + x2 + ... + xn,   lp(p12) ~ 1 + x1 + x2 + ... + xn
## # lp(p21) ~ 1 + x1 + x2 + ... + xn,   lp(p22) ~ 1 + x1 + x2 + ... + xn
## # Here I am specifying lp(p12) and lp(p22); the remaining elements
## # lp(p11) and lp(p21) are fixed at zero.
## # nrow = numRegimes, ncol = numRegimes*(numCovariates+1)
## 
## regimes <- prep.regimes(
##   values = matrix(c(0, 0, -1, 1.5,
##                     0, 0, -1, 1.5),
##                 nrow = 2, ncol = 4, byrow = T), 
##   params = matrix(c("fixed", "fixed", "int_1", "slp_1",
##                     "fixed", "fixed", "int_2", "slp_2"), 
##                 nrow = 2, ncol = 4, byrow = T), 
##   covariates = "cond")


###################################################
### code chunk number 16: NonlinearContinuousTimeModels.Rnw:383-384 (eval = FALSE)
###################################################
## file.edit(system.file("demo", "RSNonlinearODE.R", package = "dynr"))


###################################################
### code chunk number 17: NonlinearContinuousTimeModels.Rnw:387-432 (eval = FALSE)
###################################################
## # Measurement (factor loadings)
## meas <- prep.measurement(
##   values.load = diag(1, 2),
##   obs.names = c('x', 'y'),
##   state.names = c('prey', 'predator'))
## 
## # Initial conditions on the latent state and covariance
## initial <- prep.initial(
##   values.inistate = c(3, 1),
##   params.inistate = c("fixed", "fixed"),
##   values.inicov = diag(c(0.01, 0.01)), 
##   params.inicov = diag("fixed", 2),
##   values.regimep = c(.8473, 0),
##   params.regimep = c("fixed", "fixed"))
## 
## #measurement and dynamics covariances
## mdcov <- prep.noise(
##   values.latent = diag(0, 2),
##   params.latent = diag(c("fixed","fixed"), 2),
##   values.observed = diag(rep(0.5,2)),
##   params.observed = diag(rep("var_epsilon",2),2)
## )
## 
## # dynamics
## preyFormula <- prey ~ a * prey - b * prey * predator
## predFormula <- predator ~ - c * predator + d * prey * predator
## ppFormula <- list(preyFormula, predFormula)
## cPreyFormula <- prey ~ a * prey - e * prey ^ 2 - b * prey * predator
## cPredFormula <- predator ~ 
##   f * predator - c * predator ^ 2 + d * prey * predator
## cpFormula <- list(cPreyFormula, cPredFormula)
## rsFormula <- list(ppFormula, cpFormula)
## 
## dynm <- prep.formulaDynamics(formula = rsFormula,
##   startval = c(a = 2.1, c = 3, b = 1.2, d = 1.2, e = 1, f = 2),
##   isContinuousTime = TRUE)
## 
## #constraints
## tformList <- list(a ~ exp(a), b ~ exp(b), c ~ exp(c),
##     d ~ exp(d), e ~ exp(e), f ~ exp(f))
## tformInvList <- list(a ~ log(a), b ~ log(b), c ~ log(c),
##     d ~ log(d), e ~ log(e), f ~ log(f))
## trans <- prep.tfun(
##     formula.trans = tformList,
##     formula.inv = tformInvList)


###################################################
### code chunk number 18: NonlinearContinuousTimeModels.Rnw:439-462 (eval = FALSE)
###################################################
## # Cooking materials
## 
## # Put all the recipes together in a Model Specification
## model2.2 <- dynr.model(dynamics = dynm, measurement = meas,
##   noise = mdcov, initial = initial,
##   regimes = regimes, transform = trans,
##   data = data,
##   outfile = "RSNonlinearODE_1.c")
## 
## # Check the model specification using LaTeX
## printex(model2.2, ParameterAs = names(model2.2), printInit = TRUE, printRS = TRUE,
##   outFile = "RSNonlinearODE_1.tex")
## #tools::texi2pdf("RSNonlinearODE_1.tex")
## #system(paste(getOption("pdfviewer"), "RSNonlinearODE_1.pdf"))
## 
## model2.2$ub[ c("int_1", "int_2", "slp_1", "slp_2") ] <- c(0, 0, 10, 10)
## model2.2$lb[ c("int_1", "int_2", "slp_1", "slp_2") ] <- c(-10, -10, 0, 0)
## 
## # Estimate free parameters
## res2.2 <- dynr.cook(model2.2)
## 
## # Examine results
## summary(res2.2)


###################################################
### code chunk number 19: NonlinearContinuousTimeModels.Rnw:466-470 (eval = FALSE)
###################################################
## dynr.ggplot(res2.2, model2.2, style = 2, 
##   names.regime = c("Summer", "Winter"),
##   title = "", idtoPlot = 9,
##   text = element_text(size = 16))


###################################################
### code chunk number 20: NonlinearContinuousTimeModels.Rnw:492-499 (eval = FALSE)
###################################################
## plotFormula(model2.2, ParameterAs = names(model2.2)) +
##   ggtitle("(A)") +
##   theme(plot.title = element_text(hjust = 0.5, vjust = 0.01, size = 16)) 
## 
## plotFormula(model2.2, ParameterAs = coef(res2.2)) +
##   ggtitle("(B)") +
##   theme(plot.title = element_text(hjust = 0.5, vjust = 0.01, size = 16))


