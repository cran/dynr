## ----data, results="hide"-----------------------------------------------------
require(dynr)
# Data
data(Oscillator)
data <- dynr.data(Oscillator, id="id", time="times", observed="y1")

## ----measurement, results="hide"----------------------------------------------
meas <- prep.measurement(
	values.load=matrix(c(1, 0), 1, 2), 
	params.load=matrix(c('fixed', 'fixed'), 1, 2),
	state.names=c("Position","Velocity"),
	obs.names=c("y1"))

## ----noise cov, results="hide"------------------------------------------------
ecov <- prep.noise(
	values.latent=diag(c(0, 1), 2), params.latent=diag(c('fixed', 'dnoise'), 2), 
	values.observed=diag(1.5, 1), params.observed=diag('mnoise', 1)) 
dynamics <- prep.matrixDynamics(
	values.dyn=matrix(c(0, -0.1, 1, -0.2), 2, 2),
	params.dyn=matrix(c('fixed', 'spring', 'fixed', 'friction'), 2, 2), 
	isContinuousTime=TRUE)

## ----initials, results="hide"-------------------------------------------------
initial <- prep.initial(
	values.inistate=c(0, 1),
	params.inistate=c('inipos', 'fixed'), 
	values.inicov=diag(1, 2),
	params.inicov=diag('fixed', 2)) 


## ----model, results="hide"----------------------------------------------------

model <- dynr.model(dynamics=dynamics, measurement=meas, noise=ecov, initial=initial, data=data, outfile="LinearSDE.c")

## ----tex, results="hide",eval=FALSE-------------------------------------------
#  printex(model,ParameterAs=model$param.names,show=FALSE,printInit=TRUE,
#          outFile="LinearSDE.tex")
#  tools::texi2pdf("LinearSDE.tex")
#  system(paste(getOption("pdfviewer"), "LinearSDE.pdf"))

## ----cook, results="hide"-----------------------------------------------------
res <- dynr.cook(model, verbose=FALSE)

## ----serve--------------------------------------------------------------------
summary(res)

