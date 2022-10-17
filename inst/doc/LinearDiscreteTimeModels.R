### R code from vignette source 'LinearDiscreteTimeModels.Rnw'

###################################################
### code chunk number 1: LinearDiscreteTimeModels.Rnw:51-52 (eval = FALSE)
###################################################
## file.edit(system.file("demo", "RSLinearDiscreteYang.R", package = "dynr"))


###################################################
### code chunk number 2: LinearDiscreteTimeModels.Rnw:77-83 (eval = FALSE)
###################################################
## #---- Load packages ----
## require("dynr")
## #---- Read in data and create dynr data object----
## data("EMG")
## EMGdata <- dynr.data(EMG, id = 'id', time = 'time', 
##   observed = 'iEMG', covariates = 'SelfReport')


###################################################
### code chunk number 3: LinearDiscreteTimeModels.Rnw:117-122 (eval = FALSE)
###################################################
## #---- Dynamic functions ----
## recDyn <- prep.matrixDynamics(
##   values.dyn = list(matrix(.1, 1, 1), matrix(.5, 1, 1)),
##   params.dyn = list(matrix('phi_1', 1, 1), matrix('phi_2', 1, 1)),
##   isContinuousTime = FALSE)


###################################################
### code chunk number 4: LinearDiscreteTimeModels.Rnw:135-145 (eval = FALSE)
###################################################
## #---- Measurement ----
## recMeas <- prep.measurement(
##   values.load = rep(list(matrix(1, 1, 1)), 2),
##   values.int = list(matrix(4, 1, 1), matrix(3, 1, 1)),
##   params.int = list(matrix('mu_1', 1, 1), matrix('mu_2', 1, 1)),
##   values.exo = list(matrix(0, 1, 1), matrix(1, 1, 1)),
##   params.exo = list(matrix('fixed', 1, 1), matrix('beta_2', 1, 1)),
##   obs.names = c('iEMG'),
##   state.names = c('eta'),
##   exo.names = c("SelfReport"))


###################################################
### code chunk number 5: LinearDiscreteTimeModels.Rnw:153-159 (eval = FALSE)
###################################################
## #---- Dynamic and measurement noise cov structures----
## recNoise <- prep.noise(
##   values.latent = matrix(1, 1, 1),
##   params.latent = matrix('dynNoise', 1, 1),
##   values.observed = matrix(0, 1, 1),
##   params.observed = matrix('fixed', 1, 1))


###################################################
### code chunk number 6: LinearDiscreteTimeModels.Rnw:176-184 (eval = FALSE)
###################################################
## #---- Initial condition specification ----
## recIni <- prep.initial(
##   values.inistate = matrix(0, 1, 1),
##   params.inistate = matrix('fixed', 1, 1),
##   values.inicov = matrix(1, 1, 1),
##   params.inicov = matrix('fixed', 1, 1),
##   values.regimep = c(1, 0),
##   params.regimep = c('fixed', 'fixed'))


###################################################
### code chunk number 7: LinearDiscreteTimeModels.Rnw:209-213 (eval = FALSE)
###################################################
## # ---- Regimes-switching model ----
## recReg <- prep.regimes(
##   values = matrix(c(.7, -1, 0, 0), 2, 2),
##   params = matrix(c('c11', 'c21', 'fixed', 'fixed'), 2, 2))


###################################################
### code chunk number 8: LinearDiscreteTimeModels.Rnw:262-266 (eval = FALSE)
###################################################
## recReg2 <- prep.regimes(
##   values = matrix(c(.8, -1, 0, 0), 2, 2),
##   params = matrix(c('c_Delta11', 'c1', 'fixed', 'fixed'), 2, 2),
##   deviation = TRUE, refRow = 2)


###################################################
### code chunk number 9: LinearDiscreteTimeModels.Rnw:274-287 (eval = FALSE)
###################################################
## #---- Create model  ----
## 
## rsmod <- dynr.model(
##   dynamics = recDyn,
##   measurement = recMeas,
##   noise = recNoise,
##   initial = recIni,
##   regimes = recReg,
##   data = EMGdata)
## 
## #---- Create model and cook it all up  ----
## 
## yum <- dynr.cook(rsmod)


###################################################
### code chunk number 10: LinearDiscreteTimeModels.Rnw:293-295 (eval = FALSE)
###################################################
## #---- Serve it! ----
## summary(yum)


###################################################
### code chunk number 11: LinearDiscreteTimeModels.Rnw:306-307 (eval = FALSE)
###################################################
## plot(yum, dynrModel = rsmod, style = 1, textsize = 5)


###################################################
### code chunk number 12: LinearDiscreteTimeModels.Rnw:310-325 (eval = FALSE)
###################################################
## #pdf('./Figures/plotRSGG.pdf', height=7, width=12)
## dynr.ggplot(yum, dynrModel = rsmod, style = 1,
##   names.regime = c("Deactivated", "Activated"),
##   title = "(B) Results from RS-AR model", numSubjDemo = 1,
##   shape.values = c(1),
##   text = element_text(size = 24),
##   is.bw = TRUE)
## #dev.off()
## 
## autoplot(yum, dynrModel = rsmod, style = 1,
##   names.regime = c("Deactivated", "Activated"),
##   title = "(B) Results from RS-AR model", numSubjDemo = 1,
##   shape.values = c(1),
##   text = element_text(size = 16),
##   is.bw = TRUE)


###################################################
### code chunk number 13: LinearDiscreteTimeModels.Rnw:330-339 (eval = FALSE)
###################################################
## plotFormula(dynrModel = rsmod, ParameterAs = names(rsmod),
##   printDyn = TRUE, printMeas = TRUE) +
##   ggtitle("(A)") +
##   theme(plot.title = element_text(hjust = 0.5, vjust = 0.01, size = 16))
## 
## plotFormula(dynrModel = rsmod, ParameterAs = coef(yum),
##   printDyn = TRUE, printMeas = TRUE) +
##   ggtitle("(B)") +
##   theme(plot.title = element_text(hjust = 0.5, vjust = 0.01, size = 16))


###################################################
### code chunk number 14: LinearDiscreteTimeModels.Rnw:353-357 (eval = FALSE)
###################################################
## printex(rsmod,
##   ParameterAs = names(rsmod),
##   printInit = TRUE, printRS = TRUE,
##   outFile = "RSLinearDiscreteYang.tex")


###################################################
### code chunk number 15: LinearDiscreteTimeModels.Rnw:360-362 (eval = FALSE)
###################################################
## tools::texi2pdf("RSLinearDiscreteYang.tex")
## system(paste(getOption("pdfviewer"), "RSLinearDiscreteYang.pdf"))


