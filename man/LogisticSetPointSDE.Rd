% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataDoc.R
\docType{data}
\name{LogisticSetPointSDE}
\alias{LogisticSetPointSDE}
\title{Simulated time series data for a stochastic linear damped oscillator model with logistic time-varying setpoints}
\format{A data frame with 2410 rows and 6 variables}
\usage{
data(LogisticSetPointSDE)
}
\description{
A dataset simulated using a continuous-time stochastic linear damped oscillator model.
The variables are as follows:
}
\details{
\itemize{
  \item id. ID of the systems (1 to 10)
  \item times. Time index (241 time points for each system)
  \item x. Latent level variable
  \item y. Latent first derivative variable
  \item z. True values of time-varying setpoints
  \item obsy. Observed level
}
}
\examples{
# The following was used to generate the data
#---------------------------------------
\dontrun{
require(Sim.DiffProc)
freq <- -1
damp <- -.1
mu <- -2
r <- .5
b <- .1
sigma1 <- 0.1
sigma2 <- 0.1
fx <- expression(y, freq*(x-z) + damp*y, r*z*(1-b*z))
gx <- expression(0, sigma1, 0)
r3dall <- c()
for (j in 1:10){
  r3dtemp <- c(-5,0,.1)
  r3d <- r3dtemp
  for (i in seq(0.125, 30, by=0.125)){
    mod3dtemp <- snssde3d(drift=fx, diffusion=gx, M=1, t0=i-0.125,
        x0=as.numeric(r3dtemp), T=i, N=500, type="str",
        method="smilstein")
    r3dtemp <- rsde3d(mod3dtemp,at=i)
    r3d <-rbind(r3d,r3dtemp)
  }
  r3dall <- rbind(r3dall, cbind(r3d, id=j))
}

r3dall$obsy <- r3dall$x+rnorm(length(r3dall$x),0,1)
write.table(r3dall, file="LogisticSetPointSDE.txt")
}
}
\keyword{datasets}

