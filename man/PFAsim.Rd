% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataDoc.R
\docType{data}
\name{PFAsim}
\alias{PFAsim}
\title{Simulated time series data of a multisubject process factor analysis}
\format{A data frame with 2,500 rows and 10 variables}
\usage{
data(PFAsim)
}
\description{
A multiple subject dataset simulated using a two factor process factor analysis model in discrete time
with 6 observed indicators for identifying two latent factors.
The variables are as follows:
}
\details{
\itemize{
  \item ID. Person ID variable (1 to 50 because there are 50 simulated people)
  \item Time. Time ID variable (1 to 50 because there are 50 time points)
  \item V1. Noisy observed variable 1
  \item V2. Noisy observed variable 2
  \item V3. Noisy observed variable 3
  \item V4. Noisy observed variable 4
  \item V5. Noisy observed variable 5
  \item V6. Noisy observed variable 6
  \item F1. True latent variable 1 scores
  \item F2. True latent variable 2 scores
}

Variables V1, V2, and V3 load on F1, whereas variables V4, V5, V6 load on F2.  The true values of the factor loadings are 1, 2, 1, 1, 2, and 1, respectively.  The true measurement error variance is 0.5 for all variables.  The true dynamic noise covariance has F1 with a variance of 2.77, F2 with a variance of 8.40, and their covariance is 2.47.  The across-time dynamics have autoregressive effects of 0.5 for both F1 and F2 with a cross-lagged effect from F1 to F2 at 0.4.  The cross-lagged effect from F2 to F1 is zero.  The true initial latent state distribution as mean zero and a diagonal covariance matrix with var(F1) = 2 and var(F2) = 1.  The generating model is the same for all individuals.
}
\examples{
# The following was used to generate the data
\dontrun{
set.seed(12345678)
library(mvtnorm)
# setting up matrices
time      <- 50
# Occasions to throw out to wash away the effects of initial condition
npad      <- 0
np        <- 50
ne        <- 2 #Number of latent variables
ny        <- 6 #Number of manifest variables
# Residual variance-covariance matrix
psi       <- matrix(c(2.77, 2.47,
                      2.47, 8.40),
                    ncol = ne, byrow = T)
# Lambda matrix containing contemporaneous relations among
# observed variables and 2 latent variables. 
lambda    <- matrix(c(1, 0,
                      2, 0,
                      1, 0,
                      0, 1,
                      0, 2,
                      0, 1),
                    ncol = ne, byrow = TRUE)
# Measurement error variances
theta     <- diag(.5, ncol = ny, nrow = ny)
# Lagged directed relations among variables
beta      <- matrix(c(0.5, 0,
                      0.4, 0.5), 
                    ncol = ne, byrow = TRUE)
a0        <- mvtnorm::rmvnorm(1, mean = c(0, 0),
                                 sigma = matrix(c(2,0,0,1),ncol=ne))
yall <- matrix(0,nrow = time*np, ncol = ny)
eall <- matrix(0,nrow = time*np, ncol = ne)
for (p in 1:np){
  # Latent variable residuals
  zeta      <- mvtnorm::rmvnorm(time+npad, mean = c(0, 0), sigma = psi)
  # Measurement errors
  epsilon   <- rmvnorm(time, mean = c(0, 0, 0, 0, 0, 0), sigma = theta)
  # Set up matrix for contemporaneous variables
  etaC      <- matrix(0, nrow = ne, ncol = time + npad)
  # Set up matrix for lagged variables
  etaL      <- matrix(0, nrow = ne, ncol = time + npad + 1)
  
  etaL[,1]  <- a0
  etaC[,1] <- a0
  # generate factors
  for (i in 2:(time+npad)){
    etaL[ ,i] <- etaC[,i-1]
    etaC[ ,i]   <- beta \%*\% etaL[ ,i] + zeta[i, ]
  }
  etaC <- etaC[,(npad+1):(npad+time)]
  eta <- t(etaC)
  
  # generate observed series
  y   <- matrix(0, nrow = time, ncol = ny)
  for (i in 1:nrow(y)){
    y[i, ] <- lambda \%*\% eta[i, ] + epsilon[i, ]
  }
  yall[(1+(p-1)*time):(p*time),] <- y
  eall[(1+(p-1)*time):(p*time),] <- eta
}
yall <- cbind(rep(1:np,each=time),rep(1:time,np),yall)
yeall <- cbind(yall,eall)
write.table(yeall,'PFAsim.txt',row.names=FALSE,
  col.names=c("ID", "Time", paste0("V", 1:ny), paste0("F", 1:ne)))
}
}
\keyword{datasets}

