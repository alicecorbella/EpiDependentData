
# PACKAGES ----------------------------------------------------------------


if(!("mvtnorm" %in% (.packages()))){
  library(mvtnorm)
}

if(!("Rcpp" %in% (.packages()))){
  library(Rcpp)
}


# BASIC OPERATIONS --------------------------------------------------------


# file containing some of the operations that we do not do in RCPP 
# but in R, mostly copied form Functions.R

logit <- function(x){
  log(x/(1-x))
}

toNorm.mulpi <- function(x){
  u <- rep(NA, length(x))
  for (i in 1:length(x)){
    u[i] <- log(x[i]/(1-sum(x[1:i])))
  }
  return(u)
}

expit <- function(x){
  exp(x)/(1+exp(x))
}

toTh.mulpi <- function(u){
  x <- rep(NA, length(u))
  for (i in 1:length(u)){
    x[i] <- exp(u[i])/prod(exp(u[1:i])+1)
  }
  return(x)
}

toTh <- function(x){
  Th <- c(exp(x[1]), toTh.mulpi(x[2:3]), 
          expit(x[4]), expit(x[5]))
  return(Th)
}

toNorm <- function(Th){
  Norm <- c(log(Th[1]),   toNorm.mulpi(Th[2:3]), 
            logit(Th[4]), logit(Th[5]))
  return(Norm)
}


