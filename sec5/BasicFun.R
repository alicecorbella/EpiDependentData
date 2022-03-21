# PACKAGES AND DIRECTORY --------------------------------------------------
rm(list=ls())
setwd("~/EpiDependentData/sec5/")
library(Rcpp)
library(microbenchmark)

sourceCpp("BasicCppFun.cpp")
# NEGATIVE BINOMIAL -------------------------------------------------------

# density of the NEGATIVE BINOMIA r, eta
dnbinomR <- function(x, r, eta, log=FALSE){
  lfx <- (lgamma(x+r)-(lgamma(x+1)+lgamma(r)))+x*log(1-1/eta)+r*log((1/eta))
  if(log==TRUE){
    f <- lfx
  }else{
    f <- exp(lfx)
  }
  return(f)
}

# in the r parametriztion p=1/eta and size=r

# parameter values
r   <- 300
eta <- 20
# simulations
xs <- rnbinom(10, size=r, prob = 1/eta)

dnbinom(x=xs, size=r, prob = 1/eta, log=TRUE)
dnbinomR(x=xs, r=r, eta=eta, log=TRUE)
dnbinom(x=xs, size=r, prob = 1/eta)
dnbinomR(x=xs, r=r, eta=eta)
round(dnbinom(x=xs, size=r, prob = 1/eta, log=TRUE), 10)==round(dnbinomR(x=xs, r=r, eta=eta, log=TRUE), 10)
round(dnbinom(x=xs, size=r, prob = 1/eta), 10)==round(dnbinomR(x=xs, r=r, eta=eta), 10)

# mine takes slghtly longer but not much, probably because of the gamma which is in C already
microbenchmark(dnbinom(x=xs, size=r, prob = 1/eta, log=TRUE),
               dnbinomR(x=xs, r=r, eta=eta, log=TRUE),
               dnbinom(x=xs, size=r, prob = 1/eta),
               dnbinomR(x=xs, r=r, eta=eta),
               dnbinomC(x=xs, r=r, eta=eta, lg=TRUE))


sum(xs*log(1-1/rep(20, 10))*(lgamma(xs+rep(300, 10))))
(a <-dnbinomC(x=xs, r=rep(r, 10), eta=rep(eta, 10), lg=TRUE))

(a <-dnbinomC(x=xs, r=r, eta=eta, lg=TRUE))

microbenchmark(sum(dnbinom(x=xs, size=r, prob = 1/eta, log=TRUE)),
               dnbinomC(x=xs, r=r, eta=eta, lg=TRUE))

microbenchmark(dnbinom(x=xs[1], size=6, prob = 1/(296), log=TRUE)+
    dnbinom(x=xs[2], size=7, prob = 1/(297), log=TRUE)+
    dnbinom(x=xs[3], size=8, prob = 1/(298), log=TRUE)+
    dnbinom(x=xs[4], size=9, prob = 1/(299), log=TRUE)+
    dnbinom(x=xs[5], size=10, prob = 1/(300), log=TRUE)+
    dnbinom(x=xs[6], size=11, prob = 1/(301), log=TRUE)+
    dnbinom(x=xs[7], size=12, prob = 1/(302), log=TRUE)+
    dnbinom(x=xs[8], size=13, prob = 1/(303), log=TRUE)+
    dnbinom(x=xs[9], size=14, prob = 1/(304), log=TRUE)+
    dnbinom(x=xs[10], size=15, prob = 1/(305), log=TRUE))
microbenchmark(sum(dnbinom(x=xs, size=6:15, prob = 1/(296:305), log=TRUE))
)
