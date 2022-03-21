# libray and directory ----------------------------------------------------
rm(list=ls())
library(mvtnorm)
library(Rcpp)

maind <- "~/EpiDependentData/sec4"
setwd(maind)

# SOURCE and LOAD ---------------------------------------------------------

# scenrio number
task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
s              <- as.numeric(task_id_string)

load(paste(maind,"/input/SMALLsim",s,".Robj", sep=""))

DATA       <- savelist[[1]]
HIDDEN     <- savelist[[2]]
THETACONST <- savelist[[3]]
THETAvec   <- savelist[[4]]

rm(savelist)

sourceCpp(paste(maind, "/FUNCTIONS/cppfun.cpp", sep = ""))
source(paste(maind, "/FUNCTIONS/RfunALL.R", sep = ""))
source(paste(maind, "/FUNCTIONS/RfunSB.R", sep = ""))


# STARTING POINTS AND CHAIN LENGTH -----------------------------------------
start1 <- toNorm(c(THETAvec, THETACONST$thEtoH, THETACONST$thHtoIC))[1:3]
start2 <- toNorm(c(THETAvec, THETACONST$thEtoH, THETACONST$thHtoIC))[1:3]+c(0.05, 0.05, 1)/2
start3 <- toNorm(c(THETAvec, THETACONST$thEtoH, THETACONST$thHtoIC))[1:3]-c(0.05, 0.05, 1)/2

npset <- 2000

DPT.thin  <- 10
DPT.nsmpl <- 1000
SMP.thin  <- 50
SMP.nsmpl <- 1100
SMP.bin   <- 100


DPT.thinB  <- 20
DPT.nsmplB <- 1000
SMP.thinB  <- 100
SMP.nsmplB <- 1100
SMP.binB   <- 100

# RUN CHAINS --------------------------------------------------------------

dpt   <- adptMHsingle(start=start1, sdini=c(0.05, 0.05, 1), nsmpl=DPT.nsmpl, thin=DPT.thin, 
                      np=npset, stp.dpt=500, 
                      thetaconst=THETACONST,MCWM=T) 
Sys.time()
save.image(file=paste(maind,"/outputDET/SMALLsim",s,".Robj", sep=""))

ch1 <- runMHsingle (start=start1, sd=dpt[[1]], nsmpl=SMP.nsmpl, thin=SMP.thin, np=npset)
save.image(file=paste(maind,"/outputDET/SMALLsim",s,".Robj", sep=""))
Sys.time()
ch2 <- runMHsingle (start=start2, sd=dpt[[1]], nsmpl=SMP.nsmpl, thin=SMP.thin, np=npset)
save.image(file=paste(maind,"/outputDET/SMALLsim",s,".Robj", sep=""))
Sys.time()
ch3 <- runMHsingle (start=start3, sd=dpt[[1]], nsmpl=SMP.nsmpl, thin=SMP.thin, np=npset)
save.image(file=paste(maind,"/outputDET/SMALLsim",s,".Robj", sep=""))
Sys.time()

vc <- cov(rbind(ch1[[2]][(SMP.bin+1):SMP.nsmpl, ],
                ch2[[2]][(SMP.bin+1):SMP.nsmpl, ],
                ch3[[2]][(SMP.bin+1):SMP.nsmpl, ] ))

dptBLK   <- adptMHblocked(start=start1, vcini=vc, nsmpl=DPT.nsmplB, thin=DPT.thinB, 
                          np=npset, stp.dpt=100, 
                          thetaconst=THETACONST,MCWM=T)
save.image(file=paste(maind,"/outputDET/SMALLsim",s,".Robj", sep=""))
Sys.time()
ch1BLK <- runMHblocked(start=start1, vcscaled=dptBLK[[1]]*vc, nsmpl=SMP.nsmplB, thin=SMP.thinB, np=npset)
save.image(file=paste(maind,"/outputDET/SMALLsim",s,".Robj", sep=""))
Sys.time()
ch2BLK <- runMHblocked(start=start2, vcscaled=dptBLK[[1]]*vc, nsmpl=SMP.nsmplB, thin=SMP.thinB, np=npset)
save.image(file=paste(maind,"/outputDET/SMALLsim",s,".Robj", sep=""))
Sys.time()
ch3BLK <- runMHblocked(start=start3, vcscaled=dptBLK[[1]]*vc, nsmpl=SMP.nsmplB, thin=SMP.thinB, np=npset)
save.image(file=paste(maind,"/outputDET/SMALLsim",s,".Robj", sep=""))
Sys.time()
