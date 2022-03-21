rm(list=ls())
library(mvtnorm)
# library(microbenchmark)
library(TailRank)
library(Rcpp)
library(fields)
library(readxl)
library(xlsx)
library(inline)

thisdir <- "~/EpiDependentData/sec5/"
datadir <- "~/SensibleData/"
setwd(thisdir)



# DATA    --------------------------------------------------------
DATA             <- list()
thetaconst       <- list()
FOR_MODELLING    <- read_excel(paste(datadir, "USISS/USISS data for modelling.xlsx", sep = ""),sheet=1, col_names = TRUE)
data_ICU_1718    <- as.matrix((FOR_MODELLING[c(2:14, 16:35), c(4, 3)]))
data_ICU_1718    <- apply(data_ICU_1718, MARGIN =c(1,2), as.numeric)
DATA$yIC         <- as.vector(data_ICU_1718[,1])
thetaconst$detIC <- as.vector(data_ICU_1718[,2])/55268100
rm(data_ICU_1718, FOR_MODELLING)
FOR_MODELLING    <- read_excel(paste(datadir, "USISS/USISS data for modelling.xlsx", sep = ""),sheet=2, col_names = TRUE)
data_H_1718      <- as.matrix((FOR_MODELLING[c(2:14, 16:35), c(4, 3)]))
data_H_1718      <- apply(data_H_1718, MARGIN =c(1,2), as.numeric)
DATA$yH          <- as.vector(data_H_1718[,1])
thetaconst$detH  <- as.vector(data_H_1718[,2])/55268100
rm(data_H_1718, FOR_MODELLING)
# GP
rawdata <- read_excel(paste(datadir, "GP ILI/DailyGPIH_ILIByPHECByAgeGrpWinter1718.xlsx", sep = ""),
                      sheet="Data", col_names = TRUE)
dates <- (unique(rawdata$RepDate))
period <- seq(as.Date(min(dates)),as.Date(max(dates)), by=1)
newrawdata <- matrix(NA, ncol=16, nrow=length(period))
newrawdata[,1] <- period
newrawdata[,2] <- rep("allreg", length(period))
for (p in 1:length(period)){
  newrawdata[p,3:16 ] <- apply(rawdata[which(as.Date(rawdata$RepDate)==period[p]), 3:16], 2, sum) #<- as.numeric(ch)
}
pF <- apply(apply(newrawdata[, 3:9], MARGIN = c(1,2), as.numeric),
            MARGIN = 1, sum)
yF <- apply(apply(newrawdata[, 10:16], MARGIN = c(1,2), as.numeric)
            , MARGIN = 1, sum)
thetaconst$detG  <- pF/55268100
DATA$yG          <- yF
# VIROLOGY
rm(list=setdiff(ls(), c("DATA", "thetaconst", "period")))
viro.file <- paste(datadir, "Virology/PB Copy of RVU_w40tow20_1718_raw.csv", sep = "")
viro.date.fmt <- "%d-%b-%y"
regions <- c("East Midlands", "East of England", "London", "North East", "North West", "South East", "South West", "West Midlands", "Yorkshire & The Humber")## , "England")
composite.regions <- c("Midlands and East", "North", "England")
which.regions <- list("Midlands and East" = c("East Midlands", "East of England", "West Midlands"),
                      North = c("North East", "North West", "Yorkshire & The Humber"),
                      England = regions)
swab.data <- read.csv(viro.file, header = T, as.is = T)
swab.data <- swab.data[swab.data$Flu.Positive.Result %in% c("+","-"), ]
swab.data <- swab.data[swab.data$Sample.Date != "", ]
swab.data[, 1] <- as.character(swab.data[, 1])
swab.data[, 1] <- as.Date(swab.data[, 1], format = viro.date.fmt)
swab.data <- swab.data[swab.data$Sample.Date %in% period, ]
swab.data$Age.Group <- factor(swab.data$Age.Group, levels = c("<1", "1-4", "5-14", "15-24", "25-44", "45-64", "65+"), ordered = T)
swab.data$Sample.Date <- factor(swab.data$Sample.Date, levels = as.character(period))
denoms.list <- by(swab.data, swab.data[, 3], function(x) table(x$Sample.Date, x$Age.Group))
pos.list <- by(swab.data, swab.data[, 3], function(x) table(x$Sample.Date[x$Flu.Positive.Result == "+"],
                                                            x$Age.Group[x$Flu.Positive.Result == "+"]))
get.composite.regions <- function(region){
  idxs <- which(names(denoms.list) %in% which.regions[[region]])
  nc <- ncol(denoms.list[[1]])
  out <- list(denoms = NULL, pos = NULL)
  for(idx in idxs){
    out$denoms <- if(is.null(out$denoms)) denoms.list[[idx]] else out$denoms + denoms.list[[idx]]
    out$pos <- if(is.null(out$pos)) pos.list[[idx]] else out$pos + pos.list[[idx]]
  }
  out
}
add.data.list <- lapply(composite.regions, get.composite.regions)
names(add.data.list) <- composite.regions
denoms.list[names(add.data.list)] <- lapply(add.data.list, function(x) x$denoms)
pos.list[names(add.data.list)] <- lapply(add.data.list, function(x) x$pos)
denoms.list.agg <- lapply(denoms.list, function(x) matrix(apply(x, 1, sum), nrow(x), 1))
pos.list.agg <- lapply(pos.list, function(x) matrix(apply(x, 1, sum), nrow(x), 1))
length(denoms.list.agg)
length(pos.list.agg)
yV <- rep(0, length(period))
nV <- rep(0, length(period))
for (r in 1:length(denoms.list.agg)){
  yV <- yV + pos.list.agg[[r]]
  nV <- nV + denoms.list.agg[[r]]
}
thetaconst$sizeV <- as.vector(nV)
DATA$yV          <- as.vector(yV)


par(mfrow=c(4,2), mar=c(2,3,1,1))
plot(DATA$yH, type="l", main="yH")
plot(thetaconst$detH, col="gray", type="l", main="detH")
plot(DATA$yIC, type="l", main="yIC")
plot(thetaconst$detIC, col="gray", type="l", main="detIC")
plot(DATA$yG, type="l", main="yG")
plot(thetaconst$detG, col="gray", type="l", main="detIC")
plot(DATA$yV, type="l", main="yV")
lines(thetaconst$sizeV,col="gray")
plot(DATA$yV/thetaconst$sizeV, type="l", main="observed positive proportion")

# FUNCTIONS ------------------------------------------------------------

sourceCpp("CppAll.cpp")
source("Rfun.R")


# GIVEN QUANTITIES --------------------------------------------------------



thetaconst <- c(thetaconst, 
                list(N=55268100, dI=3.14, TMspd=4, TMbg=0, TMnd=231, OUTdt=7, 
                     DELstps=10, EfH=0.32, HfIC=0.4,
                     OUTdu=1, DELstpsGP=21, alEfG=3.41, btEfG=0.83))
H_tbp <- rep(0,  (thetaconst$TMnd-thetaconst$TMbg)*thetaconst$TMspd)
H_tbp[c((19*thetaconst$TMspd) :(28*thetaconst$TMspd),
        (75*thetaconst$TMspd) :(92*thetaconst$TMspd),
        (131*thetaconst$TMspd):(140*thetaconst$TMspd),
        (179*thetaconst$TMspd):(196*thetaconst$TMspd))] <- 1
thetaconst <- c(thetaconst, list(vecHtb=H_tbp))
rm(H_tbp)
thetaconst$dEtoF <-  fdisC("gamma" , c(thetaconst$alEfG, thetaconst$btEfG),  
                           thetaconst$OUTdu, thetaconst$DELstpsGP)
thetaconst$dEtoH <-  fdisC("exponential" , c(thetaconst$EfH),  
                           thetaconst$OUTdt, thetaconst$DELstps)
thetaconst$dHtoIC <-  fdisC("exponential" , c(thetaconst$HfIC),  
                            thetaconst$OUTdt, thetaconst$DELstps)


# constant things
thetaconst <- c(thetaconst, 
                list(N=55268100, dI=3.14, TMspd=4, TMbg=0, TMnd=231, OUTdt=7, DELstps=10, 
                     EfH=0.32, HfIC=0.4,  OUTdu=1, DELstpsGP=21, alEfG=3.41, btEfG=0.83))
H_tbp <- rep(0, (thetaconst$TMnd-thetaconst$TMbg)*thetaconst$TMspd)
H_tbp[c((19*thetaconst$TMspd) :(28*thetaconst$TMspd),
        (75*thetaconst$TMspd) :(92*thetaconst$TMspd),
        (131*thetaconst$TMspd):(140*thetaconst$TMspd),
        (179*thetaconst$TMspd):(196*thetaconst$TMspd))] <- 1
thetaconst       <- c(thetaconst, list(vecHtb=H_tbp))



# SET UP PARAMETERS -----------------------------------------------------------------------

DPT.thin  <- 20
DPT.nsmpl <- 500
SMP.thin  <- 100
SMP.nsmpl <- 1500
SMP.bin   <- 500

DPT.thinB  <- 20
DPT.nsmplB <- 5000
SMP.thinB  <- 40
SMP.nsmplB <- 1100
SMP.binB   <- 100

nparticle  <- 10000
# initial vectors
start1 <- toNorm(c(0.63, 0.35, 0.00005,  2,  0.1, #trans 
                   0.008,  0.05,   10,            #svrt H
                   4.6, -0.2, 1,                  #bkg
                   .005, 100,                     #svrtF
                   2, 1.8, 1.7, 1.8, 0.5, 0.1))   #dow
start2 <- toNorm(c( 0.548,0.304,0.00015,2.992,0.348,0.003,0.073,24.571, #tans svrtH
                    6.6, -0.1, 1.1,.001, 200,                           #bkg svrtF
                    2.1, 2, 1.6, 1.6, 0.01, 0.01))                      #dow
start3 <- toNorm(c( 0.596,0.384,0.00013,2.420,0.388,0.004,0.074,15.359, #tans svrtH
                    7., -0.25, 0.8,.01, 150,                            #bkg svrtF
                    1.8, 1.8, 1.8, 1.8, 0.05, 0.05))                    #dow


vecNAMES <- c(expression(beta),expression(pi),expression(iota),
              expression(d[L]),expression(kappa),
              expression(.^0*theta^{H}),expression(.^H*theta^{IC}), expression(epsilon^H),
              expression(nu[1]),expression(nu[2]),expression(nu[3]), 
              expression(iota^F),expression(epsilon^F),
              expression(omega[1],omega[2],omega[3],omega[5],omega[6],omega[7]))

# RUN ---------------------------------------------------------------------

dpt1 <- adptMHsingle(start=toNorm(c(0.63, 0.35, 0.00005,
                                    2,  0.1, 
                                    0.008,  0.05,   10,
                                    4.6, -0.2, 1,
                                    .005, 100,
                                    2, 1.8, 1.7, 1.8, 0.5, 0.1)), 
                     sdini= c(0.003, 0.01, 0.1, 0.21, 0.08, 0.08,0.1, 1.38,
                              0.01, 0.01, 0.01, 0.01, 0.01,
                              0.01)/2, 
                     nsmpl=DPT.nsmpl, thin=DPT.thin, stp.dpt=500,np=10000, THETACONST=thetaconst)
save.image("results/OUTall.RData")
dpt2 <- adptMHsingle(start=dpt1[[2]][DPT.nsmpl, ],sdini=dpt1[[1]], 
                     nsmpl=DPT.nsmpl, thin=DPT.thin, stp.dpt=500,np=10000, THETACONST=thetaconst)
save.image("results/OUTall.RData")
dpt3 <- adptMHsingle(start=dpt2[[2]][DPT.nsmpl, ],sdini=dpt2[[1]], 
                     nsmpl=DPT.nsmpl, thin=DPT.thin, stp.dpt=500,np=10000, THETACONST=thetaconst)
save.image("results/OUTall.RData")
dpt4 <- adptMHsingle(start=dpt3[[2]][DPT.nsmpl, ],sdini=dpt3[[1]], 
                     nsmpl=DPT.nsmpl, thin=DPT.thin, stp.dpt=500,np=10000, THETACONST=thetaconst)
save.image("results/OUTall.RData")


par(mfrow=c(7,3), mar=c(2, 2, 1, 0.1))
for (p in 1:19){
  plot(rbind(t(apply(dpt1[[2]], MARGIN = 1, toTh)),
             t(apply(dpt2[[2]], MARGIN = 1, toTh)),
             t(apply(dpt3[[2]], MARGIN = 1, toTh)),
             t(apply(dpt4[[2]], MARGIN = 1, toTh)))[,p], type="l", col=2, 
       main=vecNAMES[p])
}


# for parallel run this ---------------------------------------------

load("results/OUTall.RData")
library(mvtnorm)
# library(microbenchmark)
library(TailRank)
library(Rcpp)
library(fields)
library(readxl)
library(xlsx)
library(inline)
setwd(thisdir)
sourceCpp("CppAll.cpp")
source("Rfun.R")


# singlechains ------------------------------------------------------------

ch1.temp <- runMHsingle(start=start1, sd=dpt4[[1]], nsmpl=SMP.nsmpl/15, thin=SMP.thin, 
                        np=nparticle, THETACONST=thetaconst)
ch1 <- list(ch1.temp[[1]], ch1.temp[[2]], ch1.temp[[3]])
save(ch1, file = "results/ch1all.Robj")
for(split in 2:15){
  ch1.temp <- runMHsingle(start=ch1.temp[[2]][(SMP.nsmpl/15),], sd=dpt4[[1]], nsmpl=SMP.nsmpl/15, thin=SMP.thin, 
                          np=nparticle, THETACONST=thetaconst)
  ch1 <- list(rbind(ch1[[1]], ch1.temp[[1]]), 
              rbind(ch1[[2]], ch1.temp[[2]]), 
              rbind(ch1[[3]], ch1.temp[[3]]))
  save(ch1, file = "results/ch1all.Robj")
  print(paste("chunck",split, "of 15"))
}
rm(ch1.temp)

ch2.temp <- runMHsingle(start=start2, sd=dpt4[[1]], nsmpl=SMP.nsmpl/15, thin=SMP.thin, 
                        np=nparticle, THETACONST=thetaconst)
ch2 <- list(ch2.temp[[1]], ch2.temp[[2]], ch2.temp[[3]])
save(ch2, file = "results/ch2all.Robj")
for(split in 2:15){
  ch2.temp <- runMHsingle(start=ch2.temp[[2]][(SMP.nsmpl/15),], sd=dpt4[[1]], nsmpl=SMP.nsmpl/15, thin=SMP.thin, 
                          np=nparticle, THETACONST=thetaconst)
  ch2 <- list(rbind(ch2[[1]], ch2.temp[[1]]), 
              rbind(ch2[[2]], ch2.temp[[2]]), 
              rbind(ch2[[3]], ch2.temp[[3]]))
  save(ch2, file = "results/ch2all.Robj")
  print(paste("chunck",split, "of 15"))
}
rm(ch2.temp)

ch3.temp <- runMHsingle(start=start3, sd=dpt4[[1]], nsmpl=SMP.nsmpl/15, thin=SMP.thin, 
                        np=nparticle, THETACONST=thetaconst)
ch3 <- list(ch3.temp[[1]], ch3.temp[[2]], ch3.temp[[3]])
save(ch3, file = "results/ch3all.Robj")
for(split in 2:15){
  ch3.temp <- runMHsingle(start=ch3.temp[[2]][(SMP.nsmpl/15),], sd=dpt4[[1]], nsmpl=SMP.nsmpl/15, thin=SMP.thin, 
                          np=nparticle, THETACONST=thetaconst)
  ch3 <- list(rbind(ch3[[1]], ch3.temp[[1]]), 
              rbind(ch3[[2]], ch3.temp[[2]]), 
              rbind(ch3[[3]], ch3.temp[[3]]))
  save(ch3, file = "results/ch3all.Robj")
  print(paste("chunck",split, "of 15"))
}

rm(ch3.temp)

# parallel from here --------------------------------------------------

load("results/ch1all.Robj")
load("results/ch2all.Robj")
load("results/ch3all.Robj")

par(mfrow=c(7,3), mar=c(2, 2, 1, 0.1))
for (p in 1:19){
  plot(t(apply(ch1[[2]], MARGIN = 1, toTh))[,p], type="l", col=2, 
       main=vecNAMES[p],
       ylim=c(min(c(t(apply(ch1[[2]], MARGIN = 1, toTh))[,p], 
                    t(apply(ch2[[2]], MARGIN = 1, toTh))[,p], 
                    t(apply(ch3[[2]], MARGIN = 1, toTh))[,p])),
              max(c(t(apply(ch1[[2]], MARGIN = 1, toTh))[,p], 
                    t(apply(ch2[[2]], MARGIN = 1, toTh))[,p], 
                    t(apply(ch3[[2]], MARGIN = 1, toTh))[,p]))))
  lines(t(apply(ch2[[2]], MARGIN = 1, toTh))[,p], col=3)
  lines(t(apply(ch3[[2]], MARGIN = 1, toTh))[,p], col=4)
}

apply(ch1[[3]],2,sum)/(SMP.nsmpl*SMP.thin)
apply(ch2[[3]],2,sum)/(SMP.nsmpl*SMP.thin)
apply(ch3[[3]],2,sum)/(SMP.nsmpl*SMP.thin)


# vcov matrix

vcB <- cov(rbind(ch1[[2]][(SMP.bin+1):SMP.nsmpl, ],
                 ch2[[2]][(SMP.bin+1):SMP.nsmpl, ], 
                 ch3[[2]][(SMP.bin+1):SMP.nsmpl, ]))

save.image("results/OUTall.RData")

# BLOCKED CHAINS ----------------------------------------------------------


DPT.thinB  <- 20
DPT.nsmplB <- 5000
SMP.thinB  <- 200
SMP.nsmplB <- 5500
SMP.binB   <- 500

dpt1BLK1 <- adptMHblocked(start=ch1[[2]][SMP.nsmpl, ], 
                          vcini=vcB,nsmpl=DPT.nsmplB,thin=DPT.thinB, 
                          stp.dpt=500,np=10000, THETACONST=thetaconst)
save.image("results/OUTall.RData")


par(mfrow=c(7,3), mar=c(2, 2, 1, 0.1))
for (p in 1:19){
  plot(t(apply(dpt1BLK1[[2]], MARGIN = 1, toTh))[,p], type="l", col=2, 
       main=vecNAMES[p])
}

# PARALLEL CHAINS CONTINUE HERE

chBLK1.temp <- runMHblocked(start=start1, vc=dpt1BLK1[[1]]*vcB, nsmpl=SMP.nsmplB/11, 
                       thin=SMP.thinB, np=nparticle, THETACONST=thetaconst)
chBLK1 <- list(chBLK1.temp[[1]], chBLK1.temp[[2]], chBLK1.temp[[3]])
for (split in 2:11){
  chBLK1.temp <- runMHblocked(start=chBLK1.temp[[2]][(SMP.nsmplB/11),], 
                             vc=dpt1BLK1[[1]]*vcB, nsmpl=SMP.nsmplB/11, 
                       thin=SMP.thinB, np=nparticle, THETACONST=thetaconst)
  chBLK1 <- list(rbind(chBLK1[[1]], chBLK1.temp[[1]]), 
              rbind(chBLK1[[2]], chBLK1.temp[[2]]), 
              rbind(chBLK1[[3]], chBLK1.temp[[3]]))
  save(chBLK1, file = "results/ch1Ball.Robj")
  print(paste("chunck",split, "of 11"))
}
rm(chBLK1.temp)

chBLK2.temp <- runMHblocked(start=start2, vc=dpt1BLK1[[1]]*vcB, nsmpl=SMP.nsmplB/11, 
                       thin=SMP.thinB, np=nparticle, THETACONST=thetaconst)
chBLK2 <- list(chBLK2.temp[[1]], chBLK2.temp[[2]], chBLK2.temp[[3]])
for (split in 2:11){
  chBLK2.temp <- runMHblocked(start=chBLK2.temp[[2]][(SMP.nsmplB/11),], 
                             vc=dpt1BLK1[[1]]*vcB, nsmpl=SMP.nsmplB/11, 
                       thin=SMP.thinB, np=nparticle, THETACONST=thetaconst)
  chBLK2 <- list(rbind(chBLK2[[1]], chBLK2.temp[[1]]), 
              rbind(chBLK2[[2]], chBLK2.temp[[2]]), 
              rbind(chBLK2[[3]], chBLK2.temp[[3]]))
  save(chBLK2, file = "results/ch2Ball.Robj")
  print(paste("chunck",split, "of 11"))
}
rm(chBLK2.temp)

chBLK3.temp <- runMHblocked(start=start3, vc=dpt1BLK1[[1]]*vcB, nsmpl=SMP.nsmplB/11, 
                       thin=SMP.thinB, np=nparticle, THETACONST=thetaconst)
chBLK3 <- list(chBLK3.temp[[1]], chBLK3.temp[[2]], chBLK3.temp[[3]])
for (split in 2:11){
  chBLK3.temp <- runMHblocked(start=chBLK3.temp[[2]][(SMP.nsmplB/11),], 
                             vc=dpt1BLK1[[1]]*vcB, nsmpl=SMP.nsmplB/11, 
                       thin=SMP.thinB, np=nparticle, THETACONST=thetaconst)
  chBLK3 <- list(rbind(chBLK3[[1]], chBLK3.temp[[1]]), 
              rbind(chBLK3[[2]], chBLK3.temp[[2]]), 
              rbind(chBLK3[[3]], chBLK3.temp[[3]]))
  save(chBLK3, file = "results/ch3Ball.Robj")
  print(paste("chunck",split, "of 11"))
}
rm(chBLK3.temp)



# parallel, up until here

load("results/ch1Ball.Robj")
load("results/ch2Ball.Robj")
load("results/ch3Ball.Robj")
par(mfrow=c(7,3), mar=c(2, 2, 1, 0.1))
for (p in 1:19){
  plot(t(apply(chBLK1[[2]], MARGIN = 1, toTh))[,p], type="l", col=2, 
       main=vecNAMES[p],
       ylim=c(min(c(t(apply(chBLK1[[2]], MARGIN = 1, toTh))[,p], 
                    t(apply(chBLK2[[2]], MARGIN = 1, toTh))[,p], 
                    t(apply(chBLK3[[2]], MARGIN = 1, toTh))[,p])),
              max(c(t(apply(chBLK1[[2]], MARGIN = 1, toTh))[,p], 
                    t(apply(chBLK2[[2]], MARGIN = 1, toTh))[,p], 
                    t(apply(chBLK3[[2]], MARGIN = 1, toTh))[,p]))))
  lines(t(apply(chBLK2[[2]], MARGIN = 1, toTh))[,p], col=3)
  lines(t(apply(chBLK3[[2]], MARGIN = 1, toTh))[,p], col=4)
}

sum(chBLK1[[3]])/(SMP.nsmpl*SMP.thin)
sum(chBLK2[[3]])/(SMP.nsmpl*SMP.thin)
sum(chBLK3[[3]])/(SMP.nsmpl*SMP.thin)


# PLOTS SAVE --------------------------------------------------------------

pdf("OUTPUT1all.pdf")


# adaptive 
par(mfrow=c(7,3), mar=c(2, 2, 1, 0.1))
for (p in 1:19){
  plot(t(apply(dpt1[[2]], MARGIN = 1, toTh))[,p], type="l", col=2, 
       main=vecNAMES[p])
}

par(mfrow=c(7,3), mar=c(2, 2, 1, 0.1))
for (p in 1:19){
  plot(t(apply(ch1[[2]], MARGIN = 1, toTh))[,p], type="l", col=2, 
       main=vecNAMES[p],
       ylim=c(min(c(t(apply(ch1[[2]], MARGIN = 1, toTh))[,p], 
                    t(apply(ch2[[2]], MARGIN = 1, toTh))[,p], 
                    t(apply(ch3[[2]], MARGIN = 1, toTh))[,p])),
              max(c(t(apply(ch1[[2]], MARGIN = 1, toTh))[,p], 
                    t(apply(ch2[[2]], MARGIN = 1, toTh))[,p], 
                    t(apply(ch3[[2]], MARGIN = 1, toTh))[,p]))))
  lines(t(apply(ch2[[2]], MARGIN = 1, toTh))[,p], col=3)
  lines(t(apply(ch3[[2]], MARGIN = 1, toTh))[,p], col=4)
}

par(mfrow=c(7,3), mar=c(2, 2, 1, 0.1))
for (p in 1:19){
  plot(t(apply(dpt1BLK1[[2]], MARGIN = 1, toTh))[,p], type="l", col=2, 
       main=vecNAMES[p])
}



par(mfrow=c(7,3), mar=c(2, 2, 1, 0.1))
for (p in 1:19){
  plot(t(apply(chBLK1[[2]], MARGIN = 1, toTh))[,p], type="l", col=2, 
       main=vecNAMES[p],
       ylim=c(min(c(t(apply(chBLK1[[2]], MARGIN = 1, toTh))[,p], 
                    t(apply(chBLK2[[2]], MARGIN = 1, toTh))[,p], 
                    t(apply(chBLK3[[2]], MARGIN = 1, toTh))[,p])),
              max(c(t(apply(chBLK1[[2]], MARGIN = 1, toTh))[,p], 
                    t(apply(chBLK2[[2]], MARGIN = 1, toTh))[,p], 
                    t(apply(chBLK3[[2]], MARGIN = 1, toTh))[,p]))))
  lines(t(apply(chBLK2[[2]], MARGIN = 1, toTh))[,p], col=3)
  lines(t(apply(chBLK3[[2]], MARGIN = 1, toTh))[,p], col=4)
}


chTh1 <- t(apply(chBLK1[[2]], MARGIN = 1, toTh))
chTh2 <- t(apply(chBLK2[[2]], MARGIN = 1, toTh))
chTh3 <- t(apply(chBLK3[[2]], MARGIN = 1, toTh))


par(mfrow=c(13, 13), mar=rep(0.1,4))
for (i in c(1:13)){
  for (j in c(1:13)){
    if (i==j){
      plot(0,0, col="white", xaxt="n", yaxt="n")
      text(0,0,vecNAMES[j], cex=4)
    }else{plot(chTh1[seq(110,1100, by=10), j],chTh1[seq(110,1100, by=10) , i],col="red", pch=20, cex=0.1, xaxt="n", yaxt="n",
               xlim=c(min(c((chTh1)[seq(110,1100, by=10),j], (chTh2)[seq(110,1100, by=10),j], (chTh3)[seq(110,1100, by=10),j])) ,
                      max(c((chTh1)[seq(110,1100, by=10),j], (chTh2)[seq(110,1100, by=10),j], (chTh3)[seq(110,1100, by=10),j]))), 
               ylim=c(min(c((chTh1)[seq(110,1100, by=10),i], (chTh2)[seq(110,1100, by=10),i], (chTh3)[seq(110,1100, by=10),i])) , 
                      max(c((chTh1)[seq(110,1100, by=10),i], (chTh2)[seq(110,1100, by=10),i], (chTh3)[seq(110,1100, by=10),i]))))
      points(chTh2[seq(110,1100, by=10), j],chTh2[seq(110,1100, by=10) , i],col="green", pch=20, cex=0.1)
      points(chTh3[seq(110,1100, by=10), j],chTh3[seq(110,1100, by=10), i],col="blue", pch=20, cex=0.1)
    }
  }
}



par(mfrow=c(6,6), mar=rep(0.1,4))
for (i in c(9:14)){
  for (j in c(9:14)){
    if (i==j){
      plot(0,0, col="white", xaxt="n", yaxt="n")
      text(0,0,c(rep(NA, 8), expression(omega[1]),expression(omega[2]),expression(omega[3]),
                 expression(omega[5]),expression(omega[6]),expression(omega[7]))[j], cex=4)
    }else{plot(chTh1[seq(110,1100, by=10), j],chTh1[seq(110,1100, by=10) , i],col="red", pch=20, cex=0.1, xaxt="n", yaxt="n",
               xlim=c(min(c((chTh1)[seq(110,1100, by=10),j], (chTh2)[seq(110,1100, by=10),j], (chTh3)[seq(110,1100, by=10),j])) ,
                      max(c((chTh1)[seq(110,1100, by=10),j], (chTh2)[seq(110,1100, by=10),j], (chTh3)[seq(110,1100, by=10),j]))), 
               ylim=c(min(c((chTh1)[seq(110,1100, by=10),i], (chTh2)[seq(110,1100, by=10),i], (chTh3)[seq(110,1100, by=10),i])) , 
                      max(c((chTh1)[seq(110,1100, by=10),i], (chTh2)[seq(110,1100, by=10),i], (chTh3)[seq(110,1100, by=10),i]))))
      points(chTh2[seq(110,1100, by=10), j],chTh2[seq(110,1100, by=10) , i],col="green", pch=20, cex=0.1)
      points(chTh3[seq(110,1100, by=10), j],chTh3[seq(110,1100, by=10), i],col="blue", pch=20, cex=0.1)
    }
  }
}
# PRIOR to POSTERIOR DISTRIBUTION -----------------------------------------
# 
# big.sample <- rbind(chBLK1[[2]][(SMP.binB+1):SMP.nsmplB, ], 
#                     chBLK2[[2]][(SMP.binB+1):SMP.nsmplB, ], 
#                     chBLK3[[2]][(SMP.binB+1):SMP.nsmplB, ])
big.sample <- ch2[[2]][(SMP.bin+1):SMP.nsmplB, ]
th.bs <- t(apply(big.sample, 1, toTh))


par(mfrow=c(1,1), mar=rep(4,4))
# beta is a uniform between [0,4]
hist(th.bs[,1], main=expression(beta), 
     col=rgb(0/255,128/255, 128/255, 1), border="gray", lwd=2, freq=F)
curve(dunif(x, 0, 4, log=F), add=T, lwd=2, col=rgb(213/355 ,    0/255,  50/255, 1))
#pi is a beta with mean 0.375  
hist(th.bs[,2], main=expression(pi), 
     col=rgb(0/255,128/255, 128/255, 1), border="gray", lwd=2, freq=F)
curve(dbeta(x, shape1 =37.5 , shape2=62.5, log=F), add=T, lwd=2, col=rgb(213/355 ,    0/255,  50/255, 1))
# iota is a uniform between [0,0.1]
hist(th.bs[,3], main=expression(iota), 
     col=rgb(0/255,128/255, 128/255, 1), border="gray", lwd=2, freq=F)
curve(dunif(x, 0, 0.1, log=F), add=T, lwd=2, col=rgb(213/355 ,    0/255,  50/255, 1))
# dL is a LogNormal Centered in 2
hist(th.bs[,4], main=expression(d[L]), 
     col=rgb(0/255,128/255, 128/255, 1), border="gray", lwd=2, freq=F)
curve(dlnorm(x, meanlog = log(2), sd=0.5, log=F), add=T, lwd=2, col=rgb(213/355 ,    0/255,  50/255, 1))
# K is a uniform [0,2] 
hist(th.bs[,5], main=expression(kappa), 
     col=rgb(0/255,128/255, 128/255, 1), border="gray", lwd=2, freq=F)
curve( dnorm(log(x+1), 0, sd=1, log=F), add=T, lwd=2, col=rgb(213/355 ,    0/255,  50/255, 1))
# thetaH is a uniform [0,1]
hist(th.bs[,6], main=vecNAMES[6], 
     col=rgb(0/255,128/255, 128/255, 1), border="gray", lwd=2, freq=F)
curve(dunif(x, 0, 1, log=F), add=T, lwd=2, col=rgb(213/355 ,    0/255,  50/255, 1))
# HthIC is a uniform [0,1]
hist(th.bs[,7], main=vecNAMES[7], 
     col=rgb(0/255,128/255, 128/255, 1), border="gray", lwd=2, freq=F)
curve(dunif(x, 0, 1, log=F), add=T, lwd=2, col=rgb(213/355 ,    0/255,  50/255, 1))
# epsH is a is a uniform [0,100]
hist(th.bs[,8], main=vecNAMES[8], 
     col=rgb(0/255,128/255, 128/255, 1), border="gray", lwd=2, freq=F)
curve(dunif(x, 0, 100, log=F), add=T, lwd=2, col=rgb(213/355 ,    0/255,  50/255, 1))
# nus are normal with historical priors
hist(th.bs[,9], main=vecNAMES[9], 
     col=rgb(0/255,128/255, 128/255, 1), border="gray", lwd=2, freq=F)
curve(dnorm(x, mean = 6.6, sd=0.17, log=F), add=T, lwd=2, col=rgb(213/355 ,    0/255,  50/255, 1))
hist(th.bs[,10], main=vecNAMES[10], 
     col=rgb(0/255,128/255, 128/255, 1), border="gray", lwd=2, freq=F)
curve(dnorm(x, mean = -0.2, sd=0.11, log=F), add=T, lwd=2, col=rgb(213/355 ,    0/255,  50/255, 1))
hist(th.bs[,11], main=vecNAMES[11], 
     col=rgb(0/255,128/255, 128/255, 1), border="gray", lwd=2, freq=F)
curve(dnorm(x, mean = 0.99, sd=0.08, log=F), add=T, lwd=2, col=rgb(213/355 ,    0/255,  50/255, 1))

# iota is a uniform 0,1
hist(th.bs[,12], main=vecNAMES[12], 
     col=rgb(0/255,128/255, 128/255, 1), border="gray", lwd=2, freq=F)
curve(dunif(x, 0, 1, log=F), add=T, lwd=2, col=rgb(213/355 ,    0/255,  50/255, 1))
# epsF is a uniform 0, 10
hist(th.bs[,13], main=vecNAMES[13], 
     col=rgb(0/255,128/255, 128/255, 1), border="gray", lwd=2, freq=F)
curve(dunif(x, 0, 5000, log=F), add=T, lwd=2, col=rgb(213/355 ,    0/255,  50/255, 1))
# the omegas log normals
hist(th.bs[,14], main=expression(omega[1]), 
     col=rgb(0/255,128/255, 128/255, 1), border="gray", lwd=2, freq=F)
curve(dlnorm(x, meanlog = log(1), sd=0.4, log=F), add=T, lwd=2, col=rgb(213/355 ,    0/255,  50/255, 1))
hist(th.bs[,15], main=expression(omega[2]), 
     col=rgb(0/255,128/255, 128/255, 1), border="gray", lwd=2, freq=F)
curve(dlnorm(x, meanlog = log(1), sd=0.4, log=F), add=T, lwd=2, col=rgb(213/355 ,    0/255,  50/255, 1))
hist(th.bs[,16], main=expression(omega[3]), 
     col=rgb(0/255,128/255, 128/255, 1), border="gray", lwd=2, freq=F)
curve(dlnorm(x, meanlog = log(1), sd=0.4, log=F), add=T, lwd=2, col=rgb(213/355 ,    0/255,  50/255, 1))
hist(th.bs[,17], main=expression(omega[5]), 
     col=rgb(0/255,128/255, 128/255, 1), border="gray", lwd=2, freq=F)
curve(dlnorm(x, meanlog = log(1), sd=0.4, log=F), add=T, lwd=2, col=rgb(213/355 ,    0/255,  50/255, 1))
hist(th.bs[,18], main=expression(omega[6]), 
     col=rgb(0/255,128/255, 128/255, 1), border="gray", lwd=2, freq=F)
curve(dlnorm(x, meanlog = log(1), sd=0.4, log=F), add=T, lwd=2, col=rgb(213/355 ,    0/255,  50/255, 1))
hist(th.bs[,19], main=expression(omega[7]), 
     col=rgb(0/255,128/255, 128/255, 1), border="gray", lwd=2, freq=F)
curve(dlnorm(x, meanlog = log(1), sd=0.4, log=F), add=T, lwd=2, col=rgb(213/355 ,    0/255,  50/255, 1))

# OTHER TRANSFORMATION OF THE PARAMETERS ----------------------------------

# RO when kappa is not effective
R0 <- th.bs[,1]*thetaconst$dI
hist(R0, main=expression(R[0]), 
     col=rgb(0/255,128/255, 128/255, 1), border="gray", lwd=2, freq=F)
# Re when kappa is not effective
Re <- th.bs[,1]*thetaconst$dI*(1-th.bs[,2]-th.bs[,3])
hist(Re, main=expression(R[e]), breaks=30,
     col=rgb(0/255,128/255, 128/255, 1), border="gray", lwd=2, freq=F)

# RO when kappa is effective
R0 <- th.bs[,1]*thetaconst$dI*(1+th.bs[,5])
hist(R0, main=expression(paste(R[0], " when ",kappa )),
     col=rgb(0/255,128/255, 128/255, 1), border="gray", lwd=2, freq=F)
# Re when kappa is effective
Re <- th.bs[,1]*thetaconst$dI*(1-th.bs[,2]-th.bs[,3])*(1+th.bs[,5])
hist(Re, main=expression(paste(R[e], " when ",kappa )), breaks=30,
     col=rgb(0/255,128/255, 128/255, 1), border="gray", lwd=2, freq=F)

# thetha0Fprior and posterior
epsPrior <- runif(10000, 0, 5000) 
iotPrior <- runif(10000, 0, 1) 
ThGPrior <- rgamma(10000, shape= iotPrior*epsPrior, rate = epsPrior)
ThGPoste <- rgamma(dim(th.bs)[1], shape= th.bs[,12]*th.bs[,13], rate = th.bs[,13])

hist(ThGPoste, col=rgb(0/255,128/255, 128/255, 1), border="gray", lwd=2, freq=F, xlim=c(0, 0.1))
hist(ThGPrior, col=rgb(213/355 ,    0/255,  50/255, 0.4), border="gray", lwd=1, add=T, freq=F)

hist(ThGPrior, col=rgb(213/355 ,    0/255,  50/255, 0.4), border="gray", lwd=1, freq=F)

# the background:
guP <- t(apply(th.bs[,9:11], MARGIN = 1, nuTOg))
dim(guP)

MEguP <- apply(guP, MARGIN = 2, median)
LOguP <- apply(guP, MARGIN = 2, quantile, prob=0.025)
UPguP <- apply(guP, MARGIN = 2, quantile, prob=0.975)


plot(UPguP, type="l", ylim=c(0,1500), col=rgb(0/255,128/255, 128/255,1))
lines(LOguP, col=rgb(0/255,128/255, 128/255,1))
lines(MEguP, lwd=2, col=rgb(0/255,128/255, 128/255,1))
lines(nuTOg(c(6.6, -0.2,  0.99)), 
      col=rgb(213/355 ,    0/255,  50/255, 1), lwd=2)
lines(nuTOg(c(qnorm(0.025, 6.6, 0.17), qnorm(0.025, -.2, 0.11), qnorm(0.025, 0.99, 0.08))), 
      col=rgb(213/355 ,    0/255,  50/255, 1), lwd=1)
lines(nuTOg(c(qnorm(0.975, 6.6, 0.17), qnorm(0.975, -.2, 0.11), qnorm(0.975, 0.99, 0.08))), 
      col=rgb(213/355 ,    0/255,  50/255, 1), lwd=1)
lines(DATA$yG/thetaconst$detG, col=rgb(0,0, 0,.2))
legend(150, 1500, c("prior gu", "posterior gu", "data on pop scale"), 
       col=c(rgb(213/355 ,    0/255,  50/255, 1),
             rgb(0/255,128/255, 128/255,1),rgb(0,0, 0,.2)), lwd=2)


xiFuP <- matrix(NA, ncol=(thetaconst$TMbg+thetaconst$TMnd), nrow=(dim(th.bs)[1]))
# chck data fitting
for (s in 1:(dim(th.bs)[1])){
  lamu <- LamTseeiirCtv(tm_X0 = setNames(c(p_S=thetaconst$N*(1-th.bs[s,2]-th.bs[s,3]),
                                           p_E=thetaconst$N*(th.bs[s,3]*((th.bs[s,4])/(th.bs[s,4]+thetaconst$dI))),
                                           p_I=thetaconst$N*(th.bs[s,3]*((thetaconst$dI)/(th.bs[s,4]+thetaconst$dI))),
                                           p_R=thetaconst$N*(th.bs[s,2])),c("p_S","p_E","p_R","p_R")),
                        tm_Theta = c(bt=th.bs[s,1],  sgm=2/th.bs[s,4], gmm=2/thetaconst$dI),
                        tm_begin = thetaconst$TMbg, tm_end = thetaconst$TMnd,
                        tm_spd = thetaconst$TMspd, outdt = thetaconst$OUTdu,
                        tm_Tk = thetaconst$vecHtb, tm_k = th.bs[s,5])
  xiFuP[s, ] <- projCdet(probs=thetaconst$dEtoF, sizes = lamu)
}

plot(xiFuP[1, ], type="l", ylim=c(min(xiFuP), max(xiFuP)), 
     main="Some posterior smaples of the Number of New Infection")
for (s in 1:20){
  lines(xiFuP[s,])
  # Sys.sleep(1)
}
abline(v=which(thetaconst$vecHtb==1)/4, lwd=2, col=rgb(1,0,0,0.1))

# guP <- t(apply(th.bs[,6:8], MARGIN = 1, nuTOg))
# dim(guP)

rUP    <- th.bs[,12]*th.bs[,13]+guP*th.bs[,13]/xiFuP

omgUP  <- t(apply(th.bs[,14:19], MARGIN = 1, OMGfun))
dim(omgUP)

etaUP  <- matrix(NA, ncol=(thetaconst$TMbg+thetaconst$TMnd), nrow=(dim(th.bs)[1]))
for (s in 1:(dim(th.bs)[1])){
  etaUP[s,]  <- (th.bs[s,13]+thetaconst$detG*omgUP[s,]*xiFuP[s,])/th.bs[s,13]
}


alBBP  <- th.bs[,12]*th.bs[,13]
btBBP  <- (guP*th.bs[,13])/(xiFuP)



tspan <- (thetaconst$TMnd-thetaconst$TMbg)/thetaconst$OUTdu

ppyG <- matrix(NA, nrow=dim(th.bs)[1], ncol=tspan)
ppyV <- matrix(NA, nrow=dim(th.bs)[1], ncol=tspan)

for(t in 1: tspan){
  ppyG[, t] <- rnbinom((dim(th.bs)[1]), size=rUP[,t], prob = 1/etaUP[,t])
  ppyV[, t] <- rbb((dim(th.bs)[1]), N=thetaconst$sizeV[t], u=alBBP, v=btBBP[,t])
}

# summaries per each time point
SUMMyG <- matrix(NA, nrow=tspan, ncol=4)
SUMMyV <- matrix(NA, nrow=tspan, ncol=4)
colnames(SUMMyG) <- c("mean", "Me", "LoCrI", "UpCrI")
colnames(SUMMyV) <- c("mean", "Me", "LoCrI", "UpCrI")

SUMMyG[, 1] <- apply(ppyG, MARGIN = 2, mean)
SUMMyG[, 2] <- apply(ppyG, MARGIN = 2, median)
SUMMyG[, 3] <- apply(ppyG, MARGIN = 2, quantile, prob=0.025)
SUMMyG[, 4] <- apply(ppyG, MARGIN = 2, quantile, prob=0.975)
SUMMyV[, 1] <- apply(ppyV, MARGIN = 2, mean)
SUMMyV[, 2] <- apply(ppyV, MARGIN = 2, median)
SUMMyV[, 3] <- apply(ppyV, MARGIN = 2, quantile, prob=0.025)
SUMMyV[, 4] <- apply(ppyV, MARGIN = 2, quantile, prob=0.975)

plot(DATA$yG, col="white", main="Goodness of fit GP consultation data")
for (t in 1: tspan){
  points(t, SUMMyG[t, 2], col="yellowgreen")
  segments(t, SUMMyG[t, 3], t, SUMMyG[t, 4], col="yellowgreen")
  points(t, SUMMyG[t, 1], col="limegreen", pch="*")
  points(t, DATA$yG[t], col="magenta", pch=19)
}
legend(180, 2300, c("Data", "Me(CrI)", "Mean"), pch = c(19, 1, 8), 
       col=c("magenta","yellowgreen", "limegreen"), lwd=c(NA,1,NA))



plot(DATA$yG+1, col="white", log="y", main="Goodness of fit GP consultation data (log scale)")
for (t in 1: tspan){
  points(t, SUMMyG[t, 2]+1, col="yellowgreen")
  segments(t, SUMMyG[t, 3]+1, t, SUMMyG[t, 4]+1, col="yellowgreen")
  points(t, SUMMyG[t, 1]+1, col="limegreen", pch="*")
  points(t, DATA$yG[t]+1, col="magenta", pch=19)
}

legend(180, 2300, c("Data", "Me(CrI)", "Mean"), pch = c(19, 1, 8), 
       col=c("magenta","yellowgreen", "limegreen"), lwd=c(NA,1,NA))

plot(DATA$yV, col="white", main="Goodness of fit virology data")
for (t in 1: tspan){
  points(t, SUMMyV[t, 2], col="yellowgreen")
  segments(t, SUMMyV[t, 3], t, SUMMyV[t, 4], col="yellowgreen")
  points(t, SUMMyV[t, 1], col="limegreen", pch=8)
  points(t, DATA$yV[t], col="magenta", pch=19)
}
legend(180, 35, c("Data", "Me(CrI)", "Mean"), pch = c(19, 1, 8), 
       col=c("magenta","yellowgreen", "limegreen"), lwd=c(NA,1,NA))


tspanHIC <- (thetaconst$TMbg+thetaconst$TMnd)/thetaconst$OUTdt
xiH   <- matrix(NA, ncol=tspanHIC, nrow=(dim(th.bs)[1]))
xiIC  <- matrix(NA, ncol=tspanHIC, nrow=(dim(th.bs)[1]))
ppxH  <- matrix(NA, nrow=dim(th.bs)[1], ncol=tspanHIC)
ppyH  <- matrix(NA, nrow=dim(th.bs)[1], ncol=tspanHIC)
ppyIC <- matrix(NA, nrow=dim(th.bs)[1], ncol=tspanHIC)
# chck data fitting
for (s in 1:(dim(th.bs)[1])){
  lamu <- LamTseeiirCtv(tm_X0 = setNames(c(p_S=thetaconst$N*(1-th.bs[s,2]-th.bs[s,3]),
                                           p_E=thetaconst$N*(th.bs[s,3]*((th.bs[s,4])/(th.bs[s,4]+thetaconst$dI))),
                                           p_I=thetaconst$N*(th.bs[s,3]*((thetaconst$dI)/(th.bs[s,4]+thetaconst$dI))),
                                           p_R=thetaconst$N*(th.bs[s,2])),c("p_S","p_E","p_R","p_R")),
                        tm_Theta = c(bt=th.bs[s,1],  sgm=2/th.bs[s,4], gmm=2/thetaconst$dI),
                        tm_begin = thetaconst$TMbg, tm_end = thetaconst$TMnd,
                        tm_spd = thetaconst$TMspd, outdt = thetaconst$OUTdu,
                        tm_Tk = thetaconst$vecHtb, tm_k = th.bs[s,5])
  lamt  <- tolamt(lamu)
  xiH[s, ]  <- projCdet(probs=thetaconst$dEtoH, sizes = lamt)*th.bs[s,6]
  xiIC[s, ] <- projCdet(probs=thetaconst$dHtoIC, sizes = xiH[s, ])*th.bs[s,7]
  ppyIC[s, ]<- rpois(tspanHIC, as.vector(xiIC[s,]*thetaconst$detIC))
  ppxH[s, ] <- rpois(tspanHIC, as.vector(xiH[s, ]))
  ppyH[s, ] <- rbb(tspanHIC, N=ppxH[s, ], 
                   u= (thetaconst$detH/(1-thetaconst$detH))*th.bs[s, 8],th.bs[s, 8] )
}



# summaries per each time point
SUMMyH <- matrix(NA, nrow=tspanHIC, ncol=4)
SUMMyIC <- matrix(NA, nrow=tspanHIC, ncol=4)
colnames(SUMMyH) <- c("mean", "Me", "LoCrI", "UpCrI")
colnames(SUMMyIC) <- c("mean", "Me", "LoCrI", "UpCrI")

SUMMyH[, 1] <- apply(ppyH, MARGIN = 2, mean)
SUMMyH[, 2] <- apply(ppyH, MARGIN = 2, median)
SUMMyH[, 3] <- apply(ppyH, MARGIN = 2, quantile, prob=0.025)
SUMMyH[, 4] <- apply(ppyH, MARGIN = 2, quantile, prob=0.975)
SUMMyIC[, 1] <- apply(ppyIC, MARGIN = 2, mean)
SUMMyIC[, 2] <- apply(ppyIC, MARGIN = 2, median)
SUMMyIC[, 3] <- apply(ppyIC, MARGIN = 2, quantile, prob=0.025)
SUMMyIC[, 4] <- apply(ppyIC, MARGIN = 2, quantile, prob=0.975)

plot(DATA$yH, col="white", main="Goodness of fit HOSPITAL consultation data")
for (t in 1: tspanHIC){
  points(t, SUMMyH[t, 2], col="yellowgreen")
  segments(t, SUMMyH[t, 3], t, SUMMyH[t, 4], col="yellowgreen")
  points(t, SUMMyH[t, 1], col="limegreen", pch="*")
  points(t, DATA$yH[t], col="magenta", pch=19)
}
legend(180, 2300, c("Data", "Me(CrI)", "Mean"), pch = c(19, 1, 8), 
       col=c("magenta","yellowgreen", "limegreen"), lwd=c(NA,1,NA))




plot(DATA$yIC, col="white", main="Goodness of fit IC data")
for (t in 1: tspanHIC){
  points(t, SUMMyIC[t, 2], col="yellowgreen")
  segments(t, SUMMyIC[t, 3], t, SUMMyIC[t, 4], col="yellowgreen")
  points(t, SUMMyIC[t, 1], col="limegreen", pch=8)
  points(t, DATA$yIC[t], col="magenta", pch=19)
}
legend(180, 35, c("Data", "Me(CrI)", "Mean"), pch = c(19, 1, 8), 
       col=c("magenta","yellowgreen", "limegreen"), lwd=c(NA,1,NA))
dev.off()


