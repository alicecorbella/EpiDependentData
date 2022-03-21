rm(list=ls())
setwd("~/EpiDependentData/sec5/")
# install.packages("extrafont")
library("xtable")
library(extrafont)
library(Rcpp)
library(mvtnorm)
library(denstrip)
# font_import()
loadfonts()


myblue   <- rgb(0/255,128/255, 128/255, 0.95)
mygray   <- rgb(255/355 ,  192/255,  203/255, 1)
myred    <- rgb(213/355 ,    0/255,  50/255, 1)
mypurple <- rgb(119/355 ,   37/255,  131/255, 1)
myyellow <- rgb(225/255 ,  215/255,  0/255, 1)


my2blue   <- rgb(0/255,128/255, 128/255, 0.65)
my2red    <- rgb(213/355 ,    0/255,  50/255, 0.65)

load("~/EpiDependentData/sec5/results/Movable.RData")
sourceCpp("~/EpiDependentData/sec5/CppAll.cpp")
source("~/EpiDependentData/sec5/Rfun.R")

head(th.bs)


# TRANSMISSION  -----------------------------------------------------------

# fig 6.3.2: transmission parameters ---------------------------------------

pdf("fig6_3_1.pdf",  height = 4, width=5)
par( mar=c(2,4,2,1),  cex.main=2,las=1) 

# beta
xts <- seq(0,4, length.out = 10000)
plot(0,0, col=0, main=expression(beta), ylab="Density", 
     xlim=c(0.45,0.9), ylim=c(0,8))
polygon(c(xts, rev(xts)), c(dunif(xts, 0, 4), rep(0, length(xts))), 
        col=my2red, border = mygray, lwd=2)
polygon(density (th.bs[,1]), 
        col=my2blue, border = mygray, lwd=2)

# pi
xts <- seq(0,1, length.out = 10000)
plot(0,0, col=0, main=expression(pi), ylab="Density", 
     xlim=c(0.25,0.55), ylim=c(0,11))
polygon(c(xts, rev(xts)), c(dbeta(xts, shape1 =37.5 , shape2=62.5), rep(0, length(xts))), 
        col=my2red, border = mygray, lwd=2)
polygon(density (th.bs[,2]), 
        col=my2blue, border = mygray, lwd=2)

# iota
xts <- seq(0,0.1, length.out = 10000)
plot(0,0, col=0, main=expression(iota), ylab="Density", 
     xlim=c(0,0.0004), ylim=c(0,8000))
polygon(c(xts, rev(xts)), c(dunif(xts, 0, 0.1), rep(0, length(xts))), 
        col=my2red, border = NA, lwd=2)
polygon(density (th.bs[,3]), 
        col=my2blue, border = mygray, lwd=2)

# d_l
xts <- seq(0,10, length.out = 10000)
plot(0,0, col=0, main=expression(d[L]), ylab="Density", 
     xlim=c(0,8), ylim=c(0,0.5))
polygon(c(xts, rev(xts)), c(dlnorm(xts, meanlog = log(2), sd=0.5), rep(0, length(xts))), 
        col=my2red, border = mygray, lwd=2)
polygon(density (th.bs[,4]), 
        col=my2blue, border = mygray, lwd=2)

# kappa
xts <- seq(-1,4, length.out = 10000)
plot(0,0, col=0, main=expression(kappa), ylab="Density", 
     xlim=c(-0.1,0.9), ylim=c(0,5))
abline(v=0, col=mygray)
polygon(c(xts, rev(xts)), c(dnorm((log(xts+1)), 0, sd=1), rep(0, length(xts))), 
        col=my2red, border = mygray, lwd=2)
polygon(density (th.bs[,5]), 
        col=my2blue, border = mygray, lwd=2)
dev.off()


# fig 6.3.2: transmission summaries ---------------------------------------

pdf("fig6_3_2.pdf",  height = 4, width=5)
par( mar=c(2,4,2,1),  cex.main=2,las=1) 
# RO when kappa is not effective
R0  <- th.bs[,1]*thetaconst$dI
R0p <- runif(10000, 0,4)*thetaconst$dI
plot(0,0, col=0, main=expression(R[0]), ylab="Density", 
     xlim=c(1,4), ylim=c(0,2.5))
polygon(density(R0p), 
        col=my2red, border = mygray, lwd=2)
polygon(density (R0), 
        col=my2blue, border = mygray, lwd=2)

# Re when kappa is not effective
Re <- th.bs[,1]*thetaconst$dI*(1-th.bs[,2]-th.bs[,3])
Rep <- runif(10000, 0,4)*thetaconst$dI*(1-rbeta(10000, shape1 =37.5 , shape2=62.5)-runif(10000, 0,0.1))
plot(0,0, col=0, main=expression(R[e]), ylab="Density", 
     xlim=c(1,2.5), ylim=c(0,7))
polygon(density(Rep), 
        col=my2red, border = mygray, lwd=2)
polygon(density (Re), 
        col=my2blue, border = mygray, lwd=2)

# RO when kappa is effective
R0K  <- th.bs[,1]*thetaconst$dI*(1+th.bs[,5])
R0pK <- runif(10000, 0,4)*thetaconst$dI*(1+exp(rnorm(10000, 0, sd=1)-1))*(exp(rnorm(10000, mean=1, sd=1)))
# hist(R0, main=expression(R[0]), 
#      col=rgb(0/255,128/255, 128/255, 1), border="gray", lwd=2, freq=F)
plot(0,0, col=0, main=expression(kappa*R[0]), ylab="Density", 
     xlim=c(1,4), ylim=c(0,2.5))
polygon(density(R0pK), 
        col=my2red, border = mygray, lwd=2)
polygon(density (R0K), 
        col=my2blue, border = mygray, lwd=2)

# Re when kappa is not effective
ReK <- th.bs[,1]*thetaconst$dI*(1-th.bs[,2]-th.bs[,3])*(1+th.bs[,5])
RepK <- runif(10000, 0,4)*thetaconst$dI*(1-rbeta(10000, shape1 =37.5 , shape2=62.5)-runif(10000, 
        0,0.1))*(exp(rnorm(10000, mean=1, sd=1)))
plot(0,0, col=0, main=expression(kappa*R[e]), ylab="Density", 
     xlim=c(1,2.5), ylim=c(0,7))
polygon(density(RepK), 
        col=my2red, border = mygray, lwd=2)
polygon(density (ReK), 
        col=my2blue, border = mygray, lwd=2)

dev.off()

# fig 6.3.3:  number of new infections ---------------------------------------

xiFuP <- matrix(NA, ncol=(thetaconst$TMbg+thetaconst$TMnd), nrow=(dim(th.bs)[1]))

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

xiLB <- apply(xiFuP, MARGIN = 2, quantile, prob=0.025)
xiUB <- apply(xiFuP, MARGIN = 2, quantile, prob=0.975)
xiME <- apply(xiFuP, MARGIN = 2, quantile, prob=0.5)


pdf("fig6_3_3.pdf",  height = 5, width=8)
par( mar=c(4,5,1,1),  cex.main=2,las=1) 
plot(xiUB, col=0, xlab = "Time (in days)", ylab = "")
title(ylab="Number of daily new infections", line=4)
abline(v=which(thetaconst$vecHtb==1)/4, lwd=2, col=rgb(0,0,0,0.1))
polygon(x=c(1:length(xiLB),length(xiLB):1), y=c(xiLB, rev(xiUB)), 
        border=NA, col=rgb(0/255,128/255, 128/255, 0.45))
for (s in 1:20){
  lines(xiFuP[sample(1:dim(xiFuP)[1], 1),], col=rgb(0/255,128/255, 128/255, 0.25), lwd=1)
}
lines(xiME, col=myred, lwd=2)
legend(0, 300000, c("Median", "95% CrIs"), col=c(myred, rgb(0/255,128/255, 128/255, 0.45)), 
       lwd=c(2, 10))
dev.off()

#  SEVERITY PARAMTERS -----------------------------------------------------

# figure 6.3.4. severity parameters ---------------------------------------

pdf("fig6_3_4.pdf",  height = 4, width=5)
par( mar=c(2,4,2,1),  cex.main=2,las=1) 

# iotaF
xts <- seq(0,1, length.out = 10000)
plot(0,0, col=0, main=expression(iota^F), ylab="Density", 
     xlim=c(0.,0.01), ylim=c(0,800))
polygon(c(xts, rev(xts)), c(dunif(xts, 0, 1), rep(0, length(xts))), 
        col=my2red, border = my2red, lwd=2)
polygon(density (th.bs[,12]), 
        col=my2blue, border = mygray, lwd=2)

# epsF
xts <- seq(0,5000, length.out = 10000)
plot(0,0, col=0, main=expression(epsilon^F), ylab="Density", 
     xlim=c(0.,1500), ylim=c(0,0.004))
polygon(c(xts, rev(xts)), c(dunif(xts, 0, 5000), rep(0, length(xts))), 
        col=my2red, border = mygray, lwd=2)
polygon(density (th.bs[,13]), 
        col=my2blue, border = mygray, lwd=2)

# thetaH
xts <- seq(0,1, length.out = 10000)
plot(0,0, col=0, main=expression(phantom(0)^0*theta^{H}), ylab="Density", 
     xlim=c(0,0.01), ylim=c(0,700))
polygon(c(xts, rev(xts)), c(dunif(xts, 0, 1), rep(0, length(xts))), 
        col=my2red, border = my2red, lwd=2)
polygon(density (th.bs[,6]), 
        col=my2blue, border = mygray, lwd=2)
median(th.bs[,6])
# thetaIC
xts <- seq(0,1, length.out = 10000)
plot(0,0, col=0, main=expression(phantom(0)^H*theta^{IC}), ylab="Density", 
     xlim=c(0,0.15), ylim=c(0,80))
polygon(c(xts, rev(xts)), c(dunif(xts, 0, 1), rep(0, length(xts))), 
        col=my2red, border = mygray, lwd=2)
polygon(density (th.bs[,7]), 
        col=my2blue, border = mygray, lwd=2)
median(th.bs[,7])
dev.off()




# fig6.3.5 theta0F -----------------------------------------------------------------

epsPrior <- runif(10000, 0, 5000) 
iotPrior <- runif(10000, 0, 1) 
ThGPrior <- rgamma(10000, shape= iotPrior*epsPrior, rate = epsPrior)
ThGPoste <- rgamma(dim(th.bs)[1], shape= th.bs[,12]*th.bs[,13], rate = th.bs[,13])

pdf("fig6_3_5.pdf",  height = 4, width=5)
par( mar=c(2,4,2,1),  cex.main=2,las=1) 
plot(0,0, col=0, main=expression(phantom(0)^0*theta^{F}), ylab="Density", 
     xlim=c(0,0.015), ylim=c(0,280))
polygon(density (ThGPrior), 
        col=my2red, border = my2red, lwd=2)
polygon(density (ThGPoste), 
        col=my2blue, border = mygray, lwd=2)

plot(0,0, col=0, main=expression(phantom(0)^0*theta^{F}), ylab="Density", 
     xlim=c(0,1), ylim=c(0,3))
polygon(density (ThGPrior), 
        col=my2red, border = mygray, lwd=2)
dev.off()


#  BACKGROUND PARAMTERS ---------------------------------------------------


# figure 6.3.6. backround parameters ---------------------------------------

pdf("fig6_3_6.pdf",  height = 4, width=5)
par( mar=c(2,4,2,1),  cex.main=2,las=1) 

xts <- seq(-20,20, length.out = 10000)
plot(0,0, col=0, main=expression(nu[1]), ylab="Density", 
     xlim=c(6,7.2), ylim=c(0,9))
polygon(c(xts, rev(xts)), c(dnorm(xts, mean = 6.6, sd=0.17), rep(0, length(xts))), 
        col=my2red, border = mygray, lwd=2)
polygon(density (th.bs[,9]), 
        col=my2blue, border = mygray, lwd=2)

xts <- seq(-20,20, length.out = 10000)
plot(0,0, col=0, main=expression(nu[2]), ylab="Density", 
     xlim=c(-0.55,0.2), ylim=c(0,9))
polygon(c(xts, rev(xts)), c(dnorm(xts, mean = -0.2, sd=0.11), rep(0, length(xts))), 
        col=my2red, border = mygray, lwd=2)
polygon(density (th.bs[,10]), 
        col=my2blue, border = mygray, lwd=2)

xts <- seq(-20,20, length.out = 10000)
plot(0,0, col=0, main=expression(nu[3]), ylab="Density", 
     xlim=c(0.5,1.35), ylim=c(0,8))
polygon(c(xts, rev(xts)), c(dnorm(xts, mean = 0.99, sd=0.08), rep(0, length(xts))), 
        col=my2red, border = mygray, lwd=2)
polygon(density (th.bs[,11]), 
        col=my2blue, border = mygray, lwd=2)
dev.off()

# figure 6.3.7. backround ILI ---------------------------------------

guP <- t(apply(th.bs[,9:11], MARGIN = 1, nuTOg))
dim(guP)

MEguP <- apply(guP, MARGIN = 2, median)
LOguP <- apply(guP, MARGIN = 2, quantile, prob=0.025)
UPguP <- apply(guP, MARGIN = 2, quantile, prob=0.975)


pdf("fig6_3_7.pdf",  height = 5, width=8)
par( mar=c(4,5,1,1),  cex.main=2,las=1) 
plot(UPguP, type="l", ylim=c(0,600), col=0, xlab = "Time (in days)", 
     ylab = "Background ILI cases")
polygon(x=c(1:length(LOguP), length(LOguP):1), 
        y=c(nuTOg(c(qnorm(0.025, 6.6, 0.17), qnorm(0.025, -.2, 0.11), qnorm(0.025, 0.99, 0.08))),
            rev(nuTOg(c(qnorm(0.975, 6.6, 0.17), qnorm(0.975, -.2, 0.11), qnorm(0.975, 0.99, 0.08))))),
        col=rgb(213/355, 0/255,  50/255, 0.25), border=mygray)
lines(nuTOg(c(6.6, -0.2,  0.99)), col=mygray, lwd=3)
lines(nuTOg(c(6.6, -0.2,  0.99)), col=myred, lwd=2)
polygon(x=c(1:length(LOguP), length(LOguP):1), y=c(LOguP,rev(UPguP)),
        col=rgb(0/255,128/255, 128/255, 0.85), border=mygray)
lines(MEguP, col=1, lwd=3)
lines(MEguP, col=myblue, lwd=2)
legend(150, 600, c(expression(paste("prior ", b[u])), expression(paste("posterior ", b[u]))), 
       col=c(rgb(213/355, 0/255,  50/255, 0.25),
             rgb(0/255,128/255, 128/255, 0.85)), lwd=10)
dev.off()



# DAY of THE WEEK EFFECT --------------------------------------------------


# figure 6.3.8. DOW priorpost ---------------------------------------

pOmega <- exp(rmvnorm(10000, mean=rep(0,6), sigma=0.16*matrix(
  c(1, -1/6, -1/6, -1/6, -1/6, -1/6,
    -1/6, 1, -1/6, -1/6, -1/6, -1/6,
    -1/6, -1/6, 1, -1/6, -1/6, -1/6,
    -1/6, -1/6, -1/6, 1, -1/6, -1/6,
    -1/6, -1/6, -1/6, -1/6, 1, -1/6,
    -1/6, -1/6, -1/6, -1/6, -1/6, 1), byrow = T, ncol=6))
)

pdf("fig6_3_8.pdf",  height = 4, width=5)
par( mar=c(2,4,2,1),  cex.main=2,las=1) 
plot(0,0, col=0, main=expression(omega[1]), ylab="Density", 
     xlim=c(0,5), ylim=c(0,2))
polygon(density(pOmega[,1]), 
        col=my2red, border = mygray, lwd=2)
polygon(density (th.bs[,14]), 
        col=my2blue, border = mygray, lwd=2)
plot(0,0, col=0, main=expression(omega[2]), ylab="Density", 
     xlim=c(0,5), ylim=c(0,2.5))
polygon(density(pOmega[,2]), 
        col=my2red, border = mygray, lwd=2)
polygon(density (th.bs[,15]), 
        col=my2blue, border = mygray, lwd=2)
plot(0,0, col=0, main=expression(omega[3]), ylab="Density", 
     xlim=c(0,5), ylim=c(0,2.5))
polygon(density(pOmega[,3]), 
        col=my2red, border = mygray, lwd=2)
polygon(density (th.bs[,16]), 
        col=my2blue, border = mygray, lwd=2)
plot(0,0, col=0, main=expression(omega[5]), ylab="Density", 
     xlim=c(0,5), ylim=c(0,2.5))
polygon(density(pOmega[,4]), 
        col=my2red, border = mygray, lwd=2)
polygon(density (th.bs[,17]), 
        col=my2blue, border = mygray, lwd=2)
plot(0,0, col=0, main=expression(omega[6]), ylab="Density", 
     xlim=c(0,.2), ylim=c(0,80))
polygon(density(pOmega[,5]), 
        col=my2red, border = my2red, lwd=2)
polygon(density (th.bs[,18]), 
        col=my2blue, border = mygray, lwd=2)
plot(0,0, col=0, main=expression(omega[7]), ylab="Density", 
     xlim=c(0,.2), ylim=c(0,80))
polygon(density(pOmega[,6]), 
        col=my2red, border = my2red, lwd=2)
polygon(density (th.bs[,19]), 
        col=my2blue, border = mygray, lwd=2)
dev.off()


# figure 6.3.9. DOW thetarate ---------------------------------------
omega4 <- 1/(apply(th.bs[,14:19], MARGIN = 1, prod))


pdf("fig6_3_9.pdf",  height = 5, width=6)
plot(0,0, col=0, xlim=c(0.5,7.5), ylim=c(0,0.05),
     xaxt="n", xlab="week day", ylab="GP consultation rate")
axis(1, at=1:7, c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"))
denstrip(ThGPoste*th.bs[, 14], horiz=F, at=1, colmax = myblue, 
         width=0.9)
denstrip(ThGPoste*th.bs[, 15], horiz=F, at=2, colmax = myblue, 
         width=0.9)
denstrip(ThGPoste*th.bs[, 16], horiz=F, at=3, colmax = myblue, 
         width=0.9)
denstrip(ThGPoste*omega4, horiz=F, at=4, colmax = myblue, 
         width=0.9)
denstrip(ThGPoste*th.bs[, 17], horiz=F, at=5, colmax = myblue, 
         width=0.9)
denstrip(ThGPoste*th.bs[, 18], horiz=F, at=6, colmax = myblue, 
         width=0.9)
denstrip(ThGPoste*th.bs[, 19], horiz=F, at=7, colmax = myblue, 
         width=0.9)
box()
dev.off()

# DETECTION PARAMETER -----------------------------------------------------


# figure 6.3.10. epsH ---------------------------------------

pdf("fig6_3_10.pdf",  height = 4, width=5)
par( mar=c(2,4,2,1),  cex.main=2,las=1) 
xts <- seq(0,100, length.out = 10000)
plot(0,0, col=0, main=expression(epsilon^H), ylab="Density", 
     xlim=c(0,100), ylim=c(0,0.07))
polygon(c(xts, rev(xts)), c(dunif(xts, 0, 100), rep(0, length(xts))), 
        col=my2red, border = mygray, lwd=2)
polygon(density (th.bs[,8]), 
        col=my2blue, border = mygray, lwd=2)
dev.off()

zetaH <- matrix(NA, ncol=length(thetaconst$detH), nrow = length(th.bs[,8]))
for (i in 1:length(thetaconst$detH)){
  zetaH[,i] <- rbeta(length(th.bs[,8]), shape1=(thetaconst$detH[i]/(1-thetaconst$detH[i]))*th.bs[,8],
       shape2 = th.bs[,8])
}


zetaLB <- apply(zetaH, MARGIN = 2, quantile, prob=0.025)
zetaUB <- apply(zetaH, MARGIN = 2, quantile, prob=0.975)
zetaME <- apply(zetaH, MARGIN = 2, quantile, prob=0.5)


pdf("fig6_3_11.pdf",  height = 5, width=8)
par( mar=c(4,5,1,1),  cex.main=2,las=1) 
plot(zetaUB, ylim=c(0,0.5), col=0, xlab = "Time (in weeks)", ylab = "")
title(ylab="Hospitalization detection", line=4)
polygon(x=c(1:length(zetaLB),length(zetaLB):1), y=c(zetaLB, rev(zetaUB)), 
        border=NA, col=rgb(0/255,128/255, 128/255, 0.45))
for (s in 1:20){
  lines(zetaH[sample(1:dim(zetaH)[1], 1),], col=rgb(0/255,128/255, 128/255, 0.25), lwd=1)
}
lines(zetaME, col=myred, lwd=2)
legend(0, 0.5, c("Median", "95% CrIs"), col=c(myred, rgb(0/255,128/255, 128/255, 0.45)), 
       lwd=c(2, 10))
dev.off()



# TABLE -------------------------------------------------------------------

MEthPO <- round(apply(th.bs, MARGIN = 2, median), 4)
LBthPO <- round(apply(th.bs, MARGIN = 2, quantile, prob=0.025), 4)
UBthPO <- round(apply(th.bs, MARGIN = 2, quantile, prob=0.975), 4)

fortable <- paste(MEthPO, " (", LBthPO, "-", UBthPO, ")", sep="")


apply(th.bs, MARGIN = 2, median)
apply(th.bs, MARGIN = 2, quantile, prob=0.025)
apply(th.bs, MARGIN = 2, quantile, prob=0.975)


xtable(as.matrix(cbind(vecNAMES, fortable)))



# CHAINS ------------------------------------------------------------------

head(th.bs)
dim(th.bs)
vecNAMES <- c(expression(beta),expression(pi),expression(iota),
              expression(d[L]),expression(kappa),
              expression(phantom(0)^0*theta^{H}),expression(phantom(0)^H*theta^{IC}), expression(epsilon^H),
              expression(nu[1]),expression(nu[2]),expression(nu[3]), 
              expression(iota^F),expression(epsilon^F),
              expression(omega[1],omega[2],omega[3],omega[5],omega[6],omega[7]))
pdf("figChains.pdf", height=3, width=6)
par(mfrow=c(1,1), mar=c(4,2,2,1))

for (p in 1:19){
  plot(th.bs[1:5000, p], type="l", col=myblue,
       main=vecNAMES[p], ylab="", xlab = "n samples")
  lines(th.bs[5001:10000, p], col=myyellow)
  lines(th.bs[10001:15000, p], col=myred)
}
dev.off()
