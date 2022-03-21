rm(list=ls())

setwd("~/EpiDependentData/sec4/")
loadwd <- "~/EpiDependentData/sec4/results"
library(fields)
library(extrafont)
# font_import()
loadfonts()
# setting of the colors
Csus <- rgb(0/255   ,  128/255, 128/255, 1)
Cexp <- rgb(255/355 ,  192/255,  203/255, 1)
Cinf <- rgb(213/355 ,    0/255,  50/255, 1)
Crec <- rgb(119/355 ,   37/255,  131/255, 1)
Cgols <- rgb(255/255 ,   215/255,  0, 1)

CsusT <- rgb(0/255   ,  128/255, 128/255, 0.5)
CexpT <- rgb(255/355 ,  192/255,  203/255, 0.5)
CinfT <- rgb(213/355 ,    0/255,  50/255, 0.5)
CrecT <- rgb(119/355 ,   37/255,  131/255, 0.5)
colpal <- list(Csus, Cexp, Cinf, Crec, Cgols)

t.test0.pv <- function(x1n){
  n     <- length(x1n)
  hatmu <- mean(x1n)
  hatsd <- sd(x1n)
  tobs  <- hatmu/(hatsd/sqrt(n))
  pval  <- 1 - pt(tobs, df = n-1)#P(T>t.obs)
}

# sub SECTION 1 -----------------------------------------------------------

load(paste(loadwd, "OUTforTHESIS/bs3par.Robj", sep=""))

# SMALL

SMALLdep <- comp1bs3par$SMALL$dep
SMALLind <- comp1bs3par$SMALL$ind

pdf("ss_sce4.pdf",  height =4.5, width=7)
set.seed(291218)
ndx <- sample(1:500, 5)
lims <- rbind(c(0,8), c(0,8), c(0,10000))
par(mfrow=c(1, 1), mar=c(2,4,3,1))
for (p in 1:3){
  plot(density(SMALLdep[, p, ndx[1]], from=0), col=0, ylim=lims[p, ], las=1, 
       main=list(expression(paste(beta, " : SMALL dependence")),
                 expression(paste(pi, " : SMALL dependence")),
                 expression(paste(iota, " : SMALL dependence")))[[p]])
  for (i in 1:5){
    lines(density(SMALLdep[, p, ndx[i]], from=0), col=colpal[[i]], lwd=3, lty=1)
    lines(density(SMALLind[, p, ndx[i]], from=0), col=colpal[[i]], lwd=3, lty=2)
  }
  legend(0.002, lims[p, 2], c("joint dep", "miss ind") ,
         lty=c(1,2), lwd=3, bty = "n")
}

CrIRdep <- t(apply(SMALLdep, MARGIN = c(2,3), quantile, prob=0.975))-
          t(apply(SMALLdep, MARGIN = c(2,3), quantile, prob=0.025))
CrIRind <- t(apply(SMALLind, MARGIN = c(2,3), quantile, prob=0.975))-
  t(apply(SMALLind, MARGIN = c(2,3), quantile, prob=0.025))

VARdep <- t(apply(SMALLdep, MARGIN = c(2,3), var))
VARind <- t(apply(SMALLind, MARGIN = c(2,3), var))

par(mfrow=c(1, 2))
for (p in 1:3){
  par(mar=c(4, 4, 3,0.5))
  hist((VARdep[,p]- VARind[,p]), freq=F, breaks=10, ylab="Density", xlab="Variance difference", las=1,adj=1,
       main=list(expression(paste(beta, " : SMALL dependence")),
                 expression(paste(pi, " : SMALL dependence")),
                 expression(paste(iota, " : SMALL dependence")))[[p]])
  abline(v=0, lwd=3, lty=3)
  par(mar=c(4, 0.5, 3,4))
  hist((CrIRdep[,p]- CrIRind[,p]),freq=F, breaks=10, lwd=3, ylab=" ", 
       xlab="CrI Range difference",yaxt="n", adj=0,
       main=expression(paste("pairwise variability difference")))
  abline(v=0, lwd=3, lty=3)
}


VARdiff <- VARdep-VARind
apply(VARdiff<0, MARGIN = 2, sum)/500
for (i in 1:3){
  print(t.test0.pv(VARdiff[, i]))
}
# BIG

BIGdep <- comp1bs3par$BIG$dep
BIGind <- comp1bs3par$BIG$ind

set.seed(291218)
ndx <- sample(1:500, 5)
lims <- 4*rbind(c(0,8), c(0,8), c(0,10000))
limsX <- rbind(c(0.5,.75), c(0.17,.4), c(0,0.0002))
par(mfrow=c(1, 1), mar=c(2,4,3,1))
for (p in 1:3){
  plot(density(BIGdep[, p, ndx[1]]), col=0, ylim=lims[p, ],xlim=limsX[p,], las=1,
       main=list(expression(paste(beta, " : BIG dependence")),
                 expression(paste(pi, " : BIG dependence")),
                 expression(paste(iota, " : BIG dependence")))[[p]])
  for (i in 1:5){
    lines(density(BIGdep[, p, ndx[i]], from=0), col=colpal[[i]], lwd=3, lty=1)
    lines(density(BIGind[, p, ndx[i]], from=0), col=colpal[[i]], lwd=3, lty=2)
  }
  legend(limsX[p, 1], lims[p, 2], c("joint dep", "miss ind") ,
         lty=c(1,2), lwd=3, bty = "n")
}

CrIRdep <- t(apply(BIGdep, MARGIN = c(2,3), quantile, prob=0.975))-
  t(apply(BIGdep, MARGIN = c(2,3), quantile, prob=0.025))
CrIRind <- t(apply(BIGind, MARGIN = c(2,3), quantile, prob=0.975))-
  t(apply(BIGind, MARGIN = c(2,3), quantile, prob=0.025))

VARdep <- t(apply(BIGdep, MARGIN = c(2,3), var))
VARind <- t(apply(BIGind, MARGIN = c(2,3), var))

par(mfrow=c(1, 2))
for (p in 1:3){
  par(mar=c(4, 4, 3,0.5))
  hist((VARdep[,p]- VARind[,p]), freq=F, breaks=10,  ylab="Density", xlab="Variance difference", las=1,
       main=list(expression(paste(beta, " : BIG dependence")),
                 expression(paste(pi, " : BIG dependence")),
                 expression(paste(iota, " : BIG dependence")))[[p]], adj=1)
  abline(v=0, lwd=3, lty=3)
  par(mar=c(4, 0.5, 3,4))
  hist((CrIRdep[,p]- CrIRind[,p]), freq=F, breaks=10,  ylab=" ", 
       xlab="CrI Range difference",yaxt="n", adj=0,
       main=expression(paste("pairwise variability difference")))
  abline(v=0, lwd=1, lty=3)
}

VARdiff <- VARdep-VARind
apply(VARdiff<0, MARGIN = 2, sum)/500
# dev.off()

rm(list=setdiff(ls(), c("Cexp","CexpT","Cgols","Cinf",
                        "CinfT","colpal","Crec","CrecT",
                        "Csus","CsusT", "loadwd")))
# sub SECTION 2 -----------------------------------------------------------

load(paste(loadwd,"OUTforTHESIS/bs5par.Robj", sep=""))

# SMALL

SMALLdep <- comp2bs5par$SMALL$dep
SMALLind <- comp2bs5par$SMALL$ind

# pdf("fig5_7.pdf", #family="cmr10", 
#     height =4.5, width=7)
set.seed(291218)
ndx <- sample(1:500, 5)
lims <- rbind(c(0,6), c(0,6), c(0,10000), c(0,22), c(0,9) )
limsX <- rbind(c(0,1.3), c(0,.8), c(0,0.001), c(0,0.4), c(0, 0.7))
par(mfrow=c(1, 1), mar=c(2,4,3,1))
for (p in 1:5){
  plot(density(SMALLdep[, p, ndx[1]], from=0), col=0, ylim=lims[p, ], las=1, xlim=limsX[p,],
       main=list(expression(paste(beta, " : SMALL dependence")),
                 expression(paste(pi, " : SMALL dependence")),
                 expression(paste(iota, " : SMALL dependence")),
                 expression(paste(.^0, theta^H, " : SMALL dependence")),
                 expression(paste(.^H, theta^{IC}," : SMALL dependence")) )[[p]])
  for (i in 1:5){
    lines(density(SMALLdep[, p, ndx[i]], from=0), col=colpal[[i]], lwd=3, lty=1)
    lines(density(SMALLind[, p, ndx[i]], from=0), col=colpal[[i]], lwd=3, lty=2)
  }
  legend(0.002, lims[p, 2], c("joint dep", "miss ind") ,
         lty=c(1,2), lwd=3, bty = "n")
}

CrIRdep <- t(apply(SMALLdep, MARGIN = c(2,3), quantile, prob=0.975))-
  t(apply(SMALLdep, MARGIN = c(2,3), quantile, prob=0.025))
CrIRind <- t(apply(SMALLind, MARGIN = c(2,3), quantile, prob=0.975))-
  t(apply(SMALLind, MARGIN = c(2,3), quantile, prob=0.025))

VARdep <- t(apply(SMALLdep, MARGIN = c(2,3), var))
VARind <- t(apply(SMALLind, MARGIN = c(2,3), var))

par(mfrow=c(1, 2))
for (p in 1:5){
  par(mar=c(4, 4, 3,0.5))
  hist((VARdep[,p]- VARind[,p]), freq=F, breaks=10,  ylab="Density", xlab="Variance difference", las=1,
       main=list(expression(paste(beta, " : SMALL dependence")),
                 expression(paste(pi, " : SMALL dependence")),
                 expression(paste(iota, " : SMALL dependence")),
                 expression(paste(.^0, theta^H, " : SMALL dependence")),
                 expression(paste(.^H, theta^{IC}," : SMALL dependence")) )[[p]], 
       adj=1)
  abline(v=0, lwd=3, lty=3)
  par(mar=c(4, 0.5, 3,4))
  hist((CrIRdep[,p]- CrIRind[,p]), freq=F, breaks=10,  ylab=" ", 
       xlab="CrI Range difference",yaxt="n", 
       main=expression(paste("pairwise variability difference")), 
       adj=0)
  abline(v=0, lwd=3, lty=3)
}

VARdiff <- VARdep-VARind
apply(VARdiff<0, MARGIN = 2, sum)/500

# BIG

BIGdep <- comp2bs5par$BIG$dep
BIGind <- comp2bs5par$BIG$ind

set.seed(291218)
ndx <- sample(1:500, 5)
lims <- 2*rbind(c(0,4), c(0,4), c(0,13000), c(0,5.5), c(0,9) )
limsX <- rbind(c(0.3,1.3), c(0,.8), c(0,0.0002), c(0.3,0.8), c(0.75, 1))
# limsX <- rbind(c(0,.75), c(0.17,.4), c(0,0.0002), c(0,0.5), c(0.7, 1))
par(mfrow=c(1, 1), mar=c(2,4,3,1))
for (p in 1:5){
  plot(density(BIGdep[, p, ndx[1]]), col=0, ylim=lims[p, ],xlim=limsX[p,], las=1,
       main=list(expression(paste(beta, " : BIG dependence")),
                 expression(paste(pi, " : BIG dependence")),
                 expression(paste(iota, " : BIG dependence")),
                 expression(paste(.^0, theta^H, " : BIG dependence")),
                 expression(paste(.^H, theta^{IC}," : BIG dependence")) )[[p]])
  for (i in 1:5){
    lines(density(BIGdep[, p, ndx[i]], from=0), col=colpal[[i]], lwd=3, lty=1)
    lines(density(BIGind[, p, ndx[i]], from=0), col=colpal[[i]], lwd=3, lty=2)
  }
  legend(limsX[p, 1], lims[p, 2], c("joint dep", "miss ind") ,
         lty=c(1,2), lwd=3, bty = "n")
}

CrIRdep <- t(apply(BIGdep, MARGIN = c(2,3), quantile, prob=0.975))-
  t(apply(BIGdep, MARGIN = c(2,3), quantile, prob=0.025))
CrIRind <- t(apply(BIGind, MARGIN = c(2,3), quantile, prob=0.975))-
  t(apply(BIGind, MARGIN = c(2,3), quantile, prob=0.025))

VARdep <- t(apply(BIGdep, MARGIN = c(2,3), var))
VARind <- t(apply(BIGind, MARGIN = c(2,3), var))

par(mfrow=c(1, 2))
for (p in 1:5){
  par(mar=c(4, 4, 3,0.5))
  hist((VARdep[,p]- VARind[,p]), freq=F, breaks=10, ylab="Density", xlab="Variance difference", las=1, adj=1,
       main=list(expression(paste(beta, " : BIG dependence")),
                 expression(paste(pi, " : BIG dependence")),
                 expression(paste(iota, " : BIG dependence")),
                 expression(paste(.^0, theta^H, " : BIG dependence")),
                 expression(paste(.^H, theta^{IC}," : BIG dependence")) )[[p]])
  abline(v=0, lwd=3, lty=3)
  par(mar=c(4, 0.5, 3,4))
  hist((CrIRdep[,p]- CrIRind[,p]), freq=F, breaks=10,  ylab=" ", 
       xlab="CrI Range difference", 
       main=expression(paste("pairwise variability difference")), adj=0)
  abline(v=0, lwd=3, lty=3)
}


VARdiff <- VARdep-VARind
apply(VARdiff<0, MARGIN = 2, sum)/500
# dev.off()


# sub SECTION 3 -----------------------------------------------------------
rm(list=setdiff(ls(), c("Cexp","CexpT","Cgols","Cinf",
                        "CinfT","colpal","Crec","CrecT",
                        "Csus","CsusT", "loadwd")))
load(paste(loadwd, "OUTforTHESIS/SVRT.Robj", sep=""))

# ThHO

ThHOdep <- comp3SVRT$ThHO$dep
ThHOind <- comp3SVRT$ThHO$ind

# pdf("fig5_8.pdf", #family="cmr10", 
#     height =4.5, width=7)
set.seed(291218)
ndx <- sample(1:500, 5)
lims <- rbind(c(0,6), c(0,6), c(0,15000), c(0,8), c(0,20) )
limsX <- rbind(c(0.4,1.3), c(0,.8), c(0,0.0005), c(0.3,0.8), c(0, 0.25))
par(mfrow=c(1, 1), mar=c(2,4,3,1))
for (p in 1:5){
  plot(density(ThHOdep[, p, ndx[1]], from=0), col=0, ylim=lims[p, ], las=1, xlim=limsX[p,],
       main=list(expression(paste(beta, " : large ThHO")),
                 expression(paste(pi, " : large ThHO")),
                 expression(paste(iota, " : large ThHO")),
                 expression(paste(.^0, theta^H, " : large ThHO")),
                 expression(paste(.^H, theta^{IC}," : large ThHO")) )[[p]])
  for (i in 1:5){
    lines(density(ThHOdep[, p, ndx[i]], from=0), col=colpal[[i]], lwd=3, lty=1)
    lines(density(ThHOind[, p, ndx[i]], from=0), col=colpal[[i]], lwd=3, lty=2)
  }
  legend(0.002, lims[p, 2], c("joint dep", "miss ind") ,
         lty=c(1,2), lwd=3, bty = "n")
}

CrIRdep <- t(apply(ThHOdep, MARGIN = c(2,3), quantile, prob=0.975))-
  t(apply(ThHOdep, MARGIN = c(2,3), quantile, prob=0.025))
CrIRind <- t(apply(ThHOind, MARGIN = c(2,3), quantile, prob=0.975))-
  t(apply(ThHOind, MARGIN = c(2,3), quantile, prob=0.025))

VARdep <- t(apply(ThHOdep, MARGIN = c(2,3), var))
VARind <- t(apply(ThHOind, MARGIN = c(2,3), var))

par(mfrow=c(1, 2))
for (p in 1:5){
  par(mar=c(4, 4, 3,0.5))
  hist((VARdep[,p]- VARind[,p]), freq=F, breaks=10,  ylab="Density", xlab="Variance difference", las=1, adj=1,
       main=list(expression(paste(beta, " : large ThHO")),
                 expression(paste(pi, " : large ThHO")),
                 expression(paste(iota, " : large ThHO")),
                 expression(paste(.^0, theta^H, " : large ThHO")),
                 expression(paste(.^H, theta^{IC}," : large ThHO")) )[[p]])
  abline(v=0, lwd=3, lty=3)
  par(mar=c(4, 0.5, 3,4))
  hist((CrIRdep[,p]- CrIRind[,p]), freq=F, breaks=10,  ylab=" ", 
       xlab="CrI Range difference",yaxt="n", 
       main=expression(paste("pairwise variability difference")), adj=0)
  abline(v=0, lwd=3, lty=3)
}

VARdiff <- VARdep-VARind
apply(VARdiff<0, MARGIN = 2, sum)/500

# ThIC

ThICdep <- comp3SVRT$ThIC$dep
ThICind <- comp3SVRT$ThIC$ind

set.seed(291218)
ndx <- sample(1:500, 5)
lims <- rbind(c(0,6), c(0,6), c(0,23000), c(0,14.5), c(0,8) )
limsX <- rbind(c(0.3,1.3), c(0,.8), c(0,0.0006), c(0,0.3), c(0.35, 1))
# limsX <- rbind(c(0,.75), c(0.17,.4), c(0,0.0002), c(0,0.5), c(0.7, 1))
par(mfrow=c(1, 1), mar=c(2,4,3,1))
for (p in 1:5){
  plot(density(ThICdep[, p, ndx[1]]), col=0, ylim=lims[p, ],xlim=limsX[p,], las=1,
       main=list(expression(paste(beta, " : large ThIC")),
                 expression(paste(pi, " : large ThIC")),
                 expression(paste(iota, " : large ThIC")),
                 expression(paste(.^0, theta^H, " : large ThIC")),
                 expression(paste(.^H, theta^{IC}," : large ThIC")) )[[p]])
  for (i in 1:5){
    lines(density(ThICdep[, p, ndx[i]], from=0), col=colpal[[i]], lwd=3, lty=1)
    lines(density(ThICind[, p, ndx[i]], from=0), col=colpal[[i]], lwd=3, lty=2)
  }
  legend(limsX[p, 1], lims[p, 2], c("joint dep", "miss ind") ,
         lty=c(1,2), lwd=3, bty = "n")
}

CrIRdep <- t(apply(ThICdep, MARGIN = c(2,3), quantile, prob=0.975))-
  t(apply(ThICdep, MARGIN = c(2,3), quantile, prob=0.025))
CrIRind <- t(apply(ThICind, MARGIN = c(2,3), quantile, prob=0.975))-
  t(apply(ThICind, MARGIN = c(2,3), quantile, prob=0.025))

VARdep <- t(apply(ThICdep, MARGIN = c(2,3), var))
VARind <- t(apply(ThICind, MARGIN = c(2,3), var))

par(mfrow=c(1, 2))
for (p in 1:5){
  par(mar=c(4, 4, 3,0.5))
  hist((VARdep[,p]- VARind[,p]), freq=F, breaks=10,  ylab="Density", xlab="Variance difference", las=1, adj=1,
       main=list(expression(paste(beta, " : large ThIC")),
                 expression(paste(pi, " : large ThIC")),
                 expression(paste(iota, " : large ThIC")),
                 expression(paste(.^0, theta^H, " : large ThIC")),
                 expression(paste(.^H, theta^{IC}," : large ThIC")) )[[p]])
  abline(v=0, lwd=3, lty=3)
  par(mar=c(4, 0.5, 3,4))
  hist((CrIRdep[,p]- CrIRind[,p]), freq=F, breaks=10, ylab=" ", 
       xlab="CrI Range difference",yaxt="n", 
       main=expression(paste("pairwise variability difference")), adj=0)
  abline(v=0, lwd=3, lty=3)
}


VARdiff <- VARdep-VARind
apply(VARdiff<0, MARGIN = 2, sum)/500
# dev.off()



# sub SECTION 4 -----------------------------------------------------------
rm(list=setdiff(ls(), c("Cexp","CexpT","Cgols","Cinf",
                        "CinfT","colpal","Crec","CrecT",
                        "Csus","CsusT", "loadwd")))
load(paste(loadwd, "OUTforTHESIS/DTCT.Robj", sep=""))

# DtHO

DtHOdep <- comp4DTCT$DtHO$dep
DtHOind <- comp4DTCT$DtHO$ind

 # pdf("fig5_9.pdf", #family="cmr10", 
 #     height =4.5, width=7)
set.seed(291218)
ndx <- sample(1:500, 5)
lims <- rbind(c(0,6), c(0,6), c(0,15000), c(0,25), c(0,12) )
limsX <- rbind(c(0.4,1.3), c(0,.8), c(0,0.0009), c(0,0.3), c(0, 0.45))
par(mfrow=c(1, 1), mar=c(2,4,3,1))
for (p in 1:5){
  plot(density(DtHOdep[, p, ndx[1]], from=0), col=0, ylim=lims[p, ], las=1, xlim=limsX[p,],
       main=list(expression(paste(beta, " : large DtHO")),
                 expression(paste(pi, " : large DtHO")),
                 expression(paste(iota, " : large DtHO")),
                 expression(paste(.^0, theta^H, " : large DtHO")),
                 expression(paste(.^H, theta^{IC}," : large DtHO")) )[[p]])
  for (i in 1:5){
    lines(density(DtHOdep[, p, ndx[i]], from=0), col=colpal[[i]], lwd=3, lty=1)
    lines(density(DtHOind[, p, ndx[i]], from=0), col=colpal[[i]], lwd=3, lty=2)
  }
  legend(0.002, lims[p, 2], c("joint dep", "miss ind") ,
         lty=c(1,2), lwd=3, bty = "n")
}

CrIRdep <- t(apply(DtHOdep, MARGIN = c(2,3), quantile, prob=0.975))-
  t(apply(DtHOdep, MARGIN = c(2,3), quantile, prob=0.025))
CrIRind <- t(apply(DtHOind, MARGIN = c(2,3), quantile, prob=0.975))-
  t(apply(DtHOind, MARGIN = c(2,3), quantile, prob=0.025))

VARdep <- t(apply(DtHOdep, MARGIN = c(2,3), var))
VARind <- t(apply(DtHOind, MARGIN = c(2,3), var))

par(mfrow=c(1, 2))
for (p in 1:5){
  par(mar=c(4, 4, 3,0.5))
  hist((VARdep[,p]- VARind[,p]), freq=F, breaks=10,  ylab="Density", xlab="Variance difference", las=1, adj=1,
        main=list(expression(paste(beta, " : large DtHO")),
                 expression(paste(pi, " : large DtHO")),
                 expression(paste(iota, " : large DtHO")),
                 expression(paste(.^0, theta^H, " : large DtHO")),
                 expression(paste(.^H, theta^{IC}," : large DtHO")) )[[p]])
  abline(v=0, lwd=3, lty=3)
  par(mar=c(4, 0.5, 3,4))
  hist((CrIRdep[,p]- CrIRind[,p]), freq=F, breaks=10, ylab=" ", 
       xlab="CrI Range difference",yaxt="n", 
       main=expression(paste("pairwise variability difference")), adj=0)
  abline(v=0, lwd=3, lty=3)
}

VARdiff <- VARdep-VARind
apply(VARdiff<0, MARGIN = 2, sum)/500

# DtIC

DtICdep <- comp4DTCT$DtIC$dep
DtICind <- comp4DTCT$DtIC$ind

set.seed(291218)
ndx <- sample(1:500, 5)
lims <- rbind(c(0,6), c(0,6), c(0,23000), c(0,20), c(0,25) )
limsX <- rbind(c(0.3,1.3), c(0,.8), c(0,0.0006), c(0,0.3), c(0, 0.3))
# limsX <- rbind(c(0,.75), c(0.17,.4), c(0,0.0002), c(0,0.5), c(0.7, 1))
par(mfrow=c(1, 1), mar=c(2,4,3,1))
for (p in 1:5){
  plot(density(DtICdep[, p, ndx[1]]), col=0, ylim=lims[p, ],xlim=limsX[p,], las=1,
       main=list(expression(paste(beta, " : large DtIC")),
                 expression(paste(pi, " : large DtIC")),
                 expression(paste(iota, " : large DtIC")),
                 expression(paste(.^0, theta^H, " : large DtIC")),
                 expression(paste(.^H, theta^{IC}," : large DtIC")) )[[p]])
  for (i in 1:5){
    lines(density(DtICdep[, p, ndx[i]], from=0), col=colpal[[i]], lwd=3, lty=1)
    lines(density(DtICind[, p, ndx[i]], from=0), col=colpal[[i]], lwd=3, lty=2)
  }
  legend(limsX[p, 1], lims[p, 2], c("joint dep", "miss ind") ,
         lty=c(1,2), lwd=3, bty = "n")
}

CrIRdep <- t(apply(DtICdep, MARGIN = c(2,3), quantile, prob=0.975))-
  t(apply(DtICdep, MARGIN = c(2,3), quantile, prob=0.025))
CrIRind <- t(apply(DtICind, MARGIN = c(2,3), quantile, prob=0.975))-
  t(apply(DtICind, MARGIN = c(2,3), quantile, prob=0.025))

VARdep <- t(apply(DtICdep, MARGIN = c(2,3), var))
VARind <- t(apply(DtICind, MARGIN = c(2,3), var))

par(mfrow=c(1, 2))
for (p in 1:5){
  par(mar=c(4, 4, 3,0.5))
  hist((VARdep[,p]- VARind[,p]), freq=F, breaks=10,  ylab="Density", xlab="Variance difference", las=1, adj=1,
       main=list(expression(paste(beta, " : large DtIC")),
                 expression(paste(pi, " : large DtIC")),
                 expression(paste(iota, " : large DtIC")),
                 expression(paste(.^0, theta^H, " : large DtIC")),
                 expression(paste(.^H, theta^{IC}," : large DtIC")) )[[p]])
  abline(v=0, lwd=3, lty=3)
  par(mar=c(4, 0.5, 3,4))
  hist((CrIRdep[,p]- CrIRind[,p]), freq=F, breaks=10, ylab=" ", 
       xlab="CrI Range difference",yaxt="n", 
       main=expression(paste("pairwise variability difference")), adj=0)
  abline(v=0, lwd=3, lty=3)
}

VARdiff <- VARdep-VARind
apply(VARdiff<0, MARGIN = 2, sum)/500
dev.off()
