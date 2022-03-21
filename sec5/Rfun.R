# to add in the R function
nuTOg <- function(nu, tv=rep(c(40:52,1:20), each=7)){
  gw.MEAN <- exp(nu[1]+nu[2]*sin(2*pi*tv/52)+nu[3]*cos(2*pi*tv/52))/7
  return(gw.MEAN)
}

OMGfun <- function(omg, totlength=(thetaconst$TMnd-thetaconst$TMbg)/thetaconst$OUTdt){
  out <- rep(c(omg[1], omg[2], omg[3], 1/prod(omg), omg[4], omg[5], omg[6]), totlength)
}


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
  Th <- c(exp(x[1]),          # beta
          toTh.mulpi(x[2:3]), # pi, iota
          exp(x[4])+0.5,      # dL
          exp(x[5])-1,        # k
          expit(x[6]),        # 0thH
          expit(x[7]),        # HthIC
          exp(x[8]),          # epsH
          x[9],               # nu1
          x[10],              # nu2
          x[11],              # nu3
          expit(x[12]),       # iotaF
          exp(x[13]),         # epsilon F
          exp(x[14]),         # omg1
          exp(x[15]),         # omg2
          exp(x[16]),         # omg3
          exp(x[17]),         # omg5
          exp(x[18]),         # omg6
          exp(x[19]))         # omg7
  return(Th)
}

toNorm <- function(Th){
  Norm <- c(log(Th[1]),            # beta
            toNorm.mulpi(Th[2:3]), # pi, iota
            log(Th[4]-0.5),        # dL
            log(Th[5]+1),          # k
            logit(Th[6]),          # 0thH
            logit(Th[7]),          # HthIC
            log(Th[8]),            # eps
            Th[9],                 # nu1
            Th[10],                # nu2
            Th[11],                # nu3
            logit(Th[12]),         # iotaF
            log(Th[13]),           # epsF   
            log(Th[14]),           # omg1   
            log(Th[15]),           # omg2   
            log(Th[16]),           # omg3   
            log(Th[17]),           # omg5   
            log(Th[18]),           # omg6   
            log(Th[19]))           # omg7
  return(Norm)
}


# PRIORS ------------------------------------------------------------------

lP <- function(y, p){
  switch(p, 
         # beta is a uniform between [0,4]
         dunif(exp(y[1]), 0, 4, log=T)+(y[1]),
         #pi is a beta with mean 0.375 and iota is a uniform between [0,0.1]
         dbeta(toTh.mulpi(y[2:3])[1], shape1 =37.5 , shape2=62.5, log=T)+
           dunif(toTh.mulpi(y[2:3])[2], 0, 0.1, log=T)+
           ((y[2])-2*log((1+exp(y[2]))))+
           ((y[3])-(2*log(1+exp(y[3]))+log(1+exp(y[2])))) ,
         dunif(toTh.mulpi(y[2:3])[2], 0, 0.1, log=T)+
           ((y[3])-(2*log(1+exp(y[3]))+log(1+exp(y[2])))),
         # dL is a LogNormal Centered in 2
         dlnorm(exp(y[4])+0.5, meanlog = log(2), sd=0.5)+(y[4]), 
         # K is a log(-1) noramal centered in 0
         dnorm((y[5]), 0, sd=1, log=T),
         # 0thH is a uniform [0,1]
         dunif(exp(y[6])/(exp(y[6])+1), 0, 1, log=T)+((y[6])-(2*log(1+exp(y[6])))),
         # HthIC is a uniform [0,1]
         dunif(exp(y[7])/(exp(y[7])+1), 0, 1, log=T)+((y[7])-(2*log(1+exp(y[7])))),
         # epsH is a is a uniform [0,100]
         dunif(exp(y[8]), 0, 100, log=T)+(y[8]), 
         # nus are normal with historical priors
         dnorm((y[9]), mean = 6.6, sd=0.17, log=T),
         dnorm((y[10]), mean = -0.2, sd=0.11, log=T),
         dnorm((y[11]), mean = 0.99, sd=0.08, log=T),
         # iota is a uniform 0,1
         dunif(exp(y[12])/(exp(y[12])+1), 0,1, log=T)+((y[12])-(2*log(1+exp(y[12])))),
         # epsF is a uniform 0, 10
         dunif(exp(y[13]), 0, 5000, log=T)+(y[13]),
         (dmvnorm(y[14:19], mean=rep(0,6), sigma=0.16*matrix(
           c(1, -1/6, -1/6, -1/6, -1/6, -1/6,
             -1/6, 1, -1/6, -1/6, -1/6, -1/6,
             -1/6, -1/6, 1, -1/6, -1/6, -1/6,
             -1/6, -1/6, -1/6, 1, -1/6, -1/6,
             -1/6, -1/6, -1/6, -1/6, 1, -1/6,
             -1/6, -1/6, -1/6, -1/6, -1/6, 1), byrow = T, ncol=6), log=T))
         )
}

lPj <- function(y){
  lp <-  dunif(exp(y[1]), 0, 4, log=T)+(y[1])+
  #pi is a beta with mean 0.375 and iota is a uniform between [0,0.1]
  dbeta(toTh.mulpi(y[2:3])[1], shape1 =37.5 , shape2=62.5, log=T)+
    dunif(toTh.mulpi(y[2:3])[2], 0, 0.1, log=T)+
    ((y[2])-2*log((1+exp(y[2]))))+
    ((y[3])-(2*log(1+exp(y[3]))+log(1+exp(y[2])))) +
  dunif(toTh.mulpi(y[2:3])[2], 0, 0.1, log=T)+
    ((y[3])-(2*log(1+exp(y[3]))+log(1+exp(y[2]))))+
  # dL is a LogNormal Centered in 2
  dlnorm(exp(y[4])+0.5, meanlog = log(2), sd=0.5)+(y[4])+ 
  # K is a log(-1) noramal centered in 0
  dnorm((y[5]), 0, sd=1, log=T)+
  # 0thH is a uniform [0,1]
  dunif(exp(y[6])/(exp(y[6])+1), 0, 1, log=T)+((y[6])-(2*log(1+exp(y[6]))))+
  # HthIC is a uniform [0,1]
  dunif(exp(y[7])/(exp(y[7])+1), 0, 1, log=T)+((y[7])-(2*log(1+exp(y[7]))))+
  # epsH is a is a uniform [0,100]
  dunif(exp(y[8]), 0, 100, log=T)+(y[8])+ 
  # nus are normal with historical priors
  dnorm((y[9]), mean = 6.6, sd=0.17, log=T)+
  dnorm((y[10]), mean = -0.2, sd=0.11, log=T)+
  dnorm((y[11]), mean = 0.99, sd=0.08, log=T)+
  # iota is a uniform 0,1
  dunif(exp(y[12])/(exp(y[12])+1), 0,1, log=T)+((y[12])-(2*log(1+exp(y[12]))))+
  # epsF is a uniform 0, 10
  dunif(exp(y[13]), 0, 5000, log=T)+(y[13])+
  (dmvnorm(y[14:19], mean=rep(0,6), sigma=0.16*matrix(
    c(1, -1/6, -1/6, -1/6, -1/6, -1/6,
      -1/6, 1, -1/6, -1/6, -1/6, -1/6,
      -1/6, -1/6, 1, -1/6, -1/6, -1/6,
      -1/6, -1/6, -1/6, 1, -1/6, -1/6,
      -1/6, -1/6, -1/6, -1/6, 1, -1/6,
      -1/6, -1/6, -1/6, -1/6, -1/6, 1), byrow = T, ncol=6), log=T))
  return(lp)
}


# Likelihoood -------------------------------------------------------------

llyGyV <- function(xi0u    ,   # weekly number of new infections, 
                   del0toF ,   # vector of relevant delay between infecion and Hospitalization,  
                   iotaF   ,   # parametere of the Gamma prior of 0thF
                   epsF    ,   # parametere of the Gamma prior of 0thF
                   gu      ,   # background
                   zetaGu  ,   # detection of GP consultation counts
                   nV      ,   # tested samples
                   yG      ,   # data: vector of GP consultation counts
                   yV      ,   # data: vector of positive samples counts
                   omegaU ){   # vector of the day of the week effect,
  xiFu  <- projCdet(probs = del0toF, sizes = xi0u)
  rU    <- epsF*iotaF+gu*epsF/xiFu
  etaU  <- (epsF+zetaGu*omegaU*xiFu)/epsF  
  alBB  <- iotaF*epsF
  btBB  <- (gu*epsF)/(xiFu)
  llyG <- sum(dnbinom(x=yG, size=rU, prob = 1/etaU, log=T))  
  llyVi <- (dbb(x=yV, N=nV, u=alBB, v=btBB, log=T))
  llyVi[which((nV==0)&(yV==0))] <- 0
  llyV  <- sum(llyVi)
  return(llyG+llyV)
} 

# transformation dor lambdaT
tolamt <- function(lamu){
  length(lamu)/7
  lamt1 <- rep(NA, (231)/7 )
  for (i in 1:((231)/7)){
    lamt1[i] <- sum(lamu[(i*7-6):(i*7)])
  }
  return(lamt1)
}
# CHAINS ------------------------------------------------------------------


adptMHsingle <- function(start, sdini, nsmpl, thin, stp.dpt,np, THETACONST=thetaconst){
  IT         <- nsmpl*thin
  chain      <- matrix(NA, nsmpl, 19)
  old        <- start
  naccept    <- matrix(0, ncol = 14, nrow = IT)
  it         <- 1
  oldTh      <- toTh(old) 
  # incidence
  oldlEu     <- LamTseeiirCtv(tm_X0 = setNames(c(p_S=THETACONST$N*(1-oldTh[2]-oldTh[3]), 
                                                 p_E=THETACONST$N*(oldTh[3]*((oldTh[4])/(oldTh[4]+THETACONST$dI))), 
                                                 p_I=THETACONST$N*(oldTh[3]*((THETACONST$dI)/(oldTh[4]+THETACONST$dI))), 
                                                 p_R=THETACONST$N*(oldTh[2])),c("p_S","p_E","p_R","p_R")),
                              tm_Theta = c(bt=oldTh[1],  sgm=2/oldTh[4], gmm=2/THETACONST$dI),
                              tm_begin = THETACONST$TMbg, tm_end = THETACONST$TMnd, 
                              tm_spd = THETACONST$TMspd, outdt = THETACONST$OUTdu, 
                              tm_Tk = THETACONST$vecHtb, tm_k = oldTh[5])
  # incidence in weeks
  oldlEt     <- tolamt(oldlEu)
  # background
  oldGU      <- nuTOg(oldTh[9:11])
  # dayoftheweekeffect
  oldOMG     <- OMGfun(oldTh[14:19])   
  # likelihood HIC
  oldjllc <- jllC(xi0=oldlEt, del0toH=THETACONST$dEtoH, delHtoIC=THETACONST$dHtoIC, 
                  th0toH=oldTh[6], thHtoIC=oldTh[7], epsH=oldTh[8], 
                  dH=THETACONST$detH , dIC=THETACONST$detIC , 
                  yH=DATA$yH , yIC=DATA$yIC , NP=np)
  # total likelihood
  oldLl      <- llyGyV(xi0u=oldlEu, del0toF=THETACONST$dEtoF, iotaF=oldTh[12],
                       epsF=oldTh[13], gu=oldGU, zetaGu=THETACONST$detG,
                       nV=THETACONST$sizeV, yG=DATA$yG, yV=DATA$yV,
                       omegaU=oldOMG)+oldjllc
  
  for (it in 1:IT){
    for (p in 1:13){
      new        <- old
      new[p]     <- rnorm(1, mean=old[p], sdini[p])
      newTh      <- toTh(new) 
      if(p %in% c(1,2,3,4,5)){
        newlEu     <- LamTseeiirCtv(tm_X0 = setNames(c(p_S=THETACONST$N*(1-newTh[2]-newTh[3]), 
                                                       p_E=THETACONST$N*(newTh[3]*((newTh[4])/(newTh[4]+THETACONST$dI))), 
                                                       p_I=THETACONST$N*(newTh[3]*((THETACONST$dI)/(newTh[4]+THETACONST$dI))), 
                                                       p_R=THETACONST$N*(newTh[2])),c("p_S","p_E","p_R","p_R")),
                                    tm_Theta = c(bt=newTh[1],  sgm=2/newTh[4], gmm=2/THETACONST$dI),
                                    tm_begin = THETACONST$TMbg, tm_end = THETACONST$TMnd, 
                                    tm_spd = THETACONST$TMspd, outdt = THETACONST$OUTdu, 
                                    tm_Tk = THETACONST$vecHtb, tm_k = newTh[5])
        newlEt    <- tolamt(newlEu)
        newGU     <- oldGU
        newOMG    <- oldOMG
      }else if(p %in% c(9,10,11)) {
        newlEu    <- oldlEu
        newlEt    <- oldlEt
        newGU     <- nuTOg(newTh[9:11])
        newOMG    <- oldOMG
      }
      if (p %in% c(1,2,3,4,5,6,7,8)){
        newjllc <- jllC(xi0=newlEt, del0toH=THETACONST$dEtoH, delHtoIC=THETACONST$dHtoIC, 
                        th0toH=newTh[6], thHtoIC=newTh[7], epsH=newTh[8], 
                        dH=THETACONST$detH , dIC=THETACONST$detIC , 
                        yH=DATA$yH , yIC=DATA$yIC , NP=np)
      }else{
        newjllc <- oldjllc
      }
      newLl      <- llyGyV(xi0u=newlEu, del0toF=THETACONST$dEtoF, iotaF=newTh[12],
                           epsF=newTh[13], gu=newGU, zetaGu=THETACONST$detG,
                           nV=THETACONST$sizeV, yG=DATA$yG, yV=DATA$yV,
                           omegaU=newOMG)+newjllc
      
      oldlp      <- lP(y=old, p=p)
      newlp      <- lP(y=new, p=p)
      lpoldInew  <- dnorm(old[p], mean=new[p], sdini[p], log=T)
      lpnewIold  <- dnorm(new[p], mean=old[p], sdini[p], log=T)
      logrho     <- (newlp-oldlp)+(newLl-oldLl)+(lpoldInew-lpnewIold)
      if(log(runif(1))<logrho){
        old            <- new 
        oldTh          <- newTh
        oldlEu         <- newlEu 
        oldlEt         <- newlEt
        oldGU          <- newGU
        oldLl          <- newLl
        oldjllc        <- newjllc
        oldOMG         <- newOMG
        naccept[it, p] <- 1
      }
      if(it > stp.dpt){
        accp          <- sum(naccept[(it-stp.dpt):it, p])/(stp.dpt)
        if (accp<=0.2&(sdini[p]>0.000001)){
          sdini[p]  <- sdini[p]*0.999
        }
        if (accp>=0.3&(sdini[p]<2)){
          sdini[p]  <- sdini[p]*1.001
        }
      }
    }
    # sampling of the omg
    p <- 14
      new        <- old
      new[14:19]     <- rmvnorm(1, mean=old[14:19], sigma=sdini[p]*diag(1,6))
      newTh      <- toTh(new) 
      newlEu    <- oldlEu
      newlEt    <- oldlEt
      newGU     <- oldGU
      newOMG    <- OMGfun(newTh[14:19])   
      newjllc    <- oldjllc
      newLl      <- llyGyV(xi0u=newlEu, del0toF=THETACONST$dEtoF, iotaF=newTh[12],
                           epsF=newTh[13], gu=newGU, zetaGu=THETACONST$detG,
                           nV=THETACONST$sizeV, yG=DATA$yG, yV=DATA$yV,
                           omegaU=newOMG)+newjllc
      oldlp      <- lP(y=old, p=p)
      newlp      <- lP(y=new, p=p)
      lpoldInew  <- dmvnorm(old[14:19], mean=new[14:19], sigma=sdini[p]*diag(1,6), log=T)
      lpnewIold  <- dmvnorm(new[14:19], mean=old[14:19], sigma=sdini[p]*diag(1,6), log=T)
      logrho     <- (newlp-oldlp)+(newLl-oldLl)+(lpoldInew-lpnewIold)
      if(log(runif(1))<logrho){
        old            <- new 
        oldTh          <- newTh
        oldlEu         <- newlEu 
        oldlEt         <- newlEt
        oldGU          <- newGU
        oldLl          <- newLl
        oldOMG         <- newOMG
        naccept[it, p] <- 1
      }
      if(it > stp.dpt){
        accp          <- sum(naccept[(it-stp.dpt):it, p])/(stp.dpt)
        if (accp<=0.2&(sdini[p]>0.000001)){
          sdini[p]  <- sdini[p]*0.999
        }
        if (accp>=0.3&(sdini[p]<2)){
          sdini[p]  <- sdini[p]*1.001
        }
      }
      # end of the sampling of the omg
    if ((it%%thin)==0){
      chain[it/thin, ] <- old
      print(paste(it/IT*100, "%"))
      print(sdini)
    }
  }
  return(list(sdini, chain, naccept))
}


runMHsingle <- function(start, sd, nsmpl, thin, np, THETACONST=thetaconst){
  IT         <- nsmpl*thin
  chain      <- matrix(NA, nsmpl, 19)
  old        <- start
  naccept    <- matrix(0, ncol = 14, nrow = IT)
  it         <- 1
  oldTh      <- toTh(old) 
  # incidence
  oldlEu     <- LamTseeiirCtv(tm_X0 = setNames(c(p_S=THETACONST$N*(1-oldTh[2]-oldTh[3]), 
                                                 p_E=THETACONST$N*(oldTh[3]*((oldTh[4])/(oldTh[4]+THETACONST$dI))), 
                                                 p_I=THETACONST$N*(oldTh[3]*((THETACONST$dI)/(oldTh[4]+THETACONST$dI))), 
                                                 p_R=THETACONST$N*(oldTh[2])),c("p_S","p_E","p_R","p_R")),
                              tm_Theta = c(bt=oldTh[1],  sgm=2/oldTh[4], gmm=2/THETACONST$dI),
                              tm_begin = THETACONST$TMbg, tm_end = THETACONST$TMnd, 
                              tm_spd = THETACONST$TMspd, outdt = THETACONST$OUTdu, 
                              tm_Tk = THETACONST$vecHtb, tm_k = oldTh[5])
  # incidence in weeks
  oldlEt     <- tolamt(oldlEu)
  # background
  oldGU      <- nuTOg(oldTh[9:11])
  # dayoftheweekeffect
  oldOMG     <- OMGfun(oldTh[14:19])   
  # likelihood HIC
  oldjllc <- jllC(xi0=oldlEt, del0toH=THETACONST$dEtoH, delHtoIC=THETACONST$dHtoIC, 
                  th0toH=oldTh[6], thHtoIC=oldTh[7], epsH=oldTh[8], 
                  dH=THETACONST$detH , dIC=THETACONST$detIC , 
                  yH=DATA$yH , yIC=DATA$yIC , NP=np)
  # total likelihood
  oldLl      <- llyGyV(xi0u=oldlEu, del0toF=THETACONST$dEtoF, iotaF=oldTh[12],
                       epsF=oldTh[13], gu=oldGU, zetaGu=THETACONST$detG,
                       nV=THETACONST$sizeV, yG=DATA$yG, yV=DATA$yV,
                       omegaU=oldOMG)+oldjllc
    
  
  for (it in 1:IT){
    for (p in 1:13){
      new        <- old
      new[p]     <- rnorm(1, mean=old[p], sd[p])
      newTh      <- toTh(new) 
      if(p %in% c(1,2,3,4,5)){
        newlEu     <- LamTseeiirCtv(tm_X0 = setNames(c(p_S=THETACONST$N*(1-newTh[2]-newTh[3]), 
                                                       p_E=THETACONST$N*(newTh[3]*((newTh[4])/(newTh[4]+THETACONST$dI))), 
                                                       p_I=THETACONST$N*(newTh[3]*((THETACONST$dI)/(newTh[4]+THETACONST$dI))), 
                                                       p_R=THETACONST$N*(newTh[2])),c("p_S","p_E","p_R","p_R")),
                                    tm_Theta = c(bt=newTh[1],  sgm=2/newTh[4], gmm=2/THETACONST$dI),
                                    tm_begin = THETACONST$TMbg, tm_end = THETACONST$TMnd, 
                                    tm_spd = THETACONST$TMspd, outdt = THETACONST$OUTdu, 
                                    tm_Tk = THETACONST$vecHtb, tm_k = newTh[5])
        newlEt    <- tolamt(newlEu)
        newGU     <- oldGU
        newOMG    <- oldOMG
      }else if(p %in% c(9,10,11)) {
        newlEu    <- oldlEu
        newlEt    <- oldlEt
        newGU     <- nuTOg(newTh[9:11])
        newOMG    <- oldOMG
      }
      if (p %in% c(1,2,3,4,5,6,7,8)){
        newjllc <- jllC(xi0=newlEt, del0toH=THETACONST$dEtoH, delHtoIC=THETACONST$dHtoIC, 
                        th0toH=newTh[6], thHtoIC=newTh[7], epsH=newTh[8], 
                        dH=THETACONST$detH , dIC=THETACONST$detIC , 
                        yH=DATA$yH , yIC=DATA$yIC , NP=np)
      }else{
        newjllc <- oldjllc
      }
      newLl      <- llyGyV(xi0u=newlEu, del0toF=THETACONST$dEtoF, iotaF=newTh[12],
                           epsF=newTh[13], gu=newGU, zetaGu=THETACONST$detG,
                           nV=THETACONST$sizeV, yG=DATA$yG, yV=DATA$yV,
                           omegaU=newOMG)+newjllc
        
      oldlp      <- lP(y=old, p=p)
      newlp      <- lP(y=new, p=p)
      lpoldInew  <- dnorm(old[p], mean=new[p], sd[p], log=T)
      lpnewIold  <- dnorm(new[p], mean=old[p], sd[p], log=T)
      logrho     <- (newlp-oldlp)+(newLl-oldLl)+(lpoldInew-lpnewIold)
      if(log(runif(1))<logrho){
        old            <- new 
        oldTh          <- newTh
        oldlEu         <- newlEu 
        oldlEt         <- newlEt
        oldGU          <- newGU
        oldLl          <- newLl
        oldjllc        <- newjllc
        oldOMG         <- newOMG
        naccept[it, p] <- 1
      }
    }
    # sampling of the omg
    p          <- 14
    new        <- old
    new[14:19] <- rmvnorm(1, mean=old[14:19], sigma=sd[p]*diag(1,6))
    newTh      <- toTh(new) 
    newlEu     <- oldlEu
    newlEt     <- oldlEt
    newGU      <- oldGU
    newOMG     <- OMGfun(newTh[14:19])   
    newjllc    <- oldjllc
    newLl      <- llyGyV(xi0u=newlEu, del0toF=THETACONST$dEtoF, iotaF=newTh[12],
                         epsF=newTh[13], gu=newGU, zetaGu=THETACONST$detG,
                         nV=THETACONST$sizeV, yG=DATA$yG, yV=DATA$yV,
                         omegaU=newOMG)+newjllc
    oldlp      <- lP(y=old, p=p)
    newlp      <- lP(y=new, p=p)
    lpoldInew  <- dmvnorm(old[14:19], mean=new[14:19], sigma=sd[p]*diag(1,6), log=T)
    lpnewIold  <- dmvnorm(new[14:19], mean=old[14:19], sigma=sd[p]*diag(1,6), log=T)
    logrho     <- (newlp-oldlp)+(newLl-oldLl)+(lpoldInew-lpnewIold)
    if(log(runif(1))<logrho){
      old            <- new 
      oldTh          <- newTh
      oldlEu         <- newlEu 
      oldlEt         <- newlEt
      oldGU          <- newGU
      oldLl          <- newLl
      oldOMG         <- newOMG
      naccept[it, p] <- 1
    }
    # end of the sampling of the omg
    if ((it%%thin)==0){
      chain[it/thin, ] <- old
      print(paste(it/IT*100, "%"))
    }
  }
  return(list(sd, chain, naccept))
}




adptMHblocked <- function(start, vcini, nsmpl, thin, stp.dpt, np, THETACONST=thetaconst){
  IT         <- nsmpl*thin
  chain      <- matrix(NA, nsmpl, 19)
  old        <- start
  naccept    <- rep(0, IT)
  it         <- 1
  scale      <- 1
  oldTh      <- toTh(old) 
  # incidence
  oldlEu     <- LamTseeiirCtv(tm_X0 = setNames(c(p_S=THETACONST$N*(1-oldTh[2]-oldTh[3]), 
                                                 p_E=THETACONST$N*(oldTh[3]*((oldTh[4])/(oldTh[4]+THETACONST$dI))), 
                                                 p_I=THETACONST$N*(oldTh[3]*((THETACONST$dI)/(oldTh[4]+THETACONST$dI))), 
                                                 p_R=THETACONST$N*(oldTh[2])),c("p_S","p_E","p_R","p_R")),
                              tm_Theta = c(bt=oldTh[1],  sgm=2/oldTh[4], gmm=2/THETACONST$dI),
                              tm_begin = THETACONST$TMbg, tm_end = THETACONST$TMnd, 
                              tm_spd = THETACONST$TMspd, outdt = THETACONST$OUTdu, 
                              tm_Tk = THETACONST$vecHtb, tm_k = oldTh[5])
  # incidence in weeks
  oldlEt     <- tolamt(oldlEu)
  # background
  oldGU      <- nuTOg(oldTh[9:11])
  # dayoftheweekeffect
  oldOMG     <- OMGfun(oldTh[14:19])   
  # likelihood
  oldLl      <- llyGyV(xi0u=oldlEu, del0toF=THETACONST$dEtoF, iotaF=oldTh[12],
                       epsF=oldTh[13], gu=oldGU, zetaGu=THETACONST$detG,
                       nV=THETACONST$sizeV, yG=DATA$yG, yV=DATA$yV,
                       omegaU=oldOMG)+
    jllC(xi0=oldlEt, del0toH=THETACONST$dEtoH, delHtoIC=THETACONST$dHtoIC, 
         th0toH=oldTh[6], thHtoIC=oldTh[7], epsH=oldTh[8], 
         dH=THETACONST$detH , dIC=THETACONST$detIC , 
         yH=DATA$yH , yIC=DATA$yIC , NP=np)
  
  for (it in 1:IT){
    new        <- rmvnorm(1, mean=old, sigma=scale*vcini)
    newTh      <- toTh(new) 
      newlEu     <- LamTseeiirCtv(tm_X0 = setNames(c(p_S=THETACONST$N*(1-newTh[2]-newTh[3]), 
                                                     p_E=THETACONST$N*(newTh[3]*((newTh[4])/(newTh[4]+THETACONST$dI))), 
                                                     p_I=THETACONST$N*(newTh[3]*((THETACONST$dI)/(newTh[4]+THETACONST$dI))), 
                                                     p_R=THETACONST$N*(newTh[2])),c("p_S","p_E","p_R","p_R")),
                                  tm_Theta = c(bt=newTh[1],  sgm=2/newTh[4], gmm=2/THETACONST$dI),
                                  tm_begin = THETACONST$TMbg, tm_end = THETACONST$TMnd, 
                                  tm_spd = THETACONST$TMspd, outdt = THETACONST$OUTdu, 
                                  tm_Tk = THETACONST$vecHtb, tm_k = newTh[5])
      newlEt    <- tolamt(newlEu)
      newGU     <- nuTOg(newTh[9:11])
      newOMG    <- OMGfun(newTh[14:19])  

    newLl      <- llyGyV(xi0u=newlEu, del0toF=THETACONST$dEtoF, iotaF=newTh[12],
                         epsF=newTh[13], gu=newGU, zetaGu=THETACONST$detG,
                         nV=THETACONST$sizeV, yG=DATA$yG, yV=DATA$yV,
                         omegaU=newOMG)+
      jllC(xi0=newlEt, del0toH=THETACONST$dEtoH, delHtoIC=THETACONST$dHtoIC, 
           th0toH=newTh[6], thHtoIC=newTh[7], epsH=newTh[8], 
           dH=THETACONST$detH , dIC=THETACONST$detIC , 
           yH=DATA$yH , yIC=DATA$yIC , NP=np)
    oldlp      <- lPj(y=old)
    newlp      <- lPj(y=new)
    lpoldInew  <- dmvnorm(old, mean=new, sigma=scale*vcini)
    lpnewIold  <- dmvnorm(new, mean=old, sigma=scale*vcini)
    logrho     <- (newlp-oldlp)+(newLl-oldLl)+(lpoldInew-lpnewIold)
    if(log(runif(1))<logrho){
      old            <- new 
      oldTh          <- newTh
      oldlEu         <- newlEu 
      oldlEt         <- newlEt
      oldGU          <- newGU
      oldLl          <- newLl
      oldOMG         <- newOMG
      naccept[it]    <- 1
    }
    if(it > stp.dpt){
      accp          <- sum(naccept[(it-stp.dpt):it])/(stp.dpt)
      if (accp<=0.2){
        scale  <- scale*0.999
      }
      if (accp>=0.3){
        scale  <- scale*1.001
      }
    }
    if ((it%%thin)==0){
      chain[it/thin, ] <- old
      print(paste(it/IT*100, "%"))
      print(scale)
    }
  }
  return(list(scale, chain, naccept))
}



runMHblocked <- function(start, vc, nsmpl, thin, np, THETACONST=thetaconst){
  IT         <- nsmpl*thin
  chain      <- matrix(NA, nsmpl, 19)
  old        <- start
  naccept    <- rep(0, IT)
  it         <- 1
  oldTh      <- toTh(old) 
  # incidence
  oldlEu     <- LamTseeiirCtv(tm_X0 = setNames(c(p_S=THETACONST$N*(1-oldTh[2]-oldTh[3]), 
                                                 p_E=THETACONST$N*(oldTh[3]*((oldTh[4])/(oldTh[4]+THETACONST$dI))), 
                                                 p_I=THETACONST$N*(oldTh[3]*((THETACONST$dI)/(oldTh[4]+THETACONST$dI))), 
                                                 p_R=THETACONST$N*(oldTh[2])),c("p_S","p_E","p_R","p_R")),
                              tm_Theta = c(bt=oldTh[1],  sgm=2/oldTh[4], gmm=2/THETACONST$dI),
                              tm_begin = THETACONST$TMbg, tm_end = THETACONST$TMnd, 
                              tm_spd = THETACONST$TMspd, outdt = THETACONST$OUTdu, 
                              tm_Tk = THETACONST$vecHtb, tm_k = oldTh[5])
  # incidence in weeks
  oldlEt     <- tolamt(oldlEu)
  # background
  oldGU      <- nuTOg(oldTh[9:11])
  # dayoftheweekeffect
  oldOMG     <- OMGfun(oldTh[14:19])   
  # likelihood
  oldLl      <- llyGyV(xi0u=oldlEu, del0toF=THETACONST$dEtoF, iotaF=oldTh[12],
                       epsF=oldTh[13], gu=oldGU, zetaGu=THETACONST$detG,
                       nV=THETACONST$sizeV, yG=DATA$yG, yV=DATA$yV,
                       omegaU=oldOMG)+
    jllC(xi0=oldlEt, del0toH=THETACONST$dEtoH, delHtoIC=THETACONST$dHtoIC, 
         th0toH=oldTh[6], thHtoIC=oldTh[7], epsH=oldTh[8], 
         dH=THETACONST$detH , dIC=THETACONST$detIC , 
         yH=DATA$yH , yIC=DATA$yIC , NP=np)
  
  for (it in 1:IT){
    new        <- rmvnorm(1, mean=old, sigma=vc)
    newTh      <- toTh(new) 
    newlEu     <- LamTseeiirCtv(tm_X0 = setNames(c(p_S=THETACONST$N*(1-newTh[2]-newTh[3]), 
                                                   p_E=THETACONST$N*(newTh[3]*((newTh[4])/(newTh[4]+THETACONST$dI))), 
                                                   p_I=THETACONST$N*(newTh[3]*((THETACONST$dI)/(newTh[4]+THETACONST$dI))), 
                                                   p_R=THETACONST$N*(newTh[2])),c("p_S","p_E","p_R","p_R")),
                                tm_Theta = c(bt=newTh[1],  sgm=2/newTh[4], gmm=2/THETACONST$dI),
                                tm_begin = THETACONST$TMbg, tm_end = THETACONST$TMnd, 
                                tm_spd = THETACONST$TMspd, outdt = THETACONST$OUTdu, 
                                tm_Tk = THETACONST$vecHtb, tm_k = newTh[5])
    newlEt    <- tolamt(newlEu)
    newGU     <- nuTOg(newTh[9:11])
    newOMG    <- OMGfun(newTh[14:19])  
    
    newLl      <- llyGyV(xi0u=newlEu, del0toF=THETACONST$dEtoF, iotaF=newTh[12],
                         epsF=newTh[13], gu=newGU, zetaGu=THETACONST$detG,
                         nV=THETACONST$sizeV, yG=DATA$yG, yV=DATA$yV,
                         omegaU=newOMG)+
      jllC(xi0=newlEt, del0toH=THETACONST$dEtoH, delHtoIC=THETACONST$dHtoIC, 
           th0toH=newTh[6], thHtoIC=newTh[7], epsH=newTh[8], 
           dH=THETACONST$detH , dIC=THETACONST$detIC , 
           yH=DATA$yH , yIC=DATA$yIC , NP=np)
    oldlp      <- lPj(y=old)
    newlp      <- lPj(y=new)
    lpoldInew  <- dmvnorm(old, mean=new, sigma=vc)
    lpnewIold  <- dmvnorm(new, mean=old, sigma=vc)
    logrho     <- (newlp-oldlp)+(newLl-oldLl)+(lpoldInew-lpnewIold)
    if(log(runif(1))<logrho){
      old            <- new 
      oldTh          <- newTh
      oldlEu         <- newlEu 
      oldlEt         <- newlEt
      oldGU          <- newGU
      oldLl          <- newLl
      oldOMG         <- newOMG
      naccept[it]    <- 1
    }
    if ((it%%thin)==0){
      chain[it/thin, ] <- old
      print(paste(it/IT*100, "%"))
    }
  }
  return(list(vc, chain, naccept))
}

