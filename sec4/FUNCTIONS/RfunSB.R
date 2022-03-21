
# PRIORS ------------------------------------------------------------------

lP <- function(y, p, pi=THETAvec[2]){
  switch(p, 
         dunif(exp(y[1]), 0, 4, log=T)+(y[1]), 
         dlnorm(toTh.mulpi(y[2:3])[1], meanlog = log(pi), sd=(.1/pi), log=T)+
           dunif(toTh.mulpi(y[2:3])[2], 0, 0.1, log=T)+
           ((y[2])-2*log((1+exp(y[2]))))+
           ((y[3])-(2*log(1+exp(y[3]))+log(1+exp(y[2])))) ,
         dunif(toTh.mulpi(y[2:3])[2], 0, 0.1, log=T)+
           ((y[3])-(2*log(1+exp(y[3]))+log(1+exp(y[2])))))
}


lPjoint <- function(y,pi=THETAvec[2]){
  lp <- (dunif(exp(y[1]), 0, 4, log=T)+
         dlnorm(toTh.mulpi(y[2:3])[1], meanlog = log(pi), sd=(.1/pi), log=T)+
         dunif(toTh.mulpi(y[2:3])[2], 0, 0.1, log=T))+
    # log Determinant of the jacobean
    ((y[1])+
       ((y[2])-2*log((1+exp(y[2]))))+
       ((y[3])-(2*log(1+exp(y[3]))+log(1+exp(y[2])))))
}


# MH functions for the DEPENDENT model ------------------------------------

adptMHsingle <- function(start, sdini, nsmpl, thin, np, stp.dpt, 
                         thetaconst=THETACONST,MCWM=T){
  IT         <- nsmpl*thin
  chain      <- matrix(NA, nsmpl, 3)
  old        <- start
  naccept    <- matrix(0, ncol = 3, nrow = IT)
  acclike    <- matrix(NA, ncol = 3, nrow = IT)
  it         <- 1
  oldTh      <- toTh(c(old, logit(thetaconst$thEtoH), logit(thetaconst$thHtoIC))) 
  oldlE      <- LamTseirC(tm_X0 = setNames(c(p_S=thetaconst$N*(1-oldTh[2]-oldTh[3]), 
                                             p_E=thetaconst$N*(oldTh[3]*((thetaconst$dL)/(thetaconst$dL+thetaconst$dI))), 
                                             p_I=thetaconst$N*(oldTh[3]*((thetaconst$dI)/(thetaconst$dL+thetaconst$dI))), 
                                             p_R=thetaconst$N*(oldTh[2])),c("p_S","p_E","p_R","p_R")),
                          tm_Theta = c(bt=oldTh[1],  sgm=1/thetaconst$dL, gmm=1/thetaconst$dI),
                          tm_begin = thetaconst$TMbg, tm_end = thetaconst$TMnd, 
                          tm_spd = thetaconst$TMspd, outdt = thetaconst$OUTdt)
  oldLl      <- loglikesevPF(lamt_tm = oldlE,
                             yH=DATA$yH,yIC=DATA$yIC,NP=np, 
                             thEtoH=oldTh[4],thHtoIC=oldTh[5], 
                             delEtoH=thetaconst$dEtoH, delHtoIC = thetaconst$dHtoIC, 
                             detH=thetaconst$detH, detIC=thetaconst$detIC)
  for (it in 1:IT){
    for (p in 1:3){
      new        <- old
      new[p]     <- rnorm(1, mean=old[p], sdini[p])
      newTh      <- toTh(c(new, logit(thetaconst$thEtoH), logit(thetaconst$thHtoIC))) 
      newlE    <- LamTseirC(tm_X0 = setNames(c(p_S=thetaconst$N*(1-newTh[2]-newTh[3]), 
                                                 p_E=thetaconst$N*(newTh[3]*((thetaconst$dL)/(thetaconst$dL+thetaconst$dI))), 
                                                 p_I=thetaconst$N*(newTh[3]*((thetaconst$dI)/(thetaconst$dL+thetaconst$dI))), 
                                                 p_R=thetaconst$N*(newTh[2])),c("p_S","p_E","p_R","p_R")),
                              tm_Theta = c(bt=newTh[1],  sgm=1/thetaconst$dL, gmm=1/thetaconst$dI),
                              tm_begin = thetaconst$TMbg, tm_end = thetaconst$TMnd, 
                              tm_spd = thetaconst$TMspd, outdt = thetaconst$OUTdt)
      newlL      <-  loglikesevPF(lamt_tm = newlE,
                                  yH=DATA$yH,yIC=DATA$yIC,NP=np, 
                                  thEtoH=newTh[4],thHtoIC=newTh[5], 
                                  delEtoH=thetaconst$dEtoH, delHtoIC=thetaconst$dHtoIC, 
                                  detH=thetaconst$detH, detIC=thetaconst$detIC)
      if(MCWM){
        oldLl <- loglikesevPF(lamt_tm = oldlE,
                              yH=DATA$yH,yIC=DATA$yIC,NP=np, 
                              thEtoH=oldTh[4],thHtoIC=oldTh[5], 
                              delEtoH=thetaconst$dEtoH, delHtoIC=thetaconst$dHtoIC, 
                              detH=thetaconst$detH, detIC=thetaconst$detIC)
        while(oldLl==-Inf){
          oldLl <- loglikesevPF(lamt_tm = oldlE,
                                yH=DATA$yH,yIC=DATA$yIC,NP=np, 
                                thEtoH=oldTh[4],thHtoIC=oldTh[5], 
                                delEtoH=thetaconst$dEtoH, delHtoIC=thetaconst$dHtoIC, 
                                detH=thetaconst$detH, detIC=thetaconst$detIC)
          print("alert: MCWM not adequate/insufficient particles")
        }
      }
      oldlp      <- lP(y=old, p=p)
      newlp      <- lP(y=new, p=p)
      lpoldInew  <- dnorm(old[p], mean=new[p], sdini[p], log=T)
      lpnewIold  <- dnorm(new[p], mean=old[p], sdini[p], log=T)
      logrho     <- (newlp-oldlp)+(newlL-oldLl)+(lpoldInew-lpnewIold)
      if(log(runif(1))<logrho){
        old            <- new 
        oldTh          <- newTh
        oldlE          <- newlE
        oldLl          <- newlL
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
      acclike[it, p]   <- oldLl 
    }
    if ((it%%thin)==0){
      chain[it/thin, ] <- old
      print(paste(it/IT*100, "%"))
      print(sdini)
    }
  }
  return(list(sdini, chain, naccept, acclike))
}

runMHsingle <- function(start, sd, nsmpl, thin, np,  
                         thetaconst=THETACONST,MCWM=T){
  IT         <- nsmpl*thin
  chain      <- matrix(NA, nsmpl, 3)
  old        <- start
  naccept    <- matrix(0, ncol = 3, nrow = IT)
  acclike    <- matrix(NA, ncol = 3, nrow = IT)
  it         <- 1
  oldTh      <- toTh(c(old, logit(thetaconst$thEtoH), logit(thetaconst$thHtoIC))) 
  oldlE      <- LamTseirC(tm_X0 = setNames(c(p_S=thetaconst$N*(1-oldTh[2]-oldTh[3]), 
                                             p_E=thetaconst$N*(oldTh[3]*((thetaconst$dL)/(thetaconst$dL+thetaconst$dI))), 
                                             p_I=thetaconst$N*(oldTh[3]*((thetaconst$dI)/(thetaconst$dL+thetaconst$dI))), 
                                             p_R=thetaconst$N*(oldTh[2])),c("p_S","p_E","p_R","p_R")),
                          tm_Theta = c(bt=oldTh[1],  sgm=1/thetaconst$dL, gmm=1/thetaconst$dI),
                          tm_begin = thetaconst$TMbg, tm_end = thetaconst$TMnd, 
                          tm_spd = thetaconst$TMspd, outdt = thetaconst$OUTdt)
  oldLl      <- loglikesevPF(lamt_tm = oldlE,
                             yH=DATA$yH,yIC=DATA$yIC,NP=np, 
                             thEtoH=oldTh[4],thHtoIC=oldTh[5], 
                             delEtoH=thetaconst$dEtoH, delHtoIC = thetaconst$dHtoIC, 
                             detH=thetaconst$detH, detIC=thetaconst$detIC)
  for (it in 1:IT){
    for (p in 1:3){
      new        <- old
      new[p]     <- rnorm(1, mean=old[p], sd[p])
      newTh      <- toTh(c(new, logit(thetaconst$thEtoH), logit(thetaconst$thHtoIC))) 
      newlE    <- LamTseirC(tm_X0 = setNames(c(p_S=thetaconst$N*(1-newTh[2]-newTh[3]), 
                                               p_E=thetaconst$N*(newTh[3]*((thetaconst$dL)/(thetaconst$dL+thetaconst$dI))), 
                                               p_I=thetaconst$N*(newTh[3]*((thetaconst$dI)/(thetaconst$dL+thetaconst$dI))), 
                                               p_R=thetaconst$N*(newTh[2])),c("p_S","p_E","p_R","p_R")),
                            tm_Theta = c(bt=newTh[1],  sgm=1/thetaconst$dL, gmm=1/thetaconst$dI),
                            tm_begin = thetaconst$TMbg, tm_end = thetaconst$TMnd, 
                            tm_spd = thetaconst$TMspd, outdt = thetaconst$OUTdt)
      newlL      <-  loglikesevPF(lamt_tm = newlE,
                                  yH=DATA$yH,yIC=DATA$yIC,NP=np, 
                                  thEtoH=newTh[4],thHtoIC=newTh[5], 
                                  delEtoH=thetaconst$dEtoH, delHtoIC=thetaconst$dHtoIC, 
                                  detH=thetaconst$detH, detIC=thetaconst$detIC)
      if(MCWM){
        oldLl <- loglikesevPF(lamt_tm = oldlE,
                              yH=DATA$yH,yIC=DATA$yIC,NP=np, 
                              thEtoH=oldTh[4],thHtoIC=oldTh[5], 
                              delEtoH=thetaconst$dEtoH, delHtoIC=thetaconst$dHtoIC, 
                              detH=thetaconst$detH, detIC=thetaconst$detIC)
        while(oldLl==-Inf){
          oldLl <- loglikesevPF(lamt_tm = oldlE,
                                yH=DATA$yH,yIC=DATA$yIC,NP=np, 
                                thEtoH=oldTh[4],thHtoIC=oldTh[5], 
                                delEtoH=thetaconst$dEtoH, delHtoIC=thetaconst$dHtoIC, 
                                detH=thetaconst$detH, detIC=thetaconst$detIC)
          print("alert: MCWM not adequate/insufficient particles")
        }
      }
      oldlp      <- lP(y=old, p=p)
      newlp      <- lP(y=new, p=p)
      lpoldInew  <- dnorm(old[p], mean=new[p], sd[p], log=T)
      lpnewIold  <- dnorm(new[p], mean=old[p], sd[p], log=T)
      logrho     <- (newlp-oldlp)+(newlL-oldLl)+(lpoldInew-lpnewIold)
      if(log(runif(1))<logrho){
        old            <- new 
        oldTh          <- newTh
        oldlE          <- newlE
        oldLl          <- newlL
        naccept[it, p] <- 1
      }
      acclike[it, p]   <- oldLl 
    }
    if ((it%%thin)==0){
      chain[it/thin, ] <- old
      print(paste(it/IT*100, "%"))
    }
  }
  return(list(sd, chain, naccept, acclike))
}

adptMHblocked <- function(start, vcini, nsmpl, thin, np, stp.dpt, 
                         thetaconst=THETACONST,MCWM=T){
  IT         <- nsmpl*thin
  chain      <- matrix(NA, nsmpl, 3)
  old        <- start
  naccept    <- rep(0, IT)
  acclike    <- rep(NA,IT)
  it         <- 1
  scale      <- 1
  oldTh      <- toTh(c(old, logit(thetaconst$thEtoH), logit(thetaconst$thHtoIC))) 
  oldlE      <- LamTseirC(tm_X0 = setNames(c(p_S=thetaconst$N*(1-oldTh[2]-oldTh[3]), 
                                             p_E=thetaconst$N*(oldTh[3]*((thetaconst$dL)/(thetaconst$dL+thetaconst$dI))), 
                                             p_I=thetaconst$N*(oldTh[3]*((thetaconst$dI)/(thetaconst$dL+thetaconst$dI))), 
                                             p_R=thetaconst$N*(oldTh[2])),c("p_S","p_E","p_R","p_R")),
                          tm_Theta = c(bt=oldTh[1],  sgm=1/thetaconst$dL, gmm=1/thetaconst$dI),
                          tm_begin = thetaconst$TMbg, tm_end = thetaconst$TMnd, 
                          tm_spd = thetaconst$TMspd, outdt = thetaconst$OUTdt)
  oldLl      <- loglikesevPF(lamt_tm = oldlE,
                             yH=DATA$yH,yIC=DATA$yIC,NP=np, 
                             thEtoH=oldTh[4],thHtoIC=oldTh[5], 
                             delEtoH=thetaconst$dEtoH, delHtoIC = thetaconst$dHtoIC, 
                             detH=thetaconst$detH, detIC=thetaconst$detIC)
  for (it in 1:IT){
      new       <- rmvnorm(1, mean=old, sigma=scale*vcini)
      newTh     <- toTh(c(new, logit(thetaconst$thEtoH), logit(thetaconst$thHtoIC))) 
      newlE     <- LamTseirC(tm_X0 = setNames(c(p_S=thetaconst$N*(1-newTh[2]-newTh[3]), 
                                               p_E=thetaconst$N*(newTh[3]*((thetaconst$dL)/(thetaconst$dL+thetaconst$dI))), 
                                               p_I=thetaconst$N*(newTh[3]*((thetaconst$dI)/(thetaconst$dL+thetaconst$dI))), 
                                               p_R=thetaconst$N*(newTh[2])),c("p_S","p_E","p_R","p_R")),
                            tm_Theta = c(bt=newTh[1],  sgm=1/thetaconst$dL, gmm=1/thetaconst$dI),
                            tm_begin = thetaconst$TMbg, tm_end = thetaconst$TMnd, 
                            tm_spd = thetaconst$TMspd, outdt = thetaconst$OUTdt)
      newlL      <-  loglikesevPF(lamt_tm = newlE,
                                  yH=DATA$yH,yIC=DATA$yIC,NP=np, 
                                  thEtoH=newTh[4],thHtoIC=newTh[5], 
                                  delEtoH=thetaconst$dEtoH, delHtoIC=thetaconst$dHtoIC, 
                                  detH=thetaconst$detH, detIC=thetaconst$detIC)
      if(MCWM){
        oldLl <- loglikesevPF(lamt_tm = oldlE,
                              yH=DATA$yH,yIC=DATA$yIC,NP=np, 
                              thEtoH=oldTh[4],thHtoIC=oldTh[5], 
                              delEtoH=thetaconst$dEtoH, delHtoIC=thetaconst$dHtoIC, 
                              detH=thetaconst$detH, detIC=thetaconst$detIC)
        while(oldLl==-Inf){
          oldLl <- loglikesevPF(lamt_tm = oldlE,
                                yH=DATA$yH,yIC=DATA$yIC,NP=np, 
                                thEtoH=oldTh[4],thHtoIC=oldTh[5], 
                                delEtoH=thetaconst$dEtoH, delHtoIC=thetaconst$dHtoIC, 
                                detH=thetaconst$detH, detIC=thetaconst$detIC)
          print("alert: MCWM not adequate/insufficient particles")
        }
      }
      oldlp      <- lPjoint(y=old)
      newlp      <- lPjoint(y=new)
      lpoldInew  <- dmvnorm(old, mean=new, sigma=scale*vcini, log=T)
      lpnewIold  <- dmvnorm(new, mean=old, sigma=scale*vcini, log=T)
      logrho     <- (newlp-oldlp)+(newlL-oldLl)+(lpoldInew-lpnewIold)
      if(log(runif(1))<logrho){
        old            <- new 
        oldTh          <- newTh
        oldlE          <- newlE
        naccept[it]    <- 1
        oldLl          <- newlL
      }
      if(it > stp.dpt){
        accp          <- sum(naccept[(it-stp.dpt):it])/(stp.dpt)
        if (accp<=0.2&(scale>0.000001)){
          scale  <- scale*0.999
        }
        if (accp>=0.3&(scale<2)){
          scale  <- scale*1.001
        }
      }
      acclike[it]   <- oldLl 
    if ((it%%thin)==0){
      chain[it/thin, ] <- old
      print(paste(it/IT*100, "%"))
      print(scale)
    }
  }
  return(list(scale, chain, naccept, acclike))
}


runMHblocked <- function(start, vcscaled, nsmpl, thin, np, 
                          thetaconst=THETACONST,MCWM=T){
  IT         <- nsmpl*thin
  chain      <- matrix(NA, nsmpl, 3)
  old        <- start
  naccept    <- rep(0, IT)
  acclike    <- rep(NA,IT)
  it         <- 1
  oldTh      <- toTh(c(old, logit(thetaconst$thEtoH), logit(thetaconst$thHtoIC))) 
  oldlE      <- LamTseirC(tm_X0 = setNames(c(p_S=thetaconst$N*(1-oldTh[2]-oldTh[3]), 
                                             p_E=thetaconst$N*(oldTh[3]*((thetaconst$dL)/(thetaconst$dL+thetaconst$dI))), 
                                             p_I=thetaconst$N*(oldTh[3]*((thetaconst$dI)/(thetaconst$dL+thetaconst$dI))), 
                                             p_R=thetaconst$N*(oldTh[2])),c("p_S","p_E","p_R","p_R")),
                          tm_Theta = c(bt=oldTh[1],  sgm=1/thetaconst$dL, gmm=1/thetaconst$dI),
                          tm_begin = thetaconst$TMbg, tm_end = thetaconst$TMnd, 
                          tm_spd = thetaconst$TMspd, outdt = thetaconst$OUTdt)
  oldLl      <- loglikesevPF(lamt_tm = oldlE,
                             yH=DATA$yH,yIC=DATA$yIC,NP=np, 
                             thEtoH=oldTh[4],thHtoIC=oldTh[5], 
                             delEtoH=thetaconst$dEtoH, delHtoIC = thetaconst$dHtoIC, 
                             detH=thetaconst$detH, detIC=thetaconst$detIC)
  for (it in 1:IT){
    new       <- rmvnorm(1, mean=old, sigma=vcscaled)
    newTh     <- toTh(c(new, logit(thetaconst$thEtoH), logit(thetaconst$thHtoIC))) 
    newlE     <- LamTseirC(tm_X0 = setNames(c(p_S=thetaconst$N*(1-newTh[2]-newTh[3]), 
                                              p_E=thetaconst$N*(newTh[3]*((thetaconst$dL)/(thetaconst$dL+thetaconst$dI))), 
                                              p_I=thetaconst$N*(newTh[3]*((thetaconst$dI)/(thetaconst$dL+thetaconst$dI))), 
                                              p_R=thetaconst$N*(newTh[2])),c("p_S","p_E","p_R","p_R")),
                           tm_Theta = c(bt=newTh[1],  sgm=1/thetaconst$dL, gmm=1/thetaconst$dI),
                           tm_begin = thetaconst$TMbg, tm_end = thetaconst$TMnd, 
                           tm_spd = thetaconst$TMspd, outdt = thetaconst$OUTdt)
    newlL      <-  loglikesevPF(lamt_tm = newlE,
                                yH=DATA$yH,yIC=DATA$yIC,NP=np, 
                                thEtoH=newTh[4],thHtoIC=newTh[5], 
                                delEtoH=thetaconst$dEtoH, delHtoIC=thetaconst$dHtoIC, 
                                detH=thetaconst$detH, detIC=thetaconst$detIC)
    if(MCWM){
      oldLl <- loglikesevPF(lamt_tm = oldlE,
                            yH=DATA$yH,yIC=DATA$yIC,NP=np, 
                            thEtoH=oldTh[4],thHtoIC=oldTh[5], 
                            delEtoH=thetaconst$dEtoH, delHtoIC=thetaconst$dHtoIC, 
                            detH=thetaconst$detH, detIC=thetaconst$detIC)
      while(oldLl==-Inf){
        oldLl <- loglikesevPF(lamt_tm = oldlE,
                              yH=DATA$yH,yIC=DATA$yIC,NP=np, 
                              thEtoH=oldTh[4],thHtoIC=oldTh[5], 
                              delEtoH=thetaconst$dEtoH, delHtoIC=thetaconst$dHtoIC, 
                              detH=thetaconst$detH, detIC=thetaconst$detIC)
        print("alert: MCWM not adequate/insufficient particles")
      }
    }
    oldlp      <- lPjoint(y=old)
    newlp      <- lPjoint(y=new)
    lpoldInew  <- dmvnorm(old, mean=new, sigma=vcscaled, log=T)
    lpnewIold  <- dmvnorm(new, mean=old, sigma=vcscaled, log=T)
    logrho     <- (newlp-oldlp)+(newlL-oldLl)+(lpoldInew-lpnewIold)
    if(log(runif(1))<logrho){
      old            <- new 
      oldTh          <- newTh
      oldlE          <- newlE
      naccept[it]    <- 1
      oldLl          <- newlL
    }
    acclike[it]   <- oldLl 
    if ((it%%thin)==0){
      chain[it/thin, ] <- old
      print(paste(it/IT*100, "%"))
    }
  }
  return(list(vcscaled, chain, naccept, acclike))
}




# MH functions for the INDEPENDENT model ------------------------------------

INDadptMHsingle <- function(start, sdini, nsmpl, thin, stp.dpt, 
                         thetaconst=THETACONST,MCWM=T){
  IT         <- nsmpl*thin
  chain      <- matrix(NA, nsmpl, 3)
  old        <- start
  naccept    <- matrix(0, ncol = 3, nrow = IT)
  acclike    <- matrix(NA, ncol = 3, nrow = IT)
  it         <- 1
  oldTh      <- toTh(c(old, logit(thetaconst$thEtoH), logit(thetaconst$thHtoIC))) 
  oldlE      <- LamTseirC(tm_X0 = setNames(c(p_S=thetaconst$N*(1-oldTh[2]-oldTh[3]), 
                                             p_E=thetaconst$N*(oldTh[3]*((thetaconst$dL)/(thetaconst$dL+thetaconst$dI))), 
                                             p_I=thetaconst$N*(oldTh[3]*((thetaconst$dI)/(thetaconst$dL+thetaconst$dI))), 
                                             p_R=thetaconst$N*(oldTh[2])),c("p_S","p_E","p_R","p_R")),
                          tm_Theta = c(bt=oldTh[1],  sgm=1/thetaconst$dL, gmm=1/thetaconst$dI),
                          tm_begin = thetaconst$TMbg, tm_end = thetaconst$TMnd, 
                          tm_spd = thetaconst$TMspd, outdt = thetaconst$OUTdt)
  oldLl      <- loglikesevIND(lamt_tm = oldlE,
                             yH=DATA$yH,yIC=DATA$yIC,
                             thEtoH=oldTh[4],thHtoIC=oldTh[5], 
                             delEtoH=thetaconst$dEtoH, delHtoIC = thetaconst$dHtoIC, 
                             detH=thetaconst$detH, detIC=thetaconst$detIC)
  for (it in 1:IT){
    for (p in 1:3){
      new        <- old
      new[p]     <- rnorm(1, mean=old[p], sdini[p])
      newTh      <- toTh(c(new, logit(thetaconst$thEtoH), logit(thetaconst$thHtoIC))) 
      newlE    <- LamTseirC(tm_X0 = setNames(c(p_S=thetaconst$N*(1-newTh[2]-newTh[3]), 
                                               p_E=thetaconst$N*(newTh[3]*((thetaconst$dL)/(thetaconst$dL+thetaconst$dI))), 
                                               p_I=thetaconst$N*(newTh[3]*((thetaconst$dI)/(thetaconst$dL+thetaconst$dI))), 
                                               p_R=thetaconst$N*(newTh[2])),c("p_S","p_E","p_R","p_R")),
                            tm_Theta = c(bt=newTh[1],  sgm=1/thetaconst$dL, gmm=1/thetaconst$dI),
                            tm_begin = thetaconst$TMbg, tm_end = thetaconst$TMnd, 
                            tm_spd = thetaconst$TMspd, outdt = thetaconst$OUTdt)
      newlL      <-  loglikesevIND(lamt_tm = newlE,
                                  yH=DATA$yH,yIC=DATA$yIC,
                                  thEtoH=newTh[4],thHtoIC=newTh[5], 
                                  delEtoH=thetaconst$dEtoH, delHtoIC=thetaconst$dHtoIC, 
                                  detH=thetaconst$detH, detIC=thetaconst$detIC)
      if(MCWM){
        oldLl <- loglikesevIND(lamt_tm = oldlE,
                              yH=DATA$yH,yIC=DATA$yIC,
                              thEtoH=oldTh[4],thHtoIC=oldTh[5], 
                              delEtoH=thetaconst$dEtoH, delHtoIC=thetaconst$dHtoIC, 
                              detH=thetaconst$detH, detIC=thetaconst$detIC)
        while(oldLl==-Inf){
          oldLl <- loglikesevIND(lamt_tm = oldlE,
                                yH=DATA$yH,yIC=DATA$yIC,
                                thEtoH=oldTh[4],thHtoIC=oldTh[5], 
                                delEtoH=thetaconst$dEtoH, delHtoIC=thetaconst$dHtoIC, 
                                detH=thetaconst$detH, detIC=thetaconst$detIC)
          print("alert: MCWM not adequate/insufficient particles")
        }
      }
      oldlp      <- lP(y=old, p=p)
      newlp      <- lP(y=new, p=p)
      lpoldInew  <- dnorm(old[p], mean=new[p], sdini[p], log=T)
      lpnewIold  <- dnorm(new[p], mean=old[p], sdini[p], log=T)
      logrho     <- (newlp-oldlp)+(newlL-oldLl)+(lpoldInew-lpnewIold)
      if(log(runif(1))<logrho){
        old            <- new 
        oldTh          <- newTh
        oldlE          <- newlE
        oldLl          <- newlL
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
      acclike[it, p]   <- oldLl 
    }
    if ((it%%thin)==0){
      chain[it/thin, ] <- old
      # print(paste(it/IT*100, "%"))
      # print(sdini)
    }
  }
  return(list(sdini, chain, naccept, acclike))
}

INDrunMHsingle <- function(start, sd, nsmpl, thin,  
                        thetaconst=THETACONST,MCWM=T){
  IT         <- nsmpl*thin
  chain      <- matrix(NA, nsmpl, 3)
  old        <- start
  naccept    <- matrix(0, ncol = 3, nrow = IT)
  acclike    <- matrix(NA, ncol = 3, nrow = IT)
  it         <- 1
  oldTh      <- toTh(c(old, logit(thetaconst$thEtoH), logit(thetaconst$thHtoIC))) 
  oldlE      <- LamTseirC(tm_X0 = setNames(c(p_S=thetaconst$N*(1-oldTh[2]-oldTh[3]), 
                                             p_E=thetaconst$N*(oldTh[3]*((thetaconst$dL)/(thetaconst$dL+thetaconst$dI))), 
                                             p_I=thetaconst$N*(oldTh[3]*((thetaconst$dI)/(thetaconst$dL+thetaconst$dI))), 
                                             p_R=thetaconst$N*(oldTh[2])),c("p_S","p_E","p_R","p_R")),
                          tm_Theta = c(bt=oldTh[1],  sgm=1/thetaconst$dL, gmm=1/thetaconst$dI),
                          tm_begin = thetaconst$TMbg, tm_end = thetaconst$TMnd, 
                          tm_spd = thetaconst$TMspd, outdt = thetaconst$OUTdt)
  oldLl      <- loglikesevIND(lamt_tm = oldlE,
                             yH=DATA$yH,yIC=DATA$yIC,
                             thEtoH=oldTh[4],thHtoIC=oldTh[5], 
                             delEtoH=thetaconst$dEtoH, delHtoIC = thetaconst$dHtoIC, 
                             detH=thetaconst$detH, detIC=thetaconst$detIC)
  for (it in 1:IT){
    for (p in 1:3){
      new        <- old
      new[p]     <- rnorm(1, mean=old[p], sd[p])
      newTh      <- toTh(c(new, logit(thetaconst$thEtoH), logit(thetaconst$thHtoIC))) 
      newlE    <- LamTseirC(tm_X0 = setNames(c(p_S=thetaconst$N*(1-newTh[2]-newTh[3]), 
                                               p_E=thetaconst$N*(newTh[3]*((thetaconst$dL)/(thetaconst$dL+thetaconst$dI))), 
                                               p_I=thetaconst$N*(newTh[3]*((thetaconst$dI)/(thetaconst$dL+thetaconst$dI))), 
                                               p_R=thetaconst$N*(newTh[2])),c("p_S","p_E","p_R","p_R")),
                            tm_Theta = c(bt=newTh[1],  sgm=1/thetaconst$dL, gmm=1/thetaconst$dI),
                            tm_begin = thetaconst$TMbg, tm_end = thetaconst$TMnd, 
                            tm_spd = thetaconst$TMspd, outdt = thetaconst$OUTdt)
      newlL      <-  loglikesevIND(lamt_tm = newlE,
                                  yH=DATA$yH,yIC=DATA$yIC,
                                  thEtoH=newTh[4],thHtoIC=newTh[5], 
                                  delEtoH=thetaconst$dEtoH, delHtoIC=thetaconst$dHtoIC, 
                                  detH=thetaconst$detH, detIC=thetaconst$detIC)
      if(MCWM){
        oldLl <- loglikesevIND(lamt_tm = oldlE,
                              yH=DATA$yH,yIC=DATA$yIC,
                              thEtoH=oldTh[4],thHtoIC=oldTh[5], 
                              delEtoH=thetaconst$dEtoH, delHtoIC=thetaconst$dHtoIC, 
                              detH=thetaconst$detH, detIC=thetaconst$detIC)
        while(oldLl==-Inf){
          oldLl <- loglikesevIND(lamt_tm = oldlE,
                                yH=DATA$yH,yIC=DATA$yIC,
                                thEtoH=oldTh[4],thHtoIC=oldTh[5], 
                                delEtoH=thetaconst$dEtoH, delHtoIC=thetaconst$dHtoIC, 
                                detH=thetaconst$detH, detIC=thetaconst$detIC)
          print("alert: MCWM not adequate/insufficient particles")
        }
      }
      oldlp      <- lP(y=old, p=p)
      newlp      <- lP(y=new, p=p)
      lpoldInew  <- dnorm(old[p], mean=new[p], sd[p], log=T)
      lpnewIold  <- dnorm(new[p], mean=old[p], sd[p], log=T)
      logrho     <- (newlp-oldlp)+(newlL-oldLl)+(lpoldInew-lpnewIold)
      if(log(runif(1))<logrho){
        old            <- new 
        oldTh          <- newTh
        oldlE          <- newlE
        oldLl          <- newlL
        naccept[it, p] <- 1
      }
      acclike[it, p]   <- oldLl 
    }
    if ((it%%thin)==0){
      chain[it/thin, ] <- old
      # print(paste(it/IT*100, "%"))
    }
  }
  return(list(sd, chain, naccept, acclike))
}

INDadptMHblocked <- function(start, vcini, nsmpl, thin, stp.dpt, 
                          thetaconst=THETACONST,MCWM=T){
  IT         <- nsmpl*thin
  chain      <- matrix(NA, nsmpl, 3)
  old        <- start
  naccept    <- rep(0, IT)
  acclike    <- rep(NA,IT)
  it         <- 1
  scale      <- 1
  oldTh      <- toTh(c(old, logit(thetaconst$thEtoH), logit(thetaconst$thHtoIC))) 
  oldlE      <- LamTseirC(tm_X0 = setNames(c(p_S=thetaconst$N*(1-oldTh[2]-oldTh[3]), 
                                             p_E=thetaconst$N*(oldTh[3]*((thetaconst$dL)/(thetaconst$dL+thetaconst$dI))), 
                                             p_I=thetaconst$N*(oldTh[3]*((thetaconst$dI)/(thetaconst$dL+thetaconst$dI))), 
                                             p_R=thetaconst$N*(oldTh[2])),c("p_S","p_E","p_R","p_R")),
                          tm_Theta = c(bt=oldTh[1],  sgm=1/thetaconst$dL, gmm=1/thetaconst$dI),
                          tm_begin = thetaconst$TMbg, tm_end = thetaconst$TMnd, 
                          tm_spd = thetaconst$TMspd, outdt = thetaconst$OUTdt)
  oldLl      <- loglikesevIND(lamt_tm = oldlE,
                             yH=DATA$yH,yIC=DATA$yIC, 
                             thEtoH=oldTh[4],thHtoIC=oldTh[5], 
                             delEtoH=thetaconst$dEtoH, delHtoIC = thetaconst$dHtoIC, 
                             detH=thetaconst$detH, detIC=thetaconst$detIC)
  for (it in 1:IT){
    new       <- rmvnorm(1, mean=old, sigma=scale*vcini)
    newTh     <- toTh(c(new, logit(thetaconst$thEtoH), logit(thetaconst$thHtoIC))) 
    newlE     <- LamTseirC(tm_X0 = setNames(c(p_S=thetaconst$N*(1-newTh[2]-newTh[3]), 
                                              p_E=thetaconst$N*(newTh[3]*((thetaconst$dL)/(thetaconst$dL+thetaconst$dI))), 
                                              p_I=thetaconst$N*(newTh[3]*((thetaconst$dI)/(thetaconst$dL+thetaconst$dI))), 
                                              p_R=thetaconst$N*(newTh[2])),c("p_S","p_E","p_R","p_R")),
                           tm_Theta = c(bt=newTh[1],  sgm=1/thetaconst$dL, gmm=1/thetaconst$dI),
                           tm_begin = thetaconst$TMbg, tm_end = thetaconst$TMnd, 
                           tm_spd = thetaconst$TMspd, outdt = thetaconst$OUTdt)
    newlL      <-  loglikesevIND(lamt_tm = newlE,
                                yH=DATA$yH,yIC=DATA$yIC,
                                thEtoH=newTh[4],thHtoIC=newTh[5], 
                                delEtoH=thetaconst$dEtoH, delHtoIC=thetaconst$dHtoIC, 
                                detH=thetaconst$detH, detIC=thetaconst$detIC)
    if(MCWM){
      oldLl <- loglikesevIND(lamt_tm = oldlE,
                            yH=DATA$yH,yIC=DATA$yIC,
                            thEtoH=oldTh[4],thHtoIC=oldTh[5], 
                            delEtoH=thetaconst$dEtoH, delHtoIC=thetaconst$dHtoIC, 
                            detH=thetaconst$detH, detIC=thetaconst$detIC)
      while(oldLl==-Inf){
        oldLl <- loglikesevIND(lamt_tm = oldlE,
                              yH=DATA$yH,yIC=DATA$yIC,
                              thEtoH=oldTh[4],thHtoIC=oldTh[5], 
                              delEtoH=thetaconst$dEtoH, delHtoIC=thetaconst$dHtoIC, 
                              detH=thetaconst$detH, detIC=thetaconst$detIC)
        print("alert: MCWM not adequate/insufficient particles")
      }
    }
    oldlp      <- lPjoint(y=old)
    newlp      <- lPjoint(y=new)
    lpoldInew  <- dmvnorm(old, mean=new, sigma=scale*vcini, log=T)
    lpnewIold  <- dmvnorm(new, mean=old, sigma=scale*vcini, log=T)
    logrho     <- (newlp-oldlp)+(newlL-oldLl)+(lpoldInew-lpnewIold)
    if(log(runif(1))<logrho){
      old            <- new 
      oldTh          <- newTh
      oldlE          <- newlE
      naccept[it]    <- 1
      oldLl          <- newlL
    }
    if(it > stp.dpt){
      accp          <- sum(naccept[(it-stp.dpt):it])/(stp.dpt)
      if (accp<=0.2&(scale>0.000001)){
        scale  <- scale*0.999
      }
      if (accp>=0.3&(scale<2)){
        scale  <- scale*1.001
      }
    }
    acclike[it]   <- oldLl 
    if ((it%%thin)==0){
      chain[it/thin, ] <- old
      # print(paste(it/IT*100, "%"))
      # print(scale)
    }
  }
  return(list(scale, chain, naccept, acclike))
}

INDrunMHblocked <- function(start, vcscaled, nsmpl, thin, 
                         thetaconst=THETACONST,MCWM=T){
  IT         <- nsmpl*thin
  chain      <- matrix(NA, nsmpl, 3)
  old        <- start
  naccept    <- rep(0, IT)
  acclike    <- rep(NA,IT)
  it         <- 1
  oldTh      <- toTh(c(old, logit(thetaconst$thEtoH), logit(thetaconst$thHtoIC))) 
  oldlE      <- LamTseirC(tm_X0 = setNames(c(p_S=thetaconst$N*(1-oldTh[2]-oldTh[3]), 
                                             p_E=thetaconst$N*(oldTh[3]*((thetaconst$dL)/(thetaconst$dL+thetaconst$dI))), 
                                             p_I=thetaconst$N*(oldTh[3]*((thetaconst$dI)/(thetaconst$dL+thetaconst$dI))), 
                                             p_R=thetaconst$N*(oldTh[2])),c("p_S","p_E","p_R","p_R")),
                          tm_Theta = c(bt=oldTh[1],  sgm=1/thetaconst$dL, gmm=1/thetaconst$dI),
                          tm_begin = thetaconst$TMbg, tm_end = thetaconst$TMnd, 
                          tm_spd = thetaconst$TMspd, outdt = thetaconst$OUTdt)
  oldLl      <- loglikesevIND(lamt_tm = oldlE,
                             yH=DATA$yH,yIC=DATA$yIC,
                             thEtoH=oldTh[4],thHtoIC=oldTh[5], 
                             delEtoH=thetaconst$dEtoH, delHtoIC = thetaconst$dHtoIC, 
                             detH=thetaconst$detH, detIC=thetaconst$detIC)
  for (it in 1:IT){
    new       <- rmvnorm(1, mean=old, sigma=vcscaled)
    newTh     <- toTh(c(new, logit(thetaconst$thEtoH), logit(thetaconst$thHtoIC))) 
    newlE     <- LamTseirC(tm_X0 = setNames(c(p_S=thetaconst$N*(1-newTh[2]-newTh[3]), 
                                              p_E=thetaconst$N*(newTh[3]*((thetaconst$dL)/(thetaconst$dL+thetaconst$dI))), 
                                              p_I=thetaconst$N*(newTh[3]*((thetaconst$dI)/(thetaconst$dL+thetaconst$dI))), 
                                              p_R=thetaconst$N*(newTh[2])),c("p_S","p_E","p_R","p_R")),
                           tm_Theta = c(bt=newTh[1],  sgm=1/thetaconst$dL, gmm=1/thetaconst$dI),
                           tm_begin = thetaconst$TMbg, tm_end = thetaconst$TMnd, 
                           tm_spd = thetaconst$TMspd, outdt = thetaconst$OUTdt)
    newlL      <-  loglikesevIND(lamt_tm = newlE,
                                yH=DATA$yH,yIC=DATA$yIC,
                                thEtoH=newTh[4],thHtoIC=newTh[5], 
                                delEtoH=thetaconst$dEtoH, delHtoIC=thetaconst$dHtoIC, 
                                detH=thetaconst$detH, detIC=thetaconst$detIC)
    if(MCWM){
      oldLl <- loglikesevIND(lamt_tm = oldlE,
                            yH=DATA$yH,yIC=DATA$yIC,
                            thEtoH=oldTh[4],thHtoIC=oldTh[5], 
                            delEtoH=thetaconst$dEtoH, delHtoIC=thetaconst$dHtoIC, 
                            detH=thetaconst$detH, detIC=thetaconst$detIC)
      while(oldLl==-Inf){
        oldLl <- loglikesevIND(lamt_tm = oldlE,
                              yH=DATA$yH,yIC=DATA$yIC,
                              thEtoH=oldTh[4],thHtoIC=oldTh[5], 
                              delEtoH=thetaconst$dEtoH, delHtoIC=thetaconst$dHtoIC, 
                              detH=thetaconst$detH, detIC=thetaconst$detIC)
        print("alert: MCWM not adequate/insufficient particles")
      }
    }
    oldlp      <- lPjoint(y=old)
    newlp      <- lPjoint(y=new)
    lpoldInew  <- dmvnorm(old, mean=new, sigma=vcscaled, log=T)
    lpnewIold  <- dmvnorm(new, mean=old, sigma=vcscaled, log=T)
    logrho     <- (newlp-oldlp)+(newlL-oldLl)+(lpoldInew-lpnewIold)
    if(log(runif(1))<logrho){
      old            <- new 
      oldTh          <- newTh
      oldlE          <- newlE
      naccept[it]    <- 1
      oldLl          <- newlL
    }
    acclike[it]   <- oldLl 
    if ((it%%thin)==0){
      chain[it/thin, ] <- old
      # print(paste(it/IT*100, "%"))
    }
  }
  return(list(vcscaled, chain, naccept, acclike))
}