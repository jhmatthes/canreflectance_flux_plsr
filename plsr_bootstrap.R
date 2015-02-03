#This file runs PLSR models through a bootstrapping process, selecting a different random 70% of data for PLSR fitting 
#and 30% of the data for PLSR validation each iteration.
#Right now it is hard-coded to run PLSR models for predicting daily and instantaneous GPP, and daily and instantaneous CO2 flux
#and it will fail if any of these are missing from the data frame. 
#This function exports a table of fit statistics to the out.path directory to a file called table.name
#JHM, 10/22/14

plsr_bootstrap <- function(dat.all,dat.fit,plsr.samp,out.path,table.name){
  
  #calculate initial responses with 10 components for #components vs R2 plot
  plsr.gppd.init <- plsr(gppd ~ wv, 10, data = dat.all, method = "oscorespls", subset=dat.all$train, validation = "CV")
  plsr.gppi.init <- plsr(gppi ~ wv, 10, data=dat.all, method = "oscorespls", subset=dat.all$train, validation = "CV")
  plsr.gppw.init <- plsr(gppw ~ wv, 10, data=dat.all, method = "oscorespls", subset=dat.all$train, validation = "CV")
  plsr.gppm.init <- plsr(gppm ~ wv, 10, data=dat.all, method = "oscorespls", subset=dat.all$train, validation = "CV")
  plsr.wcd.init  <- plsr(wcd ~ wv, 10, data=dat.all, method = "oscorespls", subset=dat.all$train, validation = "CV")
  plsr.wci.init  <- plsr(wci ~ wv, 10, data=dat.all, method = "oscorespls", subset=dat.all$train, validation = "CV")
  plsr.wcw.init <- plsr(wcw ~ wv, 10, data=dat.all, method = "oscorespls", subset=dat.all$train, validation = "CV")
  plsr.wcm.init <- plsr(wcm ~ wv, 10, data=dat.all, method = "oscorespls", subset=dat.all$train, validation = "CV")
  
  #find the number of PLSR components to use by minimizing PRESS
  out <- list()
  out$gppd <- out$gppi <- out$gppw <- out$gppm <- out$wcm <- out$wcw <- out$wcd <- out$wci <- list()
  out$gppd$comp <- which(plsr.gppd.init$validation$PRESS==min(plsr.gppd.init$validation$PRESS))
  out$gppi$comp <- which(plsr.gppi.init$validation$PRESS==min(plsr.gppi.init$validation$PRESS))
  out$gppw$comp <- which(plsr.gppw.init$validation$PRESS==min(plsr.gppw.init$validation$PRESS))
  out$gppm$comp <- which(plsr.gppm.init$validation$PRESS==min(plsr.gppm.init$validation$PRESS))
  
  out$wcd$comp  <- which(plsr.wcd.init$validation$PRESS==min(plsr.wcd.init$validation$PRESS))
  out$wci$comp  <- which(plsr.wci.init$validation$PRESS==min(plsr.wci.init$validation$PRESS))
  out$wcw$comp  <- which(plsr.wcw.init$validation$PRESS==min(plsr.wcw.init$validation$PRESS))
  out$wcm$comp  <- which(plsr.wcm.init$validation$PRESS==min(plsr.wcm.init$validation$PRESS))
  
  #bootstrap sample storage
  out$gppd$coefs <- out$gppi$coefs <- out$wcd$coefs <- out$wci$coefs <- matrix(nrow=(length(400:900)+1),ncol=plsr.samp)
  out$gppw$coefs <- out$gppm$coefs <- out$wcw$coefs <- out$wcm$coefs <- matrix(nrow=(length(400:900)+1),ncol=plsr.samp)
  out$gppd$vip <- out$gppi$vip <- out$wcd$vip <- out$wci$vip <- matrix(nrow=(length(400:900)),ncol=plsr.samp)
  out$gppw$vip <- out$gppm$vip <- out$wcw$vip <- out$wcm$vip <- matrix(nrow=(length(400:900)),ncol=plsr.samp)
  out$gppd$r2.fit <- out$gppi$r2.fit <- out$wcd$r2.fit <- out$wci$r2.fit <- vector(length=plsr.samp)
  out$gppw$r2.fit <- out$gppm$r2.fit <- out$wcw$r2.fit <- out$wcm$r2.fit <- vector(length=plsr.samp)
  out$gppd$r2.val <- out$gppi$r2.val <- out$wcd$r2.val <- out$wci$r2.val <- vector(length=plsr.samp)
  out$gppw$r2.val <- out$gppm$r2.val <- out$wcw$r2.val <- out$wcm$r2.val <- vector(length=plsr.samp)
  out$gppd$rmse.fit <- out$gppi$rmse.fit <- out$wcd$rmse.fit <- out$wci$rmse.fit <- vector(length=plsr.samp)
  out$gppw$rmse.fit <- out$gppm$rmse.fit <- out$wcw$rmse.fit <- out$wcm$rmse.fit <- vector(length=plsr.samp)
  out$gppd$rmse.val <- out$gppi$rmse.val <- out$wcd$rmse.val <- out$wci$rmse.val <- vector(length=plsr.samp)
  out$gppw$rmse.val <- out$gppm$rmse.val <- out$wcw$rmse.val <- out$wcm$rmse.val <- vector(length=plsr.samp)
  out$gppd$bias.val <- out$gppi$bias.val <- out$wcd$bias.val <- out$wci$bias.val <- vector(length=plsr.samp)
  out$gppw$bias.val <- out$gppm$bias.val <- out$wcw$bias.val <- out$wcm$bias.val <- vector(length=plsr.samp)
  
  out$gppd$comp.var <- matrix(nrow=plsr.samp,ncol=out$gppd$comp)
  out$gppi$comp.var <- matrix(nrow=plsr.samp,ncol=out$gppi$comp)
  out$wcd$comp.var <- matrix(nrow=plsr.samp,ncol=out$wcd$comp)
  out$wci$comp.var <- matrix(nrow=plsr.samp,ncol=out$wci$comp)
  out$gppw$comp.var <- matrix(nrow=plsr.samp,ncol=out$gppw$comp)
  out$gppm$comp.var <- matrix(nrow=plsr.samp,ncol=out$gppm$comp)
  out$wcw$comp.var <- matrix(nrow=plsr.samp,ncol=out$wcw$comp)
  out$wcm$comp.var <- matrix(nrow=plsr.samp,ncol=out$wcm$comp)
  
  
  #run PLSR model through bootstrap samples
  for(n in 1:plsr.samp){
    print(n)
    
    #parse 70% data for fitting, 30% of data for validation
    sample.tmp  <- sample(1:nrow(dat.fit),.3*nrow(dat.fit),replace=FALSE) 
    train.tmp   <- !(1:nrow(dat.fit) %in% sample.tmp)
    dat.tmp     <- cbind(dat.fit,train.tmp)
    
    #run PLSR models
    plsr.gppd <- plsr(gppd ~ wv, out$gppd$comp, data=dat.tmp, method = "oscorespls", subset=dat.tmp$train.tmp, validation = "CV")
    plsr.gppi <- plsr(gppi ~ wv, out$gppi$comp, data=dat.tmp, method = "oscorespls", subset=dat.tmp$train.tmp, validation = "CV")
    plsr.gppw <- plsr(gppw ~ wv, out$gppw$comp, data=dat.tmp, method = "oscorespls", subset=dat.tmp$train.tmp, validation = "CV")
    plsr.gppm <- plsr(gppm ~ wv, out$gppm$comp, data=dat.tmp, method = "oscorespls", subset=dat.tmp$train.tmp, validation = "CV")
    plsr.wcd  <- plsr(wcd ~ wv, out$wcd$comp, data=dat.tmp, method = "oscorespls", subset=dat.tmp$train.tmp, validation = "CV")
    plsr.wci  <- plsr(wci ~ wv, out$wci$comp, data=dat.tmp, method = "oscorespls", subset=dat.tmp$train.tmp, validation = "CV")
    plsr.wcw  <- plsr(wcw ~ wv, out$wcw$comp, data=dat.tmp, method = "oscorespls", subset=dat.tmp$train.tmp, validation = "CV")
    plsr.wcm  <- plsr(wcm ~ wv, out$wcm$comp, data=dat.tmp, method = "oscorespls", subset=dat.tmp$train.tmp, validation = "CV")
    
    #save coefficients
    out$gppd$coefs[,n] <- coef(plsr.gppd,intercept=TRUE)
    out$gppi$coefs[,n] <- coef(plsr.gppi,intercept=TRUE)
    out$gppw$coefs[,n] <- coef(plsr.gppw,intercept=TRUE)
    out$gppm$coefs[,n] <- coef(plsr.gppm,intercept=TRUE)
    out$wcd$coefs[,n]  <- coef(plsr.wcd,intercept=TRUE)
    out$wci$coefs[,n]  <- coef(plsr.wci,intercept=TRUE)
    out$wcw$coefs[,n]  <- coef(plsr.wcw,intercept=TRUE)
    out$wcm$coefs[,n]  <- coef(plsr.wcm,intercept=TRUE)
    
    #save variance explained by each component
    out$gppd$comp.var[n,] <- attributes(scores(plsr.gppd))$explvar
    out$gppi$comp.var[n,] <- attributes(scores(plsr.gppi))$explvar
    out$wcd$comp.var[n,]  <- attributes(scores(plsr.wcd))$explvar
    out$wci$comp.var[n,]  <- attributes(scores(plsr.wci))$explvar
    out$gppw$comp.var[n,] <- attributes(scores(plsr.gppw))$explvar
    out$gppm$comp.var[n,] <- attributes(scores(plsr.gppm))$explvar
    out$wcw$comp.var[n,]  <- attributes(scores(plsr.wcw))$explvar
    out$wcm$comp.var[n,]  <- attributes(scores(plsr.wcm))$explvar
    
    #save loadings for each component
    if(n==1){
      out$gppd$comp.load <- loadings(plsr.gppd)
      out$gppi$comp.load <- loadings(plsr.gppi)
      out$gppw$comp.load <- loadings(plsr.gppw)
      out$gppm$comp.load <- loadings(plsr.gppm)
      out$wcd$comp.load  <- loadings(plsr.wcd)
      out$wci$comp.load  <- loadings(plsr.wci)
      out$wcw$comp.load  <- loadings(plsr.wcw)
      out$wcm$comp.load  <- loadings(plsr.wcm)
    } else{
      out$gppd$comp.load <- abind(out$gppd$comp.load,loadings(plsr.gppd),along=3)
      out$gppi$comp.load <- abind(out$gppi$comp.load,loadings(plsr.gppi),along=3)
      out$gppw$comp.load <- abind(out$gppw$comp.load,loadings(plsr.gppw),along=3)
      out$gppm$comp.load <- abind(out$gppm$comp.load,loadings(plsr.gppm),along=3)
      out$wcd$comp.load  <- abind(out$wcd$comp.load,loadings(plsr.wcd),along=3)
      out$wci$comp.load  <- abind(out$wci$comp.load,loadings(plsr.wci),along=3)
      out$wcw$comp.load  <- abind(out$wcw$comp.load,loadings(plsr.wcw),along=3)
      out$wcm$comp.load  <- abind(out$wcm$comp.load,loadings(plsr.wcm),along=3)
    }
    
    #calculate & save variable importance of projection (VIP) 
    out$gppd$vip[,n]  <- apply(VIP(plsr.gppd),2,mean)
    out$gppi$vip[,n]  <- apply(VIP(plsr.gppi),2,mean)
    out$wcd$vip[,n]  <- apply(VIP(plsr.wcd),2,mean)
    out$wci$vip[,n]  <- apply(VIP(plsr.wci),2,mean)
    out$gppw$vip[,n]  <- apply(VIP(plsr.gppw),2,mean)
    out$gppm$vip[,n]  <- apply(VIP(plsr.gppm),2,mean)
    out$wcw$vip[,n]  <- apply(VIP(plsr.wcw),2,mean)
    out$wcm$vip[,n]  <- apply(VIP(plsr.wcm),2,mean)
    
    #save stats for fitted model
    out$gppd$r2.fit[n] <- R2(plsr.gppd,estimate="train")$val[out$gppd$comp+1]
    out$gppi$r2.fit[n] <- R2(plsr.gppi,estimate="train")$val[out$gppi$comp+1]
    out$wcd$r2.fit[n]  <- R2(plsr.wcd,estimate="train")$val[out$wcd$comp+1]
    out$wci$r2.fit[n]  <- R2(plsr.wci,estimate="train")$val[out$wci$comp+1]
    out$gppw$r2.fit[n] <- R2(plsr.gppw,estimate="train")$val[out$gppw$comp+1]
    out$gppm$r2.fit[n] <- R2(plsr.gppm,estimate="train")$val[out$gppm$comp+1]
    out$wcw$r2.fit[n] <- R2(plsr.wcw,estimate="train")$val[out$wcw$comp+1]
    out$wcm$r2.fit[n] <- R2(plsr.wcm,estimate="train")$val[out$wcm$comp+1]
    
    
    out$gppd$rmse.fit[n] <- RMSEP(plsr.gppd,estimate="train")$val[out$gppd$comp+1]
    out$gppi$rmse.fit[n] <- RMSEP(plsr.gppi,estimate="train")$val[out$gppi$comp+1]
    out$wcd$rmse.fit[n]  <- RMSEP(plsr.wcd,estimate="train")$val[out$wcd$comp+1]
    out$wci$rmse.fit[n]  <- RMSEP(plsr.wci,estimate="train")$val[out$wci$comp+1]
    out$gppw$rmse.fit[n] <- RMSEP(plsr.gppw,estimate="train")$val[out$gppw$comp+1]
    out$gppm$rmse.fit[n] <- RMSEP(plsr.gppm,estimate="train")$val[out$gppm$comp+1]
    out$wcw$rmse.fit[n] <- RMSEP(plsr.wcw,estimate="train")$val[out$wcw$comp+1]
    out$wcm$rmse.fit[n] <- RMSEP(plsr.wcm,estimate="train")$val[out$wcm$comp+1]
    
    
    #calculate predictive stats
    pred.gppd <- predict(plsr.gppd, newdata=dat.tmp[!train.tmp,], ncomp=plsr.gppd$ncomp)
    pred.gppi <- predict(plsr.gppi, newdata=dat.tmp[!train.tmp,], ncomp=plsr.gppi$ncomp)
    pred.wcd  <- predict(plsr.wcd, newdata=dat.tmp[!train.tmp,], ncomp=plsr.wcd$ncomp)
    pred.wci  <- predict(plsr.wci, newdata=dat.tmp[!train.tmp,], ncomp=plsr.wci$ncomp)
    pred.gppw <- predict(plsr.gppw, newdata=dat.tmp[!train.tmp,], ncomp=plsr.gppw$ncomp)
    pred.gppm <- predict(plsr.gppm, newdata=dat.tmp[!train.tmp,], ncomp=plsr.gppm$ncomp)
    pred.wcw <- predict(plsr.wcw, newdata=dat.tmp[!train.tmp,], ncomp=plsr.wcw$ncomp)
    pred.wcm <- predict(plsr.wcm, newdata=dat.tmp[!train.tmp,], ncomp=plsr.wcm$ncomp)
    
    #save predictive fit statistics
    out$gppd$r2.val[n] <- summary(lm(pred.gppd~dat.tmp$gppd[!train.tmp]))$adj.r.squared
    out$gppi$r2.val[n] <- summary(lm(pred.gppi~dat.tmp$gppi[!train.tmp]))$adj.r.squared
    out$wcd$r2.val[n]  <- summary(lm(pred.wcd~dat.tmp$wcd[!train.tmp]))$adj.r.squared
    out$wci$r2.val[n]  <- summary(lm(pred.wci~dat.tmp$wci[!train.tmp]))$adj.r.squared
    out$gppw$r2.val[n] <- summary(lm(pred.gppw~dat.tmp$gppw[!train.tmp]))$adj.r.squared
    out$gppm$r2.val[n] <- summary(lm(pred.gppm~dat.tmp$gppm[!train.tmp]))$adj.r.squared
    out$wcw$r2.val[n] <- summary(lm(pred.wcw~dat.tmp$wcw[!train.tmp]))$adj.r.squared
    out$wcm$r2.val[n] <- summary(lm(pred.wcm~dat.tmp$wcm[!train.tmp]))$adj.r.squared
    
    out$gppd$rmse.val[n] <- summary(lm(pred.gppd~dat.tmp$gppd[!train.tmp]))$sigma
    out$gppi$rmse.val[n] <- summary(lm(pred.gppi~dat.tmp$gppi[!train.tmp]))$sigma
    out$wcd$rmse.val[n]  <- summary(lm(pred.wcd~dat.tmp$wcd[!train.tmp]))$sigma
    out$wci$rmse.val[n]  <- summary(lm(pred.wci~dat.tmp$wci[!train.tmp]))$sigma
    out$gppw$rmse.val[n] <- summary(lm(pred.gppw~dat.tmp$gppw[!train.tmp]))$sigma
    out$gppm$rmse.val[n] <- summary(lm(pred.gppm~dat.tmp$gppm[!train.tmp]))$sigma
    out$wcw$rmse.val[n] <- summary(lm(pred.wcw~dat.tmp$wcw[!train.tmp]))$sigma
    out$wcm$rmse.val[n] <- summary(lm(pred.wcm~dat.tmp$wcm[!train.tmp]))$sigma
    
  }
  
  #get mean model coefficients with 95% confidence interval
  out$gppd$coefs.mean <- apply(out$gppd$coefs[2:502,],1,"mean")
  out$gppd$coefs.quan <- apply(out$gppd$coefs[2:502,],1,"quantile",probs=c(0.05,0.95))
  out$gppi$coefs.mean <- apply(out$gppi$coefs[2:502,],1,"mean")
  out$gppi$coefs.quan <- apply(out$gppi$coefs[2:502,],1,"quantile",probs=c(0.05,0.95))
  out$wcd$coefs.mean  <- apply(out$wcd$coefs[2:502,],1,"mean")
  out$wcd$coefs.quan  <- apply(out$wcd$coefs[2:502,],1,"quantile",probs=c(0.05,0.95))
  out$wci$coefs.mean  <- apply(out$wci$coefs[2:502,],1,"mean")
  out$wci$coefs.quan  <- apply(out$wci$coefs[2:502,],1,"quantile",probs=c(0.05,0.95))
  out$gppw$coefs.mean <- apply(out$gppw$coefs[2:502,],1,"mean")
  out$gppw$coefs.quan <- apply(out$gppw$coefs[2:502,],1,"quantile",probs=c(0.05,0.95))
  out$gppm$coefs.mean <- apply(out$gppm$coefs[2:502,],1,"mean")
  out$gppm$coefs.quan <- apply(out$gppm$coefs[2:502,],1,"quantile",probs=c(0.05,0.95))
  out$wcw$coefs.mean <- apply(out$wcw$coefs[2:502,],1,"mean")
  out$wcw$coefs.quan <- apply(out$wcw$coefs[2:502,],1,"quantile",probs=c(0.05,0.95))
  out$wcm$coefs.mean <- apply(out$wcm$coefs[2:502,],1,"mean")
  out$wcm$coefs.quan <- apply(out$wcm$coefs[2:502,],1,"quantile",probs=c(0.05,0.95))
  
  out$gppd$int.mean   <- mean(out$gppd$coefs[1,])
  out$gppi$int.mean   <- mean(out$gppi$coefs[1,])
  out$wcd$int.mean    <- mean(out$wcd$coefs[1,])
  out$wci$int.mean    <- mean(out$wci$coefs[1,])
  out$gppw$int.mean   <- mean(out$gppw$coefs[1,])
  out$gppm$int.mean   <- mean(out$gppm$coefs[1,])
  out$wcw$int.mean   <- mean(out$wcw$coefs[1,])
  out$wcm$int.mean   <- mean(out$wcm$coefs[1,])
  
  #calculate mean VIP
  out$gppd$vip.mean  <- apply(out$gppd$vip,1,mean)
  out$gppi$vip.mean  <- apply(out$gppi$vip,1,mean)
  out$wcd$vip.mean  <- apply(out$wcd$vip,1,mean)
  out$wci$vip.mean  <- apply(out$wci$vip,1,mean)
  out$gppw$vip.mean  <- apply(out$gppw$vip,1,mean)
  out$gppm$vip.mean  <- apply(out$gppm$vip,1,mean)
  out$wcw$vip.mean  <- apply(out$wcw$vip,1,mean)
  out$wcm$vip.mean  <- apply(out$wcm$vip,1,mean)
  
  #calculate mean component loadings
  out$gppd$comp.mean  <- apply(out$gppd$comp.load,c(1,2),mean)
  out$gppi$comp.mean  <- apply(out$gppi$comp.load,c(1,2),mean)
  out$wcd$comp.mean   <- apply(out$wcd$comp.load,c(1,2),mean)
  out$wci$comp.mean   <- apply(out$wci$comp.load,c(1,2),mean)
  out$gppw$comp.mean  <- apply(out$gppw$comp.load,c(1,2),mean)
  out$gppm$comp.mean  <- apply(out$gppm$comp.load,c(1,2),mean)
  out$wcw$comp.mean   <- apply(out$wcw$comp.load,c(1,2),mean)
  out$wcm$comp.mean   <- apply(out$wcm$comp.load,c(1,2),mean)

  #export fit statistics to file
  fit.stats <- c(mean(out$gppd$r2.fit),mean(out$gppd$r2.val),mean(out$gppd$rmse.fit),mean(out$gppd$rmse.val),out$gppd$comp,
                 mean(out$gppi$r2.fit),mean(out$gppi$r2.val),mean(out$gppi$rmse.fit),mean(out$gppi$rmse.val),out$gppi$comp,
                 mean(out$wcd$r2.fit),mean(out$wcd$r2.val),mean(out$wcd$rmse.fit),mean(out$wcd$rmse.val),out$wcd$comp,
                 mean(out$wci$r2.fit),mean(out$wci$r2.val),mean(out$wci$rmse.fit),mean(out$wci$rmse.val),out$wci$comp,
                 mean(out$gppw$r2.fit),mean(out$gppw$r2.val),mean(out$gppw$rmse.fit),mean(out$gppw$rmse.val),out$gppw$comp,
                 mean(out$gppm$r2.fit),mean(out$gppm$r2.val),mean(out$gppm$rmse.fit),mean(out$gppm$rmse.val),out$gppm$comp,
                 mean(out$wcw$r2.fit),mean(out$wcw$r2.val),mean(out$wcw$rmse.fit),mean(out$wcw$rmse.val),out$wcw$comp,
                 mean(out$wcm$r2.fit),mean(out$wcm$r2.val),mean(out$wcm$rmse.fit),mean(out$wcm$rmse.val),out$wcm$comp)
  fit.names <- c("GPPd.r2.fit","GPPd.r2.val","GPPd.rmse.fit","GPPd.rmse.val","GPPd.comp",
                 "GPPi.r2.fit","GPPi.r2.val","GPPi.rmse.fit","GPPi.rmse.val","GPPi.comp",
                 "WCd.r2.fit","WCd.r2.val","WCd.rmse.fit","WCd.rmse.val","WCd.comp",
                 "WCi.r2.fit","WCi.r2.val","WCi.rmse.fit","WCi.rmse.val","WCi.comp",
                 "GPPw.r2.fit","GPPw.r2.val","GPPw.rmse.fit","GPPw.rmse.val","GPPw.comp",
                 "GPPm.r2.fit","GPPm.r2.val","GPPm.rmse.fit","GPPm.rmse.val","GPPm.comp",
                 "WCw.r2.fit","WCw.r2.val","WCw.rmse.fit","WCw.rmse.val","WCw.comp",
                 "WCm.r2.fit","WCm.r2.val","WCm.rmse.fit","WCm.rmse.val","WCm.comp")
  fit.dat <- rbind(fit.names,fit.stats)
  write.csv(fit.dat,file=paste(out.path,table.name,sep=""))
  
  out
  
}

