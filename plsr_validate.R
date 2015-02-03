#This code takes the output from plsr_bootstrap.R (out) and runs the mean PLSR model on an independent validation dataset (val.data)
#JHM, 10/23/14

plsr_validate <- function(val.data,out){
  
  val <- list()
  val$gppd$pred <- val$gppi$pred <- val$wcd$pred <- val$wci$pred <- vector(length=nrow(val.data))
  val$gppw$pred <- val$gppm$pred <- val$wcw$pred <- val$wcm$pred <- vector(length=nrow(val.data))
  
  for(i in 1:nrow(val.data)){
    val$gppd$pred[i] <- val.data$wv[i,1:501] %*% out$gppd$coefs.mean + out$gppd$int.mean
    val$gppi$pred[i] <- val.data$wv[i,1:501] %*% out$gppi$coefs.mean + out$gppi$int.mean
    val$wcd$pred[i]  <- val.data$wv[i,1:501] %*% out$wcd$coefs.mean + out$wcd$int.mean
    val$wci$pred[i]  <- val.data$wv[i,1:501] %*% out$wci$coefs.mean + out$wci$int.mean
    val$gppw$pred[i] <- val.data$wv[i,1:501] %*% out$gppw$coefs.mean + out$gppw$int.mean
    val$gppm$pred[i] <- val.data$wv[i,1:501] %*% out$gppm$coefs.mean + out$gppm$int.mean
    val$wcw$pred[i] <- val.data$wv[i,1:501] %*% out$wcw$coefs.mean + out$wcw$int.mean
    val$wcm$pred[i] <- val.data$wv[i,1:501] %*% out$wcm$coefs.mean + out$wcm$int.mean
  }
  print(val.data$gppd)
  
  #gppd - calculate linear fit
  dat.tmp <- data.frame(dat.gppd=val.data$gppd,pred.gppd=val$gppd$pred)
  lm.tmp  <- lm(dat.gppd ~ pred.gppd, data=dat.tmp)
  val$gppd$r2  <- summary(lm.tmp)$adj.r.squared
  val$gppd$int <- lm.tmp$coef[1]
  val$gppd$slp <- lm.tmp$coef[2]
  
  #get linear model confidence intervals
  val$gppd$prx <- seq(min(dat.tmp$pred.gppd), max(dat.tmp$pred.gppd), 0.5)
  val$gppd$ci <- predict(lm.tmp, newdata=data.frame(pred.gppd=val$gppd$prx), interval="confidence")
  val$gppd$pi <- predict(lm.tmp, newdata=data.frame(pred.gppd=val$gppd$prx),interval="prediction")
  
  #  val$gppd$ci <- predict(lm(dat.tmp$pred.gppd ~ dat.tmp$dat.gppd), new, interval="confidence")
  #  val$gppd$pi <- predict(lm.tmp, newdata=data.frame(pred.gppd=val$gppd$prx), interval="prediction")
  
  #GPPI - calculate linear fit
  dat.tmp <- data.frame(dat.gppi=val.data$gppi, pred.gppi=val$gppi$pred)
  lm.tmp  <- lm(dat.gppi ~ pred.gppi, data=dat.tmp)
  val$gppi$r2  <- summary(lm.tmp)$adj.r.squared
  val$gppi$int <- lm.tmp$coef[1]
  val$gppi$slp <- lm.tmp$coef[2]
  
  #get linear model confidence intervals
  val$gppi$prx <- seq(min(val$gppi$pred), max(val$gppi$pred), 0.5)
  val$gppi$ci <- predict(lm.tmp, newdata=data.frame(pred.gppi=val$gppi$prx), interval="confidence")
  val$gppi$pi <- predict(lm.tmp, newdata=data.frame(pred.gppi=val$gppi$prx), interval="prediction")
  
  #WCD - calculate linear fit
  dat.tmp <- data.frame(dat.wcd=val.data$wcd, pred.wcd=val$wcd$pred)
  lm.tmp  <- lm(dat.wcd ~ pred.wcd, data=dat.tmp)
  val$wcd$r2  <- summary(lm.tmp)$adj.r.squared
  val$wcd$int <- lm.tmp$coef[1]
  val$wcd$slp <- lm.tmp$coef[2]
  
  #get linear model confidence intervals
  val$wcd$prx <- seq(min(val$wcd$pred), max(val$wcd$pred), 0.5)
  val$wcd$ci <- predict(lm.tmp, newdata=data.frame(pred.wcd=val$wcd$prx), interval="confidence")
  val$wcd$pi <- predict(lm.tmp, newdata=data.frame(pred.wcd=val$wcd$prx), interval="prediction")
  
  #WCI - calculate linear fit
  dat.tmp <- data.frame(dat.wci=val.data$wci, pred.wci=val$wci$pred)
  lm.tmp  <- lm(dat.wci ~ pred.wci,data=dat.tmp)
  val$wci$r2  <- summary(lm.tmp)$adj.r.squared
  val$wci$int <- lm.tmp$coef[1]
  val$wci$slp <- lm.tmp$coef[2]
  
  #get linear model confidence intervals
  val$wci$prx <- seq(min(val$wci$pred), max(val$wci$pred), 0.5)
  val$wci$ci <- predict(lm.tmp, newdata=data.frame(pred.wci=val$wci$prx), interval="confidence")
  val$wci$pi <- predict(lm.tmp, newdata=data.frame(pred.wci=val$wci$prx), interval="prediction")
  
  #gppw - calculate linear fit
  dat.tmp <- data.frame(dat.gppw=val.data$gppw,pred.gppw=val$gppw$pred)
  lm.tmp  <- lm(dat.gppw ~ pred.gppw, data=dat.tmp)
  val$gppw$r2  <- summary(lm.tmp)$adj.r.squared
  val$gppw$int <- lm.tmp$coef[1]
  val$gppw$slp <- lm.tmp$coef[2]
  
  #get linear model confidence intervals
  val$gppw$prx <- seq(min(dat.tmp$pred.gppw), max(dat.tmp$pred.gppw), 0.5)
  val$gppw$ci <- predict(lm.tmp, newdata=data.frame(pred.gppw=val$gppw$prx), interval="confidence")
  val$gppw$pi <- predict(lm.tmp, newdata=data.frame(pred.gppw=val$gppw$prx),interval="prediction")
  
  #gppm - calculate linear fit
  dat.tmp <- data.frame(dat.gppm=val.data$gppm,pred.gppm=val$gppm$pred)
  lm.tmp  <- lm(dat.gppm ~ pred.gppm, data=dat.tmp)
  val$gppm$r2  <- summary(lm.tmp)$adj.r.squared
  val$gppm$int <- lm.tmp$coef[1]
  val$gppm$slp <- lm.tmp$coef[2]
  
  #get linear model confidence intervals
  val$gppm$prx <- seq(min(dat.tmp$pred.gppm), max(dat.tmp$pred.gppm), 0.5)
  val$gppm$ci <- predict(lm.tmp, newdata=data.frame(pred.gppm=val$gppm$prx), interval="confidence")
  val$gppm$pi <- predict(lm.tmp, newdata=data.frame(pred.gppm=val$gppm$prx),interval="prediction")
  
  #wcw - calculate linear fit
  dat.tmp <- data.frame(dat.wcw=val.data$wcw,pred.wcw=val$wcw$pred)
  lm.tmp  <- lm(dat.wcw ~ pred.wcw, data=dat.tmp)
  val$wcw$r2  <- summary(lm.tmp)$adj.r.squared
  val$wcw$int <- lm.tmp$coef[1]
  val$wcw$slp <- lm.tmp$coef[2]
  
  #get linear model confidence intervals
  val$wcw$prx <- seq(min(dat.tmp$pred.wcw), max(dat.tmp$pred.wcw), 0.5)
  val$wcw$ci <- predict(lm.tmp, newdata=data.frame(pred.wcw=val$wcw$prx), interval="confidence")
  val$wcw$pi <- predict(lm.tmp, newdata=data.frame(pred.wcw=val$wcw$prx),interval="prediction")
  
  #wcm - calculate linear fit
  dat.tmp <- data.frame(dat.wcm=val.data$wcm,pred.wcm=val$wcm$pred)
  lm.tmp  <- lm(dat.wcm ~ pred.wcm, data=dat.tmp)
  val$wcm$r2  <- summary(lm.tmp)$adj.r.squared
  val$wcm$int <- lm.tmp$coef[1]
  val$wcm$slp <- lm.tmp$coef[2]
  
  #get linear model confidence intervals
  val$wcm$prx <- seq(min(dat.tmp$pred.wcm), max(dat.tmp$pred.wcm), 0.5)
  val$wcm$ci <- predict(lm.tmp, newdata=data.frame(pred.wcm=val$wcm$prx), interval="confidence")
  val$wcm$pi <- predict(lm.tmp, newdata=data.frame(pred.wcm=val$wcm$prx),interval="prediction")
  
  val
}

