#plot the output from the plsr_validate.R script (val) against actual data (dat.val) 
#for 3 sites, 4 PLSR regression variables
#JHM, 11/3/14

plot_plsr_validate <- function(dat.val,val,plot.leg,ymin,ymax,xmin,xmax,var){

  if(var=="gppd"){
    plot(val$gppd$pred,dat.val$gppd,xlab=expression(bold("Predicted Daily GPP [g-C "*m^"-2"*d^"-1"*"]")),
         ylab=expression(bold("Observed Daily GPP [g-C "*m^"-2"*d^"-1"*"]")),ylim=c(ymin,ymax),xlim=c(ymin,ymax))
    sq.tmp <- min(floor(val$gppd$pred)):max(ceiling(val$gppd$pred))
    lines(sq.tmp,sq.tmp,lty=2)
    lines(sq.tmp,(sq.tmp*val$gppd$slp + val$gppd$int))
    #  lines(val$gppd$prx,val$gppd$ci[,1])
    lines(val$gppd$prx,val$gppd$ci[,2],lty=3)
    lines(val$gppd$prx,val$gppd$ci[,3],lty=3)
    lines(val$gppd$prx,val$gppd$pi[,1])
    lines(val$gppd$prx,val$gppd$pi[,2],col="grey")
    lines(val$gppd$prx,val$gppd$pi[,3],col="grey")
    if(plot.leg==1){
      legend(xmin,ymax,c("1:1 line","Predictive fit","95% confidence interval","95% predictive interval"),
             lty=c(2,1,3,1),col=c("black","black","black","grey"))
    }
    r2 <- bquote(bold(R^2 == .(round(val$gppd$r2,2))))
    bias <- bquote(bold(Bias == .(round(val$gppd$int,2))))
    text(xmax-4,ymin+2,r2)
    text(xmax-3.8,ymin+1,bias)
    
  } else if(var=="gppi"){
    plot(val$gppi$pred,dat.val$gppi,xlab=expression(bold("Predicted Inst. GPP ["*mu~"mol "*m^"-2"*s^"-1"*"]")),
         ylab=expression(bold("Observed Inst. GPP ["*mu~"mol "*m^"-2"*s^"-1"*"]")),ylim=c(ymin,ymax),xlim=c(ymin,ymax))
    sq.tmp <- min(floor(val$gppi$pred)):max(ceiling(val$gppi$pred))
    lines(sq.tmp,sq.tmp,lty=2)
    lines(sq.tmp,(sq.tmp*val$gppi$slp + val$gppi$int))
    lines(val$gppi$prx,val$gppi$ci[,2],lty=3)
    lines(val$gppi$prx,val$gppi$ci[,3],lty=3)
    lines(val$gppi$prx,val$gppi$pi[,1])
    lines(val$gppi$prx,val$gppi$pi[,2],col="grey")
    lines(val$gppi$prx,val$gppi$pi[,3],col="grey")
    if(plot.leg==1){
      legend(xmin,ymax,c("1:1 line","Predictive fit","95% confidence interval","95% predictive interval"),
             lty=c(2,1,3,1),col=c("black","black","black","grey"))
    }
    r2 <- bquote(bold(R^2 == .(round(val$gppi$r2,2))))
    bias <- bquote(bold(Bias == .(round(val$gppi$int,2))))
    text(xmax-4,ymin+2,r2)
    text(xmax-3.8,ymin+1,bias)
    
  } else if(var=="wcd"){
    plot(val$wcd$pred,dat.val$wcd,xlab=expression(bold("Predicted Daily CO"[2]~"Flux [g-C"*m^"-2"*d^"-1"*"]")),
         ylab=expression(bold("Observed Daily CO"[2]~"Flux [g-C"*m^"-2"*d^"-1"*"]")),ylim=c(ymin,ymax),xlim=c(ymin,ymax))
    sq.tmp <- min(floor(val$wcd$pred)):max(ceiling(val$wcd$pred))
    lines(sq.tmp,sq.tmp,lty=2)
    lines(sq.tmp,(sq.tmp*val$wcd$slp + val$wcd$int))
    #  lines(val$wcd$prx,val$wcd$ci[,1])
    lines(val$wcd$prx,val$wcd$ci[,2],lty=3)
    lines(val$wcd$prx,val$wcd$ci[,3],lty=3)
    lines(val$wcd$prx,val$wcd$pi[,1])
    lines(val$wcd$prx,val$wcd$pi[,2],col="grey")
    lines(val$wcd$prx,val$wcd$pi[,3],col="grey")
    if(plot.leg==1){
      legend(xmin,ymax,c("1:1 line","Predictive fit","95% confidence interval","95% predictive interval"),
             lty=c(2,1,3,1),col=c("black","black","black","grey"))
    }
    r2 <- bquote(bold(R^2 == .(round(val$wcd$r2,2))))
    bias <- bquote(bold(Bias == .(round(val$wcd$int,2))))
    text(xmax-4,ymin+2,r2)
    text(xmax-3.8,ymin+1,bias)
  } else if(var=="wci"){
    plot(val$wci$pred,dat.val$wci,xlab=expression(bold("Predicted Inst. CO"[2]~"Flux ["*mu~"mol "*m^"-2"*s^"-1"*"]")),
         ylab=expression(bold("Observed Inst. CO"[2]~"Flux ["*mu~"mol "*m^"-2"*s^"-1"*"]")),ylim=c(ymin,ymax),xlim=c(ymin,ymax))
    sq.tmp <- min(floor(val$wci$pred)):max(ceiling(val$wci$pred))
    lines(sq.tmp,sq.tmp,lty=2)
    lines(sq.tmp,(sq.tmp*val$wci$slp + val$wci$int))
    #  lines(val$wci$prx,val$wci$ci[,1])
    lines(val$wci$prx,val$wci$ci[,2],lty=3)
    lines(val$wci$prx,val$wci$ci[,3],lty=3)
    lines(val$wci$prx,val$wci$pi[,1])
    lines(val$wci$prx,val$wci$pi[,2],col="grey")
    lines(val$wci$prx,val$wci$pi[,3],col="grey")
    if(plot.leg==1){
      legend(xmin,ymax,c("1:1 line","Predictive fit","95% confidence interval","95% predictive interval"),
             lty=c(2,1,3,1),col=c("black","black","black","grey"))
    }
    r2 <- bquote(bold(R^2 == .(round(val$wci$r2,2))))
    bias <- bquote(bold(Bias == .(round(val$wci$int,2))))
    text(xmax-4,ymin+2.5,r2)
    text(xmax-3.8,ymin+1,bias)
  } else if(var=="gppw"){
      plot(val$gppw$pred,dat.val$gppw,xlab=expression(bold("Predicted Daily GPP [g-C "*m^"-2"*d^"-1"*"]")),
           ylab=expression(bold("Observed Daily GPP [g-C "*m^"-2"*d^"-1"*"]")),ylim=c(ymin,ymax),xlim=c(ymin,ymax))
      sq.tmp <- min(floor(val$gppw$pred)):max(ceiling(val$gppw$pred))
      lines(sq.tmp,sq.tmp,lty=2)
      lines(sq.tmp,(sq.tmp*val$gppw$slp + val$gppw$int))
      #  lines(val$gppw$prx,val$gppw$ci[,1])
      lines(val$gppw$prx,val$gppw$ci[,2],lty=3)
      lines(val$gppw$prx,val$gppw$ci[,3],lty=3)
      lines(val$gppw$prx,val$gppw$pi[,1])
      lines(val$gppw$prx,val$gppw$pi[,2],col="grey")
      lines(val$gppw$prx,val$gppw$pi[,3],col="grey")
      if(plot.leg==1){
        legend(xmin,ymax,c("1:1 line","Predictive fit","95% confidence interval","95% predictive interval"),
               lty=c(2,1,3,1),col=c("black","black","black","grey"))
      }
      r2 <- bquote(bold(R^2 == .(round(val$gppw$r2,2))))
      bias <- bquote(bold(Bias == .(round(val$gppw$int,2))))
      text(xmax-4,ymin+2,r2)
      text(xmax-3.8,ymin+1,bias)
      
  } else if(var=="gppm"){
    plot(val$gppm$pred,dat.val$gppm,xlab=expression(bold("Predicted Daily GPP [g-C "*m^"-2"*d^"-1"*"]")),
         ylab=expression(bold("Observed Daily GPP [g-C "*m^"-2"*d^"-1"*"]")),ylim=c(ymin,ymax),xlim=c(ymin,ymax))
    sq.tmp <- min(floor(val$gppm$pred)):max(ceiling(val$gppm$pred))
    lines(sq.tmp,sq.tmp,lty=2)
    lines(sq.tmp,(sq.tmp*val$gppm$slp + val$gppm$int))
    #  lines(val$gppm$prx,val$gppm$ci[,1])
    lines(val$gppm$prx,val$gppm$ci[,2],lty=3)
    lines(val$gppm$prx,val$gppm$ci[,3],lty=3)
    lines(val$gppm$prx,val$gppm$pi[,1])
    lines(val$gppm$prx,val$gppm$pi[,2],col="grey")
    lines(val$gppm$prx,val$gppm$pi[,3],col="grey")
    if(plot.leg==1){
      legend(xmin,ymax,c("1:1 line","Predictive fit","95% confidence interval","95% predictive interval"),
             lty=c(2,1,3,1),col=c("black","black","black","grey"))
    }
    r2 <- bquote(bold(R^2 == .(round(val$gppm$r2,2))))
    bias <- bquote(bold(Bias == .(round(val$gppm$int,2))))
    text(xmax-4,ymin+2,r2)
    text(xmax-3.8,ymin+1,bias)
    
  } else if(var=="wcw"){
    plot(val$wcw$pred,dat.val$wcw,xlab=expression(bold("Predicted Daily GPP [g-C "*m^"-2"*d^"-1"*"]")),
         ylab=expression(bold("Observed Daily GPP [g-C "*m^"-2"*d^"-1"*"]")),ylim=c(ymin,ymax),xlim=c(ymin,ymax))
    sq.tmp <- min(floor(val$wcw$pred)):max(ceiling(val$wcw$pred))
    lines(sq.tmp,sq.tmp,lty=2)
    lines(sq.tmp,(sq.tmp*val$wcw$slp + val$wcw$int))
    #  lines(val$wcw$prx,val$wcw$ci[,1])
    lines(val$wcw$prx,val$wcw$ci[,2],lty=3)
    lines(val$wcw$prx,val$wcw$ci[,3],lty=3)
    lines(val$wcw$prx,val$wcw$pi[,1])
    lines(val$wcw$prx,val$wcw$pi[,2],col="grey")
    lines(val$wcw$prx,val$wcw$pi[,3],col="grey")
    if(plot.leg==1){
      legend(xmin,ymax,c("1:1 line","Predictive fit","95% confidence interval","95% predictive interval"),
             lty=c(2,1,3,1),col=c("black","black","black","grey"))
    }
    r2 <- bquote(bold(R^2 == .(round(val$wcw$r2,2))))
    bias <- bquote(bold(Bias == .(round(val$wcw$int,2))))
    text(xmax-4,ymin+2,r2)
    text(xmax-3.8,ymin+1,bias)
    
  } else{
      plot(val$wcm$pred,dat.val$wcm,xlab=expression(bold("Predicted Daily GPP [g-C "*m^"-2"*d^"-1"*"]")),
           ylab=expression(bold("Observed Daily GPP [g-C "*m^"-2"*d^"-1"*"]")),ylim=c(ymin,ymax),xlim=c(ymin,ymax))
      sq.tmp <- min(floor(val$wcm$pred)):max(ceiling(val$wcm$pred))
      lines(sq.tmp,sq.tmp,lty=2)
      lines(sq.tmp,(sq.tmp*val$wcm$slp + val$wcm$int))
      #  lines(val$wcm$prx,val$wcm$ci[,1])
      lines(val$wcm$prx,val$wcm$ci[,2],lty=3)
      lines(val$wcm$prx,val$wcm$ci[,3],lty=3)
      lines(val$wcm$prx,val$wcm$pi[,1])
      lines(val$wcm$prx,val$wcm$pi[,2],col="grey")
      lines(val$wcm$prx,val$wcm$pi[,3],col="grey")
      if(plot.leg==1){
        legend(xmin,ymax,c("1:1 line","Predictive fit","95% confidence interval","95% predictive interval"),
               lty=c(2,1,3,1),col=c("black","black","black","grey"))
      }
      r2 <- bquote(bold(R^2 == .(round(val$wcm$r2,2))))
      bias <- bquote(bold(Bias == .(round(val$wcm$int,2))))
      text(xmax-4,ymin+2,r2)
      text(xmax-3.8,ymin+1,bias)
      
  }
}

