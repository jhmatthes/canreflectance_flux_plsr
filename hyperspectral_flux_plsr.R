#This code reads in the raw hyperspectral reflectance field data after they were corrected for the
#global dark reading by Sara Knox, does prelim aggregation/screening, averages daily data
#and runs PLSR models for correlation with eddy covariance flux measurements.
#
#All code is provided for fair use with no guarantees.
#Jaclyn Hatala Matthes, 2/3/15

library(pls)
library(abind)
library(R.matlab) #only necessary because flux data is in .mat format

list.dirs <- function(path=".", pattern=NULL, all.dirs=FALSE,
                      full.names=FALSE, ignore.case=FALSE) {
  
  all <- list.files(path, pattern, all.dirs,
                    full.names, recursive=FALSE, ignore.case)
  all[file.info(all)$isdir]
}

#set up directories and constants
in.path  <- "/data/" #path to hyperspectral data
fig.path <- "/figures/" #path to export figures
sites    <- c("sherman_pasture","twitchell_rice")
umol2gd <- 60*30*10^-6*12 #convert umol/(m2*s) to g/(m2*day)
plsr.samp <- 1000 #number bootstrap samples to run for PLSR

#spectrometer wavelength limits to include for analysis
wv.lm1 <- c(400,900)
wv.lm2 <- c(400,900.5)
wv.lm3 <- c(399.5,900.75)
wv.lm4 <- c(396.9,902.9)

#load rice flux data and select correct variables
#these are hard coded for our particular dataset
rice.dat <- readMat(paste(in.path,"rice_flux.mat",sep=""))
rice.dat <- rice.dat$data
rice.yr <- rice.dat[[4]]
rice.DOY <- rice.dat[[5]]
rice.time <- rice.dat[[6]]
rice.wc <- rice.dat[[103]]  #ANN method
rice.gpp <- rice.dat[[118]] #Reichstein method
rice.gpp[rice.gpp>0] <- 0

#load pasture flux data 
#these are hard coded for our particular dataset
past.dat  <- readMat(paste(in.path,"past_flux.mat",sep=""))
past.dat  <- past.dat$data
past.yr   <- past.dat[[4]]
past.DOY  <- past.dat[[5]]
past.time <- past.dat[[6]]
past.wc   <- past.dat[[113]] #ANN method
past.gpp  <- past.dat[[128]] #Reichstein method
past.gpp[past.gpp>0] <- 0

#set up storage
x <- 1 #have to do a counter since don't know a priori how many site years of data there are
wv.var <- wv.mn <- wv.ln <- matrix(nrow=1480,ncol=500) #make a huge matrix because don't know size a priori
doy <- time <- wc.daily <- wc.inst <- gpp.daily <- gpp.inst <- site.ind <- year <- vector()
gpp.weekly <- wc.weekly <- gpp.monthly <- wc.monthly <- vector()
ratio <- list()

#### LOAD, PROCESS, AND FORMAT DATA ####
#the indexing set up here will loop over the pasture, then the rice, and aggregate in that order
#the reflectance data files in this example are organized as: site > year > site_YYYY_DDD_HHMM.csv
for(s in 1:2){ 
  #loop over number of year folders in each site
  years <- as.numeric(list.files(paste(in.path,sites[s],"/",sep="")))
  for(y in 1:length(years)){
    year.files <- list.files(paste(in.path,sites[s],years[y],"/",sep=""))
    
    #loop over files with each site year
    for(f in 1:length(year.files)){ 
      
      #read in raw data file and get the DOY and time of collection
      file <- read.csv(paste(in.path,sites[s],years[y],"/",year.files[f],sep=""))
      doy[x]  <- strsplit(year.files[f],"_")[[1]][3] #get file DOY
      time[x] <- strsplit(strsplit(year.files[f],"_")[[1]][4],".csv")[[1]][1] #get file record time
      site.ind[x] <- s
      year[x] <- years[y]      
      
      #clip wavelength to the 400-900nm range due to large noise outside that range
      #changes very slightly due to differences in spectrometer calibration
      wv.tmp          <- file[,grep("wlength",colnames(file))]
      if(years[y]<2013){
        wv.clip         <- which(wv.tmp>=wv.lm1[1] & wv.tmp<=wv.lm1[2])
      } else if(years[y]==2013 & doy[x]<335){
        wv.clip         <- which(wv.tmp>=wv.lm1[1] & wv.tmp<=wv.lm1[2])
      } else if(years[y]==2013 & doy[x]>335){ 
        wv.clip         <- which(wv.tmp>=wv.lm2[1] & wv.tmp<=wv.lm2[2]) 
      } else if(years[y]==2014 & doy[x]<195) {
        wv.clip         <- which(wv.tmp>=wv.lm2[1] & wv.tmp<=wv.lm2[2]) 
      } else if(years[y]==2014 & doy[x]>195 & doy[x]<275){
        wv.clip         <- which(wv.tmp>=wv.lm3[1] & wv.tmp<=wv.lm3[2]) 
      } else if(years[y]==2014 & doy[x]>275){
        wv.clip         <- which(wv.tmp>=wv.lm4[1] & wv.tmp<=wv.lm4[2]) 
      }
      wv.ln[,x]       <- wv.tmp[wv.clip]
      
      #grab only the reflectance ratios ("ratio" columns in dataset) and 
      #normalize the measured reflectance ratio for each day to 1.0
      tmp.max       <- max(file[wv.clip,grep("ratio",colnames(file))],na.rm=TRUE)
      ratio[[x]]    <- file[wv.clip,grep("ratio",colnames(file))]/tmp.max
      
      #count how many replicates were collected
      if(is.null(ncol(ratio[[x]]))){
        repnum <- 1
      } else{
        repnum <- ncol(ratio[[x]])
      }

      if(s==1 & x==1){
        past.ratio   <- ratio[[x]]
        past.ratio.doy <- rep(doy[x],repnum)
        past.ratio.year <- rep(years[y],repnum)
      } else if(s==1 & x>1){
        past.ratio   <- cbind(past.ratio,ratio[[x]])
        past.ratio.doy <- c(past.ratio.doy,rep(doy[x],repnum))
        past.ratio.year <- c(past.ratio.year,rep(years[y],repnum))
      } else if(s==2 & !exists("rice.ratio")){
        rice.ratio  <- ratio[[x]]
        rice.ratio.doy <- rep(doy[x],repnum)
        rice.ratio.year <- rep(years[y],repnum)
      } else if(s==2 & exists("rice.ratio")){
        rice.ratio  <- cbind(rice.ratio,ratio[[x]])
        rice.ratio.doy <- c(rice.ratio.doy,rep(doy[x],repnum))
        rice.ratio.year <- c(rice.ratio.year,rep(years[y],repnum))
      }
      
      #average/get variance for ratios within each site visit
      if(length(grep("ratio",colnames(file)))>1){
        wv.var[,x]      <- apply(ratio[[x]],1,na.rm=TRUE,"sd")
        wv.mn[,x]       <- apply(ratio[[x]],1,na.rm=TRUE,"mean")
      } else{
        wv.var[,x]      <- ratio[[x]]
        wv.mn[,x]       <- ratio[[x]]
      }
      
      #grab the matching flux data for each DOY, time, and integrated flux timescale
      if(s==1){
        exact.ind <- which(past.yr==years[y] & past.DOY==as.numeric(doy)[x] & 
                             past.time==round(as.numeric(time[x])/100))
        daily.ind <- which(past.yr==years[y] & past.DOY==as.numeric(doy)[x])
        weekly.ind <- which(past.yr==years[y] & past.DOY>(as.numeric(doy)[x]-7) & past.DOY<=(as.numeric(doy)[x]))
        monthly.ind <- which(past.yr==years[y] & past.DOY>(as.numeric(doy)[x]-30) & past.DOY<=(as.numeric(doy)[x]))
        
        wc.daily[x]  <- sum(past.wc[daily.ind])*umol2gd
        wc.inst[x]   <- past.wc[exact.ind]
        wc.weekly[x] <- sum(past.wc[weekly.ind])*umol2gd
        wc.monthly[x] <- sum(past.wc[monthly.ind])*umol2gd
        gpp.daily[x] <- sum(past.gpp[daily.ind])*umol2gd
        gpp.inst[x]  <- past.gpp[exact.ind]
        gpp.weekly[x] <- sum(past.gpp[weekly.ind])*umol2gd
        gpp.monthly[x] <- sum(past.gpp[monthly.ind])*umol2gd
      } else if(s==2){
        exact.ind <- which(rice.yr==years[y] & rice.DOY==as.numeric(doy)[x] & 
                             rice.time==round(as.numeric(time[x])/100))
        daily.ind <- which(rice.yr==years[y] & rice.DOY==as.numeric(doy)[x])
        weekly.ind <- which(rice.yr==years[y] & rice.DOY>(as.numeric(doy)[x]-7) & rice.DOY<=(as.numeric(doy)[x]))
        monthly.ind <- which(rice.yr==years[y] & rice.DOY>(as.numeric(doy)[x]-30) & rice.DOY<=(as.numeric(doy)[x]))
        
        wc.daily[x]   <- sum(rice.wc[daily.ind])*umol2gd
        wc.inst[x]    <- rice.wc[exact.ind]
        wc.weekly[x]  <- sum(rice.wc[weekly.ind])*umol2gd
        wc.monthly[x] <- sum(rice.wc[monthly.ind])*umol2gd
        gpp.daily[x]  <- sum(rice.gpp[daily.ind])*umol2gd
        gpp.inst[x]   <- rice.gpp[exact.ind]
        gpp.weekly[x] <- sum(rice.gpp[weekly.ind])*umol2gd
        gpp.monthly[x] <- sum(rice.gpp[weekly.ind])*umol2gd
      }
      x <- x+1
    }
  }
}

#clean up columns full of NAs in the big matrices
x <- 1
good <- vector()
for(i in 1:ncol(wv.var)){
  na.sum <- sum(is.na(wv.var[,i])) #calculate the number of NA values in the DOY 
  #only keep the data for each DOY if less than 10% values are NA
  if(na.sum<(nrow(wv.var)/10) & !is.na(gpp.inst[i])){ 
    good[x] <- i
    x <- x+1
  }
}
wv.var <- t(wv.var[,good])
wv.mn  <- t(wv.mn[,good])
wv.ln  <- t(wv.ln[,good])
gpp.daily <- gpp.daily[good]
gpp.inst  <- gpp.inst[good]
wc.daily  <- wc.daily[good]
wc.inst   <- wc.inst[good]
wc.weekly <- wc.weekly[good]
wc.monthly <- wc.monthly[good]
gpp.weekly <- gpp.weekly[good]
gpp.monthly <- gpp.monthly[good]
site.ind  <- site.ind[good]
year      <- year[good]
doy       <- doy[good]
doy.time <- format(strptime(doy, format="%j"), format="%m-%d") 
date <- as.Date(paste(doy.time,year,sep="-"),"%m-%d-%Y")
colnames(wv.mn) <- wv.ln[1,]

#integrate mean hyperspectral reflectance to 1nm intervals
wv.int <- matrix(NA,nrow(wv.mn),length(400:900))
for(i in 1:nrow(wv.mn)){
  bad <- which(is.na(wv.mn[i,]))
  if(length(bad)==0){
    tmp <- wv.mn[i,]
  } else{
    wv.mn[i,bad] <- 0
    tmp <- wv.mn[i,]
  }
  wv.int[i,] <- spline(x=wv.ln[i,],y=tmp,xout=400:900)$y
}

#### PLSR MODELING ####
#format data frame for PLSR package
sample.val <- sample(1:length(gpp.daily),.2*length(gpp.daily),replace=FALSE) #save 20% data for validation
train   <- !(1:length(gpp.daily) %in% sample.val)

gpp.daily <- abs(gpp.daily)
gpp.inst  <- abs(gpp.inst)
gpp.weekly <- abs(gpp.weekly)
gpp.monthly <- abs(gpp.monthly)
dat.all            <- list(wv.int,gpp.daily,gpp.inst,gpp.weekly,gpp.monthly,wc.daily,wc.inst,wc.weekly,wc.monthly,train)
class(dat.all)     <- "data.frame"
names(dat.all)     <- c("wv","gppd","gppi","gppw","gppm","wcd","wci","wcw","wcm","train")
row.names(dat.all) <- 1:length(gpp.daily)

dat.val <- dat.all[dat.all$train==0,]
dat.fit <- dat.all[dat.all$train==1,]

#format pasture-only data for PLSR package
past.sample.val <- sample(1:length(gpp.daily[site.ind==1]),.2*length(gpp.daily[site.ind==1]),replace=FALSE) #save 20% data for validation
past.train   <- !(1:length(gpp.daily[site.ind==1]) %in% past.sample.val)

dat.past            <- list(wv.int[site.ind==1,],gpp.daily[site.ind==1],gpp.inst[site.ind==1],gpp.weekly[site.ind==1],gpp.monthly[site.ind==1],
                            wc.daily[site.ind==1],wc.inst[site.ind==1],wc.weekly[site.ind==1],wc.monthly[site.ind==1],past.train)
class(dat.past)     <- "data.frame"
names(dat.past)     <- c("wv","gppd","gppi","gppw","gppm","wcd","wci","wcw","wcm","train")
row.names(dat.past) <- 1:length(gpp.daily[site.ind==1])

past.dat.val <- dat.past[dat.past$train==0,]
past.dat.fit <- dat.past[dat.past$train==1,]

#format rice-only data for PLSR package
rice.sample.val <- sample(1:length(gpp.daily[site.ind==2]),.2*length(gpp.daily[site.ind==2]),replace=FALSE) #save 20% data for validation
rice.train   <- !(1:length(gpp.daily[site.ind==2]) %in% rice.sample.val)

dat.rice            <- list(wv.int[site.ind==2,],gpp.daily[site.ind==2],gpp.inst[site.ind==2],gpp.weekly[site.ind==2],gpp.monthly[site.ind==2],
                            wc.daily[site.ind==2],wc.inst[site.ind==2],wc.weekly[site.ind==2],wc.monthly[site.ind==2],rice.train)
class(dat.rice)     <- "data.frame"
names(dat.rice)     <- c("wv","gppd","gppi","gppw","gppm","wcd","wci","wcw","wcm","train")
row.names(dat.rice) <- 1:length(gpp.daily[site.ind==2])

rice.dat.val <- dat.rice[dat.rice$train==0,]
rice.dat.fit <- dat.rice[dat.rice$train==1,]

#run PLSR models for all data, pasture only, rice only 
all <- plsr_bootstrap(dat.all,dat.fit,plsr.samp,fig.path,"TABLE_PLSRstats_ALL.csv")
past <- plsr_bootstrap(dat.past,past.dat.fit,plsr.samp,fig.path,"TABLE_PLSRstats_PAST.csv")
rice <- plsr_bootstrap(dat.rice,rice.dat.fit,plsr.samp,fig.path,"TABLE_PLSRstats_RICE.csv")

#test mean PLSR model on independent validation dataset
all.val  <- plsr_validate(dat.val,all)
past.val <- plsr_validate(past.dat.val,past)
rice.val <- plsr_validate(rice.dat.val,rice)

#make validation fit stats table for export 
val.stats     <- matrix(nrow=3,ncol=16)
val.stats[1,] <- c(all.val$wci$r2,all.val$wci$int,all.val$wcd$r2,all.val$wcd$int,all.val$wcw$r2,all.val$wcw$int,
                   all.val$wcm$r2,all.val$wcm$int,
                   all.val$gppi$r2,all.val$gppi$int,all.val$gppd$r2,all.val$gppd$int,all.val$gppw$r2,all.val$gppw$int,
                   all.val$gppm$r2,all.val$gppm$int)
val.stats[2,] <- c(past.val$wci$r2,past.val$wci$int,past.val$wcd$r2,past.val$wcd$int,past.val$wcw$r2,past.val$wcw$int,
                   past.val$wcm$r2,past.val$wcm$int,
                   past.val$gppi$r2,past.val$gppi$int,past.val$gppd$r2,past.val$gppd$int,past.val$gppw$r2,past.val$gppw$int,
                   past.val$gppm$r2,past.val$gppm$int)
val.stats[3,] <- c(rice.val$wci$r2,rice.val$wci$int,rice.val$wcd$r2,rice.val$wcd$int,rice.val$wcw$r2,rice.val$wcw$int,
                   rice.val$wcm$r2,rice.val$wcm$int,
                   rice.val$gppi$r2,rice.val$gppi$int,rice.val$gppd$r2,rice.val$gppd$int,rice.val$gppw$r2,rice.val$gppw$int,
                   rice.val$gppm$r2,rice.val$gppm$int)
write.csv(val.stats,file=paste(fig.path,"PLSR_validation_stats.csv",sep=""))

####PLOT FIGURES####
#FIGURE 1: map of field sites & dates of sampling
doy.time <- format(strptime(doy, format="%j"), format="%m-%d") 
pdf(paste(fig.path,"FIG_collectiondates.pdf",sep=""),width=10,height=5)
par(mfrow=c(1,2))
plot(as.Date(doy.time[site.ind==1],"%m-%d"),year[site.ind==1],pch="x",
     ylim=c(2010,2014),ylab="Year",xlab="Collection Date",main="a) Pasture")
plot(as.Date(doy.time[site.ind==2],"%m-%d"),year[site.ind==2],pch="x",
     ylim=c(2010,2014),ylab="",xlab="Collection Date",main="b) Rice")
dev.off()

#FIGURE 2: plot gpp/wc for two sites
pdf(paste(fig.path,"FIG2_Cfluxes.pdf",sep=""),width=12,height=12)
par(mfrow=c(2,1),mar=c(5,6,4,5))
par(mar=c(3,6,2,5))
plot(date[site.ind==1],gpp.inst[site.ind==1],pch=16,ylim=c(-40,1),
     xaxt="n",xlab="",ylab=expression(bold("Instantanous GPP ("*mu~"mol "*m^"-2"*s^"-1"*")")))
points(date[site.ind==2],gpp.inst[site.ind==2],pch=3)
axis.Date(1, at = seq(date[site.ind==1][1], date[site.ind==2][63], length.out=20),
          format= "%m-%Y", las = 1)
text(min(date[site.ind==1]),1,expression(bold("a)")))
legend(min(date[site.ind==1]),-25,c("Pasture","Rice"),pch=c(16,3))

plot(date[site.ind==1],wc.daily[site.ind==1],pch=16,ylim=c(-12,8),
     xaxt="n",xlab="",ylab=expression(bold("Daily CO"[2]~"Flux (g-C "*m^"-2"*d^"-1"*")")))
points(date[site.ind==2],wc.daily[site.ind==2],pch=3)
axis.Date(1, at = seq(date[site.ind==1][1], date[site.ind==2][63], length.out=20),
          format= "%m-%Y", las = 1)
text(min(date[site.ind==1]),7,expression(bold("b)")))
legend(min(date[site.ind==1]),-5,c("Pasture","Rice"),pch=c(16,3))
dev.off()

#FIGURE 3: plot within-site variability during phenology events
pdf(paste(fig.path,"FIG3_site_phenology.pdf",sep=""),height=12,width=12)
par(mfrow=c(2,2))

#pasture 2014, DOY 100 - GREEN CANOPY
ind <- which(past.ratio.year==2014 & past.ratio.doy==100)
colors <- gray.colors(length(ind))
plot(wv.ln[1,],wv.mn[which(doy==100 & year==2014 & site.ind==1),],type="l",
     ylim=c(0,1),ylab="Reflectance",xlab="Wavelength (nm)",
     main="a) Pasture: Green pepperweed, 10 April 2014",lwd=2)
for(i in 1:length(ind)){
  lines(wv.ln[1,],past.ratio[,ind[i]],col=colors[i])
}
lines(wv.ln[1,],wv.mn[which(doy==100 & year==2014 & site.ind==1),])
legend(450,1.0,c("Daily mean",rep("Replicate",length(ind))),
       col=c(1,colors),lwd=2)

#rice 2013, DOY 212 - GREEN CANOPY
ind <- which(rice.ratio.year==2013 & rice.ratio.doy==212)
ind <- ind[1:10]
colors <- gray.colors(length(ind))
plot(wv.ln[1,],wv.mn[which(doy==212 & year==2013 & site.ind==2),],type="l",
     ylim=c(0,1),ylab="Reflectance",xlab="Wavelength (nm)",
     main="b) Rice: Green canopy, 31 July 2013",lwd=2)
for(i in 1:length(ind)){
  lines(wv.ln[1,],rice.ratio[,ind[i]],col=colors[i])
}
lines(wv.ln[1,],wv.mn[which(doy==212 & year==2013 & site.ind==2),])
legend(450,1.0,c("Daily mean",rep("Replicate",length(ind))),
       col=c(1,colors),lwd=2)

#pasture 2014, DOY 127 - WHITE FLOWERS
ind <- which(past.ratio.year==2014 & past.ratio.doy==127)
colors <- gray.colors(length(ind))
plot(wv.ln[1,],wv.mn[which(doy==127 & year==2014 & site.ind==1),],type="l",
     ylim=c(0,1),ylab="Reflectance",xlab="Wavelength (nm)",
     main="c) Pasture: White flowers, 7 May 2014",lwd=2)
for(i in 1:length(ind)){
  lines(wv.ln[1,],past.ratio[,ind[i]],col=colors[i])
}
lines(wv.ln[1,],wv.mn[which(doy==127 & year==2014 & site.ind==1),])
legend(450,1.0,c("Daily mean",rep("Replicate",length(ind))),
       col=c(1,colors),lwd=2)

#rice 2013, DOY 254 - SEEDED, DRY
ind <- which(rice.ratio.year==2013 & rice.ratio.doy==254)
ind <- ind[1:10]
colors <- gray.colors(length(ind))
plot(wv.ln[1,],wv.mn[which(doy==254 & year==2013 & site.ind==2),],type="l",
     ylim=c(0,1),ylab="Reflectance",xlab="Wavelength (nm)",
     main="d) Rice: Seeded and dry canopy, 11 September 2013",lwd=2)
for(i in 1:length(ind)){
  lines(wv.ln[1,],rice.ratio[,ind[i]],col=colors[i])
}
lines(wv.ln[1,],wv.mn[which(doy==254 & year==2013 & site.ind==2),])
legend(450,1.0,c("Daily mean",rep("Replicate",length(ind))),
       col=c(1,colors),lwd=2)
dev.off()

#FIGURE 4: plot variation in red, green, NIR reflectance for both sites over years
uniq.yrs <- unique(year)
red.ind <- which(wv.ln[1,]>688 & wv.ln[1,]<692)
grn.ind <- which(wv.ln[1,]>548 & wv.ln[1,]<552)
nir.ind <- which(wv.ln[1,]>798 & wv.ln[1,]<802)

x <- 1
grn.mn <- grn.var <- red.mn <- red.var <- nir.mn <- nir.var <- vector()
for(s in 1:2){
  for(y in 1:length(uniq.yrs)){
    ind <- which(year==uniq.yrs[y] & site.ind==s)
    for(d in 1:length(ind)){
      grn.mn[x]  <- mean(wv.mn[ind[d],grn.ind])
      grn.var[x] <- mean(wv.var[ind[d],grn.ind])
      red.mn[x]  <- mean(wv.mn[ind[d],red.ind])
      red.var[x] <- mean(wv.var[ind[d],red.ind])
      nir.mn[x]  <- mean(wv.mn[ind[d],nir.ind])
      nir.var[x] <- mean(wv.var[ind[d],nir.ind])
      x <- x+1   
    }
  }
}

grn.var[is.na(grn.var)] <- mean(grn.var,na.rm=TRUE)
plot.cols <- year
gray.cols <- rev(gray.colors(length(uniq.yrs), start = 0.1, end = 0.8, gamma = 2.2, alpha = NULL))
plot.cols[plot.cols==2010] <- gray.cols[1]
plot.cols[plot.cols==2011] <- gray.cols[2]
plot.cols[plot.cols==2012] <- gray.cols[3]
plot.cols[plot.cols==2013] <- gray.cols[4]
plot.cols[plot.cols==2014] <- gray.cols[5]

pdf(paste(fig.path,"FIG4_RGIR_refl.pdf",sep=""),width=8,height=8)
par(mfrow=c(3,2))
#GREEN
par(mar=c(4,4,2,4))
plot(as.Date(doy.time[site.ind==1],"%m-%d"), grn.mn[site.ind==1],
     ylim=c(0,0.6),pch=19, ylab="Green reflectance at 550 nm",xlab="",
     main="Pasture Reflectance",col=plot.cols[site.ind==1])
arrows(as.Date(doy.time[site.ind==1],"%m-%d"), 
       grn.mn[site.ind==1]-grn.var[site.ind==1], 
       as.Date(doy.time[site.ind==1],"%m-%d"), 
       grn.mn[site.ind==1]+grn.var[site.ind==1], length=0.05, angle=90, code=3,col=plot.cols[site.ind==1])
legend(min(as.Date(doy.time[site.ind==2],"%m-%d")),0.6,uniq.yrs,gray.cols)
text(min(as.Date(doy.time[site.ind==1],"%m-%d")),0.58,"a)",cex=1.2)

plot(as.Date(doy.time[site.ind==2],"%m-%d"), grn.mn[site.ind==2],
     ylim=c(0,0.6),pch=19, ylab="Green reflectance at 550 nm",xlab="",
     main="Rice Reflectance",col=plot.cols[site.ind==2])
arrows(as.Date(doy.time[site.ind==2],"%m-%d"), 
       grn.mn[site.ind==2]-grn.var[site.ind==2], 
       as.Date(doy.time[site.ind==2],"%m-%d"), 
       grn.mn[site.ind==2]+grn.var[site.ind==2], length=0.05, angle=90, code=3,col=plot.cols[site.ind==2])
text(min(as.Date(doy.time[site.ind==2],"%m-%d")),0.58,"b)",cex=1.2)

#RED
par(mar=c(4,4,2,4))
plot(as.Date(doy.time[site.ind==1],"%m-%d"), red.mn[site.ind==1],
     ylim=c(0,0.6),pch=19, ylab="Red reflectance at 690 nm",xlab="",
     main="",col=plot.cols[site.ind==1])
arrows(as.Date(doy.time[site.ind==1],"%m-%d"), 
       red.mn[site.ind==1]-red.var[site.ind==1], 
       as.Date(doy.time[site.ind==1],"%m-%d"), 
       red.mn[site.ind==1]+red.var[site.ind==1], length=0.05, angle=90, code=3,col=plot.cols[site.ind==1])
text(min(as.Date(doy.time[site.ind==1],"%m-%d")),0.58,"c)",cex=1.2)

par(mar=c(4,4,2,4))
plot(as.Date(doy.time[site.ind==2],"%m-%d"), red.mn[site.ind==2],
     ylim=c(0,0.6),pch=19, ylab="Red reflectance at 690 nm",xlab="",
     main="",col=plot.cols[site.ind==2])
arrows(as.Date(doy.time[site.ind==2],"%m-%d"), 
       red.mn[site.ind==2]-red.var[site.ind==2], 
       as.Date(doy.time[site.ind==2],"%m-%d"), 
       red.mn[site.ind==2]+red.var[site.ind==2], length=0.05, angle=90, code=3,col=plot.cols[site.ind==2])
text(min(as.Date(doy.time[site.ind==2],"%m-%d")),0.58,"d)",cex=1.2)

#NIR
plot(as.Date(doy.time[site.ind==1],"%m-%d"), nir.mn[site.ind==1],
     ylim=c(0,1.0),pch=19, ylab="NIR reflectance at 800 nm",xlab="",
     main="",col=plot.cols[site.ind==1])
arrows(as.Date(doy.time[site.ind==1],"%m-%d"), 
       nir.mn[site.ind==1]-nir.var[site.ind==1], 
       as.Date(doy.time[site.ind==1],"%m-%d"), 
       nir.mn[site.ind==1]+nir.var[site.ind==1], length=0.05, angle=90, code=3,col=plot.cols[site.ind==1])
text(min(as.Date(doy.time[site.ind==1],"%m-%d")),0.96,"e)",cex=1.2)

plot(as.Date(doy.time[site.ind==2],"%m-%d"), nir.mn[site.ind==2],
     ylim=c(0,1.0),pch=19, ylab="NIR reflectance at 800 nm",xlab="",
     main="",col=plot.cols[site.ind==2])
arrows(as.Date(doy.time[site.ind==2],"%m-%d"), 
       nir.mn[site.ind==2]-nir.var[site.ind==2], 
       as.Date(doy.time[site.ind==2],"%m-%d"), 
       nir.mn[site.ind==2]+nir.var[site.ind==2], length=0.05, angle=90, code=3,col=plot.cols[site.ind==2])
text(min(as.Date(doy.time[site.ind==2],"%m-%d")),0.96,"f)",cex=1.2)
dev.off()

#FIGURE 5: plot variable importance in projection (VIP)
colors <- rev(gray.colors(4, start = 0.005, end = 0.8, gamma = 1, alpha = NULL))
par(mfrow=c(2,1))
plot(400:900,all$gppi$vip.mean,type='l',ylab="VIP",
     xlab="Wavelength [nm]",lwd=1.5,col=colors[1],ylim=c(0,2))
lines(400:900,all$gppd$vip.mean,lwd=1.5,col=colors[2])
lines(400:900,all$gppw$vip.mean,lwd=1.5,col=colors[3])
lines(400:900,all$gppm$vip.mean,lwd=1.5,col=colors[4])

plot(400:900,all$wci$vip.mean,type='l',ylab="VIP",
     xlab="Wavelength [nm]",lwd=1.5,col=colors[1],ylim=c(0,2))
lines(400:900,all$wcd$vip.mean,lwd=1.5,col=colors[2])
lines(400:900,all$wcw$vip.mean,lwd=1.5,col=colors[3])
lines(400:900,all$wcm$vip.mean,lwd=1.5,col=colors[4])
legend(400,2,c("Instantaneous","Daily","Weekly","Monthly"),col=colors,lwd=1.5)

plot(400:900,past$gppi$vip.mean,type='l',ylab="VIP",
     xlab="Wavelength [nm]",lwd=1.5,col=colors[1],ylim=c(0,2))
lines(400:900,past$gppd$vip.mean,lwd=1.5,col=colors[2])
lines(400:900,past$gppw$vip.mean,lwd=1.5,col=colors[3])
lines(400:900,past$gppm$vip.mean,lwd=1.5,col=colors[4])

plot(400:900,past$wci$vip.mean,type='l',ylab="VIP",
     xlab="Wavelength [nm]",lwd=1.5,col=colors[1],ylim=c(0,2))
lines(400:900,past$wcd$vip.mean,lwd=1.5,col=colors[2])
lines(400:900,past$wcw$vip.mean,lwd=1.5,col=colors[3])
lines(400:900,past$wcm$vip.mean,lwd=1.5,col=colors[4])

plot(400:900,rice$gppi$vip.mean,type='l',ylab="VIP",
     xlab="Wavelength [nm]",lwd=1.5,col=colors[1],ylim=c(0,2))
lines(400:900,rice$gppd$vip.mean,lwd=1.5,col=colors[2])
lines(400:900,rice$gppw$vip.mean,lwd=1.5,col=colors[3])
lines(400:900,rice$gppm$vip.mean,lwd=1.5,col=colors[4])

plot(400:900,rice$wci$vip.mean,type='l',ylab="VIP",
     xlab="Wavelength [nm]",lwd=1.5,col=colors[1],ylim=c(0,2))
lines(400:900,rice$wcd$vip.mean,lwd=1.5,col=colors[2])
lines(400:900,rice$wcw$vip.mean,lwd=1.5,col=colors[3])
lines(400:900,rice$wcm$vip.mean,lwd=1.5,col=colors[4])

#FIGURE 6: testing mean PLSR models on independent validation dataset
pdf(paste(fig.path,"FIG6_independent_val.pdf",sep=""),height=16,width=12)
par(mfrow=c(4,3))
par(mar=c(4,6,3,3))
plot_plsr_validate(dat.val,all.val,1,ymin=-2,ymax=18,xmin=-2,xmax=18,"gppd")
par(mar=c(4,2,3,3))
plot_plsr_validate(past.dat.val,past.val,0,ymin=-2,ymax=12,xmin=-2,xmax=12,"gppd")
par(mar=c(4,2,3,3))
plot_plsr_validate(rice.dat.val,rice.val,0,ymin=-2,ymax=16,xmin=-2,xmax=16,"gppd")
par(mar=c(4,6,3,3))
plot_plsr_validate(dat.val,all.val,0,ymin=-2,ymax=20,xmin=-2,xmax=20,"gppi")
par(mar=c(4,2,3,3))
plot_plsr_validate(past.dat.val,past.val,0,ymin=-2,ymax=20,xmin=-2,xmax=20,"gppi")
par(mar=c(4,2,3,3))
plot_plsr_validate(rice.dat.val,rice.val,0,ymin=-2,ymax=30,xmin=-2,xmax=30,"gppi")
par(mar=c(4,6,3,3))
plot_plsr_validate(dat.val,all.val,0,ymin=-10,ymax=8,xmin=-10,xmax=8,"wcd")
par(mar=c(4,2,3,3))
plot_plsr_validate(past.dat.val,past.val,0,ymin=-5,ymax=6,xmin=-5,xmax=6,"wcd")
par(mar=c(4,2,3,3))
plot_plsr_validate(rice.dat.val,rice.val,0,ymin=-10,ymax=5,xmin=-10,xmax=5,"wcd")
par(mar=c(4,6,3,3))
plot_plsr_validate(dat.val,all.val,0,ymin=-28,ymax=10,xmin=-28,xmax=10,"wci")
par(mar=c(4,2,3,3))
plot_plsr_validate(past.dat.val,past.val,0,ymin=-16,ymax=8,xmin=-16,xmax=8,"wci")
par(mar=c(4,2,3,3))
plot_plsr_validate(rice.dat.val,rice.val,0,ymin=-28,ymax=10,xmin=-28,xmax=10,"wci")
dev.off()

#FIGURE 7: plot of change in R2 over different timescales of integrated NEE and GPP
par(mfrow=c(1,3))
par(mar=c(4,5,4,2))
plot(c(0,1,7,30),c(all.val$wci$r2,all.val$wcd$r2,all.val$wcw$r2,all.val$wcm$r2),ylim=c(0,1),
     ylab=expression(paste("R"^2," of independent validation")),xlab=expression(paste("Days of integrated CO"[2]," Flux")),cex=1.5)
points(c(0,1,7,30),c(all.val$gppi$r2,all.val$gppd$r2,all.val$gppw$r2,all.val$gppm$r2),pch=16,cex=1.5)
text(0.2,0.98,"a)",cex=1.2)
legend(20,1,c("NEE","GPP"),pch=c(1,16))

plot(c(0,1,7,30),c(past.val$wci$r2,past.val$wcd$r2,past.val$wcw$r2,past.val$wcm$r2),ylim=c(0,1),
     ylab=expression(paste("R"^2," of independent validation")),xlab=expression(paste("Days of integrated CO"[2]," Flux")),cex=1.5)
points(c(0,1,7,30),c(past.val$gppi$r2,past.val$gppd$r2,past.val$gppw$r2,past.val$gppm$r2),pch=16,cex=1.5)
text(0.2,0.98,"b)",cex=1.2)

plot(c(0,1,7,30),c(rice.val$wci$r2,rice.val$wcd$r2,rice.val$wcw$r2,rice.val$wcm$r2),ylim=c(0,1),
     ylab=expression(paste("R"^2," of independent validation")),xlab=expression(paste("Days of integrated CO"[2]," Flux")),cex=1.5)
points(c(0,1,7,30),c(rice.val$gppi$r2,rice.val$gppd$r2,rice.val$gppw$r2,rice.val$gppm$r2),pch=16,cex=1.5)
text(0.2,0.98,"c)",cex=1.2)


