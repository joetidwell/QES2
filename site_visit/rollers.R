####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ggplot2)
library(foreach)
library(doMC)
library("multicore", quietly=TRUE)
registerDoMC(max(multicore:::detectCores()-6,2)) # use all cores minus 2

source("~/ACE/global_vars.R")
source(file.path(kRootPath, "util", "load_data.R"), chdir=T)
source(file.path(kRootPath, "fitting", "interval_fitting_funcs.R"), chdir=T)
source(file.path(kRootPath, "forecast", "method_consensus_dist.R"), chdir=T)

options(stringsAsFactors = FALSE)
path.mydata <- "~/R/QES2/data"
# path.mydata <- "~/git/QES2/data"

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Get set of closed IFPS
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get closed IFPs
ifp.key <- as.data.table(read.csv(file.path(path.mydata, "ifp_key-SV.csv")))
ifp.key <- ifp.key[status=="closed",]
setkey(ifp.key, ifp_idx)
ifp.data <- as.data.table(read.csv(file.path(kRootPath, "data", "ifps", "ifps.yr4.csv")))
ifp.data <- ifp.data[, list(ifpid, q_status, q_type, date_closed)]
setkey(ifp.data, ifpid)

ifp.data <- ifp.data[ifp.key][date_closed!="NULL"]

ifp.data[,ifp_id:=paste(ifpid,q_type,sep="-")]
ifp.data <- ifp.data[q_status=="closed",]
ifp.data[, date_closed:=as.Date(date_closed)-1]

# Get resolution values
IFP.res <- data.table(read.csv(file.path(path.mydata,"lum_roll_data.csv")))
IFP.res[,IFPID:=as.character(IFPID)]
ifp.data[, ifpid:=as.character(ifpid)]
setkey(IFP.res, IFPID)
setkey(ifp.data,ifpid)

ifp.data <- IFP.res[,list(IFPID,Outcome)][ifp.data]

# Only keep IFPs which have closed, and for which we have resolution values
ifp.data <- ifp.data[Outcome!="",]

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Get all forecasts and fit distributions
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
LoadContinuous <- function(ifp.set=NULL, years = 4){

  source(file.path(kRootPath, "data","update_lumenogic_fits.R"), chdir=TRUE)

  # load forecasts
  fcasts <- data.table(read.csv(file.path(path.mydata, "fcasts.lum_roll.yr4.csv")))
  
  # rename bins
  setnames(fcasts, 
           grep("^Bin",names(fcasts),value=TRUE), 
           sub("Bin.","bin_", grep("^Bin",names(fcasts),value=TRUE)))
  # rename cutpoints
  setnames(fcasts, 
           grep("^CutPoint",names(fcasts),value=TRUE), 
           sub("CutPoint.","ctpt_", grep("^CutPoint",names(fcasts),value=TRUE)))
  # rename others
  setnames(fcasts, 
           names(fcasts)[c(1:9,23)], 
           c("lum_user_id","user_id","screen_name","type","ifp_id",
             "cifp_branch","value_type","range_low","range_high","timestamp"))

  # remove non-forecaster accounts
  fcasts <- fcasts[user_id > 20000, ]
  
  # Only keep ifps designated by ifp.set
  if(!is.null(ifp.set)) {
    fcasts <- fcasts[ifp_id %in% ifp.set,]
  }

  # add idx
  fcasts[, idx := 1:.N]
  
  # add fcast_date column from timestamps
  fcasts[, timestamp := timestamp/1000]  # convert from miliseconds to seconds
  fcasts[, date_time_gmt := as.POSIXct(timestamp, origin="1970-01-01", tz="GMT")]
  #fcasts[, fcast_date := as.Date(date_time_gmt,tz = "EST")]
  fcasts$fcast_date <- substr(as.character(fcasts$date_time_gmt),1,10)
  
  setkey(fcasts, timestamp)
  
  # Convert date-time formatted date responses to integers (days since 1970-01-01 00:00:00 GMT)
  TimeToInt <- function(date.time.tz){
    date.time.tz[is.na(date.time.tz)|date.time.tz == ""] <- NA
    date.time.tz <- as.integer(as.Date(date.time.tz))
    date.time.tz[date.time.tz < 0] <- 0
    return(date.time.tz)
  }

  # convert the cut points and ranges for date IFPs to numeric values
  cutpt.cols <- grep("^range|^ctpt",names(fcasts),value=TRUE)
  fcasts[value_type=="date",
         c(cutpt.cols):=lapply(.SD, function(x) as.character(TimeToInt(x))),
         .SDcols=cutpt.cols]

  # convert the cut points and ranges to numeric values
  fcasts[,c(cutpt.cols):=lapply(.SD, as.numeric),.SDcols=cutpt.cols]

  # Offset date bins by date of forecast
  l_ply(grep("^ctpt",names(fcasts),value=TRUE), function(x) {
    fcasts[value_type=="date", eval(as.name(x)) := eval(as.name(x)) - as.numeric(as.Date(fcast_date)), by=idx]
  })
  
  
  # add the q_type suffix to 4-digit IFPs 
  ifps      <- LoadIfps()
  ifp.types <- unique(ifps[, q_type, key = list(ifp_id = as.integer(substr(ifp_id, 0, 4)))])
  setkey(fcasts, ifp_id)
  fcasts    <- ifp.types[fcasts]
  fcasts[cifp_branch != "C1", q_type := as.integer(substr(cifp_branch, 2,2))]
  fcasts[, ifp_id:=as.character(ifp_id)]
  fcasts[, ifp_id := paste(ifp_id, q_type, sep="-")]

  # merge in q_status
  ifp.statuses <- ifps[,list(ifp_id, q_status)]
  fcasts <- ifp.statuses[fcasts]
  
  # add tournament year designation
  fcasts[, year := 4]
  
 # keep only IFPs requested via q.status
  # fcasts <- fcasts[q_status %in% q.status, ]
  
  ## Temporarily disabled updating csv
  # Load fitted parameters
  # fits <- UpdateContinuousFits()
  # setkey(fits, idx)
  # Merge with IFPs
  # setkey(fcasts, idx)
  # fcasts <- fits[fcasts]
  
  
  bin.names <- c("bin_0","bin_1","bin_2","bin_3","bin_4","bin_5","bin_6",
                 "bin_7","bin_8","bin_9","bin_10","bin_11","bin_12","bin_13")
  fcasts[,c(bin.names):=lapply(.SD, as.numeric),.SDcols=bin.names]

  
  # 'FIX' FORECASTS WITH 100% PROBABILITY IN 1 BIN
  fcasts[bin_0==100,c("bin_0","bin_1"):=list(99,1)]
  fcasts[bin_1==100,c("bin_0","bin_1","bin_2"):=list(.5,99,.5)]
  fcasts[bin_2==100,c("bin_1","bin_2","bin_3"):=list(.5,99,.5)]
  fcasts[bin_3==100,c("bin_2","bin_3","bin_4"):=list(.5,99,.5)]
  fcasts[bin_4==100,c("bin_3","bin_4","bin_5"):=list(.5,99,.5)]
  fcasts[bin_5==100,c("bin_4","bin_5","bin_6"):=list(.5,99,.5)]
  fcasts[bin_6==100,c("bin_5","bin_6","bin_7"):=list(.5,99,.5)]
  fcasts[bin_7==100,c("bin_6","bin_7","bin_8"):=list(.5,99,.5)]
  fcasts[bin_8==100,c("bin_7","bin_8","bin_9"):=list(.5,99,.5)]
  fcasts[bin_9==100,c("bin_8","bin_9","bin_10"):=list(.5,99,.5)]
  fcasts[bin_10==100,c("bin_9","bin_10","bin_11"):=list(.5,99,.5)]
  fcasts[bin_11==100,c("bin_10","bin_11","bin_12"):=list(.5,99,.5)]
  fcasts[bin_12==100,c("bin_11","bin_12","bin_13"):=list(.5,99,.5)]
  fcasts[bin_13==100,c("bin_12","bin_13"):=list(1,99)]


  # Replaces UpdateContinuousFits()
  # Refits all forecasts every time called
  fcasts <-  GetContinuousFits(fcasts)
  
  # SaveCache(fcasts, fn.name="LoadContinuousFcasts", args = args )
  return(fcasts)
}


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Generic Functions
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Switchable generic cumulative distribution function
CDF <- function(dist.name, ...) {
  return(
    switch(dist.name,
           normal = pnorm(...),
           beta   = pbeta(...),
           gamma  = pgamma(...)
    )
  )
}

# Switchable generic density function
PDF <- function(dist.name, ...) {
  return(
    switch(dist.name,
           normal = dnorm(...),
           beta   = dbeta(...),
           gamma  = dgamma(...)
    )
  )
}


# Generic Inverse Cumulative Distribution Function
QDF <- function(dist.name, ...){
  switch(dist.name,
         normal = qnorm(...),
         beta   = qbeta(...),
         gamma  = qgamma(...)
  )
}

# Generic Sampler
RDF <- function(dist.name, ...){
  switch(dist.name,
         normal = rnorm(...),
         beta   = rbeta(...),
         gamma  = rgamma(...)
  )
}

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data Wrangling
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load data
load(file.path(path.mydata,"fcasts.Rdata"))
load(file.path(path.mydata,"ifp.data.Rdata"))

load("data/fcasts.Rdata")
# Only roller conditions
rollers  <- fcasts[type%in%c("rolling","rollingCtl")] 

# merge with ifp.data
ifp.data <- ifp.data[IFPID!="1419"]
setkey(rollers, ifp_idx)
setkey(ifp.data, IFPID)
rollers <- ifp.data[rollers]

names.bins <- grep("bin_", names(fcasts), fixed=TRUE, value=TRUE)
names.ctpts <- grep("ctpt_", names(fcasts), fixed=TRUE, value=TRUE)

# Convert date outcomes to numeric
rollers[dist=="gamma",Outcome:=as.numeric(as.Date(Outcome))]
rollers[,unique(Outcome),by=c("ifp_id")]

# Remove any foracsts without specified outcomes
rollers <- rollers[!is.na(Outcome),]
rollers[,Outcome:=as.numeric(Outcome)]

# Check outcome verses medians as sanity check
rollers[,est.median:=QDF(dist,.5,par.1,par.2),by=1:nrow(rollers)]
rollers[,list(median(est.median),Outcome[1]),by=ifp_id][4:36,round(.SD,3),.SDcols=c("V1","V2")]

rollers[,list(est=round(median(est.median),digits=10),
             out=round(median(Outcome), digits=2)),by=ifp_id]


# Error
rollers[,SS.tot:=2*apply(rollers[,.SD,.SDcols=names.ctpts],1,var,na.rm=TRUE)]
rollers[, SS.fit := SS.tot - sse]
rollers[, R.sq := 1-(sse/SS.tot)]

# Drop 'unfittable'
rollers <- rollers[sse!=42,]
# Drop neg Rsq
rollers <- rollers[R.sq>0,]

# Adjust outcome (true.score) for date IFPs by date forecast made, and betas
rollers[value_type=="date",true.score:=Outcome-as.numeric(as.Date(fcast_date))]
rollers[value_type!="date",true.score:=Outcome]
rollers[dist=="beta",true.score:=Outcome/100]

# Drop forecasts made after close
rollers <- rollers[as.Date(date_closed)>=as.Date(fcast_date)]

# FOR SOME REASON MANY DATE FORECASTS HAVE NEGATIVE CUTPOINTS... WTF? FIND OUT WHY
rollers <- rollers[!(dist=="gamma" & ctpt_0<0)]

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### CRPS
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

scaleUnit <- function(x) {
  (x-min(x, na.rm=TRUE))/(max(x, na.rm=TRUE)-min(x, na.rm=TRUE))
}

GneitingCRPS <- function(distribution,par.1,par.2,true.score,n=1e6) {
  X <- RDF(distribution,n,par.1,par.2)
  X.prime <- RDF(distribution,n,par.1,par.2)
  mean(abs(X-true.score))-.5*mean(abs(X-X.prime))
}

out <- foreach(i=1:nrow(rollers)) %dopar% {
  rollers[i,GneitingCRPS(dist,par.1,par.2,true.score,n=1e6)]
}
rollers[,S.R:=unlist(out)]
rollers[is.na(S.R),.N,by=c("type")]
rollers[, S.R.r := scaleUnit(rank(S.R)),by=ifp_id]

save.image()

tmp <- rollers[,mean(S.R),by=c("screen_name","ifp_id","type","dist")]
tmp[,y:=scaleUnit(rank(V1)),by=ifp_id]

library(lme4)
mod <- lmer(y~factor(type) + (1|ifp_id) + (1|dist), data=tmp)
summary(lm(y~type,tmp))

tmp[,mean(y),by=c("type","ifp_id")]


tmp[,median(y),by=c("type","ifp_id")]
tmp[,median(y),by=c("type")]


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Consensus
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# First get 'Daily' forecasts
getDailySet <- function(fcasts, 
                        ifp,
                        decay.recent=.5, 
                        decay.min.days=3) {
  myifp <- fcasts[ifp_id==ifp,]
  bin.id <- paste0("bin_",1:13)
  ctpt.id <- paste0("ctpt_",1:12)

  min.date <- myifp[,min(as.Date(fcast_date))]
  max.date <- myifp[,as.Date(date_closed[1])]


  # Get distributions for every live day of the forecast
  roll <- foreach(d=seq.Date(min.date,max.date,"day")) %dopar% {
    # filter to forecasts made on or before d
    tmp <- myifp[as.Date(fcast_date)<=d,]

    # Keep only the most recent forecast by each user
    setkey(tmp, fcast_date)
    tmp <- tmp[idx %in% tmp[, list(last_fcast = tail(idx,n = 1)),
               by=user_id]$last_fcast]

    # get N by condition
    N <- tmp[,.N,by=type]
    # Keep only the newest forecasts
    # Either by using the decay.recent parameter and keeping the most recent X% of forecasts
    N[, n.keep.decay:=as.integer(N*decay.recent)]
    # or the number of forecasts made in the last N days as specificed by decay.min.days
    N[, n.keep.daymin:= tmp[,sum(as.Date(fcast_date)>(d-decay.min.days)),by=type]$V1]
    # whichever is larger
    N[, n.keep := max(n.keep.decay, n.keep.daymin),by=type]
    # unless N <= 5
    N[, n.keep:=ifelse(N<=5,N,n.keep)]
    N[, first:=(N-n.keep)+1]
    N[, last:=N]
    setkey(N,type)
    # print(N)
    tmp2 <- tmp[0,]
    for(condition in N$type) {
        tmp2 <- rbind(tmp2,tmp[type==condition,][N[condition,]$first:N[condition,]$last,])    
    }
    tmp2[,.N,by=type]
    tmp <- tmp2
    tmp[,roll.date:=d]
  }
  do.call(rbind,roll)
}

daily.set <- foreach(ifp=unique(rollers$ifp_id), .combine="rbind") %do% {
  getDailySet(rollers, ifp)
}

daily.set[,fcast_date:=as.Date(fcast_date)]

updateGammas <- function(daily.set, n.bins=20) {
  # Subset to forecasts that need to be updated
  # forecasts made on the same day as the roll date are fine
  # and only gammas are updated
  updates <- daily.set[dist=="gamma" & (fcast_date!=roll.date),]

  # Define refitting cutpoints for each forecast
  ctpts <- data.table(matrix(NA_real_, nrow=nrow(updates), ncol=n.bins))
  bins <- copy(ctpts)
  bin.id <- paste0("bin_",1:20)
  ctpt.id <- paste0("ctpt_",1:20)
  setnames(ctpts, ctpt.id)
  setnames(bins, bin.id)
  updates <- cbind(updates[,list(idx,ifp_id,dist,par.1,par.2,fcast_date, roll.date, type)],bins,ctpts)  
  updates[, row.id:=1:.N]  

  # Obtain probabilities for each bin
  updates[,c(bin.id):=as.list((pgamma(seq(roll.date-fcast_date, 
                                       qgamma(.9999,par.1,par.2), 
                                       length=20),par.1,par.2) - 
                                pgamma(roll.date-fcast_date, par.1, par.2)) / 
                                as.numeric(1-pgamma(roll.date-fcast_date, par.1, par.2))),
       by=row.id]

  # Convert to interval probabilities     
  updates[,c(bin.id):=as.list(updates[,bin.id,with=FALSE] - 
                              cbind(0,updates[,bin.id[-length(bin.id)],with=FALSE]))]

  # Get cutpoints
  updates[,c(ctpt.id):=as.list(seq(as.numeric(roll.date-fcast_date), 
                                   qgamma(.9999,par.1,par.2), 
                                   length=20)-as.numeric(roll.date-fcast_date)),
           by=row.id]

  # Drop unnecessary 1st ctpt/bin
  updates[,bin_1:=NULL]
  updates[,ctpt_1:=NULL]

  # refit parameters to new, renormed probabilities
  new.fits <- foreach(i=1:nrow(updates)) %dopar% {
    bin.vals <- unlist(updates[i, .SD, .SDcols=bin.id[-1]])
    ctpt.vals <- unlist(updates[i, .SD, .SDcols=ctpt.id[-1]])
    dist.name <- updates[i,]$dist
    
    tryCatch({
      if(abs(round(sum(bin.vals),2)-1)<.01) {
        fit <- FitFnInterval(probs  = bin.vals,     
                             quants = ctpt.vals,
                             dist.name = dist.name)
        return(data.table(par.1 = fit[1],
                          par.2 = fit[2],
                          sse  = fit[3]))      
      } else {
        return(data.table(par.1 = NA_real_,
                          par.2 = NA_real_,
                          sse  = NA_real_))
      }
    },
    error=function(cond){
      return(data.table(par.1 = NA_real_,
                        par.2 = NA_real_,
                        sse  = NA_real_))
    })
  }
  new.fits <- do.call(rbind,new.fits)
  daily.set[dist=="gamma" & (fcast_date!=roll.date), 
            c("par.1","par.2","sse"):=new.fits]
  return(daily.set)
}

daily.set <- updateGammas(daily.set)

#### Get Rolling Consensus
# Theta-M

ThetaM <- daily.set[,list(par.1=median(par.1, na.rm=TRUE),
                          par.2=median(par.2, na.rm=TRUE),
                          true.score=Outcome[1],
                          dist=dist[1],
                          ctpts=ctpts[1],
                          N=.N),
                    by=c("roll.date","ifp_id","type")]
ThetaM[dist=="gamma", true.score:=true.score-as.numeric(roll.date)]

# Score ThetaM
out <- foreach(i=1:nrow(ThetaM)) %dopar% {
  ThetaM[i,GneitingCRPS(dist,par.1,par.2,true.score,n=1e6)]
}

ThetaM[,S.R:=unlist(out)]
ThetaM[is.na(S.R),.N,]
ThetaM[, S.R.r := rank(S.R),by=c("ifp_id","roll.date")]

ThetaM[,sum(S.R.r==1)/(.N),by=c("type")]
ThetaM[,sum(S.R.r==2)/(.N),by=c("type")]
ThetaM[,sum(S.R.r==3)/(.N),by=c("type")]
ThetaM[,mean(S.R.r),by=c("type","ifp_id")]

ThetaM.MD <- ThetaM[,list(S.R=mean(S.R)),by=c("ifp_id","type")]
ThetaM.MD[,S.R.r:=rank(S.R),by=c("ifp_id")]

fm1 <- clm(ordered(S.R.r)~factor(type, levels=c("fixed","random","user")), data=ThetaM.MD)
summary(fm1)

# fixed < user
# random ~< user

save.image()


ThetaM.MD[,sd(S.R.r)/sqrt(.N),by=type]
ThetaM.MD[,mean(S.R.r),by=type]
