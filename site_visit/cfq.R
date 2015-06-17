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
fcasts  <- fcasts[!(type%in%c("rolling","rollingCtl"))] 

# merge with ifp.data
ifp.data <- ifp.data[IFPID!="1419"]
setkey(fcasts, ifp_idx)
setkey(ifp.data, IFPID)
fcasts <- ifp.data[fcasts]

names.bins <- grep("bin_", names(fcasts), fixed=TRUE, value=TRUE)
names.ctpts <- grep("ctpt_", names(fcasts), fixed=TRUE, value=TRUE)

# Convert date outcomes to numeric
fcasts[dist=="gamma",Outcome:=as.numeric(as.Date(Outcome))]
fcasts[,unique(Outcome),by=c("ifp_id")]

# Remove any foracsts without specified outcomes
fcasts <- fcasts[!is.na(Outcome),]
fcasts[,Outcome:=as.numeric(Outcome)]

# Check outcome verses medians as sanity check
fcasts[,est.median:=QDF(dist,.5,par.1,par.2),by=1:nrow(fcasts)]
fcasts[,list(median(est.median),Outcome[1]),by=ifp_id][4:36,round(.SD,3),.SDcols=c("V1","V2")]

fcasts[,list(est=round(median(est.median),digits=10),
             out=round(median(Outcome), digits=2)),by=ifp_id]


# Error
fcasts[,SS.tot:=2*apply(fcasts[,.SD,.SDcols=names.ctpts],1,var,na.rm=TRUE)]
fcasts[, SS.fit := SS.tot - sse]
fcasts[, R.sq := 1-(sse/SS.tot)]

# Drop 'unfittable'
fcasts <- fcasts[sse!=42,]
# Drop neg Rsq
fcasts <- fcasts[R.sq>0,]

# Adjust outcome (true.score) for date IFPs by date forecast made, and betas
fcasts[value_type=="date",true.score:=Outcome-as.numeric(as.Date(fcast_date))]
fcasts[value_type!="date",true.score:=Outcome]
fcasts[dist=="beta",true.score:=Outcome/100]

# Drop forecasts made after close
fcasts <- fcasts[as.Date(date_closed)>=as.Date(fcast_date)]

# FOR SOME REASON MANY DATE FORECASTS HAVE NEGATIVE CUTPOINTS... WTF? FIND OUT WHY
fcasts <- fcasts[!(dist=="gamma" & ctpt_0<0)]

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

out <- foreach(i=1:nrow(fcasts)) %dopar% {
  fcasts[i,GneitingCRPS(dist,par.1,par.2,true.score,n=1e6)]
}
fcasts[,S.R:=unlist(out)]
fcasts[is.na(S.R),.N,by=c("type")]
fcasts[, S.R.r := scaleUnit(rank(S.R)),by=ifp_id]

save.image()

tmp <- fcasts[,mean(S.R),by=c("screen_name","ifp_id","type","dist")]
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

daily.set <- foreach(ifp=unique(fcasts$ifp_id), .combine="rbind") %do% {
  getDailySet(fcasts, ifp)
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
ThetaM[dist=="beta", true.score:=true.score/100]
GneitingCRPS("beta",40,73,.447,n=1e6)

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

fcasts[ifp_id=="1410-0",]
ThetaM[ifp_id=="1410-0",]

ThetaM.MD <- ThetaM[,list(S.R=mean(S.R)),by=c("ifp_id","type")]
ThetaM.MD[,S.R.r:=rank(S.R),by=c("ifp_id")]

fm1 <- clm(ordered(S.R.r)~factor(type, levels=c("fixed","random","user")), data=ThetaM.MD)
summary(fm1)

# fixed < user
# random ~< user

save.image()


ThetaM.MD[,sd(S.R.r)/sqrt(.N),by=type]
ThetaM.MD[,mean(S.R.r),by=type]
setkey(fcasts,ifp_id)
setkey(ThetaM,ifp_id)

1-mean(ThetaM[type=="random",mean(S.R),by=ifp_id]$V1/fcasts[type=="random",mean(S.R),by=ifp_id]$V1)
1-mean(ThetaM[type=="user",mean(S.R),by=ifp_id]$V1/fcasts[type=="user",mean(S.R),by=ifp_id]$V1)
1-mean(ThetaM[type=="fixed",mean(S.R),by=ifp_id]$V1/fcasts[type=="fixed",mean(S.R),by=ifp_id]$V1)
sd(ThetaM[type=="random",mean(S.R),by=ifp_id]$V1/fcasts[type=="random",mean(S.R),by=ifp_id]$V1)/sqrt(36)
sd(ThetaM[type=="user",mean(S.R),by=ifp_id]$V1/fcasts[type=="user",mean(S.R),by=ifp_id]$V1)/sqrt(36)
sd(ThetaM[type=="fixed",mean(S.R),by=ifp_id]$V1/fcasts[type=="fixed",mean(S.R),by=ifp_id]$V1)/sqrt(36)

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Brier Scores
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
extremize <- function(dist,par.1,par.2,k=1) {
  switch(dist,
         gamma = {
           mu <- par.1/par.2
           v <- k*par.1/par.2^2
           beta <- mu/v
           alpha <- v*beta^2
           list(alpha,beta)
         },
         normal = {
           alpha <- par.1
           beta <- par.2*k
           list(alpha,beta)
         },
         beta = {
           mu <- par.1/(par.1+par.2)
           v <- k*(par.1*par.2)/((par.1+par.2)^2*(par.1+par.2+1))
           alpha <- (((1-mu)/v)-(1/mu))*mu^2
           beta <- (((1-mu)/v)-(1/mu))*mu*(1-mu)
           list(alpha,beta)
         })
}

blah <- function(dist,par.1,par.2) {
  switch(dist,
         gamma = {
           mu <- par.1/par.2
           v <- par.1/par.2^2
           v
         },
         normal = {
          par.2
         },
         beta = {
           mu <- par.1/(par.1+par.2)
           v <- (par.1*par.2)/((par.1+par.2)^2*(par.1+par.2+1))
           v
         })
}


BS <- function(dist, par.1, par.2, ctpt, true.score) {
  round(2*(CDF(dist, ctpt, par.1, par.2)-(true.score<ctpt))^2,4)
}
vec.BS <- Vectorize(BS)

# Fix date formats
ThetaM[,ctpts:=gsub("/","-",ctpts)]
#
ThetaM[dist=="beta",true.score:=true.score/100]
setkey(ThetaM,ifp_id)

# get GJP cutpoints
id <- unique(ThetaM$ifp_id)[2]
id <- "1489-6"
ThetaM[ifp_id==id,ctpts:="2015-03-01|2015-05-31"]

id <- "1493-6"

ThetaM[,unique(ctpts),by=ifp_id]


MDBS <- foreach(id=unique(ThetaM$ifp_id)) %do% {
  tryCatch({
    score.cut <- unlist(ThetaM[id,][1,strsplit(ctpts,"|",fixed=TRUE)])
    if(ThetaM[id,]$dist[1]=="gamma") {
      if(grepl("^..-",score.cut[1])) {
        score.cut <- sapply(score.cut, function(s) {
          tmp <- unlist(strsplit(s,"-",fixed="TRUE"))
          paste(paste0("20",tmp[3]),
                ifelse(nchar(tmp[1])==1,paste0("0",tmp[1]),tmp[1]),
                ifelse(nchar(tmp[2])==1,paste0("0",tmp[2]),tmp[2]),
                sep="-")
        })
      }
      score.cut <- unlist(sapply(score.cut, as.Date))
      score.cut <- foreach(sc=score.cut, .combine="cbind") %do% {
        matrix(sc-as.numeric(ThetaM[id,]$roll.date), ncol=1)
      }
    } else {
      score.cut <- matrix(as.numeric(rep(score.cut, times=nrow(ThetaM[id,])), ncol=1))
    }

    out <- apply(score.cut,2, function(ct) {
      ThetaM[id, list(dist=dist, 
                      type=type,
                      BS=vec.BS(dist, par.1, par.2, ct, true.score),
                      N=N)]
    })
    out <- out[[1]]
    out[,list(BS=round(mean(BS, na.rm=TRUE),4),
              N=mean(N)),by=c("type","ifp_id","dist")] 
  }, error = function(cond) {
    data.table(type=NA,ifp_id=id,dist=ThetaM[id,]$dist[1],BS=NA,N=NA)
  })
}
MDBS <- do.call(rbind,MDBS)
MDBS
save.image(file="roll.Rdata")


MDBS[ifp_id=="1415-0",mean(BS, na.rm=TRUE),by=c("ifp_id","type")]

ThetaM[,row.id:=1:nrow(ThetaM)]

ThetaM[,v:=blah(dist,par.1,par.2),by=row.id]
ThetaM[,sd(v),by=ifp_id]

ThetaM[,c("alpha","beta"):=extremize(dist,par.1,par.2,5),by=row.id]
MDBS.ex <- foreach(id=unique(ThetaM$ifp_id)) %do% {
  tryCatch({
    score.cut <- unlist(ThetaM[id,][1,strsplit(ctpts,"|",fixed=TRUE)])
    if(ThetaM[id,]$dist[1]=="gamma") {
      if(grepl("^..-",score.cut[1])) {
        score.cut <- sapply(score.cut, function(s) {
          tmp <- unlist(strsplit(s,"-",fixed="TRUE"))
          paste(paste0("20",tmp[3]),
                ifelse(nchar(tmp[1])==1,paste0("0",tmp[1]),tmp[1]),
                ifelse(nchar(tmp[2])==1,paste0("0",tmp[2]),tmp[2]),
                sep="-")
        })
      }
      score.cut <- unlist(sapply(score.cut, as.Date))
      score.cut <- foreach(sc=score.cut, .combine="cbind") %do% {
        matrix(sc-as.numeric(ThetaM[id,]$roll.date), ncol=1)
      }
    } else {
      score.cut <- matrix(as.numeric(rep(score.cut, times=nrow(ThetaM[id,])), ncol=1))
    }

    out <- apply(score.cut,2, function(ct) {
      ThetaM[id, list(dist=dist, 
                      type=type,
                      BS=vec.BS(dist, alpha, beta, ct, true.score),
                      N=N)]
    })
    out <- out[[1]]
    out[,list(BS=round(mean(BS),4),
              N=mean(N)),by=c("type","ifp_id","dist")] 
  }, error = function(cond) {
    data.table(type=NA,ifp_id=id,dist=ThetaM[id,]$dist[1],BS=NA,N=NA)
  })
}
MDBS.ex <- do.call(rbind,MDBS.ex)
# MDBS.ex
# round(cbind(MDBS.ex[,mean(BS),by=ifp_id]$V1,
# MDBS[,mean(BS),by=ifp_id]$V1),2)
# blah <- cbind(MDBS[,mean(BS),by=ifp_id]$V1,MDBS.ex[,mean(BS),by=ifp_id]$V1)
blah <- cbind(blah,MDBS.ex[,mean(BS),by=ifp_id]$V1)


id <- "1415-0"
setkey(daily.set,ifp_id)


MM <- foreach(id=unique(daily.set$ifp_id)) %do% {
  tryCatch({
    score.cut <- unlist(daily.set[id,][1,strsplit(ctpts,"|",fixed=TRUE)])
    if(daily.set[id,]$dist[1]=="gamma") {
      if(grepl("^..-",score.cut[1])) {
        score.cut <- sapply(score.cut, function(s) {
          tmp <- unlist(strsplit(s,"-",fixed="TRUE"))
          paste(paste0("20",tmp[3]),
                ifelse(nchar(tmp[1])==1,paste0("0",tmp[1]),tmp[1]),
                ifelse(nchar(tmp[2])==1,paste0("0",tmp[2]),tmp[2]),
                sep="-")
        })
      }
      score.cut <- unlist(sapply(score.cut, as.Date))
      score.cut <- foreach(sc=score.cut, .combine="cbind") %do% {
        matrix(sc-as.numeric(ThetaM[id,]$roll.date), ncol=1)
      }
    } else {
      score.cut <- matrix(as.numeric(rep(score.cut, times=nrow(daily.set[id,])), ncol=1))
    }

    out <- foreach(i=1:nrow(score.cut)) %dopar% {
      daily.set[id,][i,CDF(dist,score.cut[i,],par.1,par.2)]
    }

    # probabilities
    probs <- do.call(rbind,out)

    # get mean prob by day
    tmp.p <- cbind(daily.set[id,list(roll.date)],probs)
    tmp.p <- tmp.p[,lapply(.SD,mean),by=roll.date,.SDcols=grep("^V",names(tmp.p),value=TRUE)]$V1

    # get ctpt by day
    tmp.c <- cbind(daily.set[id,list(roll.date)],score.cut)[,V1[1],by=roll.date]$V1

    # outcome by day
    heaviside <- tmp.c < daily.set[id,true.score[1],by=roll.date]$V1

    mean(2*(tmp.p - heaviside)^2)



    tmp <- cbind(daily.set[id,list(type,roll.date)],out)
    cols <- grep("^V",names(tmp),value=TRUE)
    tmp <- tmp[,lapply(.SD,mean),.SDcols=cols,by=c("roll.date")]


    daily.set[i,true.score] < score.cut

    daily.set[i,true.score] < score.cut

    colMeans(do.call(rbind,out))

    cbind(daily.set[id,list(ifp_idx,type)],out)[,mean(out),by=type]



  }, error = function(cond) {
    data.table(type=NA,ifp_id=id,dist=ThetaM[id,]$dist[1],BS=NA,N=NA)
  })
}
MM <- do.call(rbind,MDBS.ex)



blah <- data.table(blah)
setnames(blah,c("k=1","k=1.5","k=3","k=5"))
blah[,ifp_id:=unique(MDBS.ex$ifp_id)]
# X11(width=6,height=6)


# p <- ggplot(data=MDBS, aes(x=ifp_id, y=BS, color=type)) +
#   geom_point() +
#   xlim(c(0,2))



ThetaM[id,qbeta(.5,par.1,par.2)]
ThetaM[id,par.1]
ThetaM[id,par.2]

#Add a "too good to be true" slide
#1412 is a great example
# have a horizontl line for GJP threshold
# line fro tru outcome
# curev w/ 95% CI for daily consensus

ifp.data



####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Fitted versus raw forecasts
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fitraw <- copy(rollers)


out <- foreach(i=1:nrow(fitraw)) %do% {
  quants <- unlist(fitraw[i,c(paste0("ctpt_",0:12)),with=FALSE])
  probs <-  unlist(fitraw[i,c(paste0("bin_",0:13)),with=FALSE])/100
  quants <- quants[1:(min(which(is.na(quants)))-1)]
  probs <- probs[1:(min(which(is.na(probs)))-1)]
  probs <- cumsum(probs)[-length(probs)]

  ts <- fitraw[i,true.score]

  # Empirical
  BS.emp <- 2*(c(probs,1-probs[length(probs)-1]) - c(ts < quants, ts>quants[length(quants)] ))^2

  # From fitted
  par.1 <- fitraw[i,]$par.1
  par.2 <- fitraw[i,]$par.2
  dist <- fitraw[i,]$dist
  f.probs <- CDF(dist, quants, par.1, par.2)
  BS.fit <- 2*(c(f.probs,1-f.probs[length(f.probs-1)]) - c(ts < quants, ts>quants[length(quants)] ))^2
  mean(BS.fit - BS.emp)
}

Reduce('+',out)/length(out)
out <- unlist(out)
mean(out)
sd(out)/sqrt(length(out))

setkey(fcasts, ifp_id)
setkey(ThetaM, ifp_id)

ThetaM[,median(S.R),by=ifp_id]$V1<fcasts[,median(S.R),by=ifp_id]$V1

i <- 1
out <- foreach(i=1:nrow(fitraw)) %do% {
  quants <- unlist(fitraw[i,c(paste0("ctpt_",0:12)),with=FALSE])
  probs <-  unlist(fitraw[i,c(paste0("bin_",0:13)),with=FALSE])/100
  quants <- quants[1:(min(which(is.na(quants)))-1)]
  probs <- probs[1:(min(which(is.na(probs)))-1)]
  probs <- cumsum(probs)[-length(probs)]

  ts <- fitraw[i,true.score]

  # Empirical
  # BS.emp <- 2*(c(probs,1-sum(probs)) - c(ts < quants, ts>quants[length(quants)] ))^2

  # From fitted
  par.1 <- fitraw[i,]$par.1
  par.2 <- fitraw[i,]$par.2
  dist <- fitraw[i,]$dist
  f.probs <- CDF(dist, quants, par.1, par.2)
  BS.fit <- 2*(c(f.probs,1-f.probs[length(f.probs)]) - c(ts < quants, ts>quants[length(quants)] ))^2
  mean(BS.fit)
}
 
out <- unlist(out)
mean(out)
fitraw

tmp <- data.table(ifp_id=fitraw$ifp_id,BS=out)
tmp <- tmp[,mean(BS),by=ifp_id]
tmp <- rbind(tmp,tmp)
setkey(tmp,ifp_id)

setkey(MDBS,ifp_id)
MDBS[,blah:=BS-tmp$V1]
MDBS[,mean(blah),by=type]
MDBS[,sd(blah)/sqrt(.N),by=type]


MDBS[,mean(BS,na.rm=TRUE)]


1-mean(ThetaM[type=="rolling",mean(S.R),by=ifp_id]$V1/fitraw[type=="rolling",mean(S.R),by=ifp_id]$V1)
1-mean(ThetaM[type=="rollingCtl",mean(S.R),by=ifp_id]$V1/fitraw[type=="rollingCtl",mean(S.R),by=ifp_id]$V1)

(Reduce('+',out)/length(out))

MDBS
do.call(c,out)

fcasts

length(out)

fitraw[,MDBS:=do.call(c,out)]
tmp <- fitraw[,mean(MDBS),by=ifp_id]


mean(MDBS[,mean(BS, na.rm=TRUE),by=ifp_id]$V1/fitraw[,mean(MDBS),by=ifp_id]$V1, na.rm=TRUE)


fitraw[MDBS>2,][1,]

(2*(.15-0)^2+
2*(.65-0)^2+
2*(1-0)^2+
2*(0-1)^2)/4

which(fitraw$idx==512)


fcasts[ifp_id=="1415-0",median((qgamma(.5,par.1,par.2)+as.numeric(as.Date(fcast_date))-Outcome))]
ThetaM[ifp_id=="1415-0",median((qgamma(.5,par.1,par.2)+as.numeric(as.Date(roll.date))-16307))]

ThetaM[ifp_id=="1415-0",]
fcasts[ifp_id=="1415-0",median(S.R),by=type]

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Lets look at 1415-0 to see wtf!
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tmp <- daily.set[ifp_id=="1415-0",]
gjp <- as.numeric(as.Date("2014-10-01"))
tmp[,gjp:=gjp-as.numeric(as.Date(roll.date))]
tmp[,row.id:=1:nrow(tmp)]

extremize <- function(dist,par.1,par.2,k=1) {
  switch(dist,
         gamma = {
           mu <- par.1/par.2
           v <- k*par.1/par.2^2
           beta <- mu/v
           alpha <- v*beta^2
           list(alpha,beta)
         })
}

k <- 3
tmp[,c("alpha","beta"):=extremize(dist, par.1,par.2,k=k), by=row.id]
avg.fcast <- tmp[,2*(mean(pgamma(gjp,par.1,par.2))-1)^2,by="roll.date"]
con.fcast <- tmp[,2*(pgamma(gjp[1],median(par.1),median(par.2))-1)^2,by="roll.date"]
con.e.fcast <- tmp[,2*(pgamma(gjp[1],median(alpha),median(beta))-1)^2,by="roll.date"]

mean(avg.fcast$V1)
mean(con.fcast$V1)
mean(con.e.fcast$V1)




####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Look at a good one to see how much worse it makes it
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tmp <- daily.set[ifp_id=="1416-0",]
gjp <- as.numeric(as.Date("2015-01-01"))
tmp[,gjp:=gjp-as.numeric(as.Date(roll.date))]
tmp[,row.id:=1:nrow(tmp)]

k <- 3
tmp[,c("alpha","beta"):=extremize(dist, par.1,par.2,k=k), by=row.id]
avg.fcast <- tmp[,2*(mean(pgamma(gjp,par.1,par.2))-1)^2,by="roll.date"]
con.fcast <- tmp[,2*(pgamma(gjp[1],median(par.1),median(par.2))-1)^2,by="roll.date"]
con.e.fcast <- tmp[,2*(pgamma(gjp[1],median(alpha),median(beta))-1)^2,by="roll.date"]

mean(avg.fcast$V1,na.rm=TRUE)
mean(con.fcast$V1, na.rm=TRUE)
mean(con.e.fcast$V1, na.rm=TRUE)

tmp2<-ThetaM[ifp_id=="1416-0",]
tmp2[,mean2*(pgamma(gjp-as.numeric(as.Date(roll.date)),par.1,par.2))^2]