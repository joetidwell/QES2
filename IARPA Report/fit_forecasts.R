####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ggplot2)
library(foreach)
library(doMC)
library("multicore", quietly=TRUE)
# registerDoMC(max(multicore:::detectCores()-2,2)) # use all cores minus 2
registerDoMC(28) # use all cores minus 2


source("~/ACE/global_vars.R")
source(file.path(kRootPath, "util", "load_data.R"), chdir=T)
source(file.path(kRootPath, "fitting", "interval_fitting_funcs.R"), chdir=T)
source(file.path(kRootPath, "forecast", "method_consensus_dist.R"), chdir=T)

options(stringsAsFactors = FALSE)
path.mydata <- "~/R/QES2/data"

# - Get closed IFPs
# - Get raw forecast data
# - Fit distributions
#   - for each ifp  
#     - get all forecasts from opening day through closing
#     - fit forecaster level distributions


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
IFP.res <- data.table(read.csv(file.path(path.mydata,"lum_ifp_data.csv")))
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
  fcasts <- data.table(read.csv(file.path(kDataPath, "lumenogic_cont", "fcasts.lum_cont.yr4.csv")))
  
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

fcasts <- LoadContinuous(ifp.set=ifp.data$IFPID)

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

save(fcasts, file=file.path(path.mydata,"fcasts.Rdata"))
save(ifp.data, file=file.path(path.mydata,"ifp.data.Rdata"))
