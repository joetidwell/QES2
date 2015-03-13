library(foreach)
library(doMC)
library("multicore", quietly=TRUE)
registerDoMC(max(multicore:::detectCores()-2,2)) # use all cores minus 2

source("~/ACE/global_vars.R")
source(file.path(kRootPath, "util", "load_data.R"), chdir=T)
source(file.path(kRootPath, "fitting", "interval_fitting_funcs.R"), chdir=T)
source(file.path(kRootPath, "fitting", "continuous_scoring_funcs.R"), chdir=T)
source(file.path(kRootPath, "forecast", "method_consensus_dist.R"), chdir=T)
options(stringsAsFactors = FALSE)
path.mydata <- "~/R/QES2/data"





# fcasts <- LoadContinuousFcasts(q.status="closed")
# unique(fcasts$ifp_id)

# save(fcasts, file="~/R/fcasts.Rdata")
# load("~/R/fcasts.Rdata")


out <- apply(ifp.data, 1, function(ifp) {
  fcasts[ifp_id==ifp['ifp_id'],FcastDateConsensus(fcasts=.SD, fcast.date=as.Date(ifp['date_closed']), use.closed=TRUE), by=type]
})

tmp <- names(out[[1]])
for(i in 1:length(out)) {
  setnames(out[[i]], tmp)
}
out <- do.call(rbind,out)  
out <- out[!is.na(method),]

# write.csv(out, file="~/R/QES2/data/closed_ifps.csv")


ifp.key <- as.data.table(read.csv(file.path(path.mydata, "ifp_key.csv")))
setkey(ifp.key, ifp_idx)
ifp.data <- as.data.table(read.csv(file.path(kRootPath, "data", "ifps", "ifps.yr4.csv")))
ifp.data <- ifp.data[, list(ifpid, q_status, q_type, date_closed)]
setkey(ifp.data, ifpid)
ifp.data <- ifp.data[ifp.key][date_closed!="NULL"]
ifp.data[,ifp_id:=paste(ifpid,q_type,sep="-")]
ifp.data <- ifp.data[q_status=="closed",]
ifp.data[, date_closed:=as.Date(date_closed)-1]
IFP.res <- data.table(read.csv(file.path(path.mydata,"ifp_resolution.csv")))
setnames(IFP.res, c("ifp_id", "value", "data_type"))
IFP.res[,ifp_id:=paste0(ifp_id,"-0")]
IFP.res <- IFP.res[value!="",]
IFP.res <- IFP.res[value!="VOIDED",]

setkey(IFP.res, ifp_id)
setkey(out, ifp_id)

data.ifp <- IFP.res[out]
data.ifp <- data.ifp[!is.na(data_type)]

ifp.key[,ifp_idx:=paste0(ifp_idx,"-0")]
setkey(ifp.key, ifp_idx)

data.ifp <- ifp.key[data.ifp]

data.ifp[data_type=="date", val.num := as.numeric(as.Date(value))]
data.ifp[data_type!="date", val.num := as.numeric(value)]
data.ifp[data_type=="date", cut.num := as.numeric(as.Date(ctpts, format="%m/%d/%y"))]
data.ifp[data_type!="date", cut.num := as.numeric(ctpts)]
data.ifp[,val.date := as.numeric(as.Date(fcast.date))]
data.ifp[data_type=="date", val.num := val.num-val.date]

data.ifp[,row.id:=1:.N]

#### Adjust values
# data.ifp[ifp_idx=="1410-0", CDF(dist, cut.num, par.1, par.2),by=row.id]
data.ifp[,prob:=round(CDF(dist, cut.num, par.1, par.2),2),by=row.id]


data.ifp[,list(ifp_idx,dist,type,par.1,par.2,val.num,cut.num,prob,fcast.date)]
data.ifp[ifp_idx=="1410-0",list(ifp_idx,dist,type,par.1,par.2,val.num,cut.num,prob,fcast.date)]



data.ifp[,BS:=round((CDF(dist, cut.num, par.1, par.2)-(val.num<cut.num)^2),2),by=row.id]



#### History

con.hist <- FcastHistConsensus(use.closed = TRUE,
                               q.status   = "closed",
                               end.date   = Sys.Date())

con.hist <- consensus.hist[!is.na(fcast.date),]

write.csv(con.hist, file.path(path.mydata,"con.hist.csv"))
