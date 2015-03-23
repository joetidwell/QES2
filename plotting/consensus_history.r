####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plots for historical consensus distributions...
#### Probably only works for binary IFPs (On the GJP side) at the moment
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Dependencies
packages <- c("data.table", "lme4", "arm", "ggplot2", "foreach",
              "foreach", "doMC", "xtable", "reshape2", "plyr", "dplyr")
lapply(packages, require, character.only=TRUE)
options(stringsAsFactors = FALSE)
theme_set(theme_classic())


# Globals
path.data <- "~/git/QES2/data"
path.work <- "~/git/QES2/plotting"
source("~/ACE/global_vars.R")

# Util funcs
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
# Generic Inverse Cumulative Distribution Function
QDF <- function(dist.name, ...){
  switch(dist.name,
         normal = qnorm(...),
         beta   = qbeta(...),
         gamma  = qgamma(...)
  )
}


# Load history forecasts
con.hist <- data.table(read.csv(file.path(path.data,"con.hist.csv")))
con.hist <- con.hist[,list(ifp_id,fcast.date,dist,par.1,par.2)]
setkey(con.hist, ifp_id)

# Load IFP cutpoint data
ifp.key <- as.data.table(read.csv(file.path(path.data, "ifp_key.csv")))
setkey(ifp.key, ifp_idx)

# Load IFP data
ifp.data <- as.data.table(read.csv(file.path(kRootPath, "data", "ifps", "ifps.yr4.csv")))
ifp.data <- ifp.data[, list(ifpid, q_status, q_type, date_closed)]
setkey(ifp.data, ifpid)
ifp.data <- ifp.data[ifp.key][date_closed!="NULL"]
ifp.data[,ifp_id:=paste(ifpid,q_type,sep="-")]
ifp.data <- ifp.data[q_status=="closed",]
ifp.data[, date_closed:=as.Date(date_closed)-1]
setkey(ifp.data, ifp_id)
# Load resolution values
IFP.res <- data.table(read.csv(file.path(path.data,"ifp_resolution.csv")))
setnames(IFP.res, c("ifp_id", "value", "data_type"))
IFP.res[,ifp_id:=paste0(ifp_id,"-0")]
IFP.res <- IFP.res[value!="",]
IFP.res <- IFP.res[value!="VOIDED",]
setkey(IFP.res, ifp_id)

# Merge IFP data with consensus values
con.hist <- IFP.res[con.hist]
con.hist <- con.hist[!is.na(data_type)]

# Add cutpoints
ifp.key[,ifp_idx:=paste0(ifp_idx,"-0")]
ifp.key[,ifp_id:=ifp_idx]
ifp.key[,ifp_idx:=NULL]
setkey(ifp.key, ifp_id)
con.hist <- ifp.key[,list(ctpts, ifp_id)][con.hist]

unique(con.hist$ifp_id)

# Get data type right
con.hist[data_type=="date", val.num := as.numeric(as.Date(value))]
con.hist[data_type!="date", val.num := as.numeric(value)]
con.hist[data_type=="date", cut.num := as.numeric(as.Date(ctpts, format="%m/%d/%y"))]
con.hist[data_type!="date", cut.num := as.numeric(ctpts)]
con.hist[,val.date := as.numeric(fcast.date)]

con.hist[data_type=="date", val.num := val.num-val.date]
con.hist[data_type=="date", cut.num := cut.num-val.date]
con.hist[,row.id:=1:.N]

# Add closing date
con.hist <- ifp.data[,list(ifp_id, date_closed)][con.hist]
con.hist[,date_closed:=as.numeric(as.Date(date_closed))]

# Remove forecasts after closing date
con.hist <- con.hist[(fcast.date-date_closed)<0,]


con.hist[,prob:=round(CDF(dist, cut.num, par.1, par.2),2),by=row.id]
con.hist[,BS:=2*(CDF(dist, cut.num, par.1, par.2) - (val.num<cut.num))^2,
          by=row.id]

unique(con.hist$ifp_id)

tmp <- con.hist[ifp_id=="1494-0",]
tmp
tmp[,date_closed-fcast.date]
tmp[,pgamma(val.num+2,par.1, par.2)]
tmp[,round(pgamma(cut.num, par.1, par.2),3),by=fcast.date]

as.numeric(as.Date("2015-02-03"))

con.hist <- con.hist[ifp_id!="1412-0",]

con.hist[dist=="normal",list(unique(dist),mean(BS)),by=ifp_id]

con.hist[ifp_id=="1424-0",QDF(dist, .5, par.1, par.1), by=row.id]
con.hist[ifp_id=="1424-0",]

con.hist[ifp_id=="1454-0" & fcast.date < as.numeric(as.Date("2014-11-13")),mean(BS)]



round(pgamma(.01,9,.025),2)

qgamma(.5,3,.03)

as.Date((16597+16309)/2,origin="1970-01-01")

con.hist[dist=="gamma",unique(ifp_id)]

con.norm <- con.hist[dist=="normal"]
bs.methods <- data.table(read.csv("~/Downloads/bs_methods_ifp.yr4.csv"))
setkey(bs.methods, ifp_id)

bs.methods <- bs.methods[ifp_id%in%unique(con.norm$ifp_id),]
con.norm <- con.norm[,list(score=mean(BS)),by=ifp_id]
con.norm[,method:="Consensus.ThetaM"]

norm.out <- rbind(bs.methods, con.norm)
norm.out[method!="Consensus.ThetaM", method:="GJP Tournament"]

ggplot(norm.out, aes(x=ifp_id,y=score,color=method, size=method, shape=method)) +
  geom_point() +
  scale_size_manual(name="Methods", values=c(4,2)) +
  scale_shape_manual(name="Methods", values=c(16,1)) +
  scale_color_manual(name="Methods", values=c("#666666","#BBBBBB")) +
  labs(x="IFP", y="Mean Daily Brier Score")


  lum.data <- data.table(read.csv("~/Downloads/lum.csv"))

  lum.data


closed.ifps <- data.table(read.csv(file.path(path.data,"closed_ifps.csv")))
closed.ifps

