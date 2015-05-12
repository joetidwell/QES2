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

IFP.res <- data.table(read.csv(file.path(path.mydata,"ifp_resolution.csv")))
setnames(IFP.res, c("ifp_id", "value", "data_type"))
IFP.res[,ifp_id:=paste0(ifp_id,"-0")]
IFP.res <- IFP.res[value!="",]
IFP.res <- IFP.res[value!="VOIDED",]
ifp.key[,ifp_idx:=paste0(ifp_idx,"-0")]
setkey(ifp.key, ifp_idx)