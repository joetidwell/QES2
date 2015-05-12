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