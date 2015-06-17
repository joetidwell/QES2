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
# path.mydata <- "~/R/QES2/data"
path.mydata <- "~/git/QES2/data"


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
ifp.key[,ifp_idx:=paste0(ifp_idx,"-0")]
setkey(ifp.key, ifp_idx)

# Add 
IFP.res
IFP.res <- rbind(IFP.res,list("1413-0","2100-01-01","date"))
IFP.res <- rbind(IFP.res,list("1428-0","2100-01-01","date"))
IFP.res <- rbind(IFP.res,list("1416-0","2100-01-01","date"))
IFP.res <- rbind(IFP.res,list("1449-0","2100-01-01","date"))
IFP.res <- rbind(IFP.res,list("1450-0","2100-01-01","date"))
IFP.res <- rbind(IFP.res,list("1478-0","2100-01-01","date"))
IFP.res <- rbind(IFP.res,list("1434-0","2100-01-01","date"))
IFP.res <- rbind(IFP.res,list("1461-0","2100-01-01","date"))
IFP.res <- rbind(IFP.res,list("1457-0","2100-01-01","date"))
IFP.res <- rbind(IFP.res,list("1457-0","2100-01-01","date"))


# fcasts <- LoadContinuousFcasts(q.status=c("active","closed"))

# save(fcasts, file=file.path(path.mydata,"fcasts.Rdata"))
load(file.path(path.mydata,"fcasts.Rdata"))
fcasts <- fcasts[q_status=="closed",]

#### Pick one gamma IFP, and do everthing through consensus...
num <- 41
unique(fcasts$ifp_id)[num]
ifp <- unique(fcasts$ifp_id)[num]
IFP.res[ifp_id==ifp,]
ctpt.GJP <- "2015-01-01"
# runAnIFP(fcasts, ifp, ctpt.GJP,"2014-09-16")
date.from <- as.Date("2014-08-21")
date.to <- as.Date("2014-08-25")

# table(fcasts[ifp_id==ifp,]$fcast_date)
out <- foreach(day=seq(date.from, date.to, by=1)) %do% {
  runAnIFP(fcasts, ifp, ctpt.GJP, day)
}
  tmp <- do.call(rbind,out)
  ggplot(tmp, aes(x=roll.date, y=prob, color=type, group=type)) +
    geom_point() +
    geom_line()

# con.hist <- tmp
con.hist <- rbind(con.hist,tmp)
# con.hist <- con.hist[1:(nrow(con.hist)-nrow(tmp)),]
unique(con.hist$ifp_id)

ggplot(con.hist, aes(x=roll.date, y=BS, color=type, group=type)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ifp_id, scale="free_x") +
  theme_classic()
ggsave("~/Desktop/consensus_history.pdf",width=11,height=8.5)

# save(con.hist, file=file.path(path.mydata,"con.hist.Rdata"))
load(file.path(path.mydata,"con.hist.Rdata"))



exclude <- c(1498,1484,1506,1507,1469,1416)
exclude <- paste0(exclude,"-0")
#### Get mean daily brier scores

unique(con.hist$ifp_id)%in%exclude

con.hist[ifp_id!="1415-0" & ifp_id!="1416-0",mean(BS),by=("type","ifp_id")]
con.hist[ifp_id=="1425-0",BS:=2*(prob-1)^2]

MDBS <- con.hist[,list(BS=mean(BS)),by=c("ifp_id","type")]
p <- ggplot(MDBS[type %in% c("fixed","random","user"),], aes(x=type, y=BS)) +
  geom_boxplot(notch=TRUE) +
  theme_classic()
gg <- ggplot_build(p)

pdata.MDBS <- MDBS[type %in% c("fixed","random","user"),list(median=median(BS), se=1.253*sd(BS)/sqrt(.N)),by=type]

pdata.MDBS <- MDBS[type %in% c("fixed","random","user"),list(median=mean(BS), se=sd(BS)/sqrt(.N)),by=type]


pdata.MDBS[,type:=factor(type, levels=c("fixed","random","user"), labels=c("Fixed","Random","User"))]

ggplot(pdata.MDBS, aes(x=type, y=median)) +
  geom_point() +
  geom_errorbar(aes(ymax=median+se, ymin=median-se), width=.2) +
  theme_classic() +
  # ylim(c(.15,.6)) +
  labs(x="Experimental Condition", y="Median BS (w/ SE)")


bs.methods <- data.table(read.csv("~/ACE/data/scoring/bs_methods_ifp.yr4.csv"))
setkey(bs.methods, ifp_id)
bs.methods <- bs.methods[ifp_id%in%unique(con.hist$ifp_id),]
bs.methods[,type:="GJP Aggregation"]

MDBS[,score:=BS]
MDBS[,method:=type]

# pdata <- rbind(MDBS[type %in% c("fixed","random","user"),list(ifp_id,type,score,method)], bs.methods[,list(ifp_id,type,score,method)])
# pdata[, type:=factor(type, levels=c("fixed","random","user","GJP Aggregation"), labels=c("Theta-M: Fixed", "Theta-M: Random", "Theta-M: User","GJP Aggregation"))]

pdata <- rbind(MDBS[,list(ifp_id,type,score,method)], bs.methods[,list(ifp_id,type,score,method)])
pdata[, type:=factor(type, levels=c("fixed","random","user","rolling","rollingCtl","GJP Aggregation"), labels=c("Theta-M: Fixed", "Theta-M: Random", "Theta-M: User", "Theta-M: Rolling", "Theta-M: Rolling Ctl", "GJP Aggregation"))]

unique(pdata$type)

pdata <- pdata[rev(order(type)),]

ggplot(pdata, aes(x=ifp_id,y=score, color=type, size=type, shape=type)) +
  geom_point() +
  scale_size_manual(name="Methods", values=c(4,4,4,4,4,2)) +
  scale_shape_manual(name="Methods", values=c(16,16,16,16,16,1)) +
  scale_color_manual(name="Methods", values=c("firebrick","steelblue","forestgreen","darkviolet","tomato2","#BBBBBB")) +
  labs(x="IFP", y="Mean Daily Brier Score") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.justification=c(1,1), legend.position=c(1,1))

pdata2 <- pdata[,list(score=median(score)),by=c("method","type")]
ggplot(pdata2, aes(x=type, y=score, color=type, size=type, shape=type)) +
  geom_point()

ggplot(bs.methods, aes(x=score)) +
  geom_histogram(fill="#CCCCCC",color="black") +
  coord_cartesian(xlim=c(0,1),ylim=c(0,135)) +
  geom_vline(xintercept=mean(con.hist$BS), color="firebrick", size=2)

ggplot(bs.methods, aes(x=score)) +
  geom_histogram(fill="#CCCCCC",color="black") +
  coord_cartesian(xlim=c(0,1),ylim=c(0,135)) +
  geom_vline(xintercept=mean(con.hist[,mean(BS),by=type]$V1), color="firebrick", size=2)
unique(bs.methods$method)

mean(bs.methods$score)

ggplot(bs.methods, aes(x=score)) +
  geom_vline(data=con.hist[,list(condition=mean(BS)),by=type],
             aes(xintercept=condition, color=type)) +
  stat_ecdf() +
  coord_cartesian(xlim=c(0,1.3),ylim=c(0,1.1)) 

ggplot(bs.methods, aes(x=score)) +
  geom_vline(xintercept=mean(con.hist$BS)) +
  stat_ecdf() +
  coord_cartesian(xlim=c(0,1.3),ylim=c(0,1.1)) 


roll.date <- date.from + 10
runAnIFP <- function(fcasts, ifp, ctpt.GJP, roll.date) {
  tmp.date <- as.numeric(as.Date(ctpt.GJP))
  #   fcasts <- fcasts[ifp_id==ifp,]
  bin.id <- paste0("bin_",1:13)
  ctpt.id <- paste0("ctpt_",1:12)
#   f.tmp <- fcasts[5,]
#   with(f.tmp, {
#     pgamma(ctpt_2, par.1, par.2)-pgamma(ctpt_1, par.1, par.2) - bin_2/100
#   })
  decay.recent    <- .5
  decay.min.days  <- 3
  outcome <- as.numeric(as.Date(IFP.res[ifp_id==ifp,]$value))
  closed <- as.numeric(as.Date(ifp.data[ifp_id==ifp,]$date_closed))
  roll.date <- as.numeric(as.Date(roll.date))
  
  myifp <- fcasts[ifp_id==ifp,]
  
  
  # Get initial 'Day 0' parameters
  myifp[,fcast_date:=as.numeric(as.Date(fcast_date))]
  
#   hist(myifp$fcast_date)
  
  # Eliminate forecasts made prior to roll.date
  myifp <- myifp[fcast_date<=roll.date,]
  
  # Keep only the most recent forecast by each user
  setkey(myifp, fcast_date)
  myifp <- myifp[idx %in%
               myifp[, list(last_fcast = tail(idx,n = 1)),
                   by=user_id]$last_fcast]
  
  #### Test BS, no fitting
  
#   myifp[,days.till:=outcome-fcast_date]
#   hist(myifp$days.till)
  
#   OrdBS <- function(q1,q2,q3,outcome,dist,par.1,par.2) {
#     sum((CDF(dist, c(q1,q2,q3), par.1,par.2) - (outcome < c(q1, q2, q3)))^2)/3
#   }
  
  BS <- function(prob,ctpt,outcome) {
    2*(prob-(outcome<ctpt))^2
  }
  
#   EmpBS <- function(q1,q2,q3,p1,p2,p3,p4,outcome) {
#     probs <- c(p1,p2,p3)/100
#     probs <- c(cumsum(probs),1)
#     cuts <- c(q1,q2,q3,Inf)
#     print(probs)
#     print(cuts)
#     print(outcome)
#     mean((probs-(outcome<cuts))^2)
#     #   sum((c(p1,p2,p3,1-(sum(p1,p2,p3))) - (outcome < c(q1, q2, q3, Inf)))^2)/3
#   }
#   
#   myifp[,OBS:=EmpBS(ctpt_0,ctpt_1,ctpt_2,bin_0,bin_1,bin_2,bin_3,days.till),by=idx]
#   myifp[,BS:=BS(pgamma(tmp.date-fcast_date, par.1, par.2),tmp.date-fcast_date,days.till),by=idx]
  
#   
#   hist(myifp$OBS)
#   hist(myifp$BS)
  
  
  ## Looks good so far...
#   myifp[,round(median(pgamma(tmp.date-fcast_date, par.1, par.2)),2)]
  
  
  # get N by condition
  N <- myifp[,.N,by=type]
  # Keep only the newest forecasts
  # Either by using the decay.recent parameter and keeping the most recent X% of forecasts
  N[, n.keep.decay:=as.integer(N*decay.recent)]
  # or the number of forecasts made in the last N days as specificed by decay.min.days
  N[, n.keep.daymin:=myifp[,sum(fcast_date > (outcome - decay.min.days)), by=type]$V1]
  # whichever is larger
  N[, n.keep := max(n.keep.decay, n.keep.daymin),by=type]
  # unless N <= 5
  N[, n.keep:=ifelse(N<=5,N,n.keep)]
  N[, first:=(N-n.keep)+1]
  N[, last:=N]
  setkey(N,type)
  
  tmp <- myifp[0,]
  for(condition in N$type) {
      tmp <- rbind(tmp,myifp[type==condition,][N[condition,]$first:N[condition,]$last,])    
  }
  myifp <- tmp
  rm(tmp)
  
#   f.tmp <- myifp[2,]
#   with(f.tmp, {
#     round(pgamma(ctpt_2, par.1, par.2)-pgamma(ctpt_1, par.1, par.2) - bin_2/100,4)
#   })
  
  # Define refitting cutpoints for each forecast
  n.bins <- 20
  ctpts <- data.table(matrix(NA_real_, nrow=nrow(myifp), ncol=n.bins))
  bins <- copy(ctpts)
  bin.id <- paste0("bin_",1:20)
  ctpt.id <- paste0("ctpt_",1:20)
  setnames(ctpts, ctpt.id)
  setnames(bins, bin.id)
  myifp <- cbind(myifp[,list(idx,dist,par.1,par.2,fcast_date, type)],bins,ctpts)
  
  # Obtain probabilities for each bin
  myifp[,c(ctpt.id):=as.list(seq(roll.date-fcast_date, 
                               qgamma(.9999,par.1,par.2), 
                               length=20)),
       by=idx]
  myifp[,c(bin.id):=as.list(pgamma(seq(roll.date-fcast_date, 
                                     qgamma(.9999,par.1,par.2), 
                                     length=20),par.1,par.2)),
       by=idx]
  myifp[,c(bin.id):=as.list(myifp[,bin.id,with=FALSE] - cbind(0,myifp[,bin.id[-length(bin.id)],with=FALSE]))]
  myifp[,c(ctpt.id):=as.list(myifp[,ctpt.id,with=FALSE] - myifp$ctpt_1)]
  
  # Renorm probabilities
  myifp[,denom:=1-pgamma(roll.date-fcast_date,par.1,par.2)]
  myifp[,c(paste0("bin_",1:20),"denom"):=as.list(.SD/denom),.SDcols=c(paste0("bin_",1:20),"denom"),by=idx]
  
  new.fits <- foreach(id=myifp$idx) %dopar% {
    bin.vals <- unlist(myifp[idx==id, .SD, .SDcols=bin.id])
    ctpt.vals <- unlist(myifp[idx==id, .SD, .SDcols=ctpt.id])
    dist.name <- myifp[idx==id,]$dist
  #   if(dist.name=="beta") {
  #     ctpts <- na.omit(ctpts/100)
  #   }
    
    fit <- FitFnInterval(probs  = bin.vals,
                         quants = ctpt.vals,
                         dist.name = dist.name)
    
    return(data.table(idx = id,
                      par.1 = fit[1],
                      par.2 = fit[2],
                      sse  = fit[3]))
  }
  
#   new.fits
  fits <- data.table(do.call(rbind, new.fits))
  setkey(fits, idx)
  setkey(myifp,idx)
  
#   myifp[,list(median(par.1),median(par.2)),by=type][,qgamma(.5,V1,V2)]
#   myifp[,round(qgamma(.5,par.1,par.2),2),by=type]
  myifp.rolled <- myifp[,list(idx,type,fcast_date)][fits]
  myifp.rolled.con <- myifp.rolled[,list(par.1=median(par.1),par.2=median(par.2)),by=type]
#   myifp.rolled.con[,qgamma(.5,par.1,par.2)]
#   myifp.rolled.con[,pgamma(outcome-roll.date,par.1,par.2)]
  
#   myifp.rolled.con[,prob:=round(median(pgamma(tmp.date-closed, par.1, par.2)),2),by=type]
  myifp.rolled.con[,prob:=pgamma(tmp.date-roll.date, par.1, par.2),by=type]
  myifp.rolled.con[,ifp_id:=ifp]
  myifp.rolled.con[,BS:=BS(prob,tmp.date-roll.date,outcome-roll.date)]
  myifp.rolled.con[,roll.date:=roll.date]
  setkey(myifp.rolled.con,type)  
  setkey(N,type)
  myifp.rolled.con[N[,list(type,n.keep)]]
}




closed

####
f.tmp <- fcasts[ifp_id=="1429-0",][1,]
with(f.tmp, {
  c(pgamma(outcome-as.numeric(as.Date(fcast_date)), par.1, par.2),
    pgamma(ctpt_0, par.1, par.2)-bin_0/100)
})

outcome

myifp.rolled.con[,round(pgamma(outcome-roll.date, par.1, par.2),2)]

tmp.par <- unlist(fits[,list(median(par.1),median(par.2))])

qgamma(.5, tmp.par[1], tmp.par[2])

2*(pgamma(outcome-roll.date, tmp.par[1], tmp.par[2])-1)^2

myifp[,round(pgamma(outcome-roll.date, par.1, par.2),2)]


#### Stop


# Apply fitting function to each forecast
fits <- foreach(i = 1:nrow(myifp)) %dopar% {
  fcasts.ifp[i , FitWrapper(bins  = c(bin_0,  bin_1,  bin_2,  bin_3,  bin_4),
                            ctpts = c(ctpt_0, ctpt_1, ctpt_2, ctpt_3),
                            type  = value_type,
                            date.fcast.made  = fcast_date), 
             by=c("idx")]
}
fits <- data.table(do.call(rbind,fits))
setkey(fits, idx)
setkey(fcasts.ifp, idx)
tmp <- fcasts.ifp[fits]

# Define refitting cutpoints for each forecast
n.bins <- 20
ctpts <- data.table(matrix(rep(as.numeric(NA), nrow(tmp)*n.bins), ncol=n.bins))
bins <- copy(ctpts)
bin.id <- paste0("bin_",1:20)
ctpt.id <- paste0("ctpt_",1:20)
setnames(ctpts, ctpt.id)
setnames(bins, bin.id)
new <- cbind(tmp[,list(idx,par1,par2,fcast_date)],bins,ctpts)

# Obtain probabilities for each bin
new[,c(ctpt.id):=as.list(seq(roll.date-fcast_date, 
                             qgamma(.9999,par1,par2), 
                             length=20)),
    by=idx]
new[,c(bin.id):=as.list(pgamma(seq(roll.date-fcast_date, 
                                   qgamma(.9999,par1,par2), 
                                   length=20),par1,par2)),
    by=idx]
new[,c(bin.id):=as.list(new[,bin.id,with=FALSE] - cbind(0,new[,bin.id[-length(bin.id)],with=FALSE]))]

# Renorm probabilities
new[,denom:=1-pgamma(roll.date-fcast_date,par1,par2)]
new[,c(paste0("bin_",1:20),"denom"):=as.list(.SD/denom),.SDcols=c(paste0("bin_",1:20),"denom"),by=idx]

# Refit
out <- foreach(i=1:nrow(tmp)) %dopar% {
  FitFnInterval(unlist(new[i,bin.id,with=FALSE]), unlist(new[i,ctpt.id,with=FALSE]), "gamma")
}
out <- do.call(rbind,out)
data.table(ifp_id     = unique(fcasts.ifp$ifp_id),
           fcast.date = roll.date,
           dist       = unique(fcasts.ifp$dist),
           par.1      = median(out[,1]),
           par.2      = median(out[,2]))



