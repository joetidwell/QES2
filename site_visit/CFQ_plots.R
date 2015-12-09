library(ggplot2)
library(foreach)
library(doMC)
registerDoMC()

source("~/ACE/global_vars.R")
source(file.path(kRootPath, "util", "load_data.R"), chdir=T)
source(file.path(kRootPath, "fitting", "interval_fitting_funcs.R"), chdir=T)
source(file.path(kRootPath, "forecast", "method_consensus_dist.R"), chdir=T)

theme_set(theme_classic())
options(stringsAsFactors = FALSE)
path.mydata <- "~/git/QES2/data"

# Load data from grofo
load(file.path(path.mydata,"tmp.RData"))

p <- ggplot(data=ThetaM, aes(x=type, y=S.R.r)) +
  geom_boxplot(aes(color=dist), outlier.size=0, notch=TRUE) +
  facet_wrap(~ifp_id)
p


tmp <- MDBS[,mean(BS,na.rm=TRUE),by=ifp_id]
setkey(tmp,V1)
MDBS[,ifp_id:=factor(ifp_id,levels=tmp$ifp_id)]
ggplot(data=MDBS[is.finite(BS)], aes(x=ifp_id, y=BS, color=type)) +
  geom_hline(yintercept=.5) +
  geom_point(size=5.4, color="black") +
  geom_point(size=5) +
  ylim(0,2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x="IFP", y="Mean Daily Brier Score") +
  scale_color_discrete(name="Condition")


tmp <- MDBS[,mean(BS,na.rm=TRUE),by=ifp_id]
setkey(tmp,V1)
MDBS[,ifp_id:=factor(ifp_id,levels=tmp$ifp_id)]
ggplot(data=MDBS[is.finite(BS),list(BS=mean(BS)),by=c("ifp_id","dist")], aes(x=ifp_id, y=BS, color=dist)) +
  geom_hline(yintercept=.5) +
  geom_point(size=5.4, color="black") +
  geom_point(size=5) +
  ylim(0,2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x="IFP", y="Mean Daily Brier Score") +
  scale_color_discrete(name="Distribution")

cond1h <- c(NaN,NaN,.159,.714,.202,.227,.451,.168,.323,.382,.188,.137,.221,.099,.886,.160,.546,.164,.182,.688,.345,.192,.236,.323,.277,.289,.427,.127,.210,.190,.768,.658,.502,.768,.924,.462)

tmp <- MDBS[,mean(BS,na.rm=TRUE),by=ifp_id]
setkey(tmp,V1)

MDBS[,ifp_id:=factor(ifp_id,levels=tmp$ifp_id)]
blerg <- MDBS[is.finite(BS),list(BS=mean(BS)),by=c("ifp_id","dist")]
tmp[,V1:=cond1h]
tmp[,BS:=V1]
tmp[,V1:=NULL]
tmp$dist <- "cond1h"
blerg <- rbind(blerg,tmp)
blerg[,dist:=factor(dist,levels=c("beta","gamma","normal",""))]
blerg<- blerg[!is.nan(BS)]
blerg[dist!="",dist:="Consensus"]
blerg[is.na(dist),dist:="ULinOP"]
blerg[,dist:=ordered(dist,levels=c("UlinOP","Consensus"))]

ggplot(data=blerg, aes(x=ifp_id, y=BS, color=dist)) +
  geom_hline(yintercept=.5) +
  geom_point(size=5.4, color="black") +
  geom_point(data=blerg[dist==""], aes(x=ifp_id, y=BS), color="steelblue",size=5) +
  geom_point(data=blerg[dist!=""], aes(x=ifp_id, y=BS), color="firebrick",size=5) +
  ylim(0,2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x="Forecasting Question", y="Mean Daily Brier Score") 
  # scale_color_manual(name="Blah", values=c("firebrick","firebrick","firebrick","white"))

ggplot(data=blerg, aes(x=ifp_id, y=BS, color=dist)) +
  geom_hline(yintercept=.5) +
  geom_point(size=5.4, color="black") +
  # geom_point(data=blerg[dist==""], aes(x=ifp_id, y=BS), color="steelblue",size=5) +
  # geom_point(data=blerg[dist!=""], aes(x=ifp_id, y=BS), color="firebrick",size=5) +
  geom_point(size=5) +
  ylim(0,2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x="Forecasting Question", y="Mean Daily Brier Score") +
  scale_color_manual(name="Source",values=c("firebrick","steelblue"))



tmp <- MDBS.ex[,mean(BS,na.rm=TRUE),by=ifp_id]
setkey(tmp,V1)
MDBS.ex[,ifp_id:=factor(ifp_id,levels=tmp$ifp_id)]
ggplot(data=MDBS.ex[is.finite(BS),list(BS=mean(BS)),by=c("ifp_id","dist")], aes(x=ifp_id, y=BS, color=dist)) +
  geom_hline(yintercept=.5) +
  geom_point(size=5.4, color="black") +
  geom_point(size=5) +
  ylim(0,2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x="IFP", y="Mean Daily Brier Score") +
  scale_color_discrete(name="Distribution")

blah
library(reshape2)
blah <- melt(blah)
tmp <- blah[,mean(value, na.rm=TRUE),by=ifp_id]
blerg

setkey(tmp,V1)
blah[,ifp_id:=factor(ifp_id,levels=tmp$ifp_id)]
blah[,col:=as.numeric(variable)]
blah[col==2,col:=1.5]
blah[col==3,col:=3]
blah[col==4,col:=5]

blerg2 <- blerg[dist=="ULinOP"]
blerg2[,value:=BS]
blerg2[,variable:=dist]

ggplot(data=blah[is.finite(value) & col!=1.5,], aes(x=ifp_id, y=value, color=variable)) +
  geom_hline(yintercept=.5) +
  # geom_point(size=5.4, color="black") +
  geom_point(data=blerg2,size=5,color="black") +
  geom_point(data=blerg2,size=4.6,color="white") +
  geom_point(size=5) +
  # geom_line(size=1.25) +
  ylim(0,2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x="IFP", y="Mean Daily Brier Score") +
  scale_color_manual(name="Variance\nInflation\nFactor", values=c("mistyrose","lightpink","firebrick"))
# scale_color_manual(name="Variance\nInflation\nFactor", values=c("#deebf7","#9ecae1","#3182bd"))


blah


MDBS[ifp_id=="1432-0"]

tmp <- ThetaM[ifp_id=="1442-0" & type=="random",]
tmp.ln <- data.table(value=c(as.numeric(as.Date("2014-10-09")),
                             as.numeric(as.Date("2015-06-01"))),
                     threshold=c("Outcome","GJP Cutpoint"))
pdata <- tmp[,list(median.date=roll.date+qgamma(.5,par.1,par.2),
          roll.date=roll.date,
          type=type),]
pdata[,threshold:="Forecast"]
pdata[,c("ymin","ymax"):=tmp[,list(roll.date+qgamma(.01,par.1,par.2),
                      roll.date+qgamma(.99,par.1,par.2))]]
pdata[,linetype:=c("dashed")]
pdata[,fill:=c("grey")]

p <- ggplot(pdata, aes(x=roll.date,y=median.date)) +
  geom_ribbon(aes(ymin=ymin, ymax=ymax,linetype=linetype, fill=fill), alpha=.2, size=.3, show.guide=FALSE) +  geom_line(size=1) +
  geom_hline(data=tmp.ln, aes(yintercept=value, color=threshold), size=1) +
  scale_color_manual(name="", values=c("firebrick","steelblue")) +
  scale_fill_manual(values="black", guide=FALSE) +
  scale_linetype_manual(values="dashed", guide=FALSE) +
  labs(x="Forecast Date", y="Predicted Date") +
  coord_cartesian(xlim=c(as.Date("2014-08-20"),as.Date("2014-08-25")))
p

tmp <- ThetaM[ifp_id=="1415-0" & type=="random",]
tmp.ln <- data.table(value=c(as.numeric(as.Date("2014-08-26")),
                             as.numeric(as.Date("2014-10-01"))),
                     threshold=c("Outcome","GJP Cutpoint"))
pdata <- tmp[,list(median.date=roll.date+qgamma(.5,par.1,par.2),
          roll.date=roll.date,
          type=type),]
pdata[,threshold:="Forecast"]
pdata[,c("ymin","ymax"):=tmp[,list(roll.date+qgamma(.01,par.1,par.2),
                      roll.date+qgamma(.99,par.1,par.2))]]
pdata[,linetype:=c("dashed")]
pdata[,fill:=c("grey")]

ggplot(pdata, aes(x=roll.date,y=median.date, color=threshold)) +
  geom_ribbon(aes(ymin=ymin, ymax=ymax,linetype=linetype,fill=fill), alpha=.2, size=.3, show.guide=FALSE) +  geom_line(size=1) +
  geom_hline(data=tmp.ln, aes(yintercept=value, color=threshold), size=1) +
  scale_color_manual(name="", values=c("black","firebrick","steelblue")) +
  scale_fill_manual(values="black", guide=FALSE) +
  scale_linetype_manual(values="dashed", guide=FALSE) +
  labs(x="Forecast Date", y="Predicted Date") +
  coord_cartesian(xlim=c(as.Date("2014-08-20"),as.Date("2014-08-25"))))

p


  geom_hline(yintercept=as.numeric(as.Date("2014-10-09"))) +
  geom_hline(yintercept=as.numeric(as.Date("2015-06-01")))


load("roll.Rdata")
tmp <- MDBS[,mean(BS,na.rm=TRUE),by=ifp_id]
setkey(tmp,V1)
MDBS[,ifp_id:=factor(ifp_id,levels=tmp$ifp_id)]

ggplot(data=MDBS[is.finite(BS)], aes(x=ifp_id, y=BS, color=type)) +
  geom_hline(yintercept=.5) +
  geom_point(size=5.4, color="black") +
  geom_point(size=5) +
  ylim(0,2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x="IFP", y="Mean Daily Brier Score") +
  scale_color_discrete(name="Condition")

pdata <- rbind(MDBS,
              data.table(type="Inkling Mean",
                         ifp_id=unique(MDBS$ifp_id),
                         dist=NA,
                         BS=c(.025, .067, .195, .043, .491, .564, .092, 1.246, .411, .287, .156),
                         N=NA))

ggplot(data=pdata[is.finite(BS)], aes(x=ifp_id, y=BS, color=type, group=type, alpha=type, shape=type)) +
  geom_hline(yintercept=.5) +
  geom_point(size=5.4, color="black") +
  geom_point(size=5) +
  # geom_line(size=1.25) +
  ylim(0,2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x="IFP", y="Mean Daily Brier Score") +
  scale_color_manual(name="Method", values=c("black","#d95f02","#1b9e77")) +
  scale_alpha_manual(name="Method",values=c(.5,1,1)) +
  scale_shape_manual(name="Method",values=c(1,16,16))

pdata[type=="rolling", type:="Rolling"]
pdata[type=="rollingCtl", type:="Rolling Control"]


x <- 1:100
ymin <- rep(1,100)
ymax <- rep(2,100)
y <- rep(1.5,100)
fill=rep("Forecast w/ 95% CI",100)
col = rep(c("IFP Cutpoint","Outcome"),times=c(20,80))
blah <- data.table(x,y,ymin,ymax,fill,col)
ggplot(blah, aes(x=x,y=y,ymin=ymin,ymax=ymax, color=col)) +
  # geom_ribbon() +
  geom_line(size=1) +
    # scale_fill_manual(values="grey") +
  # scale_linetype_manual(values=c("solid","dashed")) 
  scale_color_manual(values=c("firebrick","steelblue"))


ggplot(pdata, aes(x=roll.date,y=median.date)) +
  geom_ribbon(aes(ymin=ymin, ymax=ymax,linetype=linetype, fill=fill), alpha=.2, size=.3, show.guide=FALSE) +  geom_line(size=1) +
  geom_hline(data=tmp.ln, aes(yintercept=value, color=threshold), size=1) +
  scale_color_manual(name="", values=c("firebrick","steelblue")) +
  scale_fill_manual(values="black", guide=FALSE) +
  scale_linetype_manual(values="dashed", guide=FALSE) +
  labs(x="Forecast Date", y="Predicted Date") +
  coord_cartesian(xlim=c(as.Date("2014-08-20"),as.Date("2014-08-25")))

a <- c(.1,.2,.4,.2,.1)

a <- a[-1]
a[length(a)] <- a[length(a)]/2
a <- c(a,a[length(a)])
a <- a/sum(a) 

1/2^c(1:5)