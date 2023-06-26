###############################################################################
## This script uses the best available distance model and corresponding density estimates
## in conjunction with weigth data from port samples AND
## estimated weight of ROV fish to produce several biomass estimates:
##	1) Biomass from port weights
##	2) Biomass from estimated ROV weights
##	3) Biomass from estimated ROV weights from fish greater than a prescribed cutoff length
## ROV weights and the proportion greater than a certain length are estimated with a bootstrap
##
## August 20201
## Phil Joy
##
################################################################################
## Calculate Density and Biomass from best model...
###############################################################################
library(boot)
library(RColorBrewer)
library(tidyverse)

setwd("D:/Groundfish Biometrics/Yelloweye/YE Code/")

##Function for estimating weight from a given length using a given linear model of logW on logL... 
w.from.l<-function (LM, length){
	syx<-summary(LM)$sigma
	cf<-exp((syx^2)/2)
	pred.log<-predict(LM,data.frame(logL=log(length)), interval = "c")
	bias.pred.orig <- exp(pred.log) 
	pred.orig<-cf*bias.pred.orig
	pred.orig[1]
}

## Function for getting bootstrap mean and variance of that mean for avg weights

meanfun <- function(x, d) {
  return(mean(x[d]))
}


## Load your list of models that were produced in "DSR_ROV_SEAK_base_distance_modelling.R"
G<-read.csv("NSEO_18_Bootstrap_Model_List_Trunc.csv")

## Pull out the best model...
Best<-G[1,]
Best

## Select the appropriate area for the subdistrict being examined  (NSEO = 442, SSEO = 1056, CSEO = 1661, EYKT = 739)
HA<-442			

### load cleaned LW YE data created with "DSR_ROV_SEAK_LWexam.R" , and pull out the appropriate data
LW<-read.csv("NSEO_YE_LW_rel.csv")
str(LW)
LW<-LW[LW$Groundfish.Management.Area.Code == "NSEO" & LW$Year > 2018,]
nrow(LW)

## create a linear model for estimating ROV weights... 
LW.lm<-lm(logW ~ logL, data=LW) 
summary(LW.lm)

##### State what you want to use as a lower bound cutoff length... 
CUTOFF<-270	#cutoff length

## Load ROV lengths and get rid of the junk measurements... 
ROV<-read_csv("Data/SSEO_2020_species.csv") %>% filter(Species == 145)
ROV.adult <- ROV %>% filter(Stage != "JV")
## get rid of unmeasured fish	#str(ROV.adult[,36])
ROV.L<-as.data.frame(ROV.adult[,c(1:2,4:34)])
ROV.L<-ROV.L[complete.cases(ROV.L[,"Length..mm."]),]
ROV.L$Length<-ROV.L$Length..mm.  
ROV.L$L.prec<-ROV.L$Precision..mm.
ROV.L$rel.Prec<-ROV.L$L.prec/ROV.L$Length
## get rid of imprecise measurements
ROV.goodL<-ROV.L[ROV.L$rel.Prec < 0.5,]	
nrow(ROV.L)
nrow(ROV.goodL)
hist(ROV.goodL$rel.Prec) 

################################################################################################################
## use bootstrap to get mean weights of ROV fish and to get the proportion below your cutoff length

nboot<-1000
O<-data.frame(NA)
pb = txtProgressBar(min = 0, max = length(nboot), initial = 0, style=3) 

for (i in 1:nboot){	#i<-1
	#randomly sample ROV lengths with replacement... 
	ROV.rand<-sample(ROV.goodL$Length, size=length(ROV.goodL$Length),replace=TRUE)
	#get the mean length from that sample and record mean and var
	O$mean.l.all<-mean(ROV.rand)
	O$var.l.all<-var(ROV.rand)

	#get and record the proportion of lengths greater than CUTOFF value... 
	Prop.gt<-length(ROV.rand[ROV.rand>CUTOFF])/length(ROV.rand)
	O$Prop.gt<-Prop.gt
	O$var.Prop<-Prop.gt*(1-Prop.gt)/length(ROV.rand)
	
	#get mean length of fish greater than CUTOFF value for biomass calculations... 
	ROV.gt.co<-ROV.rand[ROV.rand>CUTOFF]
	O$mean.l.gtco<-mean(ROV.gt.co)
	O$var.l.gtco<-var(ROV.gt.co)

	#Weight the ROV fish with the w.from.l function that uses your prescribed linear model... 
	Ws<-sapply(ROV.rand[ROV.rand > 0], function(x)w.from.l(LM=LW.lm, length=x))	#plot(Ws~ROV.rand[ROV.rand>CUTOFF])
	O$mean.W<-mean(Ws)					
		mW.b<-boot(Ws,statistic=meanfun, R=500)				
	O$mean.W.var<-sd(mW.b$t)^2		#var(Ws)

	Ws.gtco<-sapply(ROV.gt.co, function(x)w.from.l(LM=LW.lm, length=x))	#plot(Ws~ROV.rand[ROV.rand>CUTOFF])
	O$mean.W.gtco<-mean(Ws.gtco)					
		mW.b2<-boot(Ws.gtco,statistic=meanfun, R=500)					
	O$mean.W.gtco.var<-sd(mW.b2$t)^2			#var(Ws.gtco)

	if (i == 1) {
	ROV.boot<-O
	} else {
	ROV.boot<-rbind(ROV.boot,O)
	}
setTxtProgressBar(pb,i/nboot)
}

str(ROV.boot)

###############################################################################################
## Take a look at the bootstrap results:

COLS<-brewer.pal(4,'Set1')
par(mfrow = c(3,1),mar=c(4,4,4,1)+0.1,oma=c(2,1,1,1))
plot(density(ROV.goodL$Length), xlab="Length", col="goldenrod", main="ROV lengths",
		xlim=c(min(ROV.goodL$Length),max(ROV.goodL$Length)), 
		ylim=c(0,0.03))
	polygon(density(ROV.goodL$Length),col=adjustcolor(COLS[1],alpha.f=0.25))
	#polygon(density(ROV.L$Length),col=adjustcolor(COLS[2],alpha.f=0.25))
	polygon(density(ROV.boot$mean.l.all), col=adjustcolor(COLS[2],alpha.f=0.25))
	polygon(density(ROV.boot$mean.l.gtco), col=adjustcolor(COLS[3],alpha.f=0.25))
	abline(v=mean(ROV.goodL$Length),col=COLS[1], lwd=2)
	abline(v=mean(ROV.boot$mean.l.all),col=COLS[2])
	abline(v=mean(ROV.boot$mean.l.gtco),col=COLS[3])
	mtext("ROV lengths",side=3, adj=0, padj=1.5, col=COLS[1])
	mtext("mean for all ROV",side=3, adj=0, padj=3, col=COLS[2])
	mtext("mean for ROV gt270",side=3, adj=0, padj=4.5, col=COLS[3])
plot(density(ROV.boot$Prop.gt), xlab="Proportion > 270", col="goldenrod", main="mean prop. ROV YE gt 270mm")
	polygon(density(ROV.boot$Prop.gt), col="goldenrod")
	abline(v=mean(ROV.boot$Prop.gt), col="yellow", lty=1)
hist(ROV.boot$Prop.gt, col="goldenrod", breaks=6)
	abline(v=mean(ROV.boot$Prop.gt), col="yellow", lty=1)

###############################################################
## Retrieve the values needed for estimating biomass: 
##mean proportion of YE less than cutoff length (270)

pgt270<-mean(ROV.boot$Prop.gt)
pgt270.var<-var(ROV.boot$Prop.gt)

##mean weights of individual fish from ROV
ROV.w<-mean(ROV.boot$mean.W)
ROV.w.var<-var(ROV.boot$mean.W)
sqrt(ROV.w.var)
mean(sqrt(ROV.boot$mean.W.var))

ROV.gtco.w<-mean(ROV.boot$mean.W.gtco)
ROV.gtco.w.var<-var(ROV.boot$mean.W.gtco)
sqrt(ROV.gtco.w.var)
mean(sqrt(ROV.boot$mean.W.gtco.var))

######################################################################
## Calcualte 3 biomass estimates.
## These calculations should align with the spreadsheet... 

###### Biomass with ROV (average of ROV weights), all fish

Biomass.ROV.a<-Best$Dhat*HA*pgt270*ROV.w

var.Biomass.ROV.a<-((Best$Dhat.se^2+Best$Dhat^2)*(pgt270.var+pgt270^2)*(ROV.w.var+ROV.w^2)-
	(Best$Dhat^2)*(pgt270^2)*(ROV.w^2))*HA^2
se.Biomass.ROV.a<-sqrt(var.Biomass.ROV.a)

ROV.a.lowCI<-(Biomass.ROV.a-1.645*se.Biomass.ROV.a)/1000
ROV.a.BM<-Biomass.ROV.a/1000						
ROV.a.hiCI<-(Biomass.ROV.a+1.645*se.Biomass.ROV.a)/1000

###### Biomass with ROV (average of ROV weights), GT 270mm

Biomass.ROV.gt<-Best$Dhat*HA*pgt270*ROV.gtco.w

var.Biomass.ROV.gt<-((Best$Dhat.se^2+Best$Dhat^2)*(pgt270.var+pgt270^2)*(ROV.gtco.w.var+ROV.gtco.w^2)-
	(Best$Dhat^2)*(pgt270^2)*(ROV.gtco.w^2))*HA^2
se.Biomass.ROV.gt<-sqrt(var.Biomass.ROV.gt)

ROV.gt.lowCI<-(Biomass.ROV.gt-1.645*se.Biomass.ROV.gt)/1000
ROV.gt.BM<-Biomass.ROV.gt/1000						#good 7/22
ROV.gt.hiCI<-(Biomass.ROV.gt+1.645*se.Biomass.ROV.gt)/1000

###### Biomass from CF avg weight

## first select the port data you want to us
YE.latest<-LW[LW$Year >= 2019 & LW$Groundfish.Management.Area.Code == "NSEO",]
str(YE.latest)
## get the mean weight from the selected data... 
Mean.CF.weight<-mean(YE.latest$Weight.Kilograms )

##calculate the variance of the mean using a bootstrap... 
wboot<-boot(YE.latest$Weight.Kilograms, statistic=meanfun, R=1000)
var.CF.weight<-(sd(wboot$t))^2

## Calculate the biomass: 
Biomass.CFW<-Best$Dhat*HA*pgt270*Mean.CF.weight

var.Biomass.CFW<-((Best$Dhat.se^2+Best$Dhat^2)*(pgt270.var+pgt270^2)*(var.CF.weight+Mean.CF.weight^2)-
	(Best$Dhat^2)*(pgt270^2)*(Mean.CF.weight^2))*HA^2
se.Biomass.CFW<-sqrt(var.Biomass.CFW)

CFW.lowCI<-(Biomass.CFW-1.645*se.Biomass.CFW)/1000
CFW.BM<-Biomass.CFW/1000
CFW.hiCI<-(Biomass.CFW+1.645*se.Biomass.CFW)/1000

######################################################################################