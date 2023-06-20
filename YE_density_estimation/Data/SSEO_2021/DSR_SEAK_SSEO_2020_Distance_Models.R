###################################################################
## Code for Distance Analysis and Models for
## DSR SSEO Rockfish from ROV surveys
## July 2021 - Phil Joy
###################################################################
################################################################################
#EXPLORE THE DATA..
################################
#View Summary of Data
library(tictoc)
library(RColorBrewer)
library(tidyverse)
library(truncnorm)
library(boot)

##Set working directory
setwd("D:/Groundfish Biometrics/Yelloweye/YE Code")

######################################################################
##Load cleaned data from ROV (file DSR_SEAK_SSEO_2020_pjj.R)

SSEO_distance<-read.csv("Data/SSEO_distance_data_GIStran_for_analysis.csv")

str(SSEO_distance)
summary(SSEO_distance$distance)
quantile(SSEO_distance$distance, probs=c(0.9, 0.95, 0.975),na.rm=TRUE)

DAT<-SSEO_distance

TRUNC1<-quantile(SSEO_distance$distance, probs=c(0.95),na.rm=TRUE)

with(DAT, table(Stage,Sample.Label))

################################################################
## Load ROV dada for weight lenth bootstrapping below 

ROV<-read_csv("Data/SSEO_2020_species.csv") %>% filter(Species == 145)
ROV.adult <- ROV %>% filter(Stage != "JV")
## get rid of unmeasured fish	#str(ROV.adult[,36])
ROV.L<-as.data.frame(ROV.adult[,c(1:2,4:34)])
ROV.L<-ROV.L[complete.cases(ROV.L[,"Length..mm."]),]
ROV.L$Length<-ROV.L$Length..mm.  
ROV.L$L.prec<-ROV.L$Precision..mm.
ROV.L$L.cv<-ROV.L$L.prec/ROV.L$Length

##############################################################
## Initial exam of Distance Data

#View Historgram of perpendicular distance from transect line 
hist(DAT$distance, xlab = "Distance (m)", breaks=50)
#maybe a little over-dispersion in the data...
#may consider examining clustering of fish in original data...?  
hist(DAT$Fish.L, xlab = "Fish.L (mm)")
hist(DAT$Depth, xlab = "Depth (m)")

## examine COLINEARITY in covariates... Pearson R2 over 0.15 is no bueno... 
cor(x=DAT$Fish.L, y=DAT$Depth, use="complete.obs", method=c("pearson"))
cor(x=DAT$Depth, y=DAT$Fish.L, use="complete.obs", method=c("pearson"))^2 #
summary(lm(DAT$Fish.L ~ DAT$Depth))
plot(DAT$Fish.L ~ DAT$Depth)
abline(lm(DAT$Fish.L ~ DAT$Depth))

summary(lm(DAT$Fish.L ~ factor(DAT$Stage)))
boxplot(DAT$Fish.L ~ DAT$Stage)

summary(lm(DAT$Depth ~ DAT$Stage, na.rm=T))
boxplot(DAT$Depth ~ DAT$Stage)

## box plot relationships between D and covariates... 
##check out maturity and detection...
ggplot(droplevels(DAT), aes(x = Stage, y = distance)) +
  geom_boxplot() + labs(x = "YE Maturity",
                        y = "distance (m)")
##check out depth and distance
ggplot(DAT, aes(x = Depth, y = distance)) +
  geom_point() + labs(x = "survey depth",
                        y = "distance (m)")
plot(DAT$distance~DAT$Depth)
abline(lm(DAT$distance~DAT$Depth)) #?lm
summary(lm(DAT$distance~DAT$Depth))

Shallow<-DAT[DAT$Depth <150,]
plot(Shallow$distance~Shallow$Depth)
abline(lm(Shallow$distance~Shallow$Depth)) #?lm
summary(lm(Shallow$distance~Shallow$Depth))

##2020 SSEO looks like depth and stage are uncorrelated with distance
## run models below for shits and giggles but don't include in bootstrap...

##check out Fish Length
ggplot(DAT, aes(x = Fish.L, y = distance)) +
  geom_point(alpha = 0.25, size = 1.6) + labs(x = "YE Length (mm)",
                                              y = "Distance (mm)")
plot(DAT$distance~DAT$Fish.L)
abline(lm(DAT$distance~DAT$Fish.L))
summary(lm(DAT$distance~DAT$Fish.L))

##slight tendency for bigger fish to be further away
##not every fish measured so would be complicated to use length as a covariate
##Lenth and stage are strongly related... could model with stage? 

######################################################################
## Check and set Units
######################################################################

units_table()
CU<-convert_units(distance_units="Meter", effort_units="Meter", 
	area_units="Square kilometer")

#######################################################################
## Run Some models...
## hn models, with and without 5% truncation
## hr models, with and without 5% truncation
## uniform models with and without 5% truncation
#######################################################################
## HALF NORMAL MODELS
####################################################################

YE.hn <- ds(DAT, key = "hn", adjustment = NULL,transect="line",
                  convert.units = CU)
summary(YE.hn)
summary(YE.hn$ddf)

plot(YE.hn, nc = 10)

gof_ds(YE.hn)

######################################################################
YE.hn.tr5 <- ds(DAT, key = "hn", adjustment = NULL,transect="line",
                  convert.units = CU, truncation=c(right="5%"))
summary(YE.hn.tr5)
summary(YE.hn.tr5$ddf)

plot(YE.hn.tr5, nc = 10)

gof_ds(YE.hn.tr5)
####################################################################

YE.hn.cos <- ds(DAT, key = "hn", adjustment = "cos", transect="line",
                  convert.units = CU)
summary(YE.hn.cos)
summary(YE.hn.cos$ddf)

plot(YE.hn.cos, nc = 10)

gof_ds(YE.hn.cos)

####################################################################

YE.hn.herm <- ds(DAT, key = "hn", adjustment = "herm", transect="line",
                  convert.units = CU)
summary(YE.hn.herm )
summary(YE.hn.herm $ddf)

plot(YE.hn.herm , nc = 10)

gof_ds(YE.hn.herm )

####################################################################

YE.hn.poly <- ds(DAT, key = "hn", adjustment = "poly", transect="line",
                  convert.units = CU)
summary(YE.hn.poly )
summary(YE.hn.poly $ddf)

plot(YE.hn.poly , nc = 10)

gof_ds(YE.hn.poly )

####################################################################

YE.hn.tr5.cos <- ds(DAT, key = "hn", adjustment = "cos", transect="line",
                  convert.units = CU, truncation=c(right="5%"))
summary(YE.hn.tr5.cos)
summary(YE.hn.tr5.cos$ddf)

plot(YE.hn.tr5.cos, nc = 10)

gof_ds(YE.hn.tr5.cos)

####################################################################

YE.hn.tr5.herm <- ds(DAT, key = "hn", adjustment = "herm", transect="line",
                  convert.units = CU, truncation=c(right="5%"))
summary(YE.hn.tr5.herm )
summary(YE.hn.tr5.herm $ddf)

plot(YE.hn.tr5.herm , nc = 10)

gof_ds(YE.hn.tr5.herm )

####################################################################

YE.hn.tr5.poly <- ds(DAT, key = "hn", adjustment = "poly", transect="line",
                  convert.units = CU, truncation=c(right="5%"))
summary(YE.hn.tr5.poly )
summary(YE.hn.tr5.poly $ddf)

plot(YE.hn.tr5.poly , nc = 10)

gof_ds(YE.hn.tr5.poly )

#####################################################################
YE.hn.Depth <-ds(DAT, key = "hn", transect="line", formula=~Depth,
                  convert.units = CU)
plot(YE.hn.Depth)
gof_ds(YE.hn.Depth)
#####################################################################
YE.hn.Stage <-ds(DAT, key = "hn", transect="line", formula=~Stage,
                  convert.units = CU)
plot(YE.hn.Stage)
gof_ds(YE.hn.Stage)
#####################################################################
YE.hn.Depth.Stage <-ds(DAT, key = "hn", transect="line", formula=~Depth + Stage,
                  convert.units = CU)
plot(YE.hn.Depth.Stage)
gof_ds(YE.hn.Depth.Stage)

#####################################################################
YE.hn.tr5.Depth <-ds(DAT, key = "hn", transect="line", formula=~Depth,
                  convert.units = CU, truncation=c(right="5%"))
plot(YE.hn.tr5.Depth)
gof_ds(YE.hn.tr5.Depth)
#####################################################################
YE.hn.tr5.Stage <-ds(DAT, key = "hn", transect="line", formula=~Stage,
                  convert.units = CU, truncation=c(right="5%"))
plot(YE.hn.tr5.Stage)
gof_ds(YE.hn.tr5.Stage)
#####################################################################
YE.hn.tr5.Depth.Stage <-ds(DAT, key = "hn", transect="line", formula=~Depth + Stage,
                  convert.units = CU, truncation=c(right="5%"))
plot(YE.hn.tr5.Depth.Stage)
gof_ds(YE.hn.tr5.Depth.Stage)
summary(YE.hn.tr5.Depth.Stage)
#####################################################################
## HAZARD RATE MODELS
######################################################################

YE.hr <- ds(DAT, key = "hr", adjustment = NULL, transect="line",
                  convert.units = CU)	#, monotonicity = "strict"

summary(YE.hr)
summary(YE.hr$ddf)

plot(YE.hr, nc = 10)

gof_ds(YE.hr)
#####################################################################
YE.hr.tr5 <- ds(DAT, key = "hr", adjustment = NULL,transect="line",
                  convert.units = CU, truncation=c(right="5%"))
summary(YE.hr.tr5)
summary(YE.hr.tr5$ddf)

plot(YE.hr.tr5, nc = 10)

gof_ds(YE.hr.tr5)
#####################################################################

YE.hr.cos <- ds(DAT, key = "hr", adjustment = "cos", transect="line",
                  convert.units = CU)
summary(YE.hr.cos)
summary(YE.hr.cos$ddf)

plot(YE.hr.cos, nc = 10)

gof_ds(YE.hr.cos)

####################################################################

YE.hr.herm <- ds(DAT, key = "hr", adjustment = "herm", transect="line",
                  convert.units = CU)
summary(YE.hr.herm )
summary(YE.hr.herm $ddf)

plot(YE.hr.herm , nc = 10)

gof_ds(YE.hr.herm )

####################################################################

YE.hr.poly <- ds(DAT, key = "hr", adjustment = "poly", transect="line",
                  convert.units = CU)
summary(YE.hr.poly )
summary(YE.hr.poly $ddf)

plot(YE.hr.poly , nc = 10)

gof_ds(YE.hr.poly )

####################################################################

YE.hr.tr5.cos <- ds(DAT, key = "hr", adjustment = "cos", transect="line",
                  convert.units = CU, truncation=c(right="5%"))
summary(YE.hr.tr5.cos)
summary(YE.hr.tr5.cos$ddf)

plot(YE.hr.tr5.cos, nc = 10)

gof_ds(YE.hr.tr5.cos)

####################################################################

YE.hr.tr5.herm <- ds(DAT, key = "hr", adjustment = "herm", transect="line",
                  convert.units = CU, truncation=c(right="5%"))
summary(YE.hr.tr5.herm )
summary(YE.hr.tr5.herm $ddf)

plot(YE.hr.tr5.herm , nc = 10)

gof_ds(YE.hr.tr5.herm )

####################################################################

YE.hr.tr5.poly <- ds(DAT, key = "hr", adjustment = "poly", transect="line",
                  convert.units = CU, truncation=c(right="5%"))
summary(YE.hr.tr5.poly )
summary(YE.hr.tr5.poly $ddf)

plot(YE.hr.tr5.poly , nc = 10)

gof_ds(YE.hr.tr5.poly )

#####################################################################
YE.hr.Depth <-ds(DAT, key = "hr", transect="line", formula=~Depth,
                  convert.units = CU)
plot(YE.hr.Depth)
gof_ds(YE.hr.Depth)
#####################################################################
YE.hr.Stage <-ds(DAT, key = "hr", transect="line", formula=~Stage,
                  convert.units = CU)
plot(YE.hr.Stage)
gof_ds(YE.hr.Stage)
#####################################################################
YE.hr.Depth.Stage <-ds(DAT, key = "hr", transect="line", formula=~Depth + Stage,
                  convert.units = CU)
plot(YE.hr.Depth.Stage)
gof_ds(YE.hr.Depth.Stage)

#####################################################################
YE.hr.tr5.Depth <-ds(DAT, key = "hr", transect="line", formula=~Depth,
                  convert.units = CU, truncation=c(right="5%"))
plot(YE.hr.tr5.Depth)
gof_ds(YE.hr.tr5.Depth)
#####################################################################
YE.hr.tr5.Stage <-ds(DAT, key = "hr", transect="line", formula=~Stage,
                  convert.units = CU, truncation=c(right="5%"))
plot(YE.hr.tr5.Stage)
gof_ds(YE.hr.tr5.Stage)
#####################################################################
YE.hr.tr5.Depth.Stage <-ds(DAT, key = "hr", transect="line", formula=~Depth + Stage,
                  convert.units = CU, truncation=c(right="5%"))
plot(YE.hr.tr5.Depth.Stage)
gof_ds(YE.hr.tr5.Depth.Stage)
summary(YE.hr.tr5.Depth.Stage)

#####################################################################
## Check Uniform models with adjustments...
#######################################################################

YE.unif.cos <- ds(DAT, key = "unif", adjustment = "cos", transect="line",
                  convert.units = CU)
summary(YE.unif.cos)
summary(YE.unif.cos$ddf)

plot(YE.unif.cos, nc = 10)

gof_ds(YE.unif.cos)

####################################################################

YE.unif.herm <- ds(DAT, key = "unif", adjustment = "herm", transect="line",
                  convert.units = CU)
summary(YE.unif.herm )
summary(YE.unif.herm $ddf)

plot(YE.unif.herm , nc = 10)

gof_ds(YE.unif.herm )

####################################################################

YE.unif.poly <- ds(DAT, key = "unif", adjustment = "poly", transect="line",
                  convert.units = CU)
summary(YE.unif.poly )
summary(YE.unif.poly $ddf)

plot(YE.unif.poly , nc = 10)

gof_ds(YE.unif.poly )

####################################################################

YE.unif.tr5.cos <- ds(DAT, key = "unif", adjustment = "cos", transect="line",
                  convert.units = CU, truncation=c(right="5%"))
summary(YE.unif.tr5.cos)
summary(YE.unif.tr5.cos$ddf)

plot(YE.unif.tr5.cos, nc = 10)

gof_ds(YE.unif.tr5.cos)

####################################################################

YE.unif.tr5.herm <- ds(DAT, key = "unif", adjustment = "herm", transect="line",
                  convert.units = CU, truncation=c(right="5%"))
summary(YE.unif.tr5.herm )
summary(YE.unif.tr5.herm $ddf)

plot(YE.unif.tr5.herm , nc = 10)

gof_ds(YE.unif.tr5.herm )

####################################################################

YE.unif.tr5.poly <- ds(DAT, key = "unif", adjustment = "poly", transect="line",
                  convert.units = CU, truncation=c(right="5%"))
summary(YE.unif.tr5.poly )
summary(YE.unif.tr5.poly $ddf)

plot(YE.unif.tr5.poly , nc = 10)

gof_ds(YE.unif.tr5.poly )

#####################################################################
#####################################################################
## RANK AND COMPARE MODELS
#####################################################################
## non-truncated models:

summarize_ds_models(YE.hn,YE.hn.cos,YE.hn.herm,YE.hn.poly,
	YE.hn.Depth,YE.hn.Stage,YE.hn.Depth.Stage,
	YE.hr,YE.hr.cos,YE.hr.herm,YE.hr.poly,
	YE.hr.Depth,YE.hr.Stage,YE.hr.Depth.Stage,
	YE.unif.cos,YE.unif.herm,YE.unif.poly) 

##above analysis showed much better fit and model stability with truncation
##truncated models

G<-summarize_ds_models(YE.hn.tr5, YE.hn.tr5.cos,YE.hn.tr5.herm, YE.hn.tr5.poly, 
	YE.hn.tr5.Depth,YE.hn.tr5.Stage,YE.hn.tr5.Depth.Stage, 
	YE.hr.tr5, YE.hr.tr5.cos, YE.hr.tr5.herm, YE.hr.tr5.poly, 
	YE.hr.tr5.Depth, YE.hr.tr5.Stage, YE.hr.tr5.Depth.Stage, 
	YE.unif.tr5.cos,YE.unif.tr5.herm,YE.unif.tr5.poly) 
names(G)[names(G) == "Model"]<-"Full.Model.Code"
G$Model<-substr(gsub('[\\{}]','',G$Full.Model.Code),7,50)
names(G)[names(G) == "C-vM p-value"]<-"CvM_pvalue"
names(G)[names(G) == "$\\hat{P_a}$"]<-"P_a.hat"
names(G)[names(G) == "se($\\hat{P_a}$)"]<-"se.P_a.hat"
names(G)[names(G) == "$\\Delta$AIC"]<-"Delta.AIC"
names(G)[names(G) == "Key function"]<-"Key.function"

G$Model<-as.factor(G$Model)

view(G)

##uniform models with different adjustments did well
## hn and hr models did not improve with adjustments (they revert to base hr and hn models)
##scrap them from list
G<-G[-c(4:6,11:13),]

##uniform models with different adjustments did well

##ordered the list of mdels basedon above output, G.  Did by hand, 
## can't figure out how to order the list of models based on AIC.. Argh!

Mods<-list(YE.unif.tr5.poly,YE.unif.tr5.herm,
		YE.hn.tr5,	#YE.hn.tr5.cos,YE.hn.tr5.herm,YE.hn.tr5.poly,
		YE.unif.tr5.cos,
		YE.hn.tr5.Stage, YE.hn.tr5.Depth, 
		YE.hr.tr5, #YE.hr.tr5.cos, YE.hr.tr5.herm, YE.hr.tr5.poly, 
		YE.hn.tr5.Depth.Stage, 
		YE.hr.tr5.Depth, YE.hr.tr5.Stage, YE.hr.tr5.Depth.Stage)

	G$Det.Prob<-sapply(Mods, function(x) summary(x)$dht$individuals$average.p)
	G$Nhat<-sapply(Mods, function(x) summary(x)$dht$individuals$N$Estimate)
	G$Nhat.se<-sapply(Mods, function(x) summary(x)$dht$individuals$N$se)
	G$Nhat.cv<-sapply(Mods, function(x) summary(x)$dht$individuals$N$cv)
	G$Nhat.lcl<-sapply(Mods, function(x) summary(x)$dht$individuals$N$lcl)
	G$Nhat.ucl<-sapply(Mods, function(x) summary(x)$dht$individuals$N$ucl)

	G$Dhat<-sapply(Mods, function(x) summary(x)$dht$individuals$D$Estimate)
	G$Dhat.se<-sapply(Mods, function(x) summary(x)$dht$individuals$D$se)
	G$Dhat.cv<-sapply(Mods, function(x) summary(x)$dht$individuals$D$cv)
	G$Dhat.lcl<-sapply(Mods, function(x) summary(x)$dht$individuals$D$lcl)
	G$Dhat.ucl<-sapply(Mods, function(x) summary(x)$dht$individuals$D$ucl)

### Export this file for bootstrapping...
write.csv(G, file = "Data/Bootstrap_Model_List.csv")

minAIC<-min(as.numeric(sapply(Mods,function(x)AIC(x)[2])))

AIC(Mods[[1]])

##plot detection plots
par(mfrow = c(4,3),mar=c(4,4,4,1)+0.1,oma=c(5,5,5,1))

for (i in 1:length(Mods)) {		#i<-1
	plot(Mods[[i]],cex.axis=1.3, cex.main=50, ylab="",
		xlab=paste("deltaAIC=",round(AIC(Mods[[i]])[2]-minAIC,2)))
	mtext(paste(Mods[[i]]$call$key,Mods[[i]]$call$adjustment, sep=","),side=3,adj=-0.1,cex=1)
	mtext(Mods[[i]]$call$formula,side=3,adj=1,cex=1)
}
mtext("Distance (m)", side=1,outer=TRUE,cex=2)
mtext("Frequency", side=2,outer=TRUE,cex=2)
mtext("Detection Functions for 5% Truncated Models", side=3,outer=TRUE,cex=2)

##plot GOF plots... 
par(mfrow = c(4,3),mar=c(4,4,4,1)+0.1,oma=c(5,5,5,1))

for (i in 1:length(Mods)) {		#i<-1
	gof_ds(Mods[[i]])[1]$dsgof$CvM$p
	mtext(paste(Mods[[i]]$call$key,Mods[[i]]$call$adjustment, sep=","),side=3,adj=-0.1,cex=1)
	mtext(Mods[[i]]$call$formula,side=3,adj=1,cex=1)
#	mtext(Mods[[i]])[1]$dsgof$CvM$p
}
mtext("GOF for 5% Truncated Models", side=3,outer=TRUE,cex=2)

###############################################################################
# PLOT Dhats from different models...

GG<-G[rev(order(G$Delta.AIC, G$Model)),]

ggplot(GG, aes(x=reorder(Model,desc(Model)), y=Dhat))+geom_point()+
	geom_errorbar(aes(ymin=Dhat.lcl, ymax=Dhat.ucl),width=0.2)+
	xlab("Model")+
	ylab("SSEO Yelloweye Density")+
	theme(axis.text.x = element_text(angle=90, hjust=1))+
	scale_y_continuous(name="SSEO Yelloweye Density", breaks=seq(from=0, to=5000, by=500))#, 
		#labels, limits)	

################################################################################
## Calculate Density and Biomass from best model...
###############################################################################
G<-read.csv("Data/Bootstrap_Model_List.csv")

Best<-G[1,]
HA<-1056			#Habitat area for SSEO

#####   WL regression from CF data; LW data from LW exam from script: DSR_SEAK_SSEO_2020_LW_Exam.R
w.from.l<-function (LM, length){
	syx<-summary(LM)$sigma
	cf<-exp((syx^2)/2)
	pred.log<-predict(LM,data.frame(logL=log(length)), interval = "c")
	bias.pred.orig <- exp(pred.log) 
	pred.orig<-cf*bias.pred.orig
	pred.orig[1]
}

### load cleaned LW YE data, and pull out the SSEO data
LW<-read.csv("Data/YE_LW_rel.csv")
str(LW)
LW<-LW[LW$Groundfish.Management.Area.Code == "SSEO",]
nrow(LW)

LW.lm<-lm(logW ~ logL, data=LW) 
summary(LW.lm)

##### ROV lengths bootstrepped to deal with uncertainty in ROV lengths
CUTOFF<-270	#cutoff length

#L<-which( colnames(ROV.L)=="Length" )
#Lp<-which( colnames(ROV.L)=="L.prec" )

#get rid of less precise length measurements...

ROV.L$rel.Prec<-ROV.L$L.prec/ROV.L$Length
## get rid of imprecise measurements
	#ROV.L$badL<-ROV.L$Length-3*ROV.L$L.prec
ROV.goodL<-ROV.L[ROV.L$rel.Prec < 0.4,]	#ROV.good<-ROV.L[ROV.L$L.prec > 0,]
nrow(ROV.L)
nrow(ROV.goodL)
hist(ROV.goodL$L.cv) 

nboot<-1000
O<-data.frame(NA)
pb = txtProgressBar(min = 0, max = length(nboot), initial = 0, style=3) 
for (i in 1:nboot){	#i<-1
	#get lengths by randomly sampling from mean and sd of measured ROV fish...
	#ROV.goodL$randL<-rtruncnorm(n=nrow(ROV.goodL),a=0,b=Inf,mean=ROV.goodL$Length, sd=ROV.goodL$L.prec)		#head(ROV.goodL,10)
	#Next take random sample of those lengths... so now we've accounted for inaccuracy of measurements and randomness of sampling!
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

	#now lets get 2 things : the weight of the average fish and the average or the converted weight of each fish... (very different!)
	#first lets get the weight of the average fish... 
	#O$W.mean.L<-w.from.l(LM=LW.lm , length=mean(ROV.rand))
	#O$W.mean.L.gtco<-w.from.l(LM=LW.lm , length=mean(ROV.gt.co))
	#now weight each ROV fish gt cutoff...
	Ws<-sapply(ROV.rand[ROV.rand > 0], function(x)w.from.l(LM=LW.lm, length=x))	#plot(Ws~ROV.rand[ROV.rand>CUTOFF])
	mean(Ws)
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

## Take a look at the bootstrap results

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

str(ROV.boot)
mean(ROV.boot$Prop.gt)
var(ROV.boot$Prop.gt)
as.numeric(mean(ROV.boot$var.Prop))

mean(ROV.boot$mean.W)
sqrt(mean(ROV.boot$mean.W.var))
mean(ROV.boot$mean.W.gtco)
sqrt(mean(ROV.boot$mean.W.gtco.var))
plot(density(ROV.boot$mean.W.var))

##mean proportion of YE less than cutoff length (270)

pgt270<-mean(ROV.boot$Prop.gt)
pgt270.var<-var(ROV.boot$Prop.gt)

##weights of individual fish from ROV
ROV.w<-mean(ROV.boot$mean.W)
ROV.w.var<-var(ROV.boot$mean.W)
sqrt(ROV.w.var)
mean(sqrt(ROV.boot$mean.W.var))

ROV.gtco.w<-mean(ROV.boot$mean.W.gtco)
ROV.gtco.w.var<-var(ROV.boot$mean.W.gtco)
sqrt(ROV.gtco.w.var)
mean(sqrt(ROV.boot$mean.W.gtco.var))

##average weight from average length ROV fish

mean(ROV.boot$W.mean.L)
sqrt(var(ROV.boot$W.mean.L))

###### Biomass with ROV1 (average of ROV weights), all fish

Biomass.ROV.a<-Best$Dhat*HA*pgt270*ROV.w

var.Biomass.ROV.a<-((Best$Dhat.se^2+Best$Dhat^2)*(pgt270.var+pgt270^2)*(ROV.w.var+ROV.w^2)-
	(Best$Dhat^2)*(pgt270^2)*(ROV.w^2))*HA^2
se.Biomass.ROV.a<-sqrt(var.Biomass.ROV.a)

ROV.a.lowCI<-(Biomass.ROV.a-1.645*se.Biomass.ROV.a)/1000
ROV.a.BM<-Biomass.ROV.a/1000						#not good!
ROV.a.hiCI<-(Biomass.ROV.a+1.645*se.Biomass.ROV.a)/1000

###### Biomass with ROV1 (average of ROV weights), GT 270mm

Biomass.ROV.gt<-Best$Dhat*HA*pgt270*ROV.gtco.w

var.Biomass.ROV.gt<-((Best$Dhat.se^2+Best$Dhat^2)*(pgt270.var+pgt270^2)*(ROV.gtco.w.var+ROV.gtco.w^2)-
	(Best$Dhat^2)*(pgt270^2)*(ROV.gtco.w^2))*HA^2
se.Biomass.ROV.gt<-sqrt(var.Biomass.ROV.gt)

ROV.gt.lowCI<-(Biomass.ROV.gt-1.645*se.Biomass.ROV.gt)/1000
ROV.gt.BM<-Biomass.ROV.gt/1000						#good 7/22
ROV.gt.hiCI<-(Biomass.ROV.gt+1.645*se.Biomass.ROV.gt)/1000

###### Biomass from CF avg weight

YE.LW<-read.csv("Data/YE_LW_rel.csv")
	##most recent data for select area (2021 = SSEO)
YE.latest<-YE.LW[YE.LW$Year >= 2019 & YE.LW$Groundfish.Management.Area.Code == "SSEO",]
str(YE.latest)
Mean.CF.weight<-mean(YE.latest$Weight.Kilograms )

##var.CF.weight<-var(YE.latest$Weight.Kilograms ) ##use boot to get var of estimate of the mean...
meanfun <- function(x, d) {
  return(mean(x[d]))
}
 
wboot<-boot(YE.latest$Weight.Kilograms, statistic=meanfun, R=1000)
var.CF.weight<-(sd(wboot$t))^2

Biomass.CFW<-Best$Dhat*HA*pgt270*Mean.CF.weight

var.Biomass.CFW<-((Best$Dhat.se^2+Best$Dhat^2)*(pgt270.var+pgt270^2)*(var.CF.weight+Mean.CF.weight^2)-
	(Best$Dhat^2)*(pgt270^2)*(Mean.CF.weight^2))*HA^2
se.Biomass.CFW<-sqrt(var.Biomass.CFW)

CFW.lowCI<-(Biomass.CFW-1.645*se.Biomass.CFW)/1000
CFW.BM<-Biomass.CFW/1000
CFW.hiCI<-(Biomass.CFW+1.645*se.Biomass.CFW)/1000

#######################################################################
##Go to bootstrap proceedure after selecting best models
## Bootstrap will allow model averaging and bootstrapping mean weights, lengths, etc...
## DSR_SEAK_SSEO_2020_Distance_Bootstrap_wBiomass.R
########################################################################################



















