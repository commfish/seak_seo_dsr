###################################################################
## Code for Distance Analysis and Models for
## DSR Rockfish from ROV surveys
## This script provides code for basic exploration of the process ROV distance files
## and running a base set of models for consideration
## August 2021 - Phil Joy
###################################################################
################################################################################
#EXPLORE THE DATA..
################################
#View Summary of Data
{library(tictoc)
library(RColorBrewer)
library(tidyverse)
library(truncnorm)
library(boot)
library(Distance)}


YEAR<-2023
Subd<-"EYKT"
#surveyed area (NSEO = 442, SSEO = 1056, CSEO = 1661, EYKT = 739)
surveyed_area<-739

##Set working directory
#setwd("D:/Groundfish Biometrics/Yelloweye/YE Code")

######################################################################
##Load cleaned data from ROV created with "DSR_ROV_SEAK_ROVprocessing.R"
# DAT<-distance
DAT<-read.csv(paste0("YE_density_estimation/Data/",Subd,"_",YEAR,"/",Subd,"_",YEAR,"_distance_data_GIStran_for_analysis.csv"))

DAT <- DAT[-c(1, 2, 3, 142), ]


str(DAT)
summary(DAT$distance)
quantile(DAT$distance, probs=c(0.9, 0.95, 0.975),na.rm=TRUE)

TRUNC1<-quantile(DAT$distance, probs=c(0.95),na.rm=TRUE)

with(DAT, table(Stage,Sample.Label))

##############################################################
## Examine the distance data for distributions and relationships to consider in modelling... 

#View Historgram of perpendicular distance from transect line 
hist(DAT$distance, xlab = "Distance (m)", breaks=50)
#maybe a little over-dispersion in the data...
#may consider examining clustering of fish in original data...?  
hist(DAT$Fish.L.mm, xlab = "Fish.L (mm)")
hist(DAT$Depth, xlab = "Depth (m)")

## examine COLINEARITY in covariates... Pearson R2 over 0.15 is no bueno... 
cor(x=DAT$Fish.L.mm, y=DAT$Depth, use="complete.obs", method=c("pearson")) #0.126188
cor(x=DAT$Depth, y=DAT$Fish.L.mm, use="complete.obs", method=c("pearson"))^2 #0.01592341
summary(lm(DAT$Fish.L.mm ~ DAT$Depth))
plot(DAT$Fish.L.mm ~ DAT$Depth)
abline(lm(DAT$Fish.L.mm ~ DAT$Depth))

# So when you're running models you want to avoid using two variables in the same 
# model if they are highly correlated. I was thinking about using fish length as 
# a covariate but you can't actually do that because you don't have a fish length 
# for every fish. When you get down to fitting the distance models you'll see that 
# the only covariates I used were life stage (adult/subadult) and depth. 

summary(lm(DAT$Fish.L.mm ~ factor(DAT$Stage)))
boxplot(DAT$Fish.L.mm ~ DAT$Stage)
#2023: R^2 is .1776 - fish length has a moderate effect by stage

summary(lm(DAT$Depth ~ DAT$Stage, na.rm=T))
boxplot(DAT$Depth ~ DAT$Stage)
#2023: stage doesn't seems like its a strong predictor of depth

## box plot relationships between D and covariates... 
##check out maturity and detection...
ggplot(droplevels(DAT), aes(x = Stage, y = distance)) +
  geom_boxplot() + labs(x = "YE Maturity",
                        y = "distance (m)")

summary(lm(DAT$distance~DAT$Stage))
#2023: maturity is a significant factor but does not explain the variability in distance

##check out depth and distance
ggplot(DAT, aes(x = Depth, y = distance)) +
  geom_point() + labs(x = "survey depth",
                        y = "distance (m)")
plot(DAT$distance~DAT$Depth)
abline(lm(DAT$distance~DAT$Depth)) #?lm
summary(lm(DAT$distance~DAT$Depth))
#2023: the effect of depth on distance is not significant

Shallow<-DAT[DAT$Depth <150,]
plot(Shallow$distance~Shallow$Depth)
abline(lm(Shallow$distance~Shallow$Depth)) #?lm
summary(lm(Shallow$distance~Shallow$Depth))
#in 2023 the deepest transect was at 137m

##2020 SSEO looks like depth and stage are uncorrelated with distance

##check out Fish Length
ggplot(DAT, aes(x = Fish.L.mm, y = distance)) +
  geom_point(alpha = 0.25, size = 1.6) + labs(x = "YE Length (mm)",
                                              y = "Distance (mm)")
plot(DAT$distance~DAT$Fish.L.mm)
abline(lm(DAT$distance~DAT$Fish.L.mm))
summary(lm(DAT$distance~DAT$Fish.L.mm))

##in previous years, slight tendency for bigger fish to be further away
##not every fish measured so would be complicated to use length as a covariate
##Length and stage are strongly related... could model with stage? 

######################################################################
## Check and set Conversion Units
######################################################################
## Check Effort, Area and distance units before proceeding
## SSEO 2021 Effort was in meters
## NSEO 2018 Effort already converted to km, but distance in meters... 
head(DAT)

units_table()
CU<-convert_units(distance_units="Meter", effort_units="Meter", 
	area_units="Square kilometer")
str(CU)
#######################################################################
## Run Some models...
## hn models, with and without 5% truncation
## hr models, with and without 5% truncation
## uniform models with and without 5% truncation
#######################################################################
## HALF NORMAL MODELS
####################################################################

YE.hn <- ds(DAT, key = "hn", adjustment = NULL, transect="line",
                  convert_units = CU)
summary(YE.hn)
summary(YE.hn$ddf)

plot(YE.hn, nc = 10)

gof_ds(YE.hn)

######################################################################
YE.hn.tr5 <- ds(DAT, key = "hn", adjustment = NULL,transect="line",
                convert_units = CU, truncation=c(right="5%"))
summary(YE.hn.tr5)
summary(YE.hn.tr5$ddf)

plot(YE.hn.tr5, nc = 10)

gof_ds(YE.hn.tr5)
####################################################################

YE.hn.cos <- ds(DAT, key = "hn", adjustment = "cos", transect="line",
                convert_units = CU)
summary(YE.hn.cos)
summary(YE.hn.cos$ddf)

plot(YE.hn.cos, nc = 10)

gof_ds(YE.hn.cos)

####################################################################

YE.hn.herm <- ds(DAT, key = "hn", adjustment = "herm", transect="line",
                 convert_units = CU)
summary(YE.hn.herm )
summary(YE.hn.herm $ddf)

plot(YE.hn.herm , nc = 10)

gof_ds(YE.hn.herm )

####################################################################

YE.hn.poly <- ds(DAT, key = "hn", adjustment = "poly", transect="line",
                 convert_units = CU)
summary(YE.hn.poly )
summary(YE.hn.poly $ddf)

plot(YE.hn.poly , nc = 10)

gof_ds(YE.hn.poly )

####################################################################

YE.hn.tr5.cos <- ds(DAT, key = "hn", adjustment = "cos", transect="line",
                    convert_units = CU, truncation=c(right="5%"))
summary(YE.hn.tr5.cos)
summary(YE.hn.tr5.cos$ddf)

plot(YE.hn.tr5.cos, nc = 10)

gof_ds(YE.hn.tr5.cos)

####################################################################

YE.hn.tr5.herm <- ds(DAT, key = "hn", adjustment = "herm", transect="line",
                     convert_units = CU, truncation=c(right="5%"))
summary(YE.hn.tr5.herm )
summary(YE.hn.tr5.herm $ddf)

plot(YE.hn.tr5.herm , nc = 10)

gof_ds(YE.hn.tr5.herm )

####################################################################

YE.hn.tr5.poly <- ds(DAT, key = "hn", adjustment = "poly", transect="line",
                     convert_units = CU, truncation=c(right="5%"))
summary(YE.hn.tr5.poly )
summary(YE.hn.tr5.poly $ddf)

plot(YE.hn.tr5.poly , nc = 10)

gof_ds(YE.hn.tr5.poly )

#####################################################################
YE.hn.Depth <-ds(DAT, key = "hn", transect="line", formula=~Depth,
                 convert_units = CU)
plot(YE.hn.Depth)
gof_ds(YE.hn.Depth)
#####################################################################
YE.hn.Stage <-ds(DAT, key = "hn", transect="line", formula=~Stage,
                 convert_units = CU)
plot(YE.hn.Stage)
gof_ds(YE.hn.Stage)
#####################################################################
YE.hn.Depth.Stage <-ds(DAT, key = "hn", transect="line", formula=~Depth + Stage,
                       convert_units = CU)
plot(YE.hn.Depth.Stage)
gof_ds(YE.hn.Depth.Stage)

#####################################################################
YE.hn.tr5.Depth <-ds(DAT, key = "hn", transect="line", formula=~Depth,
                     convert_units = CU, truncation=c(right="5%"))
plot(YE.hn.tr5.Depth)
gof_ds(YE.hn.tr5.Depth)
#####################################################################
YE.hn.tr5.Stage <-ds(DAT, key = "hn", transect="line", formula=~Stage,
                     convert_units = CU, truncation=c(right="5%"))
plot(YE.hn.tr5.Stage)
gof_ds(YE.hn.tr5.Stage)
#####################################################################
YE.hn.tr5.Depth.Stage <-ds(DAT, key = "hn", transect="line", formula=~Depth + Stage,
                           convert_units = CU, truncation=c(right="5%"))
plot(YE.hn.tr5.Depth.Stage)
gof_ds(YE.hn.tr5.Depth.Stage)
summary(YE.hn.tr5.Depth.Stage)
#####################################################################
## HAZARD RATE MODELS
######################################################################

YE.hr <- ds(DAT, key = "hr", adjustment = NULL, transect="line",
            convert_units = CU)	#, monotonicity = "strict"

summary(YE.hr)
summary(YE.hr$ddf)

plot(YE.hr, nc = 10)

gof_ds(YE.hr)
#####################################################################
YE.hr.tr5 <- ds(DAT, key = "hr", adjustment = NULL,transect="line",
                convert_units = CU, truncation=c(right="5%"))
summary(YE.hr.tr5)
summary(YE.hr.tr5$ddf)

plot(YE.hr.tr5, nc = 10)

gof_ds(YE.hr.tr5)
#####################################################################

YE.hr.cos <- ds(DAT, key = "hr", adjustment = "cos", transect="line",
                  convert_units = CU)
summary(YE.hr.cos)
summary(YE.hr.cos$ddf)

plot(YE.hr.cos, nc = 10)

gof_ds(YE.hr.cos)

####################################################################

YE.hr.herm <- ds(DAT, key = "hr", adjustment = "herm", transect="line",
                  convert_units = CU)
summary(YE.hr.herm )
summary(YE.hr.herm $ddf)

plot(YE.hr.herm , nc = 10)

gof_ds(YE.hr.herm )

####################################################################

YE.hr.poly <- ds(DAT, key = "hr", adjustment = "poly", transect="line",
                  convert_units = CU)
summary(YE.hr.poly )
summary(YE.hr.poly $ddf)

plot(YE.hr.poly , nc = 10)

gof_ds(YE.hr.poly )

####################################################################

YE.hr.tr5.cos <- ds(DAT, key = "hr", adjustment = "cos", transect="line",
                  convert_units = CU, truncation=c(right="5%"))
summary(YE.hr.tr5.cos)
summary(YE.hr.tr5.cos$ddf)

plot(YE.hr.tr5.cos, nc = 10)

gof_ds(YE.hr.tr5.cos)

####################################################################

YE.hr.tr5.herm <- ds(DAT, key = "hr", adjustment = "herm", transect="line",
                  convert_units = CU, truncation=c(right="5%"))
summary(YE.hr.tr5.herm )
summary(YE.hr.tr5.herm $ddf)

plot(YE.hr.tr5.herm , nc = 10)

gof_ds(YE.hr.tr5.herm )

####################################################################

YE.hr.tr5.poly <- ds(DAT, key = "hr", adjustment = "poly", transect="line",
                  convert_units = CU, truncation=c(right="5%"))
summary(YE.hr.tr5.poly )
summary(YE.hr.tr5.poly $ddf)

plot(YE.hr.tr5.poly , nc = 10)

gof_ds(YE.hr.tr5.poly )

#####################################################################
YE.hr.Depth <-ds(DAT, key = "hr", transect="line", formula=~Depth,
                  convert_units = CU)
plot(YE.hr.Depth)
gof_ds(YE.hr.Depth)
#####################################################################
YE.hr.Stage <-ds(DAT, key = "hr", transect="line", formula=~Stage,
                  convert_units = CU)
plot(YE.hr.Stage)
gof_ds(YE.hr.Stage)
#####################################################################
YE.hr.Depth.Stage <-ds(DAT, key = "hr", transect="line", formula=~Depth + Stage,
                  convert_units = CU)
plot(YE.hr.Depth.Stage)
gof_ds(YE.hr.Depth.Stage)

#####################################################################
YE.hr.tr5.Depth <-ds(DAT, key = "hr", transect="line", formula=~Depth,
                  convert_units = CU, truncation=c(right="5%"))
plot(YE.hr.tr5.Depth)
gof_ds(YE.hr.tr5.Depth)
#####################################################################
YE.hr.tr5.Stage <-ds(DAT, key = "hr", transect="line", formula=~Stage,
                  convert_units = CU, truncation=c(right="5%"))
plot(YE.hr.tr5.Stage)
gof_ds(YE.hr.tr5.Stage)
#####################################################################
YE.hr.tr5.Depth.Stage <-ds(DAT, key = "hr", transect="line", formula=~Depth + Stage,
                  convert_units = CU, truncation=c(right="5%"))
plot(YE.hr.tr5.Depth.Stage)
gof_ds(YE.hr.tr5.Depth.Stage)
summary(YE.hr.tr5.Depth.Stage)

#####################################################################
## Check Uniform models with adjustments...
#######################################################################

YE.unif.cos <- ds(DAT, key = "unif", adjustment = "cos", transect="line",
                  convert_units = CU)
summary(YE.unif.cos)
summary(YE.unif.cos$ddf)

plot(YE.unif.cos, nc = 10)

gof_ds(YE.unif.cos)

####################################################################

YE.unif.herm <- ds(DAT, key = "unif", adjustment = "herm", transect="line",
                  convert_units = CU)
summary(YE.unif.herm )
summary(YE.unif.herm $ddf)

plot(YE.unif.herm , nc = 10)

gof_ds(YE.unif.herm )

####################################################################

YE.unif.poly <- ds(DAT, key = "unif", adjustment = "poly", transect="line",
                  convert_units = CU)
summary(YE.unif.poly )
summary(YE.unif.poly $ddf)

plot(YE.unif.poly , nc = 10)

gof_ds(YE.unif.poly )

####################################################################

YE.unif.tr5.cos <- ds(DAT, key = "unif", adjustment = "cos", transect="line",
                  convert_units = CU, truncation=c(right="5%"))
summary(YE.unif.tr5.cos)
summary(YE.unif.tr5.cos$ddf)

plot(YE.unif.tr5.cos, nc = 10)

gof_ds(YE.unif.tr5.cos)

####################################################################

YE.unif.tr5.herm <- ds(DAT, key = "unif", adjustment = "herm", transect="line",
                  convert_units = CU, truncation=c(right="5%"))
summary(YE.unif.tr5.herm )
summary(YE.unif.tr5.herm $ddf)

plot(YE.unif.tr5.herm , nc = 10)

gof_ds(YE.unif.tr5.herm )

####################################################################

YE.unif.tr5.poly <- ds(DAT, key = "unif", adjustment = "poly", transect="line",
                  convert_units = CU, truncation=c(right="5%"))
summary(YE.unif.tr5.poly )
summary(YE.unif.tr5.poly $ddf)

plot(YE.unif.tr5.poly , nc = 10)

gof_ds(YE.unif.tr5.poly )

#####################################################################
#####################################################################
## RANK AND COMPARE MODELS
#####################################################################
## non-truncated models:

NT<-summarize_ds_models(YE.hn,YE.hn.cos,YE.hn.herm,YE.hn.poly,
	YE.hn.Depth,YE.hn.Stage,	YE.hn.Depth.Stage,
	YE.hr,YE.hr.cos,YE.hr.herm,YE.hr.poly,
	YE.hr.Depth,YE.hr.Stage,	YE.hr.Depth.Stage,
	YE.unif.cos,YE.unif.herm,YE.unif.poly) 
names(NT)[names(NT) == "Model"]<-"Full.Model.Code"
NT$Model<-substr(gsub('[\\{}]','',NT$Full.Model.Code),7,50)
names(NT)[names(NT) == "C-vM p-value"]<-"CvM_pvalue"
names(NT)[names(NT) == "$\\hat{P_a}$"]<-"P_a.hat"
names(NT)[names(NT) == "se($\\hat{P_a}$)"]<-"se.P_a.hat"
names(NT)[names(NT) == "$\\Delta$AIC"]<-"Delta.AIC"
names(NT)[names(NT) == "Key function"]<-"Key.function"

NT$Model<-as.factor(NT$Model)
summary(YE.hr)

##above analysis showed much better fit and model stability with truncation
##truncated models

TR<-summarize_ds_models(YE.hn.tr5, YE.hn.tr5.cos,YE.hn.tr5.herm, YE.hn.tr5.poly, 
	YE.hn.tr5.Depth,
	YE.hn.tr5.Stage,
	YE.hn.tr5.Depth.Stage, 
	YE.hr.tr5, YE.hr.tr5.cos, YE.hr.tr5.herm, YE.hr.tr5.poly, 
	YE.hr.tr5.Depth, YE.hr.tr5.Stage, YE.hr.tr5.Depth.Stage) 
	#YE.unif.tr5.cos,
	#YE.unif.tr5.herm,
	#YE.unif.tr5.poly) 
names(TR)[names(TR) == "Model"]<-"Full.Model.Code"
TR$Model<-substr(gsub('[\\{}]','',TR$Full.Model.Code),7,50)
names(TR)[names(TR) == "C-vM p-value"]<-"CvM_pvalue"
names(TR)[names(TR) == "$\\hat{P_a}$"]<-"P_a.hat"
names(TR)[names(TR) == "se($\\hat{P_a}$)"]<-"se.P_a.hat"
names(TR)[names(TR) == "$\\Delta$AIC"]<-"Delta.AIC"
names(TR)[names(TR) == "Key function"]<-"Key.function"

TR$Model<-as.factor(TR$Model)
summary(YE.hr.tr5.Stage)
view(TR)

##uniform models with different adjustments did well and needed to removed to summarize the models
## hn and hr models did not improve with adjustments (they revert to base hr and hn models)
##scrap redundant models from the list from list
#NT.red<-data.frame()
for (i in 1:nrow(NT)) { #i<-1
  m<-NT[i,]
  if (i == 1){
    NT.red<-m
    t<-2
  } else {
    if (m$Delta.AIC == NT.red[t-1,]$Delta.AIC) {} else {
      NT.red[t,]<-m
      t<-t+1
    }
  }
}


for (i in 1:nrow(TR)) { #i<-1
  m<-NT[i,]
  if (i == 1){
    TR.red<-m
    t<-2
  } else {
    if (m$Delta.AIC == TR.red[t-1,]$Delta.AIC) {} else {
      TR.red[t,]<-m
      t<-t+1
    }
  }
}
##ordered the list of mdels based on above output, TR.  Did by hand, 
## can't figure out how to order the list of models based on AIC.. Argh!
##use this so you can cut and paste and then remove quotes at beginning and end...
paste(NT.red$Model,collapse=", ")

NTMods<-list(YE.hn.Stage, YE.unif.cos, YE.hn.Depth.Stage, YE.hn, YE.unif.poly, YE.hn.Depth, YE.hr.cos, YE.hr.Stage, YE.hr, YE.hr.Depth.Stage, YE.hr.Depth, YE.unif.herm)

	NT.red$Det.Prob<-sapply(NTMods, function(x) summary(x)$dht$individuals$average.p)
	NT.red$Nhat<-sapply(NTMods, function(x) summary(x)$dht$individuals$N$Estimate)
	NT.red$Nhat.se<-sapply(NTMods, function(x) summary(x)$dht$individuals$N$se)
	NT.red$Nhat.cv<-sapply(NTMods, function(x) summary(x)$dht$individuals$N$cv)
	NT.red$Nhat.lcl<-sapply(NTMods, function(x) summary(x)$dht$individuals$N$lcl)
	NT.red$Nhat.ucl<-sapply(NTMods, function(x) summary(x)$dht$individuals$N$ucl)

	NT.red$Dhat<-sapply(NTMods, function(x) summary(x)$dht$individuals$D$Estimate)
	NT.red$Dhat.se<-sapply(NTMods, function(x) summary(x)$dht$individuals$D$se)
	NT.red$Dhat.cv<-sapply(NTMods, function(x) summary(x)$dht$individuals$D$cv)
	NT.red$Dhat.lcl<-sapply(NTMods, function(x) summary(x)$dht$individuals$D$lcl)
	NT.red$Dhat.ucl<-sapply(NTMods, function(x) summary(x)$dht$individuals$D$ucl)

paste(TR.red$Model,collapse=", ")

TrMods<-list(YE.hn.Stage, YE.unif.cos, YE.hn.Depth.Stage, YE.hn, YE.unif.poly, YE.hn.Depth, YE.hr.cos, YE.hr.Stage, YE.hr, YE.hr.Depth.Stage)

	TR.red$Det.Prob<-sapply(TrMods, function(x) summary(x)$dht$individuals$average.p)
	TR.red$Nhat<-sapply(TrMods, function(x) summary(x)$dht$individuals$N$Estimate)
	TR.red$Nhat.se<-sapply(TrMods, function(x) summary(x)$dht$individuals$N$se)
	TR.red$Nhat.cv<-sapply(TrMods, function(x) summary(x)$dht$individuals$N$cv)
	TR.red$Nhat.lcl<-sapply(TrMods, function(x) summary(x)$dht$individuals$N$lcl)
	TR.red$Nhat.ucl<-sapply(TrMods, function(x) summary(x)$dht$individuals$N$ucl)

	TR.red$Dhat<-sapply(TrMods, function(x) summary(x)$dht$individuals$D$Estimate)
	TR.red$Dhat.se<-sapply(TrMods, function(x) summary(x)$dht$individuals$D$se)
	TR.red$Dhat.cv<-sapply(TrMods, function(x) summary(x)$dht$individuals$D$cv)
	TR.red$Dhat.lcl<-sapply(TrMods, function(x) summary(x)$dht$individuals$D$lcl)
	TR.red$Dhat.ucl<-sapply(TrMods, function(x) summary(x)$dht$individuals$D$ucl)

### EXPORT MODEL LIST FOR MODEL AVERAGING AND BOOTSTRAPPING...
write.csv(NT.red, file = paste0("YE_density_estimation/Data/",Subd,"_",YEAR,"/",Subd,"_",YEAR,"_Bootstrap_Model_List_all.csv"))

write.csv(TR.red, file = paste0("YE_density_estimation/Data/",Subd,"_",YEAR,"/",Subd,"_",YEAR,"_Bootstrap_Model_List_Trunc.csv"))


#################################################################################
## Compare models visually...
###################################################
Mods<-TrMods
minAIC<-min(as.numeric(sapply(Mods,function(x)AIC(x)[2])))

AIC(Mods[[1]])

## plot detection 
par(mfrow = c(3,4),mar=c(4,4,4,1)+0.1,oma=c(5,5,5,1))

for (i in 1:length(Mods)) {		#i<-1
	plot(Mods[[i]],cex.axis=1.3, cex.main=50, ylab="",
		xlab=paste("deltaAIC=",round(AIC(Mods[[i]])[2]-minAIC,2)))
	mtext(paste(Mods[[i]]$call$key,Mods[[i]]$call$adjustment, sep=","),side=3,adj=-0.1,cex=1)
	mtext(Mods[[i]]$call$formula,side=3,adj=1,cex=1)
}
mtext("Distance (m)", side=1,outer=TRUE,cex=2)
mtext("Frequency", side=2,outer=TRUE,cex=2)
mtext("Detection Functions for 5% Truncated Models", side=3,outer=TRUE,cex=2)

##plot GOF 

par(mfrow = c(3,4),mar=c(4,4,4,1)+0.1,oma=c(5,5,5,1))

for (i in 1:length(Mods)) {		#i<-1
	gof_ds(Mods[[i]])[1]$dsgof$CvM$p
	mtext(paste(Mods[[i]]$call$key,Mods[[i]]$call$adjustment, sep=","),side=3,adj=-0.1,cex=1)
	mtext(Mods[[i]]$call$formula,side=3,adj=1,cex=1)
#	mtext(Mods[[i]])[1]$dsgof$CvM$p
}
mtext("GOF for 5% Truncated Models", side=3,outer=TRUE,cex=2)

## PLOT Dhats from different models...

GG<-NT.red[rev(order(NT.red$Delta.AIC, NT.red$Model)),]
GG<-TR.red[rev(order(TR.red$Delta.AIC, TR.red$Model)),]

ggplot(GG, aes(x=reorder(Model,Delta.AIC), y=Dhat))+geom_point()+
	geom_errorbar(aes(ymin=Dhat.lcl, ymax=Dhat.ucl),width=0.2)+
	xlab("Model")+
	ylab("NSEO Yelloweye Density")+
	theme(axis.text.x = element_text(angle=90, hjust=1))+
	scale_y_continuous(name="SSEO Yelloweye Density", breaks=seq(from=0, to=5000, by=100))#, 
		#labels, limits)	



######################################################################################
## Go to DSR_ROV_SEAK_basic_biomass_calc.R to calculate biomass from the best model
## Go to DSR_ROV_SEAK_dist_model_avg.R to get Dhat from model averaging the best models
## Go to DSR_ROV_SEAK_Full_bootstrap_biomass_estimation.R to model average and calculate biomass in bootstrap 
########################################################################################

