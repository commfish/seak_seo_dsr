#################################################################
## Examination of length and weight data for converting
## YE abundance estimates from Nhat out of distance analysis
## to biomass estimates
## July 20201, Phil Joy
#################################################################
library(tidyverse)
library(FSA)
library(fishmethods)
library(RColorBrewer)
library(boot)

setwd("D:/Groundfish Biometrics/Yelloweye/YE Code")

##not using survey data...
##Surv<-read.csv("Data/survey_bio_data.csv")
Port1<-read.csv("Data/port_sampling_bio_data.csv")
Port2<-read.csv("Data/port_sampling_bio_data2.csv")
ROV<-read_csv("Data/SSEO_2020_species.csv") %>% filter(Species == 145)

#str(Surv)
#unique(Surv[,1])
#unique(Surv$ï..Year)
#Surv$Year<-as.integer(Surv$ï..Year )
#unique(Surv$Year)
#hist(Surv$Length.Millimeters)
#min(Surv$Length.Millimeters, na.rm=TRUE)
#hist(Surv$Weight.Kilograms)
#plot(Surv$Weight.Kilograms ~ Surv$Length.Millimeters)

#unique(Surv$Length.Type)
#unique(Surv$Length.Type.Code)
#table(Surv$Length.Type.Code)
##get rid of those 3 fish measured wrong!
#nrow(Surv)
#Surv<-Surv[Surv$Length.Type.Code ==1,]

###clean ROV data
str(ROV)
hist(ROV$Length..mm.)
#get rid of non YEs

#get rid of the juvies...
nrow(ROV)
ROV.adult <- ROV %>% filter(Stage != "JV")
nrow(ROV.adult)

str(ROV.adult)

## get rid of unmeasured fish	#str(ROV.adult[,36])
ROV.L<-as.data.frame(ROV.adult[,c(1:2,4:34)])
ROV.L<-ROV.L[complete.cases(ROV.L[,"Length..mm."]),]
str(ROV.L)	#
hist(ROV.L$Horz..Dir...deg)

ROV.L$Length<-ROV.L$Length..mm.  
ROV.L$L.prec<-ROV.L$Precision..mm.
ROV.L$L.cv<-ROV.L$L.prec/ROV.L$Length
hist(ROV.L$L.cv)	#will remove big cv lengths for analysis...

##get rid of crappy measurements
ROV.L$badL<-ROV.L$Length-3*ROV.L$L.prec
ROV.good<-ROV.L[ROV.L$badL > 0,]

plot(ROV.L$Precision..mm. ~ ROV.L$Length..mm.)
hist(ROV.L$Precision..mm.)

nrow(ROV.L[ROV.L$Length..mm.<= 270,])/nrow(ROV.L)


###Clean Port
str(Port1)
str(Port2)
Port1$Year<-as.integer(as.character(Port1$ï..Year))
Port2$Year<-as.integer(as.character(Port2$ï..Year))

Port1<-subset(Port1, !is.na(Year))
Port2<-subset(Port2, !is.na(Year))
unique(Port1$Year)
unique(Port2$Year)
Port<-rbind(Port1,Port2)
unique(Port$Year)

str(Port)

#get rid of unpaired length and weight data
Port <- subset(Port,!is.na(Weight.Kilograms) & !is.na(Length.Millimeters))

##get rid of that weird system administrator stuff... 
#table(Port$Year)
#sum(is.na(Port$Year))
#nrow(Port)
#Port<-subset(Port, !is.na(Year))

## Select data from same management area, in this case SSEO
unique(Port$Groundfish.Management.Area.Code)
Port.SSEO<-Port[Port$Groundfish.Management.Area.Code == "SSEO",]
Port.CSEO<-Port[Port$Groundfish.Management.Area.Code == "CSEO",]
#Port.NSEI
#Port.SSEI
Port.NSEO<-Port[Port$Groundfish.Management.Area.Code == "NSEO",]
Port.EYKT<-Port[Port$Groundfish.Management.Area.Code == "EYKT",]

plot(density(Port.SSEO$Length.Millimeters, na.rm=TRUE),xlim=c(200,1000),
	xlab="Length", 
	col="goldenrod", main="Length Dist")
	polygon(density(Port.SSEO$Length.Millimeters, na.rm=TRUE), 
		col=adjustcolor("goldenrod",alpha.f=0.25), border="goldenrod")
	polygon(density(Port.CSEO$Length.Millimeters, na.rm=TRUE), 
		col=adjustcolor("orange",alpha.f=0.25), border="orange")
	polygon(density(Port.NSEO$Length.Millimeters, na.rm=TRUE), 
		col=adjustcolor("blue",alpha.f=0.25), border="blue")
	polygon(density(Port.EYKT$Length.Millimeters, na.rm=TRUE), 
		col=adjustcolor("forestgreen",alpha.f=0.25), border="forestgreen")
	polygon(density(ROV.L$Length..mm., na.rm=TRUE), 
		col=adjustcolor("red",alpha.f=0.25), border="red")

plot(density(Port.SSEO$Weight.Kilograms, na.rm=TRUE),#xlim=c(200,1000),
	xlab="Length", 
	col="goldenrod", main="Weight Dist")
	polygon(density(Port.SSEO$Weight.Kilograms, na.rm=TRUE), 
		col=adjustcolor("goldenrod",alpha.f=0.25), border="goldenrod")
	polygon(density(Port.CSEO$Weight.Kilograms, na.rm=TRUE), 
		col=adjustcolor("orange",alpha.f=0.25), border="orange")
	polygon(density(Port.NSEO$Weight.Kilograms, na.rm=TRUE), 
		col=adjustcolor("blue",alpha.f=0.25), border="blue")
	polygon(density(Port.EYKT$Weight.Kilograms, na.rm=TRUE), 
		col=adjustcolor("forestgreen",alpha.f=0.25), border="forestgreen")
#	polygon(density(ROV.L$Length..mm., na.rm=TRUE), 
#		col=adjustcolor("red",alpha.f=0.25), border="red")

Port21<-Port[Port$Year == 2021,]
nrow(Port21)
Port.SSEO.21<-Port.SSEO[Port.SSEO$Year == 2021,]
nrow(Port.SSEO.21)
Port.CSEO.21<-Port.CSEO[Port.CSEO$Year == 2021,]
nrow(Port.CSEO.21)
Port.NSEO.21<-Port.NSEO[Port.NSEO$Year == 2021,]
nrow(Port.NSEO.21)
Port.EYKT.21<-Port.EYKT[Port.EYKT$Year == 2021,]
nrow(Port.EYKT.21)

Port20<-Port[Port$Year == 2020,]
nrow(Port20)

Port.last3<-Port[Port$Year >= 2019,]
Port.SSEO.last3<-Port.SSEO[Port.SSEO$Year >= 2019,]
Port.NSEO.last3<-Port.NSEO[Port.NSEO$Year >= 2019,]
Port.CSEO.last3<-Port.CSEO[Port.CSEO$Year >= 2019,]
Port.EYKT.last3<-Port.EYKT[Port.EYKT$Year >= 2019,]

mean(Port.SSEO.last3$Weight.Kilograms)
sd(Port.SSEO.last3$Weight.Kilograms)
length(Port.SSEO.last3$Weight.Kilograms)

mean(Port.NSEO.last3$Weight.Kilograms)
sd(Port.NSEO.last3$Weight.Kilograms)
length(Port.NSEO.last3$Weight.Kilograms)

mean(Port.CSEO.last3$Weight.Kilograms)
sd(Port.CSEO.last3$Weight.Kilograms)
length(Port.CSEO.last3$Weight.Kilograms)

mean(Port.EYKT.last3$Weight.Kilograms)
sd(Port.EYKT.last3$Weight.Kilograms)
length(Port.EYKT.last3$Weight.Kilograms)

meanfun <- function(x, d) {
  return(mean(x[d]))
}

SSEO.w<-boot(Port.SSEO.last3$Weight.Kilograms, statistic=meanfun, R=1000)
SSEO.w
NSEO.w<-boot(Port.NSEO.last3$Weight.Kilograms, statistic=meanfun, R=1000)
NSEO.w
CSEO.w<-boot(Port.CSEO.last3$Weight.Kilograms, statistic=meanfun, R=1000)
CSEO.w
EYKT.w<-boot(Port.EYKT.last3$Weight.Kilograms, statistic=meanfun, R=1000)
EYKT.w

bo.ci<-boot.ci(bo, conf=0.95, type="bca")
summary(bo.ci)

#####################################################################
## Compare ROV lengths to Port lengths
#######################################################################
str(Port.SSEO.last3)
plot(density(Port.SSEO.last3$Length.Millimeters ), xlim=c(0,1000),
	xlab="Length", 
	col="goldenrod", main="Length Dist")
polygon(density(Port.SSEO.last3$Length.Millimeters ), col=adjustcolor("goldenrod",alpha.f=0.25))
polygon(density(ROV.L$Length), col=adjustcolor("navyblue",alpha.f=0.25))
polygon(density(ROV.good$Length), col=adjustcolor("yellow",alpha.f=0.5))

ROV.gt270<-ROV.good[ROV.good$Length >= 270,]
nrow(ROV.gt270)/nrow(ROV.good)
polygon(density(ROV.gt270$Length), col=adjustcolor("blue",alpha.f=0.25))
abline(v=mean(Port.SSEO.last3$Length.Millimeters), col="goldenrod")
abline(v=mean(ROV.L$Length), col="navyblue")
abline(v=mean(ROV.gt270$Length), col="blue")

mtext(paste("mean SSEO length = ",round(mean(Port.SSEO.last3$Length.Millimeters)),2), 
	side=3, col="goldenrod",adj=1,padj=1.5)
mtext(paste("mean ROV length = ",round(mean(ROV.L$Length)),2), 
	side=3, col="navyblue",adj=1,padj=3)
mtext(paste("mean ROV length gt 270 = ",round(mean(ROV.gt270$Length)),2), 
	side=3, col="blue",adj=1,padj=4.5)

############ For L:W Exam #############################################

min(Port$Length.Millimeters, na.rm=TRUE)
quantile(ROV$Length..mm., c(0.01,0.05,0.95,0.99), na.rm=TRUE)

####plot basic data
par(mfrow = c(2,1),mar=c(4,4,4,1)+0.1,oma=c(1,1,1,1))
plot(density(Port$Length.Millimeters, na.rm=TRUE), xlim=c(0,1000),
	xlab="Length", 
	col="goldenrod", main="Length Dist")
	polygon(density(Port$Length.Millimeters, na.rm=TRUE), 
		col=adjustcolor("goldenrod",alpha.f=0.25), border="goldenrod")
	abline(v=mean(Port$Length.Millimeters, na.rm=TRUE), col="goldenrod",lwd=2)
	polygon(density(ROV.L$Length..mm., na.rm=TRUE),
		col=adjustcolor("red",alpha.f=0.25), border="red")
	abline(v=mean(ROV.L$Length..mm., na.rm=TRUE), col="red",lwd=2)
	mtext("Port",side=3,adj=0.95 ,padj=3, col="goldenrod")
	mtext("ROV",side=3,adj=0.95 ,padj=4.5, col="red")
plot(Port$Weight.Kilograms ~ Port$Length.Millimeters, col="darkorange", pch=3, ylim=c(0,15),
	xlab="Lengh (mm)", ylab="Weight (kg)")

##Clean up and organize data for LW regressions
Port$logL <- log(Port$Length.Millimeters)
Port$logW <- log(Port$Weight.Kilograms)

Port<-Port[Port$logW > -10,]
str(Port)

#Surv.LW$logL <- log(Surv.LW$Length.Millimeters)
#Surv.LW$logW <- log(Surv.LW$Weight.Kilograms)

lm.Port<-lm(logW~logL,data=Port)
#lm.Surv<-lm(logW~logL,data=Surv)

fitPlot(lm.Port, xlab="log Length (mm)", ylab="log Weight (kg)", main ="")
fitPlot(lm.Surv, xlab="log Length (mm)", ylab="log Weight (kg)", main ="")

summary(lm.Port)
summary(lm.Surv)

######################################################################
## Calculate condition and remove outliers... 
## examine some outliers and get rid of them?? (CHECK WITH RHEA!)

#Port.LW$logW.to.logL<-Port.LW$logW/Port.LW$logL
Port$WtoL<-Port$Weight.Kilograms/(Port$Length.Millimeters^3)
quantile(Port$WtoL, c(0.001,0.01,0.99,0.999))

#Surv.LW$logW.to.logL<-Surv.LW$logW/Surv.LW$logL
#Surv.LW$WtoL<-Surv.LW$Weight.Kilograms/(Surv.LW$Length.Millimeters^3)

#hist(Port$logW.to.logL, breaks=100)
hist(Port$WtoL, breaks=100)

#hist(Surv.LW$logW.to.logL, breaks=100)
#hist(Surv.LW$WtoL, breaks=100)

## remove some outliers 
Port.noOut<-Port[Port$WtoL < quantile(Port$WtoL, c(0.9995)) 
		 & Port$WtoL > quantile(Port$WtoL, c(0.0005)),]
hist(Port.noOut$WtoL, breaks=100)
plot(Port$Weight.Kilograms ~ Port$Length.Millimeters, col="red", 
	pch=3, ylim=c(0,15), xlab="Lengh (mm)", ylab="Weight (kg)")
points(Port.noOut$Weight.Kilograms ~ Port.noOut$Length.Millimeters, col="blue", 
	pch=3)

#Surv.LW.noOut<-Surv.LW[Surv.LW$WtoL < quantile(Surv.LW$WtoL, c(0.999)) 
#		 & Surv.LW$WtoL > quantile(Surv.LW$WtoL, c(0.001)),]
#hist(Surv.LW.noOut$WtoL, breaks=100)
#plot(Surv.LW$Weight.Kilograms ~ Surv.LW$Length.Millimeters, col="red", 
#	pch=3, ylim=c(0,15), xlab="Lengh (mm)", ylab="Weight (kg)")
#points(Surv.LW.noOut$Weight.Kilograms ~ Surv.LW.noOut$Length.Millimeters, col="darkorange", 
#	pch=3, add=TRUE)

lm.Port<-lm(logW~logL,data=Port.noOut)
#lm.Surv<-lm(logW~logL,data=Surv.noOut)

summary(lm.Port)
#summary(lm.Surv)

plot(Port.noOut$logW~Port.noOut$logL, col="red")
abline(lm.Port, col="red")

#### Clean up Last 5 years data for most recent LW relationship data

Port.last5$WtoL<-Port.last5$Weight.Kilograms/(Port.last5$Length.Millimeters^3)
quantile(Port.last5$WtoL, c(0.001,0.01,0.99,0.999))

hist(Port.last5$WtoL, breaks=100)

## remove some outliers 
Port.last5.noOut<-Port.last5[Port.last5$WtoL < quantile(Port.last5$WtoL, c(0.999)) 
		 & Port.last5$WtoL > quantile(Port.last5$WtoL, c(0.001)),]
hist(Port.last5.noOut$WtoL, breaks=100)
plot(Port.last5$Weight.Kilograms ~ Port.last5$Length.Millimeters, col="red", 
	pch=3, ylim=c(0,15), xlab="Lengh (mm)", ylab="Weight (kg)")
points(Port.last5.noOut$Weight.Kilograms ~ Port.last5.noOut$Length.Millimeters, col="darkorange", 
	pch=3)

lm.Port.last5<-lm(logW~logL,data=Port.noOut)
summary(lm.Port.last5)

plot(Port.last5.noOut$logW ~ Port.last5.noOut$logL, col="red")
abline(lm.Port.last5, col="red")
nrow(Port.last5.noOut)


#######################################################################################
#### Quick analysis to see if LW relationship varies by other factors in the data...
#######################################################################################
str(Surv.LW.noOut)
unique(Surv$Year)
Surv.LW.noOut$Trip.No<-as.factor(Surv.LW.noOut$Trip.No)
boxplot(data = Surv.LW.noOut, WtoL ~ Year)
Y<-unique(Surv.LW.noOut$Year)
T<-unique(Surv.LW.noOut$Trip.No)

par(mfrow = c(2,1),mar=c(4,4,4,1)+0.1,oma=c(1,1,1,1))

plot(Surv.LW.noOut$logW~Surv.LW.noOut$logL, col="grey")
COLS<-brewer.pal(12,'Set2')
i<-1
for (y in Y){
	D<-Surv.LW.noOut[Surv.LW.noOut$Year == y,]
	points(D$logW~D$logL, col=COLS[i])
	abline(lm(D$logW~D$logL), col=COLS[i])
	i<-i+1
}
plot(Surv.LW.noOut$logW~Surv.LW.noOut$logL, col="grey")
i<-1
for (t in T){
	D<-Surv.LW.noOut[Surv.LW.noOut$Trip.No == t,]
	points(D$logW~D$logL, col=COLS[i])
	abline(lm(D$logW~D$logL), col=COLS[i])
	i<-i+1
}

unique(Surv.LW.noOut$ï..Year)
Surv.ex<-glm(data = Surv.LW.noOut, WtoL ~ Year)
summary(Surv.ex)

##PORT
str(Port.LW.noOut)
GMAC<-unique(Port.noOut$Groundfish.Management.Area.Code)
plot(Port.LW.noOut$logW~Port.LW.noOut$logL, col="grey")
i<-1
for (g in GMAC){
	D<-Port.LW.noOut[Port.LW.noOut$Groundfish.Management.Area.Code == g,]
	points(D$logW~D$logL, col=COLS[i])
	abline(lm(D$logW~D$logL), col=COLS[i])
	i<-i+1
}

PROJ<-unique(Port.noOut$Project)
plot(Port.noOut$logW~Port.noOut$logL, col="grey")
i<-1
for (p in PROJ){
	D<-Port.noOut[Port.noOut$Project == p,]
	points(D$logW~D$logL, col=COLS[i])
	abline(lm(D$logW~D$logL), col=COLS[i])
	i<-i+1
}

ST<-unique(Port.noOut$Sample.Type)
plot(Port.noOut$logW~Port.noOut$logL, col="grey")
i<-1
for (s in ST){
	D<-Port.noOut[Port.noOut$Sample.Type == s,]
	points(D$logW~D$logL, col=COLS[i])
	abline(lm(D$logW~D$logL), col=COLS[i])
	i<-i+1
}
nrow(Port.noOut[Port.noOut$Sample.Type == "Random",])
nrow(Port.noOut)
## **** Just use Random Samples **** ##
Port.noOut.Rand<-Port.noOut[Port.noOut$Sample.Type == "Random",]
##

Y<-unique(Port.noOut$Year)
COLS<-brewer.pal(length(Y),'YlOrRd')
plot(Port.noOut$logW~Port.noOut$logL, col="grey")
i<-1
for (y in Y){
	D<-Port.noOut[Port.noOut$Year == y,]
	points(D$logW~D$logL, col=COLS[i])
	abline(lm(D$logW~D$logL), col=COLS[i])
	i<-i+1
}
legend(x="bottomright",c("Tower/MR Esc. Est."),pch=c(18,19,5),
             col=c("red","blue","purple"),bty="n",text.col=c("red","blue","purple"))

Y<-unique(Port.noOut.Rand$Year)
COLS<-brewer.pal(length(Y),'YlOrRd')
plot(Port.noOut.Rand$logW~Port.noOut.Rand$logL, col="grey")
i<-1
for (y in Y){
	D<-Port.noOut.Rand[Port.noOut.Rand$Year == y,]
	points(D$logW~D$logL, col=COLS[i])
	abline(lm(D$logW~D$logL), col=COLS[i])
	i<-i+1
}


##couple of weird years... , I think bootstrapping relationship will account for uncertainty in true LW relationship...

################################################################################
#Save the Port data without outliers and reexamine avg weights from diferent areas for
#final analysis...

write.csv(Port.noOut, file = "Data/YE_LW_rel.csv")

##############################################################################
## Get weights from cleaned data 

Port.noOut.SSEO<-Port.noOut[Port.noOut$Groundfish.Management.Area.Code == "SSEO",]
Port.noOut.CSEO<-Port.noOut[Port.noOut$Groundfish.Management.Area.Code == "CSEO",]
Port.noOut.NSEO<-Port.noOut[Port.noOut$Groundfish.Management.Area.Code == "NSEO",]
Port.noOut.EYKT<-Port.noOut[Port.noOut$Groundfish.Management.Area.Code == "EYKT",]

plot(density(Port.noOut.SSEO$Length.Millimeters, na.rm=TRUE),xlim=c(200,1000),
	xlab="Length", 
	col="goldenrod", main="Length Dist")
	polygon(density(Port.noOut.SSEO$Length.Millimeters, na.rm=TRUE), 
		col=adjustcolor("goldenrod",alpha.f=0.25), border="goldenrod")
	polygon(density(Port.noOut.CSEO$Length.Millimeters, na.rm=TRUE), 
		col=adjustcolor("orange",alpha.f=0.25), border="orange")
	polygon(density(Port.noOut.NSEO$Length.Millimeters, na.rm=TRUE), 
		col=adjustcolor("blue",alpha.f=0.25), border="blue")
	polygon(density(Port.noOut.EYKT$Length.Millimeters, na.rm=TRUE), 
		col=adjustcolor("forestgreen",alpha.f=0.25), border="forestgreen")
	polygon(density(ROV.L$Length..mm., na.rm=TRUE), 
		col=adjustcolor("red",alpha.f=0.25), border="red")

plot(density(Port.noOut.SSEO$Weight.Kilograms, na.rm=TRUE),#xlim=c(200,1000),
	xlab="Length", 
	col="goldenrod", main="Weight Dist")
	polygon(density(Port.noOut.SSEO$Weight.Kilograms, na.rm=TRUE), 
		col=adjustcolor("goldenrod",alpha.f=0.25), border="goldenrod")
	polygon(density(Port.noOut.CSEO$Weight.Kilograms, na.rm=TRUE), 
		col=adjustcolor("orange",alpha.f=0.25), border="orange")
	polygon(density(Port.noOut.NSEO$Weight.Kilograms, na.rm=TRUE), 
		col=adjustcolor("blue",alpha.f=0.25), border="blue")
	polygon(density(Port.noOut.EYKT$Weight.Kilograms, na.rm=TRUE), 
		col=adjustcolor("forestgreen",alpha.f=0.25), border="forestgreen")


Port.noOut21<-Port.noOut[Port.noOut$Year == 2021,]
nrow(Port.noOut21)
Port.noOut.SSEO.21<-Port.noOut.SSEO[Port.noOut.SSEO$Year == 2021,]
nrow(Port.noOut.SSEO.21)
Port.noOut.CSEO.21<-Port.noOut.CSEO[Port.noOut.CSEO$Year == 2021,]
nrow(Port.noOut.CSEO.21)
Port.noOut.NSEO.21<-Port.noOut.NSEO[Port.noOut.NSEO$Year == 2021,]
nrow(Port.noOut.NSEO.21)
Port.noOut.EYKT.21<-Port.noOut.EYKT[Port.noOut.EYKT$Year == 2021,]
nrow(Port.noOut.EYKT.21)

Port.noOut20<-Port.noOut[Port.noOut$Year == 2020,]
nrow(Port.noOut20)

Port.noOut.last3<-Port.noOut[Port.noOut$Year >= 2019,]
Port.noOut.SSEO.last3<-Port.noOut.SSEO[Port.noOut.SSEO$Year >= 2019,]
Port.noOut.NSEO.last3<-Port.noOut.NSEO[Port.noOut.NSEO$Year >= 2019,]
Port.noOut.CSEO.last3<-Port.noOut.CSEO[Port.noOut.CSEO$Year >= 2019,]
Port.noOut.EYKT.last3<-Port.noOut.EYKT[Port.noOut.EYKT$Year >= 2019,]

mean(Port.noOut.SSEO.last3$Weight.Kilograms)
sd(Port.noOut.SSEO.last3$Weight.Kilograms)
length(Port.noOut.SSEO.last3$Weight.Kilograms)

mean(Port.noOut.NSEO.last3$Weight.Kilograms)
sd(Port.noOut.NSEO.last3$Weight.Kilograms)
length(Port.noOut.NSEO.last3$Weight.Kilograms)

mean(Port.noOut.CSEO.last3$Weight.Kilograms)
sd(Port.noOut.CSEO.last3$Weight.Kilograms)
length(Port.noOut.CSEO.last3$Weight.Kilograms)

mean(Port.noOut.EYKT.last3$Weight.Kilograms)
sd(Port.noOut.EYKT.last3$Weight.Kilograms)
length(Port.noOut.EYKT.last3$Weight.Kilograms)

meanfun <- function(x, d) {
  return(mean(x[d]))
}

SSEO.w<-boot(Port.noOut.SSEO.last3$Weight.Kilograms, statistic=meanfun, R=1000)
SSEO.w
NSEO.w<-boot(Port.noOut.NSEO.last3$Weight.Kilograms, statistic=meanfun, R=1000)
NSEO.w
CSEO.w<-boot(Port.noOut.CSEO.last3$Weight.Kilograms, statistic=meanfun, R=1000)
CSEO.w
EYKT.w<-boot(Port.noOut.EYKT.last3$Weight.Kilograms, statistic=meanfun, R=1000)
EYKT.w

bo.ci<-boot.ci(bo, conf=0.95, type="bca")
summary(bo.ci)

###############################################################################
# backtransfoming function for going from avg. length to average weight...

lm.Port<-lm(logW~logL,data=Port.noOut)
lm.Surv<-lm(logW~logL,data=Surv.LW.noOut)

#correction factor for back transforming predicted weights from given length ...
syx<-summary(lm.Port)$sigma
(cf<-exp((syx^2)/2))
#back transform 600mm YE based on Port samples...
meanL.port<-mean(Port.LW.noOut$Length.Millimeters)
(pred.log.port<-predict(lm.Port,data.frame(logL=log(meanL.port)), interval = "c"))
( bias.pred.orig.port <- exp(pred.log.port) )
(pred.orig.port<-cf*bias.pred.orig.port)

#correction factor for back transforming predicted weights from given length ...
syx<-summary(lm.Surv)$sigma
(cf<-exp((syx^2)/2))
#back transform 600mm YE based on Port samples...
meanL.surv<-mean(Surv.LW.noOut$Length.Millimeters)
(pred.log.surv<-predict(lm.Surv,data.frame(logL=log(meanL.port)), interval = "c"))
( bias.pred.orig.surv <- exp(pred.log.surv) )
(pred.orig.surv<-cf*bias.pred.orig.surv)
pred.orig[1]	#point estimate of predicted weight... 

pred.orig.port[1]/pred.orig.surv[1]
#1.6% difference in biomass between two data sources... 




























