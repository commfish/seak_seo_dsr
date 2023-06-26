#################################################################
## Examination of length and weight data from SEAK port samples and ROV samples
## most of this script is exploratory in nature
## Most important output is a csv file with Port length:weight samples with outliers removed
## August 20201, Phil Joy
#################################################################
library(tidyverse)
library(FSA)
library(fishmethods)
library(RColorBrewer)
library(boot)

#
setwd("D:/Groundfish Biometrics/Yelloweye/YE Code")

## Load the port samples that were downloaded from oceansAK

Port1<-read.csv("port_sampling_bio_data.csv")
Port2<-read.csv("port_sampling_bio_data2.csv")
                    
#Load the ROV data for YE 
ROV<-read_csv("NSEO_2018_species.csv") %>% filter(Species == 145)
		#("Data/SSEO_2020_species.csv") %>% filter(Species == 145)
###clean ROV data
#get rid of the juvies...
nrow(ROV)
ROV.adult <- ROV %>% filter(Stage != "JV")
nrow(ROV.adult)
str(ROV.adult)
## rename columns if necessary...
colnames(ROV.adult)
colnames(ROV.adult)[11]<-"Mid.X..mm."
colnames(ROV.adult)[12]<-"Mid.Y..mm."
colnames(ROV.adult)[13]<-"Mid.Z..mm."
colnames(ROV.adult)[4]<-"Fish.L"
colnames(ROV.adult)[5]<-"Precision..mm."
colnames(ROV.adult)[20]<-"Tansect.Number"

## get rid of unmeasured fish	#str(ROV.adult[,36])
ROV.L<-as.data.frame(ROV.adult[,c(1:2,4:28)])
ROV.L<-ROV.L[complete.cases(ROV.L[,"Fish.L"]),]
str(ROV.L)	#
hist(ROV.L$Horz..Dir...deg)
ROV.L$Length<-ROV.L$Fish.L			#Length..mm.  
ROV.L$L.prec<-ROV.L$Precision..mm.
ROV.L$L.cv<-ROV.L$L.prec/ROV.L$Length
hist(ROV.L$L.cv)	#will remove big cv lengths for analysis...

##get rid of crappy measurements
ROV.good<-ROV.L[ROV.L$L.cv < 0.4,]

nrow(ROV.L[ROV.L$Fish.L<= 270,])/nrow(ROV.L)

###Clean Port Files
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

#add columns with log weight and log length 
Port$logL <- log(Port$Length.Millimeters)
Port$logW <- log(Port$Weight.Kilograms)

## Select data from same management area, in this case SSEO
unique(Port$Groundfish.Management.Area.Code)
Port.SSEO<-Port[Port$Groundfish.Management.Area.Code == "SSEO",]
Port.CSEO<-Port[Port$Groundfish.Management.Area.Code == "CSEO",]
Port.NSEO<-Port[Port$Groundfish.Management.Area.Code == "NSEO",]
Port.EYKT<-Port[Port$Groundfish.Management.Area.Code == "EYKT",]

## look at length and weights from the different subdistricts... 
COLS<-brewer.pal(8,'Set1')

par(mfrow=c(2,1),mar=c(4,4,4,1),ps=12)
plot(density(Port.SSEO$Length.Millimeters, na.rm=TRUE),xlim=c(100,1000),
	xlab="Length (mm)", 
	col="goldenrod", main="Length Dist")
	polygon(density(Port.SSEO$Length.Millimeters, na.rm=TRUE), 
		col=adjustcolor(COLS[1],alpha.f=0.25), border=COLS[1])
	polygon(density(Port.CSEO$Length.Millimeters, na.rm=TRUE), 
		col=adjustcolor(COLS[2],alpha.f=0.25), border=COLS[2])
	polygon(density(Port.NSEO$Length.Millimeters, na.rm=TRUE), 
		col=adjustcolor(COLS[3],alpha.f=0.25), border=COLS[3])
	polygon(density(Port.EYKT$Length.Millimeters, na.rm=TRUE), 
		col=adjustcolor(COLS[4],alpha.f=0.25), border=COLS[4])
	polygon(density(ROV.L$Length, na.rm=TRUE), 
		col=adjustcolor(COLS[5],alpha.f=0.25), border=COLS[5])
	abline(v=mean(Port.SSEO$Length.Millimeters), col=COLS[1])
	abline(v=mean(Port.CSEO$Length.Millimeters), col=COLS[2])
	abline(v=mean(Port.NSEO$Length.Millimeters), col=COLS[3])
	abline(v=mean(Port.EYKT$Length.Millimeters), col=COLS[4])
	abline(v=mean(ROV.L$Length), col=COLS[5])
	legend(x="topright",legend=c("SSEO Port","CSEO Port","NSEO Port","EYKT Port","SSEO ROV"), 
		text.col=COLS)	#xjust=0, yjust=0)

plot(density(Port.SSEO$Weight.Kilograms, na.rm=TRUE),#xlim=c(200,1000),
	xlab="Weight (kg)", 
	col="goldenrod", main="Weight Dist")
	polygon(density(Port.SSEO$Weight.Kilograms, na.rm=TRUE), 
		col=adjustcolor(COLS[1],alpha.f=0.25), border=COLS[1])
	polygon(density(Port.CSEO$Weight.Kilograms, na.rm=TRUE), 
		col=adjustcolor(COLS[2],alpha.f=0.25), border=COLS[2])
	polygon(density(Port.NSEO$Weight.Kilograms, na.rm=TRUE), 
		col=adjustcolor(COLS[3],alpha.f=0.25), border=COLS[3])
	polygon(density(Port.EYKT$Weight.Kilograms, na.rm=TRUE), 
		col=adjustcolor(COLS[4],alpha.f=0.25), border=COLS[4])
	abline(v=mean(Port.SSEO$Weight.Kilograms), col=COLS[1])
	abline(v=mean(Port.CSEO$Weight.Kilograms), col=COLS[2])
	abline(v=mean(Port.NSEO$Weight.Kilograms), col=COLS[3])
	abline(v=mean(Port.EYKT$Weight.Kilograms), col=COLS[4])
	legend(x="topright",legend=c("SSEO Port","CSEO Port","NSEO Port","EYKT Port"), 
		text.col=COLS)	#xjust=0, yjust=0)


##consider different groupings for analysis...

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

##calcualte mean weights and the variance/sd/cv of the mean ...

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

bo.ci<-boot.ci(SSEO.w, conf=0.95, type="bca")
summary(bo.ci)

#####################################################################
## Compare ROV lengths to Port lengths
#######################################################################
str(Port.SSEO.last3)
plot(density(Port.NSEO.last3$Length.Millimeters ), xlim=c(0,1000),
	xlab="Length", 
	col="goldenrod", main="Length Dist")
polygon(density(Port.NSEO.last3$Length.Millimeters ), col=adjustcolor("goldenrod",alpha.f=0.25))
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

######################################################################
## Calculate condition and remove outliers... ; outliers = realy, really fat and really, really skinny fish
##	i.e., they seem unrealistic... 
## examine some outliers and get rid of them?? (CHECK WITH RHEA!)

# Create a column of conditions (K = W/L^3)
Port$WtoL<-Port$Weight.Kilograms/(Port$Length.Millimeters^3)
quantile(Port$WtoL, c(0.001,0.01,0.99,0.999))
hist(Port$WtoL, breaks=100)

## remove some outliers 
Port.noOut<-Port[Port$WtoL < quantile(Port$WtoL, c(0.999)) 
		 & Port$WtoL > quantile(Port$WtoL, c(0.001)),]
hist(Port.noOut$WtoL, breaks=100)

## look at the outliers and see if it makes sense... 
plot(Port$Weight.Kilograms ~ Port$Length.Millimeters, col="red", 
	pch=3, ylim=c(0,15), xlab="Lengh (mm)", ylab="Weight (kg)")
points(Port.noOut$Weight.Kilograms ~ Port.noOut$Length.Millimeters, col="blue", 
	pch=3)

lm.Port<-lm(logW~logL,data=Port.noOut)
summary(lm.Port)

plot(Port.noOut$logW~Port.noOut$logL, col="red")
abline(lm.Port, col="red")

#### Clean up Last 5 years data for most recent LW relationship data

#Port.last5$WtoL<-Port.last5$Weight.Kilograms/(Port.last5$Length.Millimeters^3)
#quantile(Port.last5$WtoL, c(0.001,0.01,0.99,0.999))

#hist(Port.last5$WtoL, breaks=100)

#######################################################################################
#### Quick analysis to see if LW relationship varies by other factors in the data...
#######################################################################################
##PORT
str(Port.noOut)

par(mfrow=c(3,1),mar=c(6,6,4,1))

GMAC<-unique(Port.noOut$Groundfish.Management.Area.Code)	#GMAC[3]
GMAC<-GMAC[c(4:7)]
plot(Port.noOut$logW~Port.noOut$logL, col="white", xlab="log Length", ylab="log Weight", 
	main = "Port Samples", cex.lab=1.8, cex.axis=1.5)
i<-1
for (g in GMAC){
	D<-Port.noOut[Port.noOut$Groundfish.Management.Area.Code == g,]
	points(D$logW~D$logL, col=COLS[i], pch=20, cex=2)
	abline(lm(D$logW~D$logL), col=COLS[i])
	i<-i+1
}
legend(x="bottomright",legend=GMAC, cex=1.5,
		text.col=COLS)	

PROJ<-unique(Port.noOut$Project)
plot(Port.noOut$logW~Port.noOut$logL, col="grey", xlab="log Length", ylab="log Weight", 
	main="Survey Type", cex.lab=1.8, cex.axis=1.5)
i<-1
for (p in PROJ){
	D<-Port.noOut[Port.noOut$Project == p,]
	points(D$logW~D$logL, col=COLS[i], pch=20, cex=2)
	abline(lm(D$logW~D$logL), col=COLS[i])
	i<-i+1
}
legend(x="bottomright",legend=PROJ,cex=1.25, bg=adjustcolor("grey",alpha.f=0.25),
		text.col=COLS)

ST<-unique(Port.noOut$Sample.Type)
plot(Port.noOut$logW~Port.noOut$logL, xlab="log Length", ylab="log Weight", col="grey", 
	main = "Sample Type", cex.lab=1.8, cex.axis=1.5)
i<-1
for (s in ST){
	D<-Port.noOut[Port.noOut$Sample.Type == s,]
	points(D$logW~D$logL, col=COLS[i], pch=20, cex=2)
	abline(lm(D$logW~D$logL), col=COLS[i])
	i<-i+1
}
legend(x="bottomright",legend=ST, cex=1.5,
		text.col=COLS)



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

write.csv(Port.noOut, file = "NSEO_YE_LW_rel.csv")

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
	polygon(density(ROV.L$Length, na.rm=TRUE), 
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
#            



























