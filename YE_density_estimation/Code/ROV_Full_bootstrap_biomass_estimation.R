################################################################################################################
## Full Bootstrap Analysis to Estimate Biomass from best Distance models (AIC based) 
## Resamples distance transects, fits models and picks best model in each bootstrap sample
## Also resamples the ROV lengths 
##	Estimates proportion greater than cutoff length that are detected in ROV surveys
## Also resamples the CF weights
## For each bootstrap sample best model used to estimate density, then
##	biomass is calculated by combining density estimate with average weight (from ROV and Port samples), area
##	covered, and proportion of the cutoff length
##  In 2021, this proceedure did not reveal any differences based on simpler analysis
##	However, the proportion estimate was not normally distributed and this proceedure may be useful in the 
##	future if that proportion value becomes significant in future surveys
##  August 2021
##  Phil Joy
#################################################################################################################
library(tictoc)
library(tidyverse)
library(Distance)
library(boot)
library(RColorBrewer)

##Set working directory
setwd("D:/Groundfish Biometrics/Yelloweye/YE Code")

################################################################################################
## 2 Functions you will need to run this:

## 1) this function is for backtransforming lengths (length) to weights using a linear model (LM) that comes from
## any combination of log weights~log lengths.
 
w.from.l<-function (LM, length){
	syx<-summary(LM)$sigma
	cf<-exp((syx^2)/2)
	pred.log<-predict(LM,data.frame(logL=log(length)), interval = "c")
	bias.pred.orig <- exp(pred.log) 
	pred.orig<-cf*bias.pred.orig
	pred.orig[1]
}

## 2) function for getting bootstrap mean and variance of that mean for avg weights

meanfun <- function(x, d) {
  return(mean(x[d]))
}

##############################################################################################
## Load and adjust files that were created with other scripts in this project: 
##############################################################################################

## Load the processed distance data file created in "DSR_ROV_SEAK_ROV_processing.R"

DAT<-read.csv("Data/SSEO_distance_data_GIStran_for_analysis.csv")

## Get ROV length data for estimating weights for biomass calculations
## Load the ROV species data and get rid of the juvenile fish that were detected
 
ROV<-read_csv("Data/SSEO_2020_species.csv") %>% filter(Species == 145)
ROV.adult <- ROV %>% filter(Stage != "JV")
## get rid of unmeasured fish	as well as columns that are just in the way...
ROV.L<-as.data.frame(ROV.adult[,c(1:2,4:34)])
ROV.L<-ROV.L[complete.cases(ROV.L[,"Length..mm."]),]
##make some simpler column names
ROV.L$Length<-ROV.L$Length..mm. 
ROV.L$L.prec<-ROV.L$Precision..mm. 
## get rid of imprecise measurements and make a new data frame with the good measurements...
ROV.L$rel.Prec<-ROV.L$L.prec/ROV.L$Length
ROV.good<-ROV.L[ROV.L$rel.Prec < 0.4,]	
nrow(ROV.L)
nrow(ROV.good)

## get summarized port sampling data created with script "DSR_ROV_SEAK_LWexam.R"
YE.LW<-read.csv("Data/YE_LW_rel.csv")
str(YE.LW)
## select the data you want to use (years, or multiple years) and the subdistrict you want to use...
## most recent data for select area (2021 = SSEO)
YE.latest<-YE.LW[YE.LW$Year >= 2019 & YE.LW$Groundfish.Management.Area.Code == "SSEO",]
nrow(YE.latest)		#colnames(YE.latest)
## get rid of some superfluous columns
YE.latest<-YE.latest[,-(c(1:16,18:21,23:36))]

### get the list of your models created in "DSR_ROV_SEAK_base_distance_modelling.R"
## these are the models that the bootstrap will consider with each replication..

G<-read.csv("Data/Bootstrap_Model_List.csv")

##########################################################################################
## Here are the objects that the bootstrap will call on: 
##########################################################################################

## Create you W:L linear model for estiming ROV weights from lengths...
LW.lm<-lm(logW~logL,data=YE.latest)

## Determine what cutoff ROV length you want to consider for estimating biomass
## here it is set to exclude ROV fish that are less than 270 mm 
CUTOFF<-270

## Establish conversion units for going from meters to square kilometers... 
CU<-convert_units(distance_units="Meter", effort_units="Meter", 
	area_units="Square kilometer")

## set habitat area size (NSEO = 442, SSEO = 1056, CSEO = 1661, EYKT = 739)
HA<-1056	
		
## From your list of models, pick the ones you want to average in the bootstrap
## Pick the models with deltaAIC less than four.  

dAIC_lt4<-G[G$Delta.AIC <= 4,]
nrow(dAIC_lt4)

## state how many replicates you want to run
nboot<-5000

## determine how many transects there are for resampling in the bootstrap...
smpls<-unique(DAT$Sample.Label)

#######################################################################################
## Run the bootstrap! 
#######################################################################################

Model.Boot<-data.frame(NA)									#output file

tic("boot")												#timer
pb = txtProgressBar(min = 0, max = length(nboot), initial = 0, style=3) 	#progress bar

for (r in 1:nboot){										#r<-1
## Resample the transects with replacement...
	for (s in 1:length(smpls)){
		rn<-sample(unique(DAT$Sample.Label), size=1,replace=TRUE)
		R1<-DAT[DAT$Sample.Label == rn,]
		R1$Orig.Sample.Label<-R1$Sample.Label
		R1$Sample.Label<-s
		if (s == 1) {
		R2<-R1
		} else {
		R2<-rbind(R2,R1)
		}
	}	#one replicated data set complete...
	
##Fit candidate models to the replicated data set... 
	Mlist<-list()
	for (i in 1:nrow(dAIC_lt4)){								#i<-1
## for each candidate model, extract the KEY, FORMULA and ADJUSTMENT terms to run the distance model
		AA<-dAIC_lt4[i,]									#str(AA)
		if (AA$Key.function == "Half-normal") {
			KEY<-"hn"
		} else if (AA$Key.function == "Hazard-rate") {
			KEY<-"hr"
		} else {KEY<-"unif"}

		if (KEY == "unif") {
			FORM<-~1
		} else {
			FORM<-AA$Formula						#str(FORM)
		}

		AA$Model<-as.character(AA$Model)
		HERM<-grepl("herm",AA$Model)
		POLY<-grepl("poly",AA$Model)
		COS<-grepl("cos",AA$Model)
		if (POLY == TRUE){
			ADJ<-"poly"
		} else if (COS == TRUE){
			ADJ<-"cos"
		} else if (HERM == TRUE){
			ADJ<-"herm"
		} else {
			ADJ<-NULL }
		
## Depending on what covariates are included, run the appropriate model, M
	## invisible(capture.output(...)) is just a wrapper to restrain printing output in the loop
	## try(...) allows the for loop to keep running if you get a bad model that doesn't converge...
	
		if (FORM == "~Depth") {
			try(invisible(capture.output(M <- ds(R2, key = KEY, adjustment = ADJ, transect="line", formula=~Depth,
                  convert.units = CU, truncation=c(right="5%"), quiet=TRUE),type="message")),silent=TRUE)
		} else if (FORM == "~Stage") {
			try(invisible(capture.output(M <- ds(R2, key = KEY, adjustment = ADJ, transect="line", formula=~Stage,
                  convert.units = CU, truncation=c(right="5%"), quiet=TRUE),type="message")),silent=TRUE)
		} else if (FORM == "~Depth + Stage") {
			try(invisible(capture.output(M <- ds(R2, key = KEY, adjustment = ADJ, transect="line", formula=~Depth + Stage,
                  convert.units = CU, truncation=c(right="5%"), quiet=TRUE),type="message")),silent=TRUE)
		} else {
			try(invisible(capture.output(M <- ds(R2, key = KEY, adjustment = ADJ, transect="line", formula=FORM,
                  convert.units = CU, truncation=c(right="5%"), quiet=TRUE),type="message")),silent=TRUE)
		}
		Mlist[[i]]<-M		
		assign(paste0("M",i),M)
	}	#finish loop running all cadidate models on replicated data set r

## select top model from list of candadate models	
	list.c <- sapply(Mlist, function(x) AIC(x)[2]> 0)
	out.list  <- Mlist[list.c]

	mAIC<-sapply(out.list,function(x) AIC(x)[2])
	AIC.vect<-unlist(mAIC)
	
	list.c.2 <- sapply(out.list, function(x) AIC(x)[2]==min(AIC.vect))
	out.list.2  <- out.list[list.c.2]

## extract the data we want from the best model	
	try(RANK<-summarize_ds_models(out.list.2[[1]]))	#this is ussually used to rank a series of model but I used it here 
									#for the top model as it was convenient for the code I had

	if (length(RANK)<= 1) { 				#"if else" to skip if RANK produces nothing because of bad model
	} else {
		#if RANK is good, extract all the data we want out of the bootstrap/// 
		names(RANK)[names(RANK) == "Model"]<-"Full.Model.Code"
		RANK$Model<-substr(gsub('[\\{}]','',RANK$Full.Model.Code),7,50)
		names(RANK)[names(RANK) == "C-vM p-value"]<-"CvM_pvalue"
		names(RANK)[names(RANK) == "$\\hat{P_a}$"]<-"P_a.hat"
		names(RANK)[names(RANK) == "se($\\hat{P_a}$)"]<-"se.P_a.hat"
		names(RANK)[names(RANK) == "$\\Delta$AIC"]<-"Delta.AIC"
		names(RANK)[names(RANK) == "Key function"]<-"Key_function"

		RANK$CoveredArea<-sapply(out.list.2[1], function(x) summary(x)$dht$individuals$summary$CoveredArea)
		RANK$Det.Prob<-sapply(out.list.2[1], function(x) summary(x)$dht$individuals$average.p)
		
		RANK$Dhat<-sapply(out.list.2[1], function(x) summary(x)$dht$individuals$D$Estimate)
		RANK$Dhat.se<-sapply(out.list.2[1], function(x) summary(x)$dht$individuals$D$se)
		RANK$Dhat.cv<-sapply(out.list.2[1], function(x) summary(x)$dht$individuals$D$cv)
		RANK$Dhat.lcl<-sapply(out.list.2[1], function(x) summary(x)$dht$individuals$D$lcl)
		RANK$Dhat.ucl<-sapply(out.list.2[1], function(x) summary(x)$dht$individuals$D$ucl)

### Calculate Biomass from the best model in this bootstrap sample... 
	# Resample ROV lengths with replacement...
		ROV.rand<-sample(ROV.good$Length, size=length(ROV.good$Length),replace=TRUE)
	#get the mean length from that sample and record mean and var
		RANK$mean.l.all<-mean(ROV.rand)
		RANK$var.l.all<-var(ROV.rand)
	#get and record the proportion of lengths greater than CUTOFF value... 
		Prop.gt<-length(ROV.rand[ROV.rand>CUTOFF])/length(ROV.rand)
		RANK$Prop.gt<-Prop.gt
		Prop.var<-Prop.gt*(1-Prop.gt)/length(ROV.rand)		#sqrt(Prop.gt*(1-Prop.gt)/length(ROV.rand))
		RANK$Prop.var<-Prop.var
	#get mean length of fish greater than CUTOFF value for biomass calculations... 
		ROV.gt.co<-ROV.rand[ROV.rand>CUTOFF]
		RANK$mean.l.gtco<-mean(ROV.gt.co)
		RANK$var.l.gtco<-var(ROV.gt.co)
	#weigh each ROV fish using the conversion function we set up at the beginning of this script... 
		Ws<-sapply(ROV.rand[ROV.rand > 0], function(x)w.from.l(LM=LW.lm, length=x))	#plot(Ws~ROV.rand[ROV.rand>CUTOFF])
		Ws.gtco<-sapply(ROV.gt.co, function(x)w.from.l(LM=LW.lm, length=x))	#plot(Ws~ROV.rand[ROV.rand>CUTOFF])
	#Get and record the mean weight of ROV fish, both all fish and just fish greater than your cutoff length... 
	# also we will get the variance of that mean using a small bootstrap... 
		RANK$mean.W<-mean(Ws)	
			mW.b<-boot(Ws,statistic=meanfun, R=500)				
		RANK$mean.W.var<-sd(mW.b$t)^2		
		
		RANK$mean.W.gtco<-mean(Ws.gtco)
			mW.b2<-boot(Ws.gtco,statistic=meanfun, R=500)					
		RANK$mean.W.gtco.var<-sd(mW.b2$t)^2			

	### Resample latest port data with replacement and get and record the mean weight of those fish...  		
		CFsample<-sample(YE.latest$Weight.Kilograms, size=length(YE.latest$Weight.Kilograms ), replace=TRUE)	
		CF.avgW<-mean(CFsample)
			mW.b3<-boot(CFsample,statistic=meanfun, R=500)	
		var.CF.avgW<-sd(mW.b3$t)^2		#var(CFsample)
		RANK$CF.mean.W<-CF.avgW
		RANK$CF.mean.W.var<-var.CF.avgW

	### Calculate biomass, its variance, and lower 90% CI for this bootstrap sample and record it... 
		RANK$Biomass.ROV.a<-RANK$Dhat*HA*Prop.gt*RANK$mean.W
		RANK$Biomass.ROV.gt<-RANK$Dhat*HA*Prop.gt*RANK$mean.W.gtco
		RANK$Biomass.CFW<-RANK$Dhat*HA*Prop.gt*CF.avgW

		RANK$var.Biomass.ROV.a<-((RANK$Dhat.se^2+RANK$Dhat^2)*(Prop.var+Prop.gt^2)*(RANK$mean.W.var+RANK$mean.W^2)-
			(RANK$Dhat^2)*(Prop.gt^2)*(RANK$mean.W^2))*HA^2
		RANK$ROV.a.lowCI<-(RANK$Biomass.ROV.a-1.645*sqrt(RANK$var.Biomass.ROV.a))

		RANK$var.Biomass.ROV.gt<-((RANK$Dhat.se^2+RANK$Dhat^2)*(Prop.var+Prop.gt^2)*(RANK$mean.W.gtco.var+RANK$mean.W.gtco^2)-
			(RANK$Dhat^2)*(Prop.gt^2)*(RANK$mean.W.gtco^2))*HA^2
		RANK$ROV.gt.lowCI<-(RANK$Biomass.ROV.gt-1.645*sqrt(RANK$var.Biomass.ROV.gt))

		RANK$var.Biomass.CFW<-((RANK$Dhat.se^2+RANK$Dhat^2)*(Prop.var+Prop.gt^2)*(var.CF.avgW+CF.avgW^2)-
			(RANK$Dhat^2)*(Prop.gt^2)*(CF.avgW^2))*HA^2
		RANK$CF.lowCI<-(RANK$Biomass.CFW-1.645*sqrt(RANK$var.Biomass.CFW))
	## convert kg to metric tons...
		RANK$Biomass.ROV.all<-RANK$Biomass.ROV.a/1000
		RANK$Biomass.ROV.gt270<-RANK$Biomass.ROV.gt/1000
		RANK$Biomass.CF<-RANK$Biomass.CFW/1000

		RANK$Biomass.ROV.all.loCI<-RANK$ROV.a.lowCI/1000
		RANK$Biomass.ROV.gt270.loCI<-RANK$ROV.gt.lowCI/1000
		RANK$Biomass.CF.loCI<-RANK$CF.lowCI/1000

	## Add this bootstrap sample to your list of bootstrap samples, recorded as Model.Boot...
		if(r == 1){
		Model.Boot<-RANK[1,]
		} else {
		Model.Boot<-rbind(Model.Boot,RANK[1,])
		}
	}

# Erase RANK and Mlist to avoid confusion when troubleshooting 
	RANK<-NA
	Mlist<-NA
	M<-NA
setTxtProgressBar(pb,r/nboot)
}

## how long it took to run (divide by 3600 to get hours...)
toc()	

##create model names that are useful for plotting... 
Model.Boot$Unique.Mods<-paste(Model.Boot$Key_function, Model.Boot$Formula)

## rename output for saving and plotting results
MB<-Model.Boot
write.csv(MB, file = "Output/Model_Avg_biomass_relP_p4_5K.csv")

######################################################################################
## Examine bootstrap
######################################################################################
## Reload your saved bootstrap output if you didn't just run it... 
#	MB<-read.csv("Output/Model_Avg_biomass_simpleboot_GIStran_5K.csv")
nrow(MB)
head(MB)

mean(MB$Dhat)
se(MB$Dhat)
var(MB$Dhat)
sqrt(var(MB$Dhat))

## Look at GOF results and get rid of worst fits... 
min(MB$CvM_pvalue)
MB<-MB[MB$CvM_pvalue > 0.1,]
nrow(MB)

## Make a table of Model frequency in the bootstrap (this is effectively your weights)
TAB<-table(MB$Unique.Mods)
TABDF<-as.data.frame(TAB)
colnames(TABDF)[1]<-"Model.long"
	TABDF$Model<-str_remove_all(TABDF$Model, "adjustment term of order")
	TABDF$Model<-str_remove_all(TABDF$Model, "adjustment terms of order")
	TABDF$Model<-str_remove_all(TABDF$Model, " NA")
TABDF$Rel.Freq<-TABDF$Freq/sum(TABDF$Freq)

##########################################################################################################
## Make some plots and examine output... 
##########################################################################################################

m<-rbind(c(1,1,2,2,7,7),
	c(1,1,2,2,7,7),
	c(3,3,4,4,7,7),
	c(3,3,4,4,7,7),
	c(5,5,5,6,6,8),
	c(5,5,5,6,6,8)
	)
layout(m)
layout.show(8)
par(mar=c(4,5,5,3), oma=c(1,1,1,1), bty="n")

#par(mfrow = c(3,3),mar=c(4,4,4,1)+0.1,oma=c(12,1,1,1))
COLS<-brewer.pal(8,'Set1')	#'Dark2'
COLS<-COLS[3:5]

plot(density(MB$P_a.hat), xlab="detection probability", col="goldenrod", main="Bootstrapped Detection Probability")
	polygon(density(MB$P_a.hat), col="goldenrod")
	abline(v=mean(MB$P_a.hat), col="red", lty=1)
	mtext(paste("mean =",round(mean(MB$P_a.hat),2)),side=3,adj=0.95 ,padj=1.5, col="black")
	mtext("A",side=3,adj=-0.1,padj=-2)

plot(density(MB$Dhat), xlab="Dhat", col="goldenrod2", main="Bootstrapped Density (YE/km2)")
	polygon(density(MB$Dhat), col="goldenrod2")
	abline(v=mean(MB$Dhat), col="red", lty=1)
		sd(MB$Dhat)/mean(MB$Dhat)
		mean(MB$Dhat.cv)
	#abline(v=quantile(MB$Dhat, c(0.1,0.9)), col="darkorange")
	mtext(paste("mean D =",round(mean(MB$Dhat),0)),side=3,adj=0.95 ,padj=1.5, col="black")
	mtext("B",side=3,adj=-0.1,padj=-2)

mean(MB$Dhat)
var(MB$Dhat)
sqrt(var(MB$Dhat))/mean(MB$Dhat)

#plot(density(MB$Prop.gt),xlab="Proportion > 270 mm",col="orange",main="Proportion ROV YE gt 270 mm")
#	polygon(density(MB$Prop.gt), col="orange")
hist(MB$Prop.gt, breaks=8, col="yellow", main="Proportion ROV > 270 mm", xlab="Proportion ROV YE > 270 mm")
	abline(v=mean(MB$Prop.gt), col="red", lty=1)
	mtext(paste("mean =",round(mean(MB$Prop.gt),3)), side=3, adj=0,padj=1.5, 
		col="black", cex=1)
	mtext("C",side=3,adj=-0.1,padj=-2)	

plot(density(MB$mean.W), xlab="mean Weight(kg)", col=COLS[1], main="Bootstrapped Mean YE Weight (kg)",
	#xlim=c(2,6)), 
	ylim=c(0,7))
	polygon(density(MB$mean.W), col=adjustcolor(COLS[1],alpha.f=0.25), border=COLS[1])
	polygon(density(MB$mean.W.gtco), col=adjustcolor(COLS[2],alpha.f=0.25), border=COLS[2])
	polygon(density(MB$CF.mean.W), col=adjustcolor(COLS[3],alpha.f=0.25), border=COLS[3])
	abline(v=mean(MB$mean.W), col=COLS[1], lty=1)
	abline(v=mean(MB$mean.W.gtco), col=COLS[2], lty=2)
	abline(v=mean(MB$CF.mean.W),col=COLS[3], lty=1)

	temp<-legend(x="topright", bty="n", 	
		legend=c("","",""), xjust=1, yjust=1)
	text(temp$rect$left + temp$rect$w, temp$text$y, 
		c(paste("ROV:",round(mean(MB$mean.W),2),"kg"),
			paste("ROVgt270:",round(mean(MB$mean.W.gtco),2),"kg"),
			paste("Port:",round(mean(MB$CF.mean.W),2),"kg")),
		col=c(COLS[1],COLS[2],COLS[3]),bty="n", pos=2, cex=1.4)

	mtext("D",side=3,adj=-0.1,padj=-2)

mean(MB$mean.W)
var(MB$mean.W)
mean(MB$mean.W.gtco)
var(MB$mean.W.gtco)
mean(MB$CF.mean.W)
var(MB$CF.mean.W)

quantile(MB$Biomass.CF,c(0.05,0.1,0.5,0.9,0.95))
quantile(MB$Biomass.ROV.all,c(0.05,0.1,0.5,0.9,0.95))
quantile(MB$Biomass.ROV.gt270,c(0.05,0.1,0.5,0.9,0.95))

plot(density(MB$Biomass.ROV.all), xlab="Biomass (mt)", col=COLS[1],
	main="YE Biomass (mt); dist. of lower 90% CI and mean", 
	ylim=c(0,0.00045), xlim=c(3000,14000))

	polygon(density(MB$Biomass.ROV.all), col=adjustcolor(COLS[1],alpha.f=0.25), border=COLS[1])
	polygon(density(MB$Biomass.ROV.gt270), col=adjustcolor(COLS[2],alpha.f=0.25), border=COLS[2])
	polygon(density(MB$Biomass.CF), col=adjustcolor(COLS[3],alpha.f=0.25), border=COLS[3])

	polygon(density(MB$Biomass.ROV.all.loCI), col=adjustcolor(COLS[1],alpha.f=0.25), border=COLS[1])
	polygon(density(MB$Biomass.ROV.gt270.loCI), col=adjustcolor(COLS[2],alpha.f=0.25), border=COLS[2])
	polygon(density(MB$Biomass.CF.loCI), col=adjustcolor(COLS[3],alpha.f=0.25), border=COLS[3])

	abline(v=mean(MB$Biomass.ROV.all), col=COLS[1])
	abline(v=mean(MB$Biomass.ROV.gt270), col=COLS[2])
	abline(v=mean(MB$Biomass.CF), col=COLS[3])

	abline(v=mean(MB$Biomass.ROV.all.loCI), col=COLS[1])
	abline(v=mean(MB$Biomass.ROV.gt270.loCI), col=COLS[2])
	abline(v=mean(MB$Biomass.CF.loCI), col=COLS[3])
	
	
	temp<-legend(x="topright", cex=1, bty="n",  # inset=0,	#pch=c(),
		legend=c("","","","","",""), xjust=1, yjust=1)
	
	text(temp$rect$left + temp$rect$w, temp$text$y, 
		c(paste("Low 90% all ROV",round(mean(MB$Biomass.ROV.all.loCI),0),"mt"),
			paste("Low 90% ROV gt270",round(mean(MB$Biomass.ROV.gt270.loCI),0),"mt"),
			paste("Low 90% CF wts",round(mean(MB$Biomass.CF.loCI),0),"mt"),
			paste("Mean all ROV",round(mean(MB$Biomass.ROV.all),0),"mt"),
			paste("Mean ROV gt270",round(mean(MB$Biomass.ROV.gt270),0),"mt"),
			paste("Mean CF wts",round(mean(MB$Biomass.CF),0),"mt")
			),
		col=c(COLS[1],COLS[2],COLS[3]),bty="n", pos=2, cex=1.4
		)


	mtext("E",side=3,adj=-0.1,padj=-2)

mean(MB$Biomass.ROV.all)
var(MB$Biomass.ROV.all)
sqrt(var(MB$Biomass.ROV.all))/mean(MB$Biomass.ROV.all)
mean(MB$Biomass.ROV.gt270)
var(MB$Biomass.ROV.gt270)
sqrt(var(MB$Biomass.ROV.gt270))/mean(MB$Biomass.ROV.gt270)
mean(MB$Biomass.CF)
var(MB$Biomass.CF)
sqrt(var(MB$Biomass.CF))/mean(MB$Biomass.CF)

hist(MB$CvM_pvalue, xlab="CvM P values", col="grey", main="GOF P-value", breaks=20)
	abline(v=mean(MB$CvM_pvalue), col="red", lty=1)
	#abline(v=getmode(MB$CvM_pvalue), col="red", lty=2)
	mtext("F",side=3,adj=-0.,padj=-1.8)
barplot(TABDF$Rel.Freq~TABDF$Model, xlab="", ylab="Proportion of replicates with lowest AIC", #xaxt="n")
	las=2, main="Model frequency in bootstrap sample", #cex=1.2,
	col=c("darkorange","darkorange","goldenrod4","goldenrod4","goldenrod4","goldenrod4","goldenrod4"
		,"goldenrod4","goldenrod4","goldenrod4","goldenrod4"))
	mtext("G",side=3,adj=-0.1,padj=-2)

######################################################################################################################
######################################################################################################################


