##############################################################################
## Bootstrap and model averaging for YE ROV surveys
## following course guidelines and bootstrapping transects,
## this bootstrap include 
## 1) Model (detection function) uncertainty, and averages across models
## 2) Uncertainty in ROV length; i.e., mean lenth to determine biomass from Nhat
##	a) resamples lengths, with each length sampled from the mean and SE precision
## 3) Uncertainty in proportion of ROV fish less than 270
## 4) Uncertainty in LW relationship... bootstrapped 1k random draws from LW data from port and survey data
#########################################################################
library(tictoc)

##Set working directory
setwd("D:/Groundfish Biometrics/Yelloweye/YE Code")

##function for backtransforming lengths to weights 
w.from.l<-function (LM, length){
	syx<-summary(LM)$sigma
	cf<-exp((syx^2)/2)
	pred.log<-predict(LM,data.frame(logL=log(length)), interval = "c")
	bias.pred.orig <- exp(pred.log) 
	pred.orig<-cf*bias.pred.orig
	pred.orig[1]
}

## Main Distance data... 
DAT<-SSEO_distance
##OR
#DAT<-read.csv("Data/SSEO_distance_data_for_analysis.csv")
DAT<-read.csv("Data/SSEO_distance_data_GIStran_for_analysis.csv")

## Get ROV length data for bootstrapping average length... 
ROV<-read_csv("Data/SSEO_2020_species.csv") %>% filter(Species == 145)
ROV.adult <- ROV %>% filter(Stage != "JV")
## get rid of unmeasured fish	#str(ROV.adult[,36])
ROV.L<-as.data.frame(ROV.adult[,c(1:2,4:34)])
ROV.L<-ROV.L[complete.cases(ROV.L[,"Length..mm."]),]
ROV.L$Length<-ROV.L$Length..mm. 
ROV.L$L.prec<-ROV.L$Precision..mm. 
## get rid of imprecise measurements
ROV.L$badL<-ROV.L$Length-3*ROV.L$L.prec
ROV.good<-ROV.L[ROV.L$badL > 0,]

## get summarized port and survey data for LW
YE.LW<-read.csv("Data/YE_LW_rel.csv")
str(YE.LW)
##most recent data for select area (2021 = SSEO)
YE.latest<-YE.LW[YE.LW$Year >= 2019 & YE.LW$Groundfish.Management.Area.Code == "SSEO",]
nrow(YE.latest)		#colnames(YE.latest)
YE.latest<-YE.latest[,-(c(1:16,18:21,23:36))]

###CUTOFF LENGTH from ROV SURVEY TO BE CONSIDERED 
CUTOFF<-270

### get model list created in DSR_SEAK_SSEO_2020_Distance_Models.R
G<-read.csv("Data/Bootstrap_Model_List.csv")

##need to use above file to run models and generate a list of candidate models
## for bootstrapping; done in R file "DSR_SEAK_SSEO_2020_Distance_Models.R"
##below, "G" is your list of candidate models from that file... 

## Pick the models with deltaAIC less than four.  These are the candidate models to average
## in the bootstrap
dAIC_lt4<-G[G$Delta.AIC <= 4,]
nrow(dAIC_lt4)

nboot<-10000	#number of replications
smpls<-unique(DAT$Sample.Label)
head(DAT)

Model.Boot<-data.frame(NA)

tic("boot")
pb = txtProgressBar(min = 0, max = length(nboot), initial = 0, style=3) 
for (r in 1:nboot){										#r<-1
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
	for (i in 1:nrow(dAIC_lt4)){								#i<-10
		## for each candidate model, extract the KEY, FORMULA and ADJUSTMENT terms to run...
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
		
		## Depending on what covariates are included, run the appropriate model
		## invisible(capture.output... is just a wrapper to restrain printing output in the loop
		## try(... allows the for loop to keep running if you get a bad model that doesn't work...
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
	
	try(RANK<-summarize_ds_models(out.list.2[[1]]))	#this is ussually used to rank a series of model but I used it here 
									#for the top model as it was convenient for the code I had
	if (length(RANK)<= 1) { 				#if else to skip in RANK produces nothing because of bad model
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
		#RANK$Nhat<-sapply(out.list.2[1], function(x) summary(x)$dht$individuals$N$Estimate)
		#RANK$Nhat.se<-sapply(out.list.2[1], function(x) summary(x)$dht$individuals$N$se)
		#RANK$Nhat.cv<-sapply(out.list.2[1], function(x) summary(x)$dht$individuals$N$cv)
		#RANK$Nhat.lcl<-sapply(out.list.2[1], function(x) summary(x)$dht$individuals$N$lcl)
		#RANK$Nhat.ucl<-sapply(out.list.2[1], function(x) summary(x)$dht$individuals$N$ucl)

		RANK$Dhat<-sapply(out.list.2[1], function(x) summary(x)$dht$individuals$D$Estimate)
		RANK$Dhat.se<-sapply(out.list.2[1], function(x) summary(x)$dht$individuals$D$se)
		RANK$Dhat.cv<-sapply(out.list.2[1], function(x) summary(x)$dht$individuals$D$cv)
		RANK$Dhat.lcl<-sapply(out.list.2[1], function(x) summary(x)$dht$individuals$D$lcl)
		RANK$Dhat.ucl<-sapply(out.list.2[1], function(x) summary(x)$dht$individuals$D$ucl)

	### Resample LW data and get LM
		Ref.rand<-round(runif(n=nrow(YE.latest), min=1, max=nrow(YE.latest)),0)
		LW.rand<-YE.latest[Ref.rand,]		#hist(Ref.rand)	#str(LW.rand)	#nrow(LW.rand)
		LWreg<-lm(logW~logL,data=LW.rand)
	#get lengths by randomly sampling from mean and sd of measured ROV fish...
		ROV.goodL$randL<-rtruncnorm(n=nrow(ROV.goodL),a=0,b=Inf,mean=ROV.goodL$Length, sd=ROV.goodL$L.prec)		#head(ROV.goodL,10)
		#Next take random sample of those lengths... so now we've accounted for inaccuracy of measurements and randomness of sampling!
		ROV.rand<-sample(ROV.goodL$randL, size=length(ROV.goodL$randL),replace=TRUE)
		#get the mean length from that sample and record mean and var
		RANK$mean.l.all<-mean(ROV.rand)
		RANK$var.l.all<-var(ROV.rand)

	#get and record the proportion of lengths greater than CUTOFF value... 
		Prop.gt<-length(ROV.rand[ROV.rand>CUTOFF])/length(ROV.rand)
		RANK$Prop.gt<-Prop.gt
		Prop.var<-sqrt(Prop.gt*(1-Prop.gt)/length(ROV.rand))
	
	#get mean length of fish greater than CUTOFF value for biomass calculations... 
		ROV.gt.co<-ROV.rand[ROV.rand>CUTOFF]
		RANK$mean.l.gtco<-mean(ROV.gt.co)
		RANK$var.l.gtco<-var(ROV.gt.co)

	#now lets get 2 things : the weight of the average fish and the average or the converted weight of each fish... (very different!)
	#first lets get the weight of the average fish... 
		RANK$W.mean.L<-w.from.l(LM=LW.lm , length=mean(ROV.rand))
		RANK$W.mean.L.gtco<-w.from.l(LM=LW.lm , length=mean(ROV.gt.co))
	#now weight each ROV fish gt cutoff...
		Ws<-sapply(ROV.rand[ROV.rand > 0], function(x)w.from.l(LM=LW.lm, length=x))	#plot(Ws~ROV.rand[ROV.rand>CUTOFF])
		RANK$mean.W<-mean(Ws)					
		RANK$mean.W.var<-var(Ws)
		Ws.gtco<-sapply(ROV.gt.co, function(x)w.from.l(LM=LW.lm, length=x))	#plot(Ws~ROV.rand[ROV.rand>CUTOFF])
		RANK$mean.W.gtco<-mean(Ws.gtco)					
		RANK$mean.W.gtco.var<-var(Ws.gtco)

	### Resample latest survey data for mean length #str(YE.latest)
		CFsample<-sample(YE.latest$Weight.Kilograms, size=length(YE.latest$Weight.Kilograms ), replace=TRUE)	#head(YE.latest)
		CF.avgW<-mean(CFsample)
		RANK$CF.mean.W<-CF.avgW

	### Biomass (NSEO = 442, SSEO = 1056, CSEO = 1661, EYKT = 739)
		RANK$Biomass.ROV.a<-RANK$Dhat*HA*Prop.gt*RANK$mean.W
		RANK$Biomass.ROV.gt<-RANK$Dhat*HA*Prop.gt*RANK$mean.W.gtco
		RANK$Biomass.CFW<-RANK$Dhat*HA*Prop.gt*CF.avgW

		if(r == 1){
		Model.Boot<-RANK[1,]
		} else {
		Model.Boot<-rbind(Model.Boot,RANK[1,])
		}
	}

#erasing RANK and Mlist to avoid confusion... 
	RANK<-NA
	Mlist<-NA
	M<-NA
setTxtProgressBar(pb,r/nboot)
}
toc()	#divide by 3600 to get hours...
Model.Boot$Unique.Mods<-paste(Model.Boot$Key_function, Model.Boot$Formula)
MB<-Model.Boot
write.csv(MB, file = "Output/Model_Avg_biomass_boot_GIStran_10K.csv")
######################################################################################
## Examine bootstrap
######################################################################################
MB<-read.csv("Output/Model_Avg_biomass_boot_5000_test.csv")
nrow(MB)
head(MB)

mean(MB$Dhat)
se(MB$Dhat)
var(MB$Dhat)
sqrt(var(MB$Dhat))

#mode calc function if needed
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

#convert biomass from kg to metric tons
MB$Biomass.ROV.all<-MB$Biomass.ROV.a/1000
MB$Biomass.ROV.gt270<-MB$Biomass.ROV.gt/1000
MB$Biomass.CF<-MB$Biomass.CFW/1000

#GOF tests
min(MB$CvM_pvalue)

#could cull some of worse fits? 
MB<-MB[MB$CvM_pvalue > 0.1,]
nrow(MB)

## Model frequency (weight)
TAB<-table(MB$Unique.Mods)
TABDF<-as.data.frame(TAB)
colnames(TABDF)[1]<-"Model.long"
	TABDF$Model<-str_remove_all(TABDF$Model, "adjustment term of order")
	TABDF$Model<-str_remove_all(TABDF$Model, "adjustment terms of order")
	TABDF$Model<-str_remove_all(TABDF$Model, " NA")
TABDF$Rel.Freq<-TABDF$Freq/sum(TABDF$Freq)

##########################################################################################################
## Pretty summary plots...
##########################################################################################################

par(mfrow = c(4,2),mar=c(4,4,4,1)+0.1,oma=c(12,1,1,1))
COLS<-brewer.pal(3,'Dark2')

plot(density(MB$P_a.hat), xlab="detection probability", col="goldenrod", main="Bootstrapped Detection Probability")
	polygon(density(MB$P_a.hat), col="goldenrod")
	abline(v=mean(MB$P_a.hat), col="red", lty=1)
	mtext(paste("mean =",round(mean(MB$P_a.hat),2)),side=3,adj=0.95 ,padj=1.5, col="red")

plot(density(MB$Dhat), xlab="Dhat", col="goldenrod2", main="Bootstrapped Dhat (YE/km2)")
	polygon(density(MB$Dhat), col="goldenrod2")
	abline(v=mean(MB$Dhat), col="red", lty=1)
		sd(MB$Dhat)/mean(MB$Dhat)
		mean(MB$Dhat.cv)
	abline(v=quantile(MB$Dhat, c(0.1,0.9)), col="darkorange")
	mtext(paste("mean D =",round(mean(MB$Dhat),2)),side=3,adj=0.95 ,padj=1.5, col="darkorange")
mean(MB$Dhat)
var(MB$Dhat)
sqrt(var(MB$Dhat))/mean(MB$Dhat)

plot(density(MB$Prop.gt),xlab="Proportion > 270 mm",col="orange",main="Proportion ROV YE gt 270 mm")
	polygon(density(MB$Prop.gt), col="orange")
	abline(v=mean(MB$Prop.gt), col="yellow", lty=1)
	mtext(paste("mean ROV prop. gt 270 mm =",round(mean(MB$Prop.gt),2)), side=3, adj=0,padj=1.5, col="darkorange")	

plot(density(MB$mean.W), xlab="mean Weight(kg)", col=COLS[1], main="Bootstrapped Mean YE Weight (kg)",
	xlim=c(2,6), ylim=c(0,6))
	polygon(density(MB$mean.W), col=adjustcolor(COLS[1],alpha.f=0.25))
	polygon(density(MB$mean.W.gtco), col=adjustcolor(COLS[2],alpha.f=0.25))
	polygon(density(MB$CF.mean.W), col=adjustcolor(COLS[3],alpha.f=0.25))
	abline(v=mean(MB$mean.W), col=COLS[1], lty=1)
	abline(v=mean(MB$mean.W.gtco), col=COLS[2], lty=2)
	abline(v=mean(MB$CF.mean.W),col=COLS[3], lty=1)
	mtext(paste("mean ROV wt. ",round(mean(MB$mean.W),2),"kg, var=",round(var(MB$mean.W),3)),
		side=3,adj=0.95 ,padj=1.5, col=COLS[1])
	mtext(paste("mean ROV wt. gt270 ",round(mean(MB$mean.W.gtco),2),"kg, var=",round(var(MB$mean.W.gtco),3)),
		side=3,adj=0.95 ,padj=3, col=COLS[2])
	mtext(paste("mean CF wt. ",round(mean(MB$CF.mean.W),2),"kg, var=",round(var(MB$CF.mean.W),3)),
		side=3,adj=0.95 ,padj=4.5, col=COLS[3])
mean(MB$mean.W)
var(MB$mean.W)
mean(MB$mean.W.gtco)
var(MB$mean.W.gtco)
mean(MB$CF.mean.W)
var(MB$CF.mean.W)

plot(density(MB$Biomass.ROV.all), xlab="Biomass (mt)", col=COLS[1],
	main="Bootstrapped YE Biomass (mt)", ylim=c(0,0.00045))
	polygon(density(MB$Biomass.ROV.all), col=adjustcolor(COLS[1],alpha.f=0.25))
	polygon(density(MB$Biomass.ROV.gt270), col=adjustcolor(COLS[2],alpha.f=0.25))
	polygon(density(MB$Biomass.CF), col=adjustcolor(COLS[3],alpha.f=0.25))
	abline(v=quantile(MB$Biomass.ROV.all, c(0.1)), col=COLS[1])
	abline(v=quantile(MB$Biomass.ROV.gt270, c(0.1)), col=COLS[2])
	abline(v=quantile(MB$Biomass.CF, c(0.1)), col=COLS[3])
	mtext(paste("10% CI all ROV",round(quantile(MB$Biomass.ROV.all, c(0.1)),0),"mt"),
		side=3,adj=0.95 ,padj=1.5, col=COLS[1])
	mtext(paste("10% CI ROV gt270",round(quantile(MB$Biomass.ROV.gt270, c(0.1)),0),"mt"),
		side=3,adj=0.95 ,padj=3, col=COLS[2])
	mtext(paste("10% CI CF weights",round(quantile(MB$Biomass.CF, c(0.1)),0),"mt"),
		side=3,adj=0.95 ,padj=4.5, col=COLS[3])
mean(MB$Biomass.ROV.all)
var(MB$Biomass.ROV.all)
sqrt(var(MB$Biomass.ROV.all))/mean(MB$Biomass.ROV.all)
mean(MB$Biomass.ROV.gt270)
var(MB$Biomass.ROV.gt270)
sqrt(var(MB$Biomass.ROV.gt270))/mean(MB$Biomass.ROV.gt270)
mean(MB$Biomass.CF)
var(MB$Biomass.CF)
sqrt(var(MB$Biomass.CF))/mean(MB$Biomass.CF)

hist(MB$CvM_pvalue, xlab="CvM P values", col="yellow3", main="GOF P-value (> 0.1 means good!)", breaks=20)
	abline(v=mean(MB$CvM_pvalue), col="red", lty=1)
	#abline(v=getmode(MB$CvM_pvalue), col="red", lty=2)
barplot(TABDF$Rel.Freq~TABDF$Model, xlab="", ylab="Proportion of replicates with lowest AIC", #xaxt="n")
	las=2, main="Model frequency in bootstrap sample",
	col=c("darkorange","darkorange","goldenrod4","goldenrod4","goldenrod4","goldenrod4","goldenrod4"
		,"goldenrod4","goldenrod4","goldenrod4","goldenrod4"))

######################################################################################################################
	
##bottstrapped detection probability
hist(MB$P_a.hat)
mean(MB$P_a.hat)
getmode(MB$P_a.hat)
sd(MB$P_a.hat)
quantile(MB$P_a.hat, c(0.05,0.95))
plot(density(MB$P_a.hat))

##bootstrap N
hist(MB$Nhat)
mean(MB$Nhat)	#mean of the mean... 
getmode(MB$Nhat)
sd(MB$Nhat)
quantile(MB$Nhat, c(0.05,0.1,0.9,0.95))
plot(density(MB$Nhat))

mean(MB$Nhat.lcl)
mean(MB$Nhat.ucl)
##Bootstrap D
hist(MB$Dhat)
mean(MB$Dhat)
getmode(MB$Dhat)
sd(MB$Dhat)
quantile(MB$Dhat, c(0.05,0.95))


##?table
Mods<-unique(Model.Boot$Unique.Mods)
plot(table(Model.Boot$Unique.Mods))
TAB<-table(Model.Boot$Unique.Mods)
TABDF<-as.data.frame(TAB)
colnames(TABDF)[1]<-"Model"
TABDF$Rel.Freq<-TABDF$Freq/sum(TABDF$Freq)

################################
