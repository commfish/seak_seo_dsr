##############################################################################
## Bootstrap and model averaging for YE ROV surveys
## This script is used to model average difference models using the recommended procedures 
## where transects are the samples that are bootstrapped and within each replicate
## the candidate models are fit to the data and the top ranked model from each iteration is used
## to produce density estimates and other outputs from the distance model
## LAST UPDATED 8/27/24 - LSC
#########################################################################
{library(tictoc)
library(Distance)
library(boot)
library(tidyverse)}

YEAR<-2023
Subd<-"EYKT"
#surveyed area (NSEO = 442, SSEO = 1056, CSEO = 1661, EYKT = 739)
surveyed_area<-739

## Load the processed ROV data produced in the "DSR_ROV_SEAK_ROV_Processing.R"
#DAT<-read.csv(paste0("YE_density_estimation/Data/",Subd,"_",YEAR,"/",Subd,"_",YEAR,"_distance_data_rtran_for_analysis.csv"))
DAT<-read.csv(paste0("YE_density_estimation/Data/",Subd,"_",YEAR,"/",Subd,"_",YEAR,"_distance_data_GIStran_for_analysis.csv"))
head(DAT)
DAT <- DAT[-c(1, 2, 3, 142), ]
### get model list created in "DSR_ROV_SEAK_base_distance_modelling.R"
G<-read.csv(paste0("YE_density_estimation/Data/",Subd,"_",YEAR,"/",Subd,"_",YEAR,"_Bootstrap_Model_List_all.csv"))
G<-read.csv(paste0("YE_density_estimation/Data/",Subd,"_",YEAR,"/",Subd,"_",YEAR,"_Bootstrap_Model_List_Trunc.csv"))

## Pick the models with deltaAIC less than four.  These are the candidate models to average
## in the bootstrap
dAIC_lt4<-G[G$Delta.AIC <= 4,]
nrow(dAIC_lt4)

## Establish conversion units for going from meters to square kilometers... 
CU<-convert_units(distance_units="Meter", effort_units="Meter", 
	area_units="Square kilometer")

## Are models truncated?
TRUNC<-"Y"		##"Y"

#################################################################
## Run the bootstrap...

nboot<-1000					#number of replications
smpls<-unique(DAT$Sample.Label)	#unique transects...
head(DAT)

Model.Boot<-data.frame(NA)		#results data frame 

tic("boot")					#time the bootstrap for reference
pb = txtProgressBar(min = 0, max = length(nboot), initial = 0, style=3) #progress bar to monitor... 

for (r in 1:nboot){										r<-1
### randomly sample the transects with replacement...
	for (s in 1:1){ #length(smpls)
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
	for (i in 1:nrow(dAIC_lt4)){								i<-1
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
		
## Depending on what covariates are included, run the appropriate model
		## invisible(capture.output... is just a wrapper to restrain printing output in the loop
		## try(... allows the for loop to keep running if you get a bad model that doesn't work...
		if (TRUNC == "Y") {
			if (FORM == "~Depth") {
				try(invisible(capture.output(M <- ds(R2, key = KEY, adjustment = ADJ, transect="line", formula=~Depth,
            	      convert_units = CU, truncation=c(right="5%"), quiet=TRUE),type="message")),silent=TRUE)
			} else if (FORM == "~Stage") {
				try(invisible(capture.output(M <- ds(R2, key = KEY, adjustment = ADJ, transect="line", formula=~Stage,
            	      convert_units = CU, truncation=c(right="5%"), quiet=TRUE),type="message")),silent=TRUE)
			} else if (FORM == "~Depth + Stage") {
				try(invisible(capture.output(M <- ds(R2, key = KEY, adjustment = ADJ, transect="line", formula=~Depth + Stage,
            	      convert_units = CU, truncation=c(right="5%"), quiet=TRUE),type="message")),silent=TRUE)
			} else {
				try(invisible(capture.output(M <- ds(R2, key = KEY, adjustment = ADJ, transect="line", formula=~1,
            	      convert_units = CU, truncation=c(right="5%"), quiet=TRUE),type="message")),silent=TRUE)
			}
		} else {
			if (FORM == "~Depth") {
				try(invisible(capture.output(M <- ds(R2, key = KEY, adjustment = ADJ, transect="line", formula=~Depth,
            	      convert_units = CU, quiet=TRUE),type="message")),silent=TRUE)
			} else if (FORM == "~Stage") {
				try(invisible(capture_output(M <- ds(R2, key = KEY, adjustment = ADJ, transect="line", formula=~Stage,
            	      convert.units = CU, quiet=TRUE),type="message")),silent=TRUE)
			} else if (FORM == "~Depth + Stage") {
				try(invisible(capture.output(M <- ds(R2, key = KEY, adjustment = ADJ, transect="line", formula=~Depth + Stage,
            	      convert_units = CU, quiet=TRUE),type="message")),silent=TRUE)
			} else {
				try(invisible(capture.output(M <- ds(R2, key = KEY, adjustment = ADJ, transect="line", formula=~1,
            	      convert_units = CU, quiet=TRUE),type="message")),silent=TRUE)
			}
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

##extract top model from list of candidates... 
	try(RANK<-summarize_ds_models(out.list.2[[1]]))	#this is usually used to rank a series of model but I used it here 
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
		RANK$Nhat<-sapply(out.list.2[1], function(x) summary(x)$dht$individuals$N$Estimate)
		RANK$Nhat.se<-sapply(out.list.2[1], function(x) summary(x)$dht$individuals$N$se)
		RANK$Nhat.cv<-sapply(out.list.2[1], function(x) summary(x)$dht$individuals$N$cv)
		RANK$Nhat.lcl<-sapply(out.list.2[1], function(x) summary(x)$dht$individuals$N$lcl)
		RANK$Nhat.ucl<-sapply(out.list.2[1], function(x) summary(x)$dht$individuals$N$ucl)

		RANK$Dhat<-sapply(out.list.2[1], function(x) summary(x)$dht$individuals$D$Estimate)
		RANK$Dhat.se<-sapply(out.list.2[1], function(x) summary(x)$dht$individuals$D$se)
		RANK$Dhat.cv<-sapply(out.list.2[1], function(x) summary(x)$dht$individuals$D$cv)
		RANK$Dhat.lcl<-sapply(out.list.2[1], function(x) summary(x)$dht$individuals$D$lcl)
		RANK$Dhat.ucl<-sapply(out.list.2[1], function(x) summary(x)$dht$individuals$D$ucl)

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
##see how long it took... 
toc()	#divide by 3600 to get hours...

## rename output for saving and analysis...
Model.Boot$Unique.Mods<-paste(Model.Boot$Key_function, Model.Boot$Formula)
MB<-Model.Boot
write.csv(MB, file = paste0("YE_density_estimation/Output/",Subd,"_",YEAR,"Trunc_Model_Avg_boot_1000.csv"))

######################################################################################
## Examine bootstrap
######################################################################################
nrow(MB)
head(MB)
mean(MB$Dhat)
sqrt(var(MB$Dhat))/mean(MB$Dhat)

#mode calc function if needed
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

#GOF tests
min(MB$CvM_pvalue)

#could cull some of worse fits? 
nrow(MB)
nrow(MB[MB$CvM_pvalue > 0.1,])
MB<-MB[MB$CvM_pvalue > 0.1,]


## Model frequency (weight)
TAB<-table(MB$Unique.Mods)
TABDF<-as.data.frame(TAB)
colnames(TABDF)[1]<-"Model.long"
	TABDF$Model<-str_remove_all(TABDF$Model, "adjustment term of order")
	TABDF$Model<-str_remove_all(TABDF$Model, "adjustment terms of order")
	TABDF$Model<-str_remove_all(TABDF$Model, " NA")
TABDF$Rel.Freq<-TABDF$Freq/sum(TABDF$Freq)

##
par(mfrow = c(3,1),mar=c(4,4,4,1)+0.1,oma=c(15,1,1,1))
hist(MB$P_a.hat, xlab="detection probability", breaks=20, col="goldenrod", main="Bootstrapped Detection Probability")
	abline(v=mean(MB$P_a.hat), col="red", lty=1)
	abline(v=getmode(MB$P_a.hat), col="red", lty=2)
hist(MB$Nhat, xlab="SSEO Nhat", breaks=20, col="goldenrod2", main="Bootstrapped SSEO Abundance")
	abline(v=mean(MB$Nhat), col="red", lty=1)
	abline(v=getmode(MB$Nhat), col="red", lty=2)
barplot(TABDF$Rel.Freq~TABDF$Model, xlab="", ylab="Proportion of replicates with lowest AIC", #xaxt="n")
	las=2, main="Model frequency in bootstrap sample",
	col=c("darkorange","darkorange","darkorange","red","goldenrod4","goldenrod4","goldenrod4"
		,"goldenrod4","goldenrod4","goldenrod4","goldenrod4"))
##or
par(mfrow = c(4,1),mar=c(4,4,4,1)+0.1,oma=c(12,1,1,1))
plot(density(MB$P_a.hat), xlab="detection probability", col="goldenrod", main="Bootstrapped Detection Probability")
	polygon(density(MB$P_a.hat), col="goldenrod")
	abline(v=mean(MB$P_a.hat), col="red", lty=1)
	abline(v=median(MB$P_a.hat), col="red", lty=2)
	mtext("mean (solid)",side=3,adj=0.95 ,padj=1.5, col="red")
	mtext("mode (dashed)",side=3,adj=0.95 ,padj=3, col="red")
plot(density(MB$Dhat), xlab="Dhat", col="goldenrod2", main="Bootstrapped SSEO Yelloweye Density")
	polygon(density(MB$Dhat), col="goldenrod2")
	abline(v=mean(MB$Dhat), col="red", lty=1)
	abline(v=median(MB$Dhat), col="red", lty=2)
	abline(v=quantile(MB$Dhat, c(0.1,0.9)), col="darkorange")
	mtext("90% CI",side=3,adj=0.95 ,padj=1.5, col="darkorange")
hist(MB$CvM_pvalue, xlab="CvM P values", col="goldenrod3", main="GOF P-value (> 0.1 means good!)", breaks=20)
	abline(v=mean(MB$CvM_pvalue), col="red", lty=1)
	#abline(v=getmode(MB$CvM_pvalue), col="red", lty=2)
barplot(TABDF$Rel.Freq~TABDF$Model, xlab="", ylab="Proportion of replicates with lowest AIC", #xaxt="n")
	las=2, main="Model frequency in bootstrap sample",
	col=c("darkorange","darkorange","darkorange","red","goldenrod4","goldenrod4","goldenrod4"
		,"goldenrod4","goldenrod4","goldenrod4","goldenrod4"))
	
##bottstrapped detection probability
hist(MB$P_a.hat)
mean(MB$P_a.hat)
getmode(MB$P_a.hat)
sd(MB$P_a.hat)
quantile(MB$P_a.hat, c(0.05,0.95))
plot(density(MB$P_a.hat))

##bootstrap N
hist(MB$Dhat)
mean(MB$Dhat)	#mean of the mean... 
getmode(MB$Nhat)
sd(MB$Dhat)
sd(MB$Dhat)/mean(MB$Dhat)

quantile(MB$Nhat, c(0.05,0.1,0.9,0.95))
plot(density(MB$Nhat))

mean(MB$Nhat.lcl)
mean(MB$Nhat.ucl)
##Bootstrap D
hist(MB$Dhat)
mean(MB$Dhat)
median(MB$Dhat)
getmode(MB$Dhat)
sd(MB$Dhat)
sd(MB$Dhat)/median(MB$Dhat)

quantile(MB$Dhat, c(0.05,0.95))


##?table
Mods<-unique(Model.Boot$Unique.Mods)
plot(table(Model.Boot$Unique.Mods))
TAB<-table(Model.Boot$Unique.Mods)
TABDF<-as.data.frame(TAB)
colnames(TABDF)[1]<-"Model"
TABDF$Rel.Freq<-TABDF$Freq/sum(TABDF$Freq)

################################
#UPDATE DENSITY DATA SHEET FOR ASSESSMENT
# YE_Density_SEO_subdistricts.csv

dens<-read_csv(paste0("Data_processing/Data/YE_Density_SEOsubdistricts.csv"))
str(dens)
max(dens$Year)
dens %>% filter(Year == YEAR)

dens$Density[dens$Subdistrict == Subd & dens$Year == YEAR] <- round(median(MB$Dhat),0)
dens$CV[dens$Subdistrict == Subd & dens$Year == YEAR] <- sd(MB$Dhat)/median(MB$Dhat)
dens$YE_abund[dens$Subdistrict == Subd & dens$Year == YEAR] <- round(median(MB$Nhat.lcl),0)

#check... do we need to insert another year?
ny<-data.frame(matrix(ncol=ncol(dens))); i<-1
subs<-unique(dens$Subdistrict)
colnames(ny)<-colnames(dens)
for (s in subs) { #s <-subs[1]
  ny[i,"Subdistrict"] <- s
  ny[i,"Year"] <- YEAR+1
  ny[i,"Density"] <- NA
  ny[i,"CV"] <- 1
  ny[i,"Area_km2"] <- unique(dens$Area_km2[dens$Subdistrict == s])
  ny[i,"YE_abund"] <- NA
  i<-i+1
}

dens$Density[dens$Subdistrict == Subd & dens$Year == YEAR] <- round(median(MB$Dhat),0)
dens$CV[dens$Subdistrict == Subd & dens$Year == YEAR] <- sd(MB$Dhat)/median(MB$Dhat)
dens$YE_abund[dens$Subdistrict == Subd & dens$Year == YEAR] <- round(median(MB$Nhat.lcl),0)

dens<-rbind(dens,ny)

write.csv(dens,file = paste0("Data_processing/Data/YE_Density_SEOsubdistricts.csv"))

# Save density details for table 14.7 in the SAFE report

dens_sums<-read_csv(paste0("Data_processing/Data/seo_dsr_density_summary_stats.csv"))
str(dens_sums)

newdat<-data.frame(matrix(ncol=ncol(dens_sums))); i<-1
colnames(newdat)<-colnames(dens_sums)
newdat$Area <- Subd
newdat$Year <- YEAR
newdat$`Number transects` <- length(unique(DAT$Sample.Label))
newdat$`Number YE`<- nrow(DAT[!is.na(DAT$Stage),])
newdat$`Meters surveyed`<-round(sum(unique(DAT$Effort)),0)
newdat$`Encounter rate (YE/m)`<- round(nrow(DAT[!is.na(DAT$Stage),])/sum(unique(DAT$Effort)),3)
newdat$Density_YE_km2 <- round(median(MB$Dhat),0)
newdat$`Lower CI`<- round(quantile(MB$Dhat, c(0.05)),0)
newdat$`Upper CI`<- round(quantile(MB$Dhat, c(0.95)),0)
newdat$CV <- sd(MB$Dhat)/median(MB$Dhat)

dens_sums <- rbind(dens_sums,newdat)

dens_sums <- dens_sums[order(dens_sums$Area,dens_sums$Year),]
view(dens_sums)

write.csv(dens_sums,file = paste0("Data_processing/Data/seo_dsr_density_summary_stats.csv"))
