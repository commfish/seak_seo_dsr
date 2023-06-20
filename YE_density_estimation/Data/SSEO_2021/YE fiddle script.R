hist(sseo_nav$DIVE_NO)
unique(sseo_nav$DIVE_NO)
?substr
?mutate

unique(sseo_nav$Dive)
?full_join

str(new_transect_qc)
tst<-new_transect_qc[new_transect_qc$DIVE_NO == 2,]

ggplot(tst, aes(ROV_X, ROV_Y)) + geom_point(aes(colour = factor(Family))) +
  facet_wrap(~DIVE_NO, scales = "free") +
  theme(axis.text.x = element_text(angle = 90))

Xr<-max(tst$ROV_X)-min(tst$ROV_X)
Yr<-max(tst$ROV_Y)-min(tst$ROV_Y)

div<-min(Xr,Yr)
mx<-max(Xr,Yr)

tst$y.adj<-tst$ROV_Y/div
min(tst$y.adj)
max
?levels
#########################################################################
for (i in 1:length (levels (new_transect_qc$DIVE_NO))) {		#i<-1	#loops through each dive, 
  sT <- new_transect_qc [new_transect_qc$DIVE_NO == levels (new_transect_qc$DIVE_NO)[i],]	#view(sT)
  
  tX <- smooth.spline(sT$Seconds, sT$ROV_X, spar = 0.7)	#?smooth.spline
  tY <- smooth.spline(sT$Seconds, sT$ROV_Y, spar = 0.7)
  pX <- predict (tX) #gives the predicted values for lat by seconds for each observation at the level 
  pY <- predict (tY) #gives the predicted values for long by seconds for each observation at the level 
  
  prSp <- data.frame (pX, Y = pY$y) # creates the data frame with the pX values of seconds=x and y=lat and the pY values of long=Y
  names (prSp) <- c("Sec", "X", "Y") #renames seconds=Sec, y=X, and Y=Y
  #view(prSp)
  
  #Calculates difference between lat and long points		#A LITTLE CONFUSED ON THIS; ARE THESE D FOR MODEL OR DIFFERENCE BETWEEN TRANSECT POINTS???
  #Lat and longs are in the UTM coordinate system (Universal Transverse Mercator coordinate system)
  lon.diff<-diff(prSp$Y) 	#?diff
  lat.diff<-diff(prSp$X) 
  dist=sqrt(lon.diff^2 + lat.diff^2)
  dist <- c(dist,0)                 #double check this code #NOT SURE WHY WE ADDED A 0 AT THE END...
  prSp$dist <- dist			#sum of these should give transect length? 
  prSp$Dive <- i				#added so that dive number is also recorded
  
  plot (sT$ROV_X, sT$ROV_Y, main = levels (new_transect_qc$DIVE_NO)[i],asp = 1, ylab = "Y", xlab = "X") #plots the observed unsmoothed points for lat (x) vs.long (y)
  lines (prSp$X, prSp$Y, lwd = 2, col = 5)  #draws the predicted line for lat vs. long (2=red and 5=blue) 
  
  #Output
  if (i == 1){
    outPut <- prSp
  } else {
    outPut <- rbind (outPut, prSp) 
  }
}

view(outPut)
D<-unique(outPut$Dive)
Tran.Length<-matrix(NA, ncol=2)
	dimnames(Tran.Length)<-list(NA,c("Dive","Length"))
	Tran.Length<-as.data.frame(Tran.Length)
i<-1
for (d in D){
	T<-outPut[outPut$Dive == d,]	
	Tran.Length[i,"Dive"]<-d
	Tran.Length[i,"Length"]<-sum(T$dist)
	i<-i+1
}
##################################################################################
plyr::rename(transect_summary, replace = c("DIVE_NO" = "Dive")) -> transect_summary
view(transect_summary)

plyr::rename(sseo_transects, replace = c("DIVE_NO" = "Dive")) -> transect_summary	#not necessary because my file has it listed as Dive...
view(transect_summary)

SSEO_survey <- full_join(sseo_transects, ye_adult_NoAttracted, by = "Dive") %>% 
  mutate(mgt_area = "SSEO", Area = 1056, distance = abs(`Mid.X..mm.` * 0.001))

view(SSEO_survey)

##replace "transect_length_m" with "Length" to match my  new code... 

SSEO_distance <- SSEO_survey %>% select(Year, mgt_area, Area, Dive, transect_length_m, distance) %>%
  mutate(YEAR = replace_na(Year, 2020)) 

plyr::rename(SSEO_distance, replace = c("mgt_area" = "Region.Label", "Dive" = "Sample.Label",
                                        "transect_length_m" = "Effort" )) -> SSEO_distance
view(SSEO_distance)
#Data has to be in a data frame in order to work in distance
as.data.frame(SSEO_distance) -> SSEO_distance

##################################################################################
unique(SSEO_survey$Number)
#################################################################################
summary(TS)
library(Distance)
tidy(TS$dht$individuals$N)
library(dsm)
dsm.cor(TS)

get(paste0((TS)))
get(TS,
AIC(TS)[2]
N(TS)
print.dht.result(TS$dht,report="both",groups=FALSE)
print.summary.dsmodel(TS)

summary(Mlist[[12]])
summary(Mlist[[13]])
dAIC_lt4[13,]

####################################################################
## extract min AIC from list of models
###**********##############
## This is it.. don't loose!
list.c <- sapply(Mlist, function(x) AIC(x)[2]> 0)
out.list  <- Mlist[list.c]

mAIC<-sapply(out.list,function(x) AIC(x)[2])
AIC.vect<-unlist(mAIC)
min(AIC.vect)

list.c.2 <- sapply(out.list, function(x) AIC(x)[2]==min(AIC.vect))
out.list.2  <- out.list[list.c.2]
summary(out.list.2[[1]])

RANK<-summarize_ds_models(out.list.2[[1]])
#####condense code

AIC.vect<-unlist(sapply(Mlist[sapply(Mlist, function(x) AIC(x)[2]> 0)],function(x) AIC(x)[2]))

list.condition.2 <- sapply(Mlist[sapply(Mlist, function(x) AIC(x)[2]> 0)], function(x) AIC(x)[2]==min(AIC.vect))
lowAIC  <- Mlist[sapply(Mlist, function(x) AIC(x)[2]> 0)][sapply(Mlist[sapply(Mlist, function(x) AIC(x)[2]> 0)], function(x) AIC(x)[2]==min(AIC.vect))]
summary(lowAIC[[1]])

RANK<-summarize_ds_models(lowAIC[[1]])

###############################################################################################
MT<-invisible(capture.output(M <- ds(R2, key = KEY, adjustment = ADJ, transect="line", formula=FORM,
                  convert.units = CU, truncation=c(right="5%"), quiet=TRUE),type="message"))
MT<-invisible(MT <- ds(R2, key = KEY, adjustment = ADJ, transect="line", formula=FORM,
                  convert.units = CU, truncation=c(right="5%"), quiet=TRUE))
invisible(capture.output(MG <- ds(R2, key = KEY, adjustment = ADJ, transect="line", formula=FORM,
                  convert.units = CU, truncation=c(right="5%"), quiet=TRUE),type="message"))
str(MT)
summary(MG)$dht$individuals$average.p

#######################################################################
## Model Average all models within delta AIC of 4 using bootstrap methods... (William and Thomas)
## THIS MODEL BOOTSTRAPPED WRONG! RESAMPLED DISTANCES AND NOT TRANSECTS! NEED TO DO TRANSECTS (MAYBE BOTH
###########################
#1) resample distances with replacement 1000 times
#2) For each resampled data set, refit the candidate models
#3) Use the model with the lowest AIC to...
#	a) calculate average p
#	b) N in covered region
#4) for 1000 resamples, calculate
#	a) the proportion of repl sets for each model
#	b) the average "average p" selected by each of the models
#5) Use a weighted average of those p's using the proportion of times that model was selected as the weight

dAIC_lt4<-G[G$Delta.AIC <= 4,]
nrow(dAIC_lt4)
##got rid of model 9 because if continually barfs in the bootstrap...
dAIC_lt4<-dAIC_lt4[-9,]


nboot<-5	#number of replications
Strat<-unique(DAT$Sample.Label)
head(DAT)

Model.Boot<-data.frame(matrix(NA, ncol=21))

tic("boot")
pb = txtProgressBar(min = 0, max = length(nboot), initial = 0, style=3) 
for (p in 1:nboot){										#p<-1
	for (s in 1:length(Strat)){
		rn<-sample(unique(DAT$Sample.Label), size=1,replace=TRUE)
		G1<-DAT[DAT$Sample.Label == rn,]
		G1$Orig.Sample.Label<-G1$Sample.Label
		G1$Sample.Label<-s
		if (s == 1)  {
		G2<-G1
		} else {
		G2<-rbind(G2,G1)
		}
	}	#one replicated data set complete...
	G2$distance<-G2$Rep.dist	#replace distance measurements with replicated data

	##Fit models... 
	Hlist<-list()
	for (i in 1:nrow(dAIC_lt4)){								#i<-9
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
		
		if (FORM == "~Depth") {
			invisible(capture.output(H <- ds(G2, key = KEY, adjustment = ADJ, transect="line", formula=~Depth,
                  convert.units = CU, truncation=c(right="5%"), quiet=TRUE),type="message"))
		} else if (FORM == "~Stage") {
			invisible(capture.output(H <- ds(G2, key = KEY, adjustment = ADJ, transect="line", formula=~Stage,
                  convert.units = CU, truncation=c(right="5%"), quiet=TRUE),type="message"))
		} else if (FORM == "~Depth + Stage") {
			invisible(capture.output(H <- ds(G2, key = KEY, adjustment = ADJ, transect="line", formula=~Depth + Stage,
                  convert.units = CU, truncation=c(right="5%"), quiet=TRUE),type="message"))
		} else {
			invisible(capture.output(H <- ds(G2, key = KEY, adjustment = ADJ, transect="line", formula=FORM,
                  convert.units = CU, truncation=c(right="5%"), quiet=TRUE),type="message"))
		}
		Hlist[[i]]<-H		#summary(M)$dht$individuals$average.p
		assign(paste0("H",i),H)
	}
	## select top model and note in output...; MODIFY DEPENDING ON HOW MANY MODELS IN dAIC_lt4!!!!!!!
	
	list.d <- sapply(Hlist, function(x) AIC(x)[2]> 0)
	out.listH  <- Hlist[list.d]

	HmAIC<-sapply(out.listH,function(x) AIC(x)[2])
	AIC.vectH<-unlist(HmAIC)
	#min(AIC.vect)

	list.c.2H <- sapply(out.listH, function(x) AIC(x)[2]==min(AIC.vectH))
	out.list.2H  <- out.list[list.c.2H]
	#summary(out.list.2[[1]])

	RANKH<-summarize_ds_models(out.list.2H[[1]])
	length(RANKH)
	names(RANK)[names(RANK) == "Model"]<-"Full.Model.Code"
	RANK$Model<-substr(gsub('[\\{}]','',RANK$Full.Model.Code),7,50)
	names(RANK)[names(RANK) == "C-vM p-value"]<-"CvM_pvalue"
	names(RANK)[names(RANK) == "$\\hat{P_a}$"]<-"P_a.hat"
	names(RANK)[names(RANK) == "se($\\hat{P_a}$)"]<-"se.P_a.hat"
	names(RANK)[names(RANK) == "$\\Delta$AIC"]<-"Delta.AIC"
	names(RANK)[names(RANK) == "Key function"]<-"Key_function"

	RANK$CoveredArea<-sapply(out.list.2[1], function(x) summary(x)$dht$individuals$summary$CoveredArea)
	RANK$Adjustment<-sapply(out.list.2[1], function(x) summary(x)$ddf$ds$aux$ddfobj$adjustment$series)
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

	## record average p, N incovered area, D and N stats... 
	if(r == 1){
	Model.Boot<-RANK[1,]
	} else {
	Model.Boot<-rbind(Model.Boot,RANK[1,])
	}
	
	#erasing RANK and Mlist to avoid confusion... 
	RANK<-NA
	Mlist<-NA
	M<-NA
setTxtProgressBar(pb,r/nboot)
}
toc()	#divide by 3600 to get hours...
Model.Boot$Unique.Mods<-paste(Model.Boot$Key_function, Model.Boot$Formula)

write.csv(Model.Boot, file = "Output/Model_Avg_boot_1k.csv")
Test.Boot<-Model.Boot

try(invisible(capture.output(M <- ds(R2, key = KEY, adjustment = ADJ, transect="line", formula=FORM,
                  convert.units = CU, truncation=c(right="5%"), quiet=TRUE),type="message")),silent=TRUE)
summary(M)

for (i in 1:nrow(dAIC_lt4)){								#i<-10
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
		Mlist[[i]]<-M		#summary(M)$dht$individuals$average.p
		assign(paste0("M",i),M)
	}
summary(Mlist[[11]])


