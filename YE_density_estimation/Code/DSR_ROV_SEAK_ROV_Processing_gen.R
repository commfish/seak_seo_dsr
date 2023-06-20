#################################################################################################
## Modified from DSR_SEAK_SSEO_2020.R
##
## This script is used to clean and organize the raw ROV data
## minor modifications from DSR_SEAK_SSEO_2020.R script
## included in this project for completion
## minimal modifications made for consistency and ease of use in this project... 
## Aug 2021
## Phil Joy
#################################################################################################

##REFERENCE INSTALLED PACKAGES
{library(tidyverse)
library(lubridate)
library(zoo)
library(Distance)
library(broom)
library(ggplot2)
library(arsenal)
library(chron)}

#Estimate Line Lengths from Sub/ROV DSR surveys
#Clean up line length data and trim out bad line lengths
#install.packages("broom", type="binary")

##Set working directory
#setwd("D:/Groundfish Biometrics/Yelloweye/NSEO_2018/")

##IMPORT NAVIGATION, QUALITY CONTROL, AND SPECIMEN DATA#########################
nav <-read.csv("Data/NSEO_2022/2022_NSEO_NAV_Final.csv")		#("SSEO_2020_Nav.csv")
qc <-read.csv("Data/NSEO_2022/QC_NSEO_2022_summary.csv")		#("SSEO_2020_QC.csv")
species <-read.csv("Data/NSEO_2022/SPECIES_NSEO_2022_summary.csv")		#("SSEO_2020_Species.csv")

head(nav)
	str(nav)
head(qc)
	str(qc)
view(qc)
	str(qc)
head(species)
#Plot ROV transects#

idseps<-max(stringr::str_count(nav$Event_ID, "_"))+1
nav %>% separate(Event_ID, into = paste0("ide",1:idseps),sep = "_") %>%
  mutate(DIVE_NO = ide4) %>%
  select(-ide1,ide2,ide3,ide4) -> nav

#ggplot(nav, aes(ROV_X, ROV_Y)) + geom_point() +
ggplot(nav, aes(x, y)) + geom_point() +
  facet_wrap(~DIVE_NO, scales = "free")

##TRANSECT LINE ESTIMATION######################################################
#Need to get both tables in similar format before joining. This selects for just the dive numbers in the nav table, not text (Column DIVE_NO)
## If GIS processing has been done, skip down to  LOAD PROCESSED GIS DATA
################################################################################
unique(nav$DIVE_NO)
d17<-nav[nav$DIVE_NO == "17",]
d17b<-nav[nav$DIVE_NO == "17b",]
summary(comparedf(d17,d17b))  #17 and 17b look different?

#nav %>% mutate(Dive = substr(DIVE_NO, 6, 7)) -> subd_nav
nav %>% mutate(Dive = as.numeric(DIVE_NO),
               ROV_X = x, ROV_Y = y) -> subd_nav
unique(subd_nav$Dive)
head(subd_nav)
#Converts dive column from text to numeric format

#Rename column names (Typically from SECONDS to Seconds, but this was already "Seconds")
plyr::rename(subd_nav, replace = c("Sec" = "Seconds"))-> subd_nav

view(subd_nav)
head(qc)

idseps<-max(stringr::str_count(qc$Filename, "_"))+1
qc %>% separate(Filename, into = paste0("ide",1:idseps),sep = "_") %>%
  mutate(DIVE_NO = ide5,
         Dive = as.numeric(DIVE_NO),
         Seconds = hours(times(Time..HMS.))*3600+minutes(times(Time..HMS.))*60+seconds(times(Time..HMS.))) %>%
  select(-c(ide1,ide2,ide3,ide4, ide5, ide6)) -> qc

unique(qc$DIVE_NO)

with(qc, table(Dive, Depth))
#qc$Dive<-as.numeric(qc$DIVE_NO)

#Join Tables using a "full join"remove.
transect <- full_join(subd_nav, qc, by = c("Dive", "Seconds"))	#SSEO2020#c("DIVE_NO", "Seconds"))
unique(transect$Dive)
view(transect %>% filter(Dive == 35))
view(transect %>% filter(Dive == 1))
#Need to fill in missing values in nav table from the quality control so the good and bad sections have
#time assignments for the entire transect.

#fill() function automatically replaces NA values with previous value
head(transect)
transect$Depth[transect$DIVE_NO.x =="1" & !is.na(transect$Depth)]

#Use select() to only keep the necessary columns i.e. Seconds, Dive #, x, y, video quality (good/bad)
transect_qc <- transect %>% fill(Family) %>% filter(!is.na(Family)) %>% 
  select(Seconds, Dive, ROV_X, ROV_Y, Family, DIVE_NOraw=DIVE_NO.x, Depth, ALDT=ALDT) %>%
  group_by(Dive) %>%
  #mutate(Depth = replace_na(Depth,mean(Depth,na.rm=TRUE)))
  mutate(n = n(),
         depth_nas = sum(is.na(Depth)),
         prop_depth_na = depth_nas/n,
         Depth = ifelse(prop_depth_na == 1, NA, replace_na(Depth,mean(Depth,na.rm=TRUE))))

view(transect_qc)  
  
transect_qc %>% filter(Dive == 2)  
  
  
nrow(transect_qc)
hist(transect_qc$Depth)
#Check Data
head(transect_qc); unique(transect_qc$Dive); unique(transect_qc$DIVE_NOraw)
view(transect_qc %>% filter (Dive == 35))
#Four rows were empty, remove rows by using this code:
# get rid of empty rows except for dives with no depth data....

new_transect_qc <- transect_qc %>% 
  drop_na(-c("Depth"))

nrow(transect_qc)
nrow(new_transect_qc)
unique(new_transect_qc$Dive)
unique(new_transect_qc$Depth)

#new_transect_qc <- (na.omit(transect_qc))
view(new_transect_qc)

#Check that rows were omitted
view(new_transect_qc %>% filter(Dive == 2))

#Output cleaned data table (make sure to use updated transect_qc table)
write.csv(new_transect_qc, file = "Output/nseo_new_transect_qc.csv")

#Use ggplot to look at data to identify good/bad areas:
new_transect_qc <- read_csv("Output/nseo_new_transect_qc.csv")

jpeg(filename = "Figures/NSEO_rov_transects2022.jpg",
     width = 12, height = 15, units = "in", res = 50)

ggplot(new_transect_qc, aes(ROV_X, ROV_Y)) + geom_point(aes(colour = factor(Family))) +
#  ggplot(new_transect_qc, aes(x, y)) + geom_point(aes(colour = factor(Family))) +
  facet_wrap(~Dive, scales = "free") +			#SSEO20 used DIVE_NO instead of Dive	
  theme(axis.text.x = element_text(angle = 90))

dev.off()

#Check line transects in ArcGIS to ensure transects follow a straight path
#ggplot2() graphs can make transects look zigzaggy when they are actually straight due to scaling issues

##SMOOTHING LINE TRANSECT DATA##################################################

#Convert "Dive" or "DIVE_NO" column from numeric value to a factor:
new_transect_qc$DIVE_NO <-factor(new_transect_qc$Dive)		#SSEO20 used DIVE_NO in place of Dive
head(new_transect_qc)
str(new_transect_qc)
dim(new_transect_qc)
#new_transect_qc$DIVE_NO <- factor(new_transect_qc$DIVE_NO)
#levels(new_transect_qc$DIVE_NO)
#is.factor(new_transect_qc$DIVE_NO)
#is.numeric(new_transect_qc$DIVE_NO)

##remove bad points here?  PJ added this to remove Bad points here instead of in GIS.  However, 
## d
#new_transect_qc<-new_transect_qc[new_transect_qc$Family == "Good",]
new_transect_qc$Depth[new_transect_qc$Dive =="1"]
view(new_transect_qc)

#Verify that "Dive" column was converted correctly:
glimpse(new_transect_qc)

#Save output as a pdf so it can be reviewed easily:
pdf("Figures/2020_SSEO_smoothed_transects.jpg")

#Set up graph window as a 2X2 frame
par(mfrow = c(2,2))

dev.off()

#SMOOTHING FUNCTION#############################################################
#Smoothing loop that also calculates distance between points 
head(new_transect_qc)
par(mfrow = c(3,3))

new_transect_qc <- new_transect_qc %>% filter(Dive != 35)
unique(new_transect_qc$Dive)

Cull<-FALSE

for (i in 1:length (unique (new_transect_qc$DIVE_NO))) {  #i<-2
  sT <- new_transect_qc [new_transect_qc$DIVE_NO == levels (new_transect_qc$DIVE_NO)[i],]
  
  tX <- smooth.spline(sT$Seconds, sT$ROV_X, spar = 0.7)
  tY <- smooth.spline(sT$Seconds, sT$ROV_Y, spar = 0.7)
  pX <- predict (tX) #gives the predicted values for lat by seconds for each observation at the level 
  pY <- predict (tY) #gives the predicted values for long by seconds for each observation at the level 
  
  prSp <- data.frame (pX, Y = pY$y) # creates the data frame with the pX values of seconds=x and y=lat and the pY values of long=Y
  names (prSp) <- c("Sec", "X", "Y") #renames seconds=Sec, y=X, and Y=Y
  
  #Calculates difference between lat and long points
  #Lat and longs are in the UTM coordinate system (Universal Transverse Mercator coordinate system)
  lon.diff<-diff(prSp$Y) 
  lat.diff<-diff(prSp$X) 
  dist=sqrt(lon.diff^2 + lat.diff^2)
  dist <- c(dist,0)                 #double check this code
  prSp$dist <- dist
  prSp$Dive <- i				#added so that dive number is also recorded

  #CULLING CODE: BLOCK this out the 1st time you run to identify outliers.
  if (Cull == TRUE) {
    if (max(prSp$dist) > 18*max(prSp$dist[prSp$dist < max(prSp$dist)])) {
      prSp<-prSp[prSp$dist < max(prSp$dist),]
    }
    #transect 4 starts tracking while they are moving into position it looks like
    #get rid of the moving in quick stuff I guess... 
    if (i == 4){
      prSp<-prSp[400:max(nrow(prSp)),]
    }
  }
  #END CULLLING CODE
    
  plot (sT$ROV_X, sT$ROV_Y, main = levels (new_transect_qc$DIVE_NO)[i],asp = 1, ylab = "Y", xlab = "X") #plots the observed unsmoothed points for lat (x) vs.long (y)
  lines (prSp$X, prSp$Y, lwd = 2, col = 5)  #draws the predicted line for lat vs. long (2=red and 5=blue) 
   
  #Output
  if (i == 1){
    outPut <- prSp
  } else {
    outPut <- rbind (outPut, prSp) 
  }
}

dev.off()

if (Cull == TRUE) {
  write.csv(outPut, file = "Output/CSEO22_smooth_transect_output_culled.csv")
} else {
  write.csv(outPut, file = "Output/NSEO22_smooth_transect_output_raw2.csv") 
}

##Check transect lengths...
D<-unique(outPut$Dive)
Tran.Length<-data.frame()

i<-1
for (d in D){
  T<-outPut[outPut$Dive == d,]	
  T<-unique(T)
  Tran.Length[i,"Dive"]<-d
  Tran.Length[i,"transect_length_m"]<-sum(T$dist)
  i<-i+1
}; Tran.Length

nrow(new_transect_qc)
nrow(outPut)
unique(outPut$Dive)
unique(new_transect_qc$Dive)

#for NSEO 2022 get rid of dive 35
new_transect_qc <- new_transect_qc %>% filter(Dive != 35)
##check for bad line lengths that are out of line with expectations of ~1000 m transects... 
#2022 CSEO, lots.  Lets package this bad stuff up for GIS measured length and then
# do our best to adjust line lengths with R
transect_pred <- cbind(new_transect_qc, predX = outPut$X, predY = outPut$Y, Dist = outPut$dist)

#Use this output for ArcGIS to determine length
write.csv(transect_pred, file = "Output/2022_NSEO_smooth_predict_GISreview.csv") #This file will be used in ArcGIS

#-------------------------------------------------------------------------------
# Lets fix transect length here
## Scrap code to dig around in the outPut file 
# identify bad transects in Tran.Length.  A lot of transects there is one egement that
# is totally wacky and easy to get rid of.  Another CSEO 2022 transect seemed to 
# indicate that the ROV motored up to the start of the transect... 

par(mfrow=c(1,1))
d4<-outPut[outPut$Dive == 4,]
d4raw<-new_transect_qc[new_transect_qc$Dive == 4,]
plot (d4raw$ROV_X, d4raw$ROV_Y,  ylab = "Y", xlab = "X", xlim=c(min(d4raw$ROV_X),max(d4raw$ROV_X))) #plots the observed unsmoothed points for lat (x) vs.long (y)
lines (d4$X, d4$Y, lwd = 2, col = 5) 
quantile(d4$dist,c(0.025,0.5,0.999))

d4raw<-d17raw[d17raw$ROV_X < max(d17raw$ROV_X),]
plot (d4raw$ROV_X, d4raw$ROV_Y,  ylab = "Y", xlab = "X", xlim=c(min(d4raw$ROV_X),max(d4raw$ROV_X)))

max(d4$dist)
max(d4$dist[d4$dist < max(d4$dist)])
plot(d4$dist[1:500])

head(prSp)

d5<-outPut[outPut$Dive == 5,]
sum(d5$dist)

d24<-outPut[outPut$Dive == 24,]
sum(d24$dist)

head(d4); head(d5); head(d24)
par(mfrow=c(3,1))
hist(d17$dist, breaks=100); hist(d5$dist, breaks=100); hist(d24$dist, breaks=100)

quantile(d5$dist,c(0.025,0.5,0.975))
quantile(d24$dist,c(0.025,0.5,0.998))

d24out<-d24[d24$dist > quantile(d24$dist,c(0.998)),]
d24in<-d24[d24$dist < quantile(d24$dist,c(0.998)),]
sum(d24in$dist); sum(d24out$dist)
hist(d24in$dist, breaks=100)

## NOW turn on the culling option and modify the culling section to clean it up
Cull<-TRUE

for (i in 1:length (levels (new_transect_qc$DIVE_NO))) {  #i<-6
  sT <- new_transect_qc [new_transect_qc$DIVE_NO == levels (new_transect_qc$DIVE_NO)[i],]
  
  tX <- smooth.spline(sT$Seconds, sT$ROV_X, spar = 0.7)
  tY <- smooth.spline(sT$Seconds, sT$ROV_Y, spar = 0.7)
  pX <- predict (tX) #gives the predicted values for lat by seconds for each observation at the level 
  pY <- predict (tY) #gives the predicted values for long by seconds for each observation at the level 
  
  prSp <- data.frame (pX, Y = pY$y) # creates the data frame with the pX values of seconds=x and y=lat and the pY values of long=Y
  names (prSp) <- c("Sec", "X", "Y") #renames seconds=Sec, y=X, and Y=Y
  
  #Calculates difference between lat and long points
  #Lat and longs are in the UTM coordinate system (Universal Transverse Mercator coordinate system)
  lon.diff<-diff(prSp$Y) 
  lat.diff<-diff(prSp$X) 
  dist=sqrt(lon.diff^2 + lat.diff^2)
  dist <- c(dist,0)                 #double check this code
  prSp$dist <- dist
  prSp$Dive <- i				#added so that dive number is also recorded
  prSp$avg.Depth <- mean(sT$Depth)
  
  #CULLING CODE: BLOCK this out the 1st time you run to identify outliers.
  if (Cull == TRUE) {
    if (max(prSp$dist) > 18*max(prSp$dist[prSp$dist < max(prSp$dist)])) {
      prSp<-prSp[prSp$dist < max(prSp$dist),]
    }
    #transect 4 starts tracking while they are moving into position it looks like
    #get rid of the moving in quick stuff I guess... 
    if (i == 4){
      prSp<-prSp[400:max(nrow(prSp)),]
    }
  }
  #END CULLLING CODE
  
  plot (sT$ROV_X, sT$ROV_Y, main = levels (new_transect_qc$DIVE_NO)[i],asp = 1, ylab = "Y", xlab = "X") #plots the observed unsmoothed points for lat (x) vs.long (y)
  lines (prSp$X, prSp$Y, lwd = 2, col = 5)  #draws the predicted line for lat vs. long (2=red and 5=blue) 
  
  #Output
  if (i == 1){
    outPut <- prSp
  } else {
    outPut <- rbind (outPut, prSp) 
  }
}
head(outPut)
D<-unique(outPut$Dive)
Tran.Length<-data.frame()

i<-1
for (d in D){
  T<-outPut[outPut$Dive == d,]	
  T<-unique(T)
  Tran.Length[i,"Dive"]<-d
  Tran.Length[i,"transect_length_m"]<-sum(T$dist)
  Tran.Length[i,"avg.depth"]<-mean(T$avg.Depth)
  i<-i+1
}; Tran.Length

##check for bad line lengths that are out of line with expectations of ~1000 m transects... 
#2022 CSEO, lots.  Lets package this bad stuff up for GIS measured length and then
# do our best to adjust line lengths with R

#transect_pred <- cbind(new_transect_qc, predX = outPut$X, predY = outPut$Y, Dist = outPut$dist)

#Tran.Length2<-Tran.Length

write.csv(Tran.Length, file = "Output/CSEO_2022_smooth_predict_lengths_Rfix.csv") 

### If using R measured transect lengths, turn Tran.Length into transects to proceed...
transects<-Tran.Length

#Combines your original dataset with the predicted output from smoothing function
#need to do this without the culling code I added above to get rid of outliers.



##CHECK TRANSECT LENGTHS IN ARCGIS##############################################

#Import smoothed data into ArcGIS (if needed)
#Convert the 2020_sseo_smoothed_predict.csv into a feature class in ArcGIS
#Create feature class from XY table, using Pred_X and Pred_Y and project the feature class
#as WGS84 UTM 8N (7N for EYKT). 
#Delete Bad Points
#Summarize Dist by Dive Number
#Export as .csv file


#If, for some reason this file does not output correctly to determine
#transect lengths in ArcGIS:
##Import "smooth_transect_output" .csv and "transect_qc" .csv into ArcGIS. 
##Export XY points from both files
##spatially join the "transect_qc" shapefile to the "smooth_transect_output shapefile.
##Open "Points to Line" tool and use joined shapefile for input.
##Export the attribute table as a .txt file - which will be your "smoothed_transect_lengths" file 
##that will need to be converted to a .csv file for the next step.

#This is the output from the smoothed transect made in ArcGIS (or rename Dist as Shape_Length)

#####################################################################################
#### LOAD PROCESSED GIS DATA
######################################################################################
##GIS TRANSECT LENGTHS:
r_transects<-Tran.Length
transects <- read.csv("Data/CSEO_2022/CSEO_2022_Smooth_Predict_Lengths_gis.csv")		#This was created in ArcGIS
#R CALC TRANSECT LENGTHS
#transects<-read.csv("Output/2022_CSEO_smooth_predict_lengths_Rfix.csv") 
str(transects)
nrow(r_transects); nrow(transects)

#add in depth column I coded in r transects... 
left_join(transects,r_transects,by="Dive") %>% 
  mutate(transect_length_m = transect_length_m.x) %>% 
  select(Dive,transect_length_m,avg.depth) -> transects
#transects$transect_length_m<-transects$Dist				#transects$GIS : from SSEO 2020

#transects<-transects[-34,]								#these 3 lines needed for SSEO 2020 but not NSEO 2018
#colnames(transects)[1]<-"Dive"
#transects$Dive<-as.integer(sseo_transects$Dive)

#sseo_transects <- read.csv("SSEO_2020_smooth_predict_lengths.csv") #This was created in above in R
##sseo_transects <-Tran.Length ##for R transect length calculation...

#transect_summary <- sseo_transects %>% group_by(DIVE_NO) %>%  					##NOT WORKING WITH NEW FILE; MAY NEED TO ADJUST BELOW...
#
#  summarise(transect_length_m = sum(Shape_Length, na.rm = TRUE))

transect_summary <-transects										##ADDED to ease changes below...
str(transect_summary)
#Verify dive transect lengths
View(transect_summary)

with(transect_summary, table(Dive,transect_length_m))
##DISTANCE ANALYSIS#############################################################

#Import ROV specimen data and filter for YE only

species <- read_csv("Data/CSEO_2022/SPECIES_CSEO_2022.csv")%>%filter(Species == 145)			#("Data/SSEO_2020_species.csv") %>% filter(Species == 145)
view(species)
str(species)

#quick sloppy check on number of fish/dive
Dives<-unique(species$Dive)

species %>% group_by(Dive) %>%
  summarize(n=n()) -> raw.count
View(raw.count)

#NA's are dive 18 in CSEO 2022 transect files... 

#need to fix dive numbers for CSEO 2022
colnames(species)[21]<-"DIVE_NOraw"
species %>% mutate(Dive = ifelse(DIVE_NOraw %in% c("17b"),17,
                                 ifelse(is.na(DIVE_NOraw),18,
                                        as.numeric(DIVE_NOraw)))) -> species



unique(species$Dive)
unique(species$DIVE_NOraw)

#colnames(species)[11]
#colnames(species)[11]<-"Mid.X..mm."
#colnames(species)[12]<-"Mid.Y..mm."
#colnames(species)[13]<-"Mid.Z..mm."
#colnames(species)[20]<-"Transect.Number"


#For the density estimate we only want adults and subadults as these are selected for in the fishery
#filter bio data so raw data is only adults and subadults for YE
ye_adult <- species %>% filter(Stage != "JV")
view(ye_adult)

#We also want to exclude fish that were attracted to the ROV.
ye_adult_NoAttracted <- ye_adult %>% filter(Activity != "Attracted")
view(ye_adult_NoAttracted)

#Join specimen and transect summary table together
#Columns are renamed to avoid confusion with specimen table						#NOT NECESSARY WITH PHIL CHANGES...
#plyr::rename(transect_summary, replace = c("DIVE_NO" = "Dive")) -> transect_summary
#view(transect_summary)

### Added to deal with NSEO 2018... not sure if necessary for other ones?... 
### this summarizes transect lengths for brevity; DID NOT WORK RIGHT?!?! SKIP 
Tsum<-data.frame()
T<-unique(transect_summary$Dive)
i<-1
for (t in T){
	TN<-transect_summary[transect_summary$Dive == t,]
	Tsum[i,"Dive"]<-t
	Tsum[i,"transect_length_m"]<-mean(TN$Dist)
	i<-i+1
}

#Make sure to change the Area for each surveyed area (NSEO = 442, SSEO = 1056, CSEO = 1661, EYKT = 739)
# Use Tsum or transect_summary depending on how it is set up that year
str(transect_summary)
str(ye_adult_NoAttracted)

survey <- full_join(transect_summary, ye_adult_NoAttracted, by = "Dive") %>% 	
  mutate(mgt_area = "CSEO", Area = 1661, distance = abs(`Mid.X..mm.` * 0.001))

## Find fish length column (varies year-to-year, and rename it ...
colnames(survey)[9]<-"Fish.L.mm"
#survey$Fish.L<-survey$Length..mm.								#PJJ added to rename fish length to make easier below...

str(survey)
head(survey)

##PREPARE DATA FOR DISTANCE ANALYSIS############################################

#If you have transects with zero fish observed you need to replace "NAs" with zero for a given transect
# PJJ: Added in Fish.L, Depth, Stage as covariates... 

distance <- survey %>% select(Year, mgt_area, Area, Dive, transect_length_m,
						    Fish.L.mm,Depth,Stage, distance,avg.depth) %>%
  mutate(YEAR = replace_na(Year, 2022)) %>% 
  group_by(Dive) %>% 
  mutate(Depth =replace_na(Depth,mean(avg.depth)))
  
head(distance)

plyr::rename(distance, replace = c("mgt_area" = "Region.Label", "Dive" = "Sample.Label",
                                        "transect_length_m" = "Effort" )) -> distance

#Data has to be in a data frame in order to work in distance
as.data.frame(distance) -> distance
unique(distance$Sample.Label)

view(distance)

write.csv(distance, file = "Data/CSEO_2022/CSEO_22_distance_data_GIStran_for_analysis.csv")
################################################################################
##2020 SSEO DENSITY ANALYSIS####################################################
################################################################################
#EXPLORE THE DATA..
################################
#View Summary of Data
summary(distance$distance)
head(distance)
#NEED TO FIND OUT LENGTH DIST OF COMM HARVEST AND EXCLUDE FISH <MINIMUM HARVEST SIZE

#View Historgram of perpendicular distance from transect line 
hist(distance$distance, xlab = "Distance (m)")
hist(distance$Fish.L.mm, xlab = "Fish.L (mm)")
hist(distance$Depth, xlab = "Depth (m)")

DAT<-distance

## examine COLINEARITY in covariates... Pearson R2 over 0.15 is no bueno... 
cor(x=DAT$Fish.L.mm, y=DAT$Depth, use="complete.obs", method=c("pearson"))
cor(x=DAT$Depth, y=DAT$Fish.L.mm, use="complete.obs", method=c("pearson"))^2 #
summary(lm(DAT$Fish.L.mm ~ DAT$Depth))
plot(DAT$Fish.L.mm ~ DAT$Depth)
abline(lm(DAT$Fish.L.mm ~ DAT$Depth))

summary(lm(DAT$Fish.L.mm ~ factor(DAT$Stage)))
boxplot(DAT$Fish.L.mm ~ DAT$Stage)

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

##check out Fish Length
ggplot(DAT, aes(x = Fish.L.mm, y = distance)) +
  geom_point(alpha = 0.25, size = 1.6) + labs(x = "YE Length (mm)",
                                              y = "Distance (mm)")
plot(DAT$distance~DAT$Fish.L)
abline(lm(DAT$distance~DAT$Fish.L))
summary(lm(DAT$distance~DAT$Fish.L))

#######################################################################
