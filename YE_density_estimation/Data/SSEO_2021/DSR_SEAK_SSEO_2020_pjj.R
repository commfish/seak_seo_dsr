#Estimate Line Lengths from Sub/ROV DSR surveys
#Clean up line length data and trim out bad line lengths
#install.packages("broom", type="binary")
##REFERENCE INSTALLED PACKAGES
library(tidyverse)
library(lubridate)
library(zoo)
library(Distance)

##Set working directory
setwd("D:/Groundfish Biometrics/Yelloweye/YE Code")

##IMPORT NAVIGATION, QUALITY CONTROL, AND SPECIMEN DATA#########################
sseo_nav <-read.csv("SSEO_2020_Nav.csv")
sseo_qc <-read.csv("SSEO_2020_QC.csv")
sseo_species <-read.csv("SSEO_2020_Species.csv")

head(sseo_nav)
	str(sseo_nav)
head(sseo_qc)
	str(sseo_qc)
view(sseo_qc)
	str(sseo_qc)
head(sseo_species)
#Plot ROV transects#

ggplot(sseo_nav, aes(ROV_X, ROV_Y)) + geom_point() +
  facet_wrap(~DIVE_NO, scales = "free")

##TRANSECT LINE ESTIMATION######################################################
#Need to get both tables in similar format before joining. This selects for just the dive numbers in the nav table, not text (Column DIVE_NO)
sseo_nav %>% mutate(Dive = substr(DIVE_NO, 6, 7)) -> sseo_nav

#Converts dive column from text to numeric format
sseo_nav$Dive <- as.numeric(sseo_nav$Dive)

#Rename column names (Typically from SECONDS to Seconds, but this was already "Seconds")
plyr::rename(sseo_nav, replace = c("Seconds" = "Seconds"))-> sseo_nav

view(sseo_nav)
view(sseo_qc)
#Join Tables using a "full join"remove.
transect <- full_join(sseo_nav, sseo_qc, by = c("DIVE_NO", "Seconds"))
view(transect)

#Need to fill in missing values in nav table from the quality control so the good and bad sections have
#time assignments for the entire transect.

#fill() function automatically replaces NA values with previous value

#Use select() to only keep the necessary columns i.e. Seconds, Dive #, x, y, video quality (good/bad)
transect_qc <- transect %>% fill(Family) %>% filter(!is.na(Family)) %>% 
  select(Seconds, DIVE_NO, ROV_X, ROV_Y, Family) 

#Check Data
view(transect_qc)

#Four rows were empty, remove rows by using this code:
new_transect_qc <- (na.omit(transect_qc))

#Check that rows were omitted
view(new_transect_qc)

#Output cleaned data table (make sure to use updated transect_qc table)
write.csv(new_transect_qc, file = "Output/new_transect_qc.csv")

#Use ggplot to look at data to identify good/bad areas:
new_transect_qc <- read_csv("Output/new_transect_qc.csv")

jpeg(filename = "Figures/SSEO_rov_transects2020.jpg",
     width = 12, height = 15, units = "in", res = 50)

ggplot(new_transect_qc, aes(ROV_X, ROV_Y)) + geom_point(aes(colour = factor(Family))) +
  facet_wrap(~DIVE_NO, scales = "free") +
  theme(axis.text.x = element_text(angle = 90))

dev.off()

#Check line transects in ArcGIS to ensure transects follow a straight path
#ggplot2() graphs can make transects look zigzaggy when they are actually straight due to scaling issues

##SMOOTHING LINE TRANSECT DATA##################################################

#Convert "Dive" or "DIVE_NO" column from numeric value to a factor:
new_transect_qc$DIVE_NO <-factor(new_transect_qc$DIVE_NO)
head(new_transect_qc)
str(new_transect_qc)
dim(new_transect_qc)
new_transect_qc$DIVE_NO <- factor(new_transect_qc$DIVE_NO)
levels(new_transect_qc$DIVE_NO)
is.factor(new_transect_qc$DIVE_NO)
is.numeric(new_transect_qc$DIVE_NO)

##remove bad points here?  PJJ Added this line to do so...						##PJ added this to remove Bad points here instead of in GIS
new_transect_qc<-new_transect_qc[new_transect_qc$Family == "Good",]

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
par(mfrow = c(6,6))
for (i in 1:length (levels (new_transect_qc$DIVE_NO))) {
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

write.csv(outPut, file = "Output/SSEO_smooth_transect_output.csv") 
View(outPut)

###Added by Phil 7/8 <-summarize transect lengths...#####
### Phil's R code to calculate transect length in R...
D<-unique(outPut$Dive)
Tran.Length<-data.frame()

i<-1
for (d in D){
	T<-outPut[outPut$Dive == d,]	
	Tran.Length[i,"Dive"]<-d
	Tran.Length[i,"transect_length_m"]<-sum(T$dist)
	i<-i+1
}

write.csv(Tran.Length, file = "Output/SSEO_2020_smooth_predict_lengths.csv") 
write.csv(Tran.Length, file = "SSEO_2020_smooth_predict_lengths.csv") #for ease of coding below

#Combines your original dataset with the predicted output from smoothing function

transect_pred <- cbind(new_transect_qc, predX = outPut$X, predY = outPut$Y, Dist = outPut$dist)

#Check and export data to be used in ArcGIS
View(transect_pred)

#Use this output for ArcGIS to determine length

write.csv(transect_pred, file = "Output/2020_SSEO_smooth_predict.csv") #This file will be used in ArcGIS

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

sseo_transects <- read.csv("Data/SSEO_2020_smooth_predict_lengths_GIS.csv") #This was created in ArcGIS
sseo_transects$transect_length_m<-sseo_transects$GIS
sseo_transects<-sseo_transects[-34,]
colnames(sseo_transects)[1]<-"Dive"
sseo_transects$Dive<-as.integer(sseo_transects$Dive)

#sseo_transects <- read.csv("SSEO_2020_smooth_predict_lengths.csv") #This was created in above in R
##sseo_transects <-Tran.Length ##for R transect length calculation...

#transect_summary <- sseo_transects %>% group_by(DIVE_NO) %>%  					##NOT WORKING WITH NEW FILE; MAY NEED TO ADJUST BELOW...
#
#  summarise(transect_length_m = sum(Shape_Length, na.rm = TRUE))

transect_summary <-sseo_transects										##ADDED to ease changes below...

#Verify dive transect lengths
View(transect_summary)

##DISTANCE ANALYSIS#############################################################

#Import ROV specimen data and filter for YE only

sseo_species <- read_csv("Data/SSEO_2020_species.csv") %>% filter(Species == 145)
view(sseo_species)
str(sseo_species)
unique(sseo_species$Dive)
unique(sseo_species$Transect.Number)

#For the density estimate we only want adults and subadults as these are selected for in the fishery
#filter bio data so raw data is only adults and subadults for YE
ye_adult <- sseo_species %>% filter(Stage != "JV")
view(ye_adult)

#We also want to exclude fish that were attracted to the ROV.
ye_adult_NoAttracted <- ye_adult %>% filter(Activity != "Attracted")
view(ye_adult_NoAttracted)

#Join specimen and transect summary table together
#Columns are renamed to avoid confusion with specimen table						#NOT NECESSARY WITH PHIL CHANGES...
#plyr::rename(transect_summary, replace = c("DIVE_NO" = "Dive")) -> transect_summary
#view(transect_summary)

#Make sure to change the Area for each surveyed area (NSEO = 442, SSEO = 1056, CSEO = 1661, EYKT = 739)
SSEO_survey <- full_join(transect_summary, ye_adult_NoAttracted, by = "Dive") %>% 
  mutate(mgt_area = "SSEO", Area = 1056, distance = abs(`Mid.X..mm.` * 0.001))

SSEO_survey$Fish.L<-SSEO_survey$Length..mm.								#PJJ added to rename fish length to make easier below...

view(SSEO_survey)

##PREPARE DATA FOR DISTANCE ANALYSIS############################################

#If you have transects with zero fish observed you need to replace "NAs" with zero for a given transect
# PJJ: Added in Fish.L, Depth, Stage as covariates... 

SSEO_distance <- SSEO_survey %>% select(Year, mgt_area, Area, Dive, transect_length_m,
						    Fish.L,Depth,Stage, distance) %>%
  mutate(YEAR = replace_na(Year, 2020)) 

plyr::rename(SSEO_distance, replace = c("mgt_area" = "Region.Label", "Dive" = "Sample.Label",
                                        "transect_length_m" = "Effort" )) -> SSEO_distance

#Data has to be in a data frame in order to work in distance
as.data.frame(SSEO_distance) -> SSEO_distance


view(SSEO_distance)

write.csv(SSEO_distance, file = "Data/SSEO_distance_data_GIStran_for_analysis.csv")
################################################################################
##2020 SSEO DENSITY ANALYSIS####################################################
################################################################################
#EXPLORE THE DATA..
################################
#View Summary of Data
summary(SSEO_distance$distance)

#NEED TO FIND OUT LENGTH DIST OF COMM HARVEST AND EXCLUDE FISH <MINIMUM HARVEST SIZE

#View Historgram of perpendicular distance from transect line 
hist(SSEO_distance$distance, xlab = "Distance (m)")
hist(SSEO_distance$Fish.L, xlab = "Fish.L (mm)")
hist(SSEO_distance$Depth, xlab = "Depth (m)")

DAT<-SSEO_distance

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

##check out Fish Length
ggplot(DAT, aes(x = Fish.L, y = distance)) +
  geom_point(alpha = 0.25, size = 1.6) + labs(x = "YE Length (mm)",
                                              y = "Distance (mm)")
plot(DAT$distance~DAT$Fish.L)
abline(lm(DAT$distance~DAT$Fish.L))
summary(lm(DAT$distance~DAT$Fish.L))

#######################################################################
## Run Some models...
##
###*Model 1 - Distance Model Fitting*###########################################
SSEO.model1 <- ds(SSEO_distance, key = "hn", adjustment = NULL,
                  convert.units = 0.000001)

summary(SSEO.model1$ddf)

plot(SSEO.model1, nc = 10)



###*Model 2 - Cosine Adjustment*################################################
SSEO.model2 <- ds(SSEO_distance, key = "hn", adjustment = "cos",
                  convert.units = 0.000001)

summary(SSEO.model2)

plot(SSEO.model2)

###*Model 3 - Cosine Adjustment with hazard rate key function*##################
SSEO.model3 <- ds(SSEO_distance, key = "hr", adjustment = "cos",
                  convert.units = 0.000001)

summary(SSEO.model3)

plot(SSEO.model3)

###*Model 4 - Hazard key function with Hermite polynomial adjustment*###########
SSEO.model4 <-ds(SSEO_distance, key = "hr", adjustment = "herm",
                 convert.units = 0.000001)

summary(SSEO.model4)

plot(SSEO.model4)


###Goodness of Fit Test#########################################################
gof_ds(SSEO.model1)

gof_ds(SSEO.model2)

gof_ds(SSEO.model3)

gof_ds(SSEO.model4)

str(SSEO.model1$dht$individuals, max = 1)

SSEO.model1$dht$individuals$summary

SSEO.model2$dht$individuals$summary

SSEO.model3$sht$individuals$summary

SSEO.model4$sht$individuals$summary

###Abundance Estimate###########################################################
SSEO.model1$dht$individuals$N

###Density Estimate#############################################################
SSEO.model1$dht$individuals$D
