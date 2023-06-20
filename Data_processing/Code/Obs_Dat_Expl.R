
wd="C:/Users/pjjoy/Documents/Groundfish Biometrics/UW Bayesian Class/Project/Observer_Data/"
setwd(wd)
getwd()

####################################################################################
# Jane's notes: 1) Goundfish Total Discards by Species.csv (CAS DATA): Species specific total discards 
# and discard rate by regulatory area. 
# Source: NMFS AKRO BLEND/Catch Accounting System (CAS). THESE DATA ARE CONFIDENTIAL. 
# Query filters: species yelloweye and halibut, nmfs areas 650 and 659, gear HAL (hook and line), years 2013-present (CAS methods have changed since the observer program restructure, so the earlier data are not consistent). 
####################################################################################
OBS3<-read.csv("Groundfish Total Discards by Species.csv")
str(OBS3)

#Jane's code...
libs <- c("tidyverse", "janitor")
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {
  install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)

theme_set(theme_minimal())

str(OBS3)
# you'll need to check out these parsing errors...
ye <- OBS3 %>% #read_csv('data/Groundfish Total Discards by Species.csv') %>% 
  janitor::clean_names() %>% 
  mutate(species = ifelse(species_code == 145, 'yelloweye', 'halibut'),
         area = ifelse(nmfs_area == 650, 'SEO', 'SEI'))
str(ye)
colnames(ye)[1]<-"year"
unique(ye$year)

rates <- ye %>% 
  group_by(year, area, species) %>% 
  dplyr::summarise(total_catch_mt = sum(total_catch_mt)) %>% 
  pivot_wider(id_cols = c(area, year), names_from = species, values_from = total_catch_mt) %>% 
  mutate(ye_bycatch = yelloweye / halibut)
str(rates)

##Data for model: 
rates$yelloweye
rates$halibut
rates$area

rates %>% 
  ggplot(aes(x = factor(year), y = ye_bycatch, col = area, shape = area)) +
  geom_point() +
  geom_line(aes(group = area)) +
  labs(x = 'year', y = 'bycatch rate',
       title = 'total catch of yelloweye/halibut')


## Not sure where the fuck this data came from?  Jesus man, keep it together!!! 

OBS<-read.csv("YE_Discards_ObserverDat.csv")
str(OBS)

unique(OBS$Year)
unique(OBS$FMP.Area)
unique(OBS$FMP.Subarea)
unique(OBS$Species.Common.Name)
unique(OBS$Species.Code)

YE.OBS<-OBS[OBS$Species.Code == 145,]
nrow(YE.OBS)
head(YE.OBS)
hist(YE.OBS$Discard.Rate)
plot(YE.OBS$Discard.Rate~YE.OBS$Year)

#just look at one vessel and one year and see what the data looks like... 

#Permits<-unique(OBS$Processor.Federal.Permit)
Vessels<-unique(OBS$Vessel.Name)

EX<-OBS[OBS$Vessel.Name == Vessels[10] & OBS$Year == 2013,]
nrow(EX)
EX.YE<-EX[EX$Species.Code == 145,]
EX.HA<-EX[EX$Species.Code == 200,]
plot(EX.YE$Discard.Rate ~ EX.YE$Year)
plot(EX.YE$Total.Catch..mt. ~ EX.YE$Year, ylim=c(0,max(EX.HA$Total.Catch..mt.)))
points(EX.HA$Total.Catch..mt. ~ EX.HA$Year, col="red")

## Vessel 10 had two landings, one inside and one outside; lets look at the outside trip
EX<-EX[EX$FMP.Area == "GOA",]

## so Halibut bm: YE bm ratio would be...
EX$Total.Catch..mt.[2]/EX$Total.Catch..mt.[1] #<-this would be how we would estimate YE bycatch with just HA data
                                              #this is one data point (one vessel in one year)
                                              #probably modeled as a binomial...

#####################################################################
# Jane's Notes: norpac_haul_hook_count_report.csv (RAW OBSERVER DATA): A combination of 
# the observer Haul, Species Composition, and Haul Hook Count sources that reports observed 
# hauls and a flattened view of longline hook count data. 
# Source: NMFS AFSC FMA Observed Debriefed Haul, Species Composition, and Haul Hook Count 
# tables that combine older Domestic (pre 2008) and newer Debriefed (2008 onwards) sources. 
# Haul and Species Composition available 1986 to present. Haul Hook Count available 2009 to present.. 
# THESE DATA ARE CONFIDENTIAL. 
# Query filters: species yelloweye and halibut, nmfs areas 650 and 659, gear longliner (hook and line), 
# years 2008-present. 

OBS2<-read.csv("norpac_haul_hook_count_report.csv", skip=8, header=T)
nrow(OBS2)
str(OBS2)
head(OBS2)
unique(OBS2$Year)
unique(OBS2$FMP.Area)
unique(OBS2$FMP.Subarea)
unique(OBS2$NMFS.Area)
length(unique(OBS2$Cruise))

unique(OBS2$Species)
unique(OBS2$Species.Name)
#lets check OBS2 data for similar
EX2.YE<-OBS2[OBS2$Species == 322,]
unique(EX2.YE$Species.Name)
head(EX2.YE)

##OK, lets start by looking at one haul on one cruise in just outside waters 
Cruises<-unique(OBS2$Cruise)

EX2<-OBS2[OBS2$Cruise == Cruises[10] & OBS2$FMP.Area == "GOA",]
Days<-unique(EX2$Haul.Date)
EX2.d1<-EX2[EX2$Haul.Date == Days[1],]

# so for this example looks like 3 actual hauls based on duplicates of $Official.Total.Catch..mt.
#to get unique haul data, look at unique values of Official.total.catch
Hauls<-unique(EX2.h1$Official.Total.Catch..mt.)
EX2.d1.h1<-EX2.d1[EX2.d1$Official.Total.Catch..mt. == Hauls[1],]

#so, presumably our YE:HA ratio would be 
EX2.d1.h1$Extrapolated.Weight..kg.[2]/EX2.d1.h1$Extrapolated.Weight..kg.[1]


