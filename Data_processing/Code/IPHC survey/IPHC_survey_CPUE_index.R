# ******************************************************************************
# Title: IPHC_survey_CPUE_index ----
# Purpose: process and organize the IPHC survey data for use as an index of 
# yelloweye rockfish abundance in Southeast Outside
# Authors: Phil Joy, Rhea Ehresmann, Caitlin Stern 
# Dates modified: Oct. 2023, Oct. 2024, June 2025 
# Notes: methods include the original bootstrap methods used in the 2022 assessment
# and the new Tweedie estimator developed after the 2023 CIE review, which now includes
# all stations even if they have 0s 
# ******************************************************************************

# ******************************************************************************
# load packages, set years ----
# ******************************************************************************

library(dplyr)
library(boot)
library(ggplot2)
library(RColorBrewer)
library(sf)
library(readr)
library(vroom)
library(scales)
library(ggpubr)
library(mgcv)
library(MuMIn)
library(GGally)
library(mgcViz)
library(wesanderson)
library(RColorBrewer)
library(tweedie)

cur.yr <- 2025
pr.yr <- 2024

source(file = paste0(here::here(), "/r_helper/Port_bio_function_rke.R")) 
## this may not work - test 

# before running the code, download the following and save them in the folder
# seak_seo_dsr/Data_Processing/Data/IPHC_raw 

# 1) Set and Pacific halibut data
# Download from https://www.iphc.int/data/fiss-survey-raw-survey-data/
# Select year range 1998-present. From drop-down menu, select areas 2C and 3A.
# From species drop-down menu, select Yelloweye Rockfish. 
# From download menu, select crosstab, select CSV format, then select the sheet
# called "Set and Pacific halibut data"
# File name should be "Set and Pacific halibut data YEAR.csv" where YEAR is the previous year

# 2) Non-Pacific halibut data
# Download from https://www.iphc.int/data/fiss-survey-raw-survey-data/
# Select year range 1998-present. From drop-down menu, select areas 2C and 3A.
# From species drop-down menu, select Yelloweye Rockfish. 
# From download menu, select crosstab, select CSV format, then select the sheet
# called "Non-Pacific halibut data"
# File name should be "Non-Pacific halibut data YEAR.csv" where YEAR is the previous year

# 3) IPHC Fishery-Independent Setline Survey hook adjustment factors
# Download from https://www.iphc.int/data/fiss-survey-raw-survey-data/
# This will download as .xlsx so you need to open it and save as .csv
# File name should be "iphc-YEAR-fiss-hadj.csv" where YEAR is the previous year

# ******************************************************************************
# read in data ----
# ******************************************************************************

BUT2C3A <- read.csv(paste0(here::here(), "/Data_Processing/Data/IPHC_raw/Set and Pacific halibut data ", pr.yr, ".csv"), skipNul = T, fileEncoding = "Latin1", sep = "\t", header = TRUE) %>%
  rename(Row.number = 1) %>%
  select(-Row.number)

YE2C3A <- read.csv(paste0(here::here(), "/Data_Processing/Data/IPHC_raw/Non-Pacific halibut data ", pr.yr, ".csv"), skipNul = T, fileEncoding = "Latin1", sep = "\t", header = TRUE) %>%
  rename(Year = 1)
  
# the hook adjustment factor is used to control for hook saturation when
# calculating CPUE
hadj <- read_csv(paste0(here::here(), "/Data_Processing/Data/IPHC_raw/iphc-", pr.yr, "-fiss-hadj.csv")) %>% 
  rename(IPHC.Reg.Area = `IPHC Reg Area`, 
         IPHC.StatArea = `IPHC Stat Area`, 
         AdjHooks.Observed = `Hooks Observed`,
         Stlkey = stlkey) %>%
  filter(IPHC.Reg.Area %in% c("2C", "3A"), IPHC.StatArea <= 200)  %>%
  mutate(Year = as.integer(Year)) %>%
  select(Stlkey, Year, Station, AdjHooks.Observed, h.adj)

str(hadj)

# ******************************************************************************  
# combine data into one data frame, called Survey ----
# note: in 2024, Rhea corrected the EYKT bounds to 137-140 longitude; previously,
# those bounds were specified as 137-139, which excluded part of the area.
# ******************************************************************************  

Survey_prep <- BUT2C3A %>%
  full_join(YE2C3A, by=c("Stlkey","Station")) %>% 
  # filter for correct area, valid sets
  filter(IPHC.Stat.Area <= 200, -MidLon.fished <= 140, Eff == "Y") %>%
  # assign stat areas to management districts
  mutate(SEdist = case_when(
    IPHC.Stat.Area %in% c(182:184,171,173,174,161:163) ~ "NSEI",
    IPHC.Stat.Area %in% c(142:144,152,153) ~ "SSEI",
    IPHC.Stat.Area %in% c(142:144,152,153,182:184,171,173,174,161:163) == FALSE & -MidLon.fished >= 137 & -MidLon.fished <=140 ~ "EYKT", 
    IPHC.Stat.Area %in% c(142:144,152,153,182:184,171,173,174,161:163) == FALSE & -MidLon.fished <= 137 & MidLat.fished >= 57.5 ~ "NSEO",
    IPHC.Stat.Area %in% c(142:144,152,153,182:184,171,173,174,161:163) == FALSE & -MidLon.fished <= 137 & MidLat.fished <57.5 & MidLat.fished >= 56 ~ "CSEO",
    IPHC.Stat.Area %in% c(142:144,152,153,182:184,171,173,174,161:163) == FALSE & -MidLon.fished <= 137 & MidLat.fished < 56 ~ "SSEO",
    TRUE ~ NA
  )) %>%
  # assign management districts to Southeast Outside or Southeast Inside
  mutate(In.Out = case_when(
    SEdist %in% c("NSEO","CSEO","SSEO","EYKT") ~ "SEO",
    SEdist %in% c("NSEI","SSEI") ~ "SEI",
    TRUE ~ NA
  )) %>%
  # rename columns
  rename(YE.obs = Number.Observed, AvgDepth.fm = AvgDepth..fm., Year = Year.x) %>%
  # create depth bins in fathoms and meters
  mutate(depth_bin = cut(AvgDepth.fm, breaks = seq(0,400,50), labels = paste (seq(50,400,50))),
         depth_m = AvgDepth.fm * 1.8288,
         depth_bin_m = cut(depth_m, breaks = seq(0,800,50), labels = paste (seq(50,800,50)))) %>%
  mutate(O32.Pacific.halibut.weight..net.lb. = as.numeric(gsub(",", "", O32.Pacific.halibut.weight..net.lb.)),
         U32.Pacific.halibut.weight..net.lb. = as.numeric(gsub(",", "", U32.Pacific.halibut.weight..net.lb.)))

# join survey data set to hook adjustment factor data set
# note: in 2024, Rhea changed the code to use AdjHooks.Observed from the hadj file to fill in missing 
# hook observed fields instead of the static "140" that had been used previously, since  
# hooks observed from the hadj file is IPHC data for that Stlkey/Station

Survey <- left_join(Survey_prep, hadj, by = c("Year","Stlkey","Station")) %>% 
  mutate(YE.obs = case_when(is.na(YE.obs) == TRUE ~ 0,
                              is.na(YE.obs) == FALSE ~ YE.obs)) %>%
  mutate(HooksObserved = case_when(
    is.na(HooksObserved) == TRUE ~ AdjHooks.Observed,
    is.na(HooksObserved) == FALSE ~ HooksObserved
  )) %>%
  mutate(HooksRetrieved = as.numeric(HooksRetrieved), 
         HooksObserved = as.numeric(HooksObserved), 
         YE.exp = HooksRetrieved/HooksObserved)

str(Survey)
head(Survey, 10)
mean(Survey$AvgDepth.fm)

# what is this for?
y <- sample(unique(Survey$Year),1)
d <- sample(unique(Survey$depth_bin[Survey$Year==y]),1)
s <- sample(unique(Survey$Station[Survey$Year==y & Survey$depth_bin==d]),1)
Dat <- Survey[Survey$Year == y & Survey$depth_bin == d & Survey$Station == s,]
unique(Dat)
nrow(Dat)

# write csv of cleaned up survey data with hook adjustment factor added
write.csv(Survey, paste0(here::here(), "/Data_processing/Data/IPHC_survey_1998-", pr.yr, "_hadj.csv"))

# ****************************************************************************** 
# Look at depth, identify stations with > 0 yelloweye versus 0 yelloweye ----
# ****************************************************************************** 

Station.sum <- Survey %>% group_by(Station) %>%
  dplyr::summarise(Lat = mean(MidLat.fished, na.rm=T),
                   Lon = mean(MidLon.fished, na.rm=T),
                   Depth = mean(AvgDepth.fm, na.rm=T),
                   obs = n(),
                   mean.hal.O32.count = mean(O32.Pacific.halibut.count, na.rm=T),
                   mean.hal.O32.wt = mean(O32.Pacific.halibut.weight..net.lb., na.rm=T),
                   mean.hal.U32.count = mean(U32.Pacific.halibut.count, na.rm=T),
                   mean.hal.U32.wt = mean(U32.Pacific.halibut.weight..net.lb., na.rm=T),
                   mean.skates = mean(No..skates.set, na.rm=T),
                   mean.hooks.ob = mean(HooksObserved, na.rm=T),
                   mean.hooks.ret = mean(HooksRetrieved, na.rm=T),
                   mean.YE = mean(YE.obs, na.rm=T),
                   var.YE = var(YE.obs, na.rm=T),
                   mean.YEcpue = mean(YE.obs/HooksObserved, na.rm=T),
                   var.YEcpue = var(YE.obs/HooksObserved, na.rm=T),
                   noYE.count = sum(YE.obs == 0),
                   prop.0 = noYE.count/obs)

ye_caught_prop <- percent(1-(nrow(Station.sum[Station.sum$mean.YE == 0,])/nrow(Station.sum)))

blanks.stations<-Station.sum$Station[Station.sum$mean.YEcpue == 0]
YE.stations_non0<-Station.sum$Station[Station.sum$mean.YEcpue != 0]  #use this for boot strap method - all stations that have seen YE 
YE.stations_10p<-Station.sum$Station[Station.sum$prop.0 < 0.9]
YE.stations_20p<-Station.sum$Station[Station.sum$prop.0 < 0.8]
YE.stations_25p<-Station.sum$Station[Station.sum$prop.0 < 0.75]
YE.stations_40p<-Station.sum$Station[Station.sum$prop.0 < 0.6]
YE.allstations <- Station.sum$Station  # this is what we want to use for Tweedie model, include all stations even with 0's

length(YE.stations_non0)
length(YE.stations_10p)
length(YE.stations_20p)
length(YE.stations_25p)
length(YE.stations_40p)
length(YE.allstations)

# ****************************************************************************** 
# load the port samples that were downloaded from OceanAK ---- 
# need yelloweye weights to get WCPUE estimates
# ****************************************************************************** 

# load stat areas information
statareas <- read.csv("Data_processing/Data/g_stat_area.csv")

# define sets of stat areas for each management area
SSc<-unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "SSEO"])
CSc<-unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "CSEO"])
NSc<-unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "NSEO"])
EYc<-unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "EYKT"])

SSIc <- unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "SSEI"])
NSIc <- unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "NSEI"])

# load biological data from port sampling
Port1<-read.csv("Data_processing/Data/SEO_YE_port_sampling_bio_data_1980-1989.csv")
Port2<-read.csv("Data_processing/Data/SEO_YE_port_sampling_bio_data_1990-1999.csv")
Port3<-read.csv("Data_processing/Data/SEO_YE_port_sampling_bio_data_2000-2009.csv")
Port4<-read.csv("Data_processing/Data/SEO_YE_port_sampling_bio_data_2010-2019.csv")
Port5<-read.csv(paste0("Data_processing/Data/SEO_YE_port_sampling_bio_data_2020-", pr.yr,".csv",sep=""))

Port <- rbind(Port1,Port2,Port3, Port4, Port5) %>%
  # correct any "EYAK" codes to "EYKT"
  mutate(Groundfish.Management.Area.Code = case_when(
    Groundfish.Management.Area.Code == "EYAK" ~ "EYKT",
    .default = Groundfish.Management.Area.Code)) %>%
  # assign correct management area codes to each stat area
  mutate(GFMU = case_when(
    Groundfish.Stat.Area %in% c(SSc) | substr(Groundfish.Stat.Area.Group,1,4) %in% c("SSEO") ~ "SSEO",
    Groundfish.Stat.Area %in% c(CSc) | substr(Groundfish.Stat.Area.Group,1,4) %in% c("CSEO") ~ "CSEO",
    Groundfish.Stat.Area %in% c(NSc) | substr(Groundfish.Stat.Area.Group,1,4) %in% c("NSEO") ~ "NSEO",
    Groundfish.Stat.Area %in% c(EYc) | substr(Groundfish.Stat.Area.Group,1,4) %in% c("EYKT") ~ "EYKT",
    TRUE ~ NA
  )) %>%
  mutate(Sex = as.factor(case_when(Sex.Code == 1 ~ "Male", Sex.Code == 2 ~ "Female"))) %>%
  filter(Groundfish.Management.Area.Code != "", is.na(Weight.Kilograms) == FALSE)

str(Port)

unique(Port$Groundfish.Management.Area.Code)  #EYAK = EYKT


# Include samples that are Random (avoid the select samples for maturity/vonB that were collected)
Port.rand <- Port %>%
  filter(Sample.Type == "Random") 

# define mean weights for each year and management area  
uYEkg <- Port.rand %>% 
  group_by(Year,Groundfish.Management.Area.Code) %>%
  summarize(uYEkg = mean(Weight.Kilograms, na.rm=TRUE),
              vYEkg = var(Weight.Kilograms, na.rm=TRUE),
              NyeKG = length(Weight.Kilograms)) %>%
  as.data.frame()

# add columns to Survey data frame  
Survey$mean.YE.kg<-NA
Survey$var.YE.kg<-NA
Survey$N.YE.kg<-NA
  
  
# Get best data for average weights and add columns to survey data
# In 2024, Rhea and Laura updated this so that we get 300 samples per year per management area 

Years<-unique(Survey$Year)
GFMA<-unique(Survey$SEdist)

  { 
  for (i in Years){    #i<-Years[1]
    for (j in GFMA){   #j<-GFMA[1]
      
      P<-Port.rand[Port.rand$Groundfish.Management.Area.Code == j,]
      #reach back from year i and see if you can get 150 weights...Rhea bumped this to 200
      Nw<-nrow(P[P$Year == i,])

      k<-0
      while (Nw < 300 & i-k >= min(Port.rand$Year) ){         #go back as far as needed to get 300 weights... 
        k<-k+1
        Nw<-nrow(P[P$Year >=i-k & P$Year <= i,])   #P[P$Year == 1990:2002,]
      }
      m<-0
      while (Nw < 300 & i+m <= max(Port.rand$Year)) {                        #go forward if you failed to find enough weights... 
        m<-m+1                                                        
        Nw<-nrow(P[P$Year >=i-k & P$Year <= i+m,])   #P[P$Year == 1986,]
      }
      #get average weights of YE...
      Sample<-P$Weight.Kilograms[P$Year >= i-k & P$Year <= i+m]
      
      #PHIL'S CODE FAILS HERE DUE TO WEIRD SUBSETTING CODE - FIXED IN FOLLOWING 3 LINES
      # Survey[,"mean.YE.kg"][Survey$Year == i & Survey$SEdist == j]<-mean(Sample)
      # Survey[,"var.YE.kg"][Survey$Year == i & Survey$SEdist == j]<-var(Sample)
      # Survey[,"N.YE.kg"][Survey$Year == i & Survey$SEdist == j]<-length(Sample)
      Survey$mean.YE.kg[Survey$Year == i & Survey$SEdist == j]<-mean(Sample)
      Survey$var.YE.kg[Survey$Year == i & Survey$SEdist == j]<-var(Sample)
      Survey$N.YE.kg[Survey$Year == i & Survey$SEdist == j]<-length(Sample)
    }
  }
}

head(Survey$mean.YE.kg)

# plot mean weight by district over time

Survey %>% group_by(SEdist,Year) %>%
  dplyr::summarise(mean_wt = mean(mean.YE.kg)) %>%
  ggplot() + 
  geom_point(aes(Year,mean_wt,col=SEdist)) +
  geom_line(aes(Year,mean_wt,col=SEdist)) 


# ******************************************************************************
# Function for generating bootstrap estimates of CPUE from IPHC surveys ----
# ******************************************************************************

YEHA.fxn <- function(Survey=Survey, Area="SEdist", Deep=250, Shallow=0,  nboot=1000){
  col.ref<-which(colnames(Survey)==Area)
  
  IPHC.cpue <- data.frame()
  Subs <- unique(Survey$SEdist) 
  Years <- unique(Survey$Year)
  
  j<-1
  for (y in Years) {  #y<-Years[26]
    for (s in Subs){  #s<-Subs[3]
      Dat<-Survey[Survey$Year == y & 
                    Survey[,col.ref] == s &
                    Survey$AvgDepth.fm > Shallow & 
                    Survey$AvgDepth.fm < Deep &
                    Survey$Eff == "Y",]
      YE_w <-unique(Dat$mean.YE.kg)
      if (nrow(Dat)>0){
        Stations<-unique(Dat$Station)
        CPUEi<-vector()
        WCPUEi<-vector()
        i<-1
        for (st in Stations){    #st<-Stations[1]   length(Stations)  st<-1
          Stat.Dat<-Dat[Dat$Station == st,] #; Stat.Dat
          #debug
          #if (nrow(Stat.Dat) > 1){aaa} else {}
          CPUE<-mean(Stat.Dat$YE.obs/Stat.Dat$HooksObserved)
          #hook adjustment factor
          CPUE<-CPUE*unique(Stat.Dat$h.adj)
          WCPUE<-CPUE*YE_w
          
          if (CPUE == 0){
            C<-0
          } else {
            C<-CPUE*Stat.Dat$HooksRetrieved
          }
          
          CPUEi[i]<-CPUE
          WCPUEi[i]<-WCPUE
          i<-i+1
        }
        
        Out<-data.frame()
        #CPUE.out<-data.frame()
        for (ii in 1:nboot){ #nboot<-1000  ii<-1
          Resamp3<-sample(CPUEi,length(CPUEi),replace=T)
          Out[ii,"CPUE"]<-mean(Resamp3, na.rm=T)
          Out[ii,"CPUE.var"]<-var(Resamp3, na.rm=T)
          Resamp4<-sample(WCPUEi,length(WCPUEi),replace=T)
          Out[ii,"WCPUE"]<-mean(Resamp4, na.rm=T)
          Out[ii,"WCPUE.var"]<-var(Resamp4, na.rm=T)
        }
        
        #hist(Out$KgHa, breaks = 25)
        IPHC.cpue[j,"Year"]<-y
        IPHC.cpue[j,"mngmt.divisions"]<-Area
        IPHC.cpue[j,"mngmt.area"]<-s
        IPHC.cpue[j,"deep.bound"]<-Deep
        IPHC.cpue[j,"shallow.bound"]<-Shallow
        IPHC.cpue[j,"no.stations"]<-length(Stations)
        
        IPHC.cpue[j,"CPUE.mean"]<-mean(CPUEi)
        IPHC.cpue[j,"CPUE.bootmean"]<-unname(quantile(Out$CPUE,c(0.5)))  #mean WCPUE from Tribuzio
        IPHC.cpue[j,"CPUE.lo95ci"]<-unname(quantile(Out$CPUE,c(0.025)))
        IPHC.cpue[j,"CPUE.hi95ci"]<-unname(quantile(Out$CPUE,c(0.975)))
        IPHC.cpue[j,"CPUE.var"]<-var(Out$CPUE)
        IPHC.cpue[j,"CPUE.cv"]<-sd(Out$CPUE)/mean(Out$CPUE)
        #    IPHC.wcpue[j,"WCPUE32.cv"]<-sqrt(var(WCPUEi.32))/mean(WCPUEi.32)
        IPHC.cpue[j,"WCPUE.mean"]<-mean(WCPUEi)
        IPHC.cpue[j,"WCPUE.bootmean"]<-unname(quantile(Out$WCPUE,c(0.5)))  #mean WCPUE from Tribuzio
        IPHC.cpue[j,"WCPUE.lo95ci"]<-unname(quantile(Out$WCPUE,c(0.025)))
        IPHC.cpue[j,"WCPUE.hi95ci"]<-unname(quantile(Out$WCPUE,c(0.975)))
        IPHC.cpue[j,"WCPUE.var"]<-var(Out$WCPUE)
        IPHC.cpue[j,"WCPUE.cv"]<-sd(Out$WCPUE)/mean(Out$WCPUE)

        j<-j+1
      } else {}
    }
  }
  return(IPHC.cpue)
}


# ******************************************************************************
# estimate CPUE using the bootstrap method ----
# ******************************************************************************

# Note: this calculation, which was used for the 2022 assessment, included only
# survey stations at which yelloweye rockfish had been observed at least once.
# The 2023 CIE review recommended including all stations, not only those at  
# which yelloweye had been observed at least once, and using the Tweedie model
# method (see next section) to calculate CPUE. The bootstrap calculation using
# non-zero stations only is preserved here for use in bridging analyses.

Survey_non0 <- Survey %>% 
  filter(Station %in% c(YE.stations_non0) & SEdist %in% c("EYKT","NSEO","CSEO","SSEO"))

SE.subdistricts_non0 <- YEHA.fxn(Survey=Survey_non0, Area="SEdist", Deep=250, Shallow=0, nboot=1000)

str(SE.subdistricts_non0)
unique(SE.subdistricts_non0$mngmt.area)

SE.subdistricts_non0 %>% filter(Year == 2023)

IPHC.index.bootstrap <- SE.subdistricts_non0 
str(IPHC.index.bootstrap)

ggplot(IPHC.index.bootstrap, aes(x=Year)) +
  geom_point(aes(y=CPUE.mean),size=2) +
  geom_point(aes(y=WCPUE.mean),size=2,col="blue") +
  geom_errorbar(aes(ymin = CPUE.lo95ci, ymax= CPUE.hi95ci),
                col="black", alpha=0.5) +
  geom_errorbar(aes(ymin = WCPUE.lo95ci, ymax= WCPUE.hi95ci),
                col="blue", alpha=0.5) +
  facet_wrap(~mngmt.area) +
  xlab("\nYear") +
  ylab("Yelloweye CPUE / WCPUE n/hook") +
  #ylab("Density (yelloweye rockfish/kmsq)") +
  scale_y_continuous(label=comma) +
  scale_x_continuous(breaks=seq(1995,2025,5)) + 
  theme (axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
         panel.grid.minor = element_blank()) +
  labs(title =~ atop("Yelloweye CPUE in IPHC FISS",scriptstyle("Stations that encountered yelloweye at least once")))

ggsave(paste0(here::here(), "/Figures/IPHC_CPUE_non0_bootstrap.png"), dpi=300,  height=5, width=5, units="in")

write.csv(IPHC.index.bootstrap, paste0(here::here(), "/Data_processing/Data/IPHC_CPUE_non0_bootstrap.csv"))


# ******************************************************************************
# model-based CPUE estimation using Tweedie model: run and compare models ----
# ******************************************************************************

# Note: include stations that have 0 yelloweye recorded in the time series. 
# Exclude stations below 250 fathoms because this is thought to be outside the 
# depth range of yelloweye.

# filter survey data for stations within the depth range
Survey_all <- Survey %>% filter(SEdist %in% c("EYKT","NSEO","CSEO","SSEO"), AvgDepth.fm <= 250)

# create data frame to use for models
IPHC_tweed.1 <- Survey_all %>% 
  group_by(Year,Station) %>%
  mutate(CPUE = h.adj*(YE.obs/HooksObserved),
         WCPUE = CPUE*mean.YE.kg) %>% 
  select(Year,Station,SEdist,IPHC.Reg.Area, IPHC.Stat.Area, IPHC.Charter.Region,        
         Purpose.Code, Date, Eff, 
         Lat = BeginLat, Lon = BeginLon, Depth = AvgDepth.fm, 
         Soak = Soak.time..min.., Temp.C, Pres = Max.Pressure..db., 
         pH, Sal = Salinity.PSU, O2_1 = Oxygen_ml,
         O2_2 = Oxygen_umol, Oxygen_sat,
         CPUE, WCPUE) %>% 
  mutate(Year = as.factor(Year),
         SEdist = as.factor(SEdist),
         Soak = as.numeric(gsub(",", "", Soak)),
         Temp.C = as.numeric(Temp.C),
         Sal = as.numeric(Sal),
         O2_1 = as.numeric(O2_1),
         O2_2 = as.numeric(O2_2),
         Oxygen_sat = as.numeric(Oxygen_sat)) %>%
  unique() %>% data.frame()

head(data.frame(IPHC_tweed.1), 30)

# see original code for correlations - some strong correlations, almost everything correlated with depth! 
# Temp, Pres pH and Salinity all correlated with each other
# probably best to just consider lat/lon, depth, and soak 
IPHC_tweed <- IPHC_tweed.1 %>% 
  select(c(Year, Station, SEdist, IPHC.Reg.Area, IPHC.Stat.Area, IPHC.Charter.Region,        
         Purpose.Code, Date, Eff, Lat, Lon, Depth, Soak, CPUE, WCPUE))

# Look at distribution of CPUE data
ggplot(IPHC_tweed, aes(WCPUE)) + geom_density(alpha = 0.4, fill = 4)
ggplot(IPHC_tweed, aes(log(WCPUE+0.01))) + geom_density(alpha = 0.4, fill = 4)
ggplot(IPHC_tweed, aes(log(WCPUE+0.1*mean(WCPUE,na.rm=T)))) + geom_density(alpha = 0.4, fill = 4)

#need complete data sets for running models: 
nrow(IPHC_tweed)
nrow(IPHC_tweed[complete.cases(IPHC_tweed),])

#only use complete cases... 
fulldat <- IPHC_tweed[complete.cases(IPHC_tweed),]

# run models
m0 <- gam(WCPUE ~ Year * SEdist, data=fulldat, gamma=1.4, family=tw(), method = "REML")
m.depth <- gam(WCPUE ~ Year * SEdist + s(Depth, k=4), data=fulldat, gamma=1.4, family=tw(), method = "REML")
m.soak <- gam(WCPUE ~ Year * SEdist + s(Soak, k=4), data=fulldat, gamma=1.4, family=tw(), method = "REML")
m.ll <- gam(WCPUE ~ Year * SEdist + te(Lon, Lat), data=fulldat, gamma=1.4, family=tw(), method = "REML")
m.depth_soak <- gam(WCPUE ~ Year * SEdist + s(Depth, k=4) + s(Soak, k=4), 
                    data=fulldat, gamma=1.4, family=tw(), method = "REML")
m.depth_ll <- gam(WCPUE ~ Year * SEdist + s(Depth, k=4) + te(Lon, Lat), 
                  data=fulldat, gamma=1.4, family=tw(), method = "REML")
m.soak_ll <- gam(WCPUE ~ Year * SEdist + s(Soak, k=4) + te(Lon, Lat) ,
                 data=fulldat, gamma=1.4, family=tw(), method = "REML")
m.global <- gam(WCPUE ~ Year * SEdist + s(Depth, k=4) + s(Soak, k=4) + te(Lon, Lat) , 
                data=fulldat, gamma=1.4, family=tw(), method = "REML")

# compare models
model.list<-list(m0,m.depth,m.soak,m.ll,
                 m.depth_soak,m.depth_ll,m.soak_ll,m.global)
names(model.list)<-c("m0","depth","soak","latlon","depth+soak","depth+latlon",
                     "soak+latlon","global")
modsum0<-data.frame(); j<-1
for (i in model.list) {
  #mod<-i
  modsum0[j,"model"]<-names(model.list[j])
  modsum0[j,"aic"]<-AIC(i)
  modsum0[j,"bic"]<-BIC(i)
  modsum0[j,"qaic"]<-QAIC(i, chat= 1/(1-summary(i)$r.sq))
  modsum0[j,"dev"]<-summary(i)$dev.expl
  modsum0[j,"rsq"]<-summary(i)$r.sq
  modsum0[j,"dev_exp"]<-summary(i)$dev.expl-summary(m0)$dev.expl
  j<-j+1
}

modsum0 %>% arrange(aic)
modsum0 %>% arrange(qaic)
modsum0 %>% arrange(bic)
modsum0 %>% arrange(-dev)  
modsum0 %>% arrange(-rsq) 

plot(m.global, page = 1, shade = TRUE, resid = TRUE, all = TRUE)
summary(m.global)

# No residual patterns, but may be some outliers
plot(fitted(m.global), resid(m.global))
abline(h = 0, col = "red", lty = 2)

plot(m.depth_soak)
plot(fitted(m.depth_soak), resid(m.depth_soak))
# Check for outliers
which(fitted(m.global) < -1.5)   

# assign model to use
mod_std_tweed <- m.global

#cpue_dat<-cpue_nom %>% select(-c(Gear,Drifts,Stat, Depth))
cpue_dat<-IPHC_tweed[complete.cases(IPHC_tweed),]
cpue_dat <- as.data.frame(IPHC_tweed)

# ******************************************************************************
# model-based CPUE estimation using Tweedie model: produce predictions ----
# ******************************************************************************

# Data set of average variables to predict CPUE from:
std_dat_tweed <- expand.grid(Year = as.factor(unique(cpue_dat$Year)),
                             SEdist = as.factor(unique(cpue_dat$SEdist)), #table(cpue_dat$Gear)
                             Depth = mean(cpue_dat$Depth),
                             Soak = mean(cpue_dat$Soak, na.rm=T), 
                             #Lat = mean(cpue_dat$Lat),
                             #Lon = mean(cpue_dat$Lon),
                             dum = 1,
                             dumstat = 1) %>%
  mutate(Lat = case_when(SEdist == "EYKT" ~ mean(cpue_dat$Lat[cpue_dat$SEdist == "EYKT"]),
                         SEdist == "NSEO" ~ mean(cpue_dat$Lat[cpue_dat$SEdist == "NSEO"]),
                         SEdist == "CSEO" ~ mean(cpue_dat$Lat[cpue_dat$SEdist == "CSEO"]),
                         SEdist == "SSEO" ~ mean(cpue_dat$Lat[cpue_dat$SEdist == "SSEO"])),
         Lon = case_when(SEdist == "EYKT" ~ mean(cpue_dat$Lon[cpue_dat$SEdist == "EYKT"]),
                         SEdist == "NSEO" ~ mean(cpue_dat$Lon[cpue_dat$SEdist == "NSEO"]),
                         SEdist == "CSEO" ~ mean(cpue_dat$Lon[cpue_dat$SEdist == "CSEO"]),
                         SEdist == "SSEO" ~ mean(cpue_dat$Lon[cpue_dat$SEdist == "SSEO"])))

# Predict CPUE
pred_cpue <- predict(mod_std_tweed, std_dat_tweed, type = "link", se = TRUE)

#checking my code with Jane's
preds<-predict.gam(mod_std_tweed, type="response", std_dat_tweed, se = TRUE)
pred_cpue
preds

#Put the standardized CPUE and SE into the data frame and convert to
#backtransformed (bt) CPUE
alpha <- 0.05  # for a 95% confidence interval on bycatch and discard estimates
z <- qnorm(1 - alpha / 2)  # Z value for 95% CI

std_dat_tweed <- std_dat_tweed %>% 
  mutate(fit = pred_cpue$fit,
         se = pred_cpue$se.fit,
         lower_link = fit - z * se,
         upper_link = fit + z * se,
         lower = mod_std_tweed$family$linkinv(lower_link),
         upper = mod_std_tweed$family$linkinv(upper_link),
         #upper = fit + (2 * se),
         #lower = fit - (2 * se),
         bt_cpue = exp(fit),
         bt_upper = upper, #exp(upper),
         bt_lower = lower, #exp(lower),
         bt_se = (bt_upper - bt_cpue) / 2  ,
         bt_cv = bt_se/bt_cpue)

# Nominal CPUE for comparison
fsh_sum_boot <- IPHC_tweed %>%
  group_by(Year,SEdist) %>%
  do(data.frame(rbind(Hmisc::smean.cl.boot(.$WCPUE)))) %>%
  mutate(calc = "set lvl kg/hook",
         fsh_cpue = Mean,
         upper = Upper,
         lower = Lower,
         cv = (upper-lower)/1.96)

std_dat_tweed %>% 
  select(Year, cpue = bt_cpue, upper = bt_upper, lower = bt_lower, SEdist) %>% 
  mutate(CPUE = "Tweedie index") %>% data.frame() -> plot_dat

plot_dat %>% data.frame() %>%
  ggplot(aes(group = 1)) + #geom_ribbon(aes(Year, ymin = lower, ymax = upper)) +
  geom_ribbon(aes(Year, ymin = lower, ymax = upper), fill = wes_palette("Darjeeling1")[c(5)],
              colour = NA, alpha = 0.5) +
  geom_point(aes(Year, cpue, colour = CPUE, shape = CPUE), size = 2, show.legend = F) +
  geom_line(aes(Year, cpue, colour = CPUE, group = CPUE), size = 1) +
  facet_wrap(~ SEdist, scales = "free") +
  scale_colour_manual(values =  wes_palette("Darjeeling1")[c(5,4)],
                      aesthetics = c("colour","fill"), name = "IPHC CPUE") +
  labs(x = "Year", y = "FISS CPUE (round kg / hook)\n") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5)) + #c(0.85, 0.9)) +
  expand_limits(y = 0)

ggsave("output/IPHC_cpue_tweedie_rke_102824.png", dpi=300, height=6, width=6, units="in")

# ******************************************************************************
# compare the bootstrap and GAM indices ----
# ******************************************************************************

boot_index <- IPHC.index.bootstrap %>% 
  mutate(SEdist = mngmt.area, fsh_cpue = WCPUE.mean, upper = WCPUE.hi95ci,
         lower = WCPUE.lo95ci, cv = (upper-lower)/1.96) %>%
  select(Year, SEdist, fsh_cpue, lower, upper, cv)

# Compare predicted cpue from GAM to nominal cpue
names(wes_palettes)

compare <- boot_index %>%
  select(Year, cpue = fsh_cpue, upper, lower, SEdist) %>% 
  mutate(CPUE = "Bootstrap index",
         Year = as.factor(Year)) %>% 
  bind_rows(std_dat_tweed %>% 
              select(Year, cpue = bt_cpue, upper = bt_upper, lower = bt_lower, SEdist) %>% 
              mutate(CPUE = "Tweedie index")) %>% 
  mutate(Year = as.numeric(as.character(Year))) %>%
  ggplot() +
  geom_ribbon(aes(Year, ymin = lower, ymax = upper, fill = CPUE), 
              colour = NA, alpha = 0.3) +
  geom_point(aes(Year, cpue, colour = CPUE, shape = CPUE), size = 2, show.legend = F) +
  geom_line(aes(Year, cpue, colour = CPUE, group = CPUE), size = 1) +
  facet_wrap(~ SEdist, scales = "free") +
  scale_colour_manual(values =  wes_palette("Darjeeling1")[c(5,4)],
                      aesthetics = c("colour","fill"), name = "IPHC CPUE") +
  labs(x = "Year", y = "FISS CPUE (round kg / hook)\n") +
  theme(legend.position = "bottom") + #c(0.85, 0.9)) +
  expand_limits(y = 0)

ggsave(paste0(here::here(), "/Figures/IPHC_CPUE_GAM_bootstrap_compare.png"), compare, dpi=300, height=6, width=6, units="in")

# Save these values
IPHC_cpue_indices <- boot_index %>%
  select(Year, SEdist, cpue = fsh_cpue, upper, lower, cv) %>% 
  mutate(CPUE = "Bootstrap index",
         Year = as.factor(Year)) %>% 
  bind_rows(std_dat_tweed %>% 
              select(Year, SEdist, cpue = bt_cpue, 
                     upper = bt_upper, lower = bt_lower, cv = bt_cv) %>% 
              mutate(CPUE = "Tweedie index")) %>% 
  #mutate(Year = as.numeric(as.character(Year)))) %>% 
  data.frame()

IPHC_cpue_indices <- left_join(IPHC_cpue_indices, IPHC.index.bootstrap %>% 
    mutate(Year = as.factor(Year)) %>%
    select(Year,SEdist = mngmt.area,no.stations),by=c("Year","SEdist"))

# Save indices for use in assessment models
write.csv(IPHC_cpue_indices, paste0(here::here(), "/Data_processing/data/IPHC_CPUE_GAM_bootstrap.csv"))

# ******************************************************************************
# end of script ----
# ******************************************************************************