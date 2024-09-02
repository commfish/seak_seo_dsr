################################################################################
## Halibut harvest reconstruction
##
## This code is used to reconstruct the halibut harvests for the SEO management
## units so that we can make estimates about yelloweye bycatch based on IPHC
## survey data.  
## Because long term halibut data management areas do not align with ADFG
## management units, there is some extrapolation that include propagation of
## error
## This is the most latest version of the code updated last 8.30.24 - LSC
################################################################################
{library(plyr)
library(dplyr)
library(boot)
library(ggplot2)
library(scales)
  }

YEAR<-2023

{#update this from here: https://oceanak.dfg.alaska.local/analytics/saw.dll?Answers&path=%2Fshared%2FCommercial%20Fisheries%2FRegion%20I%2FGroundFish%2FUser%20Reports%2FYelloweye%20Reports%20for%20Phil%2FHalibut%20harvest%20SEO%20in%20fish%20ticket%20data%202007-2022
HA.Harv<-read.csv("Data_processing/Data/Harvests/halibut_catch_data_new071422.csv", header=T)
unique(HA.Harv$year.landed) #1975-2020
#Halibut fish ticket data:
HA.Harv.update<-read.csv("Data_processing/Data/Harvests/halibut_catch_data_8.30.24.csv") %>% 
  filter(DOL.Year>2020)
unique(HA.Harv.update$DOL.Year) #1975-2020
#Halibut by IPHC area from web source data
HA.IPHCweb<-read.csv("Data_processing/Data/Harvests/Halibut_harvests_IPHCareas_1888.csv", skip=1, header=T)
#Halibut harvest from IPHC data request 1982 - present
HA.IPHCreq<-read.csv("Data_processing/Data/Harvests/Halibut_Harvest_IPHCdatareq_1982_2022.csv")


str(HA.IPHCweb)
str(HA.IPHCreq)
str(HA.Harv)
str(HA.Harv.update)
HA.Harv.update$Year<-HA.Harv.update$DOL.Year
min(HA.Harv.update$Year)

HA.Harv %>% mutate(Year = year.landed, Mgt.Area = mgmt.area, ha.lbs = round.lbs,
                   fishery.code = permit.fishery,
                   gear = gear.description,
                   ha.mt = ha.lbs*0.00045359,
                   source = "old") %>%
  dplyr::select(Year,Mgt.Area,fishery.code,gear,ha.lbs,ha.mt,source) ->HA.Harv.old

HA.Harv.update %>% mutate(Year=DOL.Year,
                          gear = Gear.Code.and.Name,
                          ha.lbs = as.numeric(Whole.Weight..sum.),
                          ha.mt = ha.lbs*0.00045359,
                          fishery.code = CFEC.Fishery.Code,
                          source = "new") %>% 
  dplyr::select(Year,Mgt.Area,fishery.code,gear,ha.lbs,ha.mt,source) ->HA.Harv.latest

temp<-rbind(HA.Harv.old,HA.Harv.latest)

ggplot(temp,aes(Year,ha.mt,col=source)) + geom_line()+facet_wrap(~Mgt.Area)

Halibut.harv.1975<-rbind(HA.Harv.old %>% 
                           filter(Year < min(HA.Harv.latest$Year)),HA.Harv.latest) %>%
  group_by(Year,Mgt.Area) %>% 
  summarise(HA.lbs = sum(ha.lbs),
            HA.mt = sum(ha.mt))

#SAVE this data for use in estimating historical bycatch

write.csv(Halibut.harv.1975,paste0("Data_processing/Data/SE_Halibut_removals_",min(Halibut.harv.1975$Year),"-",
                                   max(Halibut.harv.1975$Year),".csv"))
}
#**********************************************************************************
#*Get SEO estimates of Halibut removals pre-1975 from RP compiled data and 
#Step 1) get the contemporary data from IPHC request, get SEO data from there: HA.IPHCreq
#         a) 2C partitioned into SEO and SEI - measure this PROPORTION (2Cprops)
#         b) 3A to SEO from FISHTICKET data proportion - measure and apply (3A prop)
# step 2) get 1929-1981 data, get SEO data plus plus small proportion of 3A
# step 3) get 1888-1928 data, partition using proportions: HA.IPHCweb
#         a) apply 2Cprop and 3Aprop

# Get FISH TICKEt halibut data
#update this from here: https://oceanak.dfg.alaska.local/analytics/saw.dll?Answers&path=%2Fshared%2FCommercial%20Fisheries%2FRegion%20I%2FGroundFish%2FUser%20Reports%2FYelloweye%20Reports%20for%20Phil%2FHalibut%20harvest%20SEO%20in%20fish%20ticket%20data%202007-2022

# HA.Harv<-read.csv("Data_processing/Data/Harvests/halibut_catch_data_new071422.csv", header=T)
# #Halibut fish ticket data:
# #HA.fishtix<-read.csv("Data/Halibut harvest from fish ticket data_071522.csv")
# HA.Harv.update<-read.csv("Data_processing/Data/Harvests/Halibut harvest from fish ticket data_071522.csv") 

#Get IPHC HALIBUT DATA: 

#Halibut harvest from IPHC data request 1982 - present
# Available from: https://www.iphc.int/data/commercial-datasets and 
# under "Pacific Halibut Directed Commercial Landings" download
# "IPHC Statistical Area and Year - head-off, dressed weight; 1991-" with the net weight in lbs
# save it to the Data_processing/Data/Harvests/ folder and call it IPHC_harv_YEAR.csv
# The spreadsheet is heavily formatted so you'll need to do some cleaning to 
# get just the header row and the data for the most recent year.  This will then be
# added to the sheet Randy came up with (HA.req) and then saved for next year 

#Data available on the web is the same as what was used for the last assessment
#I put in a request for updated numbers.
#HA.newreq<-read.csv(paste0("Data_processing/Data/Harvests/IPHC_harv_",YEAR-1,".csv"))
HA.newreq<-read.csv("Data_processing/Data/Harvests/IPHC_harv_2022.csv")

#get the request Randy put in... 
HA.req<-read.csv("Data_processing/Data/Harvests/Halibut_Harvest_IPHCdatareq_1982_2022.csv")
#Halibut harvest from IPHC 1929-1975; IPHC Scientific report 67
HA.29_75<-read.csv("Data_processing/Data/Harvests/IPHC_Halibut_harv_1929-1975.csv")
#Halibut harvest IPHC 1975-1982
Ha.75_82<-read.csv("Data_processing/Data/Harvests/IPHC_Halibut_harv_1975-1982.csv")
#Halibut by IPHC area from web source data
HA.web<-read.csv("Data_processing/Data/Harvests/Halibut_harvests_IPHCareas_1888.csv", skip=1, header=T)


#-------------------------------------------------------------------------------
# 
#This step was already done above....remove this section?
# HA.Harv %>% mutate(Year = year.landed, Mgt.Area = mgmt.area, ha.lbs = round.lbs,
#                    fishery.code = permit.fishery,
#                    gear = gear.description,
#                    ha.mt = ha.lbs*0.00045359,
#                    source = "old") %>%
#   dplyr::select(Year,Mgt.Area,fishery.code,gear,ha.lbs,ha.mt,source) ->HA.Harv.old
# 
# HA.Harv.update %>% mutate(Year=DOL.Year,
#                           gear = Gear.Code.and.Name,
#                           ha.lbs = Whole.Weight..sum.,
#                           ha.mt = ha.lbs*0.00045359,
#                           fishery.code = CFEC.Fishery.Code,
#                           source = "new") %>% 
#   dplyr::select(Year,Mgt.Area,fishery.code,gear,ha.lbs,ha.mt,source) ->HA.Harv.latest
# 
# temp<-rbind(HA.Harv.old,HA.Harv.latest)

ggplot(temp,aes(Year,ha.mt,col=source)) + geom_line()+facet_wrap(~Mgt.Area)

Halibut.harv.1975<-rbind(HA.Harv.old %>% 
                           filter(Year < min(HA.Harv.latest$Year)),HA.Harv.latest) %>%
  group_by(Year,Mgt.Area) %>% 
  summarise(HA.lbs = sum(ha.lbs),
            HA.mt = sum(ha.mt))

ch98<-HA.Harv[HA.Harv$year.landed == 1998,]

#SAVE this data for use in estimating historical bycatch

write.csv(Halibut.harv.1975,paste0("Data_processing/Data/SE_Halibut_removals_",min(Halibut.harv.1975$Year),"-",
                                   max(Halibut.harv.1975$Year),".csv"))

#STEP 1------------------------------------------------------------------------
# get the contemporary data from IPHC request, get SEO data from there: HA.IPHCreq
#         a) 2C partitioned into SEO and SEI - measure this PROPORTION (2Cprops)
#         b) 3A to SEO from FISHTICKET data proportion - measure and apply (3A prop)

IOs<-unique(HA.req %>% filter(IPHC.Regulatory.Area == "2C" | IPHC.Regulatory.Area == "3A") %>%
              select(IPHC.Statistical.Area,IPHC.Stat.Area..1929.1975,IPHC.Region.2..1929.1975))

HA.newreq %>% filter(IPHC.Regulatory.Area == "2C" | IPHC.Regulatory.Area == "3A") %>% 
  mutate(Net_lbs=Net.wt..lb./1000) %>%
#         Halibut_mt = NA) %>%
  dplyr::select(Year, 
                IPHC.Regulatory.Area,
                IPHC.Statistical.Area = IPHC..Statistical.Area,
                Net_lbs) -> HA.newreq

HA.newreq <- plyr::join(HA.newreq,IOs,by="IPHC.Statistical.Area")

head(HA.req)
head(HA.newreq)

HA.req<-rbind(HA.req %>% select(Year, IPHC.Regulatory.Area, IPHC.Statistical.Area,
                                Net_lbs = Pacific.halibut.Net.wt..000.lbs.,
                                IPHC.Stat.Area..1929.1975,
                                IPHC.Region.2..1929.1975),
              HA.newreq) %>% mutate(Net_lbs = as.numeric(Net_lbs))
str(HA.req)

HA.req %>% filter(IPHC.Regulatory.Area == "2C" | IPHC.Regulatory.Area == "3A") %>%
  group_by(IPHC.Regulatory.Area, Year, IPHC.Region.2..1929.1975) %>%
  dplyr::summarise(Halibut_lbs = sum(Net_lbs)) %>%
  mutate(Halibut_mt = Halibut_lbs*0.00045359*1000,
         IO = IPHC.Region.2..1929.1975) %>% 
  dplyr::select(Year,IPHC.regarea=IPHC.Regulatory.Area,IO,Halibut_lbs,
                Halibut_mt)-> HA.req_2C_3A

HA.req_2C_3A %>% filter(Year == 2021 | Year == 2022)

Hal.SPM<-data.frame()

years<-unique(HA.req_2C_3A$Year) 
i<-1
for (y in years) {  #y<-years[1]
  subha<-Halibut.harv.1975[Halibut.harv.1975$Year == y,]
  req<-HA.req_2C_3A[HA.req_2C_3A$Year == y,]
  Hal.SPM[i,"Year"]<-y
  Hal.SPM[i,"3A.harv"]<-sum(req$Halibut_mt[req$IPHC.regarea == "3A"])
  Hal.SPM[i,"3Ayak.harv"]<-sum(req$Halibut_mt[req$IPHC.regarea == "3A" & 
                                                req$IO == "Yakutat"])
  Hal.SPM[i,"2C.harv"]<-sum(req$Halibut_mt[req$IPHC.regarea == "2C"])
  Hal.SPM[i,"EYKT"]<-subha$HA.mt[subha$Mgt.Area == "EYKT"]
  Hal.SPM[i,"SEO2C.tix"]<-sum(subha$HA.mt[subha$Mgt.Area == "NSEO" |
                                            subha$Mgt.Area == "CSEO" |
                                            subha$Mgt.Area == "SSEO"])
  Hal.SPM[i,"SEI2C.tix"]<-sum(subha$HA.mt[subha$Mgt.Area == "NSEI" |
                                            subha$Mgt.Area == "SSEI"], na.rm=T)
  Hal.SPM[i,"SEO2C.req"]<-req$Halibut_mt[req$IO == "SE-O"]
  Hal.SPM[i,"SEI2C.req"]<-req$Halibut_mt[req$IO == "SE-I"]
  
  Hal.SPM[i,"prop3A.EYKT"]<-Hal.SPM[i,"EYKT"]/Hal.SPM[i,"3A.harv"]
  Hal.SPM[i,"prop3A.EYKT_yak"]<-Hal.SPM[i,"EYKT"]/Hal.SPM[i,"3Ayak.harv"]
  Hal.SPM[i,"prop2C.SEOtix"]<-Hal.SPM[i,"SEO2C.tix"]/Hal.SPM[i,"2C.harv"]
  Hal.SPM[i,"prop2C.SEItix"]<-Hal.SPM[i,"SEI2C.tix"]/Hal.SPM[i,"2C.harv"]
  Hal.SPM[i,"prop2C.SEOreq"]<-Hal.SPM[i,"SEO2C.req"]/Hal.SPM[i,"2C.harv"]
  Hal.SPM[i,"prop2C.SEIreq"]<-Hal.SPM[i,"SEI2C.req"]/Hal.SPM[i,"2C.harv"]
  
  Hal.SPM[i,"SEO_harvest_mt"]<- Hal.SPM[i,"SEO2C.req"]+Hal.SPM[i,"EYKT"]
  
  i<-i+1
}

plot(data=Hal.SPM, SEO2C.req~Year, ylim=c(0,3000), type="l")
lines(data=Hal.SPM,SEO2C.tix~Year,type="l",col="blue")

#for hindcasting proportions lets use pre-full retention
Pre2010<-Hal.SPM[Hal.SPM$Year < 2011,]
Prop.EYKT_3A<-mean(Pre2010$prop3A.EYKT)
var.EYKT_3A<-var(Pre2010$prop3A.EYKT)
Prop.EYKT_yak<-mean(Pre2010$prop3A.EYKT_yak)
var.EYKT_yak<-var(Pre2010$prop3A.EYKT_yak)
Prop.SEO_2C<-mean(Pre2010$prop2C.SEOreq)
var.SEO_2C<-var(Pre2010$prop2C.SEOreq)

Halibut_harvest_forSPM<- Hal.SPM %>% 
  dplyr::select(Year,SEO_harvest_mt)
Halibut_harvest_forSPM$var<-0

#STEP 2------------------------------------------------------------------------
# Data from 1929-1975
# 
str(HA.29_75)
head(Hal.SPM)
yrs29<-unique(HA.29_75$Year)

HA.29_75 %>% filter(IPHC.Region.2 == "SE-O" |
                      IPHC.Region.2 == "SE-I" |
                      IPHC.Region.2 == "Yakutat") %>%
  group_by(Year,IPHC.Region.2) %>% 
  dplyr::summarize(harv.ustons = sum(Total.Catch..US.Tons..Round.Weight.)) %>%
  mutate(harv_mt = harv.ustons*0.90718474) -> prep

out29<-data.frame()
yrs29<-unique(HA.29_75$Year)
i<-1
for (y in yrs29){
  dat<-prep %>% filter(Year == y)
  out29[i,"Year"]<-y
  out29[i,"SEO_harvest_mt"]<-dat$harv_mt[dat$IPHC.Region.2 == "SE-O"] +
    dat$harv_mt[dat$IPHC.Region.2 == "Yakutat"]*Prop.EYKT_yak
  out29[i,"var"]<-var.EYKT_yak*dat$harv_mt[dat$IPHC.Region.2 == "Yakutat"]*dat$harv_mt[dat$IPHC.Region.2 == "Yakutat"]
  i<-i+1
}
sqrt(164548)

Halibut_harvest_forSPM<-rbind(Halibut_harvest_forSPM,out29)
Halibut_harvest_forSPM<-Halibut_harvest_forSPM [order(Halibut_harvest_forSPM$Year), ]

#STEP 3-------------------------------------------------------------------------
str(Ha.75_82)
str(HA.web)
#Ha.75_82 OK for SEO but missing Yakutat numbers, so will use HA.web to est and add
Ha.75_82<-Ha.75_82 %>% mutate(harv_mt = 0.90718474*Total.Catch..US.Tons..Round.Weight.)

ystofill<-setdiff (min(Halibut_harvest_forSPM$Year):max(Halibut_harvest_forSPM$Year),
                   Halibut_harvest_forSPM$Year)

out76<-data.frame()
i<-1
for (y in ystofill){
  seo<-Ha.75_82[Ha.75_82$Year == y,]
  yakref<-HA.web[HA.web$Year == y,]
  out76[i,"Year"]<-y
  out76[i,"SEO_harvest_mt"]<-seo$harv_mt[seo$Region == "SE-O"] +
    yakref$X3A.mt*Prop.EYKT_3A
  out76[i,"var"]<-var.EYKT_3A*yakref$X3A.mt*yakref$X3A.mt
  i<-i+1
}

Halibut_harvest_forSPM<-rbind(Halibut_harvest_forSPM,out76)
Halibut_harvest_forSPM<-Halibut_harvest_forSPM [order(Halibut_harvest_forSPM$Year),]

#STEP4---------------------------------------------------------------------------
eys<-unique(HA.web$Year)[HA.web$Year < min(Halibut_harvest_forSPM$Year)]

out1888<-data.frame()
i<-1

for (y in eys){ #y<-eys[1]
  dat<-HA.web[HA.web$Year == y,]
  out1888[i,"Year"]<-y
  out1888[i,"SEO_harvest_mt"]<-dat$X2C.mt*Prop.SEO_2C +
    dat$X3A.mt*Prop.EYKT_3A
  out1888[i,"var"]<-var.EYKT_3A*dat$X3A.mt*dat$X3A.mt +
    dat$X2C.mt*dat$X2C.mt*var.SEO_2C
  i<-i+1
}

Halibut_harvest_forSPM<-rbind(Halibut_harvest_forSPM,out1888)
Halibut_harvest_forSPM<-Halibut_harvest_forSPM [order(Halibut_harvest_forSPM$Year),]
#Step5----------------------------------------------------------------------
#SAVE IT! 
str(Halibut_harvest_forSPM)

ggplot(Halibut_harvest_forSPM,aes(Year,SEO_harvest_mt,col = "black")) +
  geom_ribbon(aes(ymin=SEO_harvest_mt-sqrt(var)*1.96, ymax=SEO_harvest_mt+sqrt(var)*1.96),
              alpha=0.5, col="grey") +
  geom_line(col="black") +
  xlab("\nYear") +
  ylab("SEO Halibut harvest (t)") +
  scale_y_continuous(label=comma, breaks = seq(0,6000,1000)) +
  scale_x_continuous(breaks=seq(1880,2022,5)) + 
  theme (panel.grid.minor = element_blank(),
         axis.text.x = element_text(angle = 45, vjust=1, hjust=1))

ggsave(paste0("Figures/Halibut_Harvest_SEO_", YEAR, ".png"), dpi=300,  height=3, width=6, units="in")

write.csv(Halibut_harvest_forSPM,paste0("Data_processing/Data/SEO_Halibut_removals_",
                                        min(Halibut_harvest_forSPM$Year),"-",
                                        max(Halibut_harvest_forSPM$Year),".csv"))

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
## Data check stuff here
str(HA.Harv.IPHC)
str(Halibut.harv.1975)
HA.web %>% 
  dplyr::select(Year = Year,a2C.mt=X2C.mt,a3A.mt=X3A.mt) ->old.HA

str(HA.IPHCreq)

HA.req %>% filter(IPHC.Regulatory.Area == "2C" | IPHC.Regulatory.Area == "3A") %>%
  group_by(IPHC.Regulatory.Area,Year,IPHC.Region.2..1929.1975) %>%
  dplyr::summarise(Halibut_lbs = sum(Net_lbs)) %>%
  mutate(Halibut_mt = Halibut_lbs*0.00045359*1000,
         IO = IPHC.Region.2..1929.1975)-> HA.IPHCreq_2C_3A

IPHC_2C<-HA.IPHCreq_2C_3A %>% filter(IPHC.Regulatory.Area == "2C") 
IPHC_3A<-HA.IPHCreq_2C_3A %>% filter(IPHC.Regulatory.Area == "3A")

plot(old.HA$a2C.mt ~ old.HA$Year, type="l")
lines(IPHC_2C$Halibut_mt[IPHC_2C$IO=="SE-I"]+IPHC_2C$Halibut_mt[IPHC_2C$IO=="SE-O"] ~
        IPHC_2C$Year[IPHC_2C$IO=="SE-I"],type="l",col="orange")
plot(old.HA$a3A.mt ~ old.HA$Year, type="l")
lines(IPHC_3A$Halibut_mt[IPHC_3A$IO=="Yakutat"]+IPHC_3A$Halibut_mt[IPHC_3A$IO=="Kodiak"] ~
        IPHC_3A$Year[IPHC_3A$IO=="Yakutat"],type="l",col="orange")

overlap.HA<-old.HA %>% filter(Year > 1974)

years<-unique(overlap.HA$Year) 
i<-1
for (y in years) {  #y<-years[1]
  subha<-Halibut.harv.1975[Halibut.harv.1975$Year == y,]
  overlap.HA[i,"year.check"]<-y
  overlap.HA[i,"EYKT"]<-subha$HA.mt[subha$Mgt.Area == "EYKT"]
  overlap.HA[i,"SEO2C"]<-sum(subha$HA.mt[subha$Mgt.Area == "NSEO" |
                                           subha$Mgt.Area == "CSEO" |
                                           subha$Mgt.Area == "SSEO"])
  overlap.HA[i,"SEI2C"]<-sum(subha$HA.mt[subha$Mgt.Area == "NSEI" |
                                           subha$Mgt.Area == "SSEI"], na.rm=T)
  overlap.HA[i,"prop.3A.EYKT"]<-overlap.HA[i,"EYKT"]/overlap.HA[i,"a3A.mt"]
  overlap.HA[i,"prop.2C.SEO"]<-overlap.HA[i,"SEO2C"]/overlap.HA[i,"a2C.mt"]
  overlap.HA[i,"prop.2C.SEI"]<-overlap.HA[i,"SEI2C"]/overlap.HA[i,"a2C.mt"]
  overlap.HA[i,"2c.ItoO.ratio"]<-overlap.HA[i,"SEO2C"]/(overlap.HA[i,"SEO2C"]+overlap.HA[i,"SEI2C"])
  i<-i+1
}

plot(overlap.HA$prop.2C.SEO ~ overlap.HA$Year, ylim=c(0,0.8), type="l")
lines(overlap.HA$prop.3A.EYKT ~ overlap.HA$Year, col="blue")
lines(overlap.HA$prop.2C.SEI ~ overlap.HA$Year, col="red")

plot(overlap.HA$a2C.mt ~ overlap.HA$Year, type="l", ylim=c(0,6000))
lines(overlap.HA$SEO2C ~ overlap.HA$Year, type="l", col="blue")
lines(IPHC_2C$Halibut_mt[IPHC_2C$IO=="SE-O"] ~ 
        IPHC_2C$Year[IPHC_2C$IO=="SE-O"], 
      type="l",col="darkcyan")
lines(IPHC_2C$Halibut_mt[IPHC_2C$IO=="SE-I"] ~ 
        IPHC_2C$Year[IPHC_2C$IO=="SE-I"], 
      type="l",col="purple")
lines(IPHC_2C$Halibut_mt[IPHC_2C$IO=="SE-I"]+IPHC_2C$Halibut_mt[IPHC_2C$IO=="SE-O"] ~
        IPHC_2C$Year[IPHC_2C$IO=="SE-I"],type="l",col="orange")
##HURRAY!!! IPHC data request and outside water harvest from fishticket data match
## up ... data request and web sight data match up.  Can use the SEO proportions to hindcast... 

plot(overlap.HA$a3A.mt ~ overlap.HA$Year, type="l", ylim=c(0,20000))
lines(overlap.HA$EYKT ~ overlap.HA$Year, type="l", col="blue")
lines(IPHC_3A$Halibut_mt[IPHC_3A$IO=="Yakutat"] ~ 
        IPHC_3A$Year[IPHC_3A$IO=="Yakutat"], 
      type="l",col="darkcyan")


#-------------------------------------------------------------------------------
#get inside-outside ratios 


#*******************************************************************************

mean(overlap.HA$prop.3A.EYKT); sd(overlap.HA$prop.3A.EYKT)
mean(overlap.HA$prop.2C.SEO); sd(overlap.HA$prop.2C.SEO)


ys<-c(1975,1982,1988,1993,1996,1998,2001,2003,2005,2007,2009,2011,2013,2015,2017,2019,2021)
for (y in 1:(length(ys)-1)){
  tx<-read.csv(paste0("Data/Harvests/1975-2020/",ys[y],"-",ys[y+1]-1,".csv"),
               stringsAsFactors=T)
  if(y == 1){
    haltix<-tx
  } else {
    fishtix<-rbind(tx,haltix)
  }
}
str(haltix)

unique(haltix$Stat.Area)
unique(haltix$IPHC.Regulatory.Area)
unique(haltix$IPHC.Statistical.Area)
unique(haltix$NMFS.Area)
unique(haltix$IFQ.Halibut.Area)

Halibut.harv.1975[Halibut.harv.1975$Year == 2000,]
old.HA[old.HA$Year == 2000,]

#******************************************************************************

str(YE.Harv)
Years<-unique(YE.Harv[,1])
MA<-unique(YE.Harv[,2])
unique(YE.Harv[,3])
FISH<-unique(YE.Harv[,4])
unique(YE.Harv[,5])
YE.Harv$h.code<-substr(YE.Harv$fishery,1,1)

str(YE.Subs)

str(YE.Sport)
unique(YE.Sport$species)
YE.Sport<-YE.Sport[YE.Sport$species == 145,]


str(HA.Harv)
