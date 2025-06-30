################################################################################
## Halibut harvest reconstruction
##
## This code is used to reconstruct the halibut harvests for the SEO management
## units so that we can make estimates about yelloweye bycatch based on IPHC
## survey data.  
## Because long term halibut data management areas do not align with ADFG
## management units, there is some extrapolation that include propagation of error
## 
## This is the most latest version of the code updated last 9/23/24 - LSC
################################################################################

# set up ----
{library(plyr)
library(dplyr)
library(boot)
library(ggplot2)
library(scales)}

YEAR<-2024

#ADF&G FT Data last pulled 6/30/25 LSC

##########################################################################################
### IMPORT DATA ###
##########################################################################################


## ADF&G Halibut fish Ticket Data ########################################################

## Catch Data 1975-2020 - R output - I don't know what code this originally comes from.
HA.Harv<-read.csv("Data_processing/Data/Harvests/halibut_catch_data_new071422.csv", header=T) %>% 
  mutate(Year = year.landed, Mgt.Area = mgmt.area, ha.lbs = round.lbs,
         fishery.code = permit.fishery,
         gear = gear.description,
         ha.mt = ha.lbs*0.00045359,
         source = "old") %>%
  select(Year,Mgt.Area,fishery.code,gear,ha.lbs,ha.mt,source)

unique(HA.Harv$Mgt.Area)

## Catch Data 2021-Current Year
## https://oceanak.adfg.alaska.gov/analytics/saw.dll?PortalGo&Action=prompt&path=%2Fshared%2FCommercial%20Fisheries%2FRegion%20I%2FGroundFish%2FUser%20Reports%2FYelloweye%20Reports%20for%20Phil%2FHalibut%20harvest%20SEO%20in%20fish%20ticket%20data%202021%20-%20present
## This OceanAK report originally (Halibut harvest SEO in fish ticket data 2007-2022) 
## had data from 2007-2022 but I changed the the report to include all years 
## >= 2021 to avoid overlap with HA.Harv. I also added NSEI and SSEI because the HA.Harv has inside waters
HA.Harv.update<-read.csv("Data_processing/Data/Harvests/halibut_catch_data_6.30.25.csv") %>% 
  mutate(Year=DOL.Year,
         gear = Gear.Code.and.Name,
         ha.lbs = as.numeric(Whole.Weight..sum.),
         ha.mt = ha.lbs*0.00045359,
         fishery.code = CFEC.Fishery.Code,
         source = "new",
         Mgt.Area = Mgt.Area.District) %>% 
  select(Year,Mgt.Area,fishery.code,gear,ha.lbs,ha.mt,source)

#Data Check - combine the new and old HA harvest data
temp<-rbind(HA.Harv,HA.Harv.update)
unique(temp$Year)

ggplot(temp,aes(Year,ha.mt,col=source)) + 
  geom_line()+
  facet_wrap(~Mgt.Area)

Halibut.harv.1975<-rbind(HA.Harv %>% 
                           filter(Year < min(HA.Harv.update$Year)),HA.Harv.update) %>%
  group_by(Year,Mgt.Area) %>% 
  summarise(HA.lbs = sum(ha.lbs), HA.mt = sum(ha.mt))

#SAVE this data for use in estimating historical bycatch

write.csv(Halibut.harv.1975,paste0("Data_processing/Data/SE_Halibut_removals_",min(Halibut.harv.1975$Year),"-",
                                   max(Halibut.harv.1975$Year),".csv"))

## IPHC Survey Data #####################################################################

## Halibut harvest from IPHC 1982 -2022; Data is available from: https://www.iphc.int/data/commercial-datasets 
## and under "Pacific Halibut Directed Commercial Landings" download "IPHC Statistical 
## Area and Year - head-off, dressed weight; 1991-" with the net weight in lbs
## save it to the Data_processing/Data/Harvests/ folder and call it IPHC_harv_YEAR.csv
## The spreadsheet is heavily formatted so you'll need to do some cleaning to 
## get just the header row and the data for the most recent year.  This will then be
## added to the sheet Randy came up with (HA.req) and then saved for next year 
HA.req<-read.csv("Data_processing/Data/Harvests/Halibut_Harvest_IPHCdatareq_1982_2022.csv")

## Halibut harvest from IPHC 1929-1975; Data for 1975-1982 came from IPHC Scientific report 67 (Hoag et al. 1983)
HA.29_75<-read.csv("Data_processing/Data/Harvests/IPHC_Halibut_harv_1929-1975.csv")
unique(HA.29_75$Year)

## Halibut harvest from IPHC 1975-1982; Data for 1929-1975 came from IPHC Technical report 14 (Myhre et al. 1977) 
Ha.75_82<-read.csv("Data_processing/Data/Harvests/IPHC_Halibut_harv_1975-1982.csv")

## Halibut by IPHC area from web source data; Data came via an IPHC data request (T. Khong and I. Stewart pers. comm.)
HA.IPHCweb<-read.csv("Data_processing/Data/Harvests/Halibut_harvests_IPHCareas_1888.csv", skip=1, header=T)

################################################################################
### HALIBUT HAREVEST RECONSTRUCTION ###
################################################################################

## 2022 SAFE Report - Appendix A:
## https://apps-afsc.fisheries.noaa.gov/Plan_Team/2022/GOAdsr.pdf
## Because IPHC statistical areas do not align with ADF&G management areas it
## was necessary to reconstruct halibut harvests with some uncertainty. 
## Three of the management areas (CSEO, SSEO and NSEO) were readily identifiable 
## by IPHC regulatory area or region. However, the ADF&G EYKT management area 
## comprises only a portion of IPHC regulatory area 3A and only a portion of the
## Yakutat IPHC region. As such, hindcast harvests were constructed using the average
## proportional relationship between the 3A or Yakutat harvests and the EYKT 
## management area from 1982-2021.
## For years going back to 1929 this meant applying the proportional relationship 
## between the EYKT and the Yakutat IPHC region (which comprises a portion of 
## IPHC management area 3A). For years prior to 1929 this meant applying the 
## proportional relationship between the EYKT and the entire 3A regulatory area.

## STEP 1 ######################################################################
# get the contemporary data from IPHC request, get SEO data from there: HA.req
#         a) 2C partitioned into SEO and SEI - measure this PROPORTION (2Cprops)
#         b) 3A to SEO from FISHTICKET data proportion - measure and apply (3A prop)

HA.req_2C_3A <- HA.req %>%
  filter(IPHC.Regulatory.Area %in% c("2C", "3A")) %>%
  group_by(IPHC.Regulatory.Area, Year, IPHC.Region.2..1929.1975) %>%
  summarise(Halibut_lbs = sum(Pacific.halibut.Net.wt..000.lbs.), .groups = "drop") %>%
  mutate(Halibut_mt = Halibut_lbs * 0.00045359 * 1000,
    IO = IPHC.Region.2..1929.1975) %>%
  select(Year, IPHC.regarea = IPHC.Regulatory.Area, IO, Halibut_lbs, Halibut_mt)

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

#halibut harvest in 2C from IPHC survey
plot(data=Hal.SPM, SEO2C.req~Year, ylim=c(0,3500), type="l")

#halibut harvest in 2C from ADF&G FTs
lines(data=Hal.SPM,SEO2C.tix~Year,type="l",col="blue")

## For hindcasting proportions lets use pre-full retention
## Full retention was required of all DSR captured in groundfish 
## and halibut fisheries in federal waters in 2005
## Full retention of all DSR captured in groundfish and halibut fisheries in 
## state waters in 2009.

Pre2010<-Hal.SPM[Hal.SPM$Year < 2010,]
Prop.EYKT_3A<-mean(Pre2010$prop3A.EYKT)
var.EYKT_3A<-var(Pre2010$prop3A.EYKT)
Prop.EYKT_yak<-mean(Pre2010$prop3A.EYKT_yak)
var.EYKT_yak<-var(Pre2010$prop3A.EYKT_yak)
Prop.SEO_2C<-mean(Pre2010$prop2C.SEOreq)
var.SEO_2C<-var(Pre2010$prop2C.SEOreq)

Halibut_harvest_forSPM<- Hal.SPM %>% 
  select(Year,SEO_harvest_mt)
Halibut_harvest_forSPM$var<-0

#STEP 2 ########################################################################
## Data from 1929-1975, get SEO data plus plus small proportion of 3A

yrs29<-unique(HA.29_75$Year)

prep <- HA.29_75 %>% filter(IPHC.Region.2 == "SE-O" |
                      IPHC.Region.2 == "SE-I" |
                      IPHC.Region.2 == "Yakutat") %>%
  group_by(Year,IPHC.Region.2) %>% 
  summarize(harv.ustons = sum(Total.Catch..US.Tons..Round.Weight.)) %>%
  mutate(harv_mt = harv.ustons*0.90718474)

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

Halibut_harvest_forSPM<-rbind(Halibut_harvest_forSPM,out29)
Halibut_harvest_forSPM<-Halibut_harvest_forSPM [order(Halibut_harvest_forSPM$Year), ]

## STEP 3 ######################################################################
## Data from 1888-1928 , partition using proportions: HA.IPHCweb
##         a) apply 2Cprop and 3Aprop
##Ha.75_82 OK for SEO but missing Yakutat numbers, so will use HA.IPHC web to 
## est and add

Ha.75_82 <- Ha.75_82 %>% 
  mutate(harv_mt = 0.90718474*Total.Catch..US.Tons..Round.Weight.)

ystofill<-setdiff (min(Halibut_harvest_forSPM$Year):max(Halibut_harvest_forSPM$Year),
                   Halibut_harvest_forSPM$Year)

out76<-data.frame()
i<-1
for (y in ystofill){
  seo<-Ha.75_82[Ha.75_82$Year == y,]
  yakref<-HA.IPHCweb[HA.IPHCweb$Year == y,]
  out76[i,"Year"]<-y
  out76[i,"SEO_harvest_mt"]<-seo$harv_mt[seo$Region == "SE-O"] +
    yakref$X3A.mt*Prop.EYKT_3A
  out76[i,"var"]<-var.EYKT_3A*yakref$X3A.mt*yakref$X3A.mt
  i<-i+1
}

Halibut_harvest_forSPM<-rbind(Halibut_harvest_forSPM,out76)
Halibut_harvest_forSPM<-Halibut_harvest_forSPM [order(Halibut_harvest_forSPM$Year),]

## STEP 4 ######################################################################
eys<-unique(HA.IPHCweb$Year)[HA.IPHCweb$Year < min(Halibut_harvest_forSPM$Year)]

out1888<-data.frame()
i<-1

for (y in eys){ #y<-eys[1]
  dat<-HA.IPHCweb[HA.IPHCweb$Year == y,]
  out1888[i,"Year"]<-y
  out1888[i,"SEO_harvest_mt"]<-dat$X2C.mt*Prop.SEO_2C +
    dat$X3A.mt*Prop.EYKT_3A
  out1888[i,"var"]<-var.EYKT_3A*dat$X3A.mt*dat$X3A.mt +
    dat$X2C.mt*dat$X2C.mt*var.SEO_2C
  i<-i+1
}

Halibut_harvest_forSPM<-rbind(Halibut_harvest_forSPM,out1888)
Halibut_harvest_forSPM<-Halibut_harvest_forSPM [order(Halibut_harvest_forSPM$Year),]

## STEP 5 ######################################################################
#SAVE IT! 

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

HA.IPHCweb %>% 
  select(Year = Year,
         a2C.mt=X2C.mt,
         a3A.mt=X3A.mt) ->old.HA

IPHC_2C<-HA.req_2C_3A %>% filter(IPHC.regarea == "2C") 
IPHC_3A<-HA.req_2C_3A %>% filter(IPHC.regarea == "3A")

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
# The goal here is to see if IPHC data request and outside water harvest from 
# fishticket data match up ... data request and web sight data match up.   
# something going on since the proportion is great than 1
plot(overlap.HA$prop.2C.SEO ~ overlap.HA$Year, ylim=c(0,1), type="l")
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

plot(overlap.HA$a3A.mt ~ overlap.HA$Year, type="l", ylim=c(0,20000))
lines(overlap.HA$EYKT ~ overlap.HA$Year, type="l", col="blue")
lines(IPHC_3A$Halibut_mt[IPHC_3A$IO=="Yakutat"] ~ 
        IPHC_3A$Year[IPHC_3A$IO=="Yakutat"], 
      type="l",col="darkcyan")

