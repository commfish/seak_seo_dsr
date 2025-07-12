################################################################################
## Halibut harvest reconstruction
##
## This code is used to reconstruct the halibut harvests for the SEO management
## units so that we can make estimates about yelloweye bycatch based on IPHC
## survey data. 
##
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
## 
## This is the most latest version of the code updated last 9/23/24 - LSC
################################################################################

# set up ----

{library(plyr)
library(dplyr)
library(boot)
library(ggplot2)
library(scales)
library(tidyverse)
  library(extrafont)}

###  set plotting theme to use TNR  ###
#font_import() #remove # to run this but only do this one time - it takes a while
loadfonts(device="win")
windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_bw(base_size=18,base_family='Times New Roman')
          +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

YEAR<-2024

options(scipen = 999)

################################################################################
### IMPORT DATA ###
################################################################################

## ADF&G Halibut fish Ticket Data ##############################################

## Catch Data 1975-2023 - R output - code from Rhea Ehresmann using lbs from the
## CFEC gross earning subject area in OceanAK:
## https://oceanak.adfg.alaska.gov/analytics/saw.dll?Answers&path=%2Fshared%2FCommercial%20Fisheries%2FRegion%20I%2FGroundFish%2FUser%20Reports%2FYelloweye%20Reports%20for%20Phil%2FCFEC%20Gross%20Earnings#resultsTab197eba52b18
## I have updated this spreadsheet through 2023 but we will only be using the data
## before 2006 as the fish ticket data (below) is more accurate.We have a lot more 
## flexibility with the fish ticket data and filters. The CFEC gross earnings 
## doesn’t change if fish tickets are edited the following year, so the CFEC data aren’t updated 
## even if the fish tickets are updated. 
## This does not need to be updated in the future but in case someone does want to
## The data is stored: Data_processing/Data/Harvests/CFEC Gross Earnings Data
## The code is stored: Data_processing/Code/ye_hal_summary.R
## There is an additional copy on the GF drive in Rhea's folder.
## Previous versions of this spreadsheet are called: 
## halibut_catch_data_new071422
## halibut_catch_data
HA.Harv<-read.csv("Data_processing/Data/Harvests/halibut_catch_data_cfec.csv", header=T) %>%
  mutate(Year = year.landed, Mgt.Area = mgmt.area, ha.lbs = round.lbs,
         fishery.code = permit.fishery,
         gear = gear.description,
         ha.mt = ha.lbs*0.00045359,
         source = "old") %>%
  select(Year,Mgt.Area,fishery.code,gear,ha.lbs,ha.mt,source) %>% 
  filter(Year <= 2006)## i removed the filter from the source code for longline and other gear.

lapply(HA.Harv[c("Year", "Mgt.Area", "fishery.code", "gear")], unique)

## Catch Data 2021-Current Year

## https://oceanak.adfg.alaska.gov/analytics/saw.dll?PortalGo&Action=prompt&path=%2Fshared%2FCommercial%20Fisheries%2FRegion%20I%2FGroundFish%2FUser%20Reports%2FYelloweye%20Reports%20for%20Phil%2FHalibut%20harvest%20SEO%20in%20fish%20ticket%20data%202007%20-%20present
## OceanAK report originally (Halibut harvest SEO in fish ticket data 2007-2022)
## excluded inside waters and included only halibut trips (B permits). I changed filters to match
## Ha.Harv output. I saved this OceanAK report as "Halibut harvest SEO in fish
## ticket data 2007 - present". 
## I made an archive folder for old outputs. I think these can eventually be deleted,
## but I am a data hoarder and it's nice to rerun with old outputs to recreate outputs.
## Data pulled: 7.8.25
HA.Harv.update<-read.csv("Data_processing/Data/Harvests/halibut_catch_adfg_ft_data.csv") %>%
  mutate(Year=DOL.Year,
         gear = Gear.Code.and.Name,
         ha.lbs = as.numeric(Whole.Weight..sum.),
         ha.mt = ha.lbs*0.00045359,
         fishery.code = CFEC.Fishery.Code,
         source = "new") %>%
  select(Year,Mgt.Area,fishery.code,gear,ha.lbs,ha.mt,source) %>% 
  mutate(
    gear = case_when(
      gear == "61 - Longline"~"Longline",
      gear == "05 - Hand troll"~"Hand troll",
      gear == "25 - Dinglebar troll"~"Dinglebar troll",
      gear == "15 - Power gurdy troll"~"Power gurdy troll",
      gear == "91 - Pot"~"Pot",
      gear == "99 - Other/unspecified/missing"~"Other/unspecified/missing",
      gear == "26 - Mechanical jigs"~"Mechanical jigs",TRUE~gear))

lapply(HA.Harv.update[c("Year", "Mgt.Area", "fishery.code", "gear")], unique)

## IPHC Fishery  ###############################################################

## Halibut harvest from IPHC for 2022 and beyond
## Data is available from: https://www.iphc.int/data/commercial-datasets 
## and under "Pacific Halibut Directed Commercial Landings" download "IPHC Statistical 
## Area and Year - head-off, dressed weight; 1991-" with the net weight in lbs
## save it to the Data_processing/Data/Harvests/ folder and call it IPHC_harv_YEAR.csv
## The spreadsheet is heavily formatted so you'll need to do some cleaning to 
## get just the header row and the data for the most recent year.  This will then be
## added to the sheet Randy came up with (HA.req) and then saved for next year
## An updated file can be found under “Fishery information” at 
## https://www.iphc.int/data/time-series-datasets/, 
## specifically https://www.iphc.int/uploads/2024/09/iphc-2024-tsd-026.xlsx. These 
## data are derived from logbook data and fish tickets. This report in particular 
## lags by about a year. 2024 logs will still be collected well into 2025 and 
## IPHC will generate an update sometime in the fall.

HA.newreq<-read.csv("Data_processing/Data/Harvests/IPHC_harv_2023.csv")
unique(HA.newreq$Year) #2022 - 2023

## Halibut harvest from IPHC for 1982 - 2021; Randy Peterson originally put in this request 
HA.req<-read.csv("Data_processing/Data/Harvests/Halibut_Harvest_IPHCdatareq_1982_2021.csv")
unique(HA.req$Year) #1982-2021

## Halibut harvest from IPHC 1975-1982; Data for 1929-1975 came from IPHC Technical report 14 (Myhre et al. 1977) 
Ha.75_82<-read.csv("Data_processing/Data/Harvests/IPHC_Halibut_harv_1975-1982.csv")
unique(Ha.75_82$Year) #1971-1982

## Halibut harvest from IPHC 1929-1975; Data for 1975-1982 came from IPHC Scientific report 67 (Hoag et al. 1983)
HA.29_75<-read.csv("Data_processing/Data/Harvests/IPHC_Halibut_harv_1929-1975.csv")
unique(HA.29_75$Year) #1929-1975

## Halibut by IPHC area from web source data
HA.IPHCweb<-read.csv("Data_processing/Data/Harvests/Halibut_harvests_IPHCareas_1888.csv", skip=1, header=T)

################################################################################
### HALIBUT HAREVEST ADFG FISH TICKET DATA EXPLORATION ###
################################################################################

#Data Check - combine the new and old HA harvest data
hal_fishery<-rbind(HA.Harv,HA.Harv.update)
unique(hal_fishery$Year)
unique(hal_fishery$fishery.code)
unique(hal_fishery$gear)

ggplot(hal_fishery,aes(Year,ha.mt,col=source)) +
  geom_line()+
  facet_wrap(~Mgt.Area)

## I am trying to sort out what is bycatch and what are directed trips because here's
## the thing, there are dual target trips in sablefish, DSR, and p.cod fisheries.

## In the state sablefish fisheries, fishers with halibut IFQ in regulatory area 
## 2C and a CFEC halibut permit card MUST retain all halibut over 32 inches in 
## length, up to the amount of their IFQ.

## Halibut incidentally taken during an open commercial halibut season by power 
#and hand troll gear operated for salmon consistent with applicable state laws 
## and regulations are legally taken and possessed (5 AAC 28.133[c]). Commercial 
## halibut may be retained only by Individual Fishing Quota (IFQ) permit holders 
## during the open season for halibut. Trollers making an IFQ halibut landing of 
## 500 lb or less of IFQ weight as determined pursuant to 50 CFR 679.40(h) are exempted
## from the 3-hour prior notice of landing if landed concurrently with a legal landing of salmon
## harvested using hand troll or power troll gear (50 CFR 679.5[l][1][iv][A]).


## Assigning fishery types based on fishery codes and I have only included the codes
## in the data (n=38)

hal_directed <- c("B06B","B61B","B05B","B99B","B26B","B25B","B61Z","B06Z","B09B","B91B")
yelloweye <- c("Y06A")
sablefish <- c("C61B","C09B","C06B","C61C","C61A","C50B" )
lingcod <- c("I25B","I05B")
salmon <- c("S05B","S15B","S03A","S01A","S04D")
misc <- c("M06B","M61B","M99B","L99B","M05B","S04D","M26B","P09B","D09B","G34A", "M07B","M06G")
unknown <- c("","9998","9999")

#Website for historical CFEC fishery codes: https://www.cfec.state.ak.us/misc/FshyDesH.htm
#"M06B", "M61B" & "M06G" = longline (potential pacific cod trips)
#"M05B" = hand troll
#"M26B" = mechanical jig 
#"P09B" = shrimp pot
#"D09B" = Dungeness
#"G34A" = Herring gill net
#"M07B" = trawl
#"M99B" = Experimental/Special Permit 1975-2021
#"L99B" = herring spawn on kelp (1974-1980)
#Unknown:
#"9998" = this is an interim value that is supposed  be replace with a valid value within 72 hrs
#"9999"

pal <- c("Halibut Directed"        = "#91D5DE",
         "Yelloweye Incidental"    = "#2E8289",
         "Sablefish Incidental"    = "#B4674E",
         "Lingcod Incidental"      = "#EAAE37",
         "Salmon Incidental"       = "#682C37",
         "Misc Incidental"         = "#606060",
         "Unknown"                 = "#606090")

hal_fishery <- hal_fishery %>% 
  mutate(fishery = ifelse(fishery.code %in% yelloweye, "Yelloweye Incidental",
                          #there is one yelloweye incidental from 2025, the permit holder was allowed to retain and sell the halibut
                          #even tho there are not halibut bycatch allowances for the DSR fishery because the halibut was tagged
                          ifelse(fishery.code %in% hal_directed, "Halibut Directed",
                                 ifelse(fishery.code %in% sablefish, "Sablefish Incidental",
                                        ifelse(fishery.code %in% lingcod, "Lingcod Incidental",
                                               ifelse(fishery.code %in% salmon, "Salmon Incidental",
                                                      ifelse(fishery.code %in% misc, "Misc Incidental","Unknown"))))))) %>% 
  mutate(fishery = factor(fishery, levels = names(pal)))

check <- hal_fishery %>% 
  filter(fishery.code %in% c("9998","9999",""))

unique(check$gear)

sector_fig <- hal_fishery %>% 
  ggplot(aes(x=Year, y=ha.mt)) +
  geom_col(aes(fill = fishery), width = 0.7) +
  scale_fill_manual(values = pal) + 
  ylab("Catch (mt)\n") + xlab("Year\n") +
  coord_cartesian(ylim=c(0, 7000)) +
  scale_y_continuous(label=scales::comma, breaks = seq(0, 7000, 500)) +
  scale_x_continuous(breaks=seq(1975, 2025, 5)) +
  theme(legend.title = element_blank(),
        legend.position = c(.95, .95),                #Position the legend within the plot====         
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)) ; sector_fig

incidental <- hal_fishery %>% 
  filter(!fishery=="Halibut Directed") %>% 
  ggplot(aes(x=Year, y=ha.mt)) +
  geom_col(aes(fill = fishery), width = 0.7) +
  scale_fill_manual(values = pal) + 
  ylab("Catch (mt)\n") + xlab("Year\n") + ggtitle("'Incidental' Fisheries")+
  coord_cartesian(ylim=c(0, 1700)) +
  scale_y_continuous(label=scales::comma, breaks = seq(0, 1700, 100)) +
  scale_x_continuous(breaks=seq(1975, 2025, 5)) +
  theme(legend.title = element_blank(),
        legend.position = c(.95, .95),                #Position the legend within the plot====         
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)) ; incidental

non_other_incidental <- hal_fishery %>% 
  filter(!fishery %in% c("Halibut Directed","Unknown")) %>% 
  ggplot(aes(x=Year, y=ha.mt)) +
  geom_col(aes(fill = fishery), width = 0.7) +
  scale_fill_manual(values = pal) + 
  ylab("Catch (mt)\n") + xlab("Year\n") + ggtitle("Non-Unknown Incidental Fisheries")+
  coord_cartesian(ylim=c(0, 100)) +
  scale_y_continuous(label=scales::comma, breaks = seq(0, 100, 5)) +
  scale_x_continuous(breaks=seq(1975, 2025, 5)) +
  theme(legend.title = element_blank(),
        legend.position = c(.95, .95),                #Position the legend within the plot====         
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)) ; non_other_incidental

pal2 <- c("#91D5DE", "#2E8289", "#B4674E",  "#EAAE37", "#682C37", "#606060", "#606")
#other means the fishery code was missing (all data from cfec) or the code 9999/9998 was used (all data after 2007)
unknown <- hal_fishery %>% 
  filter(fishery=="Unknown") %>% 
  mutate(gear = if_else(gear == "61 - Longline", "Longline", 
                        if_else(gear == "15 - Power gurdy troll", "Power gurdy troll",gear))) %>% 
  ggplot(aes(x=Year, y=ha.mt)) +
  geom_col(aes(fill = gear), width = 0.7) +
  scale_fill_manual(values = pal2) + 
  ylab("Catch (t)\n") + xlab("Year\n") + ggtitle("Unknown Fishery by Gear Type") +
  coord_cartesian(ylim=c(0, 1700)) +
  scale_y_continuous(label=scales::comma, breaks = seq(0, 1700, 100)) +
  scale_x_continuous(breaks=seq(1975, 2025, 5)) +
  theme(legend.title = element_blank(),
        legend.position = c(.95, .95),                #Position the legend within the plot====         
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)); unknown

#Combine 1975-2006 to 2007-present ADF&G fish ticket data
Halibut.harv.1975 <- (hal_fishery %>%
  #Filter here to keep directed halibut trips - the B permits were introduced in
  #1977, so I am also going to assume anything that is fishery "unknown" and 
  #the gear is longline is a directed halibut trip
  filter(Year %in% c("1975","1976") & fishery == "Unknown" & gear == "Longline" |
           fishery.code %in% c("B06B", "B61B", "B05B", "B99B", "B26B", 
                               "B25B", "B61Z", "B06Z", "B09B", "B91B")) %>% 
  group_by(Year, Mgt.Area) %>%
  summarise(HA.lbs = sum(ha.lbs, na.rm = TRUE),HA.mt  = sum(ha.mt,  na.rm = TRUE),.groups = "drop"))


lapply(Halibut.harv.1975[c("Year", "Mgt.Area")], unique)

#SAVE this data for use in estimating historical bycatch

write.csv(Halibut.harv.1975, paste0("Data_processing/Data/SE_Halibut_removals_",
                                    min(Halibut.harv.1975$Year), "-",
                                    max(Halibut.harv.1975$Year), ".csv"), row.names = FALSE)


################################################################################
### HALIBUT HAREVEST RECONSTRUCTION ###
################################################################################

## STEP 1 ######################################################################
# get the contemporary data from IPHC request, get SEO data from there: HA.req
#         a) 2C partitioned into SEO and SEI - measure this PROPORTION (2Cprops)
#         b) 3A to SEO from FISHTICKET data proportion - measure and apply (3A prop)

## Stat areas are further subdivided into "inside/outside" or I/O areas when it 
## was found that outside waters has a greater proportion of older fish an inside waters.
IOs<-unique(HA.req %>% filter(IPHC.Regulatory.Area == "2C" | IPHC.Regulatory.Area == "3A") %>%
              select(IPHC.Statistical.Area,IPHC.Stat.Area..1929.1975,IPHC.Region.2..1929.1975))

HA.newreq <- HA.newreq %>% 
  mutate(IPHC.Regulatory.Area = str_trim(IPHC.Regulatory.Area)) %>% #removes the space after the letters in iphc regulatory area
  filter(IPHC.Regulatory.Area == "2C" | IPHC.Regulatory.Area == "3A") %>% 
  #Net weight (lb): head off, eviscerated, ice and slime deducted weight
  mutate(Net_lbs=Net.wt..lb./1000) %>%
  select(Year,IPHC.Regulatory.Area,IPHC.Statistical.Area = IPHC..Statistical.Area,Net_lbs)

HA.newreq <- join(HA.newreq,IOs,by="IPHC.Statistical.Area")

head(HA.req)
head(HA.newreq)

HA.req<-rbind(HA.req %>% 
                select(Year,
                       IPHC.Regulatory.Area,
                       IPHC.Statistical.Area,
                       Net_lbs = Pacific.halibut.Net.wt..000.lbs.,
                       IPHC.Stat.Area..1929.1975,
                       IPHC.Region.2..1929.1975),
              HA.newreq) %>% 
  mutate(Net_lbs = as.numeric(Net_lbs)) %>% 
  arrange(Year)

HA.req_2C_3A <- HA.req %>%
  filter(IPHC.Regulatory.Area %in% c("2C", "3A")) %>%
  group_by(IPHC.Regulatory.Area, Year, IPHC.Region.2..1929.1975) %>%
  summarise(Halibut_lbs = sum(Net_lbs)) %>%
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
  
  Hal.SPM[i,"3A.harv"]<-sum(req$Halibut_mt[req$IPHC.regarea == "3A"]) # Total 3A halibut harvest (IPHC data)
  
  Hal.SPM[i,"3Ayak.harv"]<-sum(req$Halibut_mt[req$IPHC.regarea == "3A" & 
                                                req$IO == "Yakutat"]) # Portion of 3A harvest attributed to Yakutat (IPHC data)
  
  Hal.SPM[i,"2C.harv"]<-sum(req$Halibut_mt[req$IPHC.regarea == "2C"]) # Total 2C harvest (IPHC data)
  
  Hal.SPM[i,"EYKT"]<-subha$HA.mt[subha$Mgt.Area == "EYKT"] # Total EYKT harvest (ADFG FT data)
  
  Hal.SPM[i,"SEO2C.tix"]<-sum(subha$HA.mt[subha$Mgt.Area == "NSEO" |
                                            subha$Mgt.Area == "CSEO" |
                                            subha$Mgt.Area == "SSEO"]) # Total OUTSIDE harvest (ADFG FT data)
  
  Hal.SPM[i,"SEI2C.tix"]<-sum(subha$HA.mt[subha$Mgt.Area == "NSEI" |
                                            subha$Mgt.Area == "SSEI"], na.rm=T) # Total INSIDE harvest (ADFG FT data)
  
  Hal.SPM[i,"SEO2C.req"]<-req$Halibut_mt[req$IO == "SE-O"] # Total OUTSIDE harvest (IPHC data)
  
  Hal.SPM[i,"SEI2C.req"]<-req$Halibut_mt[req$IO == "SE-I"] # Total INSIDE harvest (IPHC data)

  
  Hal.SPM[i,"prop3A.EYKT"]<-Hal.SPM[i,"EYKT"]/Hal.SPM[i,"3A.harv"] # Proportion of 3A harvest attributed to EYKT (from ADFG FT data vs IPHC 3A total).
  Hal.SPM[i,"prop3A.EYKT_yak"]<-Hal.SPM[i,"EYKT"]/Hal.SPM[i,"3Ayak.harv"] # Proportion of Yakutat portion of 3A that is EYKT.
  Hal.SPM[i,"prop2C.SEOtix"]<-Hal.SPM[i,"SEO2C.tix"]/Hal.SPM[i,"2C.harv"] # Proportion of 2C harvest that is CSEO, NSEO, SSEO, using ADFG FT data.
  Hal.SPM[i,"prop2C.SEItix"]<-Hal.SPM[i,"SEI2C.tix"]/Hal.SPM[i,"2C.harv"] # Proportion of 2C harvest that is SEI, using ADFG FT  data.
  Hal.SPM[i,"prop2C.SEOreq"]<-Hal.SPM[i,"SEO2C.req"]/Hal.SPM[i,"2C.harv"] # Proportion of 2C harvest that is SEO, from IPHC data.
  Hal.SPM[i,"prop2C.SEIreq"]<-Hal.SPM[i,"SEI2C.req"]/Hal.SPM[i,"2C.harv"] # Proportion of 2C harvest that is SEI, from IPHC data.
  
  Hal.SPM[i,"SEO_harvest_mt"]<- Hal.SPM[i,"SEO2C.req"]+Hal.SPM[i,"EYKT"] 
  #Combines SE-O (from IPHC data) and EYKT (from ADFG FT data) to give total SEO harvest for that year.
  
  i<-i+1
}

#halibut harvest in SEO 2C
plot(data=Hal.SPM, SEO2C.req~Year, ylim=c(0,3500), type="l")

#halibut harvest in SEO 2C from ADF&G FTs
lines(data=Hal.SPM,SEO2C.tix~Year,type="l",col="blue")

data_check_plot <- Hal.SPM %>%
  select(Year, SEO2C.req, SEO2C.tix) %>%
  pivot_longer(cols = c(SEO2C.req, SEO2C.tix), 
               names_to = "Source", 
               values_to = "Harvest") %>% 
  ggplot(aes(x = Year, y = Harvest, color = Source)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("SEO2C.req" = "black", "SEO2C.tix" = "blue")) +
  ylim(0, 3500) +
  labs(
    title = "Halibut Harvest in SEO 2C",
    y = "Harvest (mt)", 
    x = "Year", 
    color = "Data Source");data_check_plot



summary <- Hal.SPM %>%
  summarise(across(prop3A.EYKT:prop2C.SEIreq, list(min = min, max = max))); summary

## For hindcasting proportions lets use pre-full retention
## Full retention was required for all DSR captured in groundfish and halibut 
## fisheries in federal waters starting in 2005 Full retention of all DSR 
## captured in groundfish and halibut fisheries in state waters started in 2009.

Pre2010<-Hal.SPM[Hal.SPM$Year < 2010,]

Prop.EYKT_3A<-mean(Pre2010$prop3A.EYKT) #average proportion of Area 3A harvest that came from EYKT
var.EYKT_3A<-var(Pre2010$prop3A.EYKT)

Prop.EYKT_yak<-mean(Pre2010$prop3A.EYKT_yak) #average proportion of Yakutat (3A) harvest that was from EYKT
var.EYKT_yak<-var(Pre2010$prop3A.EYKT_yak)

Prop.SEO_2C<-mean(Pre2010$prop2C.SEOreq) #average proportion of Area 2C harvest that came from SEO
var.SEO_2C<-var(Pre2010$prop2C.SEOreq)

Halibut_harvest_forSPM <- Hal.SPM %>% 
  select(Year,SEO_harvest_mt)
Halibut_harvest_forSPM$var<-0

#STEP 2 ########################################################################
## Data from 1929-1975, get SEO data plus plus small proportion of 3A

yrs29<-unique(HA.29_75$Year); yrs29 #1929-1975

prep <- HA.29_75 %>% 
  filter(IPHC.Region.2 %in% c("SE-O", "SE-I", "Yakutat")) %>%
  group_by(Year,IPHC.Region.2) %>% 
  summarize(harv.ustons = sum(Total.Catch..US.Tons..Round.Weight.)) %>%
  mutate(harv_mt = harv.ustons*0.90718474) #converting US tons to metric tons

out29<-data.frame()
yrs29<-unique(HA.29_75$Year)
i<-1
for (y in yrs29){
  dat<-prep %>% filter(Year == y)
  out29[i,"Year"]<-y
  out29[i,"SEO_harvest_mt"]<-dat$harv_mt[dat$IPHC.Region.2 == "SE-O"] +
    dat$harv_mt[dat$IPHC.Region.2 == "Yakutat"]*Prop.EYKT_yak #Proportion of Yakutat portion of 3A that is EYKT
  out29[i,"var"]<-var.EYKT_yak*dat$harv_mt[dat$IPHC.Region.2 == "Yakutat"]*dat$harv_mt[dat$IPHC.Region.2 == "Yakutat"]
  i<-i+1
}

Halibut_harvest_forSPM<-rbind(Halibut_harvest_forSPM,out29)
Halibut_harvest_forSPM<-Halibut_harvest_forSPM [order(Halibut_harvest_forSPM$Year), ]

## STEP 3 ######################################################################
## Data from 1975-1982, get SEO data plus plus small proportion of 3A
##
## Note from Phil: Ha.75_82 OK for SEO but missing Yakutat numbers, so will use 
## HA.IPHC web to est and add

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
## Data from 1888-1928 data, partition using proportions: HA.IPHCweb
##         a) apply 2Cprop and 3Aprop

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
  scale_x_continuous(breaks=seq(1880,2025,5)) + 
  theme (panel.grid.minor = element_blank(),
         axis.text.x = element_text(angle = 45, vjust=1, hjust=1))

ggsave(paste0("Figures/Halibut_Harvest_SEO_", YEAR, ".png"), dpi=300,  height=3, width=6, units="in")

write.csv(Halibut_harvest_forSPM,paste0("Data_processing/Data/SEO_Halibut_removals_",
                                        min(Halibut_harvest_forSPM$Year),"-",
                                        max(Halibut_harvest_forSPM$Year),".csv"))
min(Halibut_harvest_forSPM$Year)
max(Halibut_harvest_forSPM$Year)


#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
## Data check stuff here
## LSC - I want to come back and investigate this more. The proportion for SEI is
## greater than 1, which doesn't seem correct.

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

overlap.HA<-old.HA %>% filter(Year > 1976)

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

