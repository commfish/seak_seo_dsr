#################################################################################################
## SEO Yelloweye Rockfish NOAA Observer Data for Pop Models
## Phil Joy, Fall 2021
## Data used in production models
## Extract both BYCATCH and DISCARD data 
## 1) Create 1 file for SEO
## 2) apply Harvests of halibut in subdistricts to estimate discard and bycatch by subdistrict
#################################################################################################

library(dplyr)
library(boot)
library(ggplot2)

wd="C:/Users/pjjoy/Documents/Groundfish Biometrics/Yelloweye_Production_Models/"
setwd(wd)
getwd()

Obs<-read.csv("yelloweye-observer_dat/yelloweye/data/Groundfish Total Discards by Species.csv", header=T)
str(Obs)
Obs<-rename(Obs, Year = ï..Year)

unique(Obs$FMP.Area)
unique(Obs$FMP.Subarea)
unique(Obs$NMFS.Area)
unique(Obs$Gear)
unique(Obs$Target)
unique(Obs$Species.Code)

Years<-unique(Obs$Year)

#=================================================================================================
# 1) 
##Bootstrap mean bycatch and discard rate each year in SEO waters
OUT<-data.frame()

nboot<-1000
i<-1
for (y in Years){   #y<-Years[1]
  Dat<-Obs[Obs$Year == y & Obs$FMP.Subarea == "SE",]
  YE<-Dat[Dat$Species.Code==145,]
  HA<-Dat[Dat$Species.Code==200,]
  OUT[i,"Year"]<-y
  OUT[i,"YE.Discards"]<-sum(YE$Discards..mt.)
  OUT[i,"YE.Bycatch"]<-sum(YE$Total.Catch..mt.)
  OUT[i,"YE.Bycatch.rate"]<-sum(YE$Total.Catch..mt.)/sum(HA$Total.Catch..mt.)
  OUT[i,"YE.Discard.rate"]<-OUT[i,"YE.Discards"]/sum(YE$Total.Catch..mt.)  #prop of YE bycatch that is discarded during full retention
  
  boot<-data.frame()
  for (n in 1:nboot){  #n<-1
    YE.boot <- YE[sample(1:nrow(YE), nrow(YE), replace=TRUE), ]
    HA.boot <- HA[sample(1:nrow(YE), nrow(YE), replace=TRUE), ]
    boot[n,"YE.Discards"]<-sum(YE.boot$Discards..mt.)
    boot[n,"YE.Bycatch"]<-sum(YE.boot$Total.Catch..mt.)
    boot[n,"YE.Bycatch.rate"]<-sum(YE.boot$Total.Catch..mt.)/sum(HA.boot$Total.Catch..mt.)
    boot[n,"YE.Discard.rate"]<-boot[n,"YE.Discards"]/sum(YE.boot$Total.Catch..mt.)
  }
  
  OUT[i,"YE.Discards.boot"]<-median(boot$YE.Discards)
  OUT[i,"YE.Bycatch.boot"]<-median(boot$YE.Bycatch)
  OUT[i,"YE.Bycatch.rate.boot"]<-median(boot$YE.Bycatch.rate)
  OUT[i,"YE.Byrate.boot.var"]<-var(boot$YE.Bycatch.rate)
  OUT[i,"YE.Byrate.boot.cv"]<-sqrt(OUT[i,"YE.Byrate.boot.var"])/OUT[i,"YE.Bycatch.rate.boot"]
  OUT[i,"YE.Discard.rate.boot"]<-median(boot$YE.Discard.rate)
  OUT[i,"YE.Disrate.boot.var"]<-var(boot$YE.Discard.rate)
  OUT[i,"YE.Disrate.boot.cv"]<-sqrt(OUT[i,"YE.Disrate.boot.var"])/OUT[i,"YE.Discard.rate.boot"]
  
  i<-i+1
}

#=============================================================================================
# Save file for modelling
write.csv(OUT,"Data/SEO_YE_NOAA_Observers.csv")

#=============================================================================================
# 2) Get bycatch and discard numbers by subdistrict by applying SEO/SEI data based on halibut harvests in those subdistricts
#load Halibut harvests as YE discard index
Hali<-read.csv("Data/halibut_catch_data.csv", header=T)
Hali<-Hali[Hali$gear.description == "Longline",]
Hali<-rename(Hali, SEdist = groundfish.mgt.area.district)
Hali<-rename(Hali, Year = year.landed)

#plot(Hali$Year,Hali$round_lbs_hali)
SEsubd.Hali<-data.frame()
Yrs<-unique(Hali$Year)
SDs<-unique(Hali$SEdist)

i<-1
for (y in Yrs){  #y<-Yrs[1]
  for (s in SDs){ #s<-SDs[1]
    Dat<-Hali[Hali$Year == y & Hali$SEdist == s,]
    SEsubd.Hali[i,"Year"]<-y; SEsubd.Hali[i,"Subdistrict"]<-s
    SEsubd.Hali[i,"round_lbs_hali"]<-sum(Dat$round_lbs_hali); i<-i+1
  }
}
SEsubd.Hali$round_mt_hali<-SEsubd.Hali$round_lbs_hali*0.45359237/1000
SEsubd.Hali[47,"Year"]<-2021; SEsubd.Hali[47,"round_mt_hali"]<-mean(SEsubd.Hali[42:46,"round_mt_hali"])
plot(SEsubd.Hali$round_mt_hali)

#might as well save halibut data for model here...
write.csv(SEsubd.Hali,"Data/SEsubd_Halibut_harvs.csv")

#Just do SSO right now
SEO.Hal<-SEsubd.Hali[SEsubd.Hali$Subdistrict == "NSEO" |
                       SEsubd.Hali$Subdistrict == "CSEO" |
                       SEsubd.Hali$Subdistrict == "SSEO" |
                       SEsubd.Hali$Subdistrict == "EYKT",]

#Load IPHC survey data
IPHC<-read.csv("Data/SEsubd_YE_IPHCbycatch.csv", header=T)
IPHC<-IPHC[IPHC$SEdist == "NSEO" |
             IPHC$SEdist == "CSEO" |
             IPHC$SEdist == "SSEO" |
             IPHC$SEdist == "EYKT",]

Obs.subd<-data.frame()
Years<-unique(OUT$Year)
SD<-unique(SEO.Hal$Subdistrict)

i<-1
for (y in Years){   #y<-Years[9]
  
  if (y < 2021){
    Dat<-OUT[OUT$Year == y,]
    Hal.Dat<-SEO.Hal[SEO.Hal$Year == y,]
    Surv.Dat<-IPHC[IPHC$Year == y,]
  } else {
    Dat<-OUT[OUT$Year == y,]
    Hal.Dat<-SEO.Hal[SEO.Hal$Year >= y-6 & SEO.Hal$Year < y,]
    Surv.Dat<-IPHC[IPHC$Year >= y-6 & IPHC$Year < y,]
  }
  
  for (s in SD){  #s<-SD[4]
    Hal.Dat2<-Hal.Dat[Hal.Dat$Subdistrict == s,]
    Surv.Dat2<-Surv.Dat[Surv.Dat$SEdist == s,]
    
    Obs.subd[i,"Year"]<-y
    Obs.subd[i,"Subdistrict"]<-s
    
    if (y < 2021){
      Hal.rat<-(Hal.Dat2$round_mt_hali/sum(Hal.Dat$round_mt_hali)) #prop halibut in subd
      Surv.rat<-(Surv.Dat2$YE_to_HA_ratio /sum(Surv.Dat$YE_to_HA_ratio)) #
      Surv.rat.TRUE<-Surv.Dat2$YE_to_HA_ratio
    } else {
      Hal.rat<-(sum(Hal.Dat2$round_mt_hali)/sum(Hal.Dat$round_mt_hali)) #prop halibut in subd
      Surv.rat<-(sum(Surv.Dat2$YE_to_HA_ratio) /sum(Surv.Dat$YE_to_HA_ratio)) #
      Surv.rat.TRUE<-mean(Surv.Dat2$YE_to_HA_ratio)
    }
    
    #Hal.rat<-(Hal.Dat2$round_mt_hali/sum(Hal.Dat$round_mt_hali)) #prop halibut in subd
    #Surv.rat<-(Surv.Dat2$YE_to_HA_ratio /sum(Surv.Dat$YE_to_HA_ratio)) #
    #Surv.rat.TRUE<-Surv.Dat2$YE_to_HA_ratio
    
    Obs.subd[i,"Hal.rate"]<-Hal.rat
    Obs.subd[i,"Surv.rate"]<-Surv.rat
    Obs.subd[i,"Survey.By"]<-Surv.rat.TRUE
    Obs.subd[i,"Adj.Surv"]<-Surv.rat.TRUE/mean(Surv.Dat$YE_to_HA_ratio)
    
    Obs.subd[i,"YE.Bycheck"]<-Dat$YE.Bycatch
    Obs.subd[i,"YE.Bycatch_hal"]<-Dat$YE.Bycatch*Hal.rat
    Obs.subd[i,"YE.Bycatch_surv"]<-Dat$YE.Bycatch*Surv.rat
    Obs.subd[i,"YE.By_schwag"]<-(4*Obs.subd[i,"YE.Bycatch_hal"]+Obs.subd[i,"YE.Bycatch_surv"])/5
    
    Obs.subd[i,"YE.Discheck"]<-Dat$YE.Discards
    Obs.subd[i,"Di_schwag"]<-(4*(Dat$YE.Discards*Hal.rat)+(Dat$YE.Discards*Surv.rat))/5
    
    i<-i+1
  }
}

#================================================================================================
# Save schwaggy observer data for subdistrict level production models

write.csv(Obs.subd,"Data/SEsubd_YE_NOAA_Observers.csv")

#=================================================================================================

head(Obs.subd)
sum(Obs.subd$YE.Bycatch1)
sum(Obs.subd$YE.Bycatch2)
sum(Obs.subd$YE.By_schwag2)
mean(Obs.subd$Survey.By)
##################################################################################################
### SCRAP
head(YE)
head(YE.boot)

head(Dat)
unique(Dat$FMP.Subarea)

tot<-sum(Dat$Total.Catch..mt.)
YE.Dat<-Dat[Dat$Species.Code == 145,]
HA.Dat<-Dat[Dat$Species.Code == 200,]

#total catches of halibut and yelloweye
tot.YE<-sum(YE.Dat$Total.Catch..mt.)
tot.HA<-sum(HA.Dat$Total.Catch..mt.)
tot.YE/tot.HA

#Discard of YE
disc<-sum(YE.Dat$Discards..mt.)
ret<-sum(YE.Dat$Retained..mt.)
tot.YE<-sum(YE.Dat$Total.Catch..mt.)
disc/tot.YE
ret/tot.YE

data_s1 <- data[sample(1:nrow(data), 3), ]





