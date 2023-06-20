#################################################################################################
## SEO Yelloweye Rockfish Harvest Reconstruction for Pop Models
## Phil Joy, July 2022
## Data used in production models
#################################################################################################

library(dplyr)
library(boot)
library(ggplot2)
library(scales)

YEAR<-2022

YE.Harv<-read.csv("Data/Harvests/yelloweye_catch_data.csv", header=T)
HA.Harv<-read.csv("Data/halibut_catch_data_new071422.csv", header=T)
YE.Subs<-read.csv("Data/Harvests/subsistence.csv", header=T)
YE.Sport<-read.csv("Data/Harvests/sport.csv", header=T)

## Recent data updated
YE.Harv.update1<-read.csv("Data/Yelloweye harvest SEO in groundfish fish ticket data 2015-current.csv")
YE.Harv.update2<-read.csv("Data/Yelloweye harvest SEO in salmon fish ticket data 2015-current.csv")

#Halibut fish ticket data:
HA.Harv.update<-read.csv("Data/Halibut harvest from fish ticket data_071522.csv")
#Halibut by IPHC area from web source data
HA.IPHCweb<-read.csv("Data/Halibut_harvests_IPHCareas_1888.csv", skip=1, header=T)
#Halibut harvest from IPHC data request 1982 - present
HA.IPHCreq<-read.csv("Data/Halibut_Harvest_IPHCdatareq_1982_2022.csv")

str(HA.Harv.IPHC)
#********************************************************************************
## Recent data updated
str(YE.Harv)
str(YE.Harv.update1)
YE.Harv.update1$Year<-YE.Harv.update1[,1]
str(YE.Harv.update2)
YE.Harv.update2$Year<-YE.Harv.update2[,1]

YE.Harv.update<-rbind(YE.Harv.update1,YE.Harv.update2)

unique(YE.Harv.update$Species.Code)
colnames(YE.Harv.update)

YE.Harv.update <- YE.Harv.update %>% 
  mutate(Whole.weight.lbs = Whole.Weight..sum.,
         Mgt.Area = Mgt.Area.District) 

#*****************************************************************************
#------------------------------------------------------------------------------
# compare old and new data if you want.  C
{str(YE.Harv.update)
YE.Harv.update %>%
  group_by(Year,Mgt.Area,CFEC.Fishery.Code) %>% 
  summarize(total.removals.lbs = sum(Whole.weight.lbs)) %>%
  mutate(total.removals.mt = 0.00045359*total.removals.lbs,
         source = "update",
         Fishery.Code = CFEC.Fishery.Code,
         gear=Gear.Name )%>%
  select(Year,Mgt.Area,source,Fishery.Code,gear,total.removals.lbs,total.removals.mt) -> Comm.harv.up
str(Comm.harv.up)
#
YE.Harv %>%
  group_by(year,mgmt_area,fishery) %>% 
  summarize(total.removals.lbs = sum(whole_pounds)) %>%
  mutate(total.removals.mt = 0.00045359*total.removals.lbs,
         Mgt.Area = mgmt_area, 
         source = "original",
         Fishery.Code = fishery) %>%
  select(Year=year,Mgt.Area,source,Fishery.Code,gear,total.removals.lbs,total.removals.mt)-> Comm.harv.orig
str(Comm.harv.orig)

harvests<-rbind(Comm.harv.up,Comm.harv.orig)

ggplot(harvests,aes(Year,total.removals.mt,col = source)) +
  geom_line() + 
  facet_wrap(~Mgt.Area)

#old data and update data matches up for commercial harvests... yeah! 
# going forward set up for adding fish ticket data
YE.comm<-rbind(Comm.harv.up,Comm.harv.orig %>% 
                 filter(Year < min(Comm.harv.up$Year)))}
#--------------------------------------------------------------------------------

str(YE.Harv)
unique(YE.Harv$year)
str(YE.Harv.update)

with(YE.Harv, table(fishery,fishery_type))

YE.Harv %>% mutate(Year = year, Mgt.Area = mgmt_area, ye.lbs = whole_pounds,
                   fishery.code = fishery,
                   ha.bycatch = ifelse(substr(fishery,1,1) %in% c("B"),"hal.by","other"),
                   ye.mt = ye.lbs*0.00045359) %>%
  dplyr::select(Year,Mgt.Area,fishery.code,gear,ye.lbs,ye.mt,ha.bycatch) ->YE.Harv.old

YE.Harv.update %>% mutate(Mgt.Area = Mgt.Area.District,
                          gear = Gear.Name,
                          ye.lbs = Whole.weight.lbs,
                          ye.mt = ye.lbs*0.00045359,
                          fishery.code = CFEC.Fishery.Code,
                          ha.bycatch = ifelse(substr(fishery.code,1,1) %in% c("B"),"hal.by","other")) %>% 
  dplyr::select(Year,Mgt.Area,fishery.code,gear,ye.lbs,ye.mt,ha.bycatch) ->YE.Harv.latest

YE.comm.rem<-rbind(YE.Harv.old %>% filter(Year < min(YE.Harv.latest$Year)),
                   YE.Harv.latest)
ggplot(YE.comm.rem,aes(Year,ye.mt,col = ha.bycatch)) +
  geom_line() + 
  facet_wrap(~Mgt.Area)

max(YE.comm.rem$Year)

write.csv(YE.comm.rem,paste0("Data/SE_YE_comm_removals",min(YE.comm.rem$Year),"-",
                             max(YE.comm.rem$Year),".csv"))

#**********************************************************************************
# add in sport fish harvests... 
str(YE.Sport); str(YE.comm.rem)
unique(YE.Sport$species)
YE.Sport<-YE.Sport[YE.Sport$species == 145,]

ex<-YE.Sport[YE.Sport$mgmt_area == "CSEO",]


YE.Sport %>% mutate(Mgt.Area = mgmt_area, ye.lbs=total_lbs,
                    ye.mt = ye.lbs*0.00045359,
                    fishery.code = "sport",
                    gear = "sport",
                    ha.bycatch = "other") %>%
  dplyr::select(Year,Mgt.Area,fishery.code,gear,ye.lbs,ye.mt,ha.bycatch) ->YE.Sport1


mus<-unique(YE.Sport1$Mgt.Area)
maxyr<-max(YE.Sport1$Year)
nymis<-YEAR-maxyr

rec.avg<-data.frame()
if(YEAR > maxyr) {
  for (g in mus){  #g<-mus[1]
    for (i in nymis){  #i<-1
      dat<-YE.Sport1 %>% filter(Mgt.Area == g)
      ix<-i+which(mus == g)-1
      
      rec.avg[ix,"Year"]<-maxyr+i
      rec.avg[ix,"Mgt.Area"]<-g
      rec.avg[ix,"fishery.code"]<-"sport"
      rec.avg[ix,"gear"]<-"sport"
      rec.avg[ix,"ye.lbs"]<-NA
      rec.avg[ix,"ye.mt"]<-mean(dat$ye.mt[dat$Year>=(maxyr-1) & dat$Year <= (maxyr)])
      rec.avg[ix,"ha.bycatch"]<-NA
    }
  }
}

nrow(YE.Sport1)
YE.Sport2<-rbind(YE.Sport1,rec.avg)
nrow(YE.Sport2)
YE.removals<-rbind(YE.comm.rem,YE.Sport2) 

#********************************************************************************
#*add in subsistence
str(YE.Subs)
range(YE.Subs$Year)

ys<-unique(YE.Subs$Year)
mgmtas<-unique(YE.removals$Mgt.Area)

YE.sub<-data.frame()
i<-1
for (y in ys){     #y<-ys[1]
  for (a in mgmtas) {  #a<-mgmtas[1]
    subh<-YE.Subs[YE.Subs$Year == y,]
    YE.sub[i,"Year"]<-y
    YE.sub[i,"Mgt.Area"]<-a
    YE.sub[i,"sd.prop"]<-ifelse(YE.sub[i,"Mgt.Area"] %in% c("SSEO"), 1056/3898, 
                                ifelse(YE.sub[i,"Mgt.Area"] %in% c("CSEO"), 1661/3898, 
                                       ifelse(YE.sub[i,"Mgt.Area"] %in% c("NSEO"), 442/3898, 
                                              ifelse(YE.sub[i,"Mgt.Area"] %in% c("EYKT"), 739/3898, 
                                                     ifelse(YE.sub[i,"Mgt.Area"] %in% c("NSEI"), 0.6, 
                                                            ifelse(YE.sub[i,"Mgt.Area"] %in% c("SSEI"), 0.4, NA))))))
    if (YE.sub[i,"Mgt.Area"] == "NSEI" | YE.sub[i,"Mgt.Area"] == "SSEI") {
      YE.sub[i,"ye.lbs"] <-subh$SEI*YE.sub[i,"sd.prop"]
    } else {
      YE.sub[i,"ye.lbs"] <-subh$SEO*YE.sub[i,"sd.prop"]
    }
    i<-i+1
  }
}

YE.sub %>% mutate(ye.mt = ye.lbs*0.00045359,
                  fishery.code = "subsistence",
                  gear = "subsistence",
                  ha.bycatch = "other") %>%
  dplyr::select(Year,Mgt.Area,fishery.code,gear,ye.lbs,ye.mt,ha.bycatch) -> YE.sub1

head(YE.sub1)

mus<-unique(YE.sub1$Mgt.Area)
maxyr<-max(YE.sub1$Year)
nymis<-YEAR-maxyr

rec.avg<-data.frame()
if(YEAR > maxyr) {
  for (g in mus){  #g<-mus[1]
    for (i in nymis){  #i<-1
      dat<-YE.sub1 %>% filter(Mgt.Area == g)
      ix<-i+which(mus == g)-1
      
      rec.avg[ix,"Year"]<-maxyr+i
      rec.avg[ix,"Mgt.Area"]<-g
      rec.avg[ix,"fishery.code"]<-"sub"
      rec.avg[ix,"gear"]<-"sub"
      rec.avg[ix,"ye.lbs"]<-NA
      rec.avg[ix,"ye.mt"]<-mean(dat$ye.mt[dat$Year>=(maxyr-1) & dat$Year <= (maxyr)])
      rec.avg[ix,"ha.bycatch"]<-NA
    }
  }
}

nrow(YE.sub1)
YE.sub2<-rbind(YE.sub1,rec.avg)
nrow(YE.sub2)

YE.removals<-rbind(YE.removals,YE.sub2) 

#******************************************************************************
#*put together removal data for use in SS-SPM
str(YE.removals)

YE.removals %>% group_by(Year,Mgt.Area) %>%
  summarise(tot.rem.mt = sum(ye.mt)) ->pt1

YE.removals %>% filter(ha.bycatch == "hal.by") %>% 
  group_by(Year,Mgt.Area) %>%
  summarise(tot.hal.by = sum(ye.mt)) ->pt2

?join
All.known.YE.removals<-full_join(pt1,pt2)

#******************************************************************************
#*SAVE REMOVAL DATA for SS-SPM
write.csv(All.known.YE.removals,paste0("Data/SE_YE_known_removals_",min(All.known.YE.removals$Year),"-",
                             max(All.known.YE.removals$Year),".csv"))

#******************************************************************************
#*****************--------**************************************************
#**************---------------***************************************************
#*********************************************************************************
#*#HALIBUT removals for bycatch estimation
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
                          ha.lbs = Whole.Weight..sum.,
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

write.csv(Halibut.harv.1975,paste0("Data/SE_Halibut_removals_",min(Halibut.harv.1975$Year),"-",
                                       max(Halibut.harv.1975$Year),".csv"))

#**********************************************************************************
#*Get SEO estimates of Halibut removals pre-1975 from RP compiled data and 
#Step 1) get the contemporary data from IPHC request, get SEO data from there: HA.IPHCreq
#         a) 2C partitioned into SEO and SEI - measure this PROPORTION (2Cprops)
#         b) 3A to SEO from FISHTICKET data proportion - measure and apply (3A prop)
# step 2) get 1929-1981 data, get SEO data plus plus small proportion of 3A
# step 3) get 1888-1928 data, partition using proportions: HA.IPHCweb
#         a) apply 2Cprop and 3Aprop
HA.Harv<-read.csv("Data/halibut_catch_data_new071422.csv", header=T)
#Halibut fish ticket data:
#HA.fishtix<-read.csv("Data/Halibut harvest from fish ticket data_071522.csv")
HA.Harv.update<-read.csv("Data/Halibut harvest from fish ticket data_071522.csv") 

Halibut.harv.1975

#Halibut harvest from IPHC data request 1982 - present
HA.req<-read.csv("Data/Halibut_Harvest_IPHCdatareq_1982_2022.csv")
#Halibut harvest from IPHC 1929-1975; IPHC Scientific report 67
HA.29_75<-read.csv("Data/IPHC_Halibut_harv_1929-1975.csv")
#Halibut harvest IPHC 1975-1982
Ha.75_82<-read.csv("Data/IPHC_Halibut_harv_1975-1982.csv")
#Halibut by IPHC area from web source data
HA.web<-read.csv("Data/Halibut_harvests_IPHCareas_1888.csv", skip=1, header=T)


#-------------------------------------------------------------------------------
# 
HA.Harv %>% mutate(Year = year.landed, Mgt.Area = mgmt.area, ha.lbs = round.lbs,
                   fishery.code = permit.fishery,
                   gear = gear.description,
                   ha.mt = ha.lbs*0.00045359,
                   source = "old") %>%
  dplyr::select(Year,Mgt.Area,fishery.code,gear,ha.lbs,ha.mt,source) ->HA.Harv.old

HA.Harv.update %>% mutate(Year=DOL.Year,
                          gear = Gear.Code.and.Name,
                          ha.lbs = Whole.Weight..sum.,
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

ch98<-HA.Harv[HA.Harv$year.landed == 1998,]

write.csv(Halibut.harv.1975,paste0("Data/SE_Halibut_removals_",min(Halibut.harv.1975$Year),"-",
                                   max(Halibut.harv.1975$Year),".csv"))

#STEP 1------------------------------------------------------------------------
HA.req %>% filter(IPHC.Regulatory.Area == "2C" | IPHC.Regulatory.Area == "3A") %>%
  group_by(IPHC.Regulatory.Area,Year,IPHC.Region.2..1929.1975) %>%
  summarise(Halibut_lbs = sum(Pacific.halibut.Net.wt..000.lbs.)) %>%
  mutate(Halibut_mt = Halibut_lbs*0.00045359*1000,
         IO = IPHC.Region.2..1929.1975) %>% 
  dplyr::select(Year,IPHC.regarea=IPHC.Regulatory.Area,IO,Halibut_lbs,
                Halibut_mt)-> HA.req_2C_3A

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

plot(data=Hal.SPM,SEO2C.req~Year, ylim=c(0,3000), type="l")
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
str(HA.29_75)
head(Hal.SPM)
yrs29<-unique(HA.29$Year)

HA.29_75 %>% filter(IPHC.Region.2 == "SE-O" |
                   IPHC.Region.2 == "SE-I" |
                   IPHC.Region.2 == "Yakutat") %>%
  group_by(Year,IPHC.Region.2) %>% 
  summarize(harv.ustons = sum(Total.Catch..US.Tons..Round.Weight.)) %>%
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

write.csv(Halibut_harvest_forSPM,paste0("Data/SEO_Halibut_removals_",
                                        min(Halibut_harvest_forSPM$Year),"-",
                                   max(Halibut_harvest_forSPM$Year),".csv"))

#--------------------------------------------------------------------------------
# plot Foreign removals nicer... 
Foreign<-read.csv("Data/Foreign_YERF_SEAK.csv", skip=1)
str(Foreign)

ggplot(Foreign,aes(Year,Estimate,col = "black")) +
  geom_ribbon(aes(ymin=Estimate-0.75*Estimate, ymax=Estimate+0.75*Estimate),
              alpha=0.5, col="grey") +
  geom_line(col="black") +
  xlab("\nYear") +
  ylab("Yelloweye removals by foreign fleet (t)") +
  #scale_y_continuous(label=comma, breaks = seq(0,6000,1000)) +
  scale_x_continuous(breaks=seq(1960,1985,5)) + 
  theme (panel.grid.minor = element_blank(),
         axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
         axis.title.y = element_text(size=9))

ggsave(paste0("Figures/YE_Harvest_Foreign_", YEAR, ".png"), dpi=300,  height=3, width=4, units="in")
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
  summarise(Halibut_lbs = sum(Pacific.halibut.Net.wt..000.lbs.)) %>%
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


full_join

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

##ID areas by inside/outside

YE.Harv<- YE.Harv %>%
  mutate(IO = ifelse(mgmt_area %in% c("SSEO","CSEO","NSEO","EYKT","SEO"),"SEO", 
                         ifelse(mgmt_area %in% c("NSEI","SSEI"),"SEI",NA))) %>% 
  mutate(sd.prop = ifelse(mgmt_area %in% c("SSEO"), 1056/3898, 
                          ifelse(mgmt_area %in% c("CSEO"), 1661/3898, 
                                 ifelse(mgmt_area %in% c("NSEO"), 442/3898, 
                                        ifelse(mgmt_area %in% c("EYKT"), 739/3898, 
                                               ifelse(mgmt_area %in% c("NSEI"), 0.6, 
                                                      ifelse(mgmt_area %in% c("SSEI"), 0.4, NA)))))))
IO<-unique(YE.Harv$IO)

YE.Sport<- YE.Sport %>%
  mutate(IO = ifelse(mgmt_area %in% c("SSEO","CSEO","NSEO","EYKT","SEO"),"SEO", 
                     ifelse(mgmt_area %in% c("NSEI","SSEI"),"SEI",NA)))

HA.Harv<- HA.Harv %>%
  mutate(IO = ifelse(groundfish.mgt.area.district %in% c("SSEO","CSEO","NSEO","EYKT","SEO"),"SEO", 
                     ifelse(groundfish.mgt.area.district %in% c("NSEI","SSEI"),"SEI",NA))) 

Out.remove<-data.frame()
a<-1

for (i in IO){     #i<-IO[1]
  for (y in Years){  #y<-Years[1]
    Com<-YE.Harv[YE.Harv$year == y & YE.Harv$IO == i,]
    Spt<-YE.Sport[YE.Sport$Year == y & YE.Sport$IO == i,]
    Sub<-YE.Subs[YE.Subs$Year == y,]
    
    Out.remove[a,"Year"]<-y
    Out.remove[a,"IO"]<-i
    Out.remove[a,"Remove.lbs"]<-sum(Com$whole_pounds)+sum(Spt$total_lbs)+
      if (i == "SEO") {Sub$SEO} else {Sub$SEI}
    Out.remove[a,"YE_fr_hal_fish.lbs"]<-sum(Com$whole_pounds[Com$h.code == "B" & Com$gear == "Longline"])
    a<-a+1
  }
}

Out.remove<-Out.remove[order(Out.remove$Year),]

SE.YE.removals<-Out.remove[Out.remove$IO == "SEO",]

#*************************************************************
#SE.YE.removals is the data for the model
#*************************************************************

write.csv(SE.YE.removals,"Data/SEO_YE_removals4.csv")

#=====================================================================================
# Get harvest Data by Subdistrict

Out.remove2<-data.frame()
SD<-unique(YE.Harv$mgmt_area)

a<-1

for (i in SD){     #i<-SD[1]
  for (y in Years){  #y<-Years[1]
    Com<-YE.Harv[YE.Harv$year == y & YE.Harv$mgmt_area == i,]
    Spt<-YE.Sport[YE.Sport$Year == y & YE.Sport$mgmt_area == i,]
    Sub<-YE.Subs[YE.Subs$Year == y,]
    subprop<-unique(Com$sd.prop)
    
    Out.remove2[a,"Year"]<-y
    Out.remove2[a,"Subdistrict"]<-i
    Out.remove2[a,"Remove.lbs"]<-sum(Com$whole_pounds)+sum(Spt$total_lbs)+
      if (i == "SEO") {subprop*Sub$SEO} else {subprop*Sub$SEI}
    Out.remove2[a,"YE_fr_hal_fish.lbs"]<-sum(Com$whole_pounds[Com$h.code == "B" & Com$gear == "Longline"])
    Out.remove2[a,"Remove.mt"]<-Out.remove2[a,"Remove.lbs"]*0.45359237/1000
    Out.remove2[a,"YE_fr_hal_fish.mt"]<-Out.remove2[a,"YE_fr_hal_fish.lbs"]*0.45359237/1000
    a<-a+1
  }
}

Out.remove2<-Out.remove2[order(Out.remove2$Year),]


write.csv(Out.remove2,"Data/SEsubdistrict_YE_removals.csv")


a<-2
b<-10
c<-12

log(2+10+12)
log(2)+log(10)+log(12)
