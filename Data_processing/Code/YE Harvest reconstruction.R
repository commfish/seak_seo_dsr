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

#Pre-2021 YE DATA:
YE.Harv<-read.csv("Data_processing/Data/Harvests/yelloweye_catch_data.csv", header=T)
YE.Subs<-read.csv("Data_processing/Data/Harvests/subsistence.csv", header=T)
YE.Sport<-read.csv("Data_processing/Data/Harvests/sport.csv", header=T)

## Recent data updated
YE.Harv.update1<-read.csv("Data_processing/Data/Harvests/Yelloweye harvest SEO in groundfish fish ticket data 2015-current.csv")
YE.Harv.update2<-read.csv("Data_processing/Data/Harvests/Yelloweye harvest SEO in salmon fish ticket data 2015-current.csv")

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

write.csv(YE.comm.rem,paste0("Data_processing/Data/SE_YE_comm_removals",min(YE.comm.rem$Year),"-",
                             max(YE.comm.rem$Year),".csv"))

##calculate and plot gear breakdown of harvest... 

#**********************************************************************************
# add in sport fish harvests... 
# get the most up to date data here:
# Sport fish estimates usually lag have been a moving target as the division tries
# to settle on a satisfactory method to estimate catch and release mortality
# The current sport fish contact for this assessment is Mike Jaenike (mike.jaenike@alaska.gov)
# Sport fish data usually lags behind by at least one year due to the statewide
# harvest survey.  As such the code will insert the average catch into the newest year(s)
# for which data is missing.

str(YE.Sport); str(YE.comm.rem)
unique(YE.Sport$species)
YE.Sport<-YE.Sport[YE.Sport$species == 145,]

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
#* It has been a while since there were reliable estimates of subsistence harvests.
#* As such, recent years use the long term average

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
#*
write.csv(All.known.YE.removals,paste0("Data_processing/Data/SE_YE_known_removals_",min(All.known.YE.removals$Year),"-",
                                       as.character(Sys.Date()),".csv"))

#******************************************************************************
#*REMOVALS from FOREIGN TRAWL FISHERY from SRI
#* Derived from Donny Olson for SRI , but methods never written up
#--------------------------------------------------------------------------------
# plot Foreign removals nicer... 
Foreign<-read.csv("Data_processing/Data/Harvests/Foreign_YERF_SEAK.csv", skip=1)
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


###########################################################################
### End on 7-6-23... below is still being vetted where to put it... 

#*******************************************************************************
#*****************--------******************************************************
#************** SCRAP ********************************************************
#*******************************************************************************
#*******************************************************************************
#*#HALIBUT removals for bycatch estimation

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
