################################################################################
## Halibut HArvest revisited July 2022
## reconciling fish ticket data from Rhea with Randy Petersen's compilation of the IPHC data
###################################################################################
library(dplyr)
RheaHal<-read.csv("Data/Harvests/halibut_catch_data.csv", header=T)
RheaHal2<-read.csv("Data/halibut_catch_data_new071422.csv")
RPHal<-read.csv("Data/halibut_SEO_catch_RPderivation.csv", header=T)

#*********************************************************************
# Look at Rhea's new data
str(RheaHal2)

#notes: look at Rhea and Randy' data to check for discrepencies... contact Rhea if
# there is something.  
#Make sure Halibut data is formatted the same for wcpue and bycatch estimation 
#function below: H.catch<-Hali[Hali$SEdist == s & Hali$Year == y,]

#Load fishticket data directly...
ys<-c(1975,1982,1988,1993,1996,1998,2001,2003,2005,2007,2009,2011,2013,2015,2017,2019,2021)
for (y in 1:(length(ys)-1)){
  tx<-read.csv(paste0("Data/Harvests/1975-2020/",ys[y],"-",ys[y+1]-1,".csv"),
               stringsAsFactors=T)
  if(y == 1){
    fishtix<-tx
  } else {
    fishtix<-rbind(tx,fishtix)
  }
}
str(fishtix)
fishtix$Year<-fishtix$ï..Year.Landed
min(fishtix$Year);max(fishtix$Year)
unique(fishtix$ADFG.Management.Area.Code)
unique(fishtix$Groundfish.Mgt.Area.District)
unique(fishtix$Permit.Fishery)

#need to get subdistrict data prior to 1991 based on stat area, then compare to
# recent years where they did record ADFG areas... 
SEOtix<-fishtix %>% filter(IPHC.Statistical.Area %in% c(140:144,150:153,160:163,170:174,181:184,190,200)) %>%
  mutate(SEdist_schwag = ifelse(IPHC.Statistical.Area %in% c(182:184,171,173,174,161:163),"NSEI", 
                                ifelse(IPHC.Statistical.Area %in% c(142:144,152,153),"SSEI",
                                       ifelse(IPHC.Statistical.Area %in% c(190),"EYKT",
                                              ifelse(IPHC.Statistical.Area %in% c(140,141,150),"SSEO",
                                                     ifelse(IPHC.Statistical.Area %in% c(151),"NSEI-NSEO-SSEO",
                                                            ifelse(IPHC.Statistical.Area %in% c(160),"CSEO-NSEI",
                                                                   ifelse(IPHC.Statistical.Area %in% c(170),"CSEO-NSEO",
                                                                          ifelse(IPHC.Statistical.Area %in% c(181),"NSEI-NSEO",
                                                                                 ifelse(IPHC.Statistical.Area %in% c(185),"EYKT-NSEO",
                                                                                        ifelse(IPHC.Statistical.Area %in% c(200),"EYKT-IBS","whocares")))))))))))

str(SEOtix)
levels(SEOtix$Groundfish.Mgt.Area.District)<-c(levels(SEOtix$Groundfish.Mgt.Area.District),"no.record")
SEOtix$Groundfish.Mgt.Area.District[SEOtix$Groundfish.Mgt.Area.District == ""] <- "no.record"#"not.recorded"
#SEOtix$Groundfish.Mgt.Area.District[is.na(SEOtix$Groundfish.Mgt.Area.District)] <- "no.record"
levels(SEOtix$Groundfish.Mgt.Area.District)
head(SEOtix)
  #filter(Groundfish.Mgt.Area.District == "EYKT" |
   #                        Groundfish.Mgt.Area.District == "NSEO"  |
   #                        Groundfish.Mgt.Area.District == "CSEO" |
  #                         Groundfish.Mgt.Area.District == "SSEO")
harv.dist<-SEOtix %>% filter(Year > 1991 & Species.Code == 200)
ex<-harv.dist[harv.dist$Year == 2000 & 
                harv.dist$IPHC.Statistical.Area == 151,]
sum(ex$CFEC.Whole.Pounds)
sum(ex$CFEC.Whole.Pounds[ex$Groundfish.Mgt.Area.District == "no.record"])


View(ex)
head(harv.dist)
unique(harv.dist$Year)
post90<-seq(1992,2020,1)
#gfu<-unique(harv.dist$Groundfish.Mgt.Area.District) #c("EYKT","NSEO","CSEO","SSEO","NSEI","SSEI","")
gfu<-unique(harv.dist$SEdist_schwag)
ISA<-unique(harv.dist$IPHC.Statistical.Area)
i<-1
prop.df<-data.frame()

for (y in post90){  #y<-post90[1]
  for (s in ISA){      #s<-ISA[1] in ISA g<-gfu[1]   #g in gfu
    Dat<-harv.dist[harv.dist$Year == y & 
                     harv.dist$IPHC.Statistical.Area == s,]
    
    if (nrow(Dat) < 1) {} else 
    #View(Dat)
    {harv<-sum(Dat$CFEC.Whole.Pounds, na.rm=T)
    gfus<-unique(Dat$Groundfish.Mgt.Area.District)
    gfus2<-unique(Dat$ADFG.Management.Area.Code)
    for (g in gfus){   #head(Dat)
      Dat2<-Dat[Dat$Groundfish.Mgt.Area.District == g,]
      prop.df[i,"Year"]<-y
      prop.df[i,"gfu"]<-g
      prop.df[i,"ISA"]<-s
      prop.df[i,"gfu.catch"]<-sum(Dat2$CFEC.Whole.Pounds)
      prop.df[i,"ISA.catch"]<-sum(Dat$CFEC.Whole.Pounds)
      prop.df[i,"prop.gfu.in.ISA"]<-prop.df[i,"gfu.catch"]/prop.df[i,"ISA.catch"]
      #comp[i,"prop"]<-
      i<-i+1
    }}
  }
}

head(prop.df)
eg<-prop.df[prop.df$gfu == "CSEO",]

seo.prop<-prop.df %>% filter(gfu=="no.record"|
                               gfu=="EYKT"|
                               gfu=="NSEO"|
                               gfu=="CSEO"|
                               gfu=="SSEO")

ggplot(seo.prop) +
  geom_point (aes(Year,prop.gfu.in.ISA,col=gfu))+
  facet_wrap(~ISA)

SEOtix<-SEOtix %>% filter(Gear.Description == "Longline") %>%
  mutate(h.code = substr(Permit.Fishery,1,1))

pre_90<-SEOtix %>% filter(Year <1991)
unique(pre_90$Species.Code)

nrow(SEOtix)
unique(SEOtix$Year)
unique(SEOtix$Stat.Area)
unique(SEOtix$Statistical.Area)
unique(SEOtix$Gear.Description)
unique(SEOtix$Permit.Fishery)
unique(SEOtix$h.code)

SEOtix %>% #filter(h.code == "B") %>% 
  group_by(Year,Groundfish.Mgt.Area.District,Species.Code) %>%
  summarise(Whole_lbs = sum(CFEC.Whole.Pounds, na.rm=T)) ->Hal.longline
unique(Hal.longline$Year)
unique(Hal.longline$Species.Code)

pre90<-SEOtix %>% filter(Year <1991)
unique(pre90$Species.Code)

Hal.catch<-Hal.longline %>% filter(Species.Code == 200)
Hal.catch$hal_mt<-Hal.catch$Whole_lbs*0.00045359
unique(Hal.catch$Year)

YE.catch<-Hal.longline %>% filter(Species.Code == 145)
YE.catch$ye_mt<-YE.catch$Whole_lbs*0.00045359

Hal.catch.SEO<-Hal.catch %>% group_by(Year) %>%
  summarize(tot.mt = sum(hal_mt))

unique(Hal.catch.SEO$Year)


#************************************************************************
  HA.Harv<-RheaHal #read.csv("Data/Harvests/halibut_catch_data.csv", header=T)
  LL<-HA.Harv[HA.Harv$gear.description == "Longline",]
  head(LL)
  colnames(LL)
  
  LL<-rename(LL, SEdist = groundfish.mgt.area.district)
  LL<-rename(LL, Year = year.landed)
  Hali<-LL
  
  plot(LL$Year,LL$round_lbs_hali)
  Hcheck<-data.frame()
  Yrs<-unique(LL$Year)
  i<-1
  for (y in Yrs){  #y<-Yrs[1]
    Dat<-LL[LL$Year == y,]; Dat<-Dat[Dat$SEdist=="SSEO"|Dat$SEdist=="CSEO"|Dat$SEdist=="NSEO"|Dat$SEdist=="EYKT",]
    Hcheck[i,"Year"]<-y; Hcheck[i,"hal"]<-sum(Dat$round_lbs_hali); i<-i+1
  }
  plot(Hcheck$Year, Hcheck$hal)
 
Hali<-LL
Hali<-Hali[Hali$SEdist != "NSEI" & Hali$SEdist != "SSEI",]
str(Hali)

Hali$mt_hali<-Hali$round_lbs_hali*0.00045359

Hcheck$hal_mt<-Hcheck$hal*0.00045359
#***********************************************************************
HA.Harv2<-RheaHal2 #read.csv("Data/Harvests/halibut_catch_data.csv", header=T)
Hali2<-HA.Harv2[HA.Harv2$gear.description == "Longline",]
head(Hali2)
Hali2$Year<-Hali2$year.landed

plot(Hali2$Year,Hali2$round.lbs)
Hcheck2<-data.frame()
Yrs<-unique(Hali2$Year)
i<-1
for (y in Yrs){  #y<-Yrs[1]
  Dat<-Hali2[Hali2$Year == y,]; Dat<-Dat[Dat$mgmt.area=="SSEO"|Dat$mgmt.area=="CSEO"|Dat$mgmt.area=="NSEO"|Dat$mgmt.area=="EYKT",]
  Hcheck2[i,"Year"]<-y; Hcheck2[i,"hal"]<-sum(Dat$round.lbs); i<-i+1
}
plot(Hcheck2$Year, Hcheck$hal)

Hali2<-Hali2[Hali2$mgmt.area != "NSEI" & Hali2$mgmt.area != "SSEI",]
str(Hali2)

Hali2$mt_hali<-Hali2$round.lbs*0.00045359

Hcheck2$hal_mt<-Hcheck2$hal*0.00045359
#**************************************************************************
str(RPHal)
RPHal$SEOmt<-RPHal$SE.O*0.90718474

plot(RPHal$SEOmt~RPHal$Year)
lines(Hcheck$hal_mt~Hcheck$Year, col="blue")

CurSEOHal<-RPHal[RPHal$Year > 1974,]

plot(CurSEOHal$SEOmt~CurSEOHal$Year, ylim=c(0,max(Hcheck$hal_mt)))
lines(Hcheck$hal_mt~Hcheck$Year, col="blue")
lines(Hcheck2$hal_mt~Hcheck2$Year, col="cyan")


sub<-unique(Hali$SEdist)
for ()
  
  lines(Hali$mt_hali[Hali$SEdist=="EYKT"]~Hali$Year[Hali$SEdist=="EYKT"],col="purple")


