################################################################################
##IPHC Survey Data exploration
##
library(dplyr)
library(boot)
library(ggplot2)

#wd="C:/Users/pjjoy/Documents/Groundfish Biometrics/Yelloweye_Production_Models/"
#setwd(wd)
#getwd()

IPHCfunction<-function(){

BUT2C<-read.csv("Data/IPHC Set and Pacific halibut data 2C.csv", header=T)
YE2C<-read.csv("Data/Non-Pacific halibut data YE 2C.csv", header=T)

HA.Harv<-read.csv("Data/halibut_catch_data.csv", header=T)

Sur2C<- BUT2C %>% full_join(YE2C, by="Stlkey")

Sur2C[Sur2C$Stlkey == 20200208,]

#unique(Sur2C$IPHC.Reg.Area)
#unique(Sur2C$IPHC.Charter.Region)
#unique(Sur2C$IPHC.Stat.Area)

Sur2C<- Sur2C %>%
  mutate(SEdist = ifelse(IPHC.Stat.Area %in% c(182:184,171,173,174,161:163),"NSEI", 
                         ifelse(IPHC.Stat.Area %in% c(142:144,152,153),"SSEI",
                                ifelse(IPHC.Stat.Area %in% c(140,141,151,150,160,170,181),
                                       ifelse(-MidLon.fished>137,"EYKT",
                                              ifelse(MidLat.fished>=57.5,"NSEO",
                                                     ifelse(MidLat.fished<57.5 & MidLat.fished>=56,"CSEO",
                                                            ifelse(MidLat.fished<56,"SSEO",NA)))),NA))))

## Lets bring in 3A and pull out surveys in the Yakutat area...
BUT3A<-read.csv("Data/IPHC Set and Pacific halibut data 3A.csv", header=T)
YE3A<-read.csv("Data/Non-Pacific halibut data YE 3A.csv", header=T)

Sur3A<- BUT3A %>% full_join(YE3A, by="Stlkey")

Sur3A<-Sur3A %>%
  mutate(SEdist = ifelse(-MidLon.fished <= 137 & MidLat.fished >= 57.5,"NSEO",
                         ifelse(-MidLon.fished > 137 & -MidLon.fished <139 & MidLat.fished > 54.5,  #Lon should be less than 138 but included here to increase sample size of the index... 
                                "EYKT","whocares")))
nrow(Sur3A)
Sur3A<-Sur3A[Sur3A$SEdist != "whocares",]

## connect the two data sets
Survey<-rbind(Sur3A,Sur2C)
#make one more column to define inside versus outside for dealing with older data sets...
Survey<-Survey %>% 
  mutate(In.Out = ifelse(SEdist %in% c("NSEO","CSEO","SSEO","EYKT"),"SEO",
                         ifelse(SEdist %in% c("NSEI","SSEI"),"SEI",NA)))

Survey<-Survey[order(Survey$Year.x),]
Survey<-Survey[,c(1:32,49:57)]
Years<-sort(unique(Survey$Year.x))
GFMA<-unique(Survey$SEdist)
#==============================================================================================
## Load the port samples that were downloaded from oceansAK; part of function

Port1<-read.csv("Data/port_sampling_bio_data.csv")
Port2<-read.csv("Data/port_sampling_bio_data2.csv")
Port<-rbind(Port1,Port2)
Port$Year<-as.integer(as.character(Port$ï..Year))
str(Port)
unique(Port$Groundfish.Management.Area.Code)  #EYAK = EYKT
Port<-Port[Port$Groundfish.Management.Area.Code != "",]
Port<-subset(Port, !is.na(Weight.Kilograms))

uYEkg<- Port %>%
  group_by(Year,Groundfish.Management.Area.Code) %>%
  summarize(uYEkg = mean(Weight.Kilograms, na.rm=TRUE),
            vYEkg = var(Weight.Kilograms, na.rm=TRUE),
            NyeKG = length(Weight.Kilograms))
uYEkg<-as.data.frame(uYEkg); uYEkg

Survey$mean.YE.kg<-NA
Survey$var.YE.kg<-NA
Survey$N.YE.kg<-NA
Survey$HooksRetrieved<-as.numeric(Survey$HooksRetrieved)
Survey$YE.exp<-Survey$HooksRetrieved/Survey$HooksObserved
Survey<-Survey %>% mutate(Number.Observed = ifelse(is.na(Number.Observed),0,Number.Observed))
Survey<-Survey %>% mutate(HooksObserved  = ifelse(is.na(HooksObserved ),140,HooksObserved ))

## Get best data for average weights and add columns to survey data
for (i in Years){    #i<-Years[1]
  for (j in GFMA){   #j<-GFMA[1]
    P<-Port[Port$Groundfish.Management.Area.Code == j,]
    #reach back from year i and see if you can get 150 weights...
    Nw<-nrow(P[P$Year == i,])
    k<-0
    while (Nw < 75 & i-k >= min(Port$Year) ){         #go back as far as needed to get 150 weights... 
      k<-k+1
      Nw<-nrow(P[P$Year >=i-k & P$Year <= i,])   #P[P$Year == 1990:2002,]
    }
    m<-0
    while (Nw < 75 & i+m <= max(Port$Year)) {                        #go forward if you failed to find enough weights... 
      m<-m+1                                                        
      Nw<-nrow(P[P$Year >=i-k & P$Year <= i+m,])   #P[P$Year == 1986,]
    }
    #get average weights of YE...
    Sample<-P$Weight.Kilograms[P$Year >=i-k & P$Year <= i+m]
    Survey[,"mean.YE.kg"][Survey$Year.x == i & Survey$SEdist == j]<-mean(Sample)
    Survey[,"var.YE.kg"][Survey$Year.x == i & Survey$SEdist == j]<-var(Sample)
    Survey[,"N.YE.kg"][Survey$Year.x == i & Survey$SEdist == j]<-length(Sample)
  }
}

return(Survey)
}

#===========================================================================

Survey<-IPHCfunction()
colnames(Survey)

head(Survey,10)

#=========================================================================================
## HALIBUT HARVESTS

Halibut<-function(){
  HA.Harv<-read.csv("Data/Harvests/halibut_catch_data.csv", header=T)
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
return(Hali)
}
Hali<-Halibut()

#=============================================================================================
##YE Harvests

YE.landed.bycatch<-function() {
YE.Harv<-read.csv("Data/Harvests/yelloweye_catch_data.csv", header=T)
YE.Harv<- YE.Harv %>%
  mutate(IO = ifelse(mgmt_area %in% c("SSEO","CSEO","NSEO","EYKT","SEO"),"SEO", 
                     ifelse(mgmt_area %in% c("NSEI","SSEI"),"SEI",NA)))
YE.Harv$h.code<-substr(YE.Harv$fishery,1,1)
str(YE.Harv)
MA<-unique(YE.Harv$mgmt_area)
Years<-unique(YE.Harv$year); Years<-Years[order(Years)]
#Get YE landed bycatch by subdistrict for bootstrap below... 
YE.lnd.by<-data.frame()
a<-1

for (i in MA){     #i<-MA[1]
  for (y in Years){  #y<-Years[1]
    Com<-YE.Harv[YE.Harv$year == y & YE.Harv$mgmt_area == i,]
    YE.lnd.by[a,"Year"]<-y
    YE.lnd.by[a,"mgmt_area"]<-i
    YE.lnd.by[a,"YE_fr_hal_fish.kgs"]<-sum(Com$whole_pounds[Com$h.code == "B" & Com$gear == "Longline"])*0.45359237
    a<-a+1
  }
}
return(YE.lnd.by)
}

YE.lnd.by<-YE.landed.bycatch()
#============================================================================================
#bootstrap estimates of yelloweye catch relative to halibut catch by subdistrict and year
#bycatch from WCPUE*Halibut landings
#discards = IPHC bycatch est - landed bycatch
Obs.Rates<-read.csv("Data/SEO_YE_NOAA_Observers.csv")

Survey.boot<-data.frame()
YE.IPHC.CPUE<-data.frame()
Subs<-unique(Survey$SEdist)
Years<-unique(Survey$Year.x)
nboot<-1000

j<-1
for (y in Years) {  #y<-Years[1]
  for (s in Subs){  #s<-Subs[4]
    Dat<-Survey[Survey$Year.x == y & Survey$SEdist == s,]
    Harv<-Hali[Hali$Year == y & Hali$SEdist == s,]
    Land.YEby<-YE.lnd.by[YE.lnd.by$Year == y & YE.lnd.by$mgmt_area == s, ]
    Stations<-unique(Dat$Station.x)
    
    if (y < 2013){
      min.Discard<-mean(Obs.Rates$YE.Discard.rate.boot)
      min.D.var<-var(Obs.Rates$YE.Disrate.boot.var)
    } else {
      min.Discard<-Obs.Rates$YE.Discard.rate.boot[Obs.Rates$Year == y]
      min.D.var<-Obs.Rates$YE.Disrate.boot.var[Obs.Rates$Year == y]
    }
    
    WCPUEi<-vector()
    CPUEi<-vector()
    i<-1
    for (st in Stations){    #i<-Stations[9]   length(Stations)
      Stat.Dat<-Dat[Dat$Station.x == st,]; Stat.Dat
      #debug
      #if (nrow(Stat.Dat) > 1){aaa} else {}
      CPUE<-mean(Stat.Dat$Number.Observed/Stat.Dat$HooksObserved)
      if (CPUE == 0){
        C<-0
      } else {
        C<-CPUE*Stat.Dat$HooksRetrieved
      }
      CPUEi[i]<-CPUE
      K<-C*Dat$mean.YE.kg[1]
      WCPUE <-K/mean(Stat.Dat$O32.Pacific.halibut.weight)
      WCPUEi[i]<-WCPUE
      i<-i+1
    }
    #bootstrap WCPUE
    Out<-data.frame()
    CPUE.out<-data.frame()
    for (ii in 1:nboot){ #nboot<-1000  i<-1
      Resamp<-sample(WCPUEi,length(WCPUEi),replace=T)
      Out[ii,"WCPUE"]<-mean(Resamp, na.rm=T)
      Out[ii,"var.WCPUE"]<-var(Resamp, na.rm=T)
      Resamp2<-sample(CPUEi,length(CPUEi),replace=T)
      CPUE.out[ii,"CPUE"]<-mean(Resamp2, na.rm=T)
      CPUE.out[ii,"CPUE.var"]<-var(Resamp2, na.rm=T)
      #bootstrap bycatch and discards 12-20-21: exact same results if you pu
      Out[ii,"YE.bycatch.kg"]<-sum(Harv$round_lbs_hali)*0.45359237*Out[ii,"WCPUE"]
      Out[ii,"raw.YE.discards.kg"]<-Out[ii,"YE.bycatch.kg"]-Land.YEby$YE_fr_hal_fish.kgs #max(1,Out[ii,"YE.bycatch.kg"]-Land.YEby$YE_fr_hal_fish.kgs, na.rm=T)
      min.D.boot<-min(1,rnorm(1,mean=min.Discard ,sd=sqrt(min.D.var))*Land.YEby$YE_fr_hal_fish.kgs)
      Out[ii,"trunc.YE.discards.kg"]<-max(min.D.boot,Out[ii,"YE.bycatch.kg"]-Land.YEby$YE_fr_hal_fish.kgs, na.rm=T)
    }
    #par(mfrow=c(2,1))
    #hist((Out$YE.bycatch.kg), breaks = 25)
    #hist((Out$raw.YE.discards.kg), breaks=25)
    
    Survey.boot[j,"Year"]<-y
    Survey.boot[j,"SEdist"]<-s
    Survey.boot[j,"YE_to_HA_ratio"]<-unname(quantile(Out$WCPUE,c(0.5), na.rm=T))  #mean WCPUE from Tribuzio
    Survey.boot[j,"ratio.var"]<-var(Out$WCPUE)
    Survey.boot[j,"ratio.cv"]<-sqrt(var(Out$WCPUE))/Survey.boot[j,"YE_to_HA_ratio"]
    
    YE.IPHC.CPUE[j,"Year"]<-y
    YE.IPHC.CPUE[j,"SEdist"]<-s
    YE.IPHC.CPUE[j,"YE_CPUE"]<-unname(quantile(CPUE.out$CPUE,c(0.5)))
    YE.IPHC.CPUE[j,"YE_CPUE.var"]<-var(CPUE.out$CPUE)
    
    ##apply to actual harvest data to get bycatch
    
    Survey.boot[j,"Halibut.landed.kg"]<-sum(Harv$round_lbs_hali)*0.45359237
    Survey.boot[j,"YE.bycatch.kg"]<-Survey.boot[j,"Halibut.landed.kg"]*Survey.boot[j,"YE_to_HA_ratio"]
    Survey.boot[j,"YE.bycatch.var"]<-Survey.boot[j,"ratio.var"]*Survey.boot[j,"Halibut.landed.kg"]^2
    Survey.boot[j,"YE.bycatch.cv"]<-sqrt(Survey.boot[j,"YE.bycatch.var"])/Survey.boot[j,"YE.bycatch.kg"]
    
    #subtract known YE bycatch landings from bycatch estimates to get discard estimates
    
    Survey.boot[j,"YE.landed.bycatch"]<-Land.YEby$YE_fr_hal_fish.kgs
    #Survey.boot[j,"YE.discards.kg"]<-Survey.boot[j,"YE.bycatch.kg"]-Land.YEby$YE_fr_hal_fish.kgs
    Survey.boot[j,"YE.discards.kg"]<-max(min.Discard,Survey.boot[j,"YE.bycatch.kg"]-Land.YEby$YE_fr_hal_fish.kgs, na.rm=T)
    Survey.boot[j,"cv.YE.discards.kg"]<-sqrt(Survey.boot[j,"YE.bycatch.var"])/Survey.boot[j,"YE.discards.kg"]
    
    Survey.boot[j,"YEbycatch.boot"]<-median(Out$YE.bycatch.kg)
    Survey.boot[j,"var.YEby.boot"]<-var(Out$YE.bycatch.kg)
    Survey.boot[j,"cv.YEby.boot"]<-sqrt(var(Out$YE.bycatch.kg))/median(Out$YE.bycatch.kg)
    
    Survey.boot[j,"raw.YEdiscard.boot"]<-median(Out$raw.YE.discards.kg)
    Survey.boot[j,"var.rawYEdisc.boot"]<-var(Out$raw.YE.discards.kg)
    Survey.boot[j,"cv.rawYEdisc.boot"]<-sqrt(var(Out$raw.YE.discards.kg))/median(Out$raw.YE.discards.kg)
    
    Survey.boot[j,"trunc.YEdiscard.boot"]<-median(Out$trunc.YE.discards.kg)
    Survey.boot[j,"var.truncYEdisc.boot"]<-var(Out$trunc.YE.discards.kg)
    Survey.boot[j,"cv.truncYEdisc.boot"]<-sqrt(var(Out$trunc.YE.discards.kg))/median(Out$trunc.YE.discards.kg)
    
    j<-j+1
  }
}

str(Survey.boot)
par(mfrow=c(3,2))
for (s in Subs) {
  dat<-Survey.boot[Survey.boot$SEdist == s,]
  plot(dat$YE.discards.kg ~ dat$Year, type="l", main=s)
  lines(dat$raw.YEdiscard.boot~ dat$Year, type="l", col="blue" )
  lines(dat$trunc.YEdiscard.boot~ dat$Year, type="l", col="darkgreen")
}

#look at wcpue and observer data
str(Obs.Rates)
Subs<-Subs[c(1,2,4,5)]
par(mfrow=c(2,2), mar=c(4,4,2,1))
for (s in Subs){  #s<-Subs[1]
  dat<-Survey.boot[Survey.boot$SEdist == s,]
  plot(dat$YE_to_HA_ratio ~ dat$Year, type="p", main=s, ylim=c(0,0.3), 
       pch=17, col="goldenrod4", ylab="WCPUE/Bycatch rate", xlab="Year")
  y0<-pmax(dat$YE_to_HA_ratio-1.68*dat$ratio.cv*dat$YE_to_HA_ratio,0)
  y1<-dat$YE_to_HA_ratio+1.68*dat$ratio.cv*dat$YE_to_HA_ratio
  colp<-adjustcolor("goldenrod", alpha.f=0.5)
  polygon(c(dat$Year,rev(dat$Year)),c(y0,rev(y1)),border=NA, col=colp)
  points(Obs.Rates$YE.Bycatch.rate ~ Obs.Rates$Year, col="forestgreen", pch=18)
  y0<-pmax(Obs.Rates$YE.Bycatch.rate-1.68*Obs.Rates$YE.Byrate.boot.cv*Obs.Rates$YE.Bycatch.rate,0)
  y1<-Obs.Rates$YE.Bycatch.rate+1.68*Obs.Rates$YE.Byrate.boot.cv*Obs.Rates$YE.Bycatch.rate
  arrows(x0=Obs.Rates$Year, y0=y0, x1=Obs.Rates$Year, y1=y1,
          code=3, col="forestgreen", lwd=1, angle=90, length=0.01)
  if (s == Subs[1]){
    legend(x="topleft",c("WCPUE","SEO Observed bycatch rate"),pch=c(17,18,18),
           col=c("goldenrod4","forestgreen"),bty="n",text.col=c("goldenrod4","forestgreen"),
           cex=0.8)
  }
}

Survey.boot$YE.bycatch.cv
Survey.boot$cv.YE.discards.kg
Survey.boot$cv.truncYEdisc.boot

hist(Survey.boot$cv.YE.discards.kg, breaks=100)
mean(Survey.boot$cv.YE.discards.kg)
min(Survey.boot$cv.YE.discards.kg)
max(Survey.boot$cv.YE.discards.kg)
median(Survey.boot$cv.YE.discards.kg)
mean(Survey.boot$cv.YE.discards.kg)

## now need to go back in time to get estimated bycatch
#1st we need to scale presurvey cvs based on scale of hindcast discards...
# cv's are inversely related to the size of the discard estimate (but none are small!)
{
not0<-Survey.boot[Survey.boot$YE.discards.kg>1,]
plot(log(not0$cv.YE.discards.kg) ~ log(not0$YE.discards.kg))
log.cvs<-log(not0$cv.YE.discards.kg)
log.disc<-log(not0$YE.discards.kg)
#cv.dic<-lm(log(not0$cv.YE.discards.kg) ~ log(not0$YE.discards.kg))
cv.dic<-lm(log.cvs ~ log.disc)
abline(cv.dic); summary(cv.dic)
pred.cv<-predict(cv.dic)
points(pred.cv ~ log(not0$YE.discards.kg),col="red")

x<-log(140000)
#test.disc<-data.frame(x=c(x))
predt<-predict(cv.dic, newdata=data.frame(log.disc=c(x)), interval = "prediction")
summary(predt)
syx<-summary(cv.dic)$sigma
cf<-exp((syx^2)/2)
cf<-exp((summary(cv.dic)$sigma^2)/2)
exp(predt)*cf

pp<-log(exp(predt)*cf); pp2<-predt
points(pp[1] ~ x,col="green"); points(pp2[1] ~ x,col="darkgreen")

plot(not0$cv.YE.discards.kg ~ not0$YE.discards.kg)
points(exp(predt[1])*cf ~ exp(x), col="red")
}
#prediction function for bootstrap
cv.pred<-function(known.cv,known.disc,hind.disc){
  log.cvs<-log(known.cv)
  log.disc<-log(known.disc)
  cv.dic<-lm(log.cvs ~ log.disc)
  x<-log(hind.disc)
  predt<-predict(cv.dic, newdata=data.frame(log.disc=c(x)), interval = "prediction")
  cf<-exp((summary(cv.dic)$sigma^2)/2)
  pred.cv<-exp(predt[1])*cf
  return(pred.cv)
}
cv.pred(not0$cv.YE.discards.kg,not0$YE.discards.kg,140000)
#============================================================================

hyears<-unique(Hali$Year)
ehyears<-hyears[1:23]
Presurvey<-data.frame()
i<-1

for (y in ehyears){        #y<-ehyears[1]
  for (s in Subs){         #s<-Subs[1]
    Ratio<-Survey.boot[Survey.boot$SEdist == s,]
    not0<-Ratio[Ratio$YE.discards.kg>1,]
    Harv<-Hali[Hali$Year == y & Hali$SEdist == s,]
    Land.YEby<-YE.lnd.by[YE.lnd.by$Year == y & YE.lnd.by$mgmt_area == s, ]
    Presurvey[i,"Year"]<-y
    Presurvey[i,"SEdist"]<-s
    Presurvey[i,"YE_to_HA_ratio"]<-mean(Ratio$YE_to_HA_ratio)
    Presurvey[i,"ratio.var"]<-max(Ratio$ratio.var)
    Presurvey[i,"ratio.cv"]<-max(Ratio$ratio.cv)
    Presurvey[i,"Halibut.landed.kg"]<-sum(Harv$round_lbs_hali)*0.45359237
    Presurvey[i,"YE.bycatch.kg"]<-Presurvey[i,"Halibut.landed.kg"]*Presurvey[i,"YE_to_HA_ratio"]
    Presurvey[i,"YE.bycatch.var"]<-Presurvey[i,"ratio.var"]*Presurvey[i,"Halibut.landed.kg"]^2
 #   Presurvey[i,"YE.bycatch.cv"]<-max(Ratio$YE.bycatch.cv, na.rm=T)#sqrt(Presurvey[i,"YE.bycatch.var"])/Presurvey[i,"YE.bycatch.kg"]
    Presurvey[i,"YE.bycatch.cv"]<-sqrt(Presurvey[i,"YE.bycatch.var"])/Presurvey[i,"YE.bycatch.kg"]
    
    if (nrow(Land.YEby) < 1) {
      Presurvey[i,"YE.landed.bycatch"]<-0
      Presurvey[i,"YE.discards.kg"]<-Presurvey[i,"YE.bycatch.kg"]
#      Presurvey[i,"cv.YE.discards.kg"]<-median(Ratio$cv.YE.discards.kg, na.rm=T)#sqrt(Presurvey[i,"YE.bycatch.var"])/Presurvey[i,"YE.discards.kg"]
      Presurvey[i,"cv.YE.discards.kg"]<-cv.pred(not0$cv.YE.discards.kg,not0$YE.discards.kg,Presurvey[i,"YE.discards.kg"])
      
    } else {
      Presurvey[i,"YE.landed.bycatch"]<-Land.YEby$YE_fr_hal_fish.kgs
      Presurvey[i,"YE.discards.kg"]<-max(Presurvey[i,"YE.bycatch.kg"]-Land.YEby$YE_fr_hal_fish.kgs,10)
 #     Presurvey[i,"cv.YE.discards.kg"]<-median(Ratio$cv.YE.discards.kg, na.rm=T)#sqrt(Presurvey[i,"YE.bycatch.var"])/Presurvey[i,"YE.discards.kg"]
      Presurvey[i,"cv.YE.discards.kg"]<-cv.pred(not0$cv.YE.discards.kg,not0$YE.discards.kg,Presurvey[i,"YE.discards.kg"])
      
   }
    Presurvey[i,"YEbycatch.boot"]<-NA#median(Out$YE.bycatch.kg)
    Presurvey[i,"var.YEby.boot"]<-NA#var(Out$YE.bycatch.kg)
    Presurvey[i,"cv.YEby.boot"]<-NA#sqrt(var(Out$YE.bycatch.kg))/median(Out$YE.bycatch.kg)
    
    Presurvey[j,"raw.YEdiscard.boot"]<-NA#median(Out$raw.YE.discards.kg)
    Presurvey[j,"var.rawYEdisc.boot"]<-NA#var(Out$raw.YE.discards.kg)
    Presurvey[j,"cv.rawYEdisc.boot"]<-NA#sqrt(var(Out$raw.YE.discards.kg))/median(Out$raw.YE.discards.kg)
    
    Presurvey[j,"trunc.YEdiscard.boot"]<-NA#median(Out$trunc.YE.discards.kg)
    Presurvey[j,"var.truncYEdisc.boot"]<-NA#var(Out$trunc.YE.discards.kg)
    Presurvey[j,"cv.truncYEdisc.boot"]<-NA#sqrt(var(Out$trunc.YE.discards.kg))/median(Out$trunc.YE.discards.kg)
    i<-i+1
  }
}

unique(Survey.boot$Year)
unique(Presurvey$Year)
head(Survey.boot,20)   
head(Presurvey,20)

Presurvey<-Presurvey[!is.na(Presurvey$Year),]
Bycatch<-rbind(Presurvey,Survey.boot)


##take a look at estimated bycatch, landed bycatch and estimated discards  
Bycatch$cv.YE.discards.kg
{
dists<-unique(Bycatch$SEdist)
dists<-dists[c(1,2,4,5)]
par(mfrow=c(2,2), mar=c(4,4,2,2))    #par(mfrow=c(1,1))
for (d in dists){ #d<-dists[1]
  Dat<-Bycatch[Bycatch$SEdist == d,]
  plot(Dat$YE.bycatch.kg~Dat$Year, type="l", main=d, col="darkgrey",
       ylim=c(0,1.5*max(Dat$YE.bycatch.kg, na.rm=T)), cex=0.8,
       ylab="kg", xlab="year")
 # abline(v=1998,col="red", lty=2)
#       ylim=c(min(Dat$YE.discards.kg, na.rm=T),max(Dat$YE.bycatch.kg, na.rm=T)))
  y0<-pmax(Dat$YE.bycatch.kg-1.68*Dat$YE.bycatch.cv*Dat$YE.bycatch.kg,0)
  y1<-Dat$YE.bycatch.kg+1.68*Dat$YE.bycatch.cv*Dat$YE.bycatch.kg
#  arrows(x0=Dat$Year, y0=y0, x1=Dat$Year, y1=y1,
#         code=3, col="darkgrey", lwd=1, angle=90, length=0.01)
  colp<-adjustcolor("lightblue", alpha.f=0.999)
  polygon(c(Dat$Year,rev(Dat$Year)),c(y0,rev(y1)),border=NA, col=colp)
  points(Dat$YE.bycatch.kg~Dat$Year, type="l", main=d, col="blue", lwd=1)
  
#  points(Dat$YEbycatch.boot~Dat$Year, type="l", col="blue")
#  y0<-Dat$YEbycatch.boot-1.68*Dat$cv.YEby.boot*Dat$YEbycatch.boot
#  y1<-Dat$YEbycatch.boot+1.68*Dat$cv.YEby.boot*Dat$YEbycatch.boot
#  arrows(x0=Dat$Year, y0=y0, x1=Dat$Year, y1=y1,
#         code=3, col="blue", lwd=1, angle=90, length=0.01)
#  polygon(c(Dat$Year,rev(Dat$Year)),c(y0,rev(y1)),border=NA,col=adjustcolor("lightblue", alpha.f=0.999))
#  points(Dat$YEbycatch.boot~Dat$Year, type="l", col="blue")

  
  
#  points(Dat$raw.YEdiscard.boot~Dat$Year, type="l", col="darkcyan")
#  y0<-Dat$raw.YEdiscard.boot-1.68*Dat$cv.rawYEdisc.boot*Dat$raw.YEdiscard.boot
#  y1<-Dat$raw.YEdiscard.boot+1.68*Dat$cv.rawYEdisc.boot*Dat$raw.YEdiscard.boot
#  arrows(x0=Dat$Year, y0=y0, x1=Dat$Year, y1=y1,
#         code=3, col="darkcyan", lwd=1, angle=90, length=0.01)
  
#  points(Dat$trunc.YEdiscard.boot~Dat$Year, type="l", col="purple")
  points(Dat$YE.discards.kg~Dat$Year, type="l", col="purple", lwd=1)
#  y0<-Dat$trunc.YEdiscard.boot-1.68*Dat$cv.truncYEdisc.boot*Dat$trunc.YEdiscard.boot
#  y1<-Dat$trunc.YEdiscard.boot+1.68*Dat$cv.truncYEdisc.boot*Dat$trunc.YEdiscard.boot
  y0<-Dat$YE.discards.kg-1.68*Dat$cv.YE.discards.kg*Dat$YE.discards.kg
  y0p<-pmax(y0,0)
  y1<-Dat$YE.discards.kg+1.68*Dat$cv.YE.discards.kg*Dat$YE.discards.kg
  colp<-adjustcolor("purple", alpha.f=0.25)
  polygon(c(Dat$Year,rev(Dat$Year)),c(y0p,rev(y1)),border=NA, col=colp)
  #arrows(x0=Dat$Year, y0=y0, x1=Dat$Year, y1=y1,
  #       code=3, col="purple", lwd=1, angle=90, length=0.01)
  points(Dat$YE.landed.bycatch~Dat$Year, type="l", col="darkgreen", lwd=1, lty=1)
  points(Dat$YE.bycatch.kg~Dat$Year, type="l", main=d, col=adjustcolor("blue", alpha.f=0.5), lwd=1)
  
  if (d == dists[1]){
    legend(x="topleft",c("Expected bycatch","Landed bycatch","Est. discards"),pch=c(18,18,18),
           col=c("blue","darkgreen","purple"),bty="n",text.col=c("blue","darkgreen","purple"),
           cex=0.8)
  }
}
}
#========================================================================================
# Save subdistrict level data 

write.csv(Bycatch,"Data/SEsubd_YE_IPHCbycatch.csv")

#=========================================================================================
# Combine subdistricts for entire SEO

SEO.Bycatch<-Bycatch[Bycatch$SEdist == "SSEO" | Bycatch$SEdist == "CSEO" |
                       Bycatch$SEdist == "NSEO" | Bycatch$SEdist == "EYKT",]
head(SEO.Bycatch)
SEO.YEbycatch<-data.frame()
Ratio<-SEO.Bycatch #[Survey.boot$SEdist == s,]
not0<-Ratio[Ratio$YE.discards.kg>1,]
i<-1
for (y in hyears){  #y<-hyears[35]
  Dat<-SEO.Bycatch[SEO.Bycatch$Year == y,]
  SEO.YEbycatch[i,"Year"] <- y
  SEO.YEbycatch[i,"mean.ratio"]<-mean(Dat$YE_to_HA_ratio)
  SEO.YEbycatch[i,"ratio.var"]<-sum(Dat$ratio.var)/(4^2)
  SEO.YEbycatch[i,"ratio.cv"]<-sqrt(SEO.YEbycatch[i,"ratio.var"])/SEO.YEbycatch[i,"mean.ratio"] 
  
  SEO.YEbycatch[i,"YE.landed.bycatch"]<-sum(Dat$YE.landed.bycatch)
  
  SEO.YEbycatch[i,"YE_bycatch.kg"] <- sum(Dat$YE.bycatch.kg)
  SEO.YEbycatch[i,"YE_bycatch_var"]<-sum(Dat$YE.bycatch.var)
  SEO.YEbycatch[i,"bycatch_cv"]<-sqrt(SEO.YEbycatch[i,"YE_bycatch_var"])/SEO.YEbycatch[i,"YE_bycatch.kg"] #max(Dat$ratio.cv)
  
  SEO.YEbycatch[i,"YE_discards.kg"]<-sum(Dat$YE.discards.kg)
  SEO.YEbycatch[i,"YE_discards_cv"]<-cv.pred(not0$cv.YE.discards.kg,not0$YE.discards.kg,SEO.YEbycatch[i,"YE_discards.kg"])
  
  SEO.YEbycatch[i,"rawYE.discards.kg"]<-sum(Dat$raw.YEdiscard.boot)
  SEO.YEbycatch[i,"var.rawDisc"]<-sum(Dat$var.rawYEdisc.boot)
  SEO.YEbycatch[i,"cv.rawDisc"]<-sqrt(SEO.YEbycatch[i,"var.rawDisc"])/SEO.YEbycatch[i,"rawYE.discards.kg"]
  
  SEO.YEbycatch[i,"truncYE.discards.kg"]<-sum(Dat$trunc.YEdiscard.boot)
  SEO.YEbycatch[i,"var.truncDisc"]<-sum(Dat$var.truncYEdisc.boot)
  SEO.YEbycatch[i,"cv.truncDisc"]<-sqrt(SEO.YEbycatch[i,"var.truncDisc"])/SEO.YEbycatch[i,"truncYE.discards.kg"]
  
  i<-i+1
}

{
par(mfrow=c(1,1))
Dat<-SEO.YEbycatch
plot(Dat$YE_bycatch.kg~Dat$Year, type="l", col="darkgrey",
     ylim=c(0,1.4*max(Dat$YE_bycatch.kg, na.rm=T)),
     ylab="kg", xlab="year")
y0<-pmax(Dat$YE_bycatch.kg-1.68*Dat$bycatch_cv*Dat$YE_bycatch.kg,0)
y1<-Dat$YE_bycatch.kg+1.68*Dat$bycatch_cv*Dat$YE_bycatch.kg
#arrows(x0=Dat$Year, y0=y0, x1=Dat$Year, y1=y1,
#       code=3, col="darkgrey", lwd=1, angle=90, length=0.01)
colp<-adjustcolor("lightblue", alpha.f=0.999)
polygon(c(Dat$Year,rev(Dat$Year)),c(y0,rev(y1)),border=NA, col=colp)

y0<-pmax(Dat$YE_discards.kg-1.68*Dat$YE_discards_cv*Dat$YE_discards.kg,0)
y1<-Dat$YE_discards.kg+1.68*Dat$YE_discards_cv*Dat$YE_discards.kg
colp<-adjustcolor("purple", alpha.f=0.25)
polygon(c(Dat$Year,rev(Dat$Year)),c(y0,rev(y1)),border=NA, col=colp)
points(Dat$YE_discards.kg~Dat$Year, type="l", col="purple")

#arrows(x0=Dat$Year, y0=y0, x1=Dat$Year, y1=y1,
#       code=3, col="darkcyan", lwd=1, angle=90, length=0.01)

abline(v=1998,col="red", lty=2)
points(Dat$YE_discards.kg~Dat$Year, type="l", col="purple", lwd=2)
points(Dat$YE_bycatch.kg~Dat$Year, type="l", main=d, col="blue", lwd=2)
points(Dat$YE.landed.bycatch~Dat$Year, type="l", col="darkgreen", lty=3, lwd=2)

legend(x="topleft",c("Expected bycatch","Landed bycatch","Estimated discards"),pch=c(18,18,18),
       col=c("blue","darkgreen","purple"),bty="n",text.col=c("blue","darkgreen","purple"),
       cex=1)
#points(Dat$truncYE.discards.kg~Dat$Year, type="l", col="purple")
#y0<-Dat$truncYE.discards.kg-1.68*Dat$cv.truncDisc*Dat$truncYE.discards.kg
#y1<-Dat$truncYE.discards.kg+1.68*Dat$cv.truncDisc*Dat$truncYE.discards.kg
#arrows(x0=Dat$Year, y0=y0, x1=Dat$Year, y1=y1,
#       code=3, col="purple", lwd=1, angle=90, length=0.03)
}

#compare to known harvests for perspective...
Removals4<-read.csv("Data/SEO_YE_removals4.csv")  #YE reconstruction from Rhea
Removals4$Remove.kg<-Removals4$Remove.lbs*0.45359237      #convert to kilograms
Removals4$Year
Removals4$YE.lnd.hal.bycatch<-Removals4$YE_fr_hal_fish.lbs*0.45359237

{
  par(mfrow=c(1,1))
  Dat<-SEO.YEbycatch
  plot(Removals4$Remove.kg~Removals4$Year, type="l", col="Black",
       ylim=c(0,max(Removals4$Remove.kg, na.rm=T)),
       ylab="kg", xlab="year", lwd=2)
  
  y0<-pmax(Dat$YE_bycatch.kg-1.68*Dat$bycatch_cv*Dat$YE_bycatch.kg,0)
  y1<-Dat$YE_bycatch.kg+1.68*Dat$bycatch_cv*Dat$YE_bycatch.kg
  
  colp<-adjustcolor("lightblue", alpha.f=0.5)
  polygon(c(Dat$Year,rev(Dat$Year)),c(y0,rev(y1)),border=NA, col=colp)
  
  y0<-pmax(Dat$YE_discards.kg-1.68*Dat$YE_discards_cv*Dat$YE_discards.kg,0)
  y1<-Dat$YE_discards.kg+1.68*Dat$YE_discards_cv*Dat$YE_discards.kg
  colp<-adjustcolor("purple", alpha.f=0.25)
  polygon(c(Dat$Year,rev(Dat$Year)),c(y0,rev(y1)),border=NA, col=colp)
  points(Dat$YE_discards.kg~Dat$Year, type="l", col="purple")
  
#  abline(v=1998,col="red", lty=2)
  points(Dat$YE_discards.kg~Dat$Year, type="l", col="purple", lwd=2)
  points(Dat$YE_bycatch.kg~Dat$Year, type="l", main=d, col="blue", lwd=2)
  points(Dat$YE.landed.bycatch~Dat$Year, type="l", col="darkgreen", lty=3, lwd=2)
  
  legend(x="topright",c("Total know catches","Expected bycatch","Landed bycatch","Estimated discards"),pch=c(18,18,18,18),
         col=c("black","blue","darkgreen","purple"),bty="n",text.col=c("black","blue","darkgreen","purple"),
         cex=1)
}

#look at bycatch rates for SEO as a whole... 
{
  Dat<-SEO.YEbycatch[SEO.YEbycatch$Year > 1979,]
  plot(Dat$mean.ratio~Dat$Year, col="goldenrod4", pch=17, ylim=c(0,0.2),
       ylab="WCPUE / Bycatch rate", xlab="Year")
  y0<-pmax(Dat$mean.ratio-1.68*Dat$ratio.cv*Dat$mean.ratio,0)
  y1<-Dat$mean.ratio+1.68*Dat$ratio.cv*Dat$mean.ratio
  colp<-adjustcolor("goldenrod", alpha.f=0.5)
  polygon(c(Dat$Year,rev(Dat$Year)),c(y0,rev(y1)),border=NA, col=colp)
  points(Obs.Rates$YE.Bycatch.rate ~ Obs.Rates$Year, col="forestgreen", pch=18)
  y0<-pmax(Obs.Rates$YE.Bycatch.rate-1.68*Obs.Rates$YE.Byrate.boot.cv*Obs.Rates$YE.Bycatch.rate,0)
  y1<-Obs.Rates$YE.Bycatch.rate+1.68*Obs.Rates$YE.Byrate.boot.cv*Obs.Rates$YE.Bycatch.rate
  arrows(x0=Obs.Rates$Year, y0=y0, x1=Obs.Rates$Year, y1=y1,
         code=3, col="forestgreen", lwd=1.5, angle=90, length=0.01)
  legend(x="topleft",c("WCPUE (IPHC survey)","Observed bycatch rate (NOAA)"),pch=c(17,18,18),
           col=c("goldenrod4","forestgreen"),bty="n",text.col=c("goldenrod4","forestgreen"),
           cex=1.3)
}

##**********************************************************************
##SEO.YEbycatch is the data for modelling...
#***********************************************************************
write.csv(SEO.YEbycatch,"Data/SEO_YE_discards.csv")

##***********************************************************************
## IPHC CPUE data for pop modelling
##***********************************************************************
str(YE.IPHC.CPUE)
SEO.CPUE<-YE.IPHC.CPUE[YE.IPHC.CPUE$SEdist == "SSEO" | YE.IPHC.CPUE$SEdist == "CSEO" |
                         YE.IPHC.CPUE$SEdist == "NSEO" | YE.IPHC.CPUE$SEdist == "EYKT",]
Ys<-unique(SEO.CPUE$Year)
IPHC.YE.CPUE<-data.frame()
i<-1
for (y in Ys){  #y<-hyears[10]
  Dat<-SEO.CPUE[SEO.CPUE$Year == y,]
  IPHC.YE.CPUE[i,"Year"] <- y
  IPHC.YE.CPUE[i,"mean.CPUE"]<-mean(Dat$YE_CPUE)
  IPHC.YE.CPUE[i,"var.CPUE"] <- sum(Dat$YE_CPUE.var)/16   #ignores covariance term.. need back in the boot?
  IPHC.YE.CPUE[i,"cv.CPUE"]<-sqrt(IPHC.YE.CPUE[i,"var.CPUE"])/mean(Dat$YE_CPUE)
  i<-i+1
}

write.csv(YE.IPHC.CPUE, "Data/SEsubd_IPHC_YE_CPUE.csv")
write.csv(IPHC.YE.CPUE,"Data/SEO_IPHC_YE_CPUE.csv")
#=====================================================================================
## Bootstrap for weights of YE and weights of HA so bayesian state-space production model can do the hindcasting...
#===========================================================================================
YE_HA_weights<-data.frame()
YE.IPHC.CPUE<-data.frame()
Subs<-unique(Survey$SEdist)
Years<-unique(Survey$Year.x)
nboot<-1000

j<-1
for (y in Years) {  #y<-Years[1]
  for (s in Subs){  #s<-Subs[4]
    Dat<-Survey[Survey$Year.x == y & Survey$SEdist == s,]
    Stations<-unique(Dat$Station.x)
    KgYEi<-vector()
    log.KgYEi<-vector()
    KgHai<-vector()
    log.KgHai<-vector()
    i<-1
    for (st in Stations){    #i<-Stations[9]   length(Stations)
      Stat.Dat<-Dat[Dat$Station.x == st,]; Stat.Dat
      #debug
      #if (nrow(Stat.Dat) > 1){aaa} else {}
      CPUE<-mean(Stat.Dat$Number.Observed/Stat.Dat$HooksObserved)
      if (CPUE == 0){
        C<-0
      } else {
        C<-CPUE*Stat.Dat$HooksRetrieved
      }
      
      Kg.YE<-C*Dat$mean.YE.kg[1]
      log.KgYE<-log(Kg.YE+1)
      
      Kg.Ha<-mean(Stat.Dat$O32.Pacific.halibut.weight)   #just taking the mean here because sometimes duplicate data 
      log.KgHA<-log(Kg.Ha+1)
      
      KgYEi[i]<-Kg.YE
      log.KgYEi[i]<-log.KgYE
      KgHai[i]<-Kg.Ha
      log.KgHai[i]<-log.KgHA
      i<-i+1
    }
    #bootstrap WCPUE
    Out<-data.frame()
    CPUE.out<-data.frame()
    for (ii in 1:nboot){ #nboot<-1000  i<-1
      Resamp<-sample(KgYEi,length(KgYEi),replace=T)
      Out[ii,"KgYE"]<-mean(Resamp, na.rm=T)
      Out[ii,"var.KgYE"]<-var(Resamp, na.rm=T)
      Resamp2<-sample(KgHai,length(KgHai),replace=T)
      Out[ii,"KgHa"]<-mean(Resamp2, na.rm=T)
      Out[ii,"var.KgHa"]<-var(Resamp2, na.rm=T)
      Resamp3<-sample(log.KgYEi,length(log.KgYEi),replace=T)
      Out[ii,"log.KgYE"]<-mean(Resamp3, na.rm=T)
      Out[ii,"var.log.KgYE"]<-var(Resamp3, na.rm=T)
      Resamp4<-sample(log.KgHai,length(log.KgHai),replace=T)
      Out[ii,"log.KgHa"]<-mean(Resamp4, na.rm=T)
      Out[ii,"var.log.KgHa"]<-var(Resamp4, na.rm=T)
    }
    
    hist(Out$log.KgYE, breaks = 25)
    YE_HA_weights[j,"Year"]<-y
    YE_HA_weights[j,"SEdist"]<-s
    YE_HA_weights[j,"Kg.YE"]<-unname(quantile(Out$KgYE,c(0.5)))  #mean WCPUE from Tribuzio
    YE_HA_weights[j,"var.KgYE"]<-unname(quantile(Out$var.KgYE,c(0.5)))
    YE_HA_weights[j,"Kg.Ha"]<-unname(quantile(Out$KgHa,c(0.5)))  #mean WCPUE from Tribuzio
    YE_HA_weights[j,"var.KgHa"]<-unname(quantile(Out$var.KgHa,c(0.5)))
    
    YE_HA_weights[j,"log.KgYE"]<-unname(quantile(Out$log.KgYE,c(0.5)))  #mean WCPUE from Tribuzio
    YE_HA_weights[j,"var.log.KgYE"]<-unname(quantile(Out$var.log.KgYE,c(0.5)))
    YE_HA_weights[j,"log.KgHa"]<-unname(quantile(Out$log.KgHa,c(0.5)))  #mean WCPUE from Tribuzio
    YE_HA_weights[j,"var.log.KgHa"]<-unname(quantile(Out$var.log.KgHa,c(0.5)))
    
    j<-j+1
  }
}
YE_HA_weights[YE_HA_weights$Year == 2020,]

SEO.YE_HA_weights<-YE_HA_weights[YE_HA_weights$SEdist == "SSEO" | YE_HA_weights$SEdist == "CSEO" |
                                   YE_HA_weights$SEdist == "NSEO" | YE_HA_weights$SEdist == "EYKT",]
years<-unique(SEO.YE_HA_weights$Year)
SEO.IPHC.weights<-data.frame()
i<-1
for (y in years){  #y<-years[1]
  Dat<-SEO.YE_HA_weights[SEO.YE_HA_weights$Year == y,]
  SEO.IPHC.weights[i,"Year"] <- y
  SEO.IPHC.weights[i,"Kg.YE"]<-mean(Dat$Kg.YE)
  SEO.IPHC.weights[i,"var.KgYE"]<-sum(Dat$var.KgYE)
  SEO.IPHC.weights[i,"cv.KgYE"]<-sqrt(SEO.IPHC.weights[i,"var.KgYE"])/SEO.IPHC.weights[i,"Kg.YE"]
  
  SEO.IPHC.weights[i,"Kg.Ha"]<-mean(Dat$Kg.Ha)
  SEO.IPHC.weights[i,"var.KgHa"]<-sum(Dat$var.KgHa)
  SEO.IPHC.weights[i,"cv.KgHa"]<-sqrt(SEO.IPHC.weights[i,"var.KgHa"])/SEO.IPHC.weights[i,"Kg.Ha"]
  
  SEO.IPHC.weights[i,"log.KgYE"]<-mean(Dat$log.KgYE)
  SEO.IPHC.weights[i,"var.log.KgYE"]<-sum(Dat$var.log.KgYE)
  SEO.IPHC.weights[i,"cv.log.KgYE"]<-sqrt(SEO.IPHC.weights[i,"var.log.KgYE"])/SEO.IPHC.weights[i,"log.KgYE"]
  
  SEO.IPHC.weights[i,"log.KgHa"]<-mean(Dat$log.KgHa)
  SEO.IPHC.weights[i,"var.log.KgHa"]<-sum(Dat$var.log.KgHa)
  SEO.IPHC.weights[i,"cv.log.KgHa"]<-sqrt(SEO.IPHC.weights[i,"var.log.KgHa"])/SEO.IPHC.weights[i,"log.KgHa"]
  
  i<-i+1
}


write.csv(SEO.IPHC.weights,"Data/SEO_IPHC_Surv_Weights.csv")
