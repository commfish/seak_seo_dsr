#################################################################################
## Calculation of expected bycatch of YE in Halibut fishery using wcpue estimates
## from Tribuzio et al. 2014.  E
## Expected bycatch to be used in SS-SPM YE models for SEO and 2022 SAFE
## June 2022
## Phil Joy
##
## Background stuff on these calculations in Discards.Rproj in Yelloweye/Unreported Discards
#################################################################################
library(dplyr)
library(boot)
library(ggplot2)
library(RColorBrewer)
library(matrixStats)

source("Code/Port_bio_function.R")

#get processed IPHC survey data built in IPHC_Survey_CPUE_index.R
Survey<-read.csv(paste0("Data_processing/Data/IPHC_survey_1998-",YEAR,".csv"))

str(Survey); head(Survey,10)
str(HA.Harv)
#================================================================================
# Look at some depth stuff... 
hist(Survey$BeginDepth..fm.)
hist(Survey$EndDepth..fm.)
hist(Survey$AvgDepth.fm)

plot(Survey$YE.obs ~ Survey$AvgDepth.fm)
plot(Survey$O32.Pacific.halibut.count ~ Survey$AvgDepth.fm)
plot(Survey$U32.Pacific.halibut.count ~ Survey$AvgDepth.fm)

#==============================================================================================
## Load the port samples that were downloaded from oceansAK; 
# Need Yelloweye weights to get wcpue estimates
{  Port<-port.bio(2022)
  
  unique(Port$Sample.Type)
  Port.rand<-Port[Port$Sample.Type=="Random",]
  Port.direct<-Port[Port$Project=="Commercial Longline Trip",]
  Port.hal<-Port[Port$Project=="Commercial Halibut Longline",]; nrow(Port.hal)
  
  Survey$mean.YE.kg<-NA
  Survey$var.YE.kg<-NA
  Survey$N.YE.kg<-NA
  
  
  ## Get best data for average weights and add columns to survey data
  Years<-sort(unique(Survey$Year))
  GFMA<-unique(Survey$SEdist)
  
  for (i in Years){    #i<-Years[1]
    for (j in GFMA){   #j<-GFMA[1]
      if (i < 2008){
        P <- Port.rand %>% filter (GFMU == j); nrow(P)
        P <- P %>% filter(!is.na(Weight.Kilograms)); nrow(P)
        #reach back from year i and see if you can get 150 weights...nrow(P)
        Nw<-nrow(P[P$Year == i  & !is.na(P$Weight.Kilograms),])
        Nw<-nrow(P %>% filter(Year == i, !is.na(Weight.Kilograms)))
        #reach back from year i and see if you can get 150 weights...
        Nw<-nrow(P[P$Year == i,])
        k<-0
        while (Nw < 75 & i-k >= min(Port.rand$Year) ){         #go back as far as needed to get 150 weights... 
          k<-k+1
          Nw<-nrow(P[P$Year >=i-k & P$Year <= i,])   #P[P$Year == 1990:2002,]
        }
        m<-0
        while (Nw < 75 & i+m <= max(Port.rand$Year)) {     #go forward if you failed to find enough weights... 
          m<-m+1                                                        
          Nw<-nrow(P[P$Year >=i-k & P$Year <= i+m,])   #P[P$Year == 1986,]
        }
        #get average weights of YE...
        Sample<-P$Weight.Kilograms[P$Year >=i-k & P$Year <= i+m]
        Survey[,"mean.YE.kg"][Survey$Year == i & Survey$SEdist == j]<-mean(Sample)
        Survey[,"var.YE.kg"][Survey$Year == i & Survey$SEdist == j]<-var(Sample)
        Survey[,"N.YE.kg"][Survey$Year == i & Survey$SEdist == j]<-length(Sample)
      } else {
        P <- Port.hal %>% filter (GFMU == j); nrow(P)
        P <- P %>% filter(!is.na(Weight.Kilograms)); nrow(P)
        #reach back from year i and see if you can get 150 weights...nrow(P)
        Nw<-nrow(P[P$Year == i  & !is.na(P$Weight.Kilograms),])
        Nw<-nrow(P %>% filter(Year == i, !is.na(Weight.Kilograms)))
        #reach back from year i and see if you can get 150 weights...
        Nw<-nrow(P[P$Year == i,])
        k<-0
        while (Nw < 75 & i-k >= min(Port.hal$Year) ){         #go back as far as needed to get 150 weights... 
          k<-k+1
          Nw<-nrow(P[P$Year >=i-k & P$Year <= i,])   #P[P$Year == 1990:2002,]
        }
        m<-0
        while (Nw < 75 & i+m <= max(Port.hal$Year)) {     #go forward if you failed to find enough weights... 
          m<-m+1                                                        
          Nw<-nrow(P[P$Year >=i-k & P$Year <= i+m,])   #P[P$Year == 1986,]
        }
        #get average weights of YE...
        Sample<-P$Weight.Kilograms[P$Year >=i-k & P$Year <= i+m]
        Survey[,"mean.YE.kg"][Survey$Year == i & Survey$SEdist == j]<-mean(Sample)
        Survey[,"var.YE.kg"][Survey$Year == i & Survey$SEdist == j]<-var(Sample)
        Survey[,"N.YE.kg"][Survey$Year == i & Survey$SEdist == j]<-length(Sample)
      }
    }
  }
}
Survey$mean.YE.kg

#=============================================================================
# Load halibut harvests by subdistrict to calculate expected bycatch values;
#Deal with Rhea's subdistrict data (halibut_catch_data.csv) and Randy's halibut harvest
# reconstruction (halibut_SEO_catch_RPderivation.csv)
SubHal<-read.csv("Data/Harvests/halibut_catch_data.csv", header=T)
SubHal2<-read.csv("Data/SE_Halibut_removals_1975-2022.csv", header=T)

str(SubHal)
str(SubHal2)

SEOHal<-read.csv("Data/SEO_Halibut_removals_1888-2021.csv", header=T)

#notes: look at Rhea and Randy' data to check for discrepencies... contact Rhea if
# there is something.  
#Make sure Halibut data is formatted the same for wcpue and bycatch estimation 
#function below: H.catch<-Hali[Hali$SEdist == s & Hali$Year == y,]

Halibut<-function(){
 # HA.Harv<-read.csv("Data/Harvests/halibut_catch_data.csv", header=T)
  HA.Harv<-read.csv("Data/SE_Halibut_removals_1975-2022.csv", header=T)
  LL<-HA.Harv #[HA.Harv$gear.description == "Longline",]
  head(LL)
  colnames(LL)
  
  LL<-HA.Harv
  #LL<-rename(LL, SEdist = groundfish.mgt.area.district)
  LL<-rename(LL, SEdist = Mgt.Area)
  #LL<-rename(LL, Year = year.landed)
  Hali<-LL; head(LL)
  
  plot(LL$Year,LL$round_lbs_hali)
  plot(LL$Year,LL$HA.mt)
  Hcheck<-data.frame()
  Yrs<-unique(LL$Year)
  i<-1
  for (y in Yrs){  #y<-Yrs[1]
    Dat<-LL[LL$Year == y,]; Dat<-Dat[Dat$SEdist=="SSEO"|Dat$SEdist=="CSEO"|Dat$SEdist=="NSEO"|Dat$SEdist=="EYKT",]
    Hcheck[i,"Year"]<-y; Hcheck[i,"hal"]<-sum(Dat$HA.mt); i<-i+1
  }
  plot(Hcheck$Year, Hcheck$hal)
  return(Hali)
}
Hali<-Halibut()
str(Hali)
#=============================================================================
# 
#to gauge relationship between depth and wcpue estimates calculate for different depth 
# bins and compare... 

#*** IFQ started in 1995; prior to that it's Derby style ***

IPHC.wcpue<-data.frame()
Subs<-unique(Survey$SEdist)
Years<-unique(Survey$Year)
nboot<-1000

j<-1
for (y in Years) {  #y<-Years[1]
  for (d in Depths){  #d<-Depths[2]
    Dat<-Survey[Survey$Year == y & Survey$depth_bin == d,]
    if (nrow(Dat)>0){
      Stations<-unique(Dat$Station.x)
      WCPUEi.32<-vector()
      WCPUEi.all<-vector()
      CPUEi<-vector()
      i<-1
      for (st in Stations){    #st<-Stations[1]   length(Stations)  st<-1
        Stat.Dat<-Dat[Dat$Station.x == st,] #; Stat.Dat
        #debug
        #if (nrow(Stat.Dat) > 1){aaa} else {}
        CPUE<-mean(Stat.Dat$YE.obs/Stat.Dat$HooksObserved)
        if (CPUE == 0){
          C<-0
        } else {
          C<-CPUE*Stat.Dat$HooksRetrieved
        }
        
        CPUEi[i]<-CPUE
        K<-C*Dat$mean.YE.kg[1]
        WCPUE.32 <-K/mean(Stat.Dat$O32.Pacific.halibut.weight)
        WCPUE.all<-K/mean(Stat.Dat$O32.Pacific.halibut.weight+Stat.Dat$U32.Pacific.halibut.weight)
        WCPUEi.32[i]<-WCPUE.32
        WCPUEi.all[i]<-WCPUE.all
        i<-i+1
      }
      #get rid of NA's where no halibut were caught!!! 
      WCPUEi.32<-WCPUEi.32[!is.na(WCPUEi.32)]
      WCPUEi.all<-WCPUEi.all[!is.na(WCPUEi.all)]
      #bootstrap WCPUE
      Out<-data.frame()
      CPUE.out<-data.frame()
      for (ii in 1:nboot){ #nboot<-1000  ii<-1
        Resamp<-sample(WCPUEi.32,length(WCPUEi.32),replace=T)
        Out[ii,"WCPUE.32"]<-mean(Resamp, na.rm=T)
        Out[ii,"var.WCPUE.32"]<-var(Resamp, na.rm=T)
        Resamp2<-sample(WCPUEi.all,length(WCPUEi.all),replace=T)
        Out[ii,"WCPUE.all"]<-mean(Resamp2, na.rm=T)
        Out[ii,"var.WCPUE.all"]<-var(Resamp2, na.rm=T)
        Resamp3<-sample(CPUEi,length(CPUEi),replace=T)
        CPUE.out[ii,"CPUE"]<-mean(Resamp3, na.rm=T)
        CPUE.out[ii,"CPUE.var"]<-var(Resamp3, na.rm=T)
        #bootstrap bycatch and discards 12-20-21: exact same results if you pu
        #Out[ii,"YE.bycatch.kg"]<-sum(Harv$round_lbs_hali)*0.45359237*Out[ii,"WCPUE"]
        #Out[ii,"raw.YE.discards.kg"]<-Out[ii,"YE.bycatch.kg"]-Land.YEby$YE_fr_hal_fish.kgs #max(1,Out[ii,"YE.bycatch.kg"]-Land.YEby$YE_fr_hal_fish.kgs, na.rm=T)
        #min.D.boot<-min(1,rnorm(1,mean=min.Discard ,sd=sqrt(min.D.var))*Land.YEby$YE_fr_hal_fish.kgs)
        #Out[ii,"trunc.YE.discards.kg"]<-max(min.D.boot,Out[ii,"YE.bycatch.kg"]-Land.YEby$YE_fr_hal_fish.kgs, na.rm=T)
      }
      
      #hist(Out$KgHa, breaks = 25)
      
      for (k in 1:2){
        if (k == 1){
          IPHC.wcpue[j,"Year"]<-y
          IPHC.wcpue[j,"depth_strat"]<-d
          IPHC.wcpue[j,"length.cat"]<-32
          IPHC.wcpue[j,"WCPUE.mean"]<-mean(WCPUEi.32)
          IPHC.wcpue[j,"WCPUE.bootmean"]<-unname(quantile(Out$WCPUE.32,c(0.5)))  #mean WCPUE from Tribuzio
          IPHC.wcpue[j,"WCPUE.lo95ci"]<-unname(quantile(Out$WCPUE.32,c(0.025)))
          IPHC.wcpue[j,"WCPUE.hi95ci"]<-unname(quantile(Out$WCPUE.32,c(0.975)))
          IPHC.wcpue[j,"WCPUE.var"]<-var(Out$WCPUE.32)
          IPHC.wcpue[j,"WCPUE.cv"]<-sd(Out$WCPUE.32)/mean(WCPUEi.32)
          j<-j+1
        } else {
          IPHC.wcpue[j,"Year"]<-y
          IPHC.wcpue[j,"depth_strat"]<-d
          IPHC.wcpue[j,"length.cat"]<-"all"
          IPHC.wcpue[j,"WCPUE.mean"]<-mean(WCPUEi.all)
          IPHC.wcpue[j,"WCPUE.bootmean"]<-unname(quantile(Out$WCPUE.all,c(0.5)))  #mean WCPUE from Tribuzio
          IPHC.wcpue[j,"WCPUE.lo95ci"]<-unname(quantile(Out$WCPUE.all,c(0.025)))
          IPHC.wcpue[j,"WCPUE.hi95ci"]<-unname(quantile(Out$WCPUE.all,c(0.975)))
          IPHC.wcpue[j,"WCPUE.var"]<-var(Out$WCPUE.all)
          IPHC.wcpue[j,"WCPUE.cv"]<-sd(Out$WCPUE.all)/mean(WCPUEi.all)
          j<-j+1
        }
      }
      
      #j<-j+1
    } else {}
  }
}

# Does wcpue vary by depth?
IPHC.wcpue[IPHC.wcpue$Year == 2020,]
str(IPHC.wcpue)
IPHC.wcpue$depth_strat<-as.factor(as.numeric(IPHC.wcpue$depth_strat))

par(mfrow=c(1,1))
plot(IPHC.wcpue$WCPUE.mean[IPHC.wcpue$length.cat == 32]~
       IPHC.wcpue$depth_strat[IPHC.wcpue$length.cat == 32], 
     xlab="Depth strata (fathoms)", 
     ylab="WCPUE.32", col="red",
     type="p", pch=18)
points(IPHC.wcpue$WCPUE.mean[IPHC.wcpue$length.cat == "all"]~
         IPHC.wcpue$depth_strat[IPHC.wcpue$length.cat == "all"], 
       xlab="Depth strata (fathoms)", 
       ylab="WCPUE.all", col="blue",
       type="p", pch=18)
simp<-lm(IPHC.wcpue$WCPUE32.mean~IPHC.wcpue$depth_strat*IPHC.wcpue$Year)
summary(simp)

ggplot(IPHC.wcpue,aes(depth_strat,WCPUE.mean, fill = length.cat))+geom_boxplot()

{
  cols<-brewer.pal(length(Depths),"Set2")
  YE.depths<-Depths[1:5]
  par(mfrow=c(2,1))
  i<-1
  for (d in YE.depths){  #d<-Depths[1]
    Dat32<-IPHC.wcpue[IPHC.wcpue$depth_strat == d & IPHC.wcpue$length.cat==32,]
    
    if (d == Depths[1]){
      plot(Dat32$WCPUE.mean ~ Dat32$Year, type="l", col=cols[1], ylim=c(0,0.25), 
           lwd=2, ylab="WCPUE.32", xlab="Year")
      Y1<-Dat32$WCPUE.lo95ci; Y2<-Dat32$WCPUE.hi95ci
      polygon(c(Years,rev(Years)),c(Y1,rev(Y2)),
              col=adjustcolor(cols[i],alpha.f=0.2),border=NA)
      i<-i+1
    } else {
      lines(Dat32$WCPUE.mean ~ Dat32$Year, col=cols[i], lwd=2)
      Y1<-Dat32$WCPUE.lo95ci; Y2<-Dat32$WCPUE.hi95ci
      polygon(c(Years,rev(Years)),c(Y1,rev(Y2)),
              col=adjustcolor(cols[i],alpha.f=0.2),border=NA)
      i<-i+1
    }
  }
  legend(x="topright", title="Depth Strata",
         legend=c(YE.depths), bty="n",
         col=cols,text.col=cols)
  
  i<-1
  for (d in YE.depths){  #d<-Depths[6]
    Datall<-IPHC.wcpue[IPHC.wcpue$depth_strat == d & IPHC.wcpue$length.cat=="all",]
    if (d == Depths[1]){
      plot(Datall$WCPUE.mean ~ Datall$Year, type="l", col=cols[1], ylim=c(0,0.25), 
           lwd=2, ylab="WCPUE.all", xlab="Year")
      Y1<-Datall$WCPUE.lo95ci; Y2<-Datall$WCPUE.hi95ci
      polygon(c(Years,rev(Years)),c(Y1,rev(Y2)),
              col=adjustcolor(cols[i],alpha.f=0.2),border=NA)
      i<-i+1
    } else {
      lines(Datall$WCPUE.mean ~ Datall$Year, col=cols[i], lwd=2)
      Y1<-Datall$WCPUE.lo95ci; Y2<-Datall$WCPUE.hi95ci
      polygon(c(Years,rev(Years)),c(Y1,rev(Y2)),
              col=adjustcolor(cols[i],alpha.f=0.2),border=NA)
      i<-i+1
    }
  }
  legend(x="topright", title="Depth Strata",
         legend=c(YE.depths), bty="n",
         col=cols,text.col=cols)
}
#I suppose we could do some fancy stats, but not really necessary.
#if so, let's get depth profile of halibut fishery
#**************************************************
#*Need to get halibut fishery depth data!!!! 
#*For now below will get by area and depth strata but we can dial this in
#*once we have that data! 

#==============================================================================
# Now, lets produce wcpue and expected bycatch for all relevant management units for all depths combined
# as well as for each depth bin so we can deal with hindcasting for the variety
# of spatial resolutions we are likely to see

str(Hali)
str(SEOHal)
str(Survey)
#need by $IPHC.Reg.Area
#need by $IPHC.Stat.Area
#need by $IPHC.Charter.Region
#need by $SEdist
#need by #In.Out

#create a function to select the management area and the depth brackets 
# to calculate WCPUE and expected bycatch.

YEHA.fxn<-function(Survey=Survey,Area="SEdist",Deep=max(Survey$AvgDepth.fm), Shallow=0, 
                   nboot=1000, Hali=Hali){   #specify which halibut data set to use!!!
  col.ref<-which(colnames(Survey)==Area)
  
  IPHC.expBy<-data.frame()
  Subs<-unique(Survey[,col.ref])
  Years<-unique(Survey$Year)
  
  
  
  j<-1
  for (y in Years) {  #y<-Years[24]
    for (s in Subs){  #s<-Subs[2]
      Dat<-Survey[Survey$Year == y & Survey[,col.ref] == s &
                    Survey$AvgDepth.fm > Shallow & Survey$AvgDepth.fm < Deep,]
      if (Area == "SEdist") {
        H.catch<-Hali[Hali$SEdist == s & Hali$Year == y,]
      }
      if (Area == "In.Out") {
        H.catch <- Hali[Hali$Year == y,]
        H.catch$HA.mt<-H.catch$SEO_harvest_mt
      }
      
      if (nrow(Dat)>0){
        Stations<-unique(Dat$Station.x)
        WCPUEi.32<-vector()
        WCPUEi.all<-vector()
        expByi.32<-vector()
        expByi.all<-vector()
        CPUEi<-vector()
        i<-1
        for (st in Stations){    #st<-Stations[7]   length(Stations)  st<-1
          Stat.Dat<-Dat[Dat$Station.x == st,] #; Stat.Dat
          #debug
          #if (nrow(Stat.Dat) > 1){aaa} else {}
          CPUE<-mean(Stat.Dat$YE.obs/Stat.Dat$HooksObserved)
          if (CPUE == 0){
            C<-0
          } else {
            C<-CPUE*Stat.Dat$HooksRetrieved
          }
          
          CPUEi[i]<-CPUE
          K<-C*Dat$mean.YE.kg[1]
          WCPUE.32 <-K/mean(Stat.Dat$O32.Pacific.halibut.weight)
          WCPUE.all<-K/mean(Stat.Dat$O32.Pacific.halibut.weight+Stat.Dat$U32.Pacific.halibut.weight)
          WCPUEi.32[i]<-WCPUE.32
          WCPUEi.all[i]<-WCPUE.all
          expBy.32<-WCPUE.32*sum(H.catch$HA.mt)
          expBy.all<-WCPUE.all*sum(H.catch$HA.mt)
          expByi.32[i]<-expBy.32
          expByi.all[i]<-expBy.all
          i<-i+1
        }
        #get rid of NA's where no halibut were caught!!! 
        WCPUEi.32<-WCPUEi.32[!is.na(WCPUEi.32)]
        WCPUEi.all<-WCPUEi.all[!is.na(WCPUEi.all)]
        expByi.32<-expByi.32[!is.na(expByi.32)]
        expByi.all<-expByi.all[!is.na(expByi.all)]
        #bootstrap WCPUE
        Out<-data.frame()
        CPUE.out<-data.frame()
        for (ii in 1:nboot){ #nboot<-1000  ii<-1
          Resamp<-sample(WCPUEi.32,length(WCPUEi.32),replace=T)
          Out[ii,"WCPUE.32"]<-mean(Resamp, na.rm=T)
          Out[ii,"var.WCPUE.32"]<-var(Resamp, na.rm=T)
          Resamp2<-sample(WCPUEi.all,length(WCPUEi.all),replace=T)
          Out[ii,"WCPUE.all"]<-mean(Resamp2, na.rm=T)
          Out[ii,"var.WCPUE.all"]<-var(Resamp2, na.rm=T)
          
          Resamp3<-sample(expByi.32,length(expByi.32),replace=T)
          Out[ii,"expBy.32"]<-mean(Resamp3, na.rm=T)
          Out[ii,"var.expBy.32"]<-var(Resamp3, na.rm=T)
          
          Resamp4<-sample(expByi.all,length(expByi.all),replace=T)
          Out[ii,"expBy.all"]<-mean(Resamp4, na.rm=T)
          Out[ii,"var.expBy.all"]<-var(Resamp4, na.rm=T)
          
          Resamp5<-sample(CPUEi,length(CPUEi),replace=T)
          CPUE.out[ii,"CPUE"]<-mean(Resamp5, na.rm=T)
          CPUE.out[ii,"CPUE.var"]<-var(Resamp5, na.rm=T)
          #bootstrap bycatch and discards 12-20-21: exact same results if you pu
          #Out[ii,"YE.bycatch.kg"]<-sum(Harv$round_lbs_hali)*0.45359237*Out[ii,"WCPUE"]
          #Out[ii,"raw.YE.discards.kg"]<-Out[ii,"YE.bycatch.kg"]-Land.YEby$YE_fr_hal_fish.kgs #max(1,Out[ii,"YE.bycatch.kg"]-Land.YEby$YE_fr_hal_fish.kgs, na.rm=T)
          #min.D.boot<-min(1,rnorm(1,mean=min.Discard ,sd=sqrt(min.D.var))*Land.YEby$YE_fr_hal_fish.kgs)
          #Out[ii,"trunc.YE.discards.kg"]<-max(min.D.boot,Out[ii,"YE.bycatch.kg"]-Land.YEby$YE_fr_hal_fish.kgs, na.rm=T)
        }
        
        #hist(Out$KgHa, breaks = 25)
        IPHC.expBy[j,"Year"]<-y
        IPHC.expBy[j,"mngmt.divisions"]<-Area
        IPHC.expBy[j,"mngmt.area"]<-s
        IPHC.expBy[j,"deep.bound"]<-Deep
        IPHC.expBy[j,"shallow.bound"]<-Shallow
        IPHC.expBy[j,"no.stations"]<-length(Stations)
        
        IPHC.expBy[j,"WCPUE32.mean"]<-mean(WCPUEi.32)
        IPHC.expBy[j,"WCPUE32.bootmean"]<-unname(quantile(Out$WCPUE.32,c(0.5)))  #mean WCPUE from Tribuzio
        IPHC.expBy[j,"WCPUE32.lo95ci"]<-unname(quantile(Out$WCPUE.32,c(0.025)))
        IPHC.expBy[j,"WCPUE32.hi95ci"]<-unname(quantile(Out$WCPUE.32,c(0.975)))
        IPHC.expBy[j,"WCPUE32.var"]<-var(Out$WCPUE.32)
        IPHC.expBy[j,"WCPUE32.cv"]<-sd(Out$WCPUE.32)/mean(Out$WCPUE.32)
        #    IPHC.expBy[j,"WCPUE32.cv"]<-sqrt(var(WCPUEi.32))/mean(WCPUEi.32)
        
        IPHC.expBy[j,"WCPUEall.mean"]<-mean(WCPUEi.all)
        IPHC.expBy[j,"WCPUEall.bootmean"]<-unname(quantile(Out$WCPUE.all,c(0.5)))  #mean WCPUE from Tribuzio
        IPHC.expBy[j,"WCPUEall.lo95ci"]<-unname(quantile(Out$WCPUE.all,c(0.025)))
        IPHC.expBy[j,"WCPUEall.hi95ci"]<-unname(quantile(Out$WCPUE.all,c(0.975)))
        IPHC.expBy[j,"WCPUEall.var"]<-var(Out$WCPUE.all)
        IPHC.expBy[j,"WCPUEall.cv"]<-sd(Out$WCPUE.all)/mean(Out$WCPUE.all)
        
        #IPHC.expBy[j,"Halibut_rnd_pnds"]<-sum(H.catch$round_lbs_hali)
        IPHC.expBy[j,"Halibut_mt"]<-sum(H.catch$HA.mt)
        
        IPHC.expBy[j,"expBy32_mt.mean"]<-mean(expByi.32)
        IPHC.expBy[j,"expBy32.bootmean"]<-unname(quantile(Out$expBy.32,c(0.5)))  #mean WCPUE from Tribuzio
        IPHC.expBy[j,"expBy32.lo95ci"]<-unname(quantile(Out$expBy.32,c(0.025)))
        IPHC.expBy[j,"expBy32.hi95ci"]<-unname(quantile(Out$expBy.32,c(0.975)))
        IPHC.expBy[j,"expBy32.var"]<-var(Out$expBy.32)
        IPHC.expBy[j,"expBy32.cv"]<-sd(Out$expBy.32)/mean(Out$expBy.32)
        
        IPHC.expBy[j,"expByall_mt.mean"]<-mean(expByi.32)
        IPHC.expBy[j,"expByall.bootmean"]<-unname(quantile(Out$expBy.32,c(0.5)))  #mean WCPUE from Tribuzio
        IPHC.expBy[j,"expByall.lo95ci"]<-unname(quantile(Out$expBy.32,c(0.025)))
        IPHC.expBy[j,"expByall.hi95ci"]<-unname(quantile(Out$expBy.32,c(0.975)))
        IPHC.expBy[j,"expByall.var"]<-var(Out$expBy.32)
        IPHC.expBy[j,"expByall.cv"]<-sd(Out$expBy.32)/mean(Out$expBy.32)
        
        j<-j+1
      } else {}
    }
  }
  
 # str(IPHC.expBy)
  
  #couple of plots in output
#  Sub.cnt<-length(Subs)
#  cols<-brewer.pal(min(Sub.cnt,8),"Dark2")
  
#  if (Sub.cnt < 8) {par(mfrow=c(1,1))}
#  if (Sub.cnt > 8 & Sub.cnt < 17) {par(mfrow=c(2,1))}
#  if (Sub.cnt > 16) {par(mfrow=c(3,1))}
  
#  nplot<-ceiling(Sub.cnt/8)
  #cols<-brewer.pal(length(Subs),"Set3")
#  for (n in 1:nplot){  #n<-1
#    i<-1
#    pSubs<-Subs[(1+(n-1)*8):(8+(n-1)*8)]
#    pSubs<-pSubs[!is.na(pSubs)]
    
#    for (s in pSubs) {   #s<-pSubs[1]
#      Pdat<-IPHC.expBy[IPHC.expBy$mngmt.area == s,]
#      if (s == pSubs[1]){
#        plot(Pdat$expBy32_rnd_pds.mean ~ Pdat$Year, col=cols[i], ylim=c(0,max(IPHC.expBy$expBy32.hi95ci)), type="l",
#             ylab="expected bycatch", xlab="Year", main=Area)
#        Y1<-Pdat$expBy32.lo95ci; Y2<-Pdat$expBy32.hi95ci
#        polygon(c(Years,rev(Years)),c(Y1,rev(Y2)),
#                col=adjustcolor(cols[i],alpha.f=0.2),border=NA)
#        i<-i+1
#      } else {
#        Y1<-Pdat$expBy32.lo95ci; Y2<-Pdat$expBy32.hi95ci
#        if (length(Y1)==1){}else {
#          polygon(c(Years,rev(Years)),c(Y1,rev(Y2)),
#                  col=adjustcolor(cols[i],alpha.f=0.2),border=NA)
#          i<-i+1
#        }
        
#      }
#    }
#    i<-1
#    for (s in pSubs) {   #s<-pSubs[1]
#      Pdat<-IPHC.expBy[IPHC.expBy$mngmt.area == s,]
#      lines(Pdat$expBy32_rnd_pds.mean ~ Pdat$Year, col=cols[i], type="l",
#            lwd=1.5)
#      i<-i+1
#    }
#   legend(x="topleft", cex=0.8,
#           legend=c(pSubs), bty="n",
   #        col=cols,text.col=cols)
  #}
  
  return(IPHC.expBy)
}


colnames(Survey)
#IPHC.reg.area<- YEHA.fxn(Area="IPHC.Reg.Area",Deep=max(Survey$AvgDepth.fm), Shallow=0, nboot=1000) 
#IPHC.stat.area<- YEHA.fxn(Area="IPHC.Stat.Area",Deep=max(Survey$AvgDepth.fm), Shallow=0, nboot=1000) 
SE.subdistricts<- YEHA.fxn(Survey=Survey,Area="SEdist",Deep=max(Survey$AvgDepth.fm), Shallow=0, Hali=Hali,nboot=1000)
SE.out<- YEHA.fxn(Survey=Survey[Survey$In.Out == "SEO",], Hali = SEOHal,
                         Area="In.Out",Deep=max(Survey$AvgDepth.fm), Shallow=0, nboot=1000)

#get just outside districts
SEO.subs<-SE.subdistricts[SE.subdistricts$mngmt.area != "NSEI" & SE.subdistricts$mngmt.area != "SSEI",]
unique(SEO.subs$mngmt.area)
SEO.subs[SEO.subs$mngmt.area == "NSEO",]
SEO.subs[SEO.subs$mngmt.area == "EYKT",]
#================================================================================
# Calculate expected bycatch for years before 1998; pre-survey years
#inflate pre-1995 (Derby years) by 0.5 to account for further uncertainty...
# that is not done here, but will be done in the DATALOAD for the model (easier to change there 
# with model development...)

pres.yrs<-seq(1980,1997,by=1)
subs<-unique(SEO.subs$mngmt.area)

presurvey.bycatch<-data.frame()
i<-1
for (y in pres.yrs) { #y<-pres.yrs[1]
  for (s in subs) {  #s<-subs[1]
    ma<-SEO.subs[SEO.subs$mngmt.area ==s,]
    hal<-Hali[Hali$Year == y & Hali$SEdist ==s,]
    H.catch<-sum(hal$HA.mt)
    mean.wcpue<-mean(ma$WCPUE32.mean)
    var.wcpue<-sum(ma$WCPUE32.var)/(length(unique(ma$Year))^2)
    
    mean.wcpueall<-mean(ma$WCPUEall.mean)
    var.wcpueall<-sum(ma$WCPUEall.var)/(length(unique(ma$Year))^2)
    
    presurvey.bycatch[i,"Year"]<-y
    presurvey.bycatch[i,"mngmt.area"]<-s
    presurvey.bycatch[i,"Halibut_mt"]<-H.catch
    presurvey.bycatch[i,"WCPUE32.mean"]<-mean.wcpue
    presurvey.bycatch[i,"WCPUE32.var"]<-var.wcpue
    presurvey.bycatch[i,"expBy32_mt.mean"]<-mean.wcpue*H.catch
    presurvey.bycatch[i,"expBy32.var"]<-(H.catch^2)*var.wcpue
    presurvey.bycatch[i,"expBy32.cv"]<-sqrt(presurvey.bycatch[i,"expBy32.var"])/
                                        presurvey.bycatch[i,"expBy32_mt.mean"]
    #want to use inflated cv for historical so also save max cv from survey years... 
    presurvey.bycatch[i,"expBy32.cv2"]<-max(ma$expBy32.cv)
    
    presurvey.bycatch[i,"WCPUEall.mean"]<-mean.wcpueall
    presurvey.bycatch[i,"WCPUEall.var"]<-var.wcpueall
    presurvey.bycatch[i,"expByall_mt.mean"]<-mean.wcpueall*H.catch
    presurvey.bycatch[i,"expByall.var"]<-(H.catch^2)*var.wcpueall
    presurvey.bycatch[i,"expByall.cv"]<-sqrt(presurvey.bycatch[i,"expByall.var"])/
      presurvey.bycatch[i,"expByall_mt.mean"]
    #want to use inflated cv for historical so also save max cv from survey years... 
    presurvey.bycatch[i,"expByall.cv2"]<-max(ma$expByall.cv)
    i<-i+1
  }
}

#weird in my head that variance stays the same.. but is correct because variance of bycatch is halibut times
# same variance calculated from averages... so uncertainty a reflection of that... 
# cv's of bycatch are small - averaging reduces the uncertainty weirdly
# may just objectively apply larger cv's to bycatch estimates in the SEO SS-SPM to
# account for uncertainty... 
head(presurvey.bycatch)
#================================================================================
# Save the data for use in SS-SPMs
colnames(SEO.subs)
colnames(presurvey.bycatch)
SEO.red<-SEO.subs[,c(1,3,19,7,11,20,24,25,13,17,26,30,31)]
SEO.red$expBy32.cv2<-SEO.red$expBy32.cv
SEO.red$expByall.cv2<-SEO.red$expByall.cv
colnames(SEO.red)

SEO.expBy<-rbind(presurvey.bycatch,SEO.red)
#SEO.expBy$expBy32_mt.mean<-SEO.expBy$expBy32_mt.mean*0.00045359
colnames(SEO.expBy)

ggplot(data = SEO.expBy, aes(x = Year)) +
  geom_line(aes(y = expBy32_mt.mean, col = mngmt.area), size = 1) +
  geom_ribbon(aes(ymin = expBy32_mt.mean - expBy32.cv2*expBy32_mt.mean, 
                  ymax = expBy32_mt.mean + expBy32.cv2*expBy32_mt.mean, 
                  fill = mngmt.area),  alpha = 0.2) +
  #geom_point(aes(y = proportion, col = Source)) +  
  #expand_limits(y = c(0.0, 1)) +
  #geom_hline(yintercept = 0.5, lty = 2, col = "grey") +
  xlab("\nYear") +
  ylab("Expected yelloweye bycatch (t)\n") 
  #scale_colour_manual(values = c("black", "grey")) +
  #scale_fill_manual(values = c("black", "grey")) +
  #theme(legend.position = "none") -> byyear_plot

#plot_grid(byage_plot, byyear_plot, align = c("h"), ncol = 1)
YEAR<-2022
ggsave(paste0("Figures/exp_bycatch_in_halibut_fishery_", YEAR, ".png"), dpi=300,  height=6, width=7, units="in")

colnames(SEO.expBy)
ggplot(data = SEO.expBy %>% filter(Year > 1997), aes(x = Year)) +
  geom_line(aes(y = WCPUE32.mean, col = mngmt.area), size = 1) +
  geom_ribbon(aes(ymin = WCPUE32.mean - WCPUE32.mean*(sqrt(WCPUE32.var)/WCPUE32.mean), 
                  ymax = WCPUE32.mean + WCPUE32.mean*(sqrt(WCPUE32.var)/WCPUE32.mean), 
                  fill = mngmt.area),  alpha = 0.2) +
  xlab("\nYear") +
  ylab("wcpue/ yelloweye bycatch rate (kg yelloweye/kg legal halibut)") 

YEAR<-2022
ggsave(paste0("Figures/wcpue_ests_", YEAR, ".png"), dpi=300,  height=6, width=7, units="in")

SEO.expBy[SEO.expBy$mngmt.area == "NSEO" & SEO.expBy$Year == 2021,]
write.csv(SEO.expBy,"Data/SEO_expBy_7.31.22.csv")

#=====================================================================
# compare to landed bycatch
Removals.subd.1<-read.csv("Data/SEsubdistrict_YE_removals.csv")  #YE reconstruction from Rhea
Removals.subd<-read.csv("Data/SE_YE_known_removals_1980-2022.csv")  #YE reconstruction from most recent
#2) pre-1980 removals from SEO as a whole... 
# Removals.pre<-read.csv("Data/XXX.csv") #Waiting on Donnie; 
# generate fake data here for model development...
Removals.subd %>% filter (Mgt.Area != "SEO" & Mgt.Area != "SSEI" & Mgt.Area != "NSEI") -> Removals.subd
#  group_by(Year) %>%
#  summarise(tot.rem = sum(tot.rem.mt),
#            tot.hal.by = sum(tot.hal.by))-> Removals.subd


str(Removals.subd.1); str(Removals.subd)

unique(Removals.subd$Year)

Removals.subd$mngmt.area<-Removals.subd$Subdistrict
Removals.subd$mngmt.area<-Removals.subd$Mgt.Area

Tst<-SEO.expBy %>% left_join(Removals.subd, by=c("Year","mngmt.area"))
str(Tst)

labels<-c(eb="Expected Bycatch",lb="Landed Bycatch")
breaks<-c("eb","lb")

ggplot(data = Tst, aes(x = Year)) +
  geom_line(aes(y = expBy32_mt.mean, color="expected bycatch"), size = 1) +
  geom_ribbon(aes(ymin = expBy32_mt.mean - expBy32.cv2*expBy32_mt.mean, 
                  ymax = expBy32_mt.mean + expBy32.cv2*expBy32_mt.mean),  alpha = 0.5) +
#  geom_point(aes(y = YE_fr_hal_fish.mt, color="landed bycatch")) +  
  geom_point(aes(y = tot.hal.by, color="landed bycatch")) +
  geom_vline(xintercept=2005, col="blue") +
  geom_vline(xintercept=2000, col="blue", linetype=3) +
  facet_wrap(~ mngmt.area) +
  #expand_limits(y = c(0.0, 1)) +
  #geom_hline(yintercept = 0.5, lty = 2, col = "grey") +
  xlab("\nYear") +
  ylab("Expected and landed yelloweye bycatch (t)\n") +
  annotate("text", x=2012, y=200, col="blue", label="full retention", size=2.5) +
  #theme_minimal()+
  scale_colour_manual(labels=labels,
                      values=c("black","darkcyan")) +
                      #values = c(eb="black", lb="darkcyan"), labels=labels, breaks=breaks) +
  scale_fill_manual(values = c(eb="black", lb="darkcyan"), labels=labels, breaks=breaks) +
  theme(legend.position = "top")
  #labs(x= "\nLanded Bycatch", y="\nExpected Bycatch")

ggsave(paste0("Figures/exp_and_landed_bycatch_in_halibut_fishery_", YEAR, "_2.png"), dpi=300,  height=6, width=7, units="in")

#*******************************************************************************
#*GET SEO scale data
#*******************************************************************************
SEOHal1<-read.csv("Data/SEO_Halibut_removals_1888-2021.csv", header=T)
SEOHal2<-read.csv("Data/SE_Halibut_removals_1975-2022.csv", header=T) %>%
  filter(Mgt.Area == "CSEO" | Mgt.Area == "SSEO" | Mgt.Area == "NSEO" | Mgt.Area == "EYKT") 

ys<-unique(SEOHal2$Year)
SEOHal3<-data.frame()
j<-1
for (i in 1:length(ys)){  #i<-1
  dat<-SEOHal2[SEOHal2$Year == ys[i],]
  SEOHal3[j,"Year"]<-ys[i]
  SEOHal3[j,"SEO_harvest_mt"]<-sum(dat$HA.mt)
  SEOHal3[j,"var"]<-0
  j<-j+1
}

SEOHal<-SEOHal3

pres.yrs<-seq(1888,1997,by=1)
pres.yrs<-seq(1975,1997,by=1)
subs<-unique(SE.out$mngmt.area)

?cov

presurvey.bycatch<-data.frame()
i<-1
for (y in pres.yrs) { #y<-pres.yrs[1]
  for (s in subs) {  #s<-subs[1]
    ma<-SE.out[SE.out$mngmt.area ==s,]
    cov32<-cov(ma$WCPUE32.mean,ma$Halibut_mt)
    covall<-cov(ma$WCPUEall.mean,ma$Halibut_mt)
    hal<-SEOHal[SEOHal$Year == y,]
    H.catch<-sum(hal$SEO_harvest_mt)
    H.var<-hal$var
    mean.wcpue<-mean(ma$WCPUE32.mean)
    var.wcpue<-sum(ma$WCPUE32.var)/(length(unique(ma$Year))^2)
    
    mean.wcpueall<-mean(ma$WCPUEall.mean)
    var.wcpueall<-sum(ma$WCPUEall.var)/(length(unique(ma$Year))^2)
    
    presurvey.bycatch[i,"Year"]<-y
    presurvey.bycatch[i,"mngmt.area"]<-s
    presurvey.bycatch[i,"Halibut_mt"]<-H.catch
    presurvey.bycatch[i,"var.Halibut_mt"]<-H.var
    presurvey.bycatch[i,"WCPUE32.mean"]<-mean.wcpue
    presurvey.bycatch[i,"WCPUE32.var"]<-var.wcpue
    presurvey.bycatch[i,"expBy32_mt.mean"]<-mean.wcpue*H.catch
    presurvey.bycatch[i,"expBy32.var"]<-2*mean.wcpue*H.catch*cov32+var.wcpue*H.catch*H.catch+H.var*mean.wcpue*mean.wcpue
      #(H.catch^2)*var.wcpue
    presurvey.bycatch[i,"expBy32.cv"]<-sqrt(presurvey.bycatch[i,"expBy32.var"])/
      presurvey.bycatch[i,"expBy32_mt.mean"]
    #want to use inflated cv for historical so also save max cv from survey years... 
    presurvey.bycatch[i,"expBy32.cv2"]<-max(ma$expBy32.cv)
    
    presurvey.bycatch[i,"WCPUEall.mean"]<-mean.wcpueall
    presurvey.bycatch[i,"WCPUEall.var"]<-var.wcpueall
    presurvey.bycatch[i,"expByall_mt.mean"]<-mean.wcpueall*H.catch
    presurvey.bycatch[i,"expByall.var"]<-2*mean.wcpueall*H.catch*covall+var.wcpueall*H.catch*H.catch+H.var*mean.wcpueall*mean.wcpueall
      #(H.catch^2)*var.wcpueall
    presurvey.bycatch[i,"expByall.cv"]<-sqrt(presurvey.bycatch[i,"expByall.var"])/
      presurvey.bycatch[i,"expByall_mt.mean"]
    #want to use inflated cv for historical so also save max cv from survey years... 
    presurvey.bycatch[i,"expByall.cv2"]<-max(ma$expByall.cv)
    i<-i+1
  }
}

#-------------------------------------------------------------------------------
colnames(SE.out)
colnames(presurvey.bycatch)
SEO.red<-SE.out[,c(1,3,19,7,11,20,24,25,13,17,26,30,31)]
SEO.red$expBy32.cv2<-SEO.red$expBy32.cv
SEO.red$expByall.cv2<-SEO.red$expByall.cv
SEO.red$var.Halibut_mt<-0
colnames(SEO.red)

head(SE.out)

SE.expBy<-rbind(presurvey.bycatch,SEO.red)
#SEO.expBy$expBy32_mt.mean<-SEO.expBy$expBy32_mt.mean*0.00045359
colnames(SE.expBy)
SE.expBy$mngmt.area

ggplot(data = SE.expBy, aes(x = Year)) +
  geom_line(aes(y = expBy32_mt.mean, col = mngmt.area), size = 1) +
  geom_ribbon(aes(ymin = expBy32_mt.mean - expBy32.cv*expBy32_mt.mean, 
                  ymax = expBy32_mt.mean + expBy32.cv*expBy32_mt.mean, 
                  fill = mngmt.area),  alpha = 0.2) +
  #geom_point(aes(y = proportion, col = Source)) +  
  #expand_limits(y = c(0.0, 1)) +
  #geom_hline(yintercept = 0.5, lty = 2, col = "grey") +
  xlab("\nYear") +
  ylab("Expected yelloweye bycatch (t)\n") +
#scale_colour_manual(values = c("black", "grey")) +
#scale_fill_manual(values = c("black", "grey")) +
theme(legend.position = "none") #-> byyear_plot

#plot_grid(byage_plot, byyear_plot, align = c("h"), ncol = 1)
YEAR<-2022
ggsave(paste0("Figures/exp_bycatch_in_halibut_fishery_", YEAR, ".png"), dpi=300,  height=6, width=7, units="in")

write.csv(SE.expBy,"Data/SE_expBy_1888_2020_7.31.22.csv")

#=====================================================================
# compare to landed bycatch
Removals.subd<-read.csv("Data/SE_YE_known_removals_1980-2022.csv")  #YE reconstruction from Rhea
str(Removals.subd); unique(Removals.subd$Year)
Removals.subd$mngmt.area<-Removals.subd$Subdistrict

SEOrem<-Removals.subd %>% filter(Mgt.Area != "NSEI" &
                                   Mgt.Area != "SSEI" &
                                   Mgt.Area != "SEO") %>%
  group_by(Year) %>% 
  summarize(tot.rem.mt = sum(tot.rem.mt),
            tot.lnd.by = sum(tot.hal.by))

Tst<-SE.expBy %>% left_join(SEOrem, by=c("Year"))

labels<-c(eb="Expected Bycatch",lb="Landed Bycatch")
breaks<-c("eb","lb")

ggplot(data = Tst, aes(x = Year)) +
  geom_line(aes(y = expBy32_mt.mean, color="expected bycatch"), size = 1) +
  geom_ribbon(aes(ymin = expBy32_mt.mean - expBy32.cv2*expBy32_mt.mean, 
                  ymax = expBy32_mt.mean + expBy32.cv2*expBy32_mt.mean),  alpha = 0.5) +
  geom_line(aes(y = tot.lnd.by, col="landed bycatch"),lwd=1.1, 
            col="darkorange", lty=1) +  
  geom_vline(xintercept=2005, col="blue") +
  geom_vline(xintercept=2000, col="blue", linetype=3) +
  #expand_limits(y = c(0.0, 1)) +
  #geom_hline(yintercept = 0.5, lty = 2, col = "grey") +
  xlab("\nYear") +
  ylab("Expected and landed yelloweye bycatch (t)\n") +
  annotate("text", x=2012, y=300, col="blue", label="full retention", size=2.7) +
  #theme_minimal()+
  scale_colour_manual(labels=c("Expected Bycatch","Landed Bycatch"),
                      values=c("black","orange")) +
  #values = c(eb="black", lb="darkcyan"), labels=labels, breaks=breaks) +
  scale_fill_manual(values = c(eb="black", lb="darkorange"), 
                    labels=c("Expected Bycatch","Landed Bycatch"), 
                    breaks=c("eb","lb")) +
  theme(legend.position = "top")
#labs(x= "\nLanded Bycatch", y="\nExpected Bycatch")

ggsave(paste0("Figures/SEO_exp_and_landed_bycatch_in_halibut_fishery_", YEAR, "_2.png"), dpi=300,  height=4, width=8, units="in")

ggplot(data = Tst, aes(x = Year)) +
  geom_line(aes(y = expBy32_mt.mean, color="expected bycatch"), size = 1) +
  geom_ribbon(aes(ymin = expBy32_mt.mean - expBy32.cv2*expBy32_mt.mean, 
                  ymax = expBy32_mt.mean + expBy32.cv2*expBy32_mt.mean),  alpha = 0.5) +
  geom_point(aes(y = tot.lnd.by, color="landed bycatch")) +  
  geom_vline(xintercept=2005, col="blue") +
  geom_vline(xintercept=2000, col="blue", linetype=3) +
  facet_wrap(~ mngmt.area) +
  #expand_limits(y = c(0.0, 1)) +
  #geom_hline(yintercept = 0.5, lty = 2, col = "grey") +
  xlab("\nYear") +
  ylab("Expected and landed yelloweye bycatch (t)\n") +
  annotate("text", x=2012, y=200, col="blue", label="full retention", size=2.5) +
  #theme_minimal()+
  scale_colour_manual(labels=labels,
                      values=c("black","darkcyan")) +
  #values = c(eb="black", lb="darkcyan"), labels=labels, breaks=breaks) +
  scale_fill_manual(values = c(eb="black", lb="darkcyan"), labels=labels, breaks=breaks) +
  theme(legend.position = "top")
#labs(x= "\nLanded Bycatch", y="\nExpected Bycatch")

#---------------------------------------------------------------------------------
## Get SEO discard estimate for SAFE report table

str(Tst)

SAFE.discards<-Tst[Tst$Year > 1979,]

SAFE.discards$Discards.mt <-SAFE.discards$expBy32_mt.mean - SAFE.discards$tot.lnd.by
SAFE.discards$min.D<-1
colnames(SAFE.discards)

SAFE.discards$Discards.mt<-rowMaxs(as.matrix(SAFE.discards[,c(19,20)]))

SAFE.disc<-SAFE.discards[,c(1,19)]

write.csv(SAFE.disc,"Data/SEO_YEdisc_in_hali_fsh_TribuzioMethods_1980_2021_10.3.22.csv")






