library(dplyr)
library(boot)
library(ggplot2)
library(RColorBrewer)
#wd="C:/Users/pjjoy/Documents/Groundfish Biometrics/Yelloweye_Production_Models/"
#setwd(wd)
#getwd()

#IPHCfunction<-function(){
{
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
                           ifelse(SEdist %in% c("NSEI","SSEI"),"SEI",NA))) %>%
    rename(YE.obs = Number.Observed,
           AvgDepth.fm = AvgDepth..fm.,
           Year = Year.x)%>%
    mutate(depth_bin = cut(AvgDepth.fm, breaks = seq(0,400,50),
                           labels = paste (seq(50,400,50))),
           YE.obs = ifelse(is.na(YE.obs),0,YE.obs),
           HooksObserved  = ifelse(is.na(HooksObserved ),140,HooksObserved ))
  
  Survey$HooksRetrieved<-as.numeric(Survey$HooksRetrieved)
  Survey$YE.exp<-Survey$HooksRetrieved/Survey$HooksObserved
  
  Survey<-Survey[order(Survey$Year),]
}
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
  
  Depths<-unique(Survey$depth_bin)
  
  #==============================================================================================
  ## Load the port samples that were downloaded from oceansAK; 
  # Need Yelloweye weights to get wcpue estimates
  {
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
  
  
  ## Get best data for average weights and add columns to survey data
  Years<-sort(unique(Survey$Year))
  GFMA<-unique(Survey$SEdist)
  
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
      while (Nw < 75 & i+m <= max(Port$Year)) {     #go forward if you failed to find enough weights... 
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
  Survey$mean.YE.kg
  
  #=============================================================================
  # to gauge relationship between depth and wcpue estimates calculate for different depth 
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
  # Now, lets produce wcpue for all relevant management units for all depths combined
  # as well as for each depth bin so we can deal with hindcasting for the variety
  # of spatial resolutions we are likely to see
  
str(Survey)
  #need by $IPHC.Reg.Area
  #need by $IPHC.Stat.Area
  #need by $IPHC.Charter.Region
  #need by $SEdist
  #need by #In.Out
  
  #create a function to select the management area and the depth brackets 
  # to calculate WCPUE for...
  
  YEHA.fxn<-function(Area="IPHC.Stat.Area",Deep=max(Survey$AvgDepth.fm), Shallow=0, nboot=1000){
    col.ref<-which(colnames(Survey)==Area)
    
    IPHC.wcpue<-data.frame()
    Subs<-unique(Survey[,col.ref])
    Years<-unique(Survey$Year)
    
    j<-1
    for (y in Years) {  #y<-Years[1]
      for (s in Subs){  #s<-Subs[1]
        Dat<-Survey[Survey$Year == y & Survey[,col.ref] == s &
                      Survey$AvgDepth.fm > Shallow & Survey$AvgDepth.fm < Deep,]
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
          IPHC.wcpue[j,"Year"]<-y
          IPHC.wcpue[j,"mngmt.divisions"]<-Area
          IPHC.wcpue[j,"mngmt.area"]<-s
          IPHC.wcpue[j,"deep.bound"]<-Deep
          IPHC.wcpue[j,"shallow.bound"]<-Shallow
          IPHC.wcpue[j,"no.stations"]<-length(Stations)
          
          IPHC.wcpue[j,"WCPUE32.mean"]<-mean(WCPUEi.32)
          IPHC.wcpue[j,"WCPUE32.bootmean"]<-unname(quantile(Out$WCPUE.32,c(0.5)))  #mean WCPUE from Tribuzio
          IPHC.wcpue[j,"WCPUE32.lo95ci"]<-unname(quantile(Out$WCPUE.32,c(0.025)))
          IPHC.wcpue[j,"WCPUE32.hi95ci"]<-unname(quantile(Out$WCPUE.32,c(0.975)))
          IPHC.wcpue[j,"WCPUE32.var"]<-var(Out$WCPUE.32)
          IPHC.wcpue[j,"WCPUE32.cv"]<-sd(Out$WCPUE.32)/mean(Out$WCPUE.32)
          #    IPHC.wcpue[j,"WCPUE32.cv"]<-sqrt(var(WCPUEi.32))/mean(WCPUEi.32)
          
          IPHC.wcpue[j,"WCPUEall.mean"]<-mean(WCPUEi.all)
          IPHC.wcpue[j,"WCPUEall.bootmean"]<-unname(quantile(Out$WCPUE.all,c(0.5)))  #mean WCPUE from Tribuzio
          IPHC.wcpue[j,"WCPUEall.lo95ci"]<-unname(quantile(Out$WCPUE.all,c(0.025)))
          IPHC.wcpue[j,"WCPUEall.hi95ci"]<-unname(quantile(Out$WCPUE.all,c(0.975)))
          IPHC.wcpue[j,"WCPUEall.var"]<-var(Out$WCPUE.all)
          IPHC.wcpue[j,"WCPUEall.cv"]<-sd(Out$WCPUE.all)/mean(Out$WCPUE.all)
          
          j<-j+1
        } else {}
      }
    }
    
    #couple of plots in output
    Sub.cnt<-length(Subs)
    cols<-brewer.pal(min(Sub.cnt,8),"Dark2")
    
    if (Sub.cnt < 8) {par(mfrow=c(1,1))}
    if (Sub.cnt > 8 & Sub.cnt < 17) {par(mfrow=c(2,1))}
    if (Sub.cnt > 16) {par(mfrow=c(3,1))}
    
    nplot<-ceiling(Sub.cnt/8)
    #cols<-brewer.pal(length(Subs),"Set3")
    for (n in 1:nplot){  #n<-3
      i<-1
      pSubs<-Subs[(1+(n-1)*8):(8+(n-1)*8)]
      pSubs<-pSubs[!is.na(pSubs)]
      
      for (s in pSubs) {   #s<-pSubs[6]
        Pdat<-IPHC.wcpue[IPHC.wcpue$mngmt.area == s,]
        if (s == pSubs[1]){
          plot(Pdat$WCPUE32.mean ~ Pdat$Year, col=cols[i], ylim=c(0,0.25), type="l",
               ylab="WCPUE", xlab="Year", main=Area)
          Y1<-Pdat$WCPUE32.lo95ci; Y2<-Pdat$WCPUE32.hi95ci
          polygon(c(Years,rev(Years)),c(Y1,rev(Y2)),
                  col=adjustcolor(cols[i],alpha.f=0.2),border=NA)
          i<-i+1
        } else {
          Y1<-Pdat$WCPUE32.lo95ci; Y2<-Pdat$WCPUE32.hi95ci
          if (length(Y1)==1){}else {
            polygon(c(Years,rev(Years)),c(Y1,rev(Y2)),
                    col=adjustcolor(cols[i],alpha.f=0.2),border=NA)
            i<-i+1
          }
          
        }
      }
      i<-1
      for (s in pSubs) {   #s<-Subs[2]
        Pdat<-IPHC.wcpue[IPHC.wcpue$mngmt.area == s,]
        lines(Pdat$WCPUE32.mean ~ Pdat$Year, col=cols[i], type="l",
              lwd=1.5)
        i<-i+1
      }
      legend(x="topleft", cex=0.8,
             legend=c(pSubs), bty="n",
             col=cols,text.col=cols)
    }
    
     return(IPHC.wcpue)
  }

  
  colnames(Survey)
 IPHC.reg.area<- YEHA.fxn(Area="IPHC.Reg.Area",Deep=max(Survey$AvgDepth.fm), Shallow=0, nboot=1000) 
 IPHC.stat.area<- YEHA.fxn(Area="IPHC.Stat.Area",Deep=max(Survey$AvgDepth.fm), Shallow=0, nboot=1000) 
 SE.subdistricts<- YEHA.fxn(Area="SEdist",Deep=max(Survey$AvgDepth.fm), Shallow=0, nboot=1000)
 SE.in_and_out<- YEHA.fxn(Area="In.Out",Deep=max(Survey$AvgDepth.fm), Shallow=0, nboot=1000)
 
 SE.wcpue<-rbind(SE.subdistricts,SE.in_and_out)
 write.csv(SE.wcpue,"C:/Users/pjjoy/Documents/Groundfish Biometrics/SRI/Discard_Estimation/SE.wcpue.ests_4.28.22.csv")
 
 #Left off Friday 4/15/22: dialing part of plot function
 # more work on calculating cv's of wcpue in bootstrap
 # cleaner plots for depth exam to share with Rhea and Randy
 
 #if really ambitious, try pulling bayesian model out... 
 #==============================================================================
 # Bayesian Discard Model using Observer Data
 #=========================================================================================
 {
   library(MCMCvis)
   library(runjags)
   library(rjags)
   library(dplyr)
   library(boot)
   library(ggplot2)
   library(coda)
   library(MASS)  
   library(R2OpenBUGS)
   library(jagsUI)
 }
 #name the model
 Mod.title<-"SEOsubd_PTmod_FULL4_yes94_unif_hB_rchange_ABC005"
 #=========================================================================================
 cat('model { 

#Bycatch Model
for (i in 1:Subd){
for (t in 1:N){
  YEHA.ratio[i,t] ~ dbeta(pB1[i],pB2[i]) #T(0.001,0.2) #hyper prior based on observed bycatch in subdistrict 
  B1[i,t]<-((1-YEHA.ratio[i,t])/((cv.wcpue[i,t]*YEHA.ratio[i,t])^2) - 1/YEHA.ratio[i,t])*(YEHA.ratio[i,t]^2)   #s[t]<-(cv.wcpue[t]*wcpue[t])^2
  B2[i,t]<-B1[i,t]*(1/YEHA.ratio[i,t]-1)
  wcpue[i,t] ~ dbeta(B1[i,t],B2[i,t])
    
  #observed bycatch from NOAA observers
  By[i,t]<-YEHA.ratio[i,t]*Hal[i,t]
  logBy[i,t]<-log(By[i,t])
  tau.log.By[i,t] <- 1 / log(cv.By[i,t]*cv.By[i,t] + 1)
  By.obs[i,t] ~ dlnorm(logBy[i,t],tau.log.By[i,t])T(0.01,)  #(50,)
  }
}

#Discard Model; Discards = Est. Bycatch - landed bycatch
for (i in 1:Subd){
for (t in 1:N){
  D[i,t]<- max(1,By[i,t]-Lnd.By[i,t])   #~dunif(0,5000)  #was 10K
  logD[i,t]<-log(D[i,t])
  
  #halibut fleet observer program discards
  D.obs[i,t] ~ dlnorm(logD[i,t], tau.log.D[i,t])T(0.0001,10000)  #Observed discard 
  tau.log.D[i,t]<-1/log(cv.D[i,t]*cv.D[i,t]+1)
  }
}

#Known Catches
for (i in 1:Subd){
for (t in 1:N){
  KnC[i,t]~dunif(0,5000)
  logKnC[i,t]<-log(KnC[i,t])
  KnC.obs[i,t] ~ dlnorm(logKnC[i,t], tau.log.KnC[i,t])  #Observed discard estimate
  tau.log.KnC[i,t]<-1/log(cv.KnC[i,t]*cv.KnC[i,t]+1)  
  
#Total catches
  C[i,t]<-D[i,t]+KnC[i,t]
}
}

for (t in 1:N){
  Cseo[t]<-sum(C[,t])
  KnCseo[t]<-sum(KnC[,t])
  Byseo[t]<-sum(By[,t])
  Dseo[t]<-sum(D[,t])
}

#subdistrict level priors 
for (i in 1:Subd){
  pB1[i] ~ dnorm(hB1,0.5)
  pB2[i] ~ dnorm(hB2,0.5)
}

#priors/hyper priors
hB1 ~ dunif(0.01,1000) #dnorm(35.35735,0.5)   #based on mean YEHA ratio in observer data for SEO
hB2 ~ dunif(0.01,1000) #dnorm(751.1689,0.5)
#fcv.B ~ dlnorm(log(0.2+Tau2),0.1)T(0.05,0.5)
}', file=Mod.title) 
 
 #===========================================================================================================
 # prep model setting with data, initial values
 
 {
   data<-list(N=N,Subd=Subd,
              KnC.obs=KnC.obs, cv.KnC=cv.KnC,
              D.obs=D.obs, cv.D=cv.D, 
              By.obs=By.obs, cv.By=cv.By,
              wcpue=wcpue,cv.wcpue=cv.wcpue,
              Lnd.By=Lnd.By, Hal=Hal)
   
   KnCinits1<-KnC.obs+runif(N*Subd,-0.5*KnC.obs,0.5*KnC.obs)
   KnCinits2<-KnC.obs+runif(N*Subd,-0.5*KnC.obs,0.5*KnC.obs)
   KnCinits3<-KnC.obs+runif(N*Subd,-0.5*KnC.obs,0.5*KnC.obs)
   KnCinits4<-KnC.obs+runif(N*Subd,-0.5*KnC.obs,0.5*KnC.obs)
   KnCinits5<-KnC.obs+runif(N*Subd,-0.5*KnC.obs,0.5*KnC.obs)
   
   PYEinits1<-rbind(runif(N,0.01,0.09),runif(N,0.01,0.09),runif(N,0.01,0.09),runif(N,0.01,0.09))
   PYEinits2<-rbind(runif(N,0.01,0.09),runif(N,0.01,0.09),runif(N,0.01,0.09),runif(N,0.01,0.09))
   PYEinits3<-rbind(runif(N,0.01,0.09),runif(N,0.01,0.09),runif(N,0.01,0.09),runif(N,0.01,0.09))
   
   inits1<-list(YEHA.ratio=PYEinits1,
                KnC=KnCinits1)
   inits2<-list(YEHA.ratio=PYEinits2,
                KnC=KnCinits2)
   inits3<-list(YEHA.ratio=PYEinits3,
                KnC=KnCinits3)
   
   #inits <- list(inits1,inits2,inits3,inits4,inits5)
   inits<-list(inits1,inits2,inits3)
   
   parameters=c(
     "D","C", "KnC", "By",
     "Cseo","Byseo","Dseo","KnCseo",
     "B1","B2","YEHA.ratio",
     "pB1","pB2",
     "hB1","hB2" #, #, "p" #, "CtoFB"
     )               
 }
 
 #=========================================================================================================
 #### run JAGS ####
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  