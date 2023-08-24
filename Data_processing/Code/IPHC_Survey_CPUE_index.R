################################################################################
## This code is used to process, cull and organize the IPHC survey data for use
## as a secondary index of yelloweye rockfish abundance in the SEO
YEAR<-2023

{library(dplyr)
library(boot)
library(ggplot2)
library(RColorBrewer)
library(sf)
library(readr)
library(vroom)
library(scales)
}
#wd="C:/Users/pjjoy/Documents/Groundfish Biometrics/Yelloweye_Production_Models/"
#setwd(wd)
#getwd()

source("r_helper/Port_bio_function.R")

#*******************************************************************************
### NEED TO DO!!!!!!!! 7-12-23
# 1) Add in hook saturation correction factor for both cpue and bycatch 
# 2) filter out stations that were not effective!!! "Eff" column... Yes = good, no = toss it
# 3) filter stations that fit in the "designated yelloweye habitat" that we use for the ROV surveys
#    Station locations move by 3nm year to year, so may be in and out of actual habitat
#*******************************************************************************

#IPHCfunction<-function(){
{
  #HA.Harv<-read.csv("Data_processing/Data/halibut_catch_data.csv", header=T)
  
  #GET Latest Data and prep to add into 1998-2020 data that script is set to handle
  But2C3A_2022<-read.csv("Data_processing/Data/IPHC_raw/Set and Pacific halibut data 2C3A 2022.csv", check.names = TRUE)
#  View(But2C3A_2022)
  
  YE2C3A_2022<-read.csv("Data_processing/Data/IPHC_raw/Non-Pacific halibut data 2C3A 2022.csv", check.names = TRUE)
  unique(YE2C3A_2022$Year)
#  View(YE2C3A_2022)
  
  Surv22<-But2C3A_2022 %>% full_join(YE2C3A_2022, by=c("Stlkey","Station")) #by=c("Stlkey","Setno","Station"))
  colnames(Surv22)
  Sur2C_22<-Surv22[Surv22$IPHC.Reg.Area == "2C",]
  Sur3A_22<-Surv22[Surv22$IPHC.Reg.Area == "3A",]
  
  But2C3A_2021<-read.csv("Data_processing/Data/IPHC_raw/IPHC Set and Pacific halibut data 2C3A 2021.csv", header=T)
  YE2C3A_2021<-read.csv("Data_processing/Data/IPHC_raw/Non-Pacific halibut data YE 2C3A 2021.csv", header=T)
  
  Surv21<-But2C3A_2021 %>% full_join(YE2C3A_2021, by=c("Stlkey","Station")) #by=c("Stlkey","Setno","Station"))
  head(Surv21)
  Sur2C_21<-Surv21[Surv21$IPHC.Reg.Area == "2C",]
  Sur3A_21<-Surv21[Surv21$IPHC.Reg.Area == "3A",]
  #******************************************************************************
  #*
  BUT2C<-read.csv("Data_processing/Data/IPHC_raw/IPHC Set and Pacific halibut data 2C.csv", header=T)
  YE2C<-read.csv("Data_processing/Data/IPHC_raw/Non-Pacific halibut data YE 2C.csv", header=T)
  
  Sur2C<- BUT2C %>% full_join(YE2C, by=c("Stlkey","Station")) #by=c("Stlkey","Setno","Station"))  
  Sur2C <- Sur2C %>% rbind(Sur2C, Sur2C_21, Sur2C_22)
  
  Sur2C[Sur2C$Stlkey == 20200208,]
  
  Sur2C<- Sur2C %>%
    mutate(SEdist = ifelse(IPHC.Stat.Area %in% c(182:184,171,173,174,161:163),"NSEI", 
                           ifelse(IPHC.Stat.Area %in% c(142:144,152,153),"SSEI",
                                  ifelse(IPHC.Stat.Area %in% c(140,141,151,150,160,170,181),
                                         ifelse(-MidLon.fished>137,"EYKT",
                                                ifelse(MidLat.fished>=57.5,"NSEO",
                                                       ifelse(MidLat.fished<57.5 & MidLat.fished>=56,"CSEO",
                                                              ifelse(MidLat.fished<56,"SSEO",NA)))),NA))))
  
  nrow(Sur2C)
  ## Lets bring in 3A and pull out surveys in the Yakutat area...
  BUT3A<-read.csv("Data_processing/Data/IPHC_raw/IPHC Set and Pacific halibut data 3A.csv", header=T)
  YE3A<-read.csv("Data_processing/Data/IPHC_raw/Non-Pacific halibut data YE 3A.csv", header=T)
  
  Sur3A<- BUT3A %>% full_join(YE3A, by=c("Stlkey","Station")) #by=c("Stlkey","Setno","Station"))
  Sur3A <- Sur3A %>% rbind(Sur3A,Sur3A_21, Sur3A_22)
  
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
  
  Depths<-unique(Survey$depth_bin)
  
  Survey$O32.Pacific.halibut.weight<-as.numeric(Survey$O32.Pacific.halibut.weight)
  Survey$U32.Pacific.halibut.weight<-as.numeric(Survey$U32.Pacific.halibut.weight)
  
  Survey<-unique(Survey %>% select(-c(Row.number.x,Row.number.y)))
  #str(Survey)
  #Bring in and join hook adjustment factor to control for hook saturation when 
  # calculating CPUE
  
  hadj<-read.csv(paste0("Data_processing/Data/IPHC_raw/iphc_",YEAR,"_fiss_hadj.csv"), check.names = TRUE) %>% 
  #hadj<-hadj %>% 
    filter(IPHC.Reg.Area == "3A " | IPHC.Reg.Area == "2C ") %>%
    select(Stlkey = stlkey,Year,Station,h.adj)
    str(hadj)
    
    
  Survey <- left_join(Survey,hadj,by=c("Year","Stlkey","Station"))  
  #str(Survey)
}

{y<-sample(unique(Survey$Year),1)
d<-sample(unique(Survey$depth_bin[Survey$Year==y]),1)
s<-sample(unique(Survey$Station[Survey$Year==y & Survey$depth_bin==d]),1)
Dat<-Survey[Survey$Year == y & Survey$depth_bin == d & Survey$Station == s,]
unique(Dat); nrow(Dat)}

str(Survey); head(Survey,10)
str(HA.Harv)
mean(Survey$AvgDepth.fm)

#Survey <-Survey %>% filter(Year < 2022)
write.csv(Survey,paste0("Data_processing/Data/IPHC_survey_1998-",YEAR,".csv"))

#================================================================================
# Look at some depth stuff and identify stations that have seen yelloweye...
hist(Survey$BeginDepth..fm.)
hist(Survey$EndDepth..fm.)
hist(Survey$AvgDepth.fm)

plot(Survey$YE.obs ~ Survey$AvgDepth.fm)
plot(Survey$O32.Pacific.halibut.count ~ Survey$AvgDepth.fm)
plot(Survey$U32.Pacific.halibut.count ~ Survey$AvgDepth.fm)

Survey %>% group_by(Station) %>%
  dplyr::summarise(Lat = mean(MidLat.fished, na.rm=T),
                   Lon = mean(MidLon.fished, na.rm=T),
                   Depth = mean(AvgDepth.fm, na.rm=T),
                   obs = n(),
                   mean.hal.O32.count = mean(O32.Pacific.halibut.count, na.rm=T),
                   mean.hal.O32.wt = mean(O32.Pacific.halibut.weight, na.rm=T),
                   mean.hal.U32.count = mean(U32.Pacific.halibut.count, na.rm=T),
                   mean.hal.U32.wt = mean(U32.Pacific.halibut.weight, na.rm=T),
                   mean.skates = mean(No..skates.set, na.rm=T),
                   mean.hooks.ob = mean(HooksObserved, na.rm=T),
                   mean.hooks.ret = mean(HooksRetrieved, na.rm=T),
                   mean.YE = mean(YE.obs, na.rm=T),
                   var.YE = var(YE.obs, na.rm=T),
                   mean.YEcpue = mean(YE.obs/HooksObserved, na.rm=T),
                   var.YEcpue = var(YE.obs/HooksObserved, na.rm=T),
                   noYE.count = sum(YE.obs == 0),
                   prop.0 = noYE.count/obs
                   ) -> Station.sum
nrow(Station.sum)

hist(Station.sum$mean.YE, breaks=25); nrow(Station.sum[Station.sum$mean.YE == 0,])/nrow(Station.sum)
hist(Station.sum$mean.YEcpue, breaks=25); nrow(Station.sum[Station.sum$mean.YEcpue == 0,])/nrow(Station.sum)
hist(Station.sum$noYE.count, breaks=25)
hist(Station.sum$prop.0, breaks=25)
hist(Station.sum$prop.0[Station.sum$prop.0<1], breaks=25)
max(Station.sum$prop.0[Station.sum$prop.0<1])

blanks.stations<-Station.sum$Station[Station.sum$mean.YEcpue == 0]
YE.stations<-Station.sum$Station[Station.sum$mean.YEcpue != 0]
YE.stations_10p<-Station.sum$Station[Station.sum$prop.0 < 0.9]
YE.stations_20p<-Station.sum$Station[Station.sum$prop.0 < 0.8]
YE.stations_25p<-Station.sum$Station[Station.sum$prop.0 < 0.75]
YE.stations_40p<-Station.sum$Station[Station.sum$prop.0 < 0.6]

length(YE.stations);length(YE.stations_10p);length(YE.stations_20p);length(YE.stations_25p);length(YE.stations_40p)

plot_stations<-st_as_sf(Station.sum,coords = c("Lon","Lat"))

plot_stations<-st_set_crs(plot_stations,crs=4326)

use_stations<-YE.stations

ggplot(plot_stations %>% filter(Station %in% c(use_stations))) + 
  geom_sf(aes(color = mean.YE, size = mean.YE)) +
  scale_color_viridis_c(option = "magma",begin = 0)

ggplot(plot_stations %>% filter(Station %in% c(use_stations))) + 
  geom_sf(aes(color = mean.YEcpue, size = mean.YEcpue)) +
  scale_color_viridis_c(option = "plasma",begin = 0)

ggplot(plot_stations %>% filter(Station %in% c(use_stations))) + 
  geom_sf(aes(color = var.YEcpue, size = var.YEcpue)) +
  scale_color_viridis_c(option = "viridis",begin = 0)

ggplot(plot_stations %>% filter(Station.x %in% c(use_stations))) + 
  geom_sf(aes(color = prop.0, size = prop.0)) +
  scale_color_viridis_c(option = "plasma",begin = 0)
#==============================================================================================
## Load the port samples that were downloaded from oceansAK; 
# Need Yelloweye weights to get wcpue estimates
{
  #Port1<-read.csv("Data/port_sampling_bio_data.csv")
  #Port2<-read.csv("Data/port_sampling_bio_data2.csv")
  #Port3<-read.csv("Data/port sampling bio data 2021-2022 071522.csv")
  
  #Port<-rbind(Port1,Port2,Port3)
  
  Port<-port.bio(2022)
  str(Port)
  Port$Year<-as.integer(Port[,1])
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
# FUNCTION for generating cpue estimates from IPHC surveys for different management areas

YEHA.fxn<-function(Survey=Survey, Area="SEdist",Deep=250, Shallow=0,  nboot=1000){
  col.ref<-which(colnames(Survey)==Area)
  
  IPHC.cpue<-data.frame()
  Subs<-unique(Survey[,col.ref])
  Years<-unique(Survey$Year)
  
  j<-1
  for (y in Years) {  #y<-Years[1]
    for (s in Subs){  #s<-Subs[1]
      Dat<-Survey[Survey$Year == y & 
                    Survey[,col.ref] == s &
                    Survey$AvgDepth.fm > Shallow & 
                    Survey$AvgDepth.fm < Deep &
                    Survey$Eff == "Y",]
      if (nrow(Dat)>0){
        Stations<-unique(Dat$Station)
        #WCPUEi.32<-vector()
        #WCPUEi.all<-vector()
        CPUEi<-vector()
        i<-1
        for (st in Stations){    #st<-Stations[1]   length(Stations)  st<-1
          Stat.Dat<-Dat[Dat$Station == st,] #; Stat.Dat
          #debug
          #if (nrow(Stat.Dat) > 1){aaa} else {}
          CPUE<-mean(Stat.Dat$YE.obs/Stat.Dat$HooksObserved)
          #hook adjustment factor
          CPUE<-CPUE*unique(Stat.Dat$h.adj)
          
          if (CPUE == 0){
            C<-0
          } else {
            C<-CPUE*Stat.Dat$HooksRetrieved
          }
          
          CPUEi[i]<-CPUE
          i<-i+1
        }
        
        Out<-data.frame()
        #CPUE.out<-data.frame()
        for (ii in 1:nboot){ #nboot<-1000  ii<-1
          Resamp3<-sample(CPUEi,length(CPUEi),replace=T)
          Out[ii,"CPUE"]<-mean(Resamp3, na.rm=T)
          Out[ii,"CPUE.var"]<-var(Resamp3, na.rm=T)
        }
        
        #hist(Out$KgHa, breaks = 25)
        IPHC.cpue[j,"Year"]<-y
        IPHC.cpue[j,"mngmt.divisions"]<-Area
        IPHC.cpue[j,"mngmt.area"]<-s
        IPHC.cpue[j,"deep.bound"]<-Deep
        IPHC.cpue[j,"shallow.bound"]<-Shallow
        IPHC.cpue[j,"no.stations"]<-length(Stations)
        
        IPHC.cpue[j,"CPUE.mean"]<-mean(CPUEi)
        IPHC.cpue[j,"CPUE.bootmean"]<-unname(quantile(Out$CPUE,c(0.5)))  #mean WCPUE from Tribuzio
        IPHC.cpue[j,"CPUE.lo95ci"]<-unname(quantile(Out$CPUE,c(0.025)))
        IPHC.cpue[j,"CPUE.hi95ci"]<-unname(quantile(Out$CPUE,c(0.975)))
        IPHC.cpue[j,"CPUE.var"]<-var(Out$CPUE)
        IPHC.cpue[j,"CPUE.cv"]<-sd(Out$CPUE)/mean(Out$CPUE)
        #    IPHC.wcpue[j,"WCPUE32.cv"]<-sqrt(var(WCPUEi.32))/mean(WCPUEi.32)
       
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
  for (n in 1:nplot){  #n<-1
    i<-1
    pSubs<-Subs[(1+(n-1)*8):(8+(n-1)*8)]
    pSubs<-pSubs[!is.na(pSubs)]
    
    for (s in pSubs) {   #s<-pSubs[1]
      Pdat<-IPHC.cpue[IPHC.cpue$mngmt.area == s,]
      if (s == pSubs[1]){
        plot(Pdat$CPUE.mean ~ Pdat$Year, col=cols[i], 
             ylim=c(0,max(IPHC.cpue$CPUE.hi95ci)), type="l",
             ylab="CPUE", xlab="Year", main=Area)
        Y1<-Pdat$CPUE.lo95ci; Y2<-Pdat$CPUE.hi95ci
        polygon(c(Years,rev(Years)),c(Y1,rev(Y2)),
                col=adjustcolor(cols[i],alpha.f=0.2),border=NA)
        i<-i+1
      } else {
        Y1<-Pdat$CPUE.lo95ci; Y2<-Pdat$CPUE.hi95ci
        if (length(Y1)==1){}else {
          polygon(c(Years,rev(Years)),c(Y1,rev(Y2)),
                  col=adjustcolor(cols[i],alpha.f=0.2),border=NA)
          i<-i+1
        }
        
      }
    }
    i<-1
    for (s in pSubs) {   #s<-Subs[2]
      Pdat<-IPHC.cpue[IPHC.cpue$mngmt.area == s,]
      lines(Pdat$CPUE.mean ~ Pdat$Year, col=cols[i], type="l",
            lwd=1.5)
      i<-i+1
    }
    legend(x="topleft", cex=0.8,
           legend=c(pSubs), bty="n",
           col=cols,text.col=cols)
  }
  
  return(IPHC.cpue)
}


colnames(Survey)

#Decision: which stations to use to calculate CPUE??? 
Survey_40p<-Survey %>% filter(Station %in% c(YE.stations_40p))
Survey_non0<-Survey %>% filter(Station %in% c(YE.stations))

#IPHC.reg.area<- YEHA.fxn(Area="IPHC.Reg.Area",Deep=max(Survey$AvgDepth.fm), Shallow=0, nboot=1000) 
#IPHC.stat.area<- YEHA.fxn(Area="IPHC.Stat.Area",Deep=max(Survey$AvgDepth.fm), Shallow=0, nboot=1000) 
SE.subdistricts<- YEHA.fxn(Survey=Survey_non0, Area="SEdist",Deep=250, Shallow=0, nboot=1000)
SE.subdistricts_40p<- YEHA.fxn(Survey=Survey_40p, Area="SEdist",Deep=250, Shallow=0, nboot=1000)

str(SE.subdistricts); unique(SE.subdistricts$mngmt.area)

IPHC.index<-SE.subdistricts[SE.subdistricts$mngmt.area != "NSEI" & 
                              SE.subdistricts$mngmt.area != "SSEI",]
str(IPHC.index)

ggplot(IPHC.index, aes(x=Year)) +
  geom_point(aes(y=CPUE.mean),size=2) +
  geom_errorbar(aes(ymin = CPUE.lo95ci, ymax= CPUE.hi95ci),
                col="black", alpha=0.5) +
  facet_wrap(~mngmt.area) +
  xlab("\nYear") +
  ylab("Yelloweye CPUE n/hook") +
  #ylab("Density (yelloweye rockfish/kmsq)") +
  scale_y_continuous(label=comma) +
  scale_x_continuous(breaks=seq(1995,2025,5)) + 
  theme (axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
         panel.grid.minor = element_blank()) +
  ggtitle("Yelloweye cpue in IPHC FISS, stations that encountered yelloweye at least once")

ggsave(paste0("Figures/YE_IPHC_cpue_non0_", YEAR, ".png"), dpi=300,  height=5, width=5, units="in")

write.csv(IPHC.index,paste0("Data_processing/Data/IPHC.cpue.SEO_non0_",YEAR,".csv"))

IPHC.index_40p<-SE.subdistricts_40p[SE.subdistricts_40p$mngmt.area != "NSEI" & 
                              SE.subdistricts_40p$mngmt.area != "SSEI",]

ggplot(IPHC.index_40p, aes(x=Year)) +
  geom_point(aes(y=CPUE.mean),size=2) +
  geom_errorbar(aes(ymin = CPUE.lo95ci, ymax= CPUE.hi95ci),
                col="black", alpha=0.5) +
  facet_wrap(~mngmt.area) +
  xlab("\nYear") +
  ylab("Yelloweye CPUE n/hook") +
  #ylab("Density (yelloweye rockfish/kmsq)") +
  scale_y_continuous(label=comma) +
  scale_x_continuous(breaks=seq(1995,2025,5)) + 
  theme (axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
         panel.grid.minor = element_blank()) +
  ggtitle("Yelloweye cpue in IPHC FISS, stations that encountered yelloweye at least 40% of years")

ggsave(paste0("Figures/YE_IPHC_cpue_min40percentYE_", YEAR, ".png"), dpi=300,  height=5, width=5, units="in")

write.csv(IPHC.index_40p,paste0("Data_processing/Data/IPHC.cpue.SEO_min40percentYE_",YEAR,".csv"))

comp<-rbind(IPHC.index %>% mutate(yelloweye_presence = "at least once"),
            IPHC.index_40p %>% mutate(yelloweye_presence = "40% of years"))

ggplot(comp, aes(x=Year)) +
  geom_point(aes(y=CPUE.mean, col = yelloweye_presence),
             size=2, position = position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin = CPUE.lo95ci, ymax= CPUE.hi95ci, col=yelloweye_presence),
                alpha=0.5,
                position = position_dodge(width=0.25)) +
  facet_wrap(~mngmt.area) +
  xlab("\nYear") +
  ylab("Yelloweye CPUE n/hook") +
  #ylab("Density (yelloweye rockfish/kmsq)") +
  scale_y_continuous(label=comma) +
  scale_x_continuous(breaks=seq(1995,2025,5)) + 
  theme (axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
         panel.grid.minor = element_blank(),
         legend.position = "top") +
  ggtitle("Yelloweye cpue in IPHC FISS")

ggsave(paste0("Figures/YE_IPHC_cpue_method_comparison_", YEAR, ".png"), dpi=300,  height=5, width=5, units="in")






















