################################################################################
## This code is used to process, cull and organize the IPHC survey data for use
## as a secondary index of yelloweye rockfish abundance in the SEO
#
# Methods include the original bootstrap methods used in the 2022 assessment
# and the new Tweedie estimator developed after the 2023 CIE review
# Phil Joy
# Oct. 2023


# Updated 9/24/24 by RKE and fixed some errors and issues 
# you will need to fix file location because I did not push all of my files to github
# and i set my folder structure up differently due to confusion with current organization 
################################################################################
YEAR<-2023

{library(dplyr)
library(boot)
library(ggplot2)
library(RColorBrewer)
library(sf)
library(readr)
library(vroom)
library(scales)
library(ggpubr)
library(mgcv)
library(MuMIn)
library(GGally)
library(mgcViz)
library(wesanderson)
library(RColorBrewer)
library(tweedie)
}

source(file = paste0(here::here(), "/r_helper/Port_bio_function.R"))

#*******************************************************************************
#*******************************************************************************
# RKE truncated this code and cleaned up so you only need to import the full data files from IPHC once and
# avoid join issues 

HA.Harv<-read_csv("./Data_processing/Data/IPHC_raw/halibut_catch_data_rke.csv")
  unique(HA.Harv$year.landed)#1975-2024
 
BUT2C3A <- read_csv("./Data_processing/Data/IPHC_raw/Set and Pacific halibut data_rke.csv") %>% 
    dplyr::rename_all(~make.names(.))

YE2C3A <- read_csv("./Data_processing/Data/IPHC_raw/Non-Pacific halibut data_rke.csv") %>% 
    dplyr::rename_all(~make.names(.))
  
# RKE corrected EYKT bounds to 137-140 Longitude as that is the management area
  # Phil had 137-139 so it was exluding 1/3 of the area. Not sure if that was intentional or an error 
  # or if this change impacts the overall results at all 
Survey_prep <- BUT2C3A %>% 
    full_join(YE2C3A, by=c("Stlkey","Station")) %>% 
    filter(IPHC.Stat.Area <= 200, 
           -MidLon.fished <= 140, 
           Eff == "Y") %>%  
    mutate(SEdist = ifelse(IPHC.Stat.Area %in% c(182:184,171,173,174,161:163),"NSEI", 
                         ifelse(IPHC.Stat.Area %in% c(142:144,152,153),"SSEI",
                                ifelse(-MidLon.fished >= 137 & -MidLon.fished <=140, "EYKT", 
                                       ifelse(-MidLon.fished <= 137 & MidLat.fished >= 57.5,"NSEO",
                                              ifelse(-MidLon.fished <= 137 & MidLat.fished <57.5 & MidLat.fished >= 56,"CSEO",
                                                            ifelse(-MidLon.fished <= 137 & MidLat.fished < 56,"SSEO",NA)))))), 
           In.Out = ifelse(SEdist %in% c("NSEO","CSEO","SSEO","EYKT"),"SEO",
                           ifelse(SEdist %in% c("NSEI","SSEI"),"SEI",NA))) %>% 
    rename(YE.obs = Number.Observed,
           AvgDepth.fm = AvgDepth..fm.,
           Year = Year.x) %>%
    mutate(depth_bin = cut(AvgDepth.fm, breaks = seq(0,400,50),
                           labels = paste (seq(50,400,50))),
           depth_m = AvgDepth.fm * 1.8288,
           depth_bin_m = cut(depth_m, breaks = seq(0,800,50),
                             labels = paste (seq(50,800,50)))) %>% 
    select(-c(Row.number.x,Row.number.y))

  #Bring in and join hook adjustment factor to control for hook saturation when calculating CPUE
  # Note from PJ: The data was already posted: https://www.iphc.int/data/fiss-survey-raw-survey-data/
  hadj<-read_csv("./Data_processing/Data/IPHC_raw/iphc-2023-fiss-hadj-20231031_rke.csv") %>% 
    dplyr::rename_all(funs(make.names(.))) %>% 
    filter(IPHC.Reg.Area %in% c("2C", "3A"), 
           IPHC.StatArea <= 200)  %>%
    select(Stlkey = stlkey, Year, Station, AdjHooks.Observed = Hooks.Observed, h.adj)
    str(hadj)

  Survey <- left_join(Survey_prep,hadj,by=c("Year","Stlkey","Station")) %>% 
    mutate(YE.obs = ifelse(is.na(YE.obs),0,YE.obs),
           # RKE used AdjHooks.Observed from hadj file to fill in missing 
           # hook observed fields instead of static "140" from Phil's code which is unclear 
           # where that came from given a different number of skates are sampled 
           # he said 20 hooks observed x 7 skates sampled but some sets are a total of 6 skates set/hauled
           # so those would be 120 hooks instead of 140. Regardless, if you use the Hooks Observed
           # from the hadj file, it's IPHC data for that Stlkey/Station! So I assume it's correct. 
           HooksObserved = ifelse(is.na(HooksObserved), AdjHooks.Observed, HooksObserved), 
           HooksRetrieved = as.numeric(HooksRetrieved), 
           HooksObserved = as.numeric(HooksObserved), 
           YE.exp = HooksRetrieved/HooksObserved)
  
{y<-sample(unique(Survey$Year),1)
d<-sample(unique(Survey$depth_bin[Survey$Year==y]),1)
s<-sample(unique(Survey$Station[Survey$Year==y & Survey$depth_bin==d]),1)
Dat<-Survey[Survey$Year == y & Survey$depth_bin == d & Survey$Station == s,]
unique(Dat); nrow(Dat)}

str(Survey); head(Survey,10)
str(HA.Harv)
mean(Survey$AvgDepth.fm)

write.csv(Survey,"./Data_processing/Data/IPHC_survey_1998-2024_rke.csv")

#================================================================================
# Look at some depth stuff and identify stations that have seen yelloweye...
hist(Survey$BeginDepth..fm.)
hist(Survey$EndDepth..fm.)
hist(Survey$AvgDepth.fm)

plot(Survey$YE.obs ~ Survey$AvgDepth.fm)
plot(Survey$O32.Pacific.halibut.count ~ Survey$AvgDepth.fm)
plot(Survey$U32.Pacific.halibut.count ~ Survey$AvgDepth.fm)

#Some data for PWS YE assessment: 
with(Survey, table(YE.obs,depth_bin_m))
Survey %>% group_by(depth_bin_m) %>% summarize(ye_obs = sum(YE.obs),
                                               n_stations = n_distinct(Station)) %>%
  mutate(YE_per_station = ye_obs / n_stations) -> ye_by_depth_m
#write.csv(ye_by_depth_m,"Data_processing/Data/SE_YE_FISS_depth_forPWS.csv")

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

hist(Station.sum$mean.YE, breaks=50); 1-(nrow(Station.sum[Station.sum$mean.YE == 0,])/nrow(Station.sum))

ggplot(Station.sum,aes(Depth, mean.YE)) +
  geom_point() + geom_smooth(linewidth = 1.5, se = FALSE, linetype=1, alpha=0.75) +
  ylab("Mean number of yelloweye encountered at FISS stations") +
  xlab("Depth (fathoms)") + theme_bw() +
  scale_x_continuous(breaks=seq(0,350,25)) 
  
ggsave(paste0("output/FISS_YE_x_depth_rke_", YEAR, ".png"), dpi=300,  height=5, width=5, units="in") 

ggplot(Station.sum %>% filter(mean.YE > 0),
       aes(Depth, mean.YE)) +
  geom_point() + geom_smooth(size = 1.5, se = FALSE, linetype=1, alpha=0.75) +
  ylab("Mean number of yelloweye encountered at FISS stations") +
  xlab("Depth (fathoms)") + theme_bw() +
  scale_x_continuous(breaks=seq(0,350,25))

ggsave(paste0("output/FISS_YE_x_depth_non0_rke_", YEAR, ".png"), dpi=300,  height=5, width=5, units="in") 

ye_caught_prop<-percent(1-(nrow(Station.sum[Station.sum$mean.YE == 0,])/nrow(Station.sum)))

ggplot(Station.sum) +
  geom_histogram(aes(mean.YE), binwidth=0.25, color="grey", fill="darkgrey") +
  annotate("text", y=40, # 
           x=12.5, color="darkorange",
           label=c(paste0("Yelloweye caught at least once at ",ye_caught_prop," of stations"))) +
  annotate("text", y=90, 
           x=12.5, color="black",fontface =2, 
           label=c(paste0("All SEO FISS stations"))) +
  ylab("Number of FISS stations") +
  xlab("") + theme_bw() -> all_staions_hist

ggplot(Station.sum %>% filter(mean.YE > 0)) +
  geom_histogram(aes(mean.YE), binwidth=0.25, color="grey", fill="darkgrey") +
  annotate("text", y=15, 
           x=12.5, color="black",fontface =2,
           label=c(paste0("SEO FISS stations that encountered yelloweye at least once"))) +
  ylab("Number of FISS stations") +
  xlab("mean number of yelloweye encountered per year") + theme_bw() -> non0_stations_hist

ggarrange(all_staions_hist, non0_stations_hist,
          nrow=2)

ggsave(paste0("output/FISS_YE_histograms_rke_", YEAR, ".png"), dpi=300,  height=5, width=5, units="in")

hist(Station.sum$mean.YEcpue, breaks=25); nrow(Station.sum[Station.sum$mean.YEcpue == 0,])/nrow(Station.sum)
hist(Station.sum$noYE.count, breaks=25)
hist(Station.sum$prop.0, breaks=25)
hist(Station.sum$prop.0[Station.sum$prop.0<1], breaks=25)
max(Station.sum$prop.0[Station.sum$prop.0<1])

blanks.stations<-Station.sum$Station[Station.sum$mean.YEcpue == 0]
All.stations <- Station.sum$Station
YE.stations<-Station.sum$Station[Station.sum$mean.YEcpue != 0]
YE.stations_10p<-Station.sum$Station[Station.sum$prop.0 < 0.9]
YE.stations_20p<-Station.sum$Station[Station.sum$prop.0 < 0.8]
YE.stations_25p<-Station.sum$Station[Station.sum$prop.0 < 0.75]
YE.stations_40p<-Station.sum$Station[Station.sum$prop.0 < 0.6]

length(All.stations);length(YE.stations);length(YE.stations_10p);length(YE.stations_20p);length(YE.stations_25p);length(YE.stations_40p)

plot_stations<-st_as_sf(Station.sum,coords = c("Lon","Lat"))

plot_stations<-plot_stations %>% st_set_crs(4326)

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

ggplot(plot_stations %>% filter(Station %in% c(use_stations))) + 
  geom_sf(aes(color = prop.0, size = prop.0)) +
  scale_color_viridis_c(option = "plasma",begin = 0)
#==============================================================================================
## Load the port samples that were downloaded from oceansAK; 
# Need Yelloweye weights to get wcpue estimates
 
  Port<-port.bio(2024)
  str(Port)
  Port$Year<-as.integer(Port[,1])
  unique(Port$Groundfish.Management.Area.Code)  #EYAK = EYKT
  Port<-Port[Port$Groundfish.Management.Area.Code != "",]
  Port<-subset(Port, !is.na(Weight.Kilograms))
  
  #Rhea edited so that samples are Random. Want to avoid the select samples for maturity/vonB that were collected 
  Port.rand <- Port %>%
    filter(Sample.Type == "Random") 
  
  uYEkg <- Port.rand %>% 
    group_by(Year,Groundfish.Management.Area.Code) %>%
    summarize(uYEkg = mean(Weight.Kilograms, na.rm=TRUE),
              vYEkg = var(Weight.Kilograms, na.rm=TRUE),
              NyeKG = length(Weight.Kilograms))
  uYEkg<-as.data.frame(uYEkg); uYEkg
  
  Survey$mean.YE.kg<-NA
  Survey$var.YE.kg<-NA
  Survey$N.YE.kg<-NA
  
  
  ## Get best data for average weights and add columns to survey data
  # Rhea updated this so that we get 200 samples per year per management area 
  # looking at recent years (2022-2024) we have very few samples. For example,
  # from EYKT there are 3 boats sampled in 2022, 3 boats in 2023, and 2 boats in 2024. I don't think 
  # these are representative. For our port sampling goals, we aim for 550 per year per management area, 
  # with 50-75 samples per boat to spread distribution of samples across time and area. 
  Years<-unique(Survey$Year)
  GFMA<-unique(Survey$SEdist)
  { 
  for (i in Years){    #i<-Years[1]
    for (j in GFMA){   #j<-GFMA[1]
      
      P<-Port.rand[Port.rand$Groundfish.Management.Area.Code == j,]
      #reach back from year i and see if you can get 150 weights...rhea bumped this to 200
      Nw<-nrow(P[P$Year == i,])

      k<-0
      while (Nw < 200 & i-k >= min(Port.rand$Year) ){         #go back as far as needed to get 200 weights... 
        k<-k+1
        Nw<-nrow(P[P$Year >=i-k & P$Year <= i,])   #P[P$Year == 1990:2002,]
      }
      m<-0
      while (Nw < 200 & i+m <= max(Port.rand$Year)) {                        #go forward if you failed to find enough weights... 
        m<-m+1                                                        
        Nw<-nrow(P[P$Year >=i-k & P$Year <= i+m,])   #P[P$Year == 1986,]
      }
      #get average weights of YE...
      Sample<-P$Weight.Kilograms[P$Year >= i-k & P$Year <= i+m]
      
      #PHIL'S CODE FAILS HERE DUE TO WEIRD SUBSETTING CODE - FIXED IN FOLLOWING 3 LINES
      # Survey[,"mean.YE.kg"][Survey$Year == i & Survey$SEdist == j]<-mean(Sample)
      # Survey[,"var.YE.kg"][Survey$Year == i & Survey$SEdist == j]<-var(Sample)
      # Survey[,"N.YE.kg"][Survey$Year == i & Survey$SEdist == j]<-length(Sample)
      Survey$mean.YE.kg[Survey$Year == i & Survey$SEdist == j]<-mean(Sample)
      Survey$var.YE.kg[Survey$Year == i & Survey$SEdist == j]<-var(Sample)
      Survey$N.YE.kg[Survey$Year == i & Survey$SEdist == j]<-length(Sample)
    }
  }
}
Survey$mean.YE.kg

Survey %>% group_by(SEdist,Year) %>%
  dplyr::summarise(mean_wt = mean(mean.YE.kg)) %>%
  ggplot() + 
  geom_point(aes(Year,mean_wt,col=SEdist)) +
  geom_line(aes(Year,mean_wt,col=SEdist)) 


#=============================================================================
# FUNCTION for generating cpue estimates from IPHC surveys for different management areas

YEHA.fxn<-function(Survey=Survey, Area="SEdist",Deep=250, Shallow=0,  nboot=1000){
  col.ref<-which(colnames(Survey)==Area)
  
  IPHC.cpue<-data.frame()
  Subs<-unique(Survey$SEdist) # FIXED FROM PHIL'S CODE - DIDN'T SPECIFY COLUMN
  Years<-unique(Survey$Year)
  
  j<-1
  for (y in Years) {  #y<-Years[26]
    for (s in Subs){  #s<-Subs[3]
      Dat<-Survey[Survey$Year == y & 
                    Survey[,col.ref] == s &
                    Survey$AvgDepth.fm > Shallow & 
                    Survey$AvgDepth.fm < Deep &
                    Survey$Eff == "Y",]
      YE_w <-unique(Dat$mean.YE.kg)
      if (nrow(Dat)>0){
        Stations<-unique(Dat$Station)
        CPUEi<-vector()
        WCPUEi<-vector()
        i<-1
        for (st in Stations){    #st<-Stations[1]   length(Stations)  st<-1
          Stat.Dat<-Dat[Dat$Station == st,] #; Stat.Dat
          #debug
          #if (nrow(Stat.Dat) > 1){aaa} else {}
          CPUE<-mean(Stat.Dat$YE.obs/Stat.Dat$HooksObserved)
          #hook adjustment factor
          CPUE<-CPUE*unique(Stat.Dat$h.adj)
          WCPUE<-CPUE*YE_w
          
          if (CPUE == 0){
            C<-0
          } else {
            C<-CPUE*Stat.Dat$HooksRetrieved
          }
          
          CPUEi[i]<-CPUE
          WCPUEi[i]<-WCPUE
          i<-i+1
        }
        
        Out<-data.frame()
        #CPUE.out<-data.frame()
        for (ii in 1:nboot){ #nboot<-1000  ii<-1
          Resamp3<-sample(CPUEi,length(CPUEi),replace=T)
          Out[ii,"CPUE"]<-mean(Resamp3, na.rm=T)
          Out[ii,"CPUE.var"]<-var(Resamp3, na.rm=T)
          Resamp4<-sample(WCPUEi,length(WCPUEi),replace=T)
          Out[ii,"WCPUE"]<-mean(Resamp4, na.rm=T)
          Out[ii,"WCPUE.var"]<-var(Resamp4, na.rm=T)
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
        IPHC.cpue[j,"WCPUE.mean"]<-mean(WCPUEi)
        IPHC.cpue[j,"WCPUE.bootmean"]<-unname(quantile(Out$WCPUE,c(0.5)))  #mean WCPUE from Tribuzio
        IPHC.cpue[j,"WCPUE.lo95ci"]<-unname(quantile(Out$WCPUE,c(0.025)))
        IPHC.cpue[j,"WCPUE.hi95ci"]<-unname(quantile(Out$WCPUE,c(0.975)))
        IPHC.cpue[j,"WCPUE.var"]<-var(Out$WCPUE)
        IPHC.cpue[j,"WCPUE.cv"]<-sd(Out$WCPUE)/mean(Out$WCPUE)

        j<-j+1
      } else {}
    }
  }
  
  #couple of plots in output
  #Sub.cnt<-length(Subs)
  #cols<-brewer.pal(min(Sub.cnt,8),"Dark2")
  
  #if (Sub.cnt < 8) {par(mfrow=c(1,1))}
  #if (Sub.cnt > 8 & Sub.cnt < 17) {par(mfrow=c(2,1))}
  #if (Sub.cnt > 16) {par(mfrow=c(3,1))}
  
  #nplot<-ceiling(Sub.cnt/8)
  #cols<-brewer.pal(length(Subs),"Set3")
  #for (n in 1:nplot){  #n<-1
  #  i<-1
  #  pSubs<-Subs[(1+(n-1)*8):(8+(n-1)*8)]
  #  pSubs<-pSubs[!is.na(pSubs)]
    
  #  for (s in pSubs) {   #s<-pSubs[1]
  #    Pdat<-IPHC.cpue[IPHC.cpue$mngmt.area == s,]
  #    if (s == pSubs[1]){
  #      plot(Pdat$CPUE.mean ~ Pdat$Year, col=cols[i], 
  #           ylim=c(0,max(IPHC.cpue$CPUE.hi95ci)), type="l",
  #           ylab="CPUE", xlab="Year", main=Area)
  #      Y1<-Pdat$CPUE.lo95ci; Y2<-Pdat$CPUE.hi95ci
  #      polygon(c(Years,rev(Years)),c(Y1,rev(Y2)),
  #              col=adjustcolor(cols[i],alpha.f=0.2),border=NA)
  #      i<-i+1
  #    } else {
  #      Y1<-Pdat$CPUE.lo95ci; Y2<-Pdat$CPUE.hi95ci
  #      if (length(Y1)==1){}else {
  #        polygon(c(Years,rev(Years)),c(Y1,rev(Y2)),
  #                col=adjustcolor(cols[i],alpha.f=0.2),border=NA)
  #        i<-i+1
  #      }
        
  #    }
  #  }
  #  i<-1
  #  for (s in pSubs) {   #s<-Subs[2]
  #    Pdat<-IPHC.cpue[IPHC.cpue$mngmt.area == s,]
  #    lines(Pdat$CPUE.mean ~ Pdat$Year, col=cols[i], type="l",
  #          lwd=1.5)
  #    i<-i+1
  #  }
  #  legend(x="topleft", cex=0.8,
  #         legend=c(pSubs), bty="n",
  #         col=cols,text.col=cols)
  #}
  
  return(IPHC.cpue)
}


colnames(Survey)

#Decision: which stations to use to calculate CPUE??? 

# Note Oct 2023: Port function not set up for inside waters yet so just do outside for now: 
#Survey_40p<-Survey %>% filter(Station %in% c(YE.stations_40p))
Survey_all<-Survey %>% filter(Station %in% c(All.stations) &
                                 SEdist %in% c("EYKT","NSEO","CSEO","SSEO"))

Survey_non0<-Survey %>% filter(Station %in% c(YE.stations) &
                                 SEdist %in% c("EYKT","NSEO","CSEO","SSEO"))

#IPHC.reg.area<- YEHA.fxn(Area="IPHC.Reg.Area",Deep=max(Survey$AvgDepth.fm), Shallow=0, nboot=1000) 
#IPHC.stat.area<- YEHA.fxn(Area="IPHC.Stat.Area",Deep=max(Survey$AvgDepth.fm), Shallow=0, nboot=1000) 
SE.subdistricts<- YEHA.fxn(Survey=Survey_non0, Area="SEdist",Deep=250, Shallow=0, nboot=1000)
#SE.subdistricts_40p<- YEHA.fxn(Survey=Survey_40p, Area="SEdist",Deep=250, Shallow=0, nboot=1000)

str(SE.subdistricts); unique(SE.subdistricts$mngmt.area)

SE.subdistricts %>% filter(Year == 2023)

IPHC.index<-SE.subdistricts #[SE.subdistricts$mngmt.area != "NSEI" & SE.subdistricts$mngmt.area != "SSEI",NA]
str(IPHC.index)

ggplot(IPHC.index, aes(x=Year)) +
  geom_point(aes(y=CPUE.mean),size=2) +
  geom_point(aes(y=WCPUE.mean),size=2,col="blue") +
  geom_errorbar(aes(ymin = CPUE.lo95ci, ymax= CPUE.hi95ci),
                col="black", alpha=0.5) +
  geom_errorbar(aes(ymin = WCPUE.lo95ci, ymax= WCPUE.hi95ci),
                col="blue", alpha=0.5) +
  facet_wrap(~mngmt.area) +
  xlab("\nYear") +
  ylab("Yelloweye CPUE / WCPUE n/hook") +
  #ylab("Density (yelloweye rockfish/kmsq)") +
  scale_y_continuous(label=comma) +
  scale_x_continuous(breaks=seq(1995,2025,5)) + 
  theme (axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
         panel.grid.minor = element_blank()) +
  labs(title =~ atop("Yelloweye cpue in IPHC FISS",scriptstyle("stations that encountered yelloweye at least once")))

ggsave(paste0("output/YE_IPHC_cpue_non0_boot_rke_", YEAR, ".png"), dpi=300,  height=5, width=5, units="in")

#This output was used in the REMA model; use the Tweedie model output now!
write.csv(IPHC.index,paste0(here::here(), "/Data_processing/Data/IPHC_raw/IPHC.cpue.SEO_non0_rke_",YEAR,".csv"))

#get average number of stations in each area 
IPHC.index %>% group_by(mngmt.area) %>%
  summarize(mean_no_stations = mean(no.stations),
            min_no_stations = min(no.stations),
            max_no_stations = max(no.stations))

# Code to check what filtering more stations does. Blocked out now and not necessary
# but saved here if someone wants to revisit this. 

#IPHC.index_40p<-SE.subdistricts_40p #[SE.subdistricts_40p$mngmt.area != "NSEI" &                               SE.subdistricts_40p$mngmt.area != "SSEI",]

#ggplot(IPHC.index_40p, aes(x=Year)) +
#  geom_point(aes(y=CPUE.mean),size=2) +
#  geom_errorbar(aes(ymin = CPUE.lo95ci, ymax= CPUE.hi95ci),
#                col="black", alpha=0.5) +
#  facet_wrap(~mngmt.area) +
#  xlab("\nYear") +
#  ylab("Yelloweye CPUE n/hook") +
  #ylab("Density (yelloweye rockfish/kmsq)") +
#  scale_y_continuous(label=comma) +
#  scale_x_continuous(breaks=seq(1995,2025,5)) + 
#  theme (axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
#         panel.grid.minor = element_blank()) +
#  ggtitle("Yelloweye cpue in IPHC FISS, stations that encountered yelloweye at least 40% of years")

#ggsave(paste0("Figures/YE_IPHC_cpue_min40percentYE_", YEAR, ".png"), dpi=300,  height=5, width=5, units="in")

#This output is used in the REMA model if being more restrictive 
#write.csv(IPHC.index_40p,paste0("Data_processing/Data/IPHC.cpue.SEO_min40percentYE_",YEAR,".csv"))

#comp<-rbind(IPHC.index %>% mutate(yelloweye_presence = "at least once"),
#            IPHC.index_40p %>% mutate(yelloweye_presence = "40% of years"))

#ggplot(comp, aes(x=Year)) +
#  geom_point(aes(y=CPUE.mean, col = yelloweye_presence),
#             size=2, position = position_dodge(width=0.25)) +
#  geom_errorbar(aes(ymin = CPUE.lo95ci, ymax= CPUE.hi95ci, col=yelloweye_presence),
#                alpha=0.5,
#                position = position_dodge(width=0.25)) +
#  facet_wrap(~mngmt.area) +
#  xlab("\nYear") +
#  ylab("Yelloweye CPUE n/hook") +
  #ylab("Density (yelloweye rockfish/kmsq)") +
#  scale_y_continuous(label=comma) +
#  scale_x_continuous(breaks=seq(1995,2025,5)) + 
#  theme (axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
#         panel.grid.minor = element_blank(),
#         legend.position = "top") +
#  ggtitle("Yelloweye cpue in IPHC FISS")

#ggsave(paste0("Figures/YE_IPHC_cpue_method_comparison_", YEAR, ".png"), dpi=300,  height=5, width=5, units="in")


################################################################################
## Model based CPUE estimate using Tweedie model
##
################################################################################
# We will still restrict stations to those that have encountered YE at least once
# in the time series and not consider stations below 250 fathoms
Survey_non0<-Survey %>% filter(Station %in% c(YE.stations) &
                                 SEdist %in% c("EYKT","NSEO","CSEO","SSEO"))

IPHC_tweed <- Survey_non0 %>% filter(AvgDepth.fm <= 250) %>% 
  group_by(Year,Station) %>%
  mutate(CPUE = h.adj*YE.obs/HooksObserved,
         WCPUE = CPUE*mean.YE.kg) %>%
  select(Year,Station,SEdist,IPHC.Reg.Area, IPHC.Stat.Area, IPHC.Charter.Region,        
         Purpose.Code, Date, Eff, 
         Lat = BeginLat, Lon = BeginLon, Depth = AvgDepth.fm, 
         Soak = Soak.time..min.., Temp.C, Pres = Max.Pressure..db., 
         pH, Sal = Salinity.PSU, O2_1 = Oxygen_ml,
         O2_2 = Oxygen_umol, Oxygen_sat,
         CPUE, WCPUE) %>% 
  mutate(Year = as.factor(Year),
         SEdist = as.factor(SEdist),
         Soak = as.numeric(Soak),
         Temp.C = as.numeric(Temp.C),
         Sal = as.numeric(Sal),
         O2_1 = as.numeric(O2_1),
         O2_2 = as.numeric(O2_2),
         Oxygen_sat = as.numeric(Oxygen_sat)) %>%
  unique() %>% data.frame()

head(data.frame(IPHC_tweed),30)
colnames(IPHC_tweed)

str(IPHC_tweed)

IPHC_tweed %>% filter(Year == 2021 & Station == 3211)

nrow(IPHC_tweed %>% filter(is.na(Soak)))
nrow(IPHC_tweed %>% filter(is.na(Depth)))
nrow(IPHC_tweed %>% filter(is.na(Temp.C)))
nrow(IPHC_tweed %>% filter(is.na(Pres)))
nrow(IPHC_tweed %>% filter(is.na(pH)))
nrow(IPHC_tweed %>% filter(is.na(Sal)))
nrow(IPHC_tweed %>% filter(is.na(O2_1)))
nrow(IPHC_tweed %>% filter(is.na(Lat)))
nrow(IPHC_tweed %>% filter(is.na(Lon)))

################################################################################
## Visually see how wcpue is related to variables:

# Depth
ggplot(IPHC_tweed, aes(Depth, WCPUE)) + geom_point(shape = 20) + 
  geom_smooth(size = 2, se = FALSE)

# Soak
ggplot(IPHC_tweed, aes(Soak, WCPUE)) + geom_point(shape = 20) + 
  geom_smooth(size = 2, se = FALSE)

# Temp
hist(IPHC_tweed$Temp.C)
ggplot(IPHC_tweed, aes(Temp.C, WCPUE)) + geom_point(shape = 20) + 
  geom_smooth(size = 2, se = FALSE)

# Pressure
ggplot(IPHC_tweed, aes(Pres, WCPUE)) + geom_point(shape = 20) + 
  geom_smooth(size = 2, se = FALSE)

# pH
ggplot(IPHC_tweed, aes(pH, WCPUE)) + geom_point(shape = 20) + 
  geom_smooth(size = 2, se = FALSE)

# Salinity
ggplot(IPHC_tweed, aes(Sal, WCPUE)) + geom_point(shape = 20) + 
  geom_smooth(size = 2, se = FALSE)

# Oxygen 1
ggplot(IPHC_tweed, aes(O2_1, WCPUE)) + geom_point(shape = 20) + 
  geom_smooth(size = 2, se = FALSE)

# Oxygen 2
ggplot(IPHC_tweed, aes(O2_2, WCPUE)) + geom_point(shape = 20) + 
  geom_smooth(size = 2, se = FALSE)

# Oxygen saturation
ggplot(IPHC_tweed, aes(Oxygen_sat, WCPUE)) + geom_point(shape = 20) + 
  geom_smooth(size = 2, se = FALSE)

# LatLon
ggplot(IPHC_tweed, aes(Lat, WCPUE)) +
  geom_smooth(method = 'loess', span = 1, se = FALSE)

ggplot(IPHC_tweed, aes(Lon, WCPUE)) +
  geom_smooth(method = 'loess', span = 1, se = FALSE)

################################################################################
## Check for colinearity

IPHC_tweed %>% 
  select(Depth, Soak, Lat, Lon) %>% 
  GGally::ggpairs() 

#some strong correlations here...
# almost everything correlated with depth! 
# Temp, Pres pH and Salinity all correlated with each other... :( 

# probably best to just consider latlon, depth, and soak  colnames(IPHC_tweed)
IPHC_tweed <- IPHC_tweed %>% 
  select(Year,Station,SEdist,IPHC.Reg.Area, IPHC.Stat.Area, IPHC.Charter.Region,        
         Purpose.Code, Date, Eff, 
         Lat, Lon, Depth, 
         Soak, 
         CPUE, WCPUE)


# Look at distribution of CPUE data
ggplot(IPHC_tweed, aes(WCPUE)) + geom_density(alpha = 0.4, fill = 4)
ggplot(IPHC_tweed, aes(log(WCPUE+0.01))) + geom_density(alpha = 0.4, fill = 4)
ggplot(IPHC_tweed, aes(log(WCPUE+0.1*mean(IPHC_tweed$WCPUE,na.rm=T)))) + geom_density(alpha = 0.4, fill = 4)

#need complete data sets for running models: 
nrow(IPHC_tweed)
nrow(IPHC_tweed[complete.cases(IPHC_tweed),])

#only use complete cases... 
fulldat<-IPHC_tweed[complete.cases(IPHC_tweed),]

m0 <- gam(WCPUE ~ Year * SEdist, data=fulldat, gamma=1.4, family=tw(), method = "REML")
m.depth <- gam(WCPUE ~ Year * SEdist + s(Depth, k=4), data=fulldat, gamma=1.4, family=tw(), method = "REML")

m.soak <- gam(WCPUE ~ Year * SEdist + s(Soak, k=4), data=fulldat, gamma=1.4, family=tw(), method = "REML")
m.ll <- gam(WCPUE ~ Year * SEdist + te(Lon, Lat), data=fulldat, gamma=1.4, family=tw(), method = "REML")
m.depth_soak <- gam(WCPUE ~ Year * SEdist + s(Depth, k=4) + s(Soak, k=4), 
                    data=fulldat, gamma=1.4, family=tw(), method = "REML")
m.depth_ll <- gam(WCPUE ~ Year * SEdist + s(Depth, k=4) + te(Lon, Lat), 
                  data=fulldat, gamma=1.4, family=tw(), method = "REML")
m.soak_ll <- gam(WCPUE ~ Year * SEdist + s(Soak, k=4) + te(Lon, Lat) ,
                 data=fulldat, gamma=1.4, family=tw(), method = "REML")
m.global <- gam(WCPUE ~ Year * SEdist + s(Depth, k=4) + s(Soak, k=4) + te(Lon, Lat) , 
                data=fulldat, gamma=1.4, family=tw(), method = "REML")

model.list<-list(m0,#m0.trip,
                 m.depth,m.soak,m.ll,
                 m.depth_soak,m.depth_ll,m.soak_ll,m.global
                 )
names(model.list)<-c("m0",#"trip",
                     "depth","soak","latlon","depth+soak","depth+latlon",
                     "soak+latlon","global")
modsum0<-data.frame(); j<-1
for (i in model.list) {
  #mod<-i
  modsum0[j,"model"]<-names(model.list[j])
  modsum0[j,"aic"]<-AIC(i)
  modsum0[j,"bic"]<-BIC(i)
  modsum0[j,"qaic"]<-QAIC(i, chat= 1/(1-summary(i)$r.sq))
  modsum0[j,"dev"]<-summary(i)$dev.expl
  modsum0[j,"rsq"]<-summary(i)$r.sq
  modsum0[j,"dev_exp"]<-summary(i)$dev.expl-summary(m0)$dev.expl
  j<-j+1
}

modsum0 %>% arrange(aic)
modsum0 %>% arrange(qaic)
modsum0 %>% arrange(bic)
modsum0 %>% arrange(-dev)  
modsum0 %>% arrange(-rsq) 

plot(m.global, page = 1, shade = TRUE, resid = TRUE, all = TRUE)
summary(m.global)

# No residual patterns, but may be some outliers
plot(fitted(m.global), resid(m.global))
abline(h = 0, col = "red", lty = 2)

plot(m.depth_soak)
plot(fitted(m.depth_soak), resid(m.depth_soak))
# Check for outliers
which(fitted(m.global) < -1.5)   

mod_std_tweed <- m.global

################################################################################
## get predicted values:

# Cant have interaction between 

#cpue_dat<-cpue_nom %>% select(-c(Gear,Drifts,Stat, Depth))
cpue_dat<-IPHC_tweed[complete.cases(IPHC_tweed),]
cpue_dat <- as.data.frame(IPHC_tweed)

# Predictions ----

# Data set of average variables to predict CPUE from:
std_dat_tweed <- expand.grid(Year = as.factor(unique(cpue_dat$Year)),
                             SEdist = as.factor(unique(cpue_dat$SEdist)), #table(cpue_dat$Gear)
                             Depth = mean(cpue_dat$Depth),
                             Soak = mean(cpue_dat$Soak, na.rm=T), 
                             #Lat = mean(cpue_dat$Lat),
                             #Lon = mean(cpue_dat$Lon),
                             dum = 1,
                             dumstat = 1) %>%
  mutate(Lat = case_when(SEdist == "EYKT" ~ mean(cpue_dat$Lat[cpue_dat$SEdist == "EYKT"]),
                         SEdist == "NSEO" ~ mean(cpue_dat$Lat[cpue_dat$SEdist == "NSEO"]),
                         SEdist == "CSEO" ~ mean(cpue_dat$Lat[cpue_dat$SEdist == "CSEO"]),
                         SEdist == "SSEO" ~ mean(cpue_dat$Lat[cpue_dat$SEdist == "SSEO"])),
         Lon = case_when(SEdist == "EYKT" ~ mean(cpue_dat$Lon[cpue_dat$SEdist == "EYKT"]),
                         SEdist == "NSEO" ~ mean(cpue_dat$Lon[cpue_dat$SEdist == "NSEO"]),
                         SEdist == "CSEO" ~ mean(cpue_dat$Lon[cpue_dat$SEdist == "CSEO"]),
                         SEdist == "SSEO" ~ mean(cpue_dat$Lon[cpue_dat$SEdist == "SSEO"])))

# Predict CPUE
pred_cpue <- predict(mod_std_tweed, std_dat_tweed, type = "link", se = TRUE)
#checking my code with Jane's... checks out :)
preds<-predict.gam(mod_std_tweed, type="response", std_dat_tweed, se = TRUE)
pred_cpue;preds

#Put the standardized CPUE and SE into the data frame and convert to
#backtransformed (bt) CPUE
alpha <- 0.05  # for a 95% confidence interval on bycatch and discard estimates
z <- qnorm(1 - alpha / 2)  # Z value for 95% CI

std_dat_tweed %>% 
  mutate(fit = pred_cpue$fit,
         se = pred_cpue$se.fit,
         lower_link = fit - z * se,
         upper_link = fit + z * se,
         lower = mod_std_tweed$family$linkinv(lower_link),
         upper = mod_std_tweed$family$linkinv(upper_link),
         #upper = fit + (2 * se),
         #lower = fit - (2 * se),
         bt_cpue = exp(fit),
         bt_upper = upper, #exp(upper),
         bt_lower = lower, #exp(lower),
         bt_se = (bt_upper - bt_cpue) / 2  ,
         bt_cv = bt_se/bt_cpue
  ) -> std_dat_tweed

# Nominal CPUE for comparison ----

IPHC_tweed %>%
  group_by(Year,SEdist) %>%
  do(data.frame(rbind(Hmisc::smean.cl.boot(.$WCPUE)))) %>%
  mutate(calc = "set lvl kg/hook",
         fsh_cpue = Mean,
         upper = Upper,
         lower = Lower,
         cv = (upper-lower)/1.96) -> fsh_sum_boot

std_dat_tweed %>% 
  select(Year, cpue = bt_cpue, upper = bt_upper, lower = bt_lower, SEdist) %>% 
  mutate(CPUE = "Tweedie index") %>% data.frame() -> plot_dat

plot_dat %>% data.frame() %>%
  ggplot(aes(group = 1)) + #geom_ribbon(aes(Year, ymin = lower, ymax = upper)) +
  geom_ribbon(aes(Year, ymin = lower, ymax = upper), fill = wes_palette("Darjeeling1")[c(5)],
              colour = NA, alpha = 0.5) +
  geom_point(aes(Year, cpue, colour = CPUE, shape = CPUE), size = 2, show.legend = F) +
  geom_line(aes(Year, cpue, colour = CPUE, group = CPUE), size = 1) +
  facet_wrap(~ SEdist, scales = "free") +
  scale_colour_manual(values =  wes_palette("Darjeeling1")[c(5,4)],
                      aesthetics = c("colour","fill"), name = "IPHC CPUE") +
  labs(x = "Year", y = "FISS CPUE (round kg / hook)\n") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5)) + #c(0.85, 0.9)) +
  expand_limits(y = 0)

ggsave(paste0("output/IPHC_cpue_tweedie_rke_",YEAR,".png"), dpi=300, height=6, width=6, units="in")

#-------------------------------------------------------------------------------
# Compare Tweedie and bootstrap index if you want... 
boot_index <- IPHC.index %>% 
  mutate(SEdist = mngmt.area, fsh_cpue = WCPUE.mean, upper = WCPUE.hi95ci,
         lower = WCPUE.lo95ci, cv = (upper-lower)/1.96) %>%
  select(Year, SEdist, fsh_cpue, lower, upper, cv)

# Compare predicted cpue from gam to nominal cpue
names(wes_palettes)

boot_index %>%
  select(Year, cpue = fsh_cpue, upper, lower, SEdist) %>% 
  mutate(CPUE = "Bootstrap index",
         Year = as.factor(Year)) %>% 
  bind_rows(std_dat_tweed %>% 
              select(Year, cpue = bt_cpue, upper = bt_upper, lower = bt_lower, SEdist) %>% 
              mutate(CPUE = "Tweedie index")) %>% 
  mutate(Year = as.numeric(as.character(Year))) %>%
  ggplot() +
  geom_ribbon(aes(Year, ymin = lower, ymax = upper, fill = CPUE), 
              colour = NA, alpha = 0.3) +
  geom_point(aes(Year, cpue, colour = CPUE, shape = CPUE), size = 2, show.legend = F) +
  geom_line(aes(Year, cpue, colour = CPUE, group = CPUE), size = 1) +
  facet_wrap(~ SEdist, scales = "free") +
  scale_colour_manual(values =  wes_palette("Darjeeling1")[c(5,4)],
                      aesthetics = c("colour","fill"), name = "IPHC CPUE") +
  labs(x = "Year", y = "FISS CPUE (round kg / hook)\n") +
  theme(legend.position = "bottom") + #c(0.85, 0.9)) +
  expand_limits(y = 0)

ggsave(paste0("output/IPHC_cpue_tweedie_boot_comp_rke_",YEAR,".png"), dpi=300, height=6, width=6, units="in")

# Save these values to combine with other indices we'll generate: 
boot_index %>%
  select(Year, SEdist, cpue = fsh_cpue, upper, lower, cv) %>% 
  mutate(CPUE = "Bootstrap index",
         Year = as.factor(Year)) %>% 
  bind_rows(std_dat_tweed %>% 
              select(Year, SEdist, cpue = bt_cpue, 
                     upper = bt_upper, lower = bt_lower, cv = bt_cv) %>% 
              mutate(CPUE = "Tweedie index")) %>% 
  #mutate(Year = as.numeric(as.character(Year)))) %>% 
  data.frame() -> IPHC_cpue_indices

IPHC_cpue_indices<-left_join(IPHC_cpue_indices,IPHC.index %>% mutate(Year = as.factor(Year)) %>%
                select(Year,SEdist = mngmt.area,no.stations),by=c("Year","SEdist"))

# SAVE all the indices for used in assessment models
write.csv(IPHC_cpue_indices,paste0("./Data_processing/Data/IPHC.cpue.SEO_tweed_boot_rke_",YEAR,".csv"))

