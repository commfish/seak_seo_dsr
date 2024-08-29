#############################################################################################
## SEO Subdistrict Biomass Recalculation
## Dec 2021
## For stock assessment
## Phil Joy
## Last Updated 8/29/24 - LSC
#############################################################################################

{library(dplyr)
library(boot)
library(ggplot2)
library(plotly)
library(scales)

library(tidyverse)
library(kableExtra)}

source("r_helper/Port_bio_function.R")

meanfun <- function(x, d) {
  return(mean(x[d]))
}

# 1) Get subdistrict density data
# 2) Get port sampling data for weights
# 3) recalculate biomass

YEAR<-2023
#===========================================================================================
# 1) Get subdistrict density data
Dens<-read.csv("Data_processing/Data/YE_Density_SEOsubdistricts.csv")
str(Dens)
Dens$Density<-as.numeric(Dens$Density); Dens$YE_abund<-as.integer(Dens$YE_abund)

# Get variance of density & abundance from reported cv
Dens$var.Dens<-(Dens$Density*Dens$CV)^2

Dens$var.Abund<-Dens$var.Dens*Dens$Area_km2^2
head(Dens)

dev.off()

ggplot(Dens, aes(x=Year)) +
  geom_point(aes(y=Density),size=2) +
  geom_errorbar(aes(ymin = Density-1.95*sqrt(var.Dens), ymax= Density+1.95*sqrt(var.Dens)),
                col="black", alpha=0.5) +
  facet_wrap(~Subdistrict) +
  xlab("\nYear") +
  ylab(bquote(Yelloweye~rockfish/km^2)) +
  #ylab("Density (yelloweye rockfish/kmsq)") +
  scale_y_continuous(label=comma) +
  scale_x_continuous(breaks=seq(1995,2025,5)) + 
  theme (axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
    panel.grid.minor = element_blank())
ggsave(paste0("Figures/YE_Density_", YEAR, ".png"), dpi=300,  height=4, width=4, units="in")
#===========================================================================================
# 2) Get port sampling data for weights
#last update of port sampling data was: 08/26/24 <-update this whenever you run this!!

{Port<-port.bio(2024)

unique(Port$Sample.Type)
Port.rand<-Port[Port$Sample.Type=="Random",]
Port.direct<-Port[Port$Project=="Commercial Longline Trip",]
Port.hal<-Port[Port$Project=="Commercial Halibut Longline",]; nrow(Port.hal)
}
#Port<-Port[!is.na(Port$GFMU),]
{Dens$mean.YE.kg_rand<-NA
Dens$var.YE.kg_rand<-NA
Dens$bootvar.YE.kg_rand<-NA
Dens$N.wts_rand<-NA
Dens$wt_yrs_rand<-NA

Dens$mean.YE.kg_dir<-NA
Dens$var.YE.kg_dir<-NA
Dens$bootvar.YE.kg_dir<-NA
Dens$N.wts_dir<-NA
Dens$wt_yrs_dir<-NA

Dens$mean.YE.kg_hal<-NA
Dens$var.YE.kg_hal<-NA
Dens$bootvar.YE.kg_hal<-NA
Dens$N.wts_hal<-NA
Dens$wt_yrs_hal<-NA}

Years<-unique(Dens$Year)
GFMA<-unique(Dens$Subdistrict)

## Get best data for average weights and add columns to survey data
for (i in Years){    #i<-Years[1]
  for (j in GFMA){   #j<-GFMA[1]
    P <- Port.rand %>% filter (GFMU == j); nrow(P)
    P <- P %>% filter(!is.na(Weight.Kilograms)); nrow(P)
    #reach back from year i and see if you can get 150 weights...nrow(P)
    Nw<-nrow(P[P$Year == i  & !is.na(P$Weight.Kilograms),])
    Nw<-nrow(P %>% filter(Year == i, !is.na(Weight.Kilograms)))
    
    #This segment of code is designed to ensure that you gather at least 75 
    #weight measurements for each combination of year and subdistrict. If there 
    #aren’t enough weights for the current year, it tries to collect additional 
    #weights by expanding the range of years both backward and forward.
    
    k<-0
    while (Nw < 75 & i-k >= min(Port.rand$Year) ){         #go back as far as needed to get 150 weights... 
      k<-k+1
      Nw<-nrow(P[P$Year >=i-k & P$Year <= i,])   #P[P$Year == 1990:2002,]
    }
    m<-0
    while (Nw < 75 & i+m <= max(Port.rand$Year)) {                        #go forward if you failed to find enough weights... 
      m<-m+1                                                        
      Nw<-nrow(P[P$Year >=i-k & P$Year <= i+m,])   #P[P$Year == 1986,]
    }
    #get average weights of YE...
    Sample<-P$Weight.Kilograms[P$Year >=i-k & P$Year <= i+m]
    Dens[,"wt_yrs_rand"][Dens$Year == i & Dens$Subdistrict == j]<-paste(i-k,"-",i+m,sep="")
    Dens[,"N.wts_rand"][Dens$Year == i & Dens$Subdistrict == j]<-Nw
    Dens[,"mean.YE.kg_rand"][Dens$Year == i & Dens$Subdistrict == j]<-mean(Sample)
    Dens[,"var.YE.kg_rand"][Dens$Year == i & Dens$Subdistrict == j]<-var(Sample)
    bootsmpl<-boot(Sample, statistic=meanfun, R=1000)
    Dens[,"bootvar.YE.kg_rand"][Dens$Year == i & Dens$Subdistrict == j]<-var(bootsmpl$t)
    #sd(bootsmpl$t)^2
    
    #directed fishery...
    P <- Port.direct %>% filter (GFMU == j); nrow(P)
    P <- P %>% filter(!is.na(Weight.Kilograms)); nrow(P)
    #reach back from year i and see if you can get 150 weights...nrow(P)
    Nw<-nrow(P[P$Year == i  & !is.na(P$Weight.Kilograms),])
    Nw<-nrow(P %>% filter(Year == i, !is.na(Weight.Kilograms)))
    
    k<-0
    while (Nw < 75 & i-k >= min(Port.direct$Year) ){         #go back as far as needed to get 150 weights... 
      k<-k+1
      Nw<-nrow(P[P$Year >=i-k & P$Year <= i,])   #P[P$Year == 1990:2002,]
    }
    m<-0
    while (Nw < 75 & i+m <= max(Port.direct$Year)) {                        #go forward if you failed to find enough weights... 
      m<-m+1                                                        
      Nw<-nrow(P[P$Year >=i-k & P$Year <= i+m,])   #P[P$Year == 1986,]
    }
    #get average weights of YE...
    Sample<-P$Weight.Kilograms[P$Year >=i-k & P$Year <= i+m]
    Dens[,"N.wts_dir"][Dens$Year == i & Dens$Subdistrict == j]<-Nw
    Dens[,"wt_yrs_dir"][Dens$Year == i & Dens$Subdistrict == j]<-paste(i-k,"-",i+m,sep="")
    Dens[,"mean.YE.kg_dir"][Dens$Year == i & Dens$Subdistrict == j]<-mean(Sample)
    Dens[,"var.YE.kg_dir"][Dens$Year == i & Dens$Subdistrict == j]<-var(Sample)
    bootsmpl<-boot(Sample, statistic=meanfun, R=1000)
    Dens[,"bootvar.YE.kg_dir"][Dens$Year == i & Dens$Subdistrict == j]<-var(bootsmpl$t)
    #
    ## Halibut bycatch
    P <- Port.hal %>% filter (GFMU == j); nrow(P)
    P <- P %>% filter(!is.na(Weight.Kilograms)); nrow(P)
    #reach back from year i and see if you can get 150 weights...nrow(P)
    Nw<-nrow(P[P$Year == i  & !is.na(P$Weight.Kilograms),])
    Nw<-nrow(P %>% filter(Year == i, !is.na(Weight.Kilograms)))
    
    k<-0
    while (Nw < 75 & i-k >= min(Port.hal$Year) ){         #go back as far as needed to get 150 weights... 
      k<-k+1
      Nw<-nrow(P[P$Year >=i-k & P$Year <= i,])   #P[P$Year == 1990:2002,]
    }
    m<-0
    while (Nw < 75 & i+m <= max(Port.hal$Year)) {                        #go forward if you failed to find enough weights... 
      m<-m+1                                                        
      Nw<-nrow(P[P$Year >=i-k & P$Year <= i+m,])   #P[P$Year == 1986,]
    }
    #get average weights of YE...
    Sample<-P$Weight.Kilograms[P$Year >=i-k & P$Year <= i+m]
    Dens[,"N.wts_hal"][Dens$Year == i & Dens$Subdistrict == j]<-Nw
    Dens[,"wt_yrs_hal"][Dens$Year == i & Dens$Subdistrict == j]<-paste(i-k,"-",i+m,sep="")
    Dens[,"mean.YE.kg_hal"][Dens$Year == i & Dens$Subdistrict == j]<-mean(Sample)
    Dens[,"var.YE.kg_hal"][Dens$Year == i & Dens$Subdistrict == j]<-var(Sample)
    bootsmpl<-boot(Sample, statistic=meanfun, R=1000)
    Dens[,"bootvar.YE.kg_hal"][Dens$Year == i & Dens$Subdistrict == j]<-var(bootsmpl$t)
  }
}
head(Dens[Dens$Subdistrict == "EYKT",],20)
Dens$var.YE.kg
Dens$bootvar.YE.kg

Dens[Dens$Subdistrict == "EYKT",]
colnames(Dens)
#Save weight data for tables in SAFE------------------------------------------
Wts<-Dens[,c(2,1,9,10,11,12,13)]; head(Wts)

#need to rerun to get back to 1984? 1980
PYears<-unique(Port.rand$Year)
GFMA<-unique(Dens$Subdistrict)
head(Port.rand)
nrow(Port.rand)
colnames(Port.rand)

Port.rand[Port.rand$GFMU == "EYKT" & Port.rand$Year == 1984,]

#These are ALL biological samples for yelloweye
unique(Port$Sample.Type)

Port %>% group_by(GFMU,Year) %>%
  filter(!is.na(GFMU),
         GFMU == "SSEO" | GFMU == "CSEO" | GFMU == "NSEO" | GFMU == "EYKT") %>% 
  summarise(mean.weight.kg = mean(Weight.Kilograms, na.rm=TRUE),
            var.weight = var(Weight.Kilograms, na.rm=TRUE),
            sd.weight = sqrt(var.weight),
            no.YE = sum(!is.na(Weight.Kilograms))) %>% 
  complete(Year = 1984:YEAR) -> Wt.smpl 

Wt.smpl$no.YE<-replace(Wt.smpl$no.YE,Wt.smpl$no.YE==0,NA)

opts <- options(knitr.kable.NA = "-")            
table<-Wt.smpl %>% #group_by(Year) %>% 
  kbl(digits=2, col.names = gsub("[.]", " ", names(head(Wt.smpl)))) %>% #kable_styling
  #kable_paper("hover",full_width=F)
  kable_classic(full_width=F, html_font = "Times", position="center")
  #kable_classic_2(full_width=F, position = "left")
webshot2::install_phantomjs()
save_kable(table,"Figures/Weight_Table_SAFE14p1.pdf")
#2024: can't get this pdf to run - need to come back to it

#These are samples that have been randomly sampled only - includes directed fishery and bycatch
unique(Port.rand$Project)

Port.rand %>% group_by(GFMU,Year) %>% 
  filter(!is.na(GFMU),
         GFMU == "SSEO" | GFMU == "CSEO" | GFMU == "NSEO" | GFMU == "EYKT") %>% 
  summarise(mean.weight.kg = mean(Weight.Kilograms, na.rm=TRUE),
            var.weight = var(Weight.Kilograms, na.rm=TRUE),
            sd.weight = sqrt(var.weight),
            no.YE = sum(!is.na(Weight.Kilograms))) %>% 
  complete(Year = 1984:YEAR) -> Wt.smpl.rand

write.csv(Wt.smpl,"Data_processing/Data/SEO_YE_mean_wts_all.csv")
write.csv(Wt.smpl.rand,"Data_processing/Data/SEO_YE_mean_wts_randomonly.csv")

ggplot(data = Dens %>% filter(Year > 2007), 
       aes(mean.YE.kg_hal,mean.YE.kg_dir)) + geom_point() +
  facet_wrap(~Subdistrict)
#NOTE - bycatch ye heavier than directed
ggplot(data = Dens %>% filter(Year > 2007), 
       aes(mean.YE.kg_hal,mean.YE.kg_rand)) + geom_point() +
  facet_wrap(~Subdistrict)

Dens$N.wts_hal[Dens$Year > 2007]
#====================================================================================
# Calculate biomass

#Calculate biomass in kilograms (yelloweye abundance * mean size)
Dens$Biomass.kg<-Dens$YE_abund*Dens$mean.YE.kg_rand

#Calculate the variance of biomass
Dens$var.Biomass.kg<-(Dens$YE_abund^2)*Dens$bootvar.YE.kg_rand+(Dens$mean.YE.kg_rand^2)*Dens$var.Abund+
  Dens$var.Abund*Dens$bootvar.YE.kg_rand

#Calculate coefficient of variation for Biomass
Dens$Biomass.cv<-sqrt(Dens$var.Biomass.kg)/Dens$Biomass.kg

#Cinvert biomass (kg) to metric tons
Dens$Biomass.mt<-Dens$Biomass.kg/1000

eg <-Dens[Dens$Subdistrict =="EYKT",]; head(eg); eg$mean.YE.kg


str(Dens)

ggplot(Dens, aes(x=Year)) +
  geom_point(aes(y=Biomass.mt),size=2) +
  geom_errorbar(aes(ymin = Biomass.mt-1.95*Biomass.cv*Biomass.mt, 
                    ymax= Biomass.mt+1.95*Biomass.cv*Biomass.mt),
                col="black", alpha=0.5) +
  facet_wrap(~Subdistrict) +
  xlab("\nYear") +
  ylab("yelloweye rockfish biomass (mt)") +
  #ylab("Density (yelloweye rockfish/kmsq)") +
  scale_y_continuous(label=comma) +
  scale_x_continuous(breaks=seq(1995,2025,5)) + 
  theme (axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
         panel.grid.minor = element_blank())
ggsave(paste0("Figures/YE_Biomass_MA_", YEAR, ".png"), dpi=300,  height=4, width=4, units="in")

#plot biomass and density in same plot...
# work in progress... difficulty with secondary axis that is not worth the time right now... 
ggplot(Dens, aes(x=Year)) +
  geom_point(aes(y=Biomass.mt),size=2, col="red") +
  geom_point(aes(y=Density),size=2, col="blue") +
  geom_errorbar(aes(ymin = Biomass.mt-1.95*Biomass.cv*Biomass.mt, 
                    ymax= Biomass.mt+1.95*Biomass.cv*Biomass.mt),
                col="red", alpha=0.5) +
  geom_errorbar(aes(ymin = Density-1.95*sqrt(var.Dens), ymax= Density+1.95*sqrt(var.Dens)),
                col="blue", alpha=0.5) +
  facet_wrap(~Subdistrict) +
  xlab("\nYear") +
  ylab("Yelloweye rockfish biomass (mt)") +
  #ylab("Density (yelloweye rockfish/kmsq)") +
  scale_y_continuous(label=comma,
                     sec.axis=sec_axis(~.*0.25, name="Density (fish/km^2)", label=comma)) +
  scale_x_continuous(breaks=seq(1995,2025,5)) + 
  theme (axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
         axis.text.y = element_text(color = "black"),
         axis.text.y.right = element_text(color = "black"),
         panel.grid.minor = element_blank())

ggsave(paste0("Figures/YE_Density_&_Biomass_MA_", YEAR, ".png"), dpi=300,  height=4, width=4, units="in")

#************************************************************************************
#*
#*
# Get "Status-quo" biomass estimates and lower 90% CI

Years2<-unique(Dens$Year)
Biomass.sq<-data.frame()
k<-1
for (i in Years2){  #i<-Years[2]
  DAT<-Dens[Dens$Year == i,]
  if (i == 1994) {} else {
    bio<-vector()
    bvar<-vector()
    jj<-1
    for (j in GFMA) {  #j<-GFMA[1]
      
      Dsub<-DAT[DAT$Subdistrict == j,]
      if(is.na(Dsub$Biomass.mt)) {
        past<-Dens[Dens$Year < i & Dens$Subdistrict == j,]
        latN<-past$YE_abund
        Ns<-latN[!is.na(latN)]
        useN<-Ns[length(Ns)]
        latV<-past$var.Abund
        Vs<-latV[!is.na(latV)]
        useV<-Vs[length(Vs)]
        
        bio[jj]<-useN*Dsub$mean.YE.kg_rand
        bvar[jj]<-(bio[jj]^2)*Dsub$bootvar.YE.kg_rand+(Dsub$mean.YE.kg_rand^2)*useV+
          useV*Dsub$bootvar.YE.kg_rand
        jj<-jj+1
      } else {
        bio[jj]<-Dsub$Biomass.kg
        bvar[jj]<-Dsub$var.Biomass.kg
        jj<-jj+1
      }
    }
    Biomass.sq[k,"Year"]<-i
    Biomass.sq[k,"Biomass"]<-sum(bio)
    Biomass.sq[k,"Biom.var_noVov"]<-sum(bvar)
    #Biomass.sq[k,"Biomass.cv"]<-Biomass.sq[k,"Biomass.var"]*Biomass.sq[k,"Biomass.var"]/
    #  Biomass.sq[k,"Biomass"]
    Biomass.sq[k,"Biomass.CS"]<-bio[1]
    Biomass.sq[k,"Biomass.EY"]<-bio[2]
    Biomass.sq[k,"Biomass.NS"]<-bio[3]
    Biomass.sq[k,"Biomass.SS"]<-bio[4]
    Biomass.sq[k,"var.BiomassCS"]<-bvar[1]
    Biomass.sq[k,"var.BiomassEY"]<-bvar[2]
    Biomass.sq[k,"var.BiomassNS"]<-bvar[3]
    Biomass.sq[k,"var.BiomassSS"]<-bvar[4]
    k<-k+1
  }
}

Biomass.cov<-2*cov(x=Biomass.sq$Biomass.CS,y=Biomass.sq$Biomass.EY)+
  2*cov(Biomass.sq$Biomass.CS,Biomass.sq$Biomass.NS)+
  2*cov(Biomass.sq$Biomass.CS,Biomass.sq$Biomass.SS)+
  2*cov(Biomass.sq$Biomass.EY,Biomass.sq$Biomass.NS)+
  2*cov(Biomass.sq$Biomass.EY,Biomass.sq$Biomass.SS)+
  2*cov(Biomass.sq$Biomass.NS,Biomass.sq$Biomass.SS)

Biomass.sq$var.Biomass<-Biomass.sq$Biom.var_noVov+Biomass.cov
  
Biomass.sq$Biomass.cv<-sqrt(Biomass.sq$var.Biomass)/Biomass.sq$Biomass
Biomass.sq$Biomass_mt<-Biomass.sq$Biomass/1000
Biomass.sq$Biomass_mt_lo90<-Biomass.sq$Biomass_mt-Biomass.sq$Biomass_mt*1.68*Biomass.sq$Biomass.cv
Biomass.sq$Biomass_mt_hi90<-Biomass.sq$Biomass_mt+Biomass.sq$Biomass_mt*1.68*Biomass.sq$Biomass.cv

#write.csv(Biomass.sq,"Data_processing/Data/SEO_YE_Biomass_08292024.csv")

#plot it for SAFE report
head(Biomass.sq)
Biomass.sq %>% mutate(Bio.CSEO.mt = Biomass.CS/1000,
                      Bio.EYKT.mt = Biomass.EY/1000,
                      Bio.NSEO.mt = Biomass.NS/1000,
                      Bio.SSEO.mt = Biomass.SS/1000,
                      Bio.CSEO.cv = sqrt(var.BiomassCS)/Biomass.CS,
                      Bio.EYKT.cv = sqrt(var.BiomassEY)/Biomass.EY,
                      Bio.NSEO.cv = sqrt(var.BiomassNS)/Biomass.NS,
                      Bio.SSEO.cv = sqrt(var.BiomassSS)/Biomass.SS,
                      Bio.CSEO.lo90 = Bio.CSEO.mt-1.68*Bio.CSEO.cv*Bio.CSEO.mt,
                      Bio.EYKT.lo90 = Bio.EYKT.mt-1.68*Bio.EYKT.cv*Bio.EYKT.mt,
                      Bio.NSEO.lo90 = Bio.NSEO.mt-1.68*Bio.NSEO.cv*Bio.NSEO.mt,
                      Bio.SSEO.lo90 = Bio.SSEO.mt-1.68*Bio.SSEO.cv*Bio.SSEO.mt,) -> Biomass.sq

#write.csv(Biomass.sq,"Data_processing/Data/SEO_YE_Biomass_subd_08292024.csv")
#===================================================================================
xl<-expression()

plot_ly(y=~c(Biomass.sq$Biomass_mt), x=~c(Biomass.sq$Year), type="line", 
     ylim=c(5000,max(Biomass.sq$Biomass_mt_hi90)),
     ylab=c("SEO yelloweye biomass (t)"),
     xlab="Year")
polygon(x=c(Biomass.sq$Year,rev(Biomass.sq$Year)),
        y=c(Biomass.sq$Biomass_mt_lo90,rev(Biomass.sq$Biomass_mt_hi90)),
        col="lightblue",border=F)
lines(Biomass.sq$Biomass_mt ~ Biomass.sq$Year, type="l",lwd=2, pch=18)
points(Biomass.sq$Biomass_mt ~ Biomass.sq$Year, type="p",cex=1.5, pch=19)

CurrentDate<-Sys.Date()
write.csv(Dens,paste0("Data_processing/Data/SEO_YE_Biomass_subdistrict_",CurrentDate,".csv"))
  
  
ggplot(Biomass.sq, aes(x=Year)) +
  geom_line(aes(y=Biomass_mt),linewidth=1) +
  geom_ribbon(aes(ymin = Biomass_mt_lo90, ymax= Biomass_mt_hi90),col="lightgrey", alpha=0.5) +
  xlab("\nYear") +
  ylab("SEO yelloweye biomass (t)\n") +
  scale_y_continuous(label=comma, breaks = c(10000,15000,20000,25000,30000,35000,40000)) +
  scale_x_continuous(breaks=seq(1995,2025,5)) + 
  theme (panel.grid.major = element_blank (),panel.grid.minor = element_blank())
  
YEAR<-2023
ggsave(paste0("Figures/SEO_status_quo_biomass_", YEAR, ".png"), dpi=300,  height=4, width=5, units="in")

##Compare to original Status-quo from the archives measurements

SQ_hist<-read.csv("Data/YE_BiomassEstimates_2022_STATUS_QUO_ORIG.csv")
str(SQ_hist); unique(SQ_hist$BiomassType)

ys<-unique(SQ_hist$Year)
SQ<-data.frame(); i<-1
for (y in ys){
  dat<-SQ_hist[SQ_hist$Year==y,]
  SQ[i,"Year"]<-y-1
  SQ[i,"Biomass_mt"]<-dat$Biomass_Estimate[dat$BiomassType=="Biomass_mt"]
  SQ[i,"Biomass_mt_lo90"]<-dat$Biomass_Estimate[dat$BiomassType=="LCI_90"]
  
  if (nrow(dat[dat$BiomassType == "UCI_90",]) == 0){
    SQ[i,"Biomass_mt_hi90"]<-dat$Biomass_Estimate[dat$BiomassType=="Biomass_mt"]+
      (dat$Biomass_Estimate[dat$BiomassType=="Biomass_mt"]-
         dat$Biomass_Estimate[dat$BiomassType=="LCI_90"])
  } else {
    SQ[i,"Biomass_mt_hi90"]<-dat$Biomass_Estimate[dat$BiomassType=="UCI_90"]
  }
  
  i<-i+1
}

oldnew<-rbind(Biomass.sq %>% mutate(Method = "revised status quo") %>% 
                dplyr::select(Year,Biomass_mt,Biomass_mt_lo90,Biomass_mt_hi90,Method=Method),
              SQ %>% mutate(Method = "original status quo") %>% 
                dplyr::select(Year,Biomass_mt,Biomass_mt_lo90,Biomass_mt_hi90,Method=Method))

ggplot(oldnew, aes(x=Year)) +
  geom_ribbon(aes(ymin = Biomass_mt_lo90, ymax= Biomass_mt_hi90, fill=Method), alpha=0.35) +
  geom_line(aes(y=Biomass_mt,col=Method),size=1) +
  xlab("\nYear") +
  ylab("SEO yelloweye biomass (t)\n") +
  scale_y_continuous(label=comma, breaks = c(10000,15000,20000,25000,30000,35000,40000,45000,50000)) +
  scale_x_continuous(breaks=seq(1995,2025,5)) + 
  theme (panel.grid.minor = element_blank (), #panel.grid.major = element_blank(),
         legend.position=c(0.8,0.8))

ggsave(paste0("Figures/SEO_status_quo_biomass_comp2_", YEAR, ".png"), dpi=300,  height=4, width=5, units="in")


