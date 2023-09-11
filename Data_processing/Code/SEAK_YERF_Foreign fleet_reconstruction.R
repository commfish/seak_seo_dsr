###############################################################################
## Estimation of yelloweye removals from the foreign fleet in southeast AK
##
## Developed by Donnie Arthur with ADF&G Sport Fish Division
## See the documents folder for background an initial write-up of methods
## Donnie's original, unmodified work can be found in the SEO_DSR/Foreign removals Donnie's raw staff folder
##
##
## Adopted here and modified by Phil Joy 9-11-23
##
## Background (Donnie)
# A Pacific ocean perch trawl fishery in the Gulf of Alaska developed in the early 
# 1960’s with large effort by the U.S.S.R and Japanese fleets. At the height of 
# the fishery in 1965, the catches of all rockfish, including Pacific ocean perch, 
# exceeded 370,000 metric tons (mt). Generally, catches declined following this 
# peak until foreign fishing was abolished from the Gulf of Alaska in 1987. During 
# the early period of this foreign fishery (1961-1974), foreign catches of rockfish 
# were often reported in crude management groups, including “Pacific ocean perch” 
# or “other rockfish”, with no differentiating between species. With implementation 
# of a fishery observer program in 1975 and 1977 onward, species composition, 
# including Pacific ocean perch and Yelloweye Rockfish, of foreign catches became available. 
#
# From 1961-1974, catches by the foreign fleet were reported only as Gulf-wide; 
# however, the spatial resolution of foreign catch reporting improved in the second 
# half of the fishery. In the later years of the Gulf of Alaska foreign fishery, 
# catches were reported to International North Pacific Fisheries Commission (INPFC) 
# areas. The Southeastern INPFCs area nearly aligns with the Southeast Outside (SEO) 
# subdistrict of the Gulf of Alaska, and for the purpose of this catch reconstruction, 
# the two areas will be treated as the same. The foreign fleet was banned from 
# the Southeastern INPFC area earlier than any of the other INPFC areas in 1982. 
# Catches were reconstructed for the Southeastern INPFC area, and by proxy, for 
# the SEO subdistrict.
#
# Methods (Donnie)
#Non-parametric bootstrapping was used to reconstruct all non-P.O.P. rockfish 
# foreign catch for 1961-1972 by resampling the percentage of all non-P.O.P rockfish 
# catches from all rockfish including P.O.P. catches from 1973-1984 and multiplied 
# by the reported all rockfish catches for each given year from 1961-1972. Next, 
# the Gulf-wide yelloweye rockfish catches were reconstructed for 1961-1974 and 
# 1976 by bootstrap sampling the percent of yelloweye Rockfish catches from the 
# all non-P.O.P. rockfish catches from 1977-1985, which was then multiplied by 
# the reconstructed all non-P.O.P. catch estimates for 1961-1972 and the reported 
# all non-P.O.P catches for 1973-1974, 1976.  Lastly, the Gulf-wide catches were 
# allocated to the Southeastern INPFC area by bootstrap sampling the percentage 
# of foreign Yelloweye rockfish catches from the Southeastern INPFC area of the 
# Gulf-wide yelloweye rockfish catch totals for 1978-1981. The mean of the bootstrap 
# samples was used as estimate of Southeastern yelloweye rockfish foreign catches 
# and the percentile method was used to construct 95% bootstrap confidence intervals. 
# A coefficient of variation of 0.75 was applied to bootstrap estimates to capture 
# the high degree of uncertainty in catch reconstructions. 
#
# Methods as published in 2022 SAFE report, written by Phil Joy
# Removals of yelloweye rockfish in the foreign fleet were estimated in the years 
# prior to accurate record keeping by applying proportional relationships derived 
# from years when data was available (Table 14.6; Berger et al. 1984, 1985, 1987; 
# Berger and Weikart 1988; Forrester et al. 1978, 1983; Nelson et al. 1983; Wall 
# et al. 1978, 1979, 1980, 1981, 1982).  This reconstruction remains a work in 
# progress.  First, non-POP harvests were estimated for 1961-1972 by applying the 
# proportional relationship between non-POP and all rockfish in 1973-1984 data 
# (Table 14.6).  Secondly, yelloweye rockfish harvests were estimated for 1961-1974 
# by applying the proportional relationship between non-POP rockfish and known 
# yelloweye rockfish in 1977-1985.  Lastly, the Southeastern INPFC harvests for 
## years prior to 1978 were estimated by applying the proportional relationship 
# between Southeastern INPFC yelloweye rockfish harvests and gulf-wide yelloweye 
# rockfish harvests in 1978-1981.  During each step, non-parametric bootstrapping 
# was used to estimate variance and cv’s and the uncertainty passed forward.  A 
# coefficient of variation of 0.75 was applied to final estimates to capture the 
# high degree of uncertainty in catch reconstructions.  
################################################################################
#### Bootstrap function ####
library(purrr)
library(ggplot2)
library(dplyr)
library(ggpubr)

bootstrap <- function(data, n) {
  resampled_data <- lapply(1:n, function(i) {
    resample <- sample(data, replace = TRUE)
  })
  return(resampled_data)
}

################################################################################
## After reviewing Donnie's code, which is preserved below,  I think the methods are
## more-or-less sound, but I am worried that something might be getting lost in the
## sequential bootstrapping approach that he used.  I am going to use the same steps
## as Donnie, but the bootstrapping procedure will cover the whole estimation
## process and all the associated steps.  
################################################################################
harvest <- read.csv("Data_processing/Data/YERF_Foreign_reconstruction.csv")
str(harvest)

# 1) All non-POP rockfish harvested between '73 and '84
RF_excPOP_GOA_1973_84 <- harvest$RF_excPOP_harvest_GOA[14:25]

# 2) All rockfish harvested between '73 and '84, INCLUDING POP
RF_GOA_1973_84 <- harvest$RF_harvest_GOA[14:25]

# 3) Rockfish harvest in the GOA for entire time series
RF_GOA <- harvest$RF_harvest_GOA

# 4) proportional relationship between non-POP and POP in harvest between '73 and '84
prop_RFexcPOP_RF <- data.frame(RF_excPOP_GOA_1973_84/RF_GOA_1973_84)
mean(prop_RFexcPOP_RF[,1])

bootstrapped_samples <- bootstrap(prop_RFexcPOP_RF$
                                    RF_excPOP_GOA_1973_84.RF_GOA_1973_84,1000)
bootstrapped_samples_means=lapply(bootstrapped_samples,mean)

str(bootstrapped_samples_means)
boot_means<-unlist(bootstrapped_samples_means)

     # boot mean versus mean to examine bias: 
mean(boot_means)
mean(prop_RFexcPOP_RF[,1])

# 5)  non-POP rockfish harvest during observed years
RF_excPOP_GOA_1977_85 <- harvest$RF_excPOP_harvest_GOA[18:26]

# 6) Yelloweye harvests during observed years
YERF_GOA_1977_85 <- harvest$YERF_harvest_GOA[18:26]

# 7) nonPOP rockfish harvest in all years
RF_excPOP_GOA <- harvest$RF_excPOP_harvest_GOA

# 8) proportional relationship between YE harvest and nonPOP harvest during observed years:
prop_YERF_RF <- data.frame(YERF_GOA_1977_85/RF_excPOP_GOA_1977_85)
mean(prop_YERF_RF[,1])

# 9) YE harvested in SEAK
YERF_SEAK_1978_81 <- harvest$YERF_Southeastern[19:22]

# 10) YE harvested in the entire GOA
YERF_GOA_1978_81 <- harvest$YERF_harvest_GOA[19:22] 

# 11) Proportional relationship between SEAK YERF and GOA YERF
prop_SE_YERF <- data.frame(YERF_SEAK_1978_81/YERF_GOA_1978_81)
mean(prop_SE_YERF[,1])

# 12 Estimation loop

years <- seq(1960,1987,1)

for_rem<-data.frame(matrix(nrow=length(years)))

for_rem[,"Year"]<-years
for_rem[,"mean_prop_nonPOP_to_allRF"] <-mean(prop_RFexcPOP_RF[,1])
for_rem[,"nonPOP_RF_harv"]<-mean(prop_RFexcPOP_RF[,1])*RF_GOA
for_rem[,"mean_prop_YERF_to_nonPOP"] <- mean(prop_YERF_RF[,1])
for_rem[,"YERF_GOA_harv"]<- mean(prop_YERF_RF[,1])*for_rem[,"nonPOP_RF_harv"]
for_rem[,"mean_prop_SEAK_to_GOA"]<-mean(prop_SE_YERF[,1])
for_rem[,"YERF_SEAK_harv_full"]<-mean(prop_SE_YERF[,1])*for_rem[,"YERF_GOA_harv"]
for_rem[c(1:13),"YERF_SEAK_harv"]<-for_rem[c(1:13),"YERF_SEAK_harv_full"]
for_rem[c(14,15,17),"YERF_SEAK_harv"]<-harvest$RF_excPOP_harvest_GOA[c(14,15,17)]*
  mean(prop_YERF_RF[,1])*mean(prop_SE_YERF[,1])
for_rem[c(16,c(18:28)),"YERF_SEAK_harv"]<-harvest$YERF_harvest_GOA[c(16,c(18:28))]*mean(prop_SE_YERF[,1])

#to get confidence intervals and cv's we will bootstrap... 
nboot<-10000

i<-1
  boot_SEAK_YE<-data.frame(matrix(nrow=length(years)))
  boot_SEAK_YE[,i]<-years
  
  boot_SEAK_YE_full<-data.frame(matrix(nrow=length(years)))
  boot_SEAK_YE_full[,i]<-years
  
  for (i in 1:nboot){ #i<-1
    #boot sample proportion nonPOP:allRF
    nonPOP_allRF<-sample(prop_RFexcPOP_RF[,1],length(prop_RFexcPOP_RF[,1]),replace=T)
    YERF_to_nonPOP<-sample(prop_YERF_RF[,1],length(prop_YERF_RF[,1]),replace=T)
    SEAK_to_GOA<-sample(prop_SE_YERF[,1],length(prop_SE_YERF[,1]),replace=T)
    
  #FLAG!! I am really diverging from Donnie's approach and using the best data for each year
    # i.e, for 1975 we have estimates of YERF in GOA so will use those numbers
    for (y in 1:length(years)){  #y<-1``
      if (y < 14) {  
        boot_SEAK_YE[y,i+1]<-RF_GOA[y]*mean(nonPOP_allRF)*mean(YERF_to_nonPOP)*mean(SEAK_to_GOA)
      }
      if (y %in% c(14,15,17)) {
        boot_SEAK_YE[y,i+1]<-harvest$RF_excPOP_harvest_GOA[y]*mean(YERF_to_nonPOP)*mean(SEAK_to_GOA)
      } 
      if (y %in% c(16,c(18:28))) {
        boot_SEAK_YE[y,i+1]<-harvest$YERF_harvest_GOA[y]*mean(SEAK_to_GOA)
      }
    }
    
    boot_SEAK_YE_full[,i+1]<-RF_GOA*mean(nonPOP_allRF)*mean(YERF_to_nonPOP)*mean(SEAK_to_GOA)
    i<-i+1
  }
  
for (y in 1:length(years)) {  #y<-1
  v<-as.numeric(as.vector(boot_SEAK_YE[y,c(2:(nboot+1))]))
  for_rem[y,"lo_95"]<-quantile(v,c(0.025))
  for_rem[y,"lo_75"]<-quantile(v,c(0.125))
  for_rem[y,"lo_50"]<-quantile(v,c(0.25))
  for_rem[y,"hi_95"]<-quantile(v,c(0.975))
  for_rem[y,"hi_75"]<-quantile(v,c(0.875))
  for_rem[y,"hi_50"]<-quantile(v,c(0.75))
  for_rem[y,"var"]<-var(v)
  for_rem[y,"cv"]<-sd(v)/for_rem[y,"YERF_SEAK_harv"]
  
  v2<-as.numeric(as.vector(boot_SEAK_YE_full[y,c(2:(nboot+1))]))
  for_rem[y,"lo_95_full"]<-quantile(v2,c(0.025))
  for_rem[y,"lo_75_full"]<-quantile(v2,c(0.125))
  for_rem[y,"lo_50_full"]<-quantile(v2,c(0.25))
  for_rem[y,"hi_95_full"]<-quantile(v2,c(0.975))
  for_rem[y,"hi_75_full"]<-quantile(v2,c(0.875))
  for_rem[y,"hi_50_full"]<-quantile(v2,c(0.75))
  for_rem[y,"var_full"]<-var(v2)
  for_rem[y,"cv_full"]<-sd(v2)/for_rem[y,"YERF_SEAK_harv"]
}
   
ggplot(for_rem) + 
  geom_ribbon(aes(Year, ymin=lo_95,ymax=hi_95),fill="blue",alpha=0.2, col=NA) + 
  geom_ribbon(aes(Year, ymin=lo_75,ymax=hi_75),fill="blue",alpha=0.2, col=NA) +
  geom_ribbon(aes(Year, ymin=lo_50,ymax=hi_50),fill="blue",alpha=0.2, col=NA) + 
  geom_point(aes(Year,YERF_SEAK_harv), col="blue") +
  geom_line(aes(Year,YERF_SEAK_harv), col="blue") +
  annotate("text", y=1300, # 
           x=1980, color="blue",
           label=c(paste0("mean annual removals = ",round(mean(for_rem$YERF_SEAK_harv),0)," (t)"))) +
  #geom_point(data=harvest,aes(Year,YERF_harvest_GOA*mean(prop_SE_YERF[,1])), col="red") +
  labs(x = "Year", y = "Yelloweye rockfish removals by the foreign trawl fleet (t)",
       title = "Calculated from best available data by year") +
  scale_x_continuous(breaks=seq(1960,1990,5)) + 
  theme (axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
         panel.grid.minor = element_blank()) +
  ylim(0,1500) -> for_removals_A

for_removals_A
ggsave(paste0("Figures/foreigh_trawl_removals_1960-1987_best_dat_byyear.png"), dpi=300, height=4, width=7, units="in")

ggplot(for_rem) + 
  #geom_ribbon(aes(Year, ymin=lo_95,ymax=hi_95),fill="blue",alpha=0.2, col=NA) + 
  #geom_ribbon(aes(Year, ymin=lo_75,ymax=hi_75),fill="blue",alpha=0.2, col=NA) +
  #geom_ribbon(aes(Year, ymin=lo_50,ymax=hi_50),fill="blue",alpha=0.2, col=NA) + 
  geom_ribbon(aes(Year, ymin=lo_95_full,ymax=hi_95_full),fill="aquamarine4",alpha=0.2, col=NA) + 
  geom_ribbon(aes(Year, ymin=lo_75_full,ymax=hi_75_full),fill="aquamarine4",alpha=0.2, col=NA) +
  geom_ribbon(aes(Year, ymin=lo_50_full,ymax=hi_50_full),fill="aquamarine4",alpha=0.2, col=NA) + 
  #geom_point(aes(Year,YERF_SEAK_harv), col="blue") +
  #geom_line(aes(Year,YERF_SEAK_harv), col="blue") +
  geom_point(aes(Year,YERF_SEAK_harv_full), col="aquamarine4") +
  geom_line(aes(Year,YERF_SEAK_harv_full), col="aquamarine4") +
  annotate("text", y=1300, # 
           x=1980, color="aquamarine4",
           label=c(paste0("mean annual removals = ",round(mean(for_rem$YERF_SEAK_harv_full),0)," (t)"))) +
  #geom_point(data=harvest,aes(Year,YERF_harvest_GOA*mean(prop_SE_YERF[,1])), col="red") +
  labs(x = "Year", y = "Yelloweye rockfish removals by the foreign trawl fleet (t)",
       title = "Calculated from all rockfish estimates in GOA") +
  scale_x_continuous(breaks=seq(1960,1990,5)) + 
  theme (axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
         panel.grid.minor = element_blank()) +
  ylim(0,1500) -> for_removals_B

for_removals_B
ggsave(paste0("Figures/foreigh_trawl_removals_1960-1987_da_mthds.png"), dpi=300, height=4, width=7, units="in")

ggarrange(for_removals_A, for_removals_B,
          nrow=2)

ggsave(paste0("Figures/foreigh_trawl_removals_1960-1987.png"), dpi=300, height=8, width=7, units="in")

for_rem_SPM<-for_rem %>% select(Year, Estimate = YERF_SEAK_harv_full, cv=cv_full)

write.csv(for_rem,"Data_processing/Data/Harvests/Foreign_YERF_SEAK_ests.csv")

write.csv(for_rem_SPM,"Data_processing/Data/Harvests/Foreign_YERF_SEAK.csv")
#write.csv(for_rem, "C:\\Users\\dearthur\\Documents\\YERF_SEAK_foreign.csv",row.names=FALSE)


################################################################################
##### Donnie's original code and methods below: 
################################################################################
#### Proportion of Non-POP RF of All RF harvest Bootstrapping ####
harvest <- read.csv("Data_processing/Data/YERF_Foreign_reconstruction.csv")
str(harvest)

# 1) All non-POP rockfish harvested between '73 and '84
RF_excPOP_GOA_1973_84 <- harvest$RF_excPOP_harvest_GOA[14:25]

# 2) All rockfish harvested between '73 and '84, INCLUDING POP
RF_GOA_1973_84 <- harvest$RF_harvest_GOA[14:25]

# 3) Rockfish harvest in the GOA for entire time series
RF_GOA <- harvest$RF_harvest_GOA

# 4) proportional relationship between non-POP and POP in harvest between '73 and '84
prop_RFexcPOP_RF <- data.frame(RF_excPOP_GOA_1973_84/RF_GOA_1973_84)
mean(prop_RFexcPOP_RF[,1])

# 5) bootstrap resamples of those proportions estimated in step 4.
bootstrapped_samples <- bootstrap(prop_RFexcPOP_RF$
    RF_excPOP_GOA_1973_84.RF_GOA_1973_84,1000)
bootstrapped_samples_means=lapply(bootstrapped_samples,mean)

str(bootstrapped_samples_means)
boot_means<-unlist(bootstrapped_samples_means)

# 6) boot mean versus mean to examine bias: 
mean(boot_means)
mean(prop_RFexcPOP_RF[,1])
# 6) mean and median of 
###################################################################

#### Non-POP RF harvest reconstructed from All RF (incl. POP) harvest #####
### Bootstrapped means muliplied by reported All RF harvest###
## Proportion non-POP applied to all rockfish = est. of non-POP harvest...

{RF_GOA_reconst_1960 <- Map('*',bootstrapped_samples_means,RF_GOA[1])
RF_GOA_reconst_1961 <- Map('*',bootstrapped_samples_means,RF_GOA[2])
RF_GOA_reconst_1962 <- Map('*',bootstrapped_samples_means,RF_GOA[3])
RF_GOA_reconst_1963 <- Map('*',bootstrapped_samples_means,RF_GOA[4])
RF_GOA_reconst_1964 <- Map('*',bootstrapped_samples_means,RF_GOA[5])
RF_GOA_reconst_1965 <- Map('*',bootstrapped_samples_means,RF_GOA[6])
RF_GOA_reconst_1966 <- Map('*',bootstrapped_samples_means,RF_GOA[7])
RF_GOA_reconst_1967 <- Map('*',bootstrapped_samples_means,RF_GOA[8])
RF_GOA_reconst_1968 <- Map('*',bootstrapped_samples_means,RF_GOA[9])
RF_GOA_reconst_1969 <- Map('*',bootstrapped_samples_means,RF_GOA[10])
RF_GOA_reconst_1970 <- Map('*',bootstrapped_samples_means,RF_GOA[11])
RF_GOA_reconst_1971 <- Map('*',bootstrapped_samples_means,RF_GOA[12])
RF_GOA_reconst_1972 <- Map('*',bootstrapped_samples_means,RF_GOA[13])
RF_GOA_reconst_1973 <- Map('*',bootstrapped_samples_means,RF_GOA[14])
RF_GOA_reconst_1974 <- Map('*',bootstrapped_samples_means,RF_GOA[15])
RF_GOA_reconst_1975 <- Map('*',bootstrapped_samples_means,RF_GOA[16])
RF_GOA_reconst_1976 <- Map('*',bootstrapped_samples_means,RF_GOA[17])
RF_GOA_reconst_1977 <- Map('*',bootstrapped_samples_means,RF_GOA[18])
RF_GOA_reconst_1978 <- Map('*',bootstrapped_samples_means,RF_GOA[19])
RF_GOA_reconst_1979 <- Map('*',bootstrapped_samples_means,RF_GOA[20])
RF_GOA_reconst_1980 <- Map('*',bootstrapped_samples_means,RF_GOA[21])
RF_GOA_reconst_1981 <- Map('*',bootstrapped_samples_means,RF_GOA[22])
RF_GOA_reconst_1982 <- Map('*',bootstrapped_samples_means,RF_GOA[23])
RF_GOA_reconst_1983 <- Map('*',bootstrapped_samples_means,RF_GOA[24])
RF_GOA_reconst_1984 <- Map('*',bootstrapped_samples_means,RF_GOA[25])
RF_GOA_reconst_1985 <- Map('*',bootstrapped_samples_means,RF_GOA[26])
RF_GOA_reconst_1986 <- Map('*',bootstrapped_samples_means,RF_GOA[27])
RF_GOA_reconst_1987 <- Map('*',bootstrapped_samples_means,RF_GOA[28])}

##########################################################################

#### Proportion of gulf-wide YERF of Non-POP RF harvest Bootstrapping ####
### Bootstrap resampling of prop. YERF of Non-POP RF harvest 1977-1985###

## OK, Donnie has the bootstrap samples of non_POP rockfish harvests above
## Here he is applying the proportion of YERF:non-POP from the observed years to 
## those estimates

## YE harvest = prop_YE*non_POP estimate

## 7)  non-POP rockfish harvest during observed years
RF_excPOP_GOA_1977_85 <- harvest$RF_excPOP_harvest_GOA[18:26]

## 8) Yelloweye harvests during observed years
YERF_GOA_1977_85 <- harvest$YERF_harvest_GOA[18:26]

## 9) nonPOP rockfish harvest in all years
RF_excPOP_GOA <- harvest$RF_excPOP_harvest_GOA

# 10) proportional relationship between YE harvest and nonPOP harvest during observed years:
prop_YERF_RF <- data.frame(YERF_GOA_1977_85/RF_excPOP_GOA_1977_85)
mean(prop_YERF_RF[,1])

# 11) bootstrap those proportions estimated in step 10
YEbootstrapped_samples <- bootstrap(prop_YERF_RF$
                                    YERF_GOA_1977_85.RF_excPOP_GOA_1977_85,1000)

boot_means_YE<-unlist(YEbootstrapped_samples)

# 12) boot mean versus mean to examine bias: 
mean(boot_means_YE)
mean(prop_YERF_RF[,1])

###Sampled proportions expanded to reconstructed NON-POP RF harvest###
# Apply the proportions from step 10-11 to get estimates of YE removals 
{YERF_GOA_reconst_1960 <- Map('*',YEbootstrapped_samples,RF_GOA_reconst_1960)
YERF_GOA_reconst_1961 <- Map('*',YEbootstrapped_samples,RF_GOA_reconst_1961)
YERF_GOA_reconst_1962 <- Map('*',YEbootstrapped_samples,RF_GOA_reconst_1962)
YERF_GOA_reconst_1963 <- Map('*',YEbootstrapped_samples,RF_GOA_reconst_1963)
YERF_GOA_reconst_1964 <- Map('*',YEbootstrapped_samples,RF_GOA_reconst_1964)
YERF_GOA_reconst_1965 <- Map('*',YEbootstrapped_samples,RF_GOA_reconst_1965)
YERF_GOA_reconst_1966 <- Map('*',YEbootstrapped_samples,RF_GOA_reconst_1966)
YERF_GOA_reconst_1967 <- Map('*',YEbootstrapped_samples,RF_GOA_reconst_1967)
YERF_GOA_reconst_1968 <- Map('*',YEbootstrapped_samples,RF_GOA_reconst_1968)
YERF_GOA_reconst_1969 <- Map('*',YEbootstrapped_samples,RF_GOA_reconst_1969)
YERF_GOA_reconst_1970 <- Map('*',YEbootstrapped_samples,RF_GOA_reconst_1970)
YERF_GOA_reconst_1971 <- Map('*',YEbootstrapped_samples,RF_GOA_reconst_1971)
YERF_GOA_reconst_1972 <- Map('*',YEbootstrapped_samples,RF_GOA_reconst_1972)
YERF_GOA_reconst_1973 <- Map('*',YEbootstrapped_samples,RF_GOA_reconst_1973)
YERF_GOA_reconst_1974 <- Map('*',YEbootstrapped_samples,RF_GOA_reconst_1974)
YERF_GOA_reconst_1975 <- Map('*',YEbootstrapped_samples,RF_GOA_reconst_1975)
YERF_GOA_reconst_1976 <- Map('*',YEbootstrapped_samples,RF_GOA_reconst_1976)
YERF_GOA_reconst_1977 <- Map('*',YEbootstrapped_samples,RF_GOA_reconst_1977)
YERF_GOA_reconst_1978 <- Map('*',YEbootstrapped_samples,RF_GOA_reconst_1978)
YERF_GOA_reconst_1979 <- Map('*',YEbootstrapped_samples,RF_GOA_reconst_1979)
YERF_GOA_reconst_1980 <- Map('*',YEbootstrapped_samples,RF_GOA_reconst_1980)
YERF_GOA_reconst_1981 <- Map('*',YEbootstrapped_samples,RF_GOA_reconst_1981)
YERF_GOA_reconst_1982 <- Map('*',YEbootstrapped_samples,RF_GOA_reconst_1982)
YERF_GOA_reconst_1983 <- Map('*',YEbootstrapped_samples,RF_GOA_reconst_1983)
YERF_GOA_reconst_1984 <- Map('*',YEbootstrapped_samples,RF_GOA_reconst_1984)
YERF_GOA_reconst_1985 <- Map('*',YEbootstrapped_samples,RF_GOA_reconst_1985)
YERF_GOA_reconst_1986 <- Map('*',YEbootstrapped_samples,RF_GOA_reconst_1986)
YERF_GOA_reconst_1987 <- Map('*',YEbootstrapped_samples,RF_GOA_reconst_1987)}

{YERF_GOA_reconst_1960_means <- lapply(YERF_GOA_reconst_1960,mean)
YERF_GOA_reconst_1961_means <- lapply(YERF_GOA_reconst_1961,mean)
YERF_GOA_reconst_1962_means <- lapply(YERF_GOA_reconst_1962,mean)
YERF_GOA_reconst_1963_means <- lapply(YERF_GOA_reconst_1963,mean)
YERF_GOA_reconst_1964_means <- lapply(YERF_GOA_reconst_1964,mean)
YERF_GOA_reconst_1965_means <- lapply(YERF_GOA_reconst_1965,mean)
YERF_GOA_reconst_1966_means <- lapply(YERF_GOA_reconst_1966,mean)
YERF_GOA_reconst_1967_means <- lapply(YERF_GOA_reconst_1967,mean)
YERF_GOA_reconst_1968_means <- lapply(YERF_GOA_reconst_1968,mean)
YERF_GOA_reconst_1969_means <- lapply(YERF_GOA_reconst_1969,mean)
YERF_GOA_reconst_1970_means <- lapply(YERF_GOA_reconst_1970,mean)
YERF_GOA_reconst_1971_means <- lapply(YERF_GOA_reconst_1971,mean)
YERF_GOA_reconst_1972_means <- lapply(YERF_GOA_reconst_1972,mean)
YERF_GOA_reconst_1973_means <- lapply(YERF_GOA_reconst_1973,mean)
YERF_GOA_reconst_1974_means <- lapply(YERF_GOA_reconst_1974,mean)
YERF_GOA_reconst_1975_means <- lapply(YERF_GOA_reconst_1975,mean)
YERF_GOA_reconst_1976_means <- lapply(YERF_GOA_reconst_1976,mean)
YERF_GOA_reconst_1977_means <- lapply(YERF_GOA_reconst_1977,mean)
YERF_GOA_reconst_1978_means <- lapply(YERF_GOA_reconst_1978,mean)
YERF_GOA_reconst_1979_means <- lapply(YERF_GOA_reconst_1979,mean)
YERF_GOA_reconst_1980_means <- lapply(YERF_GOA_reconst_1980,mean)
YERF_GOA_reconst_1981_means <- lapply(YERF_GOA_reconst_1981,mean)
YERF_GOA_reconst_1982_means <- lapply(YERF_GOA_reconst_1982,mean)
YERF_GOA_reconst_1983_means <- lapply(YERF_GOA_reconst_1983,mean)
YERF_GOA_reconst_1984_means <- lapply(YERF_GOA_reconst_1984,mean)
YERF_GOA_reconst_1985_means <- lapply(YERF_GOA_reconst_1985,mean)
YERF_GOA_reconst_1986_means <- lapply(YERF_GOA_reconst_1986,mean)
YERF_GOA_reconst_1987_means <- lapply(YERF_GOA_reconst_1987,mean)}
############################################################################
## These are GOA estimates to this point.  This step is where Donnie estimates
## the proportional relationship between GOA and SEAK to get at SEAK specific estimates
## for the assessment... 
#### Prop. Southeastern YERF harvest of gulfwide YERF Bootstrapping ####
YERF_SEAK_1978_81 <- harvest$YERF_Southeastern[19:22]
YERF_GOA_1978_81 <- harvest$YERF_harvest_GOA[19:22] 

prop_SE_YERF <- data.frame(YERF_SEAK_1978_81/YERF_GOA_1978_81)
SEAK_YERF_bootstrapped <- bootstrap(prop_SE_YERF$YERF_SEAK_1978_81.YERF_GOA_1978_81,1000)

###Prop. Southeastern YERF expanded to reconstructed gulfwide YERF###
{YERF_SEAK_reconst_1960 <- Map('*',SEAK_YERF_bootstrapped,
                              YERF_GOA_reconst_1960_means)
YERF_SEAK_reconst_1961 <- Map('*',SEAK_YERF_bootstrapped,
                              YERF_GOA_reconst_1961_means)
YERF_SEAK_reconst_1962 <- Map('*',SEAK_YERF_bootstrapped,
                              YERF_GOA_reconst_1962_means)
YERF_SEAK_reconst_1963 <- Map('*',SEAK_YERF_bootstrapped,
                              YERF_GOA_reconst_1963_means)
YERF_SEAK_reconst_1964 <- Map('*',SEAK_YERF_bootstrapped,
                              YERF_GOA_reconst_1964_means)
YERF_SEAK_reconst_1965 <- Map('*',SEAK_YERF_bootstrapped,
                              YERF_GOA_reconst_1965_means)
YERF_SEAK_reconst_1966 <- Map('*',SEAK_YERF_bootstrapped,
                              YERF_GOA_reconst_1966_means)
YERF_SEAK_reconst_1967 <- Map('*',SEAK_YERF_bootstrapped,
                              YERF_GOA_reconst_1967_means)
YERF_SEAK_reconst_1968 <- Map('*',SEAK_YERF_bootstrapped,
                              YERF_GOA_reconst_1968_means)
YERF_SEAK_reconst_1969 <- Map('*',SEAK_YERF_bootstrapped,
                              YERF_GOA_reconst_1969_means)
YERF_SEAK_reconst_1970 <- Map('*',SEAK_YERF_bootstrapped,
                              YERF_GOA_reconst_1970_means)
YERF_SEAK_reconst_1971 <- Map('*',SEAK_YERF_bootstrapped,
                              YERF_GOA_reconst_1971_means)
YERF_SEAK_reconst_1972 <- Map('*',SEAK_YERF_bootstrapped,
                              YERF_GOA_reconst_1972_means)
YERF_SEAK_reconst_1973 <- Map('*',SEAK_YERF_bootstrapped,
                              YERF_GOA_reconst_1973_means)
YERF_SEAK_reconst_1974 <- Map('*',SEAK_YERF_bootstrapped,
                              YERF_GOA_reconst_1974_means)
YERF_SEAK_reconst_1975 <- Map('*',SEAK_YERF_bootstrapped,
                              YERF_GOA_reconst_1975_means)
YERF_SEAK_reconst_1976 <- Map('*',SEAK_YERF_bootstrapped,
                              YERF_GOA_reconst_1976_means)
YERF_SEAK_reconst_1977 <- Map('*',SEAK_YERF_bootstrapped,
                              YERF_GOA_reconst_1977_means)
YERF_SEAK_reconst_1978 <- Map('*',SEAK_YERF_bootstrapped,
                              YERF_GOA_reconst_1978_means)
YERF_SEAK_reconst_1979 <- Map('*',SEAK_YERF_bootstrapped,
                              YERF_GOA_reconst_1979_means)
YERF_SEAK_reconst_1980 <- Map('*',SEAK_YERF_bootstrapped,
                              YERF_GOA_reconst_1980_means)
YERF_SEAK_reconst_1981 <- Map('*',SEAK_YERF_bootstrapped,
                              YERF_GOA_reconst_1981_means)
YERF_SEAK_reconst_1982 <- Map('*',SEAK_YERF_bootstrapped,
                              YERF_GOA_reconst_1982_means)
YERF_SEAK_reconst_1983 <- Map('*',SEAK_YERF_bootstrapped,
                              YERF_GOA_reconst_1983_means)
YERF_SEAK_reconst_1984 <- Map('*',SEAK_YERF_bootstrapped,
                              YERF_GOA_reconst_1984_means)
YERF_SEAK_reconst_1985 <- Map('*',SEAK_YERF_bootstrapped,
                              YERF_GOA_reconst_1985_means)
YERF_SEAK_reconst_1986 <- Map('*',SEAK_YERF_bootstrapped,
                              YERF_GOA_reconst_1986_means)
YERF_SEAK_reconst_1987 <- Map('*',SEAK_YERF_bootstrapped,
                              YERF_GOA_reconst_1987_means)}
#####################################################################

#### SEAK YERF reconstructed yearly means and CI ####

{YERF_SEAK_mean_1960 <- mean(unlist(lapply(YERF_SEAK_reconst_1960,mean)))
YERF_SEAK_95CI_1960 <- unname(quantile(unlist(lapply(
  YERF_SEAK_reconst_1960,mean)),probs=c(0.025,0.975)))

YERF_SEAK_mean_1961 <- mean(unlist(lapply(YERF_SEAK_reconst_1961,mean)))
YERF_SEAK_95CI_1961 <- unname(quantile(unlist(
  lapply(YERF_SEAK_reconst_1961,mean)),probs=c(0.025,0.975)))     

YERF_SEAK_mean_1962 <- mean(unlist(lapply(YERF_SEAK_reconst_1962,mean)))
YERF_SEAK_95CI_1962 <- unname(quantile(unlist(
  lapply(YERF_SEAK_reconst_1962,mean)),probs=c(0.025,0.975)))

YERF_SEAK_mean_1963 <- mean(unlist(lapply(YERF_SEAK_reconst_1963,mean)))
YERF_SEAK_95CI_1963 <- unname(quantile(unlist(
  lapply(YERF_SEAK_reconst_1963,mean)),probs=c(0.025,0.975)))     

YERF_SEAK_mean_1964 <- mean(unlist(lapply(YERF_SEAK_reconst_1964,mean)))
YERF_SEAK_95CI_1964 <- unname(quantile(unlist(
  lapply(YERF_SEAK_reconst_1964,mean)),probs=c(0.025,0.975)))     

YERF_SEAK_mean_1965 <- mean(unlist(lapply(YERF_SEAK_reconst_1965,mean)))
YERF_SEAK_95CI_1965 <- unname(quantile(unlist(
  lapply(YERF_SEAK_reconst_1965,mean)),probs=c(0.025,0.975)))     

YERF_SEAK_mean_1966 <- mean(unlist(lapply(YERF_SEAK_reconst_1966,mean)))
YERF_SEAK_95CI_1966 <- unname(quantile(unlist(
  lapply(YERF_SEAK_reconst_1966,mean)),probs=c(0.025,0.975)))     

YERF_SEAK_mean_1967 <- mean(unlist(lapply(YERF_SEAK_reconst_1967,mean)))
YERF_SEAK_95CI_1967 <- unname(quantile(unlist(
  lapply(YERF_SEAK_reconst_1967,mean)),probs=c(0.025,0.975)))     

YERF_SEAK_mean_1968 <- mean(unlist(lapply(YERF_SEAK_reconst_1968,mean)))
YERF_SEAK_95CI_1968 <- unname(quantile(unlist(
  lapply(YERF_SEAK_reconst_1968,mean)),probs=c(0.025,0.975)))     

YERF_SEAK_mean_1969 <- mean(unlist(lapply(YERF_SEAK_reconst_1969,mean)))
YERF_SEAK_95CI_1969 <- unname(quantile(unlist(
  lapply(YERF_SEAK_reconst_1969,mean)),probs=c(0.025,0.975)))     

YERF_SEAK_mean_1970 <- mean(unlist(lapply(YERF_SEAK_reconst_1970,mean)))
YERF_SEAK_95CI_1970 <- unname(quantile(unlist(
  lapply(YERF_SEAK_reconst_1970,mean)),probs=c(0.025,0.975)))     

YERF_SEAK_mean_1971 <- mean(unlist(lapply(YERF_SEAK_reconst_1971,mean)))
YERF_SEAK_95CI_1971 <- unname(quantile(unlist(
  lapply(YERF_SEAK_reconst_1971,mean)),probs=c(0.025,0.975)))     

YERF_SEAK_mean_1972 <- mean(unlist(lapply(YERF_SEAK_reconst_1972,mean)))
YERF_SEAK_95CI_1972 <- unname(quantile(unlist(
  lapply(YERF_SEAK_reconst_1972,mean)),probs=c(0.025,0.975)))     

YERF_SEAK_mean_1973 <- mean(unlist(lapply(YERF_SEAK_reconst_1973,mean)))
YERF_SEAK_95CI_1973 <- unname(quantile(unlist(
  lapply(YERF_SEAK_reconst_1973,mean)),probs=c(0.025,0.975)))     

YERF_SEAK_mean_1974 <- mean(unlist(lapply(YERF_SEAK_reconst_1974,mean)))
YERF_SEAK_95CI_1974 <- unname(quantile(unlist(
  lapply(YERF_SEAK_reconst_1974,mean)),probs=c(0.025,0.975)))     

YERF_SEAK_mean_1975 <- mean(unlist(lapply(YERF_SEAK_reconst_1975,mean)))
YERF_SEAK_95CI_1975 <- unname(quantile(unlist(
  lapply(YERF_SEAK_reconst_1975,mean)),probs=c(0.025,0.975)))     

YERF_SEAK_mean_1976 <- mean(unlist(lapply(YERF_SEAK_reconst_1976,mean)))
YERF_SEAK_95CI_1976 <- unname(quantile(unlist(
  lapply(YERF_SEAK_reconst_1976,mean)),probs=c(0.025,0.975)))     

YERF_SEAK_mean_1977 <- mean(unlist(lapply(YERF_SEAK_reconst_1977,mean)))
YERF_SEAK_95CI_1977 <- unname(quantile(unlist(
  lapply(YERF_SEAK_reconst_1977,mean)),probs=c(0.025,0.975)))     

YERF_SEAK_mean_1978 <- mean(unlist(lapply(YERF_SEAK_reconst_1978,mean)))
YERF_SEAK_95CI_1978 <- unname(quantile(unlist(
  lapply(YERF_SEAK_reconst_1978,mean)),probs=c(0.025,0.975)))     

YERF_SEAK_mean_1979 <- mean(unlist(lapply(YERF_SEAK_reconst_1979,mean)))
YERF_SEAK_95CI_1979 <- unname(quantile(unlist(
  lapply(YERF_SEAK_reconst_1979,mean)),probs=c(0.025,0.975)))     

YERF_SEAK_mean_1980 <- mean(unlist(lapply(YERF_SEAK_reconst_1980,mean)))
YERF_SEAK_95CI_1980 <- unname(quantile(unlist(
  lapply(YERF_SEAK_reconst_1980,mean)),probs=c(0.025,0.975)))     

YERF_SEAK_mean_1981 <- mean(unlist(lapply(YERF_SEAK_reconst_1981,mean)))
YERF_SEAK_95CI_1981 <- unname(quantile(unlist(
  lapply(YERF_SEAK_reconst_1981,mean)),probs=c(0.025,0.975)))     

YERF_SEAK_mean_1982 <- mean(unlist(lapply(YERF_SEAK_reconst_1982,mean)))
YERF_SEAK_95CI_1982 <- unname(quantile(unlist(
  lapply(YERF_SEAK_reconst_1982,mean)),probs=c(0.025,0.975)))     

YERF_SEAK_mean_1983 <- mean(unlist(lapply(YERF_SEAK_reconst_1983,mean)))
YERF_SEAK_95CI_1983 <- unname(quantile(unlist(
  lapply(YERF_SEAK_reconst_1983,mean)),probs=c(0.025,0.975)))     

YERF_SEAK_mean_1984 <- mean(unlist(lapply(YERF_SEAK_reconst_1984,mean)))
YERF_SEAK_95CI_1984 <- unname(quantile(unlist(
  lapply(YERF_SEAK_reconst_1984,mean)),probs=c(0.025,0.975)))     

YERF_SEAK_mean_1985 <- mean(unlist(lapply(YERF_SEAK_reconst_1985,mean)))
YERF_SEAK_95CI_1985 <- unname(quantile(unlist(
  lapply(YERF_SEAK_reconst_1985,mean)),probs=c(0.025,0.975)))     

YERF_SEAK_mean_1986 <- mean(unlist(lapply(YERF_SEAK_reconst_1986,mean)))
YERF_SEAK_95CI_1986 <- unname(quantile(unlist(
  lapply(YERF_SEAK_reconst_1986,mean)),probs=c(0.025,0.975)))     

YERF_SEAK_mean_1987 <- mean(unlist(lapply(YERF_SEAK_reconst_1987,mean)))
YERF_SEAK_95CI_1987 <- unname(quantile(unlist(
  lapply(YERF_SEAK_reconst_1987,mean)),probs=c(0.025,0.975)))  }   
#####################################################################

#### Results data frame ####

Year <- c(1960:1987)

Estimate <- c(YERF_SEAK_mean_1960,YERF_SEAK_mean_1961,YERF_SEAK_mean_1962,
           YERF_SEAK_mean_1963,YERF_SEAK_mean_1964,YERF_SEAK_mean_1965,
           YERF_SEAK_mean_1966,YERF_SEAK_mean_1967,YERF_SEAK_mean_1968,
           YERF_SEAK_mean_1969,YERF_SEAK_mean_1970,YERF_SEAK_mean_1971,
           YERF_SEAK_mean_1972,YERF_SEAK_mean_1973,YERF_SEAK_mean_1974,
           YERF_SEAK_mean_1975,YERF_SEAK_mean_1976,YERF_SEAK_mean_1977,
           YERF_SEAK_mean_1978,YERF_SEAK_mean_1979,YERF_SEAK_mean_1980,
           YERF_SEAK_mean_1981,YERF_SEAK_mean_1982,YERF_SEAK_mean_1983,
           YERF_SEAK_mean_1984,YERF_SEAK_mean_1985,YERF_SEAK_mean_1986,
           YERF_SEAK_mean_1987)

L95 <- c(YERF_SEAK_95CI_1960[1],YERF_SEAK_95CI_1961[1],YERF_SEAK_95CI_1962[1],
      YERF_SEAK_95CI_1963[1],YERF_SEAK_95CI_1964[1],YERF_SEAK_95CI_1965[1],
      YERF_SEAK_95CI_1966[1],YERF_SEAK_95CI_1967[1],YERF_SEAK_95CI_1968[1],
      YERF_SEAK_95CI_1969[1],YERF_SEAK_95CI_1970[1],YERF_SEAK_95CI_1971[1],
      YERF_SEAK_95CI_1972[1],YERF_SEAK_95CI_1973[1],YERF_SEAK_95CI_1974[1],
      YERF_SEAK_95CI_1975[1],YERF_SEAK_95CI_1976[1],YERF_SEAK_95CI_1977[1],
      YERF_SEAK_95CI_1978[1],YERF_SEAK_95CI_1979[1],YERF_SEAK_95CI_1980[1],
      YERF_SEAK_95CI_1981[1],YERF_SEAK_95CI_1982[1],YERF_SEAK_95CI_1983[1],
      YERF_SEAK_95CI_1984[1],YERF_SEAK_95CI_1985[1],YERF_SEAK_95CI_1986[1],
      YERF_SEAK_95CI_1987[1])

U95 <- c(YERF_SEAK_95CI_1960[2],YERF_SEAK_95CI_1961[2],YERF_SEAK_95CI_1962[2],
      YERF_SEAK_95CI_1963[2],YERF_SEAK_95CI_1964[2],YERF_SEAK_95CI_1965[2],
      YERF_SEAK_95CI_1966[2],YERF_SEAK_95CI_1967[2],YERF_SEAK_95CI_1968[2],
      YERF_SEAK_95CI_1969[2],YERF_SEAK_95CI_1970[2],YERF_SEAK_95CI_1971[2],
      YERF_SEAK_95CI_1972[2],YERF_SEAK_95CI_1973[2],YERF_SEAK_95CI_1974[2],
      YERF_SEAK_95CI_1975[2],YERF_SEAK_95CI_1976[2],YERF_SEAK_95CI_1977[2],
      YERF_SEAK_95CI_1978[2],YERF_SEAK_95CI_1979[2],YERF_SEAK_95CI_1980[2],
      YERF_SEAK_95CI_1981[2],YERF_SEAK_95CI_1982[2],YERF_SEAK_95CI_1983[2],
      YERF_SEAK_95CI_1984[2],YERF_SEAK_95CI_1985[2],YERF_SEAK_95CI_1986[2],
      YERF_SEAK_95CI_1987[2])

Reported <- harvest$YERF_Southeastern

SEAKYERF_reconstructed_foreign_harvest <- data.frame(Year, Estimate, L95,
                                                     U95,Reported)

############################

#### CSV file writing -> change file name as neeeded ####

write.csv(SEAKYERF_reconstructed_foreign_harvest, "C:\\Users\\dearthur\\Documents\\YERF_SEAK_foreign.csv",row.names=FALSE)

#########################################################



