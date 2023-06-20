#############################################################
## Observed bycatch estimates in halibut fishery
## for use in SS-SPM YE model
## Data from Jane Sullivan at NOAA per special request to get this info
## at the subdistrict level
## June 2022 - Phil Joy
#############################################################

library("dplyr")
library("ggplot2")

Obs_mean<-read.csv("Data/obs_mean_bycatch.csv")
Obs_sum<-read.csv("Data/obs_bycatch_sum.csv")
Obs_byhaul<-read.csv("Data/obs_bycatch_byhaul.csv")
str(Obs_byhaul)
Obs_byhaul<-Obs_byhaul[order(Obs_byhaul$adfg_area, Obs_byhaul$year),]
Obs_byhaul %>%
  mutate(adfg_area = as.factor(adfg_area)) -> Obs_byhaul
## get total bycatch for each subdistrict and year
## can apply small cv's when fitting this data in model since it is "known"  

##should examine how bycatch rates compare to WCPUE estimates... 
unique(Obs_byhaul$adfg_area)
Obs_byhaul %>%
  filter(year > 2007) %>% 
  group_by(adfg_area, year) %>%
  dplyr::summarise(hal_tot = sum(halibut),
         ye_tot=sum(yelloweye),
         bycatch_rate = ye_tot/hal_tot) %>%
  select(adfg_area, year, hal_tot, ye_tot, bycatch_rate)->Tot_bycatch

str(Tot_bycatch)

View(Tot_bycatch)

Tot_bycatch %>% 
  ggplot(aes(x=year,y=bycatch_rate)) +
  geom_line() + 
  facet_wrap(adfg_area ~ .)

Tot_bycatch %>% 
  ggplot(aes(x=year,y=ye_tot)) +
  geom_line() + 
  facet_wrap(adfg_area ~ .)

#Seems like pre-2013 data is biased low; guessing there is something going on
# there and should stick to just using data from 2013 through present...
  
SEO_bycatch<-Tot_bycatch %>% 
  filter(adfg_area == "EYKT" |
           adfg_area == "NSEO" |
           adfg_area == "CSEO" |
           adfg_area == "SSEO") %>%
  filter(year > 2012)

SEO_bycatch %>% 
  ggplot(aes(x=year,y=bycatch_rate)) +
  geom_line() + 
  facet_wrap(adfg_area ~ .)

465/22179

#check with summary data:
Obs_sum %>% filter(subdistrict != "NSEI" & subdistrict != "SSEI" &
                     subdistrict !="IBS") -> Obs_sum

Obs_sum[Obs_sum$subdistrict == "CSEO" & Obs_sum$year > 2012,]
SEO_bycatch[SEO_bycatch$adfg_area == "CSEO",]


#Check with what I did with Jane's Dec. 2021 exam...
SEO_bycatch %>% group_by(year) %>%
  dplyr::summarize(tot.ye = sum(ye_tot)) %>% 
  select (year,tot.ye) -> byhaul_quick

Obs_sum %>% filter (year > 2012) %>%
  group_by(year) %>%
  dplyr::summarize(tot.ye = sum(ye_bycatch)) %>% 
  select (year,tot.ye) -> obs_quick

Janes21dat<-read.csv("Data/SEO_YE_NOAA_Observers.csv")
str(Janes21dat)

byhaul_quick; obs_quick; Janes21dat$YE.Bycatch

conv.J<-Janes21dat[Janes21dat$Year < 2021,]

plot(byhaul_quick$tot.ye ~ conv.J$YE.Bycatch)

##SHIT! Total disagreement from two data sources... 
# double check this... 
Obs.1<-read.csv("Data/Groundfish Total Discards by Species.csv", header=T) %>%
  dplyr::rename(Year =  ?..Year)

str(Obs.1)

Obs.1 %>% filter(Species.Code == 145, FMP.Subarea == "SE") %>% 
  group_by(Year) %>%
  dplyr::summarize(YE.bycatch = sum(Total.Catch..mt.)) -> YE.bycatch

View(YE.bycatch)

YE.bycatch
byhaul_quick

#Check other data Jane sent 
View(Obs_mean)

##!!! Yuck: Data from "Groundfish Total Discards by Species" and 
## CAS data -> model
# "obs_bycatch_byhaul.csv" do not match up! Need to resolve this with
# Jane.  Will talk about it on June-14-2022 when I meet with her!! Until then
# using schwaggy data currently in DATALOAD script... 
