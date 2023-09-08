###############################################################################
## DSR directed fishery CPUE for pop model
## Deep dive of fish ticket and log book data in Fishery CPUE exam.R
## This is simplified and clean code for using Fishticket data to calculate 
## nominal CPUE values for an indices of abundance for use in population model
## This just uses fish ticket data as log book data did not add much information
## Will use lbs of YE/hook
## Could control for depth and hook spacing
## Hook size was unrelated to CPUE
##
## Phil Joy
## May 2022
###############################################################################
{library(ggplot2)
library(tidyr)
library(lubridate)
library(dplyr)
library(grid)
library(mgcv)}

#function for figuring things out for joining logbook and fish ticket data
random_check<-function(data){ #data<-tix
  #length(unique(Sable_ll_CPUE$Year))
  data<-data  #data<-whaled
  
  year_list<-unique(data$Year)
  #year_check<-sample(year_list,1)
  year_check<-year_list[sample(length(year_list),1)]
  
  data<-data %>% filter(Year == year_check)
  
  adfg_list<-unique(data$ADFG_NO)
  adfg_check<-adfg_list[sample(length(adfg_list),1)]
  
  data<-data %>% filter(ADFG_NO == adfg_check)
  
  sd_list<-unique(data$Sell_Date)
  sd_check<-sd_list[sample(length(sd_list),1)]
  
  data<-data %>% filter(Sell_Date == sd_check)
  
  gsa_list<-unique(data$G_STAT_AREA)
  gsa_check<-gsa_list[sample(length(gsa_list),1)]
  
  data<-data %>% filter(G_STAT_AREA == gsa_check)
  
  return<-list(year_check,adfg_check,sd_check,gsa_check) #rands<-return
}

# 1) Load the Data
Fishtix<-read.csv("Data_processing/Data/Harvests/Longline Hooks and Ticket Pounds.csv")
str(Fishtix)
Fishtix %>% filter(TARGET_SPECIES_CODE == 145,
                   G_MANAGEMENT_AREA_CODE != "IBS") %>% 
  mutate(Trip_Date = parse_date_time(DATE_LEFT_PORT, c("%m/%d/%Y %H:%M")),
                   Sell_Date = parse_date_time(SELL_DATE, c("%m/%d/%Y %H:%M")),
                   Jdate = as.numeric(format(Trip_Date,"%j")),
                   Year = as.factor(YEAR),
                   yrs = as.numeric(YEAR),
                   GMAC = as.factor(G_MANAGEMENT_AREA_CODE),
                   Trip_no = TRIP_NO,
                   Avg_Trip_Depth = AVG_DEPTH_FATHOMS,     
                   Spacing = HOOK_SPACING,
                   No.Hooks=HOOKS,
                   lbYE = YELLOWEYE_ROCKFISH,
                   lbYE_per_hook = YELLOWEYE_ROCKFISH/HOOKS,
                   hook_size = as.numeric(gsub("[^0-9]", "", HOOK_SIZE)),
                   std_hk_sz = ifelse (hook_size %in% c(3),16,
                                       ifelse(hook_size %in% c(4),15,
                                              ifelse(hook_size %in% c(5),14,
                                                     ifelse(hook_size %in% c(6),13,hook_size))))) -> tix

#str(tix)
#eg<-Fishtix$SELL_DATE[1]
#parse_date_time(eg, c("%m/%d/%Y %H:%M"))

lb<-read.csv("Data_processing/Data/Harvests/DSR Logbooks for Phil.csv") %>%
  mutate(Year = as.factor(Year),
         ticket_class<-as.factor(substr(CFEC.Permit.Number, 1, 1)),
         Trip_no = Trip.Number, 
         ADFG_NO = ADFG.Number,
         GMAC = Groundfish.Management.Area.Code,
         Sell_Date = parse_date_time(Sell.Date, c("%m/%d/%Y %H:%M")),
         G_STAT_AREA  = Groundfish.Stat.Area,
         Set_Depth = Average.Depth.Fathoms) %>%
  filter(Trip.Primary.Target.Species.Code == 145,
         Groundfish.Management.Area.Code != "IBS",
         Species.Code == 145) 
#M pcod+rockfish, Y = rockfish targetted
#str(lb)
hist(lb$Soak.Time.Hours)
hist(lb$Soak.Time.Hours[lb$Soak.Time.Hours>50])
quantile(lb$Soak.Time.Hours,c(0.95,0.99,0.999,0.9999),na.rm=T)
lb %>% filter(Soak.Time.Hours > 50) %>% 
  select(Year,Sell.Date,Time.Set,Time.Hauled,Soak.Time.Hours)
#data entryy errors in time hauled and set data...

colnames(tix); colnames(lb)

#x<-random_check(tix)
#tix %>% filter(Year == x[[1]][1],
#               ADFG_NO == x[[2]][1],
#               Sell_Date == x[[3]][1],
#               G_STAT_AREA == x[[4]][1]) -> egtix
#lb %>% filter(Year == x[[1]][1],
#               ADFG_NO == x[[2]][1],
#               Sell_Date == x[[3]][1],
#              G_STAT_AREA == x[[4]][1]) -> eglb

cpue_dat<-left_join(tix,lb,by = c("Year", "ADFG_NO",
                                     "Sell_Date",
                                     "G_STAT_AREA","Trip_no","GMAC"),
                    multiple="all",
                    suffix = c("_log","_ftx"))  %>% 
  group_by(Year, G_STAT_AREA, ADFG_NO, Sell_Date, Trip_no) %>% 
  mutate(No_sets = n(),
         Total_soak = sum(Soak.Time.Hours),
         #Tot_numbers = sum(Numbers),
         #Tot_pounds = sum(Pounds),
         set_prop_no = Numbers/sum(Numbers),
         set_prop_lbs = Pounds/sum(Pounds),
         set_units = ifelse(is.na(Numbers),
                          "Pounds", #"Pounds", #unique(Pounds),
                          "Numbers"),
         set_lbs = ifelse(set_units == "Numbers",
                          lbYE*set_prop_no,
                          lbYE*set_prop_lbs),
         set_no_hooks = No.Hooks/No_sets,
         set_lbYE_per_hook = set_lbs/set_no_hooks) %>% ungroup() %>% 
  select(Year, Sell_Date, Trip_Date, Jdate, Ticket, 
       Trip_no,ADFG_NO,Trip.Primary.Target.Species,Effort.Primary.Target.Species,
       GMAC,G_STAT_AREA, Effort.Number, No_sets,
       hook_size,std_hk_sz,Spacing,hps = HOOKS_PER_SKATE,No.Hooks,set_no_hooks,
       Set_soak = Soak.Time.Hours, Total_soak, 
       Avg_Trip_Depth, Set_Depth, 
       Disposition, Disposition.Code, 
       YELLOWEYE_ROCKFISH,lbYE,Numbers,Pounds, lbYE_per_hook, set_lbs, set_lbYE_per_hook) %>% data.frame() #%>% filter(!is.na(Fishery.Name))

colnames(cpue_dat)

x<-random_check(tix)

cpue_dat %>% filter(Year == x[[1]][1],
              ADFG_NO == x[[2]][1],
              Sell_Date == x[[3]][1],
              G_STAT_AREA == x[[4]][1]) 

#cpue %>% 
#  select(Year, yrs, Subd, Trip_Date, Sell_Date, Trip_no, 
#         ADFG_NO, 
#         Jdate, Depth,
#         Spacing, No.Hooks, hook_size, std_hk_sz, 
#         Gear.Class.Code, Gear.Code,
#         lbYE, lbYE_per_hook,lbYE_per_hook )->cpue#

#cpue <- cpue %>% filter (!is.na(lbYE_per_hook), !is.infinite(lbYE_per_hook))
#nrow(cpue[is.na(cpue$lbYE_per_hook),])

#quick check of NA's and Infinites of key variables
#nrow(cpue[is.na(cpue$lbYE_per_hook) | is.infinite(cpue$lbYE_per_hook),])
#nrow(cpue[is.na(cpue$Jdate) | is.infinite(cpue$Jdate),])
#nrow(cpue[is.na(cpue$Depth) | is.infinite(cpue$Depth),])
#nrow(cpue[is.na(cpue$Spacing) | is.infinite(cpue$Spacing),])
#nrow(cpue[is.na(cpue$No.Hooks) | is.infinite(cpue$No.Hooks),])
#nrow(cpue[is.na(cpue$std_hk_sz) | is.infinite(cpue$std_hk_sz),])
#nrow(cpue[is.na(cpue$hook_size) | is.infinite(cpue$hook_size),])
#Look at the data again...

colnames(cpue_dat)

ggplot(data=cpue_dat, 
       aes(Year, set_lbYE_per_hook))+
  geom_boxplot(fill="lightblue", colour="blue") +
  facet_grid(GMAC ~ ., scales = "free")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="lightblue2")) +
  xlab("Year") +
  ylab("lbs yelloweye per hook") +
  labs(title="Raw YE/hook")

ggplot(data=cpue_dat, 
       aes(Year, lbYE))+
  geom_boxplot(fill="violet", colour="purple") +
  facet_grid(GMAC ~ ., scales = "free")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="violet")) +
  xlab("Year") +
  ylab("lbs yelloweye") +
  labs(title="lbs of YE")

#Look at relationship between hook variables, depth and CPUE (YE/hook)
#unique(cpue$HOOK_SIZE)

#plot(cpue$lbYE_per_hook ~ cpue$std_hk_sz)
#lines(lm(cpue$lbYE_per_hook ~ cpue$std_hk_sz))

ggplot(data=cpue_dat, 
       aes(as.factor(std_hk_sz), lbYE_per_hook))+
  geom_boxplot(fill="grey", colour="black") +
  #facet_grid(HOOK_SIZE ~ ., scales = "free")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="grey")) +
  xlab("hook size") +
  ylab("YE/hook") 

#plot(cpue$lbYE_per_hook ~ cpue$Spacing)
#lines(lm(cpue$lbYE_per_hook ~ cpue$Spacing))
ggplot(data=cpue_dat, 
       aes(as.factor(Spacing), lbYE_per_hook))+
  geom_boxplot(fill="grey", colour="black") +
  #facet_grid(HOOK_SIZE ~ ., scales = "free")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="grey"))

ggplot(data=cpue_dat, 
       aes(Set_Depth, set_lbYE_per_hook))+
  geom_point() + #geom_line()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="grey"))

#annual trends in some variables... 
ggplot(data=cpue_dat, 
       aes(as.factor(Year), Set_Depth))+
  geom_boxplot(fill="grey", colour="black") +
  facet_grid(GMAC ~ ., scales = "free")+ #geom_line()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="grey"))

# Cull out the few extreme depths...
#cpue<-cpue[cpue$Depth < 150,]
#head(cpue, 20)
#str(cpue)
#===============================================================================
# NOMINAL CPUE for EACH AREA
#library(plyr)
cpue_dat %>% 
  group_by(Year,GMAC) %>% 
  dplyr::summarise(yrs = Year,
                   fsh_cpue = mean(lbYE_per_hook,na.rm=T),
                   sd = sd(lbYE_per_hook,na.rm=T),
                   n = length(lbYE_per_hook),
                   se = sd / (n ^ (1/2)),
                   var = var(lbYE_per_hook,na.rm=T),
                   cv = sd / fsh_cpue,
                   upper = fsh_cpue + (2 * se),
                   lower = fsh_cpue - (2 * se))  -> fsh_sum 

View(fsh_sum)
fsh_sum<-fsh_sum[complete.cases(fsh_sum),]
str(fsh_sum)

#dplyr bug leaving duplicate rows
fsh_sum2<-unique(fsh_sum)
View(fsh_sum2)

fsh_sum<-fsh_sum[order(fsh_sum$GMAC, fsh_sum$Year),]
nYear<-length(unique(fsh_sum$Year))

fsh_sum2 %>% 
  ggplot() +
  geom_ribbon(aes(Year, ymin = lower, ymax = upper, fill = GMAC), 
              colour = "white", alpha = 0.2, outline.type="full") +
  geom_point(aes(Year, fsh_cpue, colour = GMAC, shape = GMAC), size = 2) +
  geom_line(aes(Year, fsh_cpue, colour = GMAC, group = GMAC), size = 1) +
  facet_grid(GMAC ~ ., scales = "free")+
  # scale_colour_grey(name = "Standardized CPUE") +
  # scale_fill_grey(name = "Standardized CPUE") +
  #scale_colour_manual(values = c("darkcyan", "goldenrod"), name = "Standardized CPUE") +
  #scale_fill_manual(values = c("darkcyan", "goldenrod"), name = "Standardized CPUE") +
  #scale_shape_manual(values = c(19, 17), name = "Standardized CPUE") +
  #scale_x_continuous(breaks = axis$breaks, labels = axis$labels) + 
  labs(x = "", y = "Nominal Fishery CPUE (round lb/hook)\n") +
  theme(legend.position = c(0.8, 0.8)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="grey")) #+
  #expand_limits(y = 5)

View(fsh_sum)
fsh_sum<-fsh_sum2

#before saving for model will be easier to insert NA values for missing years here
#would be tedious to do it in the data loading script for the SPM model

allyrs<-seq(1986, 2022,1)
SD<-unique(fsh_sum$GMAC)
#doing it witha clumsy loop...
cpue_model<-data.frame(); i<-1

colnames(fsh_sum)

for (s in SD){
  for (y in allyrs){
    Dat<-fsh_sum[fsh_sum$yrs == y & fsh_sum$GMAC == s,]
    if (nrow(Dat) == 0) {
      cpue_model[i,"Year"]<-y
      cpue_model[i,"GMAC"]<-s
      cpue_model[i,"yrs"]<-y
      cpue_model[i,"fsh_cpue"]<-NA
      cpue_model[i,5:8]<-NA
      cpue_model[i,"cv"]<-1
      cpue_model[i,10:11]<-NA
      i<-i+1
    } else {
      cpue_model[i,"Year"]<-y
      cpue_model[i,"GMAC"]<-s
      cpue_model[i,3:11]<-Dat[,3:11]
      i<-i+1
    }
  }
}
colnames(cpue_model)<-colnames(fsh_sum)
View(cpue_model)

cpue_model %>% 
  ggplot() +
  geom_ribbon(aes(Year, ymin = lower, ymax = upper, fill = GMAC), 
              colour = NA, alpha = 0.2, outline.type="full") +
  geom_point(aes(Year, fsh_cpue, colour = GMAC, shape = GMAC), size = 2) +
  geom_line(aes(Year, fsh_cpue, colour = GMAC, group = GMAC), size = 1) +
  facet_grid(GMAC ~ ., scales = "free")+
  labs(x = "", y = "Fishery CPUE (round lb/hook)\n") +
  theme(legend.position = "right") + #c(0.8, 0.8)) +
  scale_x_continuous(breaks = cpue_model$Year, labels = cpue_model$Year) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="grey")) +
  expand_limits(y = 0)

YEAR<-2023
ggsave(paste0("Figures/nom_cpue_1986_", YEAR, ".png"), 
       dpi=400, height=5, width=7.5, units="in")

write.csv(cpue_model, "Data_processing/Data/YE_dirfishery_nom_CPUE.csv")

#================================================================================
# Normalization 

#get rid of NAs and Inf for exam
cpue.com<-cpue_dat %>% filter(lbYE_per_hook != Inf) #cpue[complete.cases(cpue),]
nrow(cpue)
nrow(cpue.com)
colnames(cpue.com)
unique(cpue.com$Disposition)

ggplot(cpue.com, aes(lbYE_per_hook)) + geom_density(alpha = 0.4, fill = 4)
ggplot(cpue.com, aes(set_lbYE_per_hook)) + geom_density(alpha = 0.4, fill = 4)

# Better, but still not normal with log transformation
ggplot(cpue.com, aes(log(lbYE_per_hook + 1))) + geom_density(alpha = 0.4, fill = 4)
ggplot(cpue.com, aes(log(set_lbYE_per_hook + 1))) + geom_density(alpha = 0.4, fill = 4)
# Following Jenny Stahl and Ben Williams' work in the SSEI, increase CPUE by 10%
# of the mean per Cambell et al 1996 and Cambell 2004. Back-transform with
# exp(cpue - mean(fsh_cpue$std_cpue) * 0.1)

mean(cpue.com$lbYE_per_hook, na.rm=T)
max(cpue.com$lbYE_per_hook, na.rm=T)
min(cpue.com$lbYE_per_hook, na.rm=T)

cpue.com %>% #filter()
  mutate(cpue = log(lbYE_per_hook + (mean(cpue.com$lbYE_per_hook, na.rm=T) * 0.1)),
         set_cpue = log(set_lbYE_per_hook + (mean(cpue.com$set_lbYE_per_hook, na.rm=T) * 0.1))) -> cpue.com

ggplot(cpue.com %>% filter(!is.na(cpue)), aes(cpue)) + geom_density(alpha = 0.4, fill = 4)
ggplot(cpue.com %>% filter(!is.na(set_cpue)), aes(set_cpue)) + geom_density(alpha = 0.4, fill = 4)

#-------------------------------------------------------------------------------
# Look at variables affecting CPUE
# Hook Size
ggplot(cpue.com, aes(as.factor(std_hk_sz), cpue)) + geom_boxplot()
ggplot(cpue.com, aes(Year, cpue, fill = as.factor(std_hk_sz))) + geom_boxplot()+
  theme(axis.text.x = element_text(size = 14, angle = 90, h = 1)) +
  labs(x = "", y = "Fishery CPUE\n")
ggplot(cpue.com, aes(GMAC, cpue, fill = as.factor(std_hk_sz))) + geom_boxplot()+
  labs(x = "\nStat area", y = "Fishery CPUE\n")

#Depth
ggplot(cpue.com, aes(Avg_Trip_Depth, cpue)) + geom_point(shape = 20) + 
  geom_smooth(size = 2, se = FALSE) 

ggplot(cpue.com, aes(Set_Depth, set_cpue)) + geom_point(shape = 20) + 
  geom_smooth(size = 2, se = FALSE) 

#Hook Spacing
ggplot(cpue.com, aes(Spacing, cpue)) + geom_point(shape = 20) + 
  geom_smooth(size = 1, se = FALSE) 

ggplot(cpue.com, aes(as.factor(Spacing), cpue)) + geom_boxplot()

#Soak time
hist(cpue.com$Total_soak)
cpue.com %>% filter(Total_soak > 100000)

ggplot(cpue.com %>% filter(Total_soak < 100), aes(Total_soak, cpue)) + geom_point(shape = 20) + 
  geom_smooth(size = 2, se = FALSE) 

ggplot(cpue.com %>% filter(Total_soak < 100 & Set_soak > 0), aes(Set_soak, set_cpue)) + geom_point(shape = 20) + 
  geom_smooth(size = 2, se = FALSE) 

#Julian date
ggplot(cpue.com, aes(Jdate, cpue)) +
  geom_point(shape = 20) + 
  geom_smooth(method = 'loess', span = 1, se = FALSE) 

ggplot(cpue.com, aes(Jdate, cpue, group = Year, colour = Year)) +
  geom_smooth(method = 'loess', span = 1, se = FALSE) 

#### Additive modeling... 
# get rid of those bad soak times... 
cpue_exam <- cpue.com %>% filter(Total_soak < 100,
                                 Set_soak > 0,
                                 !is.na(set_lbYE_per_hook )) %>%
  mutate(dum = 1, dumstat = 1, GMAC = as.factor(GMAC)) %>% data.frame()

colnames(cpue_exam)

m0 <- bam(set_cpue ~ Year + s(GMAC, bs='re', by=dumstat), data=cpue_exam, gamma=1.4)
m0.hook <- bam(set_cpue ~ Year + s(GMAC, bs='re', by=dumstat) + std_hk_sz, data=cpue_exam, gamma=1.4)
m0.spacing <- bam(set_cpue ~ Year + s(GMAC, bs='re', by=dumstat) + s(Spacing, k=4), data=cpue_exam, gamma=1.4)
m0.depth <- bam(set_cpue ~ Year + s(GMAC, bs='re', by=dumstat) + s(Set_Depth, k=4), data=cpue_exam, gamma=1.4)
m0.soak <- bam(set_cpue ~ Year + s(GMAC, bs='re', by=dumstat) + s(Set_soak, k=4), data=cpue_exam, gamma=1.4)
m0.stat <- bam(set_cpue ~ Year + s(GMAC, bs='re', by=dumstat) + s(G_STAT_AREA, bs='re', by=dumstat), data=cpue_exam, gamma=1.4)
m0.boat <- bam(set_cpue ~ Year + s(GMAC, bs='re', by=dumstat) + s(ADFG_NO, bs='re', by=dumstat), data=cpue_exam, gamma=1.4)
m0.jday <- bam(set_cpue ~ Year + s(GMAC, bs='re', by=dumstat) + s(Jdate, k=4), data=cpue_exam, gamma=1.4)

model.list<-list(m0,m0.hook,m0.spacing,m0.depth,m0.soak,m0.stat,m0.boat,
                 #0.lat_lon, m0.lat,m0.lon,
                 m0.jday)
names(model.list)<-c("m0","hook","spacing","depth","soak","stat","boat",#"lat_lon",
                     #"lat","lon",
                     "jday")
modsum0<-data.frame(); j<-1
for (i in model.list) {
  #mod<-i
  modsum0[j,"model"]<-names(model.list[j])
  modsum0[j,"aic"]<-AIC(i)
  modsum0[j,"dev"]<-summary(i)$dev.expl
  modsum0[j,"rsq"]<-summary(i)$r.sq
  modsum0[j,"dev_exp"]<-summary(i)$dev.expl-summary(m0)$dev.expl
  j<-j+1
}

modsum0 %>% arrange(aic)  
modsum0 %>% arrange(-dev)  
modsum0 %>% arrange(-rsq) 

# hook size, hook spacing and depth are the most important, but jdat and location haave an effect
# boat doesn't seem to have much pull

# i'd like to do a fuller analysis, but for simplicitiy lets standardize for everything

global <- bam(set_cpue ~ Year + s(GMAC, bs='re', by=dumstat) + std_hk_sz +
                s(Spacing, k=4) + s(Set_Depth, k=4) + s(Set_soak, k=4) +
                s(G_STAT_AREA, bs='re', by=dumstat) + s(ADFG_NO, bs='re', by=dumstat) +
                s(Jdate, k=4),
              data=cpue_exam, gamma=1.4)

plot(global, page = 1, shade = TRUE, all = TRUE)

vis.gam(global, c('Set_Depth', 'Spacing'), plot.type='contour', type='response', color='topo', too.far=0.1)
vis.gam(global, c('Set_Depth', 'std_hk_sz'), plot.type='contour', type='response', color='topo', too.far=0.1)

plot.gam(global)

str(diag(vcov.gam(global)))

#Create standard dataset to get standardized CPUE for each year
with(cpue_exam,table(G_STAT_AREA,GMAC))

set_cpue$G_STAT_AREA[cpue_exam$GMAC == "EYKT"]

MA<-unique(cpue_exam$GMAC)

for (i in 1:4){  #i<-1
  std_dat_sub <- expand.grid(year = unique(cpue_exam$Year[cpue_exam$GMAC == MA[i]]),
                         GMAC = MA[i],
                         std_hk_sz = mean(cpue_exam$std_hk_sz[cpue_exam$GMAC == MA[i]], na.rm=T),
                         Spacing = mean(cpue_exam$Spacing[cpue_exam$GMAC == MA[i]], na.rm=T),
                         Set_Depth = mean(cpue_exam$Set_Depth[cpue_exam$GMAC == MA[i]], na.rm=T), 
                         Set_soak = mean(cpue_exam$Set_soak[cpue_exam$GMAC == MA[i]]), 
                         #G_STAT_AREA = cpue_exam$G_STAT_AREA[cpue_exam$GMAC == MA[i]][1],
                         G_STAT_AREA = ifelse(MA[i] == "CSEO",365701,
                                              ifelse(MA[i] == "EYKT",385800,
                                                     ifelse(MA[i] == "NSEO",365733,335431))),
                         ADFG_NO = cpue_exam$ADFG_NO[cpue_exam$GMAC == MA[i]][1],
                         Jdate = mean(cpue_exam$Jdat[cpue_exam$GMAC == MA[i]], na.rm=T),
                         dum = 0,
                         dumstat = 0) %>% 
    mutate(Year = factor(year))
  
  if (i == 1) {std_dat <- std_dat_sub} else {std_dat <- rbind(std_dat,std_dat_sub)}
}

pred_cpue <- predict(global, std_dat, type = "link", se = TRUE)

#Put the standardized CPUE and SE into the data frame and convert to
#backtransformed (bt) CPUE
std_dat %>% 
  mutate(fit = pred_cpue$fit,
         se = pred_cpue$se.fit,
         upper = fit + (2 * se),
         lower = fit - (2 * se),
         bt_cpue = exp(fit) - (mean(cpue_exam$set_cpue) * 0.1),
         bt_upper = exp(upper) - (mean(cpue_exam$set_cpue) * 0.1),
         bt_lower = exp(lower) - (mean(cpue_exam$set_cpue) * 0.1),
         bt_se = (bt_upper - bt_cpue) / 2  #,
         #bt_cv = bt_se/bt_cpue
  ) -> std_dat

# Nominal CPUE ----

cpue_exam %>% 
  group_by(Year,GMAC) %>% 
  dplyr::summarise(fsh_cpue = mean(set_lbYE_per_hook),
                   sd = sd(set_lbYE_per_hook),
                   n = length(set_lbYE_per_hook),
                   se = sd / (n ^ (1/2)),
                   var = var(set_lbYE_per_hook),
                   cv = sd / fsh_cpue,
                   upper = fsh_cpue + (2 * se),
                   lower = fsh_cpue - (2 * se)) -> fsh_sum 

#need to fill in missing years for plotting... 
allyrs<-seq(1987, 2022,1)
SD<-unique(fsh_sum$GMAC)
#doing it witha clumsy loop...
cpue_nom<-data.frame(); i<-1
cpue_std <- cpue_nom

colnames(fsh_sum)

for (s in SD){  #s <-SD[1]
  for (y in allyrs){ #y<-allyrs[2]
    nom<-fsh_sum[fsh_sum$Year == y & fsh_sum$GMAC == s,]
    std<-std_dat[std_dat$Year == y & std_dat$GMAC == s,]
    if (nrow(nom) == 0) {
      cpue_nom[i,"Year"]<-y
      cpue_nom[i,"GMAC"]<-s
      #cpue_nom[i,"yrs"]<-y
      cpue_nom[i,"fsh_cpue"]<-NA
      cpue_nom[i,5:7]<-NA
      cpue_nom[i,"cv"]<-1
      cpue_nom[i,9:10]<-NA
      cpue_nom[i,"CPUE"]<-"raw"
      cpue_std[i,"Year"]<-y
      cpue_std[i,"GMAC"]<-s
      #cpue_std[i,"yrs"]<-y
      cpue_std[i,"fsh_cpue"]<-NA
      cpue_std[i,5:7]<-NA
      cpue_std[i,"cv"]<-1
      cpue_std[i,9:10]<-NA
      cpue_std[i,"CPUE"]<-"standardized"
      i<-i+1
    } else {
      cpue_nom[i,"Year"]<-y
      cpue_nom[i,"GMAC"]<-s
      cpue_nom[i,3:10]<-nom[,3:10]
      cpue_nom[i,"CPUE"]<-"raw"
      cpue_std[i,"Year"]<-y
      cpue_std[i,"GMAC"]<-s
      cpue_std[i,"fsh_cpue"]<-std[1,"bt_cpue"]
      cpue_std[i,"sd"]<-std[1,"se"]
      cpue_std[i,"n"]<-nom[1,"n"]
      cpue_std[i,"se"]<-std[1,"se"]
      cpue_std[i,"var"]<-std[1,"se"]^2
      cpue_std[i,"cv"]<-NA
      cpue_std[i,"upper"]<-std[1,"bt_upper"]
      cpue_std[i,"lower"]<-std[1,"bt_lower"]
      cpue_std[i,"CPUE"]<-"standardized"
      i<-i+1
    }
  }
}
#colnames(cpue_nom)<-colnames(fsh_sum)
#colnames(cpue_std)<-colnames(fsh_sum)
cpue_ests<-rbind(cpue_nom,cpue_std)
unique(cpue_ests$CPUE)
# Compare predicted cpue from gam to nominal cpue
cpue_ests %>% 
  ggplot() +
  geom_ribbon(aes(Year, ymin = lower, ymax = upper, fill = CPUE), 
              colour = NA, alpha = 0.2) +
  geom_point(aes(Year, fsh_cpue, colour = CPUE, shape = CPUE), size = 2) +
  geom_line(aes(Year, fsh_cpue, colour = CPUE, group = CPUE), size = 1) +
  facet_wrap(~ GMAC) +
  # scale_colour_grey(name = "Standardized CPUE") +
  # scale_fill_grey(name = "Standardized CPUE") +
  scale_colour_manual(values = c("darkcyan", "goldenrod"), name = "Standardized CPUE") +
  scale_fill_manual(values = c("darkcyan", "goldenrod"), name = "Standardized CPUE") +
  scale_shape_manual(values = c(19, 17), name = "Standardized CPUE") +
  labs(x = "", y = "Fishery CPUE (round lb/hook)\n") +
  theme(legend.position = "top") +
  expand_limits(y = 0) +
  xlab("Year") +
  scale_x_continuous(breaks = seq(1986,2022,2), labels = seq(1986,2022,2)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="grey")) 
  
YEAR<-2023
ggsave(paste0("figures/compare_stdcpue_llfsh_", YEAR, ".png"), dpi=300, height=7, width=7, units="in")

write.csv(cpue_ests, "Data_processing/Data/YE_dirfishery_std_CPUE.csv")









################################################################################
## SCRAP 
#################################################################################
m1 <- bam(lbYE_per_hook ~ Year + Jdate + std_hk_sz + 
            s(Avg_Trip_Depth, k=4) + s(Spacing, k=4) +
            GMAC,
          data=cpue.com, gamma=1.4)
m2 <- bam(lbYE_per_hook ~ Year + Jdate + #DROP:  std_hk_sz + 
            s(Depth, k=4) + s(Spacing, k=4)+
            Subd,
          data=cpue.com, gamma=1.4)
m3 <- bam(lbYE_per_hook ~ Year  + std_hk_sz + #DROP + Jdate
            s(Depth, k=4) + s(Spacing, k=4) +
            Subd,
          data=cpue.com, gamma=1.4)
m4 <- bam(lbYE_per_hook ~ Year + Jdate + std_hk_sz + 
            s(Spacing, k=4) +
            Subd,
          data=cpue.com, gamma=1.4)
m5 <- bam(lbYE_per_hook ~ Year + Jdate + std_hk_sz + 
            s(Depth, k=4) +
            Subd,
          data=cpue.com, gamma=1.4)
summary(m1) 
summary(m2)
summary(m3)
summary(m4)
summary(m5)

AIC(m1, m2,m3, m4,m5)

#Hook spacing and depth significant, but Rsq not very exciting
#AIC likes Jdate and hook size, although not significant
# RE of subdistrict very significant

m6 <- bam(lbYE_per_hook ~ Year + # Drop both: Jdate + std_hk_sz + 
            s(Depth, k=4) + s(Spacing, k=4) + 
            Subd,
          #s(subd, bs='re', by=dumstat),
          data=cpue.com, gamma=1.4)

summary(m6)

AIC(m1,m2,m3,m4,m5,m6)

plot(fitted(m1), resid(m1))
abline(h=0,col="red")
plot(fitted(m6), resid(m6))
abline(h=0,col="red")

plot(m1, page=1, shade=TRUE, resid=TRUE, all=TRUE)
plot(m2, page=1, shade=TRUE, resid=TRUE, all=TRUE)
plot(m6, page=1, shade=TRUE, resid=TRUE, all=TRUE)

vis.gam(m6,c('Depth','Spacing'),plot.type='contour',
        type='response', color='topo',too.far=0.1)

## Go forward with m6 to compare to nominal...
#predicted values to compare to nominal cpue
std_dat <- expand.grid(Year = unique(cpue.com$Year),
                       Depth = mean(cpue.com$Depth), 
                       Spacing = mean(cpue.com$Spacing), 
                       Subd = unique(cpue.com$Subd),
                       #julian_day = median(Complete_FT_YE$Jdat),
                       dum = 0,
                       dumstat = 0) 

pred_cpue <- predict(m6, std_dat, type = "link", se = TRUE)

#checking my code with Jane's... checks out :)
preds<-predict.bam(m6, type="response", std_dat, se = TRUE)
str(preds); head(preds)

#Put the standardized CPUE and SE into the data frame and convert to
#backtransformed (bt) CPUE

std_dat %>% 
  mutate(fit = pred_cpue$fit,
         se = pred_cpue$se.fit,
         upper = fit + (2 * se),
         lower = fit - (2 * se),
         bt_cpue = exp(fit) - (mean(cpue.com$lbYE_per_hook) * 0.1),
         bt_upper = exp(upper) - (mean(cpue.com$lbYE_per_hook) * 0.1),
         bt_lower = exp(lower) - (mean(cpue.com$lbYE_per_hook) * 0.1),
         bt_se = (bt_upper - bt_cpue) / 2  #,
         #bt_cv = bt_se/bt_cpue
  ) -> std_dat

# Compare predicted cpue from gam to nominal cpue
fsh_sum %>%
  select(Year, Subd, cpue = fsh_cpue, upper, lower) %>% 
  mutate(CPUE = "Nominal") %>% 
  bind_rows(std_dat %>% 
              #select(YEAR, cpue = bt_cpue, upper = bt_upper, lower = bt_lower) %>% 
              select(Year, Subd, cpue = fit, upper , lower) %>%
              mutate(CPUE = "GAM")) -> comp
View(comp)
##!!! bt_cpue is awful high relative to nominal... something going on need to figure out... 

comp%>% 
  ggplot() +
  geom_ribbon(aes(Year, ymin = lower, ymax = upper, fill = CPUE), 
              colour = "white", alpha = 0.2) +
  geom_point(aes(Year, cpue, colour = CPUE, shape = CPUE), size = 2) +
  geom_line(aes(Year, cpue, colour = CPUE, group = CPUE), size = 1) +
  facet_grid(Subd ~ ., scales = "free")+
  # scale_colour_grey(name = "Standardized CPUE") +
  # scale_fill_grey(name = "Standardized CPUE") +
  scale_colour_manual(values = c("darkcyan", "goldenrod"), name = "Standardized CPUE") +
  scale_fill_manual(values = c("darkcyan", "goldenrod"), name = "Standardized CPUE") +
  scale_shape_manual(values = c(19, 17), name = "Standardized CPUE") +
  #scale_x_continuous(breaks = axis$breaks, labels = axis$labels) + 
  labs(x = "", y = "Fishery CPUE (round lb/hook)\n") +
  theme(legend.position = c(0.8, 0.2)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="grey")) +
  expand_limits(y = 0)



ggsave(paste0("figures/compare_stdcpue_llfsh_2022reboot_", YEAR, ".png"), dpi=300, height=4, width=7, units="in")

