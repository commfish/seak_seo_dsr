################################################################################
## DSR Directed Fishery CPUE
## For NPFMC PT and SSC fall 2022
## Phil Joy
################################################################################
library(ggplot2)
library(tidyr)
library(lubridate)
library(dplyr)
library(grid)
library(mgcv)

# 1) Load the Data

LB<-read.csv("Data/DSR Logbooks for Phil.csv")

Fishtix<-read.csv("Data/Longline Hooks and Ticket Pounds.csv")
Fishtix %>% mutate(Trip_Date = parse_date_time(DATE_LEFT_PORT, c("%m/%d/%Y %H:%M")),
                   Sell_Date = parse_date_time(SELL_DATE, c("%m/%d/%Y %H:%M")))->Fishtix

#NOTE: Species numbers are round pounds landed! 

str(LB); nrow(LB)
unique(LB$Trip.Primary.Target.Species.Code)   #Species Code 145 = YE
unique(LB$Year)

str(Fishtix); nrow(Fishtix)
unique(Fishtix$TARGET_SPECIES_CODE)
unique(Fishtix$YEAR)

# Get Data from trips marked as targetting YE...

LB_YE<-LB[LB$Trip.Primary.Target.Species.Code == 145,]; nrow(LB_YE)
FT_YE<-Fishtix[Fishtix$TARGET_SPECIES_CODE ==145 & 
                      Fishtix$G_MANAGEMENT_AREA_CODE != "IBS",]; nrow(FT_YE)
#Fishtix_YE$YEAR <-as.factor(Fishtix_YE$YEAR)

#=================================================================================
# Examine and clean Fish tiecket Data
# Calculate YE/Hook
FT_YE$YE_per_hook <- FT_YE$YELLOWEYE_ROCKFISH/FT_YE$HOOKS
FT_YE$Jdate<-as.numeric(format(FT_YE$Trip_Date,"%j"))
str(FT_YE)

#Raw YE/hook
ggplot(data=FT_YE, 
       aes(as.factor(YEAR), YE_per_hook))+
  geom_boxplot(fill="palegreen", colour="forestgreen") +
  facet_grid(G_MANAGEMENT_AREA_CODE ~ ., scales = "free")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="lightgreen")) +
  xlab("Year") +
  ylab("YE/hook") +
  labs(title="Raw YE/hook")

#Total catch by year
FT_YE %>%
  group_by(G_MANAGEMENT_AREA_CODE, YEAR) %>%
  summarize(total.YE = sum(YELLOWEYE_ROCKFISH, na.rm=T)) -> sum_catch

ggplot(sum_catch, 
       aes(x = YEAR, y = total.YE, colour=G_MANAGEMENT_AREA_CODE)) +
  geom_point() + geom_line()

View(sum_catch[sum_catch$G_MANAGEMENT_AREA_CODE == "CSEO",])

#Look at relationship between hook variables, depth and CPUE (YE/hook)
unique(FT_YE$HOOK_SIZE)

ggplot(data=FT_YE, 
       aes(as.factor(HOOK_SIZE), YE_per_hook))+
  geom_boxplot(fill="grey", colour="black") +
  #facet_grid(HOOK_SIZE ~ ., scales = "free")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="grey")) +
  xlab("Year") +
  ylab("YE/hook") +
  labs(title="Raw YE/hook")

unique(FT_YE$HOOK_SPACING)

ggplot(data=FT_YE, 
         aes(as.factor(HOOK_SPACING), YE_per_hook))+
  geom_boxplot(fill="grey", colour="black") +
  #facet_grid(HOOK_SIZE ~ ., scales = "free")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="grey"))
str(FT_YE)
ggplot(data=FT_YE, 
       aes(AVG_DEPTH_FATHOMS, YE_per_hook))+
  geom_point() + #geom_line()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="grey"))

#annual trends in some variables... 
ggplot(data=FT_YE, 
       aes(as.factor(YEAR), AVG_DEPTH_FATHOMS))+
  geom_boxplot(fill="grey", colour="black") +
  facet_grid(G_MANAGEMENT_AREA_CODE ~ ., scales = "free")+ #geom_line()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="grey"))

## So, lets filter out depths > 150 fathoms and clean up hook size stuff
## standardize hook sizes...
with(FT_YE, table())

FT_YE %>% filter(AVG_DEPTH_FATHOMS < 150) %>%
  mutate(hook_size = as.numeric(gsub("[^0-9]", "", HOOK_SIZE)),
         std_hk_sz = ifelse (hook_size %in% c(3),16,
                             ifelse(hook_size %in% c(4),15,
                                    ifelse(hook_size %in% c(5),14,
                                           ifelse(hook_size %in% c(6),13,hook_size))))) -> FT_YE

ggplot(data=FT_YE, 
       aes(as.factor(std_hk_sz), YE_per_hook))+
  geom_boxplot(fill="grey", colour="black") +
  #facet_grid(HOOK_SIZE ~ ., scales = "free")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="grey")) +
  xlab("hook size") +
  ylab("YE/hook") 

#=================================================================================
# Examine and clean Logbook Data
str(LB_YE)
unique(LB_YE$Project.Code)
hist(LB_YE$Soak.Time.Hours, breaks=100)

#get rid of non-YE landings and super long soak times since they appear to be data/time data entry errors
# will need to get rid of 0 soak time stuff to examine.  may need to ignor because times
# not recorded in a lot of sets
LB_YE %>% filter (Species.Code == 145,Soak.Time.Hours <600,
                  Groundfish.Management.Area.Code != "WYAK") %>% 
  mutate(date_set = parse_date_time(Time.Set, c("%m/%d/%Y %H:%M")),
         date_haul = parse_date_time(Time.Hauled, c("%m/%d/%Y %H:%M")),
         haul_jdate = as.numeric(format(date_haul,"%j")),
         Date.Left.Port = parse_date_time(Date.Left.Port, c("%m/%d/%Y %H:%M")),
         Sell.Date = parse_date_time(Sell.Date, c("%m/%d/%Y %H:%M")))->LB_YE.1

#how many sets without soak time data? 
nrow(LB_YE.1[LB_YE.1$Soak.Time.Hours == 0,])/nrow(LB_YE.1)

str(LB_YE.1)

#how many NA's in Pounds
nrow(LB_YE.1[is.na(LB_YE.1$Pounds),])/nrow(LB_YE.1)
  # some most of the time they didn't weight fish (head slap)

#how many predation interactions
unique(LB_YE.1$Depredation) #not an issue
unique(LB_YE.1$Depredation.Code)
unique(LB_YE.1$Number.of.Skates.Impacted.by.Depredation)

#Effort numbers
hist(LB_YE.1$Effort.Number, breaks=25)

##
ggplot(data=LB_YE.1, 
       aes(as.factor(Year), Pounds))+
  geom_boxplot(fill="grey", colour="black") +
#  facet_grid(Year ~ ., scales = "free")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="grey"))

ggplot(data=LB_YE.1, 
       aes(as.factor(Year), Numbers))+
  geom_boxplot(fill="grey", colour="black") +
  #  facet_grid(Year ~ ., scales = "free")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="grey"))

#Lets look at Numbers... 
Nos<-nrow(LB_YE.1[is.na(LB_YE.1$Numbers),])/nrow(LB_YE.1)
Wts<-nrow(LB_YE.1[is.na(LB_YE.1$Pounds),])/nrow(LB_YE.1)
Nos+Wts
#well damn - they either entered data as numbers OR weights... 

hist(LB_YE.1$Pounds, breaks=100)
hist(LB_YE.1$Numbers, breaks=100)

#looks like either pounds/trip or numbers/trip
unique(LB_YE.1$Trip.Number)
LB_YE.1[LB_YE.1$Trip.Number == 6 & LB_YE.1$Year == 1986,]

LB_YE.1 %>% group_by(Year, Trip.Number) %>%
  mutate(Rec.unit = ifelse(is.na(Pounds),"count","weight"),
         Pnds_per_trip = sum(Pounds),
         Cnt_per_trip = sum(Numbers),
         Sets_per_trip = max(Effort.Number),
         Tot.soak_per_trip = sum(Soak.Time.Hours)) %>% 
  mutate(Pnds_per_set = sum(Pounds)/Sets_per_trip,
         Pnds_per_soak = sum(Pounds)/Tot.soak_per_trip,
         Cnt_per_set = sum(Numbers)/Sets_per_trip,
         Cnt_per_soak = sum(Numbers)/Tot.soak_per_trip)  %>% 
  ungroup() %>%
  group_by(Year) %>% 
  mutate(N_trips = sum(unique(Trip.Number)))->LB_YE.2

unique(LB_YE.2$Sets_per_trip)
unique(LB_YE.2$Cnt_per_trip)

plot(LB_YE.2$N_trips ~ LB_YE.2$Year)

Nos<-nrow(LB_YE.2[is.na(LB_YE.2$Cnt_per_trip),])/nrow(LB_YE.2)
Wts<-nrow(LB_YE.2[is.na(LB_YE.2$Pnds_per_trip),])/nrow(LB_YE.2)
Nos+Wts

#Effort plots
ggplot(data=LB_YE.2 %>% filter(Rec.unit == "weight"), 
       aes(as.factor(Year), Sets_per_trip))+
  geom_violin(fill="grey", colour="black") +
  geom_jitter(width=0.2, height=0, col="blue", size=0.7, alpha=0.15) +
  facet_grid(Groundfish.Management.Area.Code ~ ., scales = "free")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="grey"))

ggplot(data=LB_YE.2 %>% filter(Rec.unit == "weight"), 
       aes(as.factor(Year), Tot.soak_per_trip))+
  geom_violin(fill="grey", colour="black") +
  geom_jitter(width=0.2, height=0, col="green", size=0.7, alpha=0.15) +
  facet_grid(Groundfish.Management.Area.Code ~ ., scales = "free")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="grey"))

#weight plots
ggplot(data=LB_YE.2 %>% filter(Rec.unit == "weight"), 
       aes(as.factor(Year), Pnds_per_trip))+
  geom_violin(fill="grey", colour="black") +
  geom_jitter(width=0.2, height=0, col="blue", size=0.7, alpha=0.15) +
  facet_grid(Groundfish.Management.Area.Code ~ ., scales = "free")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="grey"))

ggplot(data=LB_YE.2 %>% filter(Rec.unit == "weight"), 
       aes(as.factor(Year), Pnds_per_set))+
  geom_violin(fill="grey", colour="black") +
  geom_jitter(width=0.2, height=0, col="red", size=0.7, alpha=0.15) +
  facet_grid(Groundfish.Management.Area.Code ~ ., scales = "free")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="grey"))

ggplot(data=LB_YE.2 %>% filter(Rec.unit == "weight", Pnds_per_soak >=0), 
       aes(as.factor(Year), Pnds_per_soak))+
  geom_violin(fill="grey", colour="black") +
  geom_jitter(width=0.2, height=0, col="violet", size=0.7, alpha=0.15) +
  facet_grid(Groundfish.Management.Area.Code ~ ., scales = "free")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="grey"))

View(LB_YE.2[LB_YE.2$Pnds_per_soak == -533.3333,])
hist(LB_YE.2$Pnds_per_soak)

#count plots
ggplot(data=LB_YE.2 %>% filter(Rec.unit == "count"), 
       aes(as.factor(Year), Cnt_per_trip))+
  geom_violin(fill="grey", colour="black") +
  geom_jitter(width=0.2, height=0, col="red", size=0.7, alpha=0.15) +
  facet_grid(Groundfish.Management.Area.Code ~ ., scales = "free")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="grey"))

ggplot(data=LB_YE.2 %>% filter(Rec.unit == "count"), 
       aes(as.factor(Year), Cnt_per_set))+
  geom_violin(fill="grey", colour="black") +
  geom_jitter(width=0.2, height=0, col="blue", size=0.7, alpha=0.15) +
  facet_grid(Groundfish.Management.Area.Code ~ ., scales = "free")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="grey"))

ggplot(data=LB_YE.2 %>% filter(Rec.unit == "count"), 
       aes(as.factor(Year), Cnt_per_soak))+
  geom_violin(fill="grey", colour="black") +
  geom_jitter(width=0.2, height=0, col="orange", size=0.7, alpha=0.15) +
  facet_grid(Groundfish.Management.Area.Code ~ ., scales = "free")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="grey"))

#===================================================================================
# 1) Fish Tickets: Pounds/Hook = gets rid of uncertainty about the number of sets
  # soak time can be ignored I think, but will examine Log book data
  # For Fish ticket data need to figure out effects of hook size & depth

# 2) For Log book data need to figure out effects of depth and soak time
  #if soak time is an issue, would need to merge logbook and fish ticket data to account for that

#1)...
str(FT_YE)
nrow(FT_YE[is.na(FT_YE$YE_per_hook),])/nrow(FT_YE)
nrow(FT_YE[is.na(FT_YE$YEAR),])/nrow(FT_YE)
nrow(FT_YE[is.na(FT_YE$Jdate),])/nrow(FT_YE)
nrow(FT_YE[is.na(FT_YE$std_hk_sz),])/nrow(FT_YE)
nrow(FT_YE[is.na(FT_YE$AVG_DEPTH_FATHOMS),])/nrow(FT_YE)
nrow(FT_YE[is.na(FT_YE$HOOK_SPACING),])/nrow(FT_YE)

Complete_FT_YE<-FT_YE[!is.na(FT_YE$Jdate) & !is.na(FT_YE$std_hk_sz) & 
                        !is.na(FT_YE$HOOK_SPACING) & !is.na(FT_YE$YE_per_hook) &
                        !is.infinite(FT_YE$YE_per_hook),]
nrow(Complete_FT_YE)
nrow(Complete_FT_YE[is.na(Complete_FT_YE$AVG_DEPTH_FATHOMS),])

Complete_FT_YE[is.infinite(Complete_FT_YE$YE_per_hook),]
Complete_FT_YE[is.infinite(Complete_FT_YE$HOOK_SPACING),]

Complete_FT_YE$YEAR<-as.factor(Complete_FT_YE$YEAR)

fsh_cpue<-Complete_FT_YE
fsh_cpue$subd<-as.factor(fsh_cpue$G_MANAGEMENT_AREA_CODE)
fsh_cpue$dumstat=1

m1 <- bam(YE_per_hook ~ YEAR + Jdate + std_hk_sz + 
            s(AVG_DEPTH_FATHOMS, k=4) + s(HOOK_SPACING, k=4) +
            subd,
          #s(subd, bs='re', by=dumstat),
          data=fsh_cpue, gamma=1.4)
m2 <- bam(YE_per_hook ~ YEAR + Jdate + #DROP:  std_hk_sz + 
            s(AVG_DEPTH_FATHOMS, k=4) + s(HOOK_SPACING, k=4)+
            subd,
          #s(subd, bs='re', by=dumstat),
          data=fsh_cpue, gamma=1.4)
m3 <- bam(YE_per_hook ~ YEAR  + std_hk_sz + #DROP + Jdate
            s(AVG_DEPTH_FATHOMS, k=4) + s(HOOK_SPACING, k=4) +
            subd,
            #s(subd, bs='re', by=dumstat),
          data=fsh_cpue, gamma=1.4)
m4 <- bam(YE_per_hook ~ YEAR + Jdate + std_hk_sz + 
             s(HOOK_SPACING, k=4) +
            subd,
            #s(subd, bs='re', by=dumstat), #Drop s(AVG_DEPTH_FATHOMS, k=4) +
          data=fsh_cpue, gamma=1.4)
m5 <- bam(YE_per_hook ~ YEAR + Jdate + std_hk_sz + 
            s(AVG_DEPTH_FATHOMS, k=4) +
            subd,
          #s(subd, bs='re', by=dumstat), # DROP + s(HOOK_SPACING, k=4),
          data=fsh_cpue, gamma=1.4)
summary(m1) 
summary(m2)
summary(m3)
summary(m4)
summary(m5)

AIC(m1, m2,m3, m4,m5)

#Hook spacing and depth significant, but Rsq not very exciting
#AIC likes Jdate and hook size, although not significant
# RE of subdistrict very significant

m6 <- bam(YE_per_hook ~ YEAR + # Drop both: Jdate + std_hk_sz + 
            s(AVG_DEPTH_FATHOMS, k=4) + s(HOOK_SPACING, k=4) + 
            subd,
            #s(subd, bs='re', by=dumstat),
          data=fsh_cpue, gamma=1.4)

summary(m6)

AIC(m1,m2,m3,m4,m5,m6)

plot(fitted(m1), resid(m1))
abline(h=0,col="red")

plot(m1, page=1, shade=TRUE, resid=TRUE, all=TRUE)
plot(m6, page=1, shade=TRUE, resid=TRUE, all=TRUE)
#YE_per_hook increases with depth and hook spacing
vis.gam(m6,c('AVG_DEPTH_FATHOMS','HOOK_SPACING'),plot.type='contour',
        type='response', color='topo',too.far=0.1)
vis.gam(m1,c('AVG_DEPTH_FATHOMS','HOOK_SPACING'),plot.type='contour',
        type='response', color='topo',too.far=0.1)


#predicted values to compare to nominal cpue
std_dat <- expand.grid(year = unique(fsh_cpue$YEAR),
                       AVG_DEPTH_FATHOMS = mean(fsh_cpue$AVG_DEPTH_FATHOMS), 
                       HOOK_SPACING = mean(fsh_cpue$HOOK_SPACING), 
                       subd = unique(fsh_cpue$subd),
                       #julian_day = median(Complete_FT_YE$Jdat),
                       dum = 0,
                       dumstat = 0) %>% 
  mutate(YEAR = factor(year))

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
         bt_cpue = exp(fit) - (mean(fsh_cpue$YE_per_hook) * 0.1),
         bt_upper = exp(upper) - (mean(fsh_cpue$YE_per_hook) * 0.1),
         bt_lower = exp(lower) - (mean(fsh_cpue$YE_per_hook) * 0.1),
         bt_se = (bt_upper - bt_cpue) / 2  #,
         #bt_cv = bt_se/bt_cpue
  ) -> std_dat

# Nominal CPUE ----

fsh_cpue %>% 
  group_by(YEAR,subd) %>% 
  dplyr::summarise(fsh_cpue = mean(YE_per_hook),
                   sd = sd(YE_per_hook),
                   n = length(YE_per_hook),
                   se = sd / (n ^ (1/2)),
                   var = var(YE_per_hook),
                   cv = sd / fsh_cpue,
                   upper = fsh_cpue + (2 * se),
                   lower = fsh_cpue - (2 * se)) -> fsh_sum 

# Compare predicted cpue from gam to nominal cpue
fsh_sum %>%
  select(YEAR, subd, cpue = fsh_cpue, upper, lower) %>% 
  mutate(CPUE = "Nominal") %>% 
  bind_rows(std_dat %>% 
              #select(YEAR, cpue = bt_cpue, upper = bt_upper, lower = bt_lower) %>% 
              select(YEAR, cpue = fit, upper , lower, subd ) %>%
              mutate(CPUE = "GAM")) -> comp
View(comp)
##!!! bt_cpue is awful high relative to nominal... something going on need to figure out... 

comp%>% 
  ggplot() +
  geom_ribbon(aes(YEAR, ymin = lower, ymax = upper, fill = CPUE), 
              colour = "white", alpha = 0.2) +
  geom_point(aes(YEAR, cpue, colour = CPUE, shape = CPUE), size = 2) +
  geom_line(aes(YEAR, cpue, colour = CPUE, group = CPUE), size = 1) +
  facet_grid(subd ~ ., scales = "free")+
  # scale_colour_grey(name = "Standardized CPUE") +
  # scale_fill_grey(name = "Standardized CPUE") +
  scale_colour_manual(values = c("darkcyan", "goldenrod"), name = "Standardized CPUE") +
  scale_fill_manual(values = c("darkcyan", "goldenrod"), name = "Standardized CPUE") +
  scale_shape_manual(values = c(19, 17), name = "Standardized CPUE") +
  #scale_x_continuous(breaks = axis$breaks, labels = axis$labels) + 
  labs(x = "", y = "Fishery CPUE (round lb/hook)\n") +
  theme(legend.position = c(0.8, 0.2)) +
  expand_limits(y = 0)



ggsave(paste0("figures/compare_stdcpue_llfsh_2022reboot_", YEAR, ".png"), dpi=300, height=4, width=7, units="in")

### SAVE csv of nominal CPUE... 
fsh_sum

#==============================
#2) Look at log bok stuff now
str(LB_YE.2)
lb_cpue<-LB_YE.2 %>% 
  select(Year, date_haul, Jdat = haul_jdate, subd=Groundfish.Management.Area.Code,
         depth = Average.Depth.Fathoms,
         soak.sum=Tot.soak_per_trip, 
         Pnds_per_set, Pnds_per_soak,
         Cnt_per_set, Cnt_per_soak)

nrow(lb_cpue[is.na(lb_cpue$depth) | is.infinite(lb_cpue$depth),])/nrow(lb_cpue)
nrow(lb_cpue[is.na(lb_cpue$soak.sum) | is.infinite(lb_cpue$soak.sum),])/nrow(lb_cpue)
nrow(lb_cpue[is.na(lb_cpue$Jdat) | is.infinite(lb_cpue$Jdat),])/nrow(lb_cpue)
nrow(lb_cpue[is.na(lb_cpue$Pnds_per_set) | is.infinite(lb_cpue$Pnds_per_set),])/nrow(lb_cpue)
nrow(lb_cpue[is.na(lb_cpue$Pnds_per_soak) | is.infinite(lb_cpue$Pnds_per_soak),])/nrow(lb_cpue)
nrow(lb_cpue[is.na(lb_cpue$Cnt_per_set) | is.infinite(lb_cpue$Cnt_per_set),])/nrow(lb_cpue)
nrow(lb_cpue[is.na(lb_cpue$Cnt_per_soak) | is.infinite(lb_cpue$Cnt_per_soak),])/nrow(lb_cpue)

lb_counts<-lb_cpue[!is.na(lb_cpue$depth) & !is.infinite(lb_cpue$depth) &
                     !is.na(lb_cpue$Cnt_per_set) & !is.infinite(lb_cpue$Cnt_per_set) &
                     !is.na(lb_cpue$Cnt_per_soak) & !is.infinite(lb_cpue$Cnt_per_soak),]
lb_pnds<-lb_cpue[!is.na(lb_cpue$depth) & !is.infinite(lb_cpue$depth) &
                   !is.na(lb_cpue$Pnds_per_set) & !is.infinite(lb_cpue$Pnds_per_set) &
                   !is.na(lb_cpue$Pnds_per_soak) & !is.infinite(lb_cpue$Pnds_per_soak),]
lb_counts$Year<-as.factor(lb_counts$Year)
lb_pnds$Year<-as.factor(lb_pnds$Year)
# number of YE caught per set... 
str(lb_counts)
m1 <- bam(Cnt_per_set ~ Year + Jdat + subd+ 
            s(depth, k=4) + s(soak.sum, k=4),
          data=lb_counts, gamma=1.4)
m2 <- bam(Cnt_per_set ~ Year + subd+ #drop Jdat
            s(depth, k=4) + s(soak.sum, k=4),
          data=lb_counts, gamma=1.4)
m3 <- bam(Cnt_per_set ~ Year + Jdat + subd+ 
            s(soak.sum, k=4),  #drop depth
          data=lb_counts, gamma=1.4)
m4 <- bam(Cnt_per_set ~ Year + Jdat + subd+ 
            s(depth, k=4) , #drop soak
          data=lb_counts, gamma=1.4)
m5 <- bam(Cnt_per_set ~ Year + Jdat + subd, #+ 
           # s(depth, k=4) + s(soak.sum, k=4),
          data=lb_counts, gamma=1.4)


summary(m1) 
summary(m2)
summary(m3)
summary(m4)
summary(m5)

AIC(m1, m2,m3, m4,m5)

#Hook spacing and depth significant, but Rsq not very exciting
#AIC likes Jdate and hook size, although not significant
# RE of subdistrict very significant

plot(fitted(m1), resid(m1))
abline(h=0,col="red")

plot(m1, page=1, shade=TRUE, resid=TRUE, all=TRUE)
plot(m6, page=1, shade=TRUE, resid=TRUE, all=TRUE)
#YE_per_hook increases with depth and hook spacing

vis.gam(m1,c('depth','soak.sum'),plot.type='contour',
        type='response', color='topo',too.far=0.1)
vis.gam(m1,c('depth','Jdat'),plot.type='contour',
        type='response', color='topo',too.far=0.1)
vis.gam(m1,c('soak.sum','Jdat'),plot.type='contour',
        type='response', color='topo',too.far=0.1)

#predicted values to compare to nominal cpue
std_dat <- expand.grid(Year = unique(lb_counts$Year),
                       depth = mean(lb_counts$depth), 
                       soak.sum = mean(lb_counts$soak.sum), 
                       subd = unique(lb_counts$subd),
                       Jdat = median(lb_counts$Jdat),
                       dum = 0,
                       dumstat = 0) %>% 
  mutate(Year = factor(Year))

pred_cpue <- predict(m1, std_dat, type = "link", se = TRUE)

#checking my code with Jane's... checks out :)
preds<-predict.bam(m1, type="response", std_dat, se = TRUE)
str(preds); head(preds)

#Put the standardized CPUE and SE into the data frame and convert to
#backtransformed (bt) CPUE

std_dat %>% 
  mutate(fit = pred_cpue$fit,
         se = pred_cpue$se.fit,
         upper = fit + (2 * se),
         lower = fit - (2 * se),
         bt_cpue = exp(fit) - (mean(fsh_cpue$YE_per_hook) * 0.1),
         bt_upper = exp(upper) - (mean(fsh_cpue$YE_per_hook) * 0.1),
         bt_lower = exp(lower) - (mean(fsh_cpue$YE_per_hook) * 0.1),
         bt_se = (bt_upper - bt_cpue) / 2  #,
         #bt_cv = bt_se/bt_cpue
  ) -> std_dat

# Nominal CPUE ----

lb_counts %>% 
  group_by(Year,subd) %>% 
  dplyr::summarise(fsh_cpue = mean(Cnt_per_set),
                   sd = sd(Cnt_per_set),
                   n = length(Cnt_per_set),
                   se = sd / (n ^ (1/2)),
                   var = var(Cnt_per_set),
                   cv = sd / fsh_cpue,
                   upper = fsh_cpue + (2 * se),
                   lower = fsh_cpue - (2 * se)) -> fsh_sum 

# Compare predicted cpue from gam to nominal cpue
fsh_sum %>%
  select(Year, subd, cpue = fsh_cpue, upper, lower) %>% 
  mutate(CPUE = "Nominal") %>% 
  bind_rows(std_dat %>% 
              #select(YEAR, cpue = bt_cpue, upper = bt_upper, lower = bt_lower) %>% 
              select(Year, cpue = fit, upper , lower, subd ) %>%
              mutate(CPUE = "GAM")) -> comp
View(comp)
##!!! bt_cpue is awful high relative to nominal... something going on need to figure out... 

comp%>% 
  ggplot() +
  geom_ribbon(aes(Year, ymin = lower, ymax = upper, fill = CPUE), 
              colour = "white", alpha = 0.2) +
  geom_point(aes(Year, cpue, colour = CPUE, shape = CPUE), size = 2) +
  geom_line(aes(Year, cpue, colour = CPUE, group = CPUE), size = 1) +
  facet_grid(subd ~ ., scales = "free")+
  # scale_colour_grey(name = "Standardized CPUE") +
  # scale_fill_grey(name = "Standardized CPUE") +
  scale_colour_manual(values = c("darkcyan", "goldenrod"), name = "Standardized CPUE") +
  scale_fill_manual(values = c("darkcyan", "goldenrod"), name = "Standardized CPUE") +
  scale_shape_manual(values = c(19, 17), name = "Standardized CPUE") +
  #scale_x_continuous(breaks = axis$breaks, labels = axis$labels) + 
  labs(x = "", y = "Fishery CPUE (# yelloweye/set)\n") +
  theme(legend.position = c(0.8, 0.2)) +
  expand_limits(y = 0)

#==================================================================================
# Should merge the log book and fish ticket data to do this right... (sigh)



###################################################################################
# scrap and notes
#=================================================================================
m1 <- bam(YE_per_hook ~ YEAR + Jdate + std_hk_sz + 
            s(AVG_DEPTH_FATHOMS, k=4) + s(HOOK_SPACING, k=4),
          data=Complete_FT_YE, gamma=1.4)
m1a <- bam(cpue ~ Year + Gear + s(depth, k=4) + s(soak, k=4) + #hook size removed
             s(Stat, bs='re', by=dumstat)+ s(Adfg, bs='re', by=dum), data=fsh_cpue_cl, gamma=1.4)
m2 <- bam(cpue ~ Year + Gear + Hook_size + s(depth, k=4) + s(soak, k=4) + 
            s(Adfg, bs='re', by=dum), data=fsh_cpue_cl, gamma=1.4)
m3 <- bam(cpue ~ Year + Gear + Hook_size + s(depth, k=4) + s(soak, k=4) +
            s(Stat, bs='re', by=dumstat), data=fsh_cpue_cl, gamma=1.4)
m4 <- bam(cpue ~ Year + Gear + Hook_size + s(depth, k=4) + s(soak, k=4), 
          data=fsh_cpue_cl, gamma=1.4)

summary(m1); summary(m1a) 
summary(m2)
summary(m3)
summary(m4)

AIC(m1,m1a, m2,m3, m4)




str(LB_YE.2)
str(FT_YE)
FT_YE %>% mutate(Year = YEAR,
                 Trip.Number = TRIP_NO)-> FT_YE

library(plyr)
Testjoin<-join_all(list(LB_YE.2,FT_YE,by=c("Year","Trip.Number")))




