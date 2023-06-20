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
library(ggplot2)
library(tidyr)
library(lubridate)
library(dplyr)
library(grid)
library(mgcv)

# 1) Load the Data
Fishtix<-read.csv("Data/Longline Hooks and Ticket Pounds.csv")
str(Fishtix)
Fishtix %>% filter(TARGET_SPECIES_CODE == 145,
                   G_MANAGEMENT_AREA_CODE != "IBS") %>% 
  mutate(Trip_Date = parse_date_time(DATE_LEFT_PORT, c("%m/%d/%Y %H:%M")),
                   Sell_Date = parse_date_time(SELL_DATE, c("%m/%d/%Y %H:%M")),
                   Jdate = as.numeric(format(Trip_Date,"%j")),
                   Year = as.factor(YEAR),
                   yrs = as.numeric(YEAR),
                   Subd = as.factor(G_MANAGEMENT_AREA_CODE),
                   Trip_no = TRIP_NO,
                   Depth = AVG_DEPTH_FATHOMS,     
                   Spacing = HOOK_SPACING,
                   No.Hooks=HOOKS,
                   lbYE = YELLOWEYE_ROCKFISH,
                   lbYE_per_hook = YELLOWEYE_ROCKFISH/HOOKS,
                   hook_size = as.numeric(gsub("[^0-9]", "", HOOK_SIZE)),
                   std_hk_sz = ifelse (hook_size %in% c(3),16,
                                       ifelse(hook_size %in% c(4),15,
                                              ifelse(hook_size %in% c(5),14,
                                                     ifelse(hook_size %in% c(6),13,hook_size))))) -> cpue

str(cpue)

cpue %>% 
  select(Year, yrs, Subd, Trip_Date, Sell_Date, Trip_no, 
         ADFG_NO, 
         Jdate, Depth,
         Spacing, No.Hooks, hook_size, std_hk_sz, 
         lbYE, lbYE_per_hook)->cpue

cpue <- cpue %>% filter (!is.na(lbYE_per_hook), !is.infinite(lbYE_per_hook))
nrow(cpue[is.na(cpue$lbYE_per_hook),])

#quick check of NA's and Infinites of key variables
nrow(cpue[is.na(cpue$lbYE_per_hook) | is.infinite(cpue$lbYE_per_hook),])
nrow(cpue[is.na(cpue$Jdate) | is.infinite(cpue$Jdate),])
nrow(cpue[is.na(cpue$Depth) | is.infinite(cpue$Depth),])
nrow(cpue[is.na(cpue$Spacing) | is.infinite(cpue$Spacing),])
nrow(cpue[is.na(cpue$No.Hooks) | is.infinite(cpue$No.Hooks),])
nrow(cpue[is.na(cpue$std_hk_sz) | is.infinite(cpue$std_hk_sz),])
nrow(cpue[is.na(cpue$hook_size) | is.infinite(cpue$hook_size),])
#Look at the data again...

ggplot(data=cpue, 
       aes(Year, lbYE_per_hook))+
  geom_boxplot(fill="lightblue", colour="blue") +
  facet_grid(Subd ~ ., scales = "free")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="lightblue2")) +
  xlab("Year") +
  ylab("lbs yelloweye per hook") +
  labs(title="Raw YE/hook")

ggplot(data=cpue, 
       aes(Year, lbYE))+
  geom_boxplot(fill="violet", colour="purple") +
  facet_grid(Subd ~ ., scales = "free")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="violet")) +
  xlab("Year") +
  ylab("lbs yelloweye") +
  labs(title="lbs of YE")

#Look at relationship between hook variables, depth and CPUE (YE/hook)
unique(cpue$HOOK_SIZE)

plot(cpue$lbYE_per_hook ~ cpue$std_hk_sz)
lines(lm(cpue$lbYE_per_hook ~ cpue$std_hk_sz))

ggplot(data=cpue, 
       aes(as.factor(std_hk_sz), lbYE_per_hook))+
  geom_boxplot(fill="grey", colour="black") +
  #facet_grid(HOOK_SIZE ~ ., scales = "free")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="grey")) +
  xlab("hook size") +
  ylab("YE/hook") 

plot(cpue$lbYE_per_hook ~ cpue$Spacing)
#lines(lm(cpue$lbYE_per_hook ~ cpue$Spacing))
ggplot(data=cpue, 
       aes(as.factor(Spacing), lbYE_per_hook))+
  geom_boxplot(fill="grey", colour="black") +
  #facet_grid(HOOK_SIZE ~ ., scales = "free")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="grey"))

ggplot(data=cpue, 
       aes(Depth, lbYE_per_hook))+
  geom_point() + #geom_line()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="grey"))

#annual trends in some variables... 
ggplot(data=cpue, 
       aes(as.factor(Year), Depth))+
  geom_boxplot(fill="grey", colour="black") +
  facet_grid(Subd ~ ., scales = "free")+ #geom_line()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="grey"))

# Cull out the few extreme depths...
cpue<-cpue[cpue$Depth < 150,]
head(cpue, 20)
str(cpue)
#===============================================================================
# NOMINAL CPUE for EACH AREA
#library(plyr)
cpue %>% 
  group_by(Year,Subd) %>% 
  dplyr::summarise(yrs = yrs,
                   fsh_cpue = mean(lbYE_per_hook),
                   sd = sd(lbYE_per_hook),
                   n = length(lbYE_per_hook),
                   se = sd / (n ^ (1/2)),
                   var = var(lbYE_per_hook),
                   cv = sd / fsh_cpue,
                   upper = fsh_cpue + (2 * se),
                   lower = fsh_cpue - (2 * se))  -> fsh_sum 

View(fsh_sum)
fsh_sum<-fsh_sum[complete.cases(fsh_sum),]
str(fsh_sum)

#dplyr bug leaving duplicate rows
fsh_sum2<-unique(fsh_sum)
View(fsh_sum2)

fsh_sum<-fsh_sum[order(fsh_sum$Subd, fsh_sum$Year),]
nYear<-length(unique(fsh_sum$Year))

fsh_sum2 %>% 
  ggplot() +
  geom_ribbon(aes(Year, ymin = lower, ymax = upper, fill = Subd), 
              colour = "white", alpha = 0.2, outline.type="full") +
  geom_point(aes(Year, fsh_cpue, colour = Subd, shape = Subd), size = 2) +
  geom_line(aes(Year, fsh_cpue, colour = Subd, group = Subd), size = 1) +
  facet_grid(Subd ~ ., scales = "free")+
  # scale_colour_grey(name = "Standardized CPUE") +
  # scale_fill_grey(name = "Standardized CPUE") +
  #scale_colour_manual(values = c("darkcyan", "goldenrod"), name = "Standardized CPUE") +
  #scale_fill_manual(values = c("darkcyan", "goldenrod"), name = "Standardized CPUE") +
  #scale_shape_manual(values = c(19, 17), name = "Standardized CPUE") +
  #scale_x_continuous(breaks = axis$breaks, labels = axis$labels) + 
  labs(x = "", y = "Fishery CPUE (round lb/hook)\n") +
  theme(legend.position = c(0.8, 0.8)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="grey")) +
  expand_limits(y = 0)

View(fsh_sum)
fsh_sum<-fsh_sum2
#before saving for model will be easier to insert NA values for missing years here
#would be tedious to do it in the data loading script for the SPM model

allyrs<-seq(min(fsh_sum$yrs), max(fsh_sum$yrs),1)
SD<-unique(fsh_sum$Subd)
#doing it witha clumsy loop...
cpue_model<-data.frame(); i<-1

for (s in SD){
  for (y in allyrs){
    Dat<-fsh_sum[fsh_sum$yrs == y & fsh_sum$Subd == s,]
    if (nrow(Dat) == 0) {
      cpue_model[i,"Year"]<-y
      cpue_model[i,"Subd"]<-s
      cpue_model[i,"yrs"]<-y
      cpue_model[i,"fsh_cpue"]<-NA
      cpue_model[i,5:8]<-NA
      cpue_model[i,"cv"]<-1
      cpue_model[i,10:11]<-NA
      i<-i+1
    } else {
      cpue_model[i,"Year"]<-y
      cpue_model[i,"Subd"]<-s
      cpue_model[i,3:11]<-Dat[,3:11]
      i<-i+1
    }
  }
}
colnames(cpue_model)<-colnames(fsh_sum)
View(cpue_model)

cpue_model %>% 
  ggplot() +
  geom_ribbon(aes(Year, ymin = lower, ymax = upper, fill = Subd), 
              colour = "white", alpha = 0.2, outline.type="full") +
  geom_point(aes(Year, fsh_cpue, colour = Subd, shape = Subd), size = 2) +
  geom_line(aes(Year, fsh_cpue, colour = Subd, group = Subd), size = 1) +
  facet_grid(Subd ~ ., scales = "free")+
  labs(x = "", y = "Fishery CPUE (round lb/hook)\n") +
  theme(legend.position = c(0.8, 0.8)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background=element_rect(fill="grey")) +
  expand_limits(y = 0)

write.csv(cpue_model, "Data/YE_dirfishery_nom_CPUE.csv")

#================================================================================
# Make sure we still feel good about ignoring depth, hook spacing and hook size... 

#get rid of NAs and Inf for exam
cpue.com<-cpue[complete.cases(cpue),]
nrow(cpue)
nrow(cpue.com)
colnames(cpue.com)

m1 <- bam(lbYE_per_hook ~ Year + Jdate + std_hk_sz + 
            s(Depth, k=4) + s(Spacing, k=4) +
            Subd,
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

