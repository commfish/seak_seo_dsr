################################################################################
## File for compiling simulation files
##
## 8/26/24 Phil Joy
################################################################################
## NOTES: The models need to be examined because there are occasional bad starting
##        values flagged as log(neg) warnings.I tried to change this mid-stream
##        by putting bounds on q but that didn't fix anything so need to reexamine
##        before running again. 
##    2:  Need to record process error ests and biases in sims. r is always over-
##        estimated. During the development of my original model I found similar
##        results and the more I constrained process error the lower r was and vice
##        versa. That may be the key to these models is finding the appropriate
##        constraint on process error. OR ignoring it all together?
##
################################################################################

library(stringr)
library(ggplot2)
library(dplyr)

# Gather up the simulation results and link them to the simulation statistics: 

files_list <- list.files("Production_models_STAN/Output/sim_res")

prefix1 <- "results_"
prefix2 <- "sim_stats_"
sufix1 <- "_nope"

PEmods <- subset_list <- files_list[!grepl("_nope", files_list)]
nonPE <- subset_list <- files_list[grepl("_nope", files_list)]

#look at results from PE models:

res_files <- str_remove(PEmods[grep(paste0("^", prefix1), PEmods)], prefix1)
#simstat_files <- str_remove(files_list, prefix2)
stat_files <- str_remove(res_files[grep(paste0("^", prefix1), PEmods)], "ko")
stat_files <- str_remove(stat_files[grep(paste0("^", prefix1), PEmods)], "_Best")
stat_files <- str_remove(stat_files[grep(paste0("^", prefix1), PEmods)], "_Best_Iest")
stat_files <- str_remove(stat_files[grep(paste0("^", prefix1), PEmods)], "_Iest")

stat_files <- unique(stat_files)

for (i in 1:length(res_files)) {
  res <- read.csv(paste0("Production_models_STAN/Output/sim_res/results_", 
                         res_files[i]))
  resfile <- paste0("Production_models_STAN/Output/sim_res/results_", 
                    res_files[i])
  resfile
  
  clean1 <- str_remove(resfile,"Production_models_STAN/Output/sim_res/results_ko")
  clean2 <- str_remove(clean1,"_Best")
  clean3 <- str_remove(clean2,"_Iest")
  clean3
  
  match_stats <- paste0("Production_models_STAN/Output/sim_res/sim_stats",
                        clean3)
  
  stats <- read.csv(match_stats)
  
  results <- merge(res,stats,by = c("iter","harv_hist","Hmax"), all.x = TRUE)
  str(results)
  
  if (i == 1) {
    master <- results
  } else {
    master <- rbind(master,results)
  }
}

unique(master$model)
unique(master$harv_hist)
unique(master$Hmax)

#-------------------------------------------------------------------------------
# These are the original 1st draft simulations.
# Not different from the extra sims I ran to examine the effect of dropping 
# process error. Skip down to that section ... LINE 554
#-------------------------------------------------------------------------------

# Examine the results graphically: 

str(master)
colnames(master)

master2 <- master[!is.na(master$model),]
master2 <- master2 %>%
  mutate(prop_div = div_trans / (iterations * 3),
         model = as.factor(model),
         mod_name = ifelse(model %in% c("ko"),"Index only, OE = PE",#Index only, OE = PE",
                           ifelse(model %in% c("ko_estB"),"Index & Biomass, Index OE = PE, Biomass OE independent",#"Index & Biomass, Index OE = PE, Biomass OE independent",
                                  ifelse(model %in% c("ko_estB_estI"),"Index & Biomass, OE and PE independent",#"Index & Biomass, OE and PE independent",
                                         ifelse(model %in% c("ko_estI"),"Index only, OE independent of PE",NA))))) #"Index only, OE independent of PE",NA)))))
str(master2)
unique(master2$mod_name)
max(master2$div_trans)
master2[master2$div_trans == max(master2$div_trans),]

library(wesanderson)
names(wes_palettes)

# convergence:
hist(master2$param_not_conv, breaks = 100)
with(master2, table(param_not_conv))
# get rid of 6 unconverged models
master2 <- master2 %>% filter(param_not_conv == 0)

cols <- wes_palette("IsleofDogs1")

## DIVERGENT TRANSITIONS: 

ggplot(master2) + # %>% filter(model == "ko")) +
  scale_color_manual(values = cols) +
  geom_point(aes(prop_div, rbias_1, col = harv_hist)) +
  geom_point(aes(prop_div, rbias_2, col = harv_hist)) +
  geom_point(aes(prop_div, rbias_3, col = harv_hist)) +
  facet_wrap(~mod_name) +
  geom_smooth(aes(prop_div, rbias_1, col = harv_hist),
              method = 'loess', span = 2, se = FALSE) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in r") + xlab("Proportions of iterations with divergent transitions")
  
  ggplot(master2) + # %>% filter(model == "ko")) +
    scale_color_manual(values = cols) +
    geom_point(aes(prop_div, Kbias_1, col = harv_hist)) +
    geom_point(aes(prop_div, Kbias_2, col = harv_hist)) +
    geom_point(aes(prop_div, Kbias_3, col = harv_hist)) +
    facet_wrap(~mod_name) +
    geom_smooth(aes(prop_div, Kbias_1, col = harv_hist),
                method = 'loess', span = 2, se = FALSE) +
    geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in K") + xlab("Proportions of iterations with divergent transitions")  
  
  ggplot(master2 ) + # %>% filter(model == "ko")) +
    scale_color_manual(values = cols) +
    geom_point(aes(prop_div, qbias_1, col = harv_hist)) +
    geom_point(aes(prop_div, qbias_2, col = harv_hist)) +
    geom_point(aes(prop_div, qbias_3, col = harv_hist)) +
    facet_wrap(~mod_name) +
    geom_smooth(aes(prop_div, qbias_1, col = harv_hist),
                method = 'loess', span = 2, se = FALSE) +
    geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
    ylab("Bias in q") + xlab("Proportions of iterations with divergent transitions")  

ggplot(master2 ) + # %>% filter(model == "ko")) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_tot_cont1, prop_div, col = harv_hist)) +
  geom_point(aes(true_tot_cont2, prop_div, col = harv_hist)) +
  geom_point(aes(true_tot_cont3, prop_div, col = harv_hist)) +
  facet_wrap(~mod_name) +
  ylab("Proportion of iterations with divergent transitions") + 
  xlab("Contrast in true biomass")

ggplot(master2) + # %>% filter(model == "ko")) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_cont1, prop_div, col = harv_hist)) +
  geom_point(aes(true_reb_cont2, prop_div, col = harv_hist)) +
  geom_point(aes(true_reb_cont3, prop_div, col = harv_hist)) +
  facet_wrap(~mod_name) +
  ylab("Proportion of iterations with divergent transitions") + 
  xlab("Proportion that biomass has rebounded from low point")

ggplot(master2) + # %>% filter(model == "ko")) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_len1, prop_div, col = harv_hist)) +
  geom_point(aes(true_reb_len2, prop_div, col = harv_hist)) +
  geom_point(aes(true_reb_len3, prop_div, col = harv_hist)) +
  facet_wrap(~mod_name) +
  ylab("Proportion of iterations with divergent transitions") + 
  xlab("Proportion that biomass has rebounded from low point")

## BIAS IN r ##----------
#Contrast
cols <- wes_palette("Zissou1Continuous"); cols

ggplot(master2) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_tot_cont1, rbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont2, rbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont3, rbias_3, col = harv_hist), alpha = 0.1) +
  facet_wrap(~mod_name) +
  geom_smooth(aes(true_tot_cont1, rbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  ylab("Bias in r") + xlab("Contrast in true biomass")

#Rebound
ggplot(master2) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_cont1, rbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont2, rbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont3, rbias_3, col = harv_hist), alpha = 0.1) +
  facet_wrap(~mod_name) +
  geom_smooth(aes(true_reb_cont1, rbias_1, col = harv_hist,),
              method = 'loess', span = 1, se = TRUE) +
  ylab("Bias in r") + xlab("Biomass rebound from low biomass")

#Rebound length
ggplot(master2) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_len1, rbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len2, rbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len3, rbias_3, col = harv_hist), alpha = 0.1) +
  facet_wrap(~mod_name) +
  geom_smooth(aes(true_reb_len1, rbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  ylab("Bias in r") + xlab("Time between low biomass and end of time series")



## BIAS IN K ##-------------
cols <- wes_palette("AsteroidCity2")
#Contrast
ggplot(master2) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_tot_cont1, Kbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont2, Kbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont3, Kbias_3, col = harv_hist), alpha = 0.1) +
  facet_wrap(~mod_name) +
  geom_smooth(aes(true_tot_cont1, Kbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in K") + xlab("Contrast in true biomass")

#Rebound
ggplot(master2) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_cont1, Kbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont2, Kbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont3, Kbias_3, col = harv_hist), alpha = 0.1) +
  facet_wrap(~mod_name) +
  geom_smooth(aes(true_reb_cont1, Kbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in K") + xlab("Biomass rebound from low biomass")

#Rebound length
ggplot(master2) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_len1, Kbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len2, Kbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len3, Kbias_3, col = harv_hist), alpha = 0.1) +
  facet_wrap(~mod_name) +
  geom_smooth(aes(true_reb_len1, Kbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in K") + xlab("Time between low biomass and end of time series")

## BIAS IN q ##-------------
cols <- wes_palette("IsleofDogs1")
#Contrast
ggplot(master2) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_tot_cont1, qbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont2, qbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont3, qbias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_tot_cont1, qbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in q") + xlab("Contrast in true biomass")

#Rebound
ggplot(master2) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_cont1, qbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont2, qbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont3, qbias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_cont1, qbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in q") + xlab("Biomass rebound from low biomass")

#Rebound length
ggplot(master2) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_len1, qbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len2, qbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len3, qbias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_len1, qbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in q") + xlab("Time between low biomass and end of time series")

## BIAS IN Overall PE ##----------------------------------
#Contrast

ggplot(master2) + 
  scale_color_manual(values = cols) +
  facet_wrap(~mod_name) + 
  geom_violin(aes(as.factor(harv_hist), pebias_1))

#check without extreme values:
ggplot(master2 %>% filter(abs(pebias_1) < 10)) + 
  scale_color_manual(values = cols) +
  facet_wrap(~mod_name) + 
  geom_violin(aes(as.factor(harv_hist), pebias_1)) +
  geom_boxplot(aes(as.factor(harv_hist), pebias_1),width=.1)


ggplot(master2) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_tot_cont1, pebias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont2, pebias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont3, pebias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_tot_cont1, pebias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in process error") + xlab("Contrast in true biomass") #+
#  scale_y_log10(labels = function(x) format(x, big.mark = ",",
#                                            scientific = FALSE),
#                breaks = 10^(0:10), minor_breaks = rep(1:9, 21)*(10^rep(0:10, each=9))) +
#  annotation_logticks(sides = "l") + 

#Rebound
ggplot(master2) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_cont1, pebias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont2, pebias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont3, pebias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_cont1, pebias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in process error") + xlab("Biomass rebound from low biomass")

#Rebound length
ggplot(master2) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_len1, pebias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len2, pebias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len3, pebias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_len1, pebias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in process error") + xlab("Time between low biomass and end of time series")

## BIAS IN Absolute PE; i.e. magnitude of process error ##----------------------
#Contrast

ggplot(master2) + 
  scale_color_manual(values = cols) +
  facet_wrap(~mod_name) + 
  geom_violin(aes(as.factor(harv_hist), abspebias_1)) +
  geom_boxplot(aes(as.factor(harv_hist), abspebias_1),width=.1)

ggplot(master2) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_tot_cont1, abspebias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont2, abspebias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont3, abspebias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_tot_cont1, abspebias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in absolute process error (magnitude of process error)") + xlab("Contrast in true biomass") #+
#  scale_y_log10(labels = function(x) format(x, big.mark = ",",
#                                            scientific = FALSE),
#                breaks = 10^(0:10), minor_breaks = rep(1:9, 21)*(10^rep(0:10, each=9))) +
#  annotation_logticks(sides = "l") + 

#Rebound
ggplot(master2) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_cont1, abspebias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont2, abspebias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont3, abspebias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_cont1, abspebias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in absolute process error (magnitude of process error)") + xlab("Biomass rebound from low biomass")

#Rebound length
ggplot(master2) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_len1, abspebias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len2, abspebias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len3, abspebias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_len1, abspebias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in absolute process error (magnitude of process error)") + xlab("Time between low biomass and end of time series")

## BRPS: MSY ## ----------------------------------------------------------------
#Contrast
ggplot(master2) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_tot_cont1, MSYbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont2, MSYbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont3, MSYbias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_tot_cont1, MSYbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in MSY") + xlab("Contrast in true biomass")

#Rebound
ggplot(master2) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_cont1, MSYbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont2, MSYbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont3, MSYbias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_cont1, MSYbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in MSY") + xlab("Biomass rebound from low biomass")

#Rebound length
ggplot(master2) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_len1, MSYbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len2, MSYbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len3, MSYbias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_len1, MSYbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in MSY") + xlab("Time between low biomass and end of time series")

## Bmsy ##
#Contrast
ggplot(master2) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_tot_cont1, Bmsybias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont2, Bmsybias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont3, Bmsybias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_tot_cont1, Bmsybias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in Bmsy") + xlab("Contrast in true biomass")

#Rebound
ggplot(master2) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_cont1, Bmsybias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont2, Bmsybias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont3, Bmsybias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_cont1, Bmsybias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in Bmsy") + xlab("Biomass rebound from low biomass")

#Rebound length
ggplot(master2) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_len1, Bmsybias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len2, Bmsybias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len3, Bmsybias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_len1, Bmsybias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in Bmsy") + xlab("Time between low biomass and end of time series")

## Fmsy ##
#Contrast
ggplot(master2) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_tot_cont1, Fmsybias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont2, Fmsybias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont3, Fmsybias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_tot_cont1, Fmsybias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in Fmsy") + xlab("Contrast in true biomass")

#Rebound
ggplot(master2) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_cont1, Fmsybias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont2, Fmsybias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont3, Fmsybias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_cont1, Fmsybias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in Fmsy") + xlab("Biomass rebound from low biomass")

#Rebound length
ggplot(master2) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_len1, Fmsybias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len2, Fmsybias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len3, Fmsybias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_len1, Fmsybias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in Fmsy") + xlab("Time between low biomass and end of time series")

## Stock status ##
#Contrast
ggplot(master2) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_tot_cont1, SSbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont2, SSbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont3, SSbias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_tot_cont1, SSbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in Stock Status (B/B0)") + xlab("Contrast in true biomass")

#Rebound
ggplot(master2) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_cont1, SSbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont2, SSbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont3, SSbias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_cont1, SSbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in Stock Status (B/B0)") + xlab("Biomass rebound from low biomass")

#Rebound length
ggplot(master2) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_len1, SSbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len2, SSbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len3, SSbias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_len1, SSbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in Stock Status (B/B0)") + xlab("Time between low biomass and end of time series")

#-------------------------------------------------------------------------------
# Is there a relationship between magnitude of PE bias (absolute PE) and and r bias?
ggplot(master2) +
  geom_point(aes(abspebias_1, rbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(abspebias_2, rbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(abspebias_3, rbias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(abspebias_1, rbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  ylab("Bias in r") + xlab("Bias in absolute process error (magnitude)")




################################################################################
#-------------------------------------------------------------------------------
# What happens when we ignore process error?
#-------------------------------------------------------------------------------
################################################################################

# Gather up the simulation results and link them to the simulation statistics:

nonPE <- subset_list <- files_list[grepl("_nope", files_list)]
#NOTATION: ko, ko_Best, ko_Best_Iset and ko_Iest are normal PE models
#          nope_B_I and nope_I are with no process error  
prefix3 <- "nope_"
suffix1 <- "_nope"
#look at results from PE models:

res_files2 <- str_remove(nonPE[grep(paste0("^", prefix1), nonPE)], prefix1)
#res_files2 <- str_remove(res_files2[grep(paste0("^", suffix1), res_files2)], suffix1)
#simstat_files <- str_remove(files_list, prefix2)
stat_files2 <- str_remove(res_files2[grep(paste0("^", prefix1), nonPE)], "ko")
stat_files2 <- str_remove(stat_files2[grep(paste0("^", prefix1), nonPE)], "_Best")
stat_files2 <- str_remove(stat_files2[grep(paste0("^", prefix1), nonPE)], "_Best_Iest")
stat_files2 <- str_remove(stat_files2[grep(paste0("^", prefix1), nonPE)], "_Iest")
stat_files2 <- str_remove(stat_files2[grep(paste0("^", prefix1), nonPE)], "_B_I")
stat_files2 <- str_remove(stat_files2[grep(paste0("^", prefix1), nonPE)], "_I")
stat_files2 <- unique(sub("^nope","",stat_files2))

unique(res_files2)

for (i in 1:length(res_files2)) { #i <- 29 #29
  res <- read.csv(paste0("Production_models_STAN/Output/sim_res/results_", 
                         res_files2[i]))
  resfile <- paste0("Production_models_STAN/Output/sim_res/results_", 
                    res_files2[i])
  resfile
  
  clean1 <- str_remove(resfile,"Production_models_STAN/Output/sim_res/results")
  clean1aa <- sub("^_nope","",clean1) #str_remove(clean1,"_nope")
  clean1a <- str_remove(clean1aa,"_ko")
  clean2 <- str_remove(clean1a,"_Best")
  clean3 <- str_remove(clean2,"_Iest")
  clean4 <- str_remove(clean3,"_B")
  clean5 <- str_remove(clean4,"_I")
  clean5
  
  match_stats <- paste0("Production_models_STAN/Output/sim_res/sim_stats",
                        clean5)
  
  stats <- read.csv(match_stats)
  
  results <- merge(res,stats,by = c("iter","harv_hist","Hmax"), all.x = TRUE)
  str(results)
  
  if (i == 1) {
    master3 <- results
  } else {
    master3 <- rbind(master3,results)
  }
}

unique(master3$model)
unique(master3$harv_hist)
unique(master3$Hmax)

#-------------------------------------------------------------------------------
str(master)
colnames(master)

master3 <- master3[!is.na(master3$model),]
master3 <- master3 %>%
  mutate(prop_div = div_trans / (iterations * 3),
         model = as.factor(model),
         mod_name = ifelse(model %in% c("ko"),"Index only, OE = PE",#Index only, OE = PE",
                           ifelse(model %in% c("ko_estB"),"Index & Biomass, Index OE = PE, Biomass OE independent",#"Index & Biomass, Index OE = PE, Biomass OE independent",
                                  ifelse(model %in% c("ko_estB_estI"),"Index & Biomass, OE and PE independent",#"Index & Biomass, OE and PE independent",
                                         ifelse(model %in% c("ko_estI"),"Index only, OE independent of PE",
                                                ifelse(model %in% c("fit_nope_B_I"),"Index & Biomass, NO process error",
                                                       ifelse(model %in% c("fit_nope_I"),"Index only, NO process error",NA))))))) %>%
  distinct()

#Damnit... ran nope_B_I when I thought I was running nope_I. Ignore those results for now!
# Fix code when rerun!!!!!!!!

master3 <- master3 %>% filter(model != "fit_nope_I")

# convergence:
hist(master3$param_not_conv, breaks = 100)
with(master3, table(param_not_conv))
# get rid of 6 unconverged models
master3 <- master3 %>% filter(param_not_conv == 0)

cols <- wes_palette("IsleofDogs1")

## Divergent transitions: 

ggplot(master3) + # %>% filter(model == "ko")) +
  scale_color_manual(values = cols) +
  geom_point(aes(prop_div, rbias_1, col = harv_hist)) +
  geom_point(aes(prop_div, rbias_2, col = harv_hist)) +
  geom_point(aes(prop_div, rbias_3, col = harv_hist)) +
  facet_wrap(~mod_name) +
  #geom_smooth(aes(prop_div, rbias_1, col = harv_hist),
  #            method = 'loess', span = 1, se = FALSE) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in r") + xlab("Proportions of iterations with divergent transitions")

ggplot(master3) + # %>% filter(model == "ko")) +
  scale_color_manual(values = cols) +
  geom_point(aes(prop_div, Kbias_1, col = harv_hist)) +
  geom_point(aes(prop_div, Kbias_2, col = harv_hist)) +
  geom_point(aes(prop_div, Kbias_3, col = harv_hist)) +
  facet_wrap(~mod_name) +
  #geom_smooth(aes(prop_div, Kbias_1, col = harv_hist),
  #            method = 'loess', span = 1, se = FALSE) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in K") + xlab("Proportions of iterations with divergent transitions")  

ggplot(master3 ) + # %>% filter(model == "ko")) +
  scale_color_manual(values = cols) +
  geom_point(aes(prop_div, qbias_1, col = harv_hist)) +
  geom_point(aes(prop_div, qbias_2, col = harv_hist)) +
  geom_point(aes(prop_div, qbias_3, col = harv_hist)) +
  facet_wrap(~mod_name) +
  #geom_smooth(aes(prop_div, qbias_1, col = harv_hist),
  #            method = 'loess', span = 1, se = FALSE) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in q") + xlab("Proportions of iterations with divergent transitions")  

ggplot(master3 ) + # %>% filter(model == "ko")) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_tot_cont1, prop_div, col = harv_hist)) +
  geom_point(aes(true_tot_cont2, prop_div, col = harv_hist)) +
  geom_point(aes(true_tot_cont3, prop_div, col = harv_hist)) +
  facet_wrap(~mod_name) +
  ylab("Proportion of iterations with divergent transitions") + 
  xlab("Contrast in true biomass")

ggplot(master3) + # %>% filter(model == "ko")) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_cont1, prop_div, col = harv_hist)) +
  geom_point(aes(true_reb_cont2, prop_div, col = harv_hist)) +
  geom_point(aes(true_reb_cont3, prop_div, col = harv_hist)) +
  facet_wrap(~mod_name) +
  ylab("Proportion of iterations with divergent transitions") + 
  xlab("Proportion that biomass has rebounded from low point")

ggplot(master3) + # %>% filter(model == "ko")) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_len1, prop_div, col = harv_hist)) +
  geom_point(aes(true_reb_len2, prop_div, col = harv_hist)) +
  geom_point(aes(true_reb_len3, prop_div, col = harv_hist)) +
  facet_wrap(~mod_name) +
  ylab("Proportion of iterations with divergent transitions") + 
  xlab("Proportion that biomass has rebounded from low point")

## r ##----------
#Contrast
cols <- wes_palette("Zissou1Continuous")

ggplot(master3) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_tot_cont1, rbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont2, rbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont3, rbias_3, col = harv_hist), alpha = 0.1) +
  facet_wrap(~mod_name) +
  geom_smooth(aes(true_tot_cont1, rbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  ylab("Bias in r") + xlab("Contrast in true biomass") +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) 

#Rebound
ggplot(master3) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_cont1, rbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont2, rbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont3, rbias_3, col = harv_hist), alpha = 0.1) +
  facet_wrap(~mod_name) +
  geom_smooth(aes(true_reb_cont1, rbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  ylab("Bias in r") + xlab("Biomass rebound from low biomass")+
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) 

#Rebound length
ggplot(master3) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_len1, rbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len2, rbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len3, rbias_3, col = harv_hist), alpha = 0.1) +
  facet_wrap(~mod_name) +
  geom_smooth(aes(true_reb_len1, rbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  ylab("Bias in r") + xlab("Time between low biomass and end of time series") +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) 



## K ##-------------
cols <- wes_palette("AsteroidCity2")
#Contrast
ggplot(master3) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_tot_cont1, Kbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont2, Kbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont3, Kbias_3, col = harv_hist), alpha = 0.1) +
  facet_wrap(~mod_name) +
  geom_smooth(aes(true_tot_cont1, Kbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in K") + xlab("Contrast in true biomass")

#Rebound
ggplot(master3) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_cont1, Kbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont2, Kbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont3, Kbias_3, col = harv_hist), alpha = 0.1) +
  facet_wrap(~mod_name) +
  geom_smooth(aes(true_reb_cont1, Kbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in K") + xlab("Biomass rebound from low biomass")

#Rebound length
ggplot(master3) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_len1, Kbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len2, Kbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len3, Kbias_3, col = harv_hist), alpha = 0.1) +
  facet_wrap(~mod_name) +
  geom_smooth(aes(true_reb_len1, Kbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in K") + xlab("Time between low biomass and end of time series")

## q ##-------------
cols <- wes_palette("IsleofDogs1")
#Contrast
ggplot(master3) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_tot_cont1, qbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont2, qbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont3, qbias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_tot_cont1, qbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in q") + xlab("Contrast in true biomass")

#Rebound
ggplot(master3) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_cont1, qbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont2, qbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont3, qbias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_cont1, qbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in q") + xlab("Biomass rebound from low biomass")

#Rebound length
ggplot(master3) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_len1, qbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len2, qbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len3, qbias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_len1, qbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in q") + xlab("Time between low biomass and end of time series")


## BRPS: MSY ## ----------------------------------------------------------------
#Contrast
ggplot(master3) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_tot_cont1, MSYbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont2, MSYbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont3, MSYbias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_tot_cont1, MSYbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in MSY") + xlab("Contrast in true biomass")

#Rebound
ggplot(master3) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_cont1, MSYbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont2, MSYbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont3, MSYbias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_cont1, MSYbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in MSY") + xlab("Biomass rebound from low biomass")

#Rebound length
ggplot(master3) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_len1, MSYbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len2, MSYbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len3, MSYbias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_len1, MSYbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in MSY") + xlab("Time between low biomass and end of time series")

## Bmsy ##
#Contrast
ggplot(master3) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_tot_cont1, Bmsybias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont2, Bmsybias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont3, Bmsybias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_tot_cont1, Bmsybias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in Bmsy") + xlab("Contrast in true biomass")

#Rebound
ggplot(master3) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_cont1, Bmsybias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont2, Bmsybias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont3, Bmsybias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_cont1, Bmsybias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in Bmsy") + xlab("Biomass rebound from low biomass")

#Rebound length
ggplot(master3) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_len1, Bmsybias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len2, Bmsybias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len3, Bmsybias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_len1, Bmsybias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in Bmsy") + xlab("Time between low biomass and end of time series")

## Fmsy ##
#Contrast
ggplot(master3) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_tot_cont1, Fmsybias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont2, Fmsybias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont3, Fmsybias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_tot_cont1, Fmsybias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in Fmsy") + xlab("Contrast in true biomass")

#Rebound
ggplot(master3) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_cont1, Fmsybias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont2, Fmsybias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont3, Fmsybias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_cont1, Fmsybias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in Fmsy") + xlab("Biomass rebound from low biomass")

#Rebound length
ggplot(master3) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_len1, Fmsybias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len2, Fmsybias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len3, Fmsybias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_len1, Fmsybias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in Fmsy") + xlab("Time between low biomass and end of time series")

## Stock status ##
#Contrast
ggplot(master3) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_tot_cont1, SSbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont2, SSbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont3, SSbias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_tot_cont1, SSbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in Stock Status (B/B0)") + xlab("Contrast in true biomass")

#Rebound
ggplot(master3) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_cont1, SSbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont2, SSbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont3, SSbias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_cont1, SSbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in Stock Status (B/B0)") + xlab("Biomass rebound from low biomass")

#Rebound length
ggplot(master3) +
  scale_color_manual(values = cols) +
  geom_point(aes(true_reb_len1, SSbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len2, SSbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len3, SSbias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_len1, SSbias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in Stock Status (B/B0)") + xlab("Time between low biomass and end of time series")

#-------------------------------------------------------------------------------
# ANALYSIS: Very, very preliminary look at the results. Need to think about how
# to do this properly... 

# of course the strata measurements are highly correlated, so for now will use
# averages: 
master3 <- master3 %>% 
  mutate(contrast_mean = rowMeans(select(.,true_tot_cont1,true_tot_cont2,true_tot_cont3)),
         rebound_mean = rowMeans(select(.,true_reb_cont1,true_reb_cont2,true_reb_cont3)),
         reb_length_mean = rowMeans(select(.,true_reb_len1,true_reb_len2,true_reb_len3)),
         depletion_mean = rowMeans(select(.,Depl_1,Depl_2,Depl_3)),
         rbias_mean = rowMeans(select(.,rbias_1,rbias_2,rbias_3)),
         Kbias_mean = rowMeans(select(.,Kbias_1,Kbias_2,Kbias_3)),
         qbias_mean = rowMeans(select(.,qbias_1,qbias_2,qbias_3)),
         pebias_mean = rowMeans(select(.,pebias_1,pebias_2,pebias_3)),
         abspebias_mean = rowMeans(select(.,abspebias_1,abspebias_2,abspebias_3)),
         MSYbias_mean = rowMeans(select(.,MSYbias_1,MSYbias_2,MSYbias_3)),
         Bmsybias_mean = rowMeans(select(.,Bmsybias_1,Bmsybias_2,Bmsybias_3)),
         Fmsybias_mean = rowMeans(select(.,Fmsybias_1,Fmsybias_2,Fmsybias_3)))

vars <- master3 %>% select(contrast_mean, rebound_mean, reb_length_mean,
                           depletion_mean)
cor(vars, use = "complete.obs") 
#depletion and contrast are of course correlated so just look at contrast

# Divergences
hist(master3$prop_div, breaks = 50)

library(tweedie)
library(mgcv)

div_exam <- gam(prop_div ~ harv_hist + model + 
                  s(contrast_mean, k= 4) + s(rebound_mean, k=4) + 
                  s(reb_length_mean, k=4),
                data = master3, family=tw(), method = "REML")

summary(div_exam)
plot(div_exam, page = 1)

#R bias
hist(master3$rbias_mean, breaks = 50)

r_exam <- glm(rbias_mean ~ harv_hist + model + prop_div +
                contrast_mean + rebound_mean + 
                reb_length_mean,
              data = master3)

summary(r_exam); plot(r_exam)

#Kbias
hist(master3$Kbias_mean, breaks = 50)

K_exam <- glm(Kbias_mean ~ harv_hist + model + prop_div +
                contrast_mean + rebound_mean + 
                reb_length_mean,
              data = master3)

summary(K_exam); plot(K_exam)

#qbias
hist(master3$qbias_mean, breaks = 50)

q_exam <- glm(qbias_mean ~ harv_hist + model + prop_div +
                contrast_mean + rebound_mean + 
                reb_length_mean +
                model * contrast_mean +
                model * rebound_mean +
                model * reb_length_mean,
              data = master3)

summary(q_exam); plot(q_exam)

#pebias

#MSY bias
hist(master3$MSYbias_mean, breaks = 50)

MSY_exam <- glm(MSYbias_mean ~ harv_hist + model + prop_div +
                contrast_mean + rebound_mean + 
                reb_length_mean, #+
                #model * contrast_mean +
                #model * rebound_mean +
                #model * reb_length_mean,
              data = master3)

summary(MSY_exam); plot(MSY_exam)

#Bmsy bias
hist(master3$Bmsybias_mean, breaks = 50)

Bmsy_exam <- glm(Bmsybias_mean ~ harv_hist + model + prop_div +
                  contrast_mean + rebound_mean + 
                  reb_length_mean, #+
                #model * contrast_mean +
                #model * rebound_mean +
                #model * reb_length_mean,
                data = master3)

summary(Bmsy_exam); plot(Bmsy_exam)

#Fmsy bias
hist(master3$Fmsybias_mean, breaks = 50)

Fmsy_exam <- glm(Fmsybias_mean ~ harv_hist + model + prop_div +
                   contrast_mean + rebound_mean + 
                   reb_length_mean, #+
                 #model * contrast_mean +
                 #model * rebound_mean +
                 #model * reb_length_mean,
                 data = master3)

summary(Fmsy_exam); plot(Fmsy_exam)





















