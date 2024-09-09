################################################################################
## File for compiling simulation files
##
## 8/26/24 Phil Joy
################################################################################
## NOTES: The models need to be examined because there are occasional bad starting
##        values flagged as log(neg) warnings. Thy are coming from the B values
##        which are in the production model in the transformed parameter section.
##        I could not come up with a way to control how stan draws starting values
##        for them so it is an unfortunate feature of the current model. 
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


prefix3 <- "nope_"
suffix1 <- "_nope"
#look at results from PE models:

res_files2 <- str_remove(files_list[grep(paste0("^", prefix1), files_list)], prefix1)
#res_files2 <- str_remove(res_files2[grep(paste0("^", suffix1), res_files2)], suffix1)
#simstat_files <- str_remove(files_list, prefix2)
stat_files2 <- str_remove(res_files2[grep(paste0("^", prefix1), files_list)], "ko")
stat_files2 <- str_remove(stat_files2[grep(paste0("^", prefix1), files_list)], "_Best")
stat_files2 <- str_remove(stat_files2[grep(paste0("^", prefix1), files_list)], "_Best_Iest")
stat_files2 <- str_remove(stat_files2[grep(paste0("^", prefix1), files_list)], "_Iest")
stat_files2 <- str_remove(stat_files2[grep(paste0("^", prefix1), files_list)], "_B_I")
stat_files2 <- str_remove(stat_files2[grep(paste0("^", prefix1), files_list)], "_I")
stat_files2 <- unique(sub("^nope","",stat_files2))

unique(res_files2)

for (i in 1:length(res_files2)) { #i <- 49 29 #29
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
master3 <- master3[!is.na(master3$model),]
master3 <- master3 %>%
  mutate(prop_div = div_trans / (iterations * 3),
         model = as.factor(model),
         mod_name = ifelse(model %in% c("ko"),"Index only, OE = PE",#Index only, OE = PE",
                           ifelse(model %in% c("ko_estB"),"Index & Biomass, Index OE = PE, Biomass OE independent",#"Index & Biomass, Index OE = PE, Biomass OE independent",
                                  ifelse(model %in% c("ko_estB_estI"),"Index & Biomass, OE and PE independent",#"Index & Biomass, OE and PE independent",
                                         ifelse(model %in% c("ko_estI"),"Index only, OE independent of PE",
                                                ifelse(model %in% c("fit_nope_B_I"),"No PE, Index & Biomass",
                                                       ifelse(model %in% c("fit_nope_I"),"No PE, Index only",NA))))))) %>%
  distinct()


# convergence:
hist(master3$param_not_conv, breaks = 100)
with(master3, table(param_not_conv))
# get rid of 2 unconverged models
master3 <- master3 %>% filter(param_not_conv == 0)


library(wesanderson)
cols <- wes_palette("IsleofDogs1")

#--------------------------------------------------------------------------------
## Divergent transitions: 
master3 %>% group_by(mod_name) %>%
  dplyr::summarise(iters = n(),
                   prop_w_divs = sum(prop_div > 0)/iters,
                   mean_prop_div = sum(prop_div)/iters,
                   mean_r_bias1 = sum(rbias_1)/iters,
                   var_r_bias1 = var(rbias_1),
                   prop_r_in_90CI = sum(r_wi_bnds_1 == TRUE)/iters,
                   #mean_r_bias2 = sum(rbias_2)/iters,
                   #mean_r_bias3 = sum(rbias_3)/iters,
                   mean_K_bias1 = sum(Kbias_1)/iters,
                   var_K_bias1 = var(Kbias_1),
                   prop_K_in_90CI = sum(K_wi_bnds_1 == TRUE)/iters,
                   mean_q_bias1 = sum(qbias_1)/iters,
                   var_q_bias1 = var(qbias_1),
                   prop_qin_90CI = "Bug in code!", #sum(q_wi_bnds_1 == TRUE)/iters, #FLAG!!! Bug in code and not recorded. Fixed 9/6/24pj 
                   #mean_pe_bias1 = sum(pebias_1)/iters,
                   mean_abspe_bias1 = sum(abspebias_1)/iters,
                   var_abspe_bias1 = var(abspebias_1),
                   prop_abspe_in_90CI = sum(abspe_wi_bnds_1 == TRUE)/iters,
                   mean_MSY_bias1 = sum(MSYbias_1)/iters,
                   var_MSY_bias1 = var(MSYbias_1),
                   prop_MSY_in_90CI = sum(MSYS_wi_bnds_1 == TRUE)/iters,
                   mean_Bmsy_bias1 = sum(Bmsybias_1)/iters,
                   var_Bmsy_bias1 = var(Bmsybias_1),
                   prop_Bmsy_in_90CI = sum(Bmsy_wi_bnds_1 == TRUE)/iters,
                   mean_Fmsy_bias1 = sum(Fmsybias_1)/iters,
                   var_Fmsy_bias1 = var(Fmsybias_1),
                   prop_Fmsy_in_90CI = sum(Fmsy_wi_bnds_1 == TRUE)/iters,
                   mean_SS_bias1 = sum(SSbias_1)/iters,
                   var_SS_bias1 = var(SSbias_1),
                   prop_SS_in_90CI = sum(SS_wi_bnds_1 == TRUE)/iters) -> summary

View(summary)
write.csv(summary,"Production_models_STAN/Output/Sim_test_results_Sept24.csv")

ggplot(master3) +
  geom_boxplot(aes(mod_name,prop_div), fill="white",alpha = 1,
               show.legend = FALSE) 

ggplot(master3, aes(prop_div)) + 
  geom_histogram() +
  facet_wrap(~mod_name) +
  xlab("Proportion of iterations that are divergent")

ggplot(master3) + # %>% filter(model == "ko")) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  geom_point(aes(prop_div, rbias_1, col = harv_hist)) +
  geom_point(aes(prop_div, rbias_2, col = harv_hist)) +
  geom_point(aes(prop_div, rbias_3, col = harv_hist)) +
  facet_wrap(~mod_name) +
  geom_smooth(aes(x = prop_div, y = rbias_1, col = harv_hist, fill = harv_hist),
              method = "lm", formula = y ~ x) +
  #geom_smooth(aes(prop_div, rbias_1, col = harv_hist),
  #            method = 'loess', span = 1, se = FALSE) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in r") + xlab("Proportions of iterations with divergent transitions")

ggplot(master3 %>% filter(prop_div > 0)) + # %>% filter(model == "ko")) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  geom_point(aes(prop_div, rbias_1, col = harv_hist)) +
  geom_point(aes(prop_div, rbias_2, col = harv_hist)) +
  geom_point(aes(prop_div, rbias_3, col = harv_hist)) +
  facet_wrap(~mod_name) +
  geom_smooth(aes(x = prop_div, y = rbias_1, col = harv_hist, fill = harv_hist),
              method = "lm", formula = y ~ x) +
  #geom_smooth(aes(prop_div, rbias_1, col = harv_hist),
  #            method = 'loess', span = 1, se = FALSE) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in r") + xlab("Proportions of iterations with divergent transitions")

ggplot(master3) + # %>% filter(prop_div > 0)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  geom_point(aes(prop_div, Kbias_1, col = harv_hist)) +
  geom_point(aes(prop_div, Kbias_2, col = harv_hist)) +
  geom_point(aes(prop_div, Kbias_3, col = harv_hist)) +
  facet_wrap(~mod_name) +
  geom_smooth(aes(x = prop_div, y = Kbias_1, col = harv_hist, fill = harv_hist), 
              method = "lm", formula = y ~ x) +
  #geom_smooth(aes(prop_div, Kbias_1, col = harv_hist),
  #            method = 'loess', span = 1, se = FALSE) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in K") + xlab("Proportions of iterations with divergent transitions")  

ggplot(master3 ) + # %>% filter(prop_div > 0)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  geom_point(aes(prop_div, qbias_1, col = harv_hist)) +
  geom_point(aes(prop_div, qbias_2, col = harv_hist)) +
  geom_point(aes(prop_div, qbias_3, col = harv_hist)) +
  facet_wrap(~mod_name) +
  geom_smooth(aes(x = prop_div, y = qbias_1, col = harv_hist, fill = harv_hist),
              method = "lm", formula = y ~ x) +
  #geom_smooth(aes(prop_div, qbias_1, col = harv_hist),
  #            method = 'loess', span = 1, se = FALSE) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in q") + xlab("Proportions of iterations with divergent transitions")  

ggplot(master3 ) + # %>% filter(prop_div > 0)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  geom_point(aes(true_tot_cont1, prop_div, col = harv_hist)) +
  geom_point(aes(true_tot_cont2, prop_div, col = harv_hist)) +
  geom_point(aes(true_tot_cont3, prop_div, col = harv_hist)) +
  facet_wrap(~mod_name) +
  #geom_smooth(aes(x = true_tot_cont1, y = prop_div, col = harv_hist, fill = harv_hist),
  #            method = "lm", formula = y ~ x) +
  geom_smooth(aes(true_tot_cont1, prop_div, col = harv_hist, fill = harv_hist),
              method = 'loess', span = 1, se = TRUE) +
  ylab("Proportion of iterations with divergent transitions") + 
  xlab("Contrast in true biomass") 

ggplot(master3) + # %>% filter(prop_div > 0)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  geom_point(aes(true_reb_cont1, prop_div, col = harv_hist)) +
  geom_point(aes(true_reb_cont2, prop_div, col = harv_hist)) +
  geom_point(aes(true_reb_cont3, prop_div, col = harv_hist)) +
  facet_wrap(~mod_name) +
  #geom_smooth(aes(x = true_reb_cont1, y = prop_div), method = "lm", formula = y ~ x) +
  geom_smooth(aes(true_reb_cont1, prop_div, col = harv_hist, fill = harv_hist),
              method = 'loess', span = 1, se = TRUE) +
  ylab("Proportion of iterations with divergent transitions") + 
  xlab("Proportion that biomass has rebounded from low point")

ggplot(master3) + # %>% filter(prop_div > 0)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  geom_point(aes(true_reb_len1, prop_div, col = harv_hist)) +
  geom_point(aes(true_reb_len2, prop_div, col = harv_hist)) +
  geom_point(aes(true_reb_len3, prop_div, col = harv_hist)) +
  facet_wrap(~mod_name) +
  #geom_smooth(aes(x = true_reb_len1, y = prop_div), method = "lm", formula = y ~ x) +
  geom_smooth(aes(true_reb_len1, prop_div, col = harv_hist, fill = harv_hist),
              method = 'loess', span = 1, se = TRUE) +
  ylab("Proportion of iterations with divergent transitions") + 
  xlab("Number of years between low biomass and end of time series")

## r ##----------
#Contrast
cols <- wes_palette("Zissou1Continuous")

ggplot(master3) + # %>% filter(prop_div > 0)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  geom_point(aes(true_tot_cont1, rbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont2, rbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont3, rbias_3, col = harv_hist), alpha = 0.1) +
  facet_wrap(~mod_name) +
  geom_smooth(aes(true_tot_cont1, rbias_1, col = harv_hist, fill = harv_hist),
              method = 'loess', span = 1, se = TRUE) +
  ylab("Bias in r") + xlab("Contrast in true biomass") +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) 

#Rebound
ggplot(master3) + # %>% filter(prop_div > 0)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  geom_point(aes(true_reb_cont1, rbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont2, rbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont3, rbias_3, col = harv_hist), alpha = 0.1) +
  facet_wrap(~mod_name) +
  geom_smooth(aes(true_reb_cont1, rbias_1, col = harv_hist, fill = harv_hist),
              method = 'loess', span = 1, se = TRUE) +
  ylab("Bias in r") + xlab("Biomass rebound from low biomass")+
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) 

#Rebound length
ggplot(master3) + 
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  geom_point(aes(true_reb_len1, rbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len2, rbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len3, rbias_3, col = harv_hist), alpha = 0.1) +
  facet_wrap(~mod_name) +
  geom_smooth(aes(true_reb_len1, rbias_1, col = harv_hist, fill = harv_hist),
              method = 'loess', span = 1, se = TRUE) +
  ylab("Bias in r") + xlab("Time between low biomass and end of time series") +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) 



## K ##-------------
cols <- wes_palette("AsteroidCity2")
#Contrast
ggplot(master3) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  geom_point(aes(true_tot_cont1, Kbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont2, Kbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont3, Kbias_3, col = harv_hist), alpha = 0.1) +
  facet_wrap(~mod_name) +
  geom_smooth(aes(true_tot_cont1, Kbias_1, col = harv_hist, fill = harv_hist),
              method = 'loess', span = 1, se = TRUE) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in K") + xlab("Contrast in true biomass")

#Rebound
ggplot(master3) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  geom_point(aes(true_reb_cont1, Kbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont2, Kbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont3, Kbias_3, col = harv_hist), alpha = 0.1) +
  facet_wrap(~mod_name) +
  geom_smooth(aes(true_reb_cont1, Kbias_1, col = harv_hist, fill = harv_hist),
              method = 'loess', span = 1, se = TRUE) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in K") + xlab("Biomass rebound from low biomass")

#Rebound length
ggplot(master3) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  geom_point(aes(true_reb_len1, Kbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len2, Kbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len3, Kbias_3, col = harv_hist), alpha = 0.1) +
  facet_wrap(~mod_name) +
  geom_smooth(aes(true_reb_len1, Kbias_1, col = harv_hist, fill = harv_hist),
              method = 'loess', span = 1, se = TRUE) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in K") + xlab("Time between low biomass and end of time series")

## q ##-------------
cols <- wes_palette("IsleofDogs1")
#Contrast
ggplot(master3) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  geom_point(aes(true_tot_cont1, qbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont2, qbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont3, qbias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_tot_cont1, qbias_1, col = harv_hist, fill = harv_hist),
              method = 'loess', span = 1, se = TRUE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in q") + xlab("Contrast in true biomass")

#Rebound
ggplot(master3) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  geom_point(aes(true_reb_cont1, qbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont2, qbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont3, qbias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_cont1, qbias_1, col = harv_hist, fill = harv_hist),
              method = 'loess', span = 1, se = TRUE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in q") + xlab("Biomass rebound from low biomass")

#Rebound length
ggplot(master3) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  geom_point(aes(true_reb_len1, qbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len2, qbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len3, qbias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_len1, qbias_1, col = harv_hist, fill = harv_hist),
              method = 'loess', span = 1, se = TRUE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in q") + xlab("Time between low biomass and end of time series")

## BIAS IN Overall PE ##----------------------------------
#Contrast
unique(master3$model)
ggplot(master3 %>% filter(model != "fit_nope_B_I" & model != "fit_nope_I")) + 
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  facet_wrap(~mod_name) + 
  geom_violin(aes(as.factor(harv_hist), pebias_1))

#check without extreme values:
ggplot(master3 %>% filter(model != "fit_nope_B_I" & model != "fit_nope_I")) + 
  scale_color_manual(values = cols) +
  facet_wrap(~mod_name) + 
  geom_violin(aes(as.factor(harv_hist), pebias_1)) +
  geom_boxplot(aes(as.factor(harv_hist), pebias_1),width=.1) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) 


ggplot(master3 %>% filter(model != "fit_nope_B_I" & model != "fit_nope_I")) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  geom_point(aes(true_tot_cont1, pebias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont2, pebias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont3, pebias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_tot_cont1, pebias_1, col = harv_hist, fill = harv_hist),
              method = 'loess', span = 1, se = TRUE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in process error") + xlab("Contrast in true biomass") #+
#  scale_y_log10(labels = function(x) format(x, big.mark = ",",
#                                            scientific = FALSE),
#                breaks = 10^(0:10), minor_breaks = rep(1:9, 21)*(10^rep(0:10, each=9))) +
#  annotation_logticks(sides = "l") + 

#Rebound
ggplot(master3 %>% filter(model != "fit_nope_B_I" & model != "fit_nope_I")) +
  scale_color_manual(values = cols) + scale_fill_manual(values = cols) +
  geom_point(aes(true_reb_cont1, pebias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont2, pebias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont3, pebias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_cont1, pebias_1, col = harv_hist),
              method = 'loess', span = 1, se = FALSE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in process error") + xlab("Biomass rebound from low biomass")

#Rebound length
ggplot(master3 %>% filter(model != "fit_nope_B_I" & model != "fit_nope_I")) +
  scale_color_manual(values = cols) + scale_fill_manual(values = cols) +
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

ggplot(master3 %>% filter(model != "fit_nope_B_I" & model != "fit_nope_I")) + 
  scale_color_manual(values = cols) + scale_fill_manual(values = cols) +
  facet_wrap(~mod_name) + 
  geom_violin(aes(as.factor(harv_hist), abspebias_1)) +
  geom_boxplot(aes(as.factor(harv_hist), abspebias_1),width=.1) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) 

ggplot(master3 %>% filter(model != "fit_nope_B_I" & model != "fit_nope_I")) +
  scale_color_manual(values = cols) + scale_fill_manual(values = cols) +
  geom_point(aes(true_tot_cont1, abspebias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont2, abspebias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont3, abspebias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_tot_cont1, abspebias_1, col = harv_hist, fill = harv_hist),
              method = 'loess', span = 1, se = TRUE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in absolute process error (magnitude of process error)") + xlab("Contrast in true biomass") #+
#  scale_y_log10(labels = function(x) format(x, big.mark = ",",
#                                            scientific = FALSE),
#                breaks = 10^(0:10), minor_breaks = rep(1:9, 21)*(10^rep(0:10, each=9))) +
#  annotation_logticks(sides = "l") + 

#Rebound
ggplot(master3 %>% filter(model != "fit_nope_B_I" & model != "fit_nope_I")) +
  scale_color_manual(values = cols) + scale_fill_manual(values = cols) +
  geom_point(aes(true_reb_cont1, abspebias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont2, abspebias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont3, abspebias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_cont1, abspebias_1, col = harv_hist, fill = harv_hist),
              method = 'loess', span = 1, se = TRUE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in absolute process error (magnitude of process error)") + xlab("Biomass rebound from low biomass")

#Rebound length
ggplot(master3 %>% filter(model != "fit_nope_B_I" & model != "fit_nope_I")) +
  scale_color_manual(values = cols) + scale_fill_manual(values = cols) +
  geom_point(aes(true_reb_len1, abspebias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len2, abspebias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len3, abspebias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_len1, abspebias_1, col = harv_hist, fill = harv_hist),
              method = 'loess', span = 1, se = TRUE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in absolute process error (magnitude of process error)") + xlab("Time between low biomass and end of time series")


# Is there a relationship between magnitude of PE bias (absolute PE) and and r bias?
ggplot(master3 %>% filter(model != "fit_nope_B_I" & model != "fit_nope_I")) +
  geom_point(aes(abspebias_1, rbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(abspebias_2, rbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(abspebias_3, rbias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(abspebias_1, rbias_1, col = harv_hist, fill = harv_hist),
              method = 'loess', span = 1, se = TRUE) +
  facet_wrap(~mod_name) +
  ylab("Bias in r") + xlab("Bias in absolute process error (magnitude)")
## BRPS: MSY ## ----------------------------------------------------------------
#Contrast
ggplot(master3) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  geom_point(aes(true_tot_cont1, MSYbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont2, MSYbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont3, MSYbias_3, col = harv_hist), alpha = 0.1) +
  #geom_smooth(aes(x = true_tot_cont1, y = MSYbias_1, col = harv_hist, fill = harv_hist),
  #            method = "lm", formula = y ~ x) +
  geom_smooth(aes(true_tot_cont1, MSYbias_1, col = harv_hist, fill = harv_hist),
              method = 'loess', span = 1, se = TRUE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in MSY") + xlab("Contrast in true biomass")

#Rebound
ggplot(master3) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  geom_point(aes(true_reb_cont1, MSYbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont2, MSYbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont3, MSYbias_3, col = harv_hist), alpha = 0.1) +
  #geom_smooth(aes(x = true_reb_cont1, y = MSYbias_1, col = harv_hist, fill = harv_hist), method = "lm", formula = y ~ x) +
  geom_smooth(aes(true_reb_cont1, MSYbias_1, col = harv_hist, fill = harv_hist),
              method = 'loess', span = 1, se = TRUE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in MSY") + xlab("Biomass rebound from low biomass")

#Rebound length
ggplot(master3) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  geom_point(aes(true_reb_len1, MSYbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len2, MSYbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len3, MSYbias_3, col = harv_hist), alpha = 0.1) +
  #geom_smooth(aes(x = true_reb_len1, y = MSYbias_1, col = harv_hist, fill = harv_hist), method = "lm", formula = y ~ x) +
  geom_smooth(aes(true_reb_len1, MSYbias_1, col = harv_hist, fill = harv_hist),
              method = 'loess', span = 1, se = TRUE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in MSY") + xlab("Time between low biomass and end of time series")

## Bmsy ##
#Contrast
ggplot(master3) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  geom_point(aes(true_tot_cont1, Bmsybias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont2, Bmsybias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont3, Bmsybias_3, col = harv_hist), alpha = 0.1) +
  #geom_smooth(aes(x = true_tot_cont1, y = log(Bmsybias_1), col = harv_hist, fill = harv_hist), 
  #            method = "lm", formula = y ~ x) +
  geom_smooth(aes(true_tot_cont1, Bmsybias_1, col = harv_hist, fill = harv_hist),
              method = 'loess', span = 1, se = TRUE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in Bmsy") + xlab("Contrast in true biomass")

#Rebound
ggplot(master3) +
  scale_color_manual(values = cols) + scale_fill_manual(values = cols) +
  geom_point(aes(true_reb_cont1, Bmsybias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont2, Bmsybias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont3, Bmsybias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_cont1, Bmsybias_1, col = harv_hist, fill = harv_hist),
              method = 'loess', span = 1, se = TRUE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in Bmsy") + xlab("Biomass rebound from low biomass")

#Rebound length
ggplot(master3) +
  scale_color_manual(values = cols) + scale_fill_manual(values = cols) +
  geom_point(aes(true_reb_len1, Bmsybias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len2, Bmsybias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len3, Bmsybias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_len1, Bmsybias_1, col = harv_hist, fill = harv_hist),
              method = 'loess', span = 1, se = TRUE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in Bmsy") + xlab("Time between low biomass and end of time series")

## Fmsy ##
#Contrast
ggplot(master3) +
  scale_color_manual(values = cols) + scale_fill_manual(values = cols) +
  geom_point(aes(true_tot_cont1, Fmsybias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont2, Fmsybias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont3, Fmsybias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_tot_cont1, Fmsybias_1, col = harv_hist, fill = harv_hist),
              method = 'loess', span = 1, se = TRUE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in Fmsy") + xlab("Contrast in true biomass")

#Rebound
ggplot(master3) +
  scale_color_manual(values = cols) + scale_fill_manual(values = cols) +
  geom_point(aes(true_reb_cont1, Fmsybias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont2, Fmsybias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont3, Fmsybias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_cont1, Fmsybias_1, col = harv_hist, fill = harv_hist),
              method = 'loess', span = 1, se = TRUE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in Fmsy") + xlab("Biomass rebound from low biomass")

#Rebound length
ggplot(master3) +
  scale_color_manual(values = cols) + scale_fill_manual(values = cols) +
  geom_point(aes(true_reb_len1, Fmsybias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len2, Fmsybias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len3, Fmsybias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_len1, Fmsybias_1, col = harv_hist, fill = harv_hist),
              method = 'loess', span = 1, se = TRUE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in Fmsy") + xlab("Time between low biomass and end of time series")

## Stock status ##
#Contrast
ggplot(master3) +
  scale_color_manual(values = cols) + scale_fill_manual(values = cols) +
  geom_point(aes(true_tot_cont1, SSbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont2, SSbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_tot_cont3, SSbias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_tot_cont1, SSbias_1, col = harv_hist, fill = harv_hist),
              method = 'loess', span = 1, se = TRUE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in Stock Status (B/B0)") + xlab("Contrast in true biomass")

#Rebound
ggplot(master3) +
  scale_color_manual(values = cols) + scale_fill_manual(values = cols) +
  geom_point(aes(true_reb_cont1, SSbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont2, SSbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_cont3, SSbias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_cont1, SSbias_1, col = harv_hist, fill = harv_hist),
              method = 'loess', span = 1, se = TRUE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in Stock Status (B/B0)") + xlab("Biomass rebound from low biomass")

#Rebound length
ggplot(master3) +
  scale_color_manual(values = cols) + scale_fill_manual(values = cols) +
  geom_point(aes(true_reb_len1, SSbias_1, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len2, SSbias_2, col = harv_hist), alpha = 0.1) +
  geom_point(aes(true_reb_len3, SSbias_3, col = harv_hist), alpha = 0.1) +
  geom_smooth(aes(true_reb_len1, SSbias_1, col = harv_hist, fill = harv_hist),
              method = 'loess', span = 1, se = TRUE) +
  facet_wrap(~mod_name) +
  geom_hline(yintercept = 0,linetype = "dashed", color = "red", size = 1) +
  ylab("Bias in Stock Status (B/B0)") + xlab("Time between low biomass and end of time series")

#-------------------------------------------------------------------------------
# ANALYSIS: just mucking about here. Need to think about how
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


#################################################################################
#################################################################################
#################################################################################

#PEmods <- subset_list <- files_list[!grepl("_nope", files_list)]
#nonPE <- subset_list <- files_list[grepl("_nope", files_list)]

#look at results from PE models:

res_files <- str_remove(files_list[grep(paste0("^", prefix1), files_list)], prefix1)
#simstat_files <- str_remove(files_list, prefix2)
stat_files <- str_remove(res_files[grep(paste0("^", prefix1), PEmods)], "ko")
stat_files <- str_remove(stat_files[grep(paste0("^", prefix1), PEmods)], "_Best")
stat_files <- str_remove(stat_files[grep(paste0("^", prefix1), PEmods)], "_Best_Iest")
stat_files <- str_remove(stat_files[grep(paste0("^", prefix1), PEmods)], "_Iest")
stat_files <- str_remove(stat_files[grep(paste0("^", prefix1), nonPE)], "_B_I")
stat_files <- str_remove(stat_files2[grep(paste0("^", prefix1), nonPE)], "_I")


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
# process error. Skip down to that section ... LINE 550
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

#nonPE <- subset_list <- files_list[grepl("_nope", files_list)]
#NOTATION: ko, ko_Best, ko_Best_Iset and ko_Iest are normal PE models
#          nope_B_I and nope_I are with no process error  


















