################################################################################
## REMA for SEO Yelloweye rockfish
## from Jane Sullivan's package "rema" with guidance from Jane :)
## Rerunnig for Nov. 2022 plan team meeting
################################################################################
#If REMA package need to be loaded, here it is: 
# install.packages("devtools")
#devtools::install_github("afsc-assessments/rema", dependencies = TRUE, build_vignettes = TRUE)

#From Jane Sullivan 9-1-2023:
#If you estimate additional observation errors in your model, good news! You can now visualize how the additional estimated observation error relates to the assumed observation error based on the design-based estimator.   

#Documentation:
#  >The original issue is here: https://github.com/afsc-assessments/rema/issues/14
#>A reproducible example is in this section of the "Fitting to an additional CPUE survey" vignette: https://afsc-assessments.github.io/rema/articles/ex2_cpue.html#model-2-estimating-additional-observation-error-in-the-two-surveys
#>The new functions are "tidy_extra_cv()" and "plot_extra_cv()"

#Important notes:
#  >If you want access to this new functionality, you have to update rema! 
#  devtools::install_github("afsc-assessments/rema", dependencies = TRUE, build_vignettes = TRUE) **add force=TRUE if that doesn't work.
#>This will not affect any other functionality of the rema package. Don't worry about your existing code breaking.


# Example R scripts are downloaded when `rema` is installed. Locate them on your computer by running the following commands:
(rema_path <- find.package('rema'))
(rema_examples <- file.path(rema_path, 'example_scripts'))
list.files(rema_examples)

library("rema")
library("scales")
library("viridis")
library("dplyr")
library("ggplot2")
library("ggpubr")
library("kableExtra")

YEAR<-2024

cbpalette <- c("#009E73", "#0072B2", "#E69F00", "#56B4E9", "#D55E00", "#CC79A7","#F0E442", "black", "grey")

#ROV biomass .. check date tag to get most recent
bio<-read.csv("Data_processing/Data/SEO_YE_Biomass_subdistrict_2024-10-28.csv")
#str(bio_new)
str(bio)

bio<-bio %>% select(strata = Subdistrict, year = Year, biomass = Biomass.mt,
                    cv = Biomass.cv)

#IPHC longline survey CPUE 
# using all survey stations where yelloweye have been seen at least once: 
#ind<-read.csv(paste0("Data_processing/Data/IPHC.cpue.SEO_non0_",YEAR,".csv"))
#Or.. use more restrictive... this using stations where YE seen at least 40% of the time: 
#ind<-read.csv(paste0("Data_processing/Data/IPHC.cpue.SEO_min40percentYE_",YEAR,".csv"))
# with non0 stations and bootstrap index for comparison
ind<-read.csv(paste0(here::here(), "/Data_processing/Data/IPHC.cpue.SEO_non0_2022methods_plusHADJ_kghook_lon_2024.csv"))
ind <-ind %>% select(strata = SEdist, year = Year, cpue = WCPUE.bootmean, cv = WCPUE.cv)

# use Tweedie model output following CIE review recommendations
# also use all IPHC survey stations < 250 fathoms, following CIE review recommendations
#ind<-read.csv(paste0("./Data_processing/Data/IPHC.cpue.SEO_tweed_boot_cas_2024.csv"))
#ind<-read.csv(paste0("./Data_processing/Data/IPHC.cpue.SEO_tweed_boot_rke_102824.csv"))

#ind <-ind %>% filter(CPUE == "Tweedie index") %>% select(strata = SEdist, year = Year, cpue, cv)

#allmods<-read.csv("Data/SEO_YE_Biomass_all_models.csv")# %>% 
#sq<-read.csv("Data/SEO_YE_Biomass_subd_100722.csv")
#  filter(Model_name == 21.1)

str(bio)
str(ind)
str(sq)

??prepare_rema_input
?fit_rema
?tidy_rema
?plot_rema
?compare_rema_models

#-------------------------------------------------------------------------
#Model 1: 
Scrapnotesin1<-prepare_rema_input(model_name = 'REMA_SEO_YE_strataPE',
                        multi_survey = TRUE,  #0 or 1? 
                        biomass_dat = bio,
                        cpue_dat = ind,
                        extra_biomass_cv = list(assumption="extra_cv", #NULL,
                                                pointer_extra_biomass_cv = c(1,1,1,1)), #"extra_cv",
                        extra_cpue_cv = list(assumption="extra_cv",  #NULL
                                             pointer_extra_cpue_cv = c(1,1,1,1)), #"extra_cv",
                        PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1)),
                        q_options = list(pointer_q_cpue = c(1,2,3,4)))
#------------------------------------------------------------------------
# Model 22.1 single (1) process error, unique strata q values, 
# NO extra variance term for biomass estimates
in22_1<-prepare_rema_input(model_name = '22.1',
                           multi_survey = TRUE,  #0 or 1? 
                           biomass_dat = bio,
                           cpue_dat = ind,
                           extra_biomass_cv = NULL, #list(assumption="extra_cv", #NULL,
                                                   #pointer_extra_biomass_cv = c(1,1,1,1)), #"extra_cv",
                           extra_cpue_cv = NULL, #list(assumption="extra_cv",  #NULL
                           #pointer_extra_cpue_cv = c(1,1,1,1)), #"extra_cv",
                           PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1)),
                           q_options = list(pointer_q_cpue = c(1,2,3,4)))

names(in22_1)

m22_1 <- fit_rema(in22_1)
check_convergence(m22_1)
output22_1 <- tidy_rema(m22_1)  #output from model and return estimates and clean data frames..
names(output22_1)
print(output22_1$total_predicted_biomass,n = 30)

output22_1$biomass_by_cpue_strata
# can ignore because not relevant
output22_1$parameter_estimates
?tidy_rema

plots22_1 <- plot_rema(output22_1, biomass_ylab = 'ROV biomass', cpue_ylab = 'IPHC setline survey CPUE')
plots22_1$biomass_by_strata
plots22_1$cpue_by_strata
cowplot::plot_grid(plots22_1$biomass_by_strata + facet_wrap(~strata, ncol = 1),
                   plots22_1$cpue_by_strata + facet_wrap(~strata, ncol = 1),
                   nrow = 1)
plots22_1$total_predicted_cpue # note that total cpue is not available because nominal cpue is not summable
plots22_1$total_predicted_biomass

#------------------------------------------------------------------------
# Model 22.2 ----
# single (1) process error, unique strata q values, 
# one extra variance term for biomass estimates
in22_2<-prepare_rema_input(model_name = '22.2',
                           multi_survey = TRUE,  #0 or 1? 
                           biomass_dat = bio,
                           cpue_dat = ind,
                           extra_biomass_cv = list(assumption="extra_cv", #NULL,
                                                   pointer_extra_biomass_cv = c(1,1,1,1)), #"extra_cv",
                           extra_cpue_cv = NULL, #list(assumption="extra_cv",  #NULL
                                                #pointer_extra_cpue_cv = c(1,1,1,1)), #"extra_cv",
                           PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1)),
                           q_options = list(pointer_q_cpue = c(1,2,3,4)))

names(in22_2)

m22_2 <- fit_rema(in22_2)
check_convergence(m22_2)
output22_2 <- tidy_rema(m22_2)  #output from model and return estimates and clean data frames..
names(output22_2)
print(output22_2$total_predicted_biomass,n = 31)

output22_2$biomass_by_cpue_strata
# can ignore because not relevant
output22_2$parameter_estimates

plots22_2 <- plot_rema(output22_2, biomass_ylab = 'ROV biomass', cpue_ylab = 'IPHC setline survey CPUE')
plots22_2$biomass_by_strata
plots22_2$cpue_by_strata
cowplot::plot_grid(plots22_2$biomass_by_strata + facet_wrap(~strata, ncol = 1),
                   plots22_2$cpue_by_strata + facet_wrap(~strata, ncol = 1),
                   nrow = 1)
plots22_2$total_predicted_cpue # note that total cpue is not available because nominal cpue is not summable
plots22_2$total_predicted_biomass

get_osa_residuals(m22_2)

ggsave(file.path(paste0(here::here(), "/REMA/Figures/biomass_by_strata_22_2.png")), plot = plots22_2$biomass_by_strata, height = 0.6*7, width = 7, units = "in")

ggsave(file.path(paste0(here::here(), "/REMA/Figures/cpue_by_strata_22_2.png")), plot = plots22_2$cpue_by_strata, height = 0.6*7, width = 7, units = "in")

ggsave(file.path(paste0(here::here(), "/REMA/Figures/total_predicted_biomass_22_2.png")), plot = plots22_2$total_predicted_biomass, height = 0.6*7, width = 7, units = "in")


write.csv(output22_2$total_predicted_biomass, file = paste0(here::here(), "/REMA/Output/REMA_total_predicted_biomass_22_2_2022methods_plusHADJ_kghook.csv"))

# plots for SAFE ----

tidy_22_2 <- tidy_rema(m22_2)

p1 <- plot_rema(tidy_22_2)$biomass_by_strata + 
  ggtitle(label = "Model fits to the ADF&G survey data",
          subtitle = "ROV survey biomass strata") + 
  geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci),
              col = cbpalette[4], fill = cbpalette[4], alpha = 0.4) +
  geom_line() +
  geom_point(aes(x = year, y = obs)) +
  geom_errorbar(aes(x = year, ymin = obs_lci, ymax = obs_uci)) +
  ylab("SEO yelloweye biomass (t)") +
  theme_bw()

p1

ggsave(p1, file = paste0(here::here(), "/REMA/Figures/2024/adfg_survey_fits_22_2_nomCPUE.png"), width = 6, height = 4, units = "in", bg = "white")


p2 <- plot_rema(tidy_22_2)$cpue_by_strata + 
  ggtitle(label = "Model fits to the IPHC survey data",
          subtitle = "IPHC survey strata") + 
  geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci),
              col = cbpalette[4], fill = cbpalette[4], alpha = 0.4) +
  geom_line() +
  geom_point(aes(x = year, y = obs)) +
  geom_errorbar(aes(x = year, ymin = obs_lci, ymax = obs_uci)) +
  ylab("CPUE (kg per hook)") +
  theme_bw()

p2

ggsave(p2, file = paste0(here::here(), "/REMA/Figures/2024/iphc_survey_fits_22_2_nomCPUE.png"), width = 6, height = 4, units = "in", bg = "white")

p3 <- plot_rema(tidy_22_2)$total_predicted_biomass + 
  #ggtitle(label = "Total predicted biomass") + 
  geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci),
              col = cbpalette[4], fill = cbpalette[4], alpha = 0.4) +
  geom_line() +
  #geom_point(aes(x = year, y = obs)) +
  #geom_errorbar(aes(x = year, ymin = obs_lci, ymax = obs_uci))
  xlab("Year") +
  ylab("SEO yelloweye biomass (t)") +
  scale_x_continuous(breaks = seq(1994, 2024, 4)) +
  theme_bw()

p3

ggsave(p3, file = paste0(here::here(), "/REMA/Figures/2024/est_biomass_22_2_nomCPUE.png"), width = 6, height = 4, units = "in", bg = "white")

p4 <- cowplot::plot_grid(plot_rema(tidy_22_2)$biomass_by_strata + 
                     ggtitle(label = "Model fits to the ADF&G survey biomass") +
                     geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci),
                                 col = cbpalette[4], fill = cbpalette[4], alpha = 0.4) +
                     geom_line() +
                     geom_point(aes(x = year, y = obs)) +
                     geom_errorbar(aes(x = year, ymin = obs_lci, ymax = obs_uci)) +
                     xlab("Year") +
                     ylab("Biomass (t)") +
                     theme_bw() +
                     facet_wrap(~strata, ncol = 4),
                   plot_rema(tidy_22_2)$cpue_by_strata + 
                     ggtitle(label = "Model fits to the IPHC survey CPUE") +
                     geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci),
                                 col = cbpalette[4], fill = cbpalette[4], alpha = 0.4) +
                     geom_line() +
                     geom_point(aes(x = year, y = obs)) +
                     geom_errorbar(aes(x = year, ymin = obs_lci, ymax = obs_uci)) +
                     xlab("Year") +
                     ylab("IPHC setline survey CPUE") +
                     theme_bw() + facet_wrap(~strata, ncol = 4),
                   nrow = 2)

p4

ggsave(p4, file = paste0(here::here(), "/REMA/Figures/2024/fits_22_2_nomCPUE.png"), width = 10, height = 4, units = "in", bg = "white")


# standardized index plots ----

# first standardize each survey time series to its mean
sur.1 <- in22_2$biomass_dat %>%
  rename("obs_index" = "biomass") %>%
  mutate(index = "Biomass index") %>%
  rbind(in22_2$cpue_dat %>%
          rename("obs_index" = "cpue") %>%
          mutate(index = "CPUE index")) %>%
  group_by(index, strata) %>%
  mutate(mean_index = mean(obs_index, na.rm = T), sd_index = sd(obs_index, na.rm = T)) %>%
  mutate(std_index = (obs_index - mean_index) / sd_index) %>%
  mutate(std_l95 = std_index * exp(-1.96 * sqrt(log(1 + cv^2))),
         std_u95 = std_index * exp(1.96 * sqrt(log(1 + cv^2))))

# plot surveys on one plot ----

sur.1.plot <- ggplot(sur.1, aes(x = year, y = std_index, color = index)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(x = year, ymin = std_l95, ymax = std_u95, color = index), width = 0) +
  scale_color_manual(values = c(cbpalette[2], cbpalette[3])) +
  xlab("Year") +
  ylab("Standardized index") + 
  labs(color="Index") +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  scale_x_continuous(breaks = seq(min(sur.1$year), max(sur.1$year), by = 5)) + 
  theme_bw() +
  facet_wrap(~strata, ncol = 2)

ggsave(plot = sur.1.plot, 
       filename = paste0(here::here(), "/REMA/Figures/2024/stand_indices.png"), 
       width = 7, 
       height = 7 * (4/6), units = "in",
       bg = "white")

# version with larger text for presentation

sur.1.plot.pres <- ggplot(sur.1, aes(x = year, y = std_index, color = index)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(x = year, ymin = std_l95, ymax = std_u95, color = index), width = 0) +
  scale_color_manual(values = c(cbpalette[2], cbpalette[3])) +
  xlab("Year") +
  ylab("Standardized index") + 
  labs(color="Index") +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  scale_x_continuous(breaks = seq(min(sur.1$year), max(sur.1$year), by = 5)) + 
  theme_bw() +
  theme(legend.title=element_blank()) +
  theme(legend.text = element_text(size = 20)) +
  theme(strip.text = element_text(size = 20)) +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~strata, ncol = 2)

ggsave(plot = sur.1.plot.pres, 
       filename = paste0(here::here(), "/REMA/Figures/2024/stand_indices_pres.png"), 
       width = 9, 
       height = 7, units = "in",
       bg = "white")


# compare to previous assessment's biomass trajectory ----

biomass.prev <- read.csv(paste0(here::here(), "/REMA/Output/REMA_total_predicted_biomass_22_2_2022.csv"))

compare.prev <- output22_2$total_predicted_biomass %>%
  mutate(model = "22.2 - 2024") %>%
  select(model, year, pred, pred_lci, pred_uci) %>%
  rbind(biomass.prev %>%
          mutate(model = "22.2 - 2022"))

cp.plot <- ggplot(compare.prev) +
  geom_line(aes(x = year, y = pred, group = model, color = model)) +
  geom_ribbon(aes(x = year, ymin = pred_lci, ymax = pred_uci, group = model, fill = model), alpha = 0.25) +
  scale_color_manual(values = c(cbpalette[1], cbpalette[2])) +
  scale_fill_manual(values = c(cbpalette[1], cbpalette[2])) +
  xlab("Year") +
  ylab("Estimated biomass (t)") + 
  labs(color="Model", fill = "Model") +
  scale_x_continuous(breaks = seq(min(sur.1$year), max(sur.1$year), by = 2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_bw()

ggsave(cp.plot, file = paste0(here::here(), "/REMA/Figures/2024/est_biomass_compare.png"), width = 6, height = 4, units = "in", bg = "white")

ggsave(cp.plot, file = paste0(here::here(), "/REMA/Figures/2024/est_biomass_compare_pres.png"), width = 8, height = 4, units = "in", bg = "white")

# compare the IPHC CPUE index used in the previous assessment with the index used in this assessment ----

# CV = SD / mean
# SD = CV * mean
# SE = SD / sqrt(n)

cpue.prev <- read.csv(paste0(here::here(), "/Data_processing/Data/IPHC.cpue.SEO_non0_2022.csv")) %>%
  rename("CPUE" = "Mean") %>%
  mutate(model = "22.2 - 2022") %>%
  mutate(SE = (CV * CPUE) / sqrt(no.stations)) %>%
  mutate(lower = CPUE + SE * qnorm(0.05 / 2), upper = CPUE + SE * qnorm(1 - 0.05 / 2)) %>%
  select(c(Year, Stratum, CPUE, CV, upper, lower, model))
cpue.curr <- read.csv(paste0(here::here(), "/Data_processing/Data/IPHC.cpue.SEO_non0_2022methods_2024.csv")) %>%
  select(Year, Stratum = SEdist, CPUE = CPUE.bootmean, CV = CPUE.cv, upper, lower) %>%
  #filter(CPUE == "Tweedie index") %>%
  #select(c(Year, SEdist, cpue, cv, upper, lower)) %>%
  #rename("Stratum" = "SEdist", "CPUE" = "cpue", "CV" = "cv") %>%
  mutate(model = "22.2 - 2024")

cpue.compare <- rbind(cpue.prev, cpue.curr)

cpue.plot <- ggplot(cpue.compare) +
  geom_point(aes(x = Year, y = CPUE, group = model, color = model)) +
  geom_ribbon(aes(x = Year, ymin = lower, ymax = upper, group = model, fill = model), alpha = 0.25) +
  scale_color_manual(values = c(cbpalette[1], cbpalette[2])) +
  scale_fill_manual(values = c(cbpalette[1], cbpalette[2])) +
  xlab("Year") +
  ylab("CPUE") + 
  labs(color="Model", fill = "Model") +
  theme_bw() +
  facet_wrap(~Stratum, ncol = 2)


# Calculation of OFL and ABC ----

# for Tier 5, F_OFL is set equal to M, which in the past was 0.02 but for this year is 0.044
# F_ABCmax = 0.75*M = 0.015

F_OFL <- 0.044
F_ABCmax <- 0.75*F_OFL

OFL_ABC <- output22_2$total_predicted_biomass %>% 
  filter(year == 2024) %>%
  mutate(ABCmax = F_ABCmax * pred, OFL = F_OFL * pred, F.OFL = F_OFL, F.ABCmax = F_ABCmax)

write.csv(OFL_ABC, file = paste0(here::here(), "/REMA/Output/REMA_ref_pts_22_2.csv"))

nonYE.ofl <- 26
nonYE.abc <- 20
tot.DSR.ofl <- nonYE.ofl + OFL_ABC$OFL
tot.DSR.maxabc <- nonYE.abc + OFL_ABC$ABCmax

# parameter estimates

par.est <- tidy_22_2$parameter_estimates %>%
  mutate(across(where(is.numeric), ~ round(., 8)))

write.csv(par.est, file = paste0(here::here(), "/REMA/Output/REMA_param_est_22_2.csv"))


#------------------------------------------------------------------------
# Model 22.7 ----
# single (1) process error, unique strata q values, 
#  extra cv only on CPUE index, not on biomass estimates

in22_7<-prepare_rema_input(model_name = '22.7',
                           multi_survey = TRUE,  #0 or 1? 
                           biomass_dat = bio,
                           cpue_dat = ind,
                           extra_biomass_cv = #list(assumption="extra_cv", #NULL,
                             NULL,
                                                   #pointer_extra_biomass_cv = c(1,1,1,1)), #"extra_cv",
                           extra_cpue_cv = list(assumption="extra_cv",  #NULL
                                                pointer_extra_cpue_cv = c(1,1,1,1)),
                           PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1)),
                           q_options = list(pointer_q_cpue = c(1,2,3,4)))

names(in22_7)

m22_7 <- fit_rema(in22_7)
check_convergence(m22_7)
output22_7 <- tidy_rema(m22_7)  #output from model and return estimates and clean data frames..
names(output22_7)
print(output22_7$total_predicted_biomass,n = 30)

output22_7$biomass_by_cpue_strata
# can ignore because not relevant
output22_7$parameter_estimates

plots22_7 <- plot_rema(output22_7, biomass_ylab = 'ROV biomass', cpue_ylab = 'IPHC setline survey CPUE')
plots22_7$biomass_by_strata
plots22_7$cpue_by_strata
cowplot::plot_grid(plots22_7$biomass_by_strata + facet_wrap(~strata, ncol = 1),
                   plots22_7$cpue_by_strata + facet_wrap(~strata, ncol = 1),
                   nrow = 1)
plots22_7$total_predicted_cpue # note that total cpue is not available because nominal cpue is not summable
plots22_7$total_predicted_biomass

#-----------------------------------------------------------------------------
# As per Sept 2022 Plan Team instructions, rerunning models without IPHC CPUE data...
# 22.3 is the SPM so skipping that designation
# 22.4 will be 22.1 without IPHC CPUE
# 22.5 will be 22.2 without IPHC CPUE

#--------------------------------------------------------------------------------
# 22.4  single (1) process error, NO extra variance term for biomass estimates
in22_4<-prepare_rema_input(model_name = '22.4',
                           multi_survey = FALSE,  #0 or 1? 
                           biomass_dat = bio,
                           extra_biomass_cv = NULL, #list(assumption="extra_cv", #NULL,
                           PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1)))

names(in22_4)

m22_4 <- fit_rema(in22_4)
check_convergence(m22_4)
output22_4 <- tidy_rema(m22_4)  #output from model and return estimates and clean data frames..
names(output22_4)
print(output22_4$total_predicted_biomass,n = 30)

output22_4$biomass_by_cpue_strata
# can ignore because not relevant
output22_4$parameter_estimates

plots22_4 <- plot_rema(output22_4, biomass_ylab = 'ROV biomass', cpue_ylab = 'IPHC setline survey CPUE')
plots22_4$biomass_by_strata
plots22_4$cpue_by_strata
cowplot::plot_grid(plots22_4$biomass_by_strata + facet_wrap(~strata, ncol = 1),
                   #plots22_4$cpue_by_strata + facet_wrap(~strata, ncol = 1),
                   nrow = 1)
#plots22_4$total_predicted_cpue # note that total cpue is not available because nominal cpue is not summable
plots22_4$total_predicted_biomass

print(output22_4$total_predicted_biomass, n = 31)

#------------------------------------------------------------------------
# Model 22.2 single (1) process error, one extra variance term for biomass estimates
in22_5<-prepare_rema_input(model_name = '22.5',
                           multi_survey = FALSE,  #0 or 1? 
                           biomass_dat = bio,
                           extra_biomass_cv = list(assumption="extra_cv", #NULL,
                                                   pointer_extra_biomass_cv = c(1,1,1,1)), #"extra_cv",
                           PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1)))

names(in22_5)

m22_5 <- fit_rema(in22_5)
check_convergence(m22_5)
output22_5 <- tidy_rema(m22_5)  #output from model and return estimates and clean data frames..
names(output22_5)
print(output22_5$total_predicted_biomass,n = 30)

output22_5$parameter_estimates

plots22_5 <- plot_rema(output22_5, biomass_ylab = 'ROV biomass', cpue_ylab = 'IPHC setline survey CPUE')
plots22_5$biomass_by_strata
cowplot::plot_grid(plots22_5$biomass_by_strata + facet_wrap(~strata, ncol = 1),
                   nrow = 1)
plots22_5$total_predicted_biomass

#------------------------------------------------------------------------
# Model 22.6 single (1) process error, unique strata q values, 
#  extra variance term for biomass estimates & IPHC CPUE
in22_6<-prepare_rema_input(model_name = '22.6',
                           multi_survey = TRUE,  #0 or 1? 
                           biomass_dat = bio,
                           cpue_dat = ind,
                           extra_biomass_cv = list(assumption="extra_cv", #NULL,
                                                   pointer_extra_biomass_cv = c(1,1,1,1)), #"extra_cv",
                           extra_cpue_cv = list(assumption="extra_cv",  #NULL
                                                  pointer_extra_cpue_cv = c(1,1,1,1)), #"extra_cv",
                           PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1)),
                           q_options = list(pointer_q_cpue = c(1,2,3,4)))

names(in22_6)

m22_6 <- fit_rema(in22_6)
check_convergence(m22_6)
output22_6 <- tidy_rema(m22_6)  #output from model and return estimates and clean data frames..
names(output22_6)
print(output22_6$total_predicted_biomass,n = 30)

output22_6$biomass_by_cpue_strata
# can ignore because not relevant
output22_6$parameter_estimates

plots22_6 <- plot_rema(output22_6, biomass_ylab = 'ROV biomass', cpue_ylab = 'IPHC setline survey CPUE')
plots22_6$biomass_by_strata
plots22_6$cpue_by_strata
cowplot::plot_grid(plots22_6$biomass_by_strata + facet_wrap(~strata, ncol = 1),
                   plots22_6$cpue_by_strata + facet_wrap(~strata, ncol = 1),
                   nrow = 1)
plots22_6$total_predicted_cpue # note that total cpue is not available because nominal cpue is not summable
plots22_6$total_predicted_biomass

#------------------------------------------------------------------------
# Model 22.6 single (1) process error, unique strata q values, 
#  extra variance term for biomass estimates & IPHC CPUE
in22_7<-prepare_rema_input(model_name = '22.7',
                           multi_survey = TRUE,  #0 or 1? 
                           biomass_dat = bio,
                           cpue_dat = ind,
                           extra_biomass_cv = list(assumption="extra_cv", #NULL,
                                                   pointer_extra_biomass_cv = c(1,1,1,1)), #"extra_cv",
                           extra_cpue_cv = list(assumption="extra_cv",  #NULL
                                                pointer_extra_cpue_cv = c(1,2,3,4)), #"extra_cv",
                           PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1)),
                           q_options = list(pointer_q_cpue = c(1,2,3,4)))

names(in22_7)

m22_7 <- fit_rema(in22_7)
check_convergence(m22_7)
output22_7 <- tidy_rema(m22_7)  #output from model and return estimates and clean data frames..
names(output22_7)
print(output22_7$total_predicted_biomass,n = 30)

output22_7$biomass_by_cpue_strata
# can ignore because not relevant
output22_7$parameter_estimates

plots22_7 <- plot_rema(output22_7, biomass_ylab = 'ROV biomass', cpue_ylab = 'IPHC setline survey CPUE')
plots22_7$biomass_by_strata
plots22_7$cpue_by_strata
cowplot::plot_grid(plots22_7$biomass_by_strata + facet_wrap(~strata, ncol = 1),
                   plots22_7$cpue_by_strata + facet_wrap(~strata, ncol = 1),
                   nrow = 1)
plots22_7$total_predicted_cpue # note that total cpue is not available because nominal cpue is not summable
plots22_7$total_predicted_biomass

#-------------------------------------------------------------------------------
# COMPARE MODELS...
?compare_rema_models
compare_bio_IPHC <- compare_rema_models(rema_models = list(m22_1, m22_2, m22_6,m22_7),
                               biomass_ylab = 'ROV biomass',
                               cpue_ylab = 'IPHC setline survey CPUE')

compare_bio_IPHC$aic # The single process error model has the lowest AIC and is the simplest model
names(compare_bio_IPHC$plots)

compare_bio_IPHC$plots$total_predicted_biomass

cowplot::plot_grid(compare_bio_IPHC$plots$biomass_by_strata +
                     facet_wrap(~strata, nrow = 1) +
                     theme(legend.position = 'top'),
                   compare_bio_IPHC$plots$cpue_by_strata +
                     facet_wrap(~strata, nrow = 1) +
                     theme(legend.position = 'none'),
                   nrow = 2) # use , rel_heights = c(0.52, 0.48)) if you want to get anal-retentive about it...

#biomass only models
compare_bio <- compare_rema_models(rema_models = list(m22_4, m22_5),
                                        biomass_ylab = 'ROV biomass')

compare_bio$aic # The single process error model has the lowest AIC and is the simplest model
names(compare_bio$plots)

compare_bio$plots$total_predicted_biomass

cowplot::plot_grid(compare_bio$plots$biomass_by_strata +
                     facet_wrap(~strata, nrow = 1) +
                     theme(legend.position = 'top'),
                   nrow = 1)
ggsave(paste0("REMA/Figures/REMA_22.5_22.5_MA_biomass_comp_", YEAR, ".png"), dpi=300,  height=4, width=7, units="in")
#-------------------------------------------------------------------------------
# Make plots comparing status quo and REMA models for SAFE report
#1) make data frame to plot different models...
pal<-c("mediumorchid4","dodgerblue4","chocolate1","cyan4")

bio22_1<-as.data.frame(output22_1$total_predicted_biomass)
bio22_2<-as.data.frame(output22_2$total_predicted_biomass)
bio22_4<-as.data.frame(output22_4$total_predicted_biomass)
bio22_5<-as.data.frame(output22_5$total_predicted_biomass)
bio22_6<-as.data.frame(output22_6$total_predicted_biomass)
bio22_7<-as.data.frame(output22_7$total_predicted_biomass)

# this code compared the biomass estimation methods pre-2022.  REMA estimates have
# been accepted as preferred methods and this is defunct. 
#sq.comp<-sq %>% select(year=Year, pred = Biomass_mt, cv=Biomass.cv) %>%
#  mutate(model_name = "status_quo (21.1)", 
#         variable = "tot_biomass_pred",
#         pred_lci = pred-1.96*cv*pred,
#         pred_uci = pred+1.96*cv*pred)
#str(sq.comp)

mods<-rbind(bio22_1,bio22_2,bio22_4,bio22_5,bio22_6,bio22_7)
str(mods)

ggplot(mods, aes(x=year, col=model_name, fill=model_name)) +
  geom_ribbon(aes(ymin = pred_lci, ymax= pred_uci,col=NULL, fill=model_name), alpha=0.3) +
  geom_line(aes(y=pred, col=model_name, linetype = model_name),size=1) +
  xlab("\nYear") +
  ylab("SEO yelloweye biomass (t)\n") +
  scale_y_continuous(label=comma, breaks = c(10000,15000,20000,25000,30000,35000,40000)) +
  scale_x_continuous(breaks=seq(1995,2025,5)) + 
  theme (panel.grid.minor = element_blank()) + 
  #scale_color_brewer (palette = "Set1") +
  scale_color_viridis(discrete = TRUE, option = "H", end=0.75)+
  scale_fill_viridis(discrete = TRUE, option = "H", end=0.75)

# --------------------------------------------------------------------------------
# Get table of biomass estimates
str(mods)

write.csv(mods,"REMA/Output/REAM_Biomass_ests_for_SAFE_table.csv")

mods_tab<-mods
mods_tab[,c(4:6)]<-formatC(as.numeric(unlist(mods_tab[,c(4:6)])),format="d",digits=0,big.mark=",")
  
opts <- options(knitr.kable.NA = "-")            
table<-mods_tab %>% #group_by(Year) %>% 
  select(-variable) %>%
  rename(Model = model_name, Year = year, Biomass = pred, Lower_CI = pred_lci, Upper_CI = pred_uci) %>% 
  kbl(digits=0) %>% 
  column_spec(3) %>% #,format.args = list(big.mark = ",",scientific = FALSE)) %>%
      #format.args = list(big.mark = ",",
      #                            scientific = FALSE)) %>% #, col.names = gsub("[.]", " ", names(head(mods)))) %>% #kable_styling
  #kable_paper("hover",full_width=F)
  kable_classic(full_width=F, html_font = "Times", position="center")
#kable_classic_2(full_width=F, position = "left")
save_kable(table,"REMA/Output/REMA_Biomass_Table.pdf")

#--------------------------------------------------------------------------------
# More comparson plots

NV<-ggplot(mods %>% filter(model_name == "22.1" | model_name == "22.4" ),
       aes(x=year)) +
  geom_line(aes(y=pred, col=model_name, linetype = model_name),size=0.5) +
  geom_ribbon(aes(ymin = pred_lci, ymax= pred_uci,col=NULL, fill=model_name), alpha=0.3) +
  xlab("") +
  ylab("biomass (t)") +
  labs(title = "no extra variance on biomass estimates") +
  scale_y_continuous(label=comma, breaks = c(10000,15000,20000,25000,30000,35000,40000)) +
  scale_x_continuous(breaks=seq(1995,2025,5)) + 
  theme (panel.grid.minor = element_blank()) +
  #scale_color_viridis(discrete = TRUE, option = "A", end=0.75)+
  #scale_fill_viridis(discrete = TRUE, option = "A", end=0.75)
  scale_color_manual(values = c(pal[1],pal[3],"black")) +
  scale_fill_manual(values = c(pal[1],pal[3],"black"))
#ggsave(paste0("Figures/REMA_IPHC_comp_noXtraPE_", YEAR, ".png"), dpi=300,  height=4, width=6, units="in")

XV<-ggplot(mods %>% filter(model_name == "22.2" | model_name == "22.5"),
       aes(x=year)) +
  geom_line(aes(y=pred, col=model_name, linetype = model_name),size=0.5) +
  geom_ribbon(aes(ymin = pred_lci, ymax= pred_uci,col=NULL, fill=model_name), alpha=0.25) +
  xlab("Year") +
  ylab("biomass (t)") +
  labs(title = "extra variance on biomass estimates") +
  scale_y_continuous(label=comma, breaks = c(10000,15000,20000,25000,30000,35000,40000)) +
  scale_x_continuous(breaks=seq(1995,2025,5)) + 
  theme (panel.grid.minor = element_blank()) +
  scale_color_manual(values = c(pal[2],pal[4],"black")) +
  scale_fill_manual(values = c(pal[2],pal[4],"black"))

#(paste0("Figures/REMA_IPHC_comp_XtraPE_", YEAR, ".png"), dpi=300,  height=4, width=6, units="in")
figure<-ggarrange(NV, XV, labels = c("A","B"), ncol=1, nrow=2)
figure
ggsave(paste0("REMA/Figures/REMA_totbio_modcomp_", YEAR, ".png"), dpi=300,  height=4, width=6, units="in")

#---------------------------------------------------------------------------------
bio22_1s<-as.data.frame(output22_1$biomass_by_strata)
bio22_2s<-as.data.frame(output22_2$biomass_by_strata)
bio22_4s<-as.data.frame(output22_4$biomass_by_strata)
bio22_5s<-as.data.frame(output22_5$biomass_by_strata)
bio22_6s<-as.data.frame(output22_6$biomass_by_strata)
bio22_7s<-as.data.frame(output22_7$biomass_by_strata)

mods_strata<-rbind(bio22_1s,bio22_2s,bio22_4s,bio22_5s, bio22_6s, bio22_7s)

Norm<-ggplot(mods_strata %>% filter(model_name == "22.1" | model_name == "22.4" ),
       aes(x=year)) +
  geom_line(aes(y=pred, col=model_name, linetype = model_name),size=0.5) +
  geom_ribbon(aes(ymin = pred_lci, ymax= pred_uci,col=NULL, fill=model_name), alpha=0.25) +
  geom_point(aes(year, obs),size=1) +
  geom_errorbar(aes(ymin=obs_lci, ymax=obs_uci),size=0.25) +
  facet_wrap(~strata, scales="free_y", nrow=1) +
  xlab("") +
  ylab("biomass (t)") +
  labs(title = "no extra variance on biomass estimates") +
  scale_y_continuous(label=comma) + #, breaks = c(10000,15000,20000,25000,30000,35000,40000)) +
  scale_x_continuous(breaks=seq(1995,2025,5)) + 
  theme (panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c(pal[1],pal[3])) +
  scale_fill_manual(values = c(pal[1],pal[3]))

xtraV<-ggplot(mods_strata %>% filter(model_name == "22.2" | model_name == "22.5" ),
       aes(x=year)) +
  geom_line(aes(y=pred, col=model_name, linetype = model_name),size=0.5) +
  geom_ribbon(aes(ymin = pred_lci, ymax= pred_uci,col=NULL, fill=model_name), alpha=0.25) +
  geom_point(aes(year, obs),size=1) +
  geom_errorbar(aes(ymin=obs_lci, ymax=obs_uci),size=0.25) +
  facet_wrap(~strata, scales="free_y", nrow=1) +
  xlab("Year") +
  ylab("biomass (t)") +
  labs(title = "extra variance on biomass estimates") +
  scale_y_continuous(label=comma) + #, breaks = c(10000,15000,20000,25000,30000,35000,40000)) +
  scale_x_continuous(breaks=seq(1995,2025,5)) + 
  theme (panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c(pal[2],pal[4])) +
  scale_fill_manual(values = c(pal[2],pal[4]))

figure<-ggarrange(Norm, xtraV, labels = c("A","B"), ncol=1, nrow=2)
figure

ggsave(paste0("REMA/Figures/REMA_strata_modcomp_", YEAR, ".png"), dpi=300,  height=5, width=9, units="in")

#--------------------------------------------------------------------------------
cowplot::plot_grid(compare_bio_IPHC$plots$biomass_by_strata +
                     facet_wrap(~strata, nrow = 1) +
                     theme(legend.position = 'top'),
                   compare_bio_IPHC$plots$cpue_by_strata +
                     facet_wrap(~strata, nrow = 1) +
                     theme(legend.position = 'none'),
                   nrow = 2) 

bio22_1c<-as.data.frame(output22_1$cpue_by_strata)
bio22_2c<-as.data.frame(output22_2$cpue_by_strata)
bio22_6c<-as.data.frame(output22_6$cpue_by_strata)
bio22_7c<-as.data.frame(output22_7$cpue_by_strata)
mods_cpue<-rbind(bio22_1c,bio22_2c,bio22_6c,bio22_7c)

str(mods_cpue)
str(mods_strata)

Bio_IPHC<-ggplot(mods_strata %>% filter(model_name == "22.1" | model_name == "22.2" ),
              aes(x=year,color=model_name)) +
  geom_ribbon(aes(ymin = pred_lci, ymax= pred_uci,col=NULL, fill=model_name), alpha=0.25) +
  geom_line(aes(y=pred, col=model_name, linetype = model_name),size=0.5) +
  geom_point(aes(year, obs),size=1, col="black", legend=FALSE) +
  geom_errorbar(aes(ymin=obs_lci, ymax=obs_uci),size=0.25, col="black", legend=FALSE) +
  facet_wrap(~strata, scales="free_y", nrow=1) +
  xlab("Year") +
  ylab("ROV biomass (t)") +
  labs(title = "Biomass") +
  scale_y_continuous(label=comma) + #, breaks = c(10000,15000,20000,25000,30000,35000,40000)) +
  scale_x_continuous(breaks=seq(1995,2025,5)) + 
  theme (panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c(pal[1],pal[2])) +
  scale_fill_manual(values = c(pal[1],pal[2]))

  
CPUE_IPHC<-ggplot(mods_cpue %>% filter(model_name == "22.1" | model_name == "22.2" ),
              aes(x=year, 
                  color=model_name)) +
  geom_ribbon(aes(ymin = pred_lci, ymax= pred_uci,col=NULL, fill=model_name), alpha=0.25) +
  geom_line(aes(y=pred, col=model_name, linetype = model_name),size=0.5) +
  geom_point(aes(year, obs, col="black"),size=1, col="black", legend=FALSE) +
  geom_errorbar(aes(ymin=obs_lci, ymax=obs_uci),size=0.25, col="black", legend=FALSE) +
  facet_wrap(~strata, scales="free_y", nrow=1) +
  xlab("Year") +
  ylab("IPHC setline survey CPUE") +
  labs(title = "IPHC CPUE") +
  scale_y_continuous(label=comma) + #, breaks = c(10000,15000,20000,25000,30000,35000,40000)) +
  scale_x_continuous(breaks=seq(1995,2025,5)) + 
  theme (panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c(pal[1],pal[2])) +
  scale_color_manual(values = c(pal[1],pal[2])) 


figure<-ggarrange(Bio_IPHC, CPUE_IPHC, ncol=1, nrow=2, common.legend = TRUE, legend="top")
figure
ggsave(paste0("REMA/Figures/REMA_cpue_modcomp_", YEAR, ".png"), dpi=300,  height=5, width=9, units="in")

scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73",
                                                              "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))

#----------------------------------------------------------------
Bio_IPHC<-ggplot(mods_strata %>% filter(model_name == "22.2" | model_name == "22.6" | model_name == "22.7" ),
                 aes(x=year,color=model_name)) +
  geom_ribbon(aes(ymin = pred_lci, ymax= pred_uci,col=NULL, fill=model_name), alpha=0.25) +
  geom_line(aes(y=pred, col=model_name, linetype = model_name),size=0.5) +
  geom_point(aes(year, obs),size=1, col="black", legend=FALSE) +
  geom_errorbar(aes(ymin=obs_lci, ymax=obs_uci),size=0.25, col="black", legend=FALSE) +
  facet_wrap(~strata, scales="free_y", nrow=1) +
  xlab("Year") +
  ylab("ROV biomass (t)") +
  labs(title = "Biomass") +
  scale_y_continuous(label=comma) + #, breaks = c(10000,15000,20000,25000,30000,35000,40000)) +
  scale_x_continuous(breaks=seq(1995,2025,5)) + 
  theme (panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c(pal[1],pal[2],pal[3])) +
  scale_fill_manual(values = c(pal[1],pal[2],pal[3]))


CPUE_IPHC<-ggplot(mods_cpue %>% filter(model_name == "22.2" | model_name == "22.6" | model_name == "22.7" ),
                  aes(x=year, 
                      color=model_name)) +
  geom_ribbon(aes(ymin = pred_lci, ymax= pred_uci,col=NULL, fill=model_name), alpha=0.25) +
  geom_line(aes(y=pred, col=model_name, linetype = model_name),size=0.5) +
  geom_point(aes(year, obs, col="black"),size=1, col="black", legend=FALSE) +
  geom_errorbar(aes(ymin=obs_lci, ymax=obs_uci),size=0.25, col="black", legend=FALSE) +
  facet_wrap(~strata, scales="free_y", nrow=1) +
  xlab("Year") +
  ylab("IPHC setline survey CPUE") +
  labs(title = "IPHC CPUE") +
  scale_y_continuous(label=comma) + #, breaks = c(10000,15000,20000,25000,30000,35000,40000)) +
  scale_x_continuous(breaks=seq(1995,2025,5)) + 
  theme (panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c(pal[1],pal[2],pal[3])) +
  scale_color_manual(values = c(pal[1],pal[2],pal[3])) 


figure<-ggarrange(Bio_IPHC, CPUE_IPHC, ncol=1, nrow=2, common.legend = TRUE, legend="top")
figure
ggsave(paste0("REMA/Figures/REMA_cpue_modcomp2_", YEAR, ".png"), dpi=300,  height=5, width=9, units="in")





