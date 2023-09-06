################################################################################
## REMA for SEO Yelloweye rockfish
## from Jane Sullivan's package "rema" with guidance from Jane :)
## Rerunnig for Nov. 2022 plan team meeting
################################################################################
#If REMA package need to be loaded, here it is: 
# install.packages("devtools")
devtools::install_github("afsc-assessments/rema", dependencies = TRUE, build_vignettes = TRUE)

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

{library("rema")
library("scales")
library("viridis")
library("dplyr")
library("ggplot2")
library("ggpubr")
library("kableExtra")}

YEAR<-2023

#ROV biomass .. check date tag to get most recent
bio<-read.csv("Data_processing/Data/SEO_YE_Biomass_subdistrict_2023-06-21.csv")
str(bio_new)
str(bio)

bio<-bio %>% select(strata = Subdistrict, year = Year, biomass = Biomass.mt,
                    cv = Biomass.cv)

#IPHC longline survey CPUE 
# using all survey stations where yelloweye have been seen at least once: 
ind<-read.csv(paste0("Data_processing/Data/IPHC.cpue.SEO_non0_",YEAR,".csv"))
#Or.. use more restrictive... this using stations where YE seen at least 40% of the time: 
#ind<-read.csv(paste0("Data_processing/Data/IPHC.cpue.SEO_min40percentYE_",YEAR,".csv"))

ind<-ind %>% select(strata = mngmt.area, year = Year, cpue = CPUE.mean, cv = CPUE.cv)

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
# Model 22.2 single (1) process error, unique strata q values, 
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
print(output22_2$total_predicted_biomass,n = 30)

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





