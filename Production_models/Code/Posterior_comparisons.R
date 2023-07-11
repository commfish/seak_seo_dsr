################################################################################
## SS-SPM Model Comparison
## Code for comparing posteriors of key parameters and 
## Deviance and DIC values
## June 2022
## Phil Joy
###############################################################################

{library(tidyverse)
  library(ggridges)
  library(ggplot2)
  library(MCMCvis)
  library(runjags)
  library(rjags)
  library(dplyr)
  library(boot)
  library(ggplot2)
  library(coda)
  library(MASS)  
  library(R2OpenBUGS)
  library(jagsUI)
  library(ggmcmc)
  library(IDPmisc)
  library(scales)}

source("Code/Posterior_Plotting/YE_SPM_posterior_exams_Func.R")

#function for ranking models by DIC and deviance
ranker<-function(ModList){
  Model.Rank<-data.frame(); j<-1
  for (i in 1:length(ModList)) { # i<-1
    Mod<-ModList[[i]]
    Modname<-ModList[i]
    tbl<-jags.View(Mod, title="", digits=3)
    tbl<-as.data.frame(tbl)
    
    Model.Rank[i,"Model"]<-names(Modname)
    Model.Rank[i,"full_name"]<-Mod.fullname[i]
    Model.Rank[i,"deviance"]<-tbl$mean[tbl$parameter == "deviance"]
    Model.Rank[i,"pD"]<-Mod$pD
    Model.Rank[i,"DIC"]<-Mod$DIC
  }
  Model.Rank<-Model.Rank[order(Model.Rank$DIC),]
  return(Model.Rank)
}

#Load models from appropriate folder and label them.  All will be called "post"
#when leaded
# lists from MULTIMODEL_RUN...: 
{
Mod0e<-"PT_2i_Ope_TEST1_0.3m"
Mod0a<-"PT_3i_pe_TEST1_0.5m"
Mod0<-"PT_2i_pe02_T1e_r0405_d1_nofishing_0.3m"  
Mod1<-"PT_2i_pe015_T1e_r0405_d1_nofishing_0.35m"
Mod2<-"PT_2i_pe015_T1n_r0405_d1_nofishing_0.3m"
Mod3<-"PT_2i_pe015_T12e_r0405_d1_nofishing_0.3m"
Mod4<-"PT_2i_pe01_T1e_r0405_d1_nofishing_0.3m"
#drop '94 from best model
#drop '94 from best model
Mod5<-"PT_2i_pe015_T1n_r0405_no94_d1_nofishing_0.25m"
Mod6<-"PT_2i_pe01_T1n_r0405_no94_d1_nofishing_0.25m"
#r sensitivity
Mod7<-"PT_2i_pe015_T1e_r0410_d1_nofishing_0.3m"
Mod8<-"PT_2i_pe015_T1e_r0510_d1_nofishing_0.3m"
Mod9<-"PT_2i_pe015_T1e_r0610_d1_nofishing_0.3m"
Mod10<-"PT_2i_pe015_T1e_runi_d1_nofishing_0.3m"
Mod11<-"PT_2i_pe015_T1e_rg2r50_d1_nofishing_0.35m"
Mod12<-"PT_2i_pe015_T1e_rg5r100_d1_nofishing_0.35m"
Mod13<-"PT_2i_pe015_T1e_rg9r200_d1_nofishing_0.35m"

Mod14<-"PT_2i_pe015_T1e_rgamv035_d1_nofishing_0.3m"
#Mod15<-"PT_2i_pe015_T1e_rgamv009_d1_nofishing_0.35m"
#Mod16<-"PT_2i_pe015_T1e_rgamv0034_d1_nofishing_0.35m"
Mod17<-"PT_2i_pe015_T1e_rgamv0012_d1_nofishing_0.3m"
Mod18<-"PT_2i_pe015_T1e_runif_d1_nofishing_0.3m"

Mod19<-"PT_2i_pe015_T1e_rln0405_d1_nofishing_0.3m"
Mod20<-"PT_2i_pe015_T1e_rln0410_d1_nofishing_0.3m"
#Mod21<-"PT_2i_pe015_T1e_rln0505_d1_nofishing"
Mod22<-"PT_2i_pe015_T1e_rln1010_d1_nofishing_0.3m"
}
#r comparisons Beta models starting on 7-15-22; r comp from model development
{Mod23<-"PT_2i_pe01_T1e_rbeta_B1-1_B2-1__derb_0_1000K"
#Mod23<-  "PT_2i_pe01_T1e_rbeta_B1-1_B2-1_upv-3_derb_0_1400k"
Mod24<-"PT_2i_pe01_T1e_rbeta_B1-1.43_B2-19__derb_0_1000K"
Mod25<-"PT_2i_pe01_T1e_rbeta_B1-1.52_B2-35__derb_0_1000K"
Mod26<-"PT_2i_pe01_T1e_rbeta_B1-1.61_B2-67__derb_0_1000K"

#Mod27<-"PT2i_base_PHASE3_norm-phi_B1-1_B2-1_phmu-0.7_phsig-1.2_Kmu-10.6_Ksig-0.3_derb_0_700k"
#Mod28<-"PT2i_fullcatch_B1-1.43_B2-19_derb_0_1600k"
}

#PHASE 1 MODELS
{phase1.3d0<-"PT_2i_pe01_T1e_rbeta_B1-1_B2-1_upv-3_derb_0_1400k"
phase1.3dundest<-"PT_2i_pe01_T1e_rbeta_B1-1_B2-1_upv-3_derb_0.3_1400k"
phase1.3dovest<-"PT_2i_pe01_T1e_rbeta_B1-1_B2-1_upv-3_derb_-0.3_1400k"

phase1.5d0<-"PT_2i_pe01_T1e_rbeta_B1-1_B2-1_upv-5_derb_0_1400k"
phase1.5dundest<-"PT_2i_pe01_T1e_rbeta_B1-1_B2-1_upv-5_derb_0.3_1400k"
phase1.5dovest<-"PT_2i_pe01_T1e_rbeta_B1-1_B2-1_upv-5_derb_-0.3_1400k"}

#Phase2 models
{phase2.3d0<-"PT2i_fullcatch_B1-1_B2-1_upv-3_derb_0_1600k"
  phase2.3dundest<-"PT2i_fullcatch_B1-1_B2-1_upv-3_derb_0.3_1600k"
  phase2.3dovest<-"PT2i_fullcatch_B1-1_B2-1_upv-3_derb_-0.3_1600k"
  
  phase2.5d0<-"PT2i_fullcatch_B1-1_B2-1_upv-5_derb_0_1600k"
  phase2.5dundest<-"PT2i_fullcatch_B1-1_B2-1_upv-5_derb_0.3_1600k"
  phase2.5dovest<-"PT2i_fullcatch_B1-1_B2-1_upv-5_derb_-0.3_1600k"}

#Phase3 models
{
  phase3.mod_PE.r_unif<-"PHASE3_B1-1_B2-1_upv-3_phmu-0.7_phsig-1.2_Kmu-10.7_Ksig-9.7_derb_0_1500k"
  phase3.mod_PE.r_broad<-"PHASE3_B1-1.48_B2-23_upv-3_phmu-0.7_phsig-1.2_Kmu-10.7_Ksig-9.7_derb_0_1500k"
  phase3.mod_PE.r_mod<-"PHASE3_B1-1.25_B2-31_upv-3_phmu-0.7_phsig-1.2_Kmu-10.7_Ksig-9.7_derb_0_1500k"
  phase3.mod_PE.r_narrow<-"PHASE3_B1-1.24_B2-53_upv-3_phmu-0.7_phsig-1.2_Kmu-10.7_Ksig-9.7_derb_0_1500k"
  
  phase3.min_PE.r_unif<-"PHASE3_B1-1_B2-1_upv-5_phmu-0.7_phsig-1.2_Kmu-10.6_Ksig-9.2_derb_0_1500k"
  phase3.min_PE.r_broad<-"PHASE3_B1-1.48_B2-23_upv-5_phmu-0.7_phsig-1.2_Kmu-10.6_Ksig-9.2_derb_0_1500k"
  phase3.min_PE.r_mod<-"PHASE3_B1-1.25_B2-31_upv-5_phmu-0.7_phsig-1.2_Kmu-10.6_Ksig-9.2_derb_0_1500k"
  phase3.min_PE.r_narrow<-"PHASE3_B1-1.24_B2-53_upv-5_phmu-0.7_phsig-1.2_Kmu-10.6_Ksig-9.2_derb_0_1500k"
}

#phase3 simulations
{
  phase3_real<-"PHASE3_B1-1_B2-1_upv-5_phmu-0.7_phsig-1.2_Kmu-10.6_Ksig-9.2_derb_0_1500k"
  sim1<-"postsim1"
  sim2<-"postsim2"
  sim3<-"postsim3"
  sim4<-"postsim4"
  sim5<-"postsim5"
  sim6<-"postsim6"
  sim7<-"postsim7"
  sim8<-"postsim8"
  sim9<-"postsim9"
  sim10<-"postsim10"
  sim11<-"postsim11"
  sim12<-"postsim12"
  sim13<-"postsim13"
  sim14<-"postsim14"
  sim15<-"postsim15"
}

#Risk Analysis
{
  Risk_nobias<-"RISK3_B1-1_B2-1_upv-5_derb_0_recABC_-0perc_1500k"
  Risk_biaslo<-"RISK3_B1-1_B2-1_upv-5_derb_0.3_recABC_-0perc_1500k"
  Risk_biashi<-"RISK3_B1-1_B2-1_upv-5_derb_-0.3_recABC_-0perc_1500k"
}
#=================================================================
# load models and name them

#------------------------------------------------------------
# R prior development models
Mod.fullname<-c(Mod23,Mod24,Mod25,Mod26)

load(file=paste("Model Output/",Mod23,"/post.Rdata", sep=""))
beta_uni<-post

load(file=paste("Model Output/",Mod24,"/post.Rdata", sep=""))
beta_broad<-post

load(file=paste("Model Output/",Mod25,"/post.Rdata", sep=""))
beta_mod<-post

load(file=paste("Model Output/",Mod26,"/post.Rdata", sep=""))
beta_narrow<-post

ModList.rdev<-list(
  beta_uni=beta_uni,
  beta_broad=beta_broad,
  beta_mod=beta_mod,
  beta_narrow=beta_narrow)


Modnames.rdev<-c("uniform beta(1,1)",
                 "broad beta(1.4,19)",
                 "moderate beta(1.5,35)",
                 "narrow beta(1.6,67)") 

MR.rdev<-ranker(ModList<-ModList.rdev); MR.rdev$test<-"r.dev"

Model.Rank[order(Model.Rank$deviance),]
Model.Rank[order(Model.Rank$DIC),]

write.csv(MR.rdev,
          file=paste("Model Output/Model_ranking_r.dev", Sys.Date(),".csv", sep=""))
#------------------------------------------------------------
# PHASE 1 comparisons 
Mod.fullname<-c(phase1.3d0,phase1.3dundest,phase1.3dovest,
                phase1.5d0,phase1.5dundest,phase1.5dovest)

load(file=paste("Model Output/",phase1.3d0,"/post.Rdata", sep=""))
lv3_derb0<-post

load(file=paste("Model Output/",phase1.3dundest,"/post.Rdata", sep=""))
lv3_derb3p<-post

load(file=paste("Model Output/",phase1.3dovest,"/post.Rdata", sep=""))
lv3_derb3o<-post

load(file=paste("Model Output/",phase1.5d0,"/post.Rdata", sep=""))
lv5_derb0<-post

load(file=paste("Model Output/",phase1.5dundest,"/post.Rdata", sep=""))
lv5_derb3p<-post

load(file=paste("Model Output/",phase1.5dovest,"/post.Rdata", sep=""))
lv5_derb3o<-post

ModList.st1<-list(
  lv3_derb0=lv3_derb0,
  lv3_derb3p=lv3_derb3p,
  lv3_derb3o=lv3_derb3o,
  lv5_derb0=lv5_derb0,
  lv5_derb3p=lv5_derb3p,
  lv5_derb3o=lv5_derb3o)

ModList.st1.sig<-list(
  lv3_derb0=lv3_derb0,
  lv5_derb0=lv5_derb0)

ModList.st1.db<-list(
  lv3_derb0=lv3_derb0,
  lv3_derb3p=lv3_derb3p,
  lv3_derb3o=lv3_derb3o)

sqrt(exp(-3))
sqrt(exp(-5))
Modnames.st1<-c("mod sig, derby by unbiased",
                 "mod sig, wcpue under derby by",
                "mod sig, wcpue over derby by",
                "min sig, derby by unbiased",
                "min sig, wcpue under derby by",
                "min sig, wcpue over derby by") 

Modnames.st1.sig<-c("log var ~ U(-10,-3)",
                "log var ~ U(-10,-5)")

Modnames.st1.db<-c("wcpue unbiased",
                "wcpue biased low",
                "wcpue biased high") 

MR.stage1<-ranker(ModList<-ModList.st1); MR.rdev$test<-"stage1"

write.csv(MR.stage1,
          file=paste("Model Output/Model_ranking_stage1.dev", Sys.Date(),".csv", sep=""))
#------------------------------------------------------------
# PHASE 2 COMPARISONS
Mod.fullname<-c(phase2.3d0,phase2.3dundest,phase2.3dovest,
                phase2.5d0,phase2.5dundest,phase2.5dovest)

load(file=paste("Model Output/",phase2.3d0,"/post.Rdata", sep=""))
lv3_derb0<-post

load(file=paste("Model Output/",phase2.3dundest,"/post.Rdata", sep=""))
lv3_derb3p<-post

load(file=paste("Model Output/",phase2.3dovest,"/post.Rdata", sep=""))
lv3_derb3o<-post

load(file=paste("Model Output/",phase2.5d0,"/post.Rdata", sep=""))
lv5_derb0<-post

load(file=paste("Model Output/",phase2.5dundest,"/post.Rdata", sep=""))
lv5_derb3p<-post

load(file=paste("Model Output/",phase2.5dovest,"/post.Rdata", sep=""))
lv5_derb3o<-post

ModList.st2<-list(
  lv3_derb0=lv3_derb0,
  lv3_derb3p=lv3_derb3p,
  lv3_derb3o=lv3_derb3o,
  lv5_derb0=lv5_derb0,
  lv5_derb3p=lv5_derb3p,
  lv5_derb3o=lv5_derb3o)

ModList.st2.sig<-list(
  lv3_derb0=lv3_derb0,
  lv5_derb0=lv5_derb0)

ModList.st2.db<-list(
  lv3_derb0=lv3_derb0,
  lv3_derb3p=lv3_derb3p,
  lv3_derb3o=lv3_derb3o)

sqrt(exp(-3))
sqrt(exp(-5))
Modnames.st1<-c("mod sig, derby by unbiased",
                "mod sig, wcpue under derby by",
                "mod sig, wcpue over derby by",
                "min sig, derby by unbiased",
                "min sig, wcpue under derby by",
                "min sig, wcpue over derby by") 

Modnames.st2.sig<-c("log var ~ U(-10,-3)",
                    "log var ~ U(-10,-5)")

Modnames.st2.db<-c("wcpue unbiased",
                   "wcpue biased low",
                   "wcpue biased high") 

MR.stage2<-ranker(ModList<-ModList.st2); MR.stage2$test<-"stage2"

write.csv(MR.stage2,
          file=paste("Model Output/Model_ranking_stage2.dev", Sys.Date(),".csv", sep=""))

#*******************************************************************************
#*PHASE 3 COMPARISONS
{Mod.fullname<-c(phase3.mod_PE.r_unif,phase3.mod_PE.r_broad,phase3.mod_PE.r_mod,phase3.mod_PE.r_narrow,
                phase3.min_PE.r_unif,phase3.min_PE.r_broad,phase3.min_PE.r_mod,phase3.min_PE.r_narrow)

load(file=paste("Model Output/",phase3.mod_PE.r_unif,"/post.Rdata", sep=""))
lv3_r_unif<-post

load(file=paste("Model Output/",phase3.mod_PE.r_broad,"/post.Rdata", sep=""))
lv3_r_broad<-post

load(file=paste("Model Output/",phase3.mod_PE.r_mod,"/post.Rdata", sep=""))
lv3_r_mod<-post

load(file=paste("Model Output/",phase3.mod_PE.r_narrow,"/post.Rdata", sep=""))
lv3_r_nar<-post

load(file=paste("Model Output/",phase3.min_PE.r_unif,"/post.Rdata", sep=""))
lv5_r_unif<-post

load(file=paste("Model Output/",phase3.min_PE.r_broad,"/post.Rdata", sep=""))
lv5_r_broad<-post

load(file=paste("Model Output/",phase3.min_PE.r_mod,"/post.Rdata", sep=""))
lv5_r_mod<-post

load(file=paste("Model Output/",phase3.min_PE.r_narrow,"/post.Rdata", sep=""))
lv5_r_nar<-post

ModList.st3<-list(
  lv3_r_unif=lv3_r_unif,
  lv3_r_broad=lv3_r_broad,
  lv3_r_mod=lv3_r_mod,
  lv3_r_nar=lv3_r_nar,
  lv5_r_unif=lv5_r_unif,
  lv5_r_broad=lv5_r_broad,
  lv5_r_mod=lv5_r_mod,
  lv5_r_nar=lv5_r_nar)

ModList.st3.sig<-list(
  )

ModList.st3.r<-list(
  )

Modnames.st3<-c("mod sig, unif R prior",
                "mod sig, broad R prior",
                "mod sig, mod R prior",
                "mod sig, narrow R prior",
                "min sig, unif R prior",
                "min sig, broad R prior",
                "min sig, mod R prior",
                "min sig, narrow R prior") 

Modnames.st3.sig<-c()

Modnames.st3.r<-c() 
}
MR.stage3<-ranker(ModList<-ModList.st3); MR.stage3$test<-"stage3"

write.csv(MR.stage3,
          file=paste("Model Output/Model_ranking_stage3.dev", Sys.Date(),".csv", sep=""))

#!!! CHECK table to determine interactions between sigma and r prior and then decide best way to compare
# graphically.  Hopefully there is no interactions

#****************************************************************************
#* Preferred model simulated data check
{
  phase3_real<-"PHASE3_B1-1_B2-1_upv-5_phmu-0.7_phsig-1.2_Kmu-10.6_Ksig-9.2_derb_0_1500k"
  sim1<-"postsim1"
  sim2<-"postsim2"
  sim3<-"postsim3"
  sim4<-"postsim4"
  sim5<-"postsim5"
  sim6<-"postsim6"
}

{Mod.fullname<-c(phase3_real,#sim1,
                 sim2,sim3,sim4,sim5,sim6,sim7,sim8,sim9,sim10,sim11,sim12)

load(file=paste("Model Output/",phase3_real,"/post.Rdata", sep=""))
Real_Data<-post

load(file=paste("Model Output/",phase3_real,"/postsim1.Rdata", sep=""))
Sim1<-post

load(file=paste("Model Output/",phase3_real,"/postsim2.Rdata", sep=""))
Sim2<-post

load(file=paste("Model Output/",phase3_real,"/postsim3.Rdata", sep=""))
Sim3<-post

load(file=paste("Model Output/",phase3_real,"/postsim4.Rdata", sep=""))
Sim4<-post

load(file=paste("Model Output/",phase3_real,"/postsim5.Rdata", sep=""))
Sim5<-post

load(file=paste("Model Output/",phase3_real,"/postsim6.Rdata", sep=""))
Sim6<-post
load(file=paste("Model Output/",phase3_real,"/postsim7.Rdata", sep=""))
Sim7<-post
load(file=paste("Model Output/",phase3_real,"/postsim8.Rdata", sep=""))
Sim8<-post
load(file=paste("Model Output/",phase3_real,"/postsim9.Rdata", sep=""))
Sim9<-post
load(file=paste("Model Output/",phase3_real,"/postsim10.Rdata", sep=""))
Sim10<-post
load(file=paste("Model Output/",phase3_real,"/postsim11.Rdata", sep=""))
Sim11<-post
load(file=paste("Model Output/",phase3_real,"/postsim12.Rdata", sep=""))
Sim12<-post
load(file=paste("Model Output/",phase3_real,"/postsim13.Rdata", sep=""))
Sim13<-post
load(file=paste("Model Output/",phase3_real,"/postsim14.Rdata", sep=""))
Sim14<-post
load(file=paste("Model Output/",phase3_real,"/postsim15.Rdata", sep=""))
Sim15<-post
}

ModList.sims<-list(
  Real_Data=Real_Data,
  #Sim1=Sim1,
  Sim2=Sim2,
  Sim3=Sim3,
  Sim4=Sim4,
  Sim5=Sim5,Sim6=Sim6,Sim7=Sim7,Sim8=Sim8,Sim9=Sim9,Sim10=Sim10,
  Sim11=Sim11,Sim12=Sim12,Sim13=Sim13,Sim14=Sim14,Sim15=Sim15)

Modnames.sims<-c("True Data",
                #"Sim1",
                "Sim2",
                "Sim3",
                "Sim4",
                "Sim5","Sim6","Sim7","Sim8","Sim9","Sim10",
                "Sim11","Sim12","Sim13","Sim14","Sim15"
                ) 

MR.sims<-ranker(ModList<-ModList.sims); MR.sims$test<-"sims"

write.csv(MR.stage3,
          file=paste("Model Output/Model_ranking_sims.dev", Sys.Date(),".csv", sep=""))

#*******************************************************************************
#* Derby effect in Stage 3

{
  Risk_nobias<-"RISK3_B1-1_B2-1_upv-5_derb_0_recABC_-0perc_1500k"
  Risk_biaslo<-"RISK3_B1-1_B2-1_upv-5_derb_0.3_recABC_-0perc_1500k"
  Risk_biashi<-"RISK3_B1-1_B2-1_upv-5_derb_-0.3_recABC_-0perc_1500k"
}

Mod.fullname<-c(Risk_nobias,Risk_biaslo,Risk_biashi)

load(file=paste("Model Output/",Risk_nobias,"/post.Rdata", sep=""))
Risk_nobias<-post

load(file=paste("Model Output/",Risk_biaslo,"/post.Rdata", sep=""))
Risk_biaslo<-post

load(file=paste("Model Output/",Risk_biashi,"/post.Rdata", sep=""))
Risk_biashi<-post

ModList.derby<-list(
  Risk_nobias=Risk_nobias,
  Risk_biaslo=Risk_biaslo,
  Risk_biashi=Risk_biashi)

Modnames.derby<-c("wcpue ~ bycatch",
                 "wcpue < bycatch",
                 "wcpue > bycatch") 

MR.derby<-ranker(ModList<-ModList.derby); MR.derby$test<-"derby effect"

write.csv(MR.derby,
          file=paste("Model Output/Model_ranking_derby.dev", Sys.Date(),".csv", sep=""))
#*******************************************************************************
#*Look at trace plots if convergence still not quite there

MCMCtrace(beta_uni, params=c("R.hyp" ,"r","Kseo","K","phi","phi.hyp",#"qfCPUE","qsCPUE",
                         "Tau1","Tau2","Tau3","sigma","eta","rB1","rB2"),
          ISB=TRUE, pdf=FALSE, Rhat=TRUE, file="Model Output/Param_trace.pdf")
#===============================================================================
## Compare key parameters of models:
Mods<-ModList.sims #ModList.sims     #ModList.st3.r  ModList.st3.sig
Modnames<-Modnames.sims

dev.off()
filename<-paste("Phase2_derby_bias_comp_",Sys.Date(),sep="")
filename<-paste("Phase3_sigma_comp_",Sys.Date(),sep="")
filename<-paste("Phase3_R_comp_",Sys.Date(),sep="")
filename<-paste("Phase3_allcomp_",Sys.Date(),sep="")
filename<-paste("Phase3_simulations_",Sys.Date(),sep="")

dir.create(paste("Figures/",filename,sep=""))

Kplot<-paramcomp(SEOq="FALSE",Mods=Mods,Modnames=Modnames,param="K", xlow=00,xhi=10000,10000,
          pscale=2,xlabel="K (t)")
Kplot
ggsave(paste0("Figures/",filename,"/K_comp_subs.png"), 
       dpi=300, height=8, width=8, units="in")

Kseo_plot<-paramcomp(SEOq="TRUE",Mods=Mods,Modnames=Modnames,param="Kseo", 100,100000,10000,
          pscale=1.25,xlabel="K (t)")
Kseo_plot
ggsave(paste0("Figures/",filename,"/K_comp.png"), 
       dpi=300, height=5, width=5, units="in")

rplot<-paramcomp(SEOq="FALSE",Mods=Mods,Modnames=Modnames,param="r", -0.005,0.2,10000,
          pscale=2,xlabel="R")
rplot
ggsave(paste0("Figures/",filename,"/r_comp_subs.png"), 
       dpi=300, height=8, width=8, units="in")

Rhypplot<-paramcomp(SEOq="TRUE",Mods=Mods,Modnames=Modnames,param="R.hyp", -0.005,0.2,10000,
          pscale=1.25,xlabel="hyper r")
Rhypplot
ggsave(paste0("Figures/",filename,"/r_comp.png"), 
       dpi=300, height=5, width=5, units="in")

sigmaplot<-paramcomp(SEOq="TRUE",Mods=Mods,Modnames=Modnames,param="sigma", -0.005,0.25,10000,
          pscale=1.25,xlabel="sigma")
sigmaplot
ggsave(paste0("Figures/",filename,"/sigma.png"), 
       dpi=300, height=5, width=5, units="in")

phiplot<-paramcomp(SEOq="FALSE",Mods=Mods,Modnames=Modnames,param="phi", 0,1.15,10000,
          pscale=2,xlabel="phi (1980 biomass:virgin biomass ratio)")
phiplot
ggsave(paste0("Figures/",filename,"/phi_subs.png"), 
       dpi=300, height=8, width=8, units="in")

bigplot<-paramcomp(SEOq="TRUE",Mods=Mods,Modnames=Modnames,param="bigphi", 0,1.15,10000,
                   pscale=2,xlabel="phi (1980 biomass:virgin biomass ratio)")
bigphiplot
ggsave(paste0("Figures/",filename,"/bigphi.png"), 
       dpi=300, height=8, width=8, units="in")

etaplot<-paramcomp(SEOq="TRUE",Mods=Mods,Modnames=Modnames,param="eta", -100,2000,10000,
          pscale=1.5,xlabel="eta")
etaplot
ggsave(paste0("Figures/",filename,"/eta.png"), 
       dpi=300, height=5, width=5, units="in")

paramcomp(SEOq="TRUE",Mods=Mods,Modnames=Modnames,param="rB1", 0,100,10000,
          pscale=2,xlabel="rB1")

paramcomp(SEOq="TRUE",Mods=Mods,Modnames=Modnames,param="rB2", -10,2000,10000,
          pscale=2,xlabel="rB2")

CBtoKplot<-paramcomp(SEOq="FALSE",Mods=Mods,Modnames=Modnames,param="CBtoK", 0,1,10000,
          pscale=2,xlabel="Current biomass:virgin biomass")
CBtoKplot
ggsave(paste0("Figures/",filename,"/CBtoK_subs.png"), 
       dpi=300, height=8, width=8, units="in")

CBtoKseoplot<-paramcomp(SEOq="TRUE",Mods=Mods,Modnames=Modnames,param="CBtoKseo", 0,0.9,10000,
          pscale=1.3,xlabel="Current:virgin biomass")
CBtoKseoplot
ggsave(paste0("Figures/",filename,"/CBtoK.png"), 
       dpi=300, height=5, width=5, units="in")

Fmsyplot<-paramcomp(SEOq="FALSE",Mods=Mods,Modnames=Modnames,param="Fmsy", -0.005,0.1,10000,
          pscale=2,xlabel="Fmsy")
Fmsyplot
ggsave(paste0("Figures/",filename,"/Fmsy_subs.png"), 
       dpi=300, height=8, width=8, units="in")

Fmsyseoplot<-paramcomp(SEOq="TRUE",Mods=Mods,Modnames=Modnames,param="Fmsyseo", -0.005,0.075,10000,
          pscale=1.3,xlabel="Fmsy")
Fmsyseoplot
ggsave(paste0("Figures/",filename,"/Fmsy.png"), 
       dpi=300, height=5, width=5, units="in")

Bmsyplot<-paramcomp(SEOq="FALSE",Mods=Mods,Modnames=Modnames,param="Bmsy", 0,20000,10000,
          pscale=2,xlabel="Bmsy (mt)")
Bmsyplot
ggsave(paste0("Figures/",filename,"/Bmsy_subs.png"), 
       dpi=300, height=8, width=8, units="in")

Bmsyseoplot<-paramcomp(SEOq="TRUE",Mods=Mods,Modnames=Modnames,param="Bmsyseo", 0,80000,10000,
          pscale=1.5,xlabel="Bmsy (mt)")
Bmsyseoplot
ggsave(paste0("Figures/",filename,"/Bmsy.png"), 
       dpi=300, height=5, width=5, units="in")

SSplot<-paramcomp(SEOq="FALSE",Mods=Mods,Modnames=Modnames,param="Stock.Status", 0,2.2,10000,
          pscale=2,xlabel="Stock Status (B2022:B40)")
SSplot
ggsave(paste0("Figures/",filename,"/StockStatus_subs.png"), 
       dpi=300, height=8, width=8, units="in")

SS2plot<-paramcomp(SEOq="TRUE",Mods=Mods,Modnames=Modnames,param="Stock.Status.SEO", 0,2.2,10000,
          pscale=1.25,xlabel="Stock Status (B2022:B40)")
SS2plot
ggsave(paste0("Figures/",filename,"/StockStatus.png"), 
       dpi=300, height=5, width=5, units="in")

Projplot<-paramcomp(SEOq="TRUE",Mods=Mods,Modnames=Modnames,param="Proj1_biomass", 0,2.2,10000,
          pscale=1.25,xlabel="2023 projected biomass (t)")
Projplot
ggsave(paste0("Figures/",filename,"/Proj1_biomass.png"), 
       dpi=300, height=5, width=5, units="in")

piplot<-paramcomp(SEOq="FALSE",Mods=Mods,Modnames=Modnames,param="pi", 0,1,10000,
          pscale=2,xlabel="% SEO K in management unit")
piplot
ggsave(paste0("Figures/",filename,"/pi_subs.png"), 
       dpi=300, height=8, width=8, units="in")

#Group plots
plot_grid(SS2plot,CBtoKseoplot, ncol = 2)


#=================================================================================
# Correlation exam for picking priors...

ModList.r<-list(
  base_e15_r0405=base_e15,
  base_e15_r0410=base_e15_1Tau_r0410,
  base_e15_r0510=base_e15_1Tau_r0510,
  base_e15_r0610=base_e15_1Tau_r0610,
  base_e15_runi=base_e15_1Tau_runi,
  base_e15_1Tau_rg2r50=base_e15_1Tau_rg2r50,
  base_e15_1Tau_rg5r100=base_e15_1Tau_rg5r100,
  base_e15_1Tau_rg9r200=base_e15_1Tau_rg9r200
)

post<-ModList.r[[1]]

rKpairs(ModList.r,1)

for (j in 1:length(ModList.beta)) {
  rKpairs(ModList.beta,j)
}

dev.off()

## by hand here in case want to look at other pairs...
rchains<-MCMCchains(  object=post,  params = "rB1",  excl = NULL,
             ISB = TRUE,  #ignore square brackets
             mcmc.list = FALSE,chain_num = NULL
)

Kchains<-MCMCchains(  object=post,  params = "K",  excl = NULL,
                      ISB = TRUE,  #ignore square brackets
                      mcmc.list = FALSE,chain_num = NULL
)

rK<-cbind(rchains,Kchains); head(rK)

#pairs(data.frame(rK))

betterPairs(data.frame(rK))
mtext(names(ModList.r)[[1]], side=3, padj=-4)

#-------------------------------------------------------------------------------
# parameter est. extraction for table in SAFE

Mod.fullname<-c(phase3.mod_PE.r_unif,phase3.mod_PE.r_broad,phase3.mod_PE.r_mod,phase3.mod_PE.r_narrow,
                phase3.min_PE.r_unif,phase3.min_PE.r_broad,phase3.min_PE.r_mod,phase3.min_PE.r_narrow)
Mod.fullname<-c(Risk_nobias,Risk_biaslo,Risk_biashi)

ModList.st3<-list(
  lv3_r_unif=lv3_r_unif,
  lv3_r_broad=lv3_r_broad,
  lv3_r_mod=lv3_r_mod,
  lv3_r_nar=lv3_r_nar,
  lv5_r_unif=lv5_r_unif,
  lv5_r_broad=lv5_r_broad,
  lv5_r_mod=lv5_r_mod,
  lv5_r_nar=lv5_r_nar)

ModList.derby<-list(
  Risk_nobias=Risk_nobias,
  Risk_biaslo=Risk_biaslo,
  Risk_biashi=Risk_biashi)



param.list<-c(ModList.st3,ModList.derby)

params<-data.frame()

for (m in 1:length(param.list)) {  #m<-1
  params[m,"model"]<-names(param.list[m])
  
  
}


#================================================================================
#  Simulation Exams
Kseo_plot<-paramcomp2(SEOq="TRUE",Mods=Mods,Modnames=Modnames,param="Kseo", 100,100000,10000,
                     pscale=1.25,xlabel="K (t)",vvalue=40367)
Kseo_plot
ggsave(paste0("Figures/",filename,"/K_comp.png"), 
       dpi=300, height=5, width=5, units="in")

rplot<-paramcomp2(SEOq="FALSE",Mods=Mods,Modnames=Modnames,param="r", -0.005,0.2,10000,
                 pscale=3,xlabel="R",vvalue=c(0.016,0.023,0.015,0.013))
rplot
ggsave(paste0("Figures/",filename,"/r_comp_subs.png"), 
       dpi=300, height=8, width=8, units="in")

sigmaplot<-paramcomp2(SEOq="TRUE",Mods=Mods,Modnames=Modnames,param="sigma", -0.005,0.25,10000,
                     pscale=1.25,xlabel="sigma",vvalue=0.058)
sigmaplot
ggsave(paste0("Figures/",filename,"/sigma.png"), 
       dpi=300, height=5, width=5, units="in")

phiplot<-paramcomp2(SEOq="FALSE",Mods=Mods,Modnames=Modnames,param="phi", 0,1.15,10000,
                   pscale=2,xlabel="phi (1980 biomass:virgin biomass ratio)",
                   vvalue=c(0.764,0.55,0.829,0.812))
phiplot
ggsave(paste0("Figures/",filename,"/phi_subs.png"), 
       dpi=300, height=8, width=8, units="in")

piplot<-paramcomp2(SEOq="FALSE",Mods=Mods,Modnames=Modnames,param="pi", 0,1,10000,
                  pscale=2,xlabel="% SEO K in management unit",
                  vvalue=c(0.25,0.115,0.348,0.288))
piplot
ggsave(paste0("Figures/",filename,"/pi_subs.png"), 
       dpi=300, height=8, width=8, units="in")


qplot<-paramcomp2(SEOq="FALSE",Mods=Mods,Modnames=Modnames,param="qsCPUE", 0,2,10000,
                   pscale=2,xlabel="q EYKT",
                   vvalue=c(2.165,10.997,5.068,3.596))   #c(2.165,10.997,5.068,3.596)
qplot
ggsave(paste0("Figures/",filename,"/q_subs.png"), 
       dpi=300, height=8, width=8, units="in")


####''''

CBtoKplot<-paramcomp(SEOq="FALSE",Mods=Mods,Modnames=Modnames,param="CBtoK", 0,1,10000,
                     pscale=2,xlabel="Current biomass:virgin biomass")
CBtoKplot
ggsave(paste0("Figures/",filename,"/CBtoK_subs.png"), 
       dpi=300, height=8, width=8, units="in")

CBtoKseoplot<-paramcomp(SEOq="TRUE",Mods=Mods,Modnames=Modnames,param="CBtoKseo", 0,0.9,10000,
                        pscale=1.3,xlabel="Current:virgin biomass")
CBtoKseoplot
ggsave(paste0("Figures/",filename,"/CBtoK.png"), 
       dpi=300, height=5, width=5, units="in")

Fmsyplot<-paramcomp(SEOq="FALSE",Mods=Mods,Modnames=Modnames,param="Fmsy", -0.005,0.1,10000,
                    pscale=2,xlabel="Fmsy")
Fmsyplot
ggsave(paste0("Figures/",filename,"/Fmsy_subs.png"), 
       dpi=300, height=8, width=8, units="in")

Fmsyseoplot<-paramcomp(SEOq="TRUE",Mods=Mods,Modnames=Modnames,param="Fmsyseo", -0.005,0.075,10000,
                       pscale=1.3,xlabel="Fmsy")
Fmsyseoplot
ggsave(paste0("Figures/",filename,"/Fmsy.png"), 
       dpi=300, height=5, width=5, units="in")

Bmsyplot<-paramcomp(SEOq="FALSE",Mods=Mods,Modnames=Modnames,param="Bmsy", 0,20000,10000,
                    pscale=2,xlabel="Bmsy (mt)")
Bmsyplot
ggsave(paste0("Figures/",filename,"/Bmsy_subs.png"), 
       dpi=300, height=8, width=8, units="in")

Bmsyseoplot<-paramcomp(SEOq="TRUE",Mods=Mods,Modnames=Modnames,param="Bmsyseo", 0,80000,10000,
                       pscale=1.5,xlabel="Bmsy (mt)")
Bmsyseoplot
ggsave(paste0("Figures/",filename,"/Bmsy.png"), 
       dpi=300, height=5, width=5, units="in")

SSplot<-paramcomp(SEOq="FALSE",Mods=Mods,Modnames=Modnames,param="Stock.Status", 0,2.2,10000,
                  pscale=2,xlabel="Stock Status (B2022:B40)")
SSplot
ggsave(paste0("Figures/",filename,"/StockStatus_subs.png"), 
       dpi=300, height=8, width=8, units="in")

SS2plot<-paramcomp(SEOq="TRUE",Mods=Mods,Modnames=Modnames,param="Stock.Status.SEO", 0,2.2,10000,
                   pscale=1.25,xlabel="Stock Status (B2022:B40)")
SS2plot
ggsave(paste0("Figures/",filename,"/StockStatus.png"), 
       dpi=300, height=5, width=5, units="in")

Projplot<-paramcomp(SEOq="TRUE",Mods=Mods,Modnames=Modnames,param="Proj1_biomass", 0,2.2,10000,
                    pscale=1.25,xlabel="2023 projected biomass (t)")
Projplot
ggsave(paste0("Figures/",filename,"/Proj1_biomass.png"), 
       dpi=300, height=5, width=5, units="in")






















