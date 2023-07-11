##############################################################################
## Code to examine the results of simulations from Data_simulator.R and 
## Model_simulated_data.R
## January 2023
## P. Joy
## Code order to conduct simulations
## identify which model you are working with and make sure consistent between scripts!
## 1) Simmulate data to look at and save values and plots: Data_simulator.R
## 2) Run the model on the simulated data (looooong time): Model_simulated_data.R
## 3) Examine simulation results: Sim_exam.R
##
##############################################################################
{
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
  library(scales)
  library(tidyverse)
  library(ggridges)
  
  source("Code/2022_DSR_SAFE_models/Phase1/DATALOAD_SEO_YE_SPM_Func_1980.R")
  source("Code/2022_DSR_SAFE_models/Phase3/DATAPREP_SPM_1980_PHASE3.R")
  source("Code/2022_DSR_SAFE_models/Phase1/PLOT_SPM80.R")
  source("Code/Posterior_Plotting/YE_SPM_posterior_exams_Func.R")
}

### identify which model results you're working with.  This should match the other
# 2 scripts!!!
res.to.sim<-"PHASE3_B1-1_B2-1_upv-5_phmu-0.7_phsig-1.2_Kmu-10.6_Ksig-9.2_derb_0_1500k"

### Pull up parameter values used in the simulation:

Sim.params<-read.csv(paste0("Model Output/",res.to.sim,"/simulations/simulated_values.csv"))
str(Sim.params)

##Load simulation results... 

a<-read.csv(paste0("Model Output/",res.to.sim,"/simulations/results_sims_1-50.csv"))
b<-read.csv(paste0("Model Output/",res.to.sim,"/simulations/results_sims_51-100.csv"))
sims<-rbind(a,b)
str(sims)
##Questions to answer: did model correctly estimate the parameters?!

#r---------------------------------------------------------------------------------
true.r<-Sim.params$r
sims$parameter
sim.r<-sims %>% filter(parameter == "r[1]" | parameter == "r[2]" |
                         parameter == "r[3]" | parameter == "r[4]")
rs<-unique(sim.r$parameter)

r.comps<-data.frame()

for (i in 1:length(rs)) {  #i<-1
  area<-sim.r[sim.r$parameter == rs[i],]
  true<-true.r[i]
  
  plot(density(area$mean))
  lines(density(area$X50.), col="blue")
  abline(v=true, col="red")
  
  r.comps[i,"parameter"]<-"r"
  
  r.comps[i,"true"]<-true
  r.comps[i,"CI_contains_true"]<-sum(area$X2.5.< true & area$X97.5. > true, na.rm=TRUE)/nrow(area)
  
  r.comps[i,"mean_sim.mean"]<-mean(area$mean)
  r.comps[i,"mean_sim.med"]<-mean(area$X50.)
  r.comps[i,"mean_dif"]<-r.comps[i,"mean_sim.mean"]-r.comps[i,"true"]
  r.comps[i,"med_dif"]<-r.comps[i,"mean_sim.med"]-r.comps[i,"true"]
  
  r.comps[i,"mean.bias"]<-r.comps[i,"mean_dif"]/true
  r.comps[i,"med.bias"]<-r.comps[i,"med_dif"]/true
  
}

#K -----------------------------------------------------------------------------
true.K<-exp(Sim.params$logKseo)*Sim.params$pi

sim.K<-sims %>% filter(parameter == "K[1]" | parameter == "K[2]" |
                         parameter == "K[3]" | parameter == "K[4]")

Ks<-unique(sim.K$parameter)

K.comps<-data.frame()

for (i in 1:length(Ks)) {  #i<-1
  area<-sim.K[sim.K$parameter == Ks[i],]
  true<-true.K[i]
  
  plot(density(area$mean))
  lines(density(area$X50.), col="blue")
  abline(v=true, col="red")
  
  K.comps[i,"parameter"]<-"K"
  K.comps[i,"true"]<-true
  K.comps[i,"CI_contains_true"]<-sum(area$X2.5.< true & area$X97.5. > true, na.rm=TRUE)/nrow(area)
  
  K.comps[i,"mean_sim.mean"]<-mean(area$mean)
  K.comps[i,"mean_sim.med"]<-mean(area$X50.)
  K.comps[i,"mean_dif"]<-K.comps[i,"mean_sim.mean"]-K.comps[i,"true"]
  K.comps[i,"med_dif"]<-K.comps[i,"mean_sim.med"]-K.comps[i,"true"]
  
  K.comps[i,"mean.bias"]<-K.comps[i,"mean_dif"]/true
  K.comps[i,"med.bias"]<-K.comps[i,"med_dif"]/true
  
}

#Fmsy, Bmsy, MSY --------------------------------------------------------------------------
#MSY[i]<-r[i]*K[i]/((p+1)^((p+1)/p)) #r*K/4 for Schaefer
#Bmsy[i]<-0.4*K[i]  #0.5 for Schaefer
#Fmsy[i]<-MSY[i]/Bmsy[i]

p <- 0.18815
true.MSY<-true.r*true.K/((p+1)^((p+1)/p))
true.Bmsy<-0.4*true.K
true.F<-(true.r*true.K/((p+1)^((p+1)/p)))/(0.4*true.K)

unique(sims$parameter)

sim.F<-sims %>% filter(parameter == "Fmsy[1]" | parameter == "Fmsy[2]" |
                         parameter == "Fmsy[3]" | parameter == "Fmsy[4]")

Fs<-unique(sim.F$parameter)

F.comps<-data.frame()

for (i in 1:length(Fs)) {  #i<-1
  area<-sim.F[sim.F$parameter == Fs[i],]
  true<-true.F[i]
  
  F.comps[i,"parameter"]<-"Fmsy"
  F.comps[i,"true"]<-true
  F.comps[i,"CI_contains_true"]<-sum(area$X2.5.< true & area$X97.5. > true, na.rm=TRUE)/nrow(area)
  
  F.comps[i,"mean_sim.mean"]<-mean(area$mean)
  F.comps[i,"mean_sim.med"]<-mean(area$X50.)
  F.comps[i,"mean_dif"]<-F.comps[i,"mean_sim.mean"]-F.comps[i,"true"]
  F.comps[i,"med_dif"]<-F.comps[i,"mean_sim.med"]-F.comps[i,"true"]
  
  F.comps[i,"mean.bias"]<-F.comps[i,"mean_dif"]/true
  F.comps[i,"med.bias"]<-F.comps[i,"med_dif"]/true
}

sim.B<-sims %>% filter(parameter == "Bmsy[1]" | parameter == "Bmsy[2]" |
                         parameter == "Bmsy[3]" | parameter == "Bmsy[4]")

Bs<-unique(sim.B$parameter)

B.comps<-data.frame()

for (i in 1:length(Bs)) {  #i<-1
  area<-sim.B[sim.B$parameter == Bs[i],]
  true<-true.Bmsy[i]
  
  B.comps[i,"parameter"]<-"Bmsy"
  B.comps[i,"true"]<-true
  B.comps[i,"CI_contains_true"]<-sum(area$X2.5.< true & area$X97.5. > true, na.rm=TRUE)/nrow(area)
  
  B.comps[i,"mean_sim.mean"]<-mean(area$mean)
  B.comps[i,"mean_sim.med"]<-mean(area$X50.)
  B.comps[i,"mean_dif"]<-B.comps[i,"mean_sim.mean"]-B.comps[i,"true"]
  B.comps[i,"med_dif"]<-B.comps[i,"mean_sim.med"]-B.comps[i,"true"]
  
  B.comps[i,"mean.bias"]<-B.comps[i,"mean_dif"]/true
  B.comps[i,"med.bias"]<-B.comps[i,"med_dif"]/true
}
#unique(sims$parameter)
sim.MSY<-sims %>% filter(parameter == "MSY[1]" | parameter == "MSY[2]" |
                         parameter == "MSY[3]" | parameter == "MSY[4]")

MSYs<-unique(sim.MSY$parameter)

MSY.comps<-data.frame()

for (i in 1:length(MSYs)) {  #i<-1
  area<-sim.MSY[sim.MSY$parameter == MSYs[i],]
  true<-true.MSY[i]
  
  MSY.comps[i,"parameter"]<-"MSY"
  MSY.comps[i,"true"]<-true
  MSY.comps[i,"CI_contains_true"]<-sum(area$X2.5.< true & area$X97.5. > true, na.rm=TRUE)/nrow(area)
  
  MSY.comps[i,"mean_sim.mean"]<-mean(area$mean)
  MSY.comps[i,"mean_sim.med"]<-mean(area$X50.)
  MSY.comps[i,"mean_dif"]<-MSY.comps[i,"mean_sim.mean"]-MSY.comps[i,"true"]
  MSY.comps[i,"med_dif"]<-MSY.comps[i,"mean_sim.med"]-MSY.comps[i,"true"]
  
  MSY.comps[i,"mean.bias"]<-MSY.comps[i,"mean_dif"]/true
  MSY.comps[i,"med.bias"]<-MSY.comps[i,"med_dif"]/true
}

# Phi-------------------------------------------------------------------------
true.phi<-Sim.params$phi
sims$parameter
sim.phi<-sims %>% filter(parameter == "phi[1]" | parameter == "phi[2]" |
                         parameter == "phi[3]" | parameter == "phi[4]")
phis<-unique(sim.phi$parameter)

phi.comps<-data.frame()

for (i in 1:length(phis)) {  #i<-1
  area<-sim.phi[sim.phi$parameter == phis[i],]
  true<-true.phi[i]
  
  phi.comps[i,"parameter"]<-"phi"
  
  phi.comps[i,"true"]<-true
  phi.comps[i,"CI_contains_true"]<-sum(area$X2.5.< true & area$X97.5. > true, na.rm=TRUE)/nrow(area)
  
  phi.comps[i,"mean_sim.mean"]<-mean(area$mean)
  phi.comps[i,"mean_sim.med"]<-mean(area$X50.)
  phi.comps[i,"mean_dif"]<-phi.comps[i,"mean_sim.mean"]-phi.comps[i,"true"]
  phi.comps[i,"med_dif"]<-phi.comps[i,"mean_sim.med"]-phi.comps[i,"true"]
  
  phi.comps[i,"mean.bias"]<-phi.comps[i,"mean_dif"]/true
  phi.comps[i,"med.bias"]<-phi.comps[i,"med_dif"]/true
  
}

# Pi-------------------------------------------------------------------------
true.pi<-Sim.params$pi
sims$parameter
sim.pi<-sims %>% filter(parameter == "pi[1]" | parameter == "pi[2]" |
                           parameter == "pi[3]" | parameter == "pi[4]")
pis<-unique(sim.pi$parameter)

pi.comps<-data.frame()

for (i in 1:length(pis)) {  #i<-1
  area<-sim.pi[sim.pi$parameter == pis[i],]
  true<-true.pi[i]
  
  pi.comps[i,"parameter"]<-"pi"
  
  pi.comps[i,"true"]<-true
  pi.comps[i,"CI_contains_true"]<-sum(area$X2.5.< true & area$X97.5. > true, na.rm=TRUE)/nrow(area)
  
  pi.comps[i,"mean_sim.mean"]<-mean(area$mean)
  pi.comps[i,"mean_sim.med"]<-mean(area$X50.)
  pi.comps[i,"mean_dif"]<-pi.comps[i,"mean_sim.mean"]-pi.comps[i,"true"]
  pi.comps[i,"med_dif"]<-pi.comps[i,"mean_sim.med"]-pi.comps[i,"true"]
  
  pi.comps[i,"mean.bias"]<-pi.comps[i,"mean_dif"]/true
  pi.comps[i,"med.bias"]<-pi.comps[i,"med_dif"]/true
  
}

# Tau1-------------------------------------------------------------------------
true.tau1<-Sim.params$Tau1
#sims$parameter
sim.tau1<-sims %>% filter(parameter == "Tau1") 
tau1s<-unique(sim.tau1$parameter)

tau1.comps<-data.frame()
i<-1
  area<-sim.tau1[sim.tau1$parameter == tau1s[i],]
  true<-true.tau1[i]
  
  tau1.comps[i,"parameter"]<-"Tau1"
  
  tau1.comps[i,"true"]<-true
  tau1.comps[i,"CI_contains_true"]<-sum(area$X2.5.< true & area$X97.5. > true, na.rm=TRUE)/nrow(area)
  
  tau1.comps[i,"mean_sim.mean"]<-mean(area$mean)
  tau1.comps[i,"mean_sim.med"]<-mean(area$X50.)
  tau1.comps[i,"mean_dif"]<-tau1.comps[i,"mean_sim.mean"]-tau1.comps[i,"true"]
  tau1.comps[i,"med_dif"]<-tau1.comps[i,"mean_sim.med"]-tau1.comps[i,"true"]
  
  tau1.comps[i,"mean.bias"]<-tau1.comps[i,"mean_dif"]/true
  tau1.comps[i,"med.bias"]<-tau1.comps[i,"med_dif"]/true

#-Tau3--------------------------------------------------------------------------
  true.tau3<-Sim.params$Tau3
  #sims$parameter
  sim.tau3<-sims %>% filter(parameter == "Tau3") 
  tau3s<-unique(sim.tau3$parameter)
  
  tau3.comps<-data.frame()
  i<-1
  area<-sim.tau3[sim.tau3$parameter == tau3s[i],]
  true<-true.tau3[i]
  
  tau3.comps[i,"parameter"]<-"Tau3"
  
  tau3.comps[i,"true"]<-true
  tau3.comps[i,"CI_contains_true"]<-sum(area$X2.5.< true & area$X97.5. > true, na.rm=TRUE)/nrow(area)
  
  tau3.comps[i,"mean_sim.mean"]<-mean(area$mean)
  tau3.comps[i,"mean_sim.med"]<-mean(area$X50.)
  tau3.comps[i,"mean_dif"]<-tau3.comps[i,"mean_sim.mean"]-tau3.comps[i,"true"]
  tau3.comps[i,"med_dif"]<-tau3.comps[i,"mean_sim.med"]-tau3.comps[i,"true"]
  
  tau3.comps[i,"mean.bias"]<-tau3.comps[i,"mean_dif"]/true
  tau3.comps[i,"med.bias"]<-tau3.comps[i,"med_dif"]/true
  
  #-sigma--------------------------------------------------------------------------
  true.sigma<-sqrt(exp(Sim.params$logvar))
  #sims$parameter
  sim.sigma<-sims %>% filter(parameter == "sigma") 
  sigmas<-unique(sim.sigma$parameter)
  
  sigma.comps<-data.frame()
  i<-1
  area<-sim.sigma[sim.sigma$parameter == sigmas[i],]
  true<-true.sigma[i]
  
  sigma.comps[i,"parameter"]<-"sigma"
  
  sigma.comps[i,"true"]<-true
  sigma.comps[i,"CI_contains_true"]<-sum(area$X2.5.< true & area$X97.5. > true, na.rm=TRUE)/nrow(area)
  
  sigma.comps[i,"mean_sim.mean"]<-mean(area$mean)
  sigma.comps[i,"mean_sim.med"]<-mean(area$X50.)
  sigma.comps[i,"mean_dif"]<-sigma.comps[i,"mean_sim.mean"]-sigma.comps[i,"true"]
  sigma.comps[i,"med_dif"]<-sigma.comps[i,"mean_sim.med"]-sigma.comps[i,"true"]
  
  sigma.comps[i,"mean.bias"]<-sigma.comps[i,"mean_dif"]/true
  sigma.comps[i,"med.bias"]<-sigma.comps[i,"med_dif"]/true
#------------------------------------------------------------------------------
# Put together and save results...

sim.sum<-rbind(r.comps,K.comps,F.comps,B.comps,MSY.comps, phi.comps,pi.comps,
               tau1.comps,tau3.comps, sigma.comps)
sim.sum[,2]<-as.numeric(sim.sum[,2])
sim.sum[,5]<-as.numeric(sim.sum[,5])
sim.sum[,6]<-as.numeric(sim.sum[,6])
sim.sum[,7]<-as.numeric(sim.sum[,7])

write.csv(sim.sum, file=paste0("Model Output/",res.to.sim,"/simulations/sim_results_summary.csv"))

#-------------------------------------------------------------------------------
# Any better at the level of the entire SEO?
true.Kseo<-exp(Sim.params$logKseo)
true.R<-Sim.params$R  #<- didn't model with hyper prior? maybe that's the problem







