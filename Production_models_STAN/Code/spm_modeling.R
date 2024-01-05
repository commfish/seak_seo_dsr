################################################################################
# State-space surplus production model simulations and model exams
# This code relies on a number of functions in the stan_helper.R script. This script
# provides a simulated data set of biomass estimates and cpue indices of a hypothetical
# fish population.  
# The code then provides functions for formatting the data, bundling initial values
# and prescribing prior values and specifing what data is available (the simulation
# and models provide options for 1 biomass estimate and up to 2 cpue indices). The simulations
# and models also provide for specifying more than one strata to consider. 
# The user may specify the nature of the productivity curve such that Bmsy may occur
# at any value relative to K. 
# The operating and estimation models are Pella-Tomlinson surplus production models
# that allows for estimating process error.  
# Author: Phil Joy
# Created: 1/1/24
################################################################################
{library("rstan")
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  library("bayesplot")
  library("rstanarm")
  library("ggplot2")
  library("shinystan")
  library("hexbin")
  library("dplyr")
  library("ggplot2")
  library(tidyr)
  library(zoo)
  library("wesanderson")}

source("Production_models_STAN/code/stan_helper.R")

set.seed(3480)
#--------------------------------------------------------------------------------
# Decide on the Pella tomlinson curve:
# Where does Pmsy occur relative to virgin biomass/K?
pmsy <- 0.46

pt <- pella_toml(0.46)

pte_p <- pt$pte_p
pte_m <- pt$pte_m

#-------------------------------------------------------------------------------
# Load previous simulation: 
 oldsim<-readRDS("Production_models_STAN/output/TEST.rds")
 list2env(oldsim,.GlobalEnv)
 plot_sim(oldsim)

# ... or SIMULATE a new rockfish data set... : ------------------------------------------------
# Outline of simulation: The simulation is set up with a prescribed F that increases 
# in the early part of the time series to a peak, after which it declines to a stable level.
# The F values are conditioned on r so that r can be varied.The discard schedule is set up 
# with a separate F also linked to the r value and is high and variable early in the
# time series and then decreases to nothing as the fishery becomes full-retention in the last
# portion of the time series. The ratio of discards to catches can be modified with 
# by_to_C so that the fishery could be dominated by bycatch if desired. 
# The simulation is set up for variable numbers of years and strata.
# (increasing the number of strata really slows down model convergence).
# The simulations are set up so that extra uncertainty in biomass estimates and 
# index estimates can be applied. The simulation is also set up to start the population
# at levels below K (depletion) if desired.  Lastly the simulations are set up to 
# apply process error with autocorrelation (rho).  The level of process error is 
# linked to the r value so that greater r produces greater process error.
# The second stage of the simulation (simulate2) allows the user to use only a subset of the 
# simulated biomass and index data to be more in line with what a real data set would
# look like. 

# years and strata
years <- 45
strata <- 3

r<-0.04
K<-4000
rho <- 0.75

tau_b1 <- 0 #extra variance in early part of time series for biomass survey (0 means there is none); 
tau_b2 <- 0 #extra variance in second part of time series for biomass survey 
tau_f1 <- 0 #extra variance in early part of time series for cpue index
tau_f2 <- 0 #extra variance in second part of time series for cpue index
by_to_C <- c(0.2,1)# how does bycatch compare to known catches? This is a rough approximation that scales the 
# the two catches so that you can have different dynamics for simulation testing of the model
# and the nature of the fishery.  A lot of rock fish are bycatch fisheries and this
# can be explored with simulation testing to see if the SPM can handle that scenario...
# This can also be used to increase or decrease the magnitude of the catch... 
depletion <- 0 # how depleted is the stock relative to virgin biomass in time period 1

sim <- simulate (years = years, strata=strata, r=r, K=K, rho=rho,
                     tau_b1 = tau_b1, tau_b2 = tau_b2, 
                     tau_f1 = tau_f1, tau_f2 = tau_f2,
                     depletion = depletion,
                     by_to_C = by_to_C) 
  
list2env(sim,.GlobalEnv)  

plot_sim(sim)
  
# Save the simulated data for later:
saveRDS(sim, "Production_models_STAN/output/TEST.rds")

# 1st step of simulation includes fake data for all years... which isn't likely
# So lets simulate what data is actually available...

data <- simulate2(sim = sim, 
                      bio_srv_frq = 3, # how frequent are biomass surveys? 3 = every third year. 
                  #                      Biomass surveys are staggered so only one strata sampled per year
                      bio_srv_start = 10, # 1st year that biomass survey occurs relative to catch history 
                      cpue_frq = 1, cpue_start = 10, #not staggered... cpue index all strata for each year that is happens... 
                      srv_frq = 1, srv_start = 20,
                      c_cv = 0.05)

list2env(data,.GlobalEnv)  

#-------------------------------------------------------------------------------
# Bundle the data! 

# NOTE: you can change PRIOR values in this list.

# NOTE: this function is set up to accept matrices of abundance set up by year and
# strata with NA's and the prep_data function will vectorize that data for use in 
# Stan.

stan_dat <- prep_data(data = data, # set up to accept data from simulation set up
                      bio_switch = 0, #1 for biomass data, 0 for no biomass data
                      ind1_switch = 1, #1 for cpue index data, 0 for no cpue data
                      ind2_switch = 0, #1 for 2nd cpue index (survey in simulations), 0 for none
                      bio_xV_lastY = 1, #last year that has extra variance.  Only set up for extra variance in one time period that starts in year 0 and goes to this year
                      cpue_xV_lastY = 1,
                      srv_xV_lastY = 1,
                      parameterization = "NC", #"C" for centered parameterization, "NC" for non-centered model
                      pmsy = 0.46, #peak of productivity curve
                      pe_switch = 0, #0 = no process error, 1 = est. process error
                      r_prior_switch = 0, #0=beta, 1=gamma, 2=lognormal
                      rp1 = 1, #1st parameter for r prior
                      rp2 = 1,# 2nd parameter for r prior
                      k1 = 10, #mean for lognormal K parameter
                      k2 = 2, #sd for lognormal K parameter
                      pe_bound = -4, #upper log bound of process error
                      P1_mu = log(1),  # log of the proportion of virgin biomass in year 1
                      P1_sig = 0.01, 
                      isig2_mu = 5,  #prior mu for fit of P
                      isig2_sig = 0.102)  #prior sig for fit of P

#-------------------------------------------------------------------------------
# Establish initial values: 

i1tau2 = list(rep(20,stan_dat$S), # For marginalized catchability models
              c(30,50,30,15),
              rep(30,strata),
              rep(30,strata)) #c(20,40,40,40),
i2tau2 = list(rep(20,stan_dat$S), # For marginalized catchability models
              rep(30,strata),
              rep(15,strata),
              rep(30,strata)) 
isigma2=c(100,140,80,90) # uncertainty around P (proportion of K) estimates
logvar=c(-6,-7,-5,-9) # maximum process error
fsigma2=c(100,100,100,100) # F sigma for F models (not available yet)
medP1 = c(1,1,1,1) # P in year 1.  Value between 0 and 1.

#true_inits <- list(r=r*1.1, K=K*0.9, 
#                   iq1=c(1/q1[1],1/q1[2],1/q1[3],1/q1[4]), 
#                   iq2=c(1/q2[1],1/q2[2],1/q2[3],1/q2[4]), 
#                   isigma2=50, logvar=-5, eps = rep(0,years),
#                   P=P1, C = C+0.000001, pi = c(0.25,0.25,0.25,0.25))

# Bundle initial values.  This function will also produce crude in itial values
# for the P, pi and q parameters
inits <- new_inits(stan_dat = stan_dat,
                  raw_dat = data,
                  chains = 3,
                  r_guess = 0.035,
                  K_guess = 5000,
                  i1tau2 = i1tau2,
                  i2tau2 = i2tau2,
                  isigma2 = isigma2,
                  medP1 = medP1,
                  logvar = logvar,
                  fsigma2 = fsigma2)

#-------------------------------------------------------------------------------
# Run the model!! 

# Model options
# C = centered models
# NC = non-centered model
# margQ = marginalized catchability
# xV = models that accommodate extra variance on biomass and index data
# D = models that accomodate discard data separate from the known catch data
# P1 = models that accomodate biomass being less than K in time period 1.  Note that
#      all centered models accomodate this, but the NC models required extra work
#      to add that functionality. 

# CENTERED MODELS:
# PT_C.stan = basic centered model
# PT_C_margQ.stan = centered model with marginalized catchability
# PT_C_xV.stan = centered model with extra variance capability (not coded yet 12/1/23)
# PT_C_margQ_xV.stan = centered model with extra variance and marginalized catchability 
# PT_C_xV_D.stan = centered model with extra variance and accommodates discard data
# PT_C_margQ_xV_D.stan = centered model with extra variance and marginalized catchability that allows for discards 

# NON_CENTERED MODELS:
# PT_NC.stan = non-centered model
# PT_NC_margQ.stan = non-centered model with marginalized catchability
# PT_NC_xV.stan = non-centered model with extra variance capability 
# PT_NC_margQ_xV.stan = non-centered model with extra variance and marginalized catchability 
# PT_NC_xV_D.stan = non-centered model with extra variance and accommodates discard data
# PT_NC_margQ_xV_D.stan = non-centered model with extra variance and marginalized catchability that accommodates discards 
# PT_NC_margQ_xV_P1.stan = non-centered model with extra variance and marginalized 
#                           catchability that allows biomass in time period 1 to 
#                           be less than K. (NOTE that centered models can accomodate
#                           P1 < K)
# PT_NC_xV_D_P1.stan = ... 
# PT_NC_margQ_xV_D_P1.stan = ...

model = 'PT_NC_xV_D' #'PT_C_margQ_xV.stan' # PT_NC_xV.stan' # PT_NC_margQ_xV.stan'
                       # 'PT_C_xV_D.stan' # 'PT_NC_xV_D.stan'
                       # 'PT_C_margQ_xV_D.stan' # 'PT_NC_margQ_xV_D.stan'


# set iterations and hmc specifics
iters <- 1000000
chains <- 3
inits <- inits #, inits3) #, inits3)
burnin <- 0.4 #proportion of chain length used as warmup/burnin

adapt_delta <- 0.9
stepsize <- 0.1
max_treedepth <- 15
stepsize_jitter <- 0

tstart <- Sys.time()
fit <- stan(file = paste0("Production_models_STAN/Models/",model,".stan"), data = stan_dat, 
            iter = iters, chains = chains, cores=chains,
            init=inits, warmup=burnin*iters, verbose=F, thin=(iters-burnin*iters)/500,
          control = list(adapt_delta = adapt_delta, stepsize = stepsize, 
                         max_treedepth = max_treedepth, stepsize_jitter = stepsize_jitter)) #stepsize_jitter default(0), values between 0 and 1
                                #metric (string, one of "unit_e", "diag_e", "dense_e", defaults to "diag_e")
runtime <- Sys.time() - tstart; runtime

#to save model for later and then reload
saveRDS(fit, paste0("Production_models_STAN/output/",model,"_",iters/1000,"K.rds"))
readRDS(paste0("Production_models_STAN/output/",model,"_",iters/1000,"K.rds"))

# closeAllConnections() #this command is if console stops displaying error messages... 

# quick and dirty diagnostics
check_energy(fit)
check_treedepth(fit)
check_divergences(fit)

get_num_divergent(fit)
get_divergent_iterations(fit)
get_num_max_treedepth(fit) 
get_bfmi(fit) 
get_low_bfmi_chains(fit)

launch_shinystan(fit)
# SKIP DOWN past tuning methods TO PLOT RESULTS

#-------------------------------------------------------------------------------
# Save initial values for future runs...
save_inits <- save_inits(fit = fit, parameterization = "NC", margQ = FALSE)
saveRDS(save_inits,paste0("Production_models_STAN/output/inits.rds"))

tuned_inits<-save_inits #readRDS("Production_models_STAN/output/inits.rds")
# check if you have the right number of chains
str(inits)

iters <- 1000000
# rerun with the new initial values... 
tstart <- Sys.time()
fit2 <- stan(file = paste0("Production_models_STAN/Models/",model,".stan"), data = stan_dat, 
            iter = iters, chains = chains, cores=chains,
            init=tuned_inits, warmup=burnin*iters, verbose=F, thin=(iters-burnin*iters)/500,
            control = list(adapt_delta = 0.9, stepsize = 0.1, max_treedepth = 15, stepsize_jitter = 0))
runtime <- Sys.time() - tstart; runtime

# quick and dirty diagnostics
check_energy(fit)
check_treedepth(fit)
check_divergences(fit)

get_num_divergent(fit)
get_divergent_iterations(fit)
get_num_max_treedepth(fit) 
get_bfmi(fit) 
get_low_bfmi_chains(fit)

launch_shinystan(fit)

# ------------------------------------------------------------------------------
# RESULTS:  
fit <- fit2

samples <- rstan::extract(fit, permuted = TRUE)
N<-years

results <- data.frame()
i <- 1
# model performance: This is set up this way to enable looping over many simulations
# and recording the results... 
{
results[i,"model"] <- model
results[i,"run_date"] <- format(Sys.time(), "%m/%d/%y")
results[i,"iters"] <- iters #<- 10000
#results[i,"metric"] <- mm
results[i,"adapt_delta"] <- adapt_delta
results[i,"stepsize"] <- stepsize
results[i,"max_treedepth"] <- max_treedepth
results[i,"stepsize_jitter"] <- stepsize_jitter
results[i,"isig2mu"] <- stan_dat$isig2_mu #sig_mus[s]
results[i,"isig2_sig"] <- stan_dat$isig2_sig #sig_sigs[ss]
results[i,"r_dist"] <- stan_dat$r_prior_switch #r_dists[r]
results[i,"rp1"] <- stan_dat$rp1
results[i,"rp2"] <- stan_dat$rp2
results[i,"k1"] <- stan_dat$k1
results[i,"k2"] <- stan_dat$k2

# comparison in catches and discards:
results[i,"disc_catch_ratio"] <- median(stan_dat$D_obs)/median(stan_dat$C_obs)

# contrast in biomass/indices
cons <- vector()
if (stan_dat$bio_switch == 1) {
  for (i in 1:stan_dat$S) {
    cons[i] <- max(data$B_ests[i], na.rm = T) / min(data$B_ests[i], na.rm = T)
  }
  results[i,"time_series_contrast"] <- max(cons)
} else {
  for (i in 1:stan_dat$S) {
    cons[i] <- max(data$I1_ests[i], na.rm = T) / min(data$I1_ests[i], na.rm = T)
  }
  results[i,"time_series_contrast"] <- max(cons)
}

# autocorrelation setting
results[i,"true_rho"] <- sim$rho

# model diagnostics
for (j in 1:data$strata) {
  results[i,paste0("energy",j)] <- get_bfmi(fit)[j]
}

results[i,"treedepth"] <- get_num_max_treedepth(fit)
results[i,"prop_divergent"] <- get_num_divergent(fit) / ((1 - burnin) * iters/(((iters-burnin*iters)/500) / chains))
results[i,"run_time"] <- runtime
results[i,"Rhat_gt_1.1"] <- any(summary(fit)$summary[,"Rhat"] > 1.1, na.rm=T)

# simulated values:
results[i,"true_r"] <- sim$r
results[i,"est_r"] <- quantile(samples$r,probs=,c(0.5))
results[i,"bias_r"] <- (quantile(samples$r,probs=,c(0.5)) - sim$r) /  sim$r
results[i,"within_50_r"] <- ifelse(sim$r < quantile(samples$r,probs=c(0.75)) &
                                     sim$r > quantile(samples$r,probs=c(0.25)),
                                   "TRUE","FALSE")
results[i,"within_95_r"] <- ifelse(sim$r < quantile(samples$r,probs=c(0.975)) &
                                     sim$r > quantile(samples$r,probs=c(0.025)),
                                   "TRUE","FALSE")
results[i,"true_K"] <- sim$K
results[i,"est_K"] <- quantile(samples$K,probs=,c(0.5))
results[i,"bias_K"] <- (quantile(samples$K,probs=,c(0.5)) - sim$K) /  sim$K
results[i,"within_50_K"] <- ifelse(sim$K < quantile(samples$K,probs=c(0.75)) &
                                     sim$K > quantile(samples$K,probs=c(0.25)),
                                   "TRUE","FALSE")
results[i,"within_95_K"] <- ifelse(sim$K < quantile(samples$K,probs=c(0.975)) &
                                     sim$K > quantile(samples$K,probs=c(0.025)),
                                   "TRUE","FALSE")


true_msy <- sim$r*sim$K/((sim$pte_p+1)^((sim$pte_p+1)/sim$pte_p))
results[i,"true_msy"] <- true_msy
results[i,"est_msy"] <- quantile(samples$MSY,probs=,c(0.5))
results[i,"bias_msy"] <- (quantile(samples$MSY,probs=,c(0.5)) - true_msy) /  true_msy
results[i,"within_50_msy"] <- ifelse(true_msy < quantile(samples$MSY,probs=c(0.75)) &
                                       true_msy > quantile(samples$MSY,probs=c(0.25)),
                                   "TRUE","FALSE")
results[i,"within_95_msy"] <- ifelse(true_msy < quantile(samples$MSY,probs=c(0.975)) &
                                       true_msy > quantile(samples$MSY,probs=c(0.025)),
                                   "TRUE","FALSE")

true_Bmsy <- stan_dat$msy_frac*sim$K
results[i,"true_Bmsy"] <- true_Bmsy
results[i,"est_Bmsy"] <- quantile(samples$Bmsy,probs=,c(0.5))
results[i,"bias_Bmsy"] <- (quantile(samples$Bmsy,probs=,c(0.5)) - true_Bmsy) /  true_Bmsy
results[i,"within_50_Bmsy"] <- ifelse(true_Bmsy < quantile(samples$Bmsy,probs=c(0.75)) &
                                        true_Bmsy > quantile(samples$Bmsy,probs=c(0.25)),
                                     "TRUE","FALSE")
results[i,"within_95_Bmsy"] <- ifelse(true_Bmsy < quantile(samples$Bmsy,probs=c(0.975)) &
                                        true_Bmsy > quantile(samples$Bmsy,probs=c(0.025)),
                                     "TRUE","FALSE")

true_Fmsy <- (sim$r*sim$K/((sim$pte_p+1)^((sim$pte_p+1)/sim$pte_p))) / (stan_dat$msy_frac*sim$K)
results[i,"true_Fmsy"] <- true_Fmsy
results[i,"est_Fmsy"] <- quantile(samples$Fmsy,probs=,c(0.5))
results[i,"bias_Fmsy"] <- (quantile(samples$Fmsy,probs=,c(0.5)) - true_Fmsy) /  true_Fmsy
results[i,"within_50_Fmsy"] <- ifelse(true_Fmsy < quantile(samples$Fmsy,probs=c(0.75)) &
                                        true_Fmsy > quantile(samples$Fmsy,probs=c(0.25)),
                                      "TRUE","FALSE")
results[i,"within_95_Fmsy"] <- ifelse(true_Fmsy < quantile(samples$Fmsy,probs=c(0.975)) &
                                        true_Fmsy > quantile(samples$Fmsy,probs=c(0.025)),
                                      "TRUE","FALSE")

true_stockstatus <- sum(bio[years,]) / sim$K
results[i,"true_stock_status"] <- true_stockstatus
results[i,"est_stock_status"] <- quantile(samples$Stock_status,probs=,c(0.5))
results[i,"bias_stock_status"] <- (quantile(samples$Stock_status,probs=,c(0.5)) - true_stockstatus) /  true_stockstatus
results[i,"within_50_stock_status"] <- ifelse(true_stockstatus < quantile(samples$Stock_status,probs=c(0.75)) &
                                                true_stockstatus > quantile(samples$Stock_status,probs=c(0.25)),
                                      "TRUE","FALSE")
results[i,"within_95_stock_status"] <- ifelse(true_stockstatus < quantile(samples$Stock_status,probs=c(0.975)) &
                                                true_stockstatus > quantile(samples$Stock_status,probs=c(0.025)),
                                      "TRUE","FALSE")

# strata specific comps
for (j in 1:data$strata) { #j <- 1
  results[i,paste0("true_pi",j)] <- sim$pi[j]
  results[i,paste0("est_pi",j)] <- quantile(samples$pi[,j],probs=c(0.5))
  results[i,paste0("bias_pi",j)] <- (quantile(samples$pi[,j],probs=c(0.5)) - sim$pi[j]) /  sim$pi[j]
  results[i,paste0("within_50_r_pi",j)] <- ifelse(sim$pi[j] < quantile(samples$pi[,j],probs=c(0.75)) &
                                       sim$pi[j] > quantile(samples$pi[j],probs=c(0.25)),
                                     "TRUE","FALSE")
  results[i,paste0("within_95_r_pi",j)] <- ifelse(sim$pi[j] < quantile(samples$pi[,j],probs=c(0.975)) &
                                       sim$pi[j] > quantile(samples$pi[,j],probs=c(0.025)),
                                     "TRUE","FALSE")
  
  if (stan_dat$ind1_switch == 1) {
    results[i,paste0("true_q1_",j)] <- sim$q1[j]
    results[i,paste0("est_q1_",j)] <- quantile(samples$q1[,j],probs=c(0.5))
    results[i,paste0("bias_q1_",j)] <- (quantile(samples$q1[,j],probs=c(0.5)) - sim$q1[j]) /  sim$q1[j]
    results[i,paste0("within_50_r_q1_",j)] <- ifelse(sim$q1[j] < quantile(samples$q1[,j],probs=c(0.75)) &
                                                       sim$q1[j] > quantile(samples$q1[j],probs=c(0.25)),
                                                     "TRUE","FALSE")
    results[i,paste0("within_95_r_q1_",j)] <- ifelse(sim$q1[j] < quantile(samples$q1[,j],probs=c(0.975)) &
                                                       sim$q1[j] > quantile(samples$q1[,j],probs=c(0.025)),
                                                     "TRUE","FALSE")
  }
  
  if (stan_dat$ind2_switch == 1) {
    results[i,paste0("true_q2_",j)] <- sim$q2[j]
    results[i,paste0("est_q2_",j)] <- quantile(samples$q2[,j],probs=c(0.5))
    results[i,paste0("bias_q2_",j)] <- (quantile(samples$q2[,j],probs=c(0.5)) - sim$q2[j]) /  sim$q2[j]
    results[i,paste0("within_50_r_q2_",j)] <- ifelse(sim$q2[j] < quantile(samples$q2[,j],probs=c(0.75)) &
                                                       sim$q2[j] > quantile(samples$q2[j],probs=c(0.25)),
                                                     "TRUE","FALSE")
    results[i,paste0("within_95_r_q2_",j)] <- ifelse(sim$q2[j] < quantile(samples$q2[,j],probs=c(0.975)) &
                                                       sim$q2[j] > quantile(samples$q2[,j],probs=c(0.025)),
                                                     "TRUE","FALSE")
  }
}

# process error comparisons
if (stan_dat$pe_switch == 1) {
  results[i, "process_error"] <- "Estimated"
  quants4<-array(dim=c(3,years))
  for (k in 1:years){
    quants4[,k] <- quantile(samples$eps[,k],probs=,c(0.05,0.5,0.95))
  }
  quants4<-t(quants4)
  quants4<-cbind(quants4,seq(1,years))
  colnames(quants4) <- c("lo90","median","hi90","year")
  quants4<-as.data.frame(quants4)
  
  lm <- lm(quants4$median ~ sim$pe)
  plot(quants4$median ~ sim$pe); abline(v=0,col="red"); abline(h=0,col="red"); abline(lm, col="blue")
  results[i,"pe_rsq_comp"] <- summary(lm)$adj.r.squared
  
  comp <- cbind(quants4$median,sim$pe)
  results[i,"pe_dir_agr"] <- (nrow(comp[comp[,1] > 0 & comp[,2] > 0,]) + nrow(comp[comp[,1] < 0 & comp[,2] < 0,]) )/ sim$years
} else {
  results[i, "process_error"] <- "Not Estimated"
  results[i,"pe_rsq_comp"] <- NA
  results[i,"pe_dir_agr"] <- NA
}

}

results
#-------------------------------------------------------------------------------
# PLOT results

mod_name <- paste0(model,"_",i)
  
plot_biomass(fit=fit,years=years,strata=strata,
                           save=FALSE, 
                           bio_dat = data$B_ests, bio_cv = data$B_cv,
                           sim=TRUE, 
                           sim_bio = bio, mod_name=mod_name,
             units = "mt", year1=1980)
  
plot_p(fit=fit,years=years,strata=strata,
               save=FALSE, 
               bio_dat = data$B_ests, bio_cv = data$B_cv,
               sim=TRUE, 
               sim_bio = bio, mod_name=mod_name,
               units = "mt", year1=1980)

plot_index1(fit=fit,years=years,strata=strata,
            save=FALSE, 
            cpue_dat = data$I1_ests, cpue_cv = data$I1_cv,
            sim=TRUE, 
            sim_cpue = cpue, mod_name=mod_name,
            units = "mt",
            year1 = 1980,
            cpue_name = "Fishery CPUE",
            cpue_type = "kgs/box")

plot_index2(fit=fit,years=years,strata=strata,
            save=FALSE, 
            cpue_dat = data$I2_ests, cpue_cv = data$I2_cv,
            sim=TRUE, 
            sim_cpue = cpue, mod_name=mod_name,
            units = "mt",
            year1 = 1980,
            cpue_name = "Fishery CPUE",
            cpue_type = "kgs/box")

plot_catch(fit=fit,years=years,strata=strata,
                       save=FALSE, 
                       dat = data, 
                       sim=TRUE, 
                       mod_name=mod_name,
                       units = "mt",
                       year1 = 1980) 

plot_discards(fit=fit,years=years,strata=strata,
           save=FALSE, 
           dat = data, 
           sim=TRUE, 
           mod_name=mod_name,
           units = "mt",
           year1 = 1980) 

plot_pe(fit = fit, sim_pe = pe, sim=TRUE, save=FALSE, mod_name=mod_name, year1 = 1980)
  
# POSTERIOR PREDICTIVE CHECKS: ------------------------------------------------- 

ppc_dens_overlay(stan_dat$B_obs, samples$Bio_new)
ppc_ecdf_overlay(stan_dat$B_obs, samples$Bio_new)
ppc_intervals(stan_dat$B_obs, samples$Bio_new)
ppc_scatter_avg(log(stan_dat$B_obs), log(samples$Bio_new))

ppc_dens_overlay(stan_dat$I1_obs, samples$I1_new)
ppc_ecdf_overlay(stan_dat$I1_obs, samples$I1_new)
ppc_intervals(stan_dat$I1_obs, samples$I1_new)
ppc_scatter_avg(log(stan_dat$I1_obs), log(samples$I1_new))

ppc_dens_overlay(stan_dat$I2_obs, samples$I2_new)
ppc_ecdf_overlay(stan_dat$I2_obs, samples$I2_new)
ppc_intervals(stan_dat$I2_obs, samples$I2_new)
ppc_scatter_avg(log(stan_dat$I2_obs), log(samples$I2_new))

# more fun DIAGNOSTICS!!... ----------------------------------------------------

posterior <- as.array(fit)
np <- nuts_params(fit)

params <- c("r", "K","isigma2", 
            "i1tau2[1]","i1tau2[2]","i1tau2[3]",
            "pi[1]","pi[2]","pi[3]")

mcmc_pairs(posterior, pars = params)#,"iq","isigma2"))

mcmc_pairs(
  posterior,
  pars = params,
  #transformations = list(sigma = "log"), # show log(sigma) instead of sigma
  #off_diag_fun = "hex" # use hexagonal heatmaps instead of scatterplots
  condition = pairs_condition(nuts = "accept_stat__"),
  np = np
)

bayesplot::color_scheme_set("brightblue")
mcmc_pairs(
  fit,
  pars = params,#"pi[1]","pi[2]","pi[3]","pi[4]"),
  #transformations = list(sigma = "log"),
  condition = pairs_condition(nuts = "accept_stat__"),
  off_diag_args = list(size = 3/4, alpha = 1/3), # size and transparency of scatterplot points
  np_style = pairs_style_np(div_color = "green", div_shape = 2, div_size=0.5), # color and shape of the divergences
  np = np
)

pairs(
  fit,
  pars = params,
  condition = pairs_condition(nuts = "divergent__")
)

mcmc_pairs(
  fit,
  pars = params,
  diag_fun = "dens",
  off_diag_fun = "hex",
  np_style = pairs_style_np(div_color = "green", div_shape = 2)
)

pairs(
  fit,
  pars = params,
  #diag_fun = "dens",
  #off_diag_fun = "hex",
  #np_style = pairs_style_np(div_color = "green", div_shape = 2), # color and shape of the divergences
  condition = pairs_condition(nuts = "divergent__")
)


#-------------------------------------------------------------------------------
# sigma prior associated with divergences

log(4000)
log(3500)
log(3000); log(5000)

ks <- rlnorm(100000,10,2); hist(ks[ks<10000], breaks = 1000)

par(mfrow=c(1,1))
sigs<-rgamma(100000,3.785,0.102) #; hist(sigs, breaks = 1000)
sigs0<-rgamma(100000,3.785,0.0102)
#hist(samples$isigma2, breaks = 1000)
plot(density(sigs)); lines(density(samples$isigma2), col="blue"); lines(density(sigs0), col="forestgreen")
lines(density(samples$isigma2), col = "darkcyan")
sigs2<-rgamma(100000,7,0.1); lines(density(sigs2), col="red")
sigs2a<-rgamma(100000,7,0.075); lines(density(sigs2a), col="orange")
sigs3<-rgamma(100000,5.5,0.1); lines(density(sigs3), col="violet")
sigs4<-rgamma(100000,4.5,0.1); lines(density(sigs4), col="purple")
beta<-rbeta(10000,1,1); lines(density(beta), col="aquamarine4")

#notes 11/9:
# gamma(3.785,0.102) gets us down to <5% deviants with uniform on K
# need to try lognormal on K with above gamma: GOOD! still <5% divergent
# uninformative gamma(0.5,0.001): Boo!! ~60% divergent! smaller step sizes noted
# uninformative gamma(0.5,0.001) T[0.01, 1000]: Boo!! ~60% divergent! smaller step sizes noted
# need to try better fitting to posterior: gamma(7,0.1):  Pretty good! ~10% divergent
# gamma(5.5,0.1): ~5% divergent
# gamma(4.5,0.1): ~5% divergent
# gamma(7,0.075): 30% divergent.  Blech! 
# Note: posterior always seems to shift right from the prior, even when the prior
#       is placed right on top of the previous posterior.
# No sigma prior: 67% divergent YUCK!!!
# 

# Other notes: Truncating P[] is currently [0.05,1.6].  It would seem that 
# truncating down to 0.01 would allow for a full suite of possibilities but
# model currently crashes when I try that.  Not sure why that is at this point,
# but something to work on. 

ps <- rlnorm(10000,log(0.7),0.1); plot(density(ps))
