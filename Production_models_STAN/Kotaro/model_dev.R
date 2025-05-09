rm(list=ls())
gc()

library(MASS) ## needed for the multivariate normal distribution
library(dplyr)
library(tidyverse)
library(wesanderson)

source("Production_models_STAN/Code/stan_helper.R")

Niter = 1000   # number of simulation iterations
Narea = 3      # number of sub-area


Catch_type = "F"    # "MSY", or "F"

#for (i in 1:Niter){
i = 123

### Setting up the OM 

set.seed(i)
rs = runif(Narea, 0.045, 0.055)    # variability in rs
Ks = runif(Narea, 100000, 200000)
Nyear = 40
Year_start = 1   # when we have survey/CPUE info from; ignore for now and keep at 1. Dealing with it a little different...  
p = 0.18815      # Pella Tom productivity curve exponent
var_ratio = 1  # ration proc / obs error variance
sigE = 0.05
varE = rep(sigE^2, Narea)
varO = var_ratio*varE
sigO = sqrt(varO)


B <- matrix(0, nrow=Nyear, ncol=Narea)
P <- matrix(0, nrow=Nyear, ncol=Narea)
logB <- matrix(0, nrow=Nyear, ncol=Narea)
B[1,] = Ks
P[1,] = rep(1, Narea)
MSYs =  rs * Ks/(p+1)^((p+1)/p)
HMSYs = MSYs/Ks
Hmax = 3
Catch <- matrix(0, nrow=Nyear-1, ncol=Narea)

# Simulate Catch (two-way trip right now)

  ## Catch is independent of stock size 
    if (Catch_type == "MSY") val = sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-10), seq(Hmax, 0.5, length.out=10))*MSYs[x], 0.05*MSYs[x]))
    
  ## Catch is dependent of stock size: HMSy * scaler describing fishing pressure over time... 0.5HMSY to 3*HMSYs then down again.
    if (Catch_type == "F")  H = sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-10), seq(Hmax, 0.5, length.out=10))*HMSYs[x], 0.1*HMSYs[x]))

# Need to create a time series of process error
  epsilon = mvrnorm(Nyear, -0.5*varE, diag(Narea)*varE)

# Update the population dynamics
  for (yr in 1:(Nyear-1)){
    for (sp in 1:Narea){
      if (Catch_type == "MSY") 
      {
        B[yr+1,sp] = (max( B[yr,sp] + rs[sp]*B[yr,sp]/p*(1 - (B[yr,sp]/Ks[sp])^p) - val[yr,sp], 1))*exp(epsilon[yr,sp])
        Catch[yr,sp] = val[yr,sp]
      }
      if (Catch_type == "F") 
      {
        B[yr+1,sp] = (max( B[yr,sp] + rs[sp]*B[yr,sp]/p*(1 - (B[yr,sp]/Ks[sp])^p) - B[yr,sp]*H[yr,sp], 1))*exp(epsilon[yr,sp])
        Catch[yr,sp] = B[yr,sp]*H[yr,sp]
      }
      # P[yr+1,sp] = fmax(P[yr,sp] + (rs[sp] / p) * P[yr,sp] * (1 - pow(P[yr,sp] , p)) - Catch[yr,sp] / K[sp], 0.001)*exp(epsilon[yr,sp])
    }}

  matplot(1:Nyear, B, type="l")
  

### Now trying to create some abundance indices
  qs <- runif(Narea, 0.00001, 0.0001)
 
  epsilonO = mvrnorm(Nyear, -0.5*varO, diag(Narea)*varO)

  IAs = sim_IAs = B
  for (yr in 1: (Nyear)){
    for (sp in 1:Narea){
      IAs[yr,sp] <- qs[sp]*B[yr,sp]*exp(epsilonO[yr,sp])
      sim_IAs[yr,sp] <- qs[sp]*B[yr,sp]
    }
  }

  matplot(1:Nyear, IAs, type="l")
  
  ## Add in biomass estimates:
  sigB = 0.15
  varB = rep(sigB^2, Narea)
  epsilonB = mvrnorm(Nyear, varB, diag(Narea)*varB)
  
  BEs = BElo = BEhi = B
  
  for (yr in 1: (Nyear)){
    for (sp in 1:Narea){
      BEs[yr,sp] <- B[yr,sp]*exp(epsilonB[yr,sp])
      BElo[yr,sp] <- exp(log(BEs[yr,sp]) - log(1.96 * sigB * BEs[yr,sp]) ) 
      BEhi[yr,sp] <- exp(log(BEs[yr,sp]) + log(1.96 * sigB * BEs[yr,sp]) ) 
    }
  }
  
  # But each area is only visited every third year starting in the 11th year of the fishery:
  bio_srv_frq <- 3
  bio_srv_start <- 11
  
  bio_est <- data.frame()
  for (k in 1:Nyear){
    bio_est[k,c(1:Narea)] <- NA
  }
  for (k in 1:Narea){
    bio_est[seq(bio_srv_start+k,Nyear,bio_srv_frq),k] <- 
      BEs[seq(bio_srv_start+k,Nyear,bio_srv_frq),k]
    bio_est_hi <- BEhi[seq(bio_srv_start+k,Nyear,bio_srv_frq),k]
    bio_est_lo <- BElo[seq(bio_srv_start+k,Nyear,bio_srv_frq),k]
  }
  
  B_cvs <- ifelse(is.na(bio_est),NA,sigB)
  B_ests <- bio_est
  
  #need to vectorize the biomass estimates:B_obs <- as.vector(data$B_ests[!is.na(data$B_ests)])
  B_obs <- as.vector(B_ests[!is.na(B_ests)])
  B_cv <- as.vector(B_cvs[!is.na(B_cvs)])
  N_Bobs <- sum (!is.na(B_ests))
  B_pos <- list()
  for (i in 1:Narea) {
    B_pos[[i]] <- which(!is.na(B_ests[,i]))
  }
  B_pos <- as.vector(unlist(B_pos))
  S_Bobs <- vector()
  for (i in 1:Narea) {
    S_Bobs[i] <- length(B_ests[,i][!is.na(B_ests[,i])])
  }

### Now making the abundance indices not available for the first years & calculating total catch
  Catch_all <- apply(Catch, 1, sum)
  B_all <- apply(B, 1, sum)
  K_all <- sum(Ks)
  
  #Year_start <- 1
  IA_obs <- I_obs_plot <- IAs[Year_start:Nyear,] 
  matplot(IA_obs, type="l")
  Btrue <- B[Year_start:Nyear,]
  matplot(Btrue, type="l")
  # matplot(Catch[Year_start:Nyear,], type="l")
  Btrue[1,]/Ks
  
  par(mfrow=c(4,1), mar=c(2,4,1,1))
  matplot(Btrue, type="l")
  matplot(epsilon[Year_start:Nyear,], type="l")
  matplot(IA_obs, type="l")
  matplot(Catch[Year_start:(Nyear-1),], type="l")

  #Eventually vectorize the index...
# START of CPUE DATAtry getting rid of first five years of index data:
  cpue_start <- 6 #models behaved well! 
  cpue_start <- 11
  
  cpue2 <- data.frame()
  for (i in 1:Nyear){
    cpue2[i,c(1:Narea)] <- NA
  }
  for (i in 1:Narea){
    cpue2[seq(cpue_start,Nyear,1),i] <- IA_obs[seq(cpue_start,Nyear,1),i]
  }
  
  #I_obs <- as.vector(IA_obs[!is.na(IA_obs)]) 
  I_obs <- as.vector(cpue2[!is.na(cpue2)])
  I_cv <- I_cv_plot <- ifelse(is.na(cpue2),NA,sigO) 
  I_cv <- as.vector(I_cv[!is.na(I_cv)])
  N_Iobs <- sum (!is.na(I_obs))
  I_pos <- list()
  for (i in 1:Narea) {
    I_pos[[i]] <- which(!is.na(cpue2[,i]))
  }
  I_pos <- as.vector(unlist(I_pos))
  S_Iobs <- vector()
  for (i in 1:Narea) {
    S_Iobs[i] <- length(cpue2[,i][!is.na(cpue2[,i])])
  }
  

### If using STan

library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


## Play around with informative priors, start at true, add more noisiness
if (Narea == 1){
  data = list(Npre = length(Year_start:Nyear), 
              S=Narea, 
              p = 0.18815,
              C_obs = as.vector(Catch[Year_start:(Nyear-1),]),
              I_obs=as.vector(IA_obs), 
              ratio = 1,
              N_Bobs = N_Bobs,
              B_obs = B_obs,
              B_cv = B_cv,
              S_Bobs = S_Bobs,
              B_pos = B_pos)
  
  init_create <- function(chain_id=1){
    set.seed(chain_id+123)
    inits <- list(
      # MSY = 15000*apply( data$C_obs, 1, mean)/apply( data$C_obs, 1, mean)[1],
      # MSY = MSYs,
      r = runif(data$S, 0.01, 0.1),
      # PP_init = runif(4, 0.8, 1),
      PP_init = runif(data$S, 0.99, 1), # will assume time series starts from unfished biomass
      # K = 300000 *apply( data$C_obs, 1, mean)/apply( data$C_obs, 1, mean)[1],
      K = Ks,
      # logqs = runif(4, -10, -6),
      iqs = 1/qs,
      # sigma_proc = runif(1, 0.001, 0.1)
      # sigma_proc = sigE
      isigma2 = (1/sigE)^2
    )
    # inits$PE = t(mvrnorm(19, rep(-0.5*inits$sigma_proc^2, data$S), diag(data$S)*inits$sigma_proc^2) )
    # inits$PE = matrix(0, nrow=data$S, ncol=Nyear - Year_start)
    inits$PE = rep(0, Nyear - Year_start)
    # inits$PP_init = inits$PP_init/sum(inits$PP_init)
    return(inits)
  }
  
  
}

if (Narea > 1){
  data = list(Npre = length(Year_start:Nyear), 
              S=Narea, 
              p = 0.18815,
              C_obs = t(Catch[Year_start:(Nyear-1),]),
              I_obs=t(IA_obs), 
              ratio = 1,
              N_Bobs = N_Bobs,
              B_obs = B_obs,
              B_cv = B_cv,
              S_Bobs = S_Bobs,
              B_pos = B_pos)
  
  init_create <- function(chain_id=1){
    set.seed(chain_id+123)
    inits <- list(
    # MSY = 15000*apply( data$C_obs, 1, mean)/apply( data$C_obs, 1, mean)[1],
    # MSY = MSYs,
    r = runif(data$S, 0.01, 0.1),
    # PP_init = runif(4, 0.8, 1),
    PP_init = runif(data$S, 0.99, 1),
    # K = 300000 *apply( data$C_obs, 1, mean)/apply( data$C_obs, 1, mean)[1],
    K = runif(data$S, Ks/2, Ks*2),
    # logqs = runif(4, -10, -6),
    # iqs = 1/qs,
    # sigma_proc = runif(1, 0.001, 0.1)
    # sigma_proc = sigE
    isigma2 = (1/sigE)^2
    )
    # inits$PE = t(mvrnorm(19, rep(-0.5*inits$sigma_proc^2, data$S), diag(data$S)*inits$sigma_proc^2) )
    inits$PE = matrix(0, nrow=data$S, ncol=Nyear - Year_start)
    # inits$PE = rep(0, Nyear - Year_start)
    # inits$PP_init = inits$PP_init/sum(inits$PP_init)
    inits$iqs = 1/(qs*inits$K/Ks)
    return(inits)
  }
}


init_ll <- lapply(1:chains, function(id) init_create(chain_id = id))

### Check that the inits is OK 
check_inits <- function(x, data, type="Biomass") {
  PP <- matrix(0, data$S, data$Npre)
  Ipred <- matrix(0, data$S, data$Npre)
  LL <- matrix(0, data$S, data$Npre)
  logB <- matrix(0, data$S, data$Npre)
  B <- matrix(0, data$S, data$Npre)
  PP[,1] <- x$PP_init
  logB[,1] <- log(PP[,1]) + log(x$K)
  B[,1] <- exp(logB[,1])
  if (type =="MSY") MSY <- x$MSY
  K <- x$K
  if (type !="MSY") r <- x$r
  nn <- (data$p + 1)
  bb = nn^(nn/(nn-1))
  PE <- x$PE
  for (t in 1:(data$Npre-1))
  {    
    for (i in 1: data$S)
    {
      if (type =="MSY") {
       PP[i,t+1] = max(PP[i,t] + MSY[i]/K[i]*bb*(PP[i,t]- PP[i,t]^nn)-data$C_obs[i,t]/K[i],0.000001)*exp(PE[i,t]);
       logB[i,t+1] = log(PP[i,t+1]) + log(K[i]);
      }
      if (type !="MSY") logB[i,t+1] = log(max(B[i,t] + (r[i]/data$p)*B[i,t]*(1-(B[i,t]/K[i])^data$p) - data$C_obs[i,t], 1)*exp(PE[i,t]));
      B[i,t+1] = exp(logB[i,t+1]);
    }
  }
  
  for (i in 1:data$S)
  {
    for (t in 1:data$Npre) 
    {       
      Ipred[i,t] = log(1/x$iqs[i]) + logB[i,t] -0.5*(exp(x$sigma_proc)*data$ratio)^2
      LL[i,t] = log(dlnorm(data$I_obs[i,t], Ipred[i,t], exp(x$sigma_proc)*data$ratio))
    }
  }
  
  return(list(Ipred, B, LL))
}
# check_inits(init_ll[[1]], data = data, type="Biomass")
# check_inits(init_ll[[2]], data = data, type="MSY")
# check_inits(init_ll[[3]], data = data, type="MSY")

#-----------------------------------------------------------------
iters <- 2000
chains <- 3
burnin <- 0.6 #proportion of chain length used as warmup/burnin

adapt_delta <- 0.99
stepsize <- 0.01
max_treedepth <- 15
stepsize_jitter <- 0

rs
#------------------------------------------------------------------------------
# Model running:
#------------------------------------------------------------------------------
#Original: we have index data all the way through the time series
# and OE and PE are equal

tstart <- Sys.time()
fit_ind_only <- stan(file = paste0("Production_models_STAN/Kotaro/kotaro.stan"), 
            data = data, init = init_ll, #inits, inits),
            iter = iters, chains = chains, cores=chains, seed=123,
            warmup=burnin*iters, verbose=F, thin=1,
            control = list(adapt_delta = adapt_delta, stepsize = stepsize, 
                           max_treedepth = max_treedepth, stepsize_jitter = stepsize_jitter)) #stepsize_jitter default(0), values between 0 and 1
#metric (string, one of "unit_e", "diag_e", "dense_e", defaults to "diag_e")
runtime <- Sys.time() - tstart; runtime
launch_shinystan(fit_ind_only)
#-------------------------------------------------------------------------------
# try with biomass estimates on top of complete index time series: 
fit_bio <- stan(file = paste0("Production_models_STAN/Kotaro/kotaro_estB.stan"), 
            data = data, init = init_ll, #inits, inits),
            iter = iters, chains = chains, cores=chains, seed=123,
            warmup=burnin*iters, verbose=F, thin=1,
            control = list(adapt_delta = adapt_delta, stepsize = stepsize, 
                           max_treedepth = max_treedepth, stepsize_jitter = stepsize_jitter))
# conv and no dts 7/9/24
launch_shinystan(fit_bio)
#-------------------------------------------------------------------------------
# Now the index data doesn't start until year 6 of the catch data...
data$N <- Nyear
data$Npre <- length(cpue_start:Nyear)
data$I_obs <- data$I_obs[,c((cpue_start):Nyear)]

adapt_delta <- 0.999
stepsize <- 0.001

tstart <- Sys.time()
fit_ind_only2 <- stan(file = paste0("Production_models_STAN/Models/ko.stan"), 
                     data = data, init = init_ll, #inits, inits),
                     iter = iters, chains = chains, cores=chains, seed=123,
                     warmup=burnin*iters, verbose=F, thin=1,
                     control = list(adapt_delta = adapt_delta, stepsize = stepsize, 
                                    max_treedepth = max_treedepth, stepsize_jitter = stepsize_jitter)) #stepsize_jitter default(0), values between 0 and 1
#metric (string, one of "unit_e", "diag_e", "dense_e", defaults to "diag_e")
runtime <- Sys.time() - tstart; runtime
launch_shinystan(fit_ind_only2)
# Start to get divergences here: 23, so not crazy. It doesn't like not having an index all the way thru 
# converged but estimates are off too
# halved stepsize from 0.01 to 0.005: still 23
# increased adapt delta from 0.99 to 0.995: 26 dts
# increased adapt delta from 0.99 to 0.999: 30 dts... going up!
# decreased adapt delta to 0.9: 44 dts :(
# adapt del to 0.999 and stepsize to 0.001: still 29 dts! damn! 

# Second try with narrower rs... and no dts!!! ... data quality... 
# cpue dat starting in year 11: 32 dts, need longer chains... 
#-------------------------------------------------------------------------------
# Estimating a singel r?
init_create_1r <- function(chain_id=1){
  set.seed(chain_id+123)
  inits <- list(
    r = runif(1, 0.01, 0.1),
    PP_init = runif(data1$S, 0.99, 1),
    K = runif(data1$S, Ks/2, Ks*2),
    isigma2 = (1/sigE)^2
  )
  inits$PE = matrix(0, nrow=data1$S, ncol=Nyear - 1) #inits$PE = matrix(0, nrow=data$S, ncol=Nyear - Year_start)
  inits$iqs = 1/(qs*inits$K/Ks)
  return(inits)
}
data1 <- data
init_ll_1r <- lapply(1:chains, function(id) init_create_1r(chain_id = id))

fit_ind_only_1r <- stan(file = paste0("Production_models_STAN/Models/ko_1r.stan"), 
                      data = data, init = init_ll_1r, #inits, inits),
                      iter = iters, chains = chains, cores=chains, seed=123,
                      warmup=burnin*iters, verbose=F, thin=1,
                      control = list(adapt_delta = adapt_delta, stepsize = stepsize, 
                                     max_treedepth = max_treedepth, stepsize_jitter = stepsize_jitter))
# no dts with cpue dat starting in year 6!
# cpue dat starting in year 11? 22 divergences/ 26 with clamped down step size and adapt delta
launch_shinystan(fit_ind_only_1r)

#-------------------------------------------------------------------------------
# Now adding biomass index and just estimating a single r for all 3 areas
adapt_delta <- 0.99
stepsize <- 0.01

fit_bio2 <- stan(file = paste0("Production_models_STAN/Models/ko_estB.stan"), 
                data = data, init = init_ll_1r, #inits, inits),
                iter = iters, chains = chains, cores=chains, seed=123,
                warmup=burnin*iters, verbose=F, thin=1,
                control = list(adapt_delta = adapt_delta, stepsize = stepsize, 
                               max_treedepth = max_treedepth, stepsize_jitter = stepsize_jitter))
launch_shinystan(fit_bio2)
# no dts...
# cpue dat starting in year 11? still good! ... bc is using whole cpue time series. duh. 

#------------------------------------------------------------------------------
# Try with entered index and biomass error and 
# get rid of assumption that OE and PE are equal

data2 = list(Npre = length(Year_start:Nyear), 
            N = Nyear,
            S=Narea, 
            p = 0.18815,
            C_obs = t(Catch[Year_start:(Nyear-1),]),
            N_Iobs = N_Iobs,
            I_obs=I_obs, 
            I_cv = I_cv,
            S_Iobs = S_Iobs,
            I_pos = I_pos,
            ratio = 1,
            N_Bobs = N_Bobs,
            B_obs = B_obs,
            B_cv = B_cv,
            S_Bobs = S_Bobs,
            B_pos = B_pos)

fit_bio_ind <- stan(file = paste0("Production_models_STAN/Kotaro/kotaro_estB_estI.stan"), 
                    data = data2, init = init_ll, #inits, inits),
                    iter = iters, chains = chains, cores=chains, seed=123,
                    warmup=burnin*iters, verbose=F, thin=1,
                    control = list(adapt_delta = adapt_delta, stepsize = stepsize, 
                                   max_treedepth = max_treedepth, stepsize_jitter = stepsize_jitter))
# converged and no dts 7/9
launch_shinystan(fit_bio_ind)
# cpue dat starting in year 11? still good!! 
#------------------------------------------------------------------------------
# Entered index and biomass error and OE separate from PE 
# Also, just one r

I_obs_mat <- IAs[cpue_start:Nyear,]
data3 = list(Npre = length(cpue_start:Nyear), 
             N = Nyear, 
             S=Narea, 
             p = 0.18815,
             C_obs = t(Catch[1:(Nyear-1),]),
             #I_obs_mat = I_obs_mat, # matrix of Iobs fo Kotaro's model... 
             N_Iobs = N_Iobs,
             I_obs=I_obs, 
             I_cv = I_cv,
             S_Iobs = S_Iobs,
             I_pos = I_pos,
             ratio = 1,
             N_Bobs = N_Bobs,
             B_obs = B_obs,
             B_cv = B_cv,
             S_Bobs = S_Bobs,
             B_pos = B_pos)

fit_bio_ind2 <- stan(file = paste0("Production_models_STAN/Models/ko_estB_estI.stan"), 
                    data = data3, init = init_ll_1r, #inits, inits),
                    iter = iters, chains = chains, cores=chains, seed=123,
                    warmup=burnin*iters, verbose=F, thin=1,
                    control = list(adapt_delta = adapt_delta, stepsize = stepsize, 
                                   max_treedepth = max_treedepth, stepsize_jitter = stepsize_jitter))
# converged and no dts 7/9
launch_shinystan(fit_bio_ind2)
# cpue dat starting in year 11? still good! 
#-------------------------------------------------------------------------------
# Index only with OE and PE separate:

fit_ind <- stan(file = paste0("Production_models_STAN/Models/ko_estI.stan"), 
                    data = data3, init = init_ll_1r, #inits, inits),
                    iter = iters, chains = chains, cores=chains, seed=123,
                    warmup=burnin*iters, verbose=F, thin=1,
                    control = list(adapt_delta = adapt_delta, stepsize = stepsize, 
                                   max_treedepth = max_treedepth, stepsize_jitter = stepsize_jitter))
#Divergences back with this version. Need to keep PE and OE equal with just the index.
# second try and no dts! data quality... 
launch_shinystan(fit_ind)
# cpue dat starting in year 11? 25dts
#-------------------------------------------------------------------------------

# Diagnostics: 

fit <- fit_bio
list_of_draws <- extract(fit)
print(names(list_of_draws))

fit_summary <- summary(fit)
View(fit_summary$summary)

check_energy(fit)
check_treedepth(fit)
check_divergences(fit)

get_num_divergent(fit)
get_divergent_iterations(fit)
get_num_max_treedepth(fit) 
get_bfmi(fit) 
get_low_bfmi_chains(fit)

library(shinystan)
launch_shinystan(fit_ind_only)


library(bayesplot)
mcmc_scatter(
  as.matrix(fit),
  pars = c("K[1]", "r[1]"),
  np = nuts_params(fit),
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8) )
mcmc_scatter(
  as.matrix(fit),
  pars = c("K[1]", "iqs[1]"),
  np = nuts_params(fit),
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8) )


## So starting from most informative inits then move away

#-------------------------------------------------------------------------------
plot_biomass(fit=fit,years=Nyear,strata=Narea,
             save=FALSE, 
             bio_dat = B_ests, bio_cv = B_cvs,
             sim=TRUE, 
             sim_bio = B, mod_name="Kotaro_Best",
             units = "NA", year1=1)

plot_index1(fit=fit,years=Nyear,strata=Narea,
            save=FALSE, 
            cpue_dat = I_obs_plot, cpue_cv = I_cv_plot,
            sim=TRUE, 
            sim_cpue = sim_IAs, mod_name="Kotaro_Best",
            units = "mt",
            year1 = 1980,
            cpue_name = "Fishery CPUE",
            cpue_type = "fish/effort")

plot_catch(fit=fit,years=Nyear,strata=Narea,
           save=FALSE, 
           dat = data, 
           sim=TRUE, 
           mod_name=mod_name,
           units = "mt",
           year1 = 1980) 

plot_pe(fit = fit, sim_pe = pe, sim=TRUE, save=FALSE, mod_name=mod_name, year1 = 1980)

#-------------------------------------------------------------------------------
# Try with entered index and biomass error: 
data = list(Npre = length(Year_start:Nyear), 
            S=Narea, 
            p = 0.18815,
            C_obs = t(Catch[Year_start:(Nyear-1),]),
            N_Iobs = N_Iobs,
            I_obs=I_obs, 
            I_cv = I_cv,
            S_Iobs = S_Iobs,
            I_pos = I_pos,
            ratio = 1,
            N_Bobs = N_Bobs,
            B_obs = B_obs,
            B_cv = B_cv,
            S_Bobs = S_Bobs,
            B_pos = B_pos)

fit_bio_ind <- stan(file = paste0("Production_models_STAN/Kotaro/kotaro_estB_estI.stan"), 
                data = data, init = init_ll, #inits, inits),
                iter = iters, chains = chains, cores=chains, seed=123,
                warmup=burnin*iters, verbose=F, thin=1,
                control = list(adapt_delta = adapt_delta, stepsize = stepsize, 
                               max_treedepth = max_treedepth, stepsize_jitter = stepsize_jitter))


#-------------------------------------------------------------------------------
### If using TMB
library(TMB)
library(TMBhelper)
version = paste0(getwd(), "/src/model"); dllversion = "model"
compile(paste0(version, ".cpp"))#, "-O1 -g", DLLFLAGS="")
dyn.load(dynlib(version))

set.seed(345)
tmb_data = list(Npre = length(Year_start:Nyear), 
            S=Narea, 
            p = 0.18815,
            C_obs = t(Catch[Year_start:Nyear,]),
            I_obs=t(IA_obs), ratio = 1)

tmb_params <- list(logr = log(runif(tmb_data$S, 0.04, 0.05)),
                   PP_init_logit = boot::logit(rep(0.5, tmb_data$S)),
                   logK = log(runif(tmb_data$S, 100000, 500000)),
                   log_invq = log(runif(tmb_data$S, 1000 , 20000)),
                   log_sigma_proc = log(sigE),
                   PE = matrix(0, nrow=tmb_data$S, ncol=tmb_data$Npre - 1)
) 

# check_inits(x=list(r = exp(tmb_params$logr),
#                  PP_init = boot::inv.logit(tmb_params$PP_init_logit) , 
#                  K = exp(tmb_params$logK)  , 
#                  logqs = (tmb_params$logq) , 
#                  sigma_proc = exp(tmb_params$log_sigma_proc), 
#                  PE = tmb_params$PE), 
#             data = tmb_data, 
#             type = "r")
# 

obj <- TMB::MakeADFun(data = tmb_data, parameters = tmb_params, map = list(logr = rep(factor(1),2)),
                                 random = c("PE"), DLL = dllversion, silent = TRUE)

opt <- fit_tmb(obj, lower=-15, upper=20, getsd=FALSE, bias.correct=FALSE, control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))          

opt$diagnostics
sdreport <- sdreport(obj)
qwe <- summary(sdreport, "report")
asd <- obj$env$report()

rowMeans(tmb_data$I_obs) - rowMeans(asd$logB)

g <- as.numeric(obj$gr(opt$par))
h <- stats::optimHess(opt$par, obj$fn, obj$gr)
par1 <- opt$par - solve(h, g)



