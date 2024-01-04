################################################################################
## DIVERGENT EXAM
##
## Divergent transitions are a huge pain in the ass with the HMC algorythm.
## This is a script to set up a simulated data set identical to spm_model_exploration.R
## script.  After simulating data and setting initial values the script ends
## with a section to loop through various models, priors, etc. to see if we
## can reduce divergent transitions.  
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
  library("wesanderson")}

# if needed : closeAllConnections()
#util <- new.env()
source("code/helper.R")

#setwd("code/")
set.seed(3480)
#--------------------------------------------------------------------------------
# Decide on the Pella tomlinson curve:
# Where does Pmsy occur?
pmsy <- 0.46

pt <- pella_toml(0.46)

pte_p <- pt$pte_p
pte_m <- pt$pte_m

#-------------------------------------------------------------------------------
# simulate some data
# use previous simulation:
oldsim<-readRDS("sim_data/TEST.rds")
#names(oldsim)<-c("years", "strata", "r", "K", "pi", "q1", "q2", "q3",
#                 "q_err", "d_err","b_err", "s_err", "pesd", "rho", "pe",
#                 "bio","bio_est","prop", "cpue", "surv", "C","D", "pte_m","pte_p")
#list2env(oldsim,.GlobalEnv)

# ... or Simulate a new rockfish data set... : ------------------------------------------------

#N years
years <- 45
strata <- 3

r<-0.04
K<-4000
rho <- 0.75

tau_b1 <- 0 #extra variance in early part of time series; 
tau_b2 <- 0 #extra variance in second part of time series
tau_f1 <- 0 #extra variance in early part of time series
tau_f2 <- 0 #extra variance in second part of time series
by_to_C <- c(0.2,1)# how does bycatch compare to known catches? This is a rough approximation that scales the 
# the two catches so that you can have different dynamics for simulation testing of the model
# and the nature of the fishery.  A lot of rock fish are bycatch fisheries and this
# can be explored with simulation testing to see if the SPM can handle that scenario...
# This can also be used to increase or decrease the magnitude of the catch... 
depletion <- 0 #how depleted is the stock relative to virgin biomass

sim <- simulate (years = years, strata=strata, r=r, K=K, rho=rho,
                 tau_b1 = tau_b1, tau_b2 = tau_b2, 
                 tau_f1 = tau_f1, tau_f2 = tau_f2,
                 depletion = depletion,
                 by_to_C = by_to_C) 

list2env(sim,.GlobalEnv)  

plot_sim(sim)

# Save the simulated data for later:
saveRDS(sim,
        "sim_data/TEST.rds")
saveRDS(list(years, strata, r, K, pi, q1, q2, q3, q_err, d_err, b_err, s_err, pesd, rho, pe,
             bio, bio_est, prop, cpue, surv, C, D, pte_m, pte_p),
        "last_sim_depl.rds")

################################################################################
# Format the fake data so that:
# 1) we only have some biomass estimates.  Either randomm removal or just one every 4 years or whatever you want!
# Add in some random missing values: 
bio_est2 <- apply(bio_est,2,function(x) {
  x[sample(c(1:years), floor(years/2))] <- NA
  x
})

#make it ragged so there are different number of samples per strata
bio_est2[runif(3,1,30),1]<-NA
bio_est2[runif(3,1,30),4]<-NA

# Or, lets mimic SEO YE and have no biomass estimates at the beginning of the time series until the pop is 
# already heading down and then have the biomass estimate staggered every four years
# in each of the strata with only one strata being sampled per year: 
{
  bio_est2 <- data.frame()
  for (i in 1:years){
    bio_est2[i,c(1:strata)] <- NA
  }
  for (i in 1:strata){
    bio_est2[seq(10+i,years,strata),i]<-bio_est[seq(10+i,years,strata),i]
  }
  
  B_cv2 <- ifelse(is.na(bio_est2),NA,min(b_err)) 
  
  
  N_Bobs <- sum (!is.na(bio_est2))
  
  #N_Bobs <- c(length(bio_est2[,1][!is.na(bio_est2[,1])]),
  #            length(bio_est2[,2][!is.na(bio_est2[,2])]),
  #            length(bio_est2[,3][!is.na(bio_est2[,3])]),
  #            length(bio_est2[,4][!is.na(bio_est2[,4])]))
  
  B_pos <- list()
  for (i in 1:strata) {
    B_pos[[i]] <- which(!is.na(bio_est2[,i]))
  }
  #B_pos <- list(which(!is.na(bio_est2[,1])),
  #              which(!is.na(bio_est2[,2])),
  #              which(!is.na(bio_est2[,3])),
  #              which(!is.na(bio_est2[,4])))
  
  #I guess we vectorize this to deal with the ragged arrays...
  B_obs <- as.vector(bio_est2[!is.na(bio_est2)])
  B_cv <- as.vector(B_cv2[!is.na(B_cv2)])
  N_Bobs <- N_Bobs#sum(N_Bobs) #total number of observations
  S_Bobs <- vector()
  for (i in 1:strata) {
    S_Bobs[i] <- length(bio_est2[,i][!is.na(bio_est2[,i])])
  }
  #S_Bobs <- c(length(bio_est2[,1][!is.na(bio_est2[,1])]),
  #            length(bio_est2[,2][!is.na(bio_est2[,2])]),
  #            length(bio_est2[,3][!is.na(bio_est2[,3])]),
  #            length(bio_est2[,4][!is.na(bio_est2[,4])])) # number of observations for each strata
  B_pos <- as.vector(unlist(B_pos))
  
  # 2) Extra variance to consider?
  bio_est2; B_pos
  # declare which years to apply extra variance to.  Here we apply it to all years before year 31
  bio_tau_vec<-ifelse(B_pos < 31,0,0)  #same length as B_obs and set up as time blocks for extra variance term; see S_Bobs... 
  #if you want different tau's by strata you'll need to coordinate that in this vector
  # only interested in one extra variance term across strata so that is how this is set up now
  bio_tau_no<-length(unique(bio_tau_vec[bio_tau_vec != 0]))
} #number of extra variance terms to estimate... 

# 3) CPUE index 1:
# either remove random values, or specify what values you want to get rid of
# There is a section to declare extra variance if you want to consider it... 
{
  cpue2 <- apply(cpue,2,function(x) {
    x[sample(c(1:years), floor(years/3))] <- NA
    x
  })
  
  cpue2[runif(3,1,30),1]<-NA
  cpue2[runif(3,1,30),strata]<-NA
  
  # or mimick the iphc cpue and lop off the last and the first ~12 years of samples
  cpue2<-cpue
  cpue2[c(1:14,years),]<-NA
  
  I1_cv2 <- ifelse(is.na(cpue2),NA,min(q_err))
  
  N_I1obs <- sum (!is.na(cpue2))
  #N_I1obs <- c(length(cpue2[,1][!is.na(cpue2[,1])]),
  #             length(cpue2[,2][!is.na(cpue2[,2])]),
  #             length(cpue2[,3][!is.na(cpue2[,3])]),
  #             length(cpue2[,4][!is.na(cpue2[,4])]))
  
  I1_pos <- list()
  for (i in 1:strata) {
    I1_pos[[i]] <- which(!is.na(cpue2[,i]))
  }
  #I1_pos <- list(which(!is.na(cpue2[,1])),
  #               which(!is.na(cpue2[,2])),
  #               which(!is.na(cpue2[,3])),
  #               which(!is.na(cpue2[,4])))
  
  #I guess we vectorize this to deal with the ragged arrays...
  I1_obs <- as.vector(cpue2[!is.na(cpue2)])
  I1_cv <- as.vector(I1_cv2[!is.na(I1_cv2)])
  N_I1obs <- N_I1obs #sum(N_I1obs) #total number of observations
  S_I1obs <- vector()
  for (i in 1:strata) {
    S_I1obs[i] <- length(cpue2[,i][!is.na(cpue2[,i])])
  }
  #S_I1obs <- c(length(cpue2[,1][!is.na(cpue2[,1])]),
  #             length(cpue2[,2][!is.na(cpue2[,2])]),
  #             length(cpue2[,3][!is.na(cpue2[,3])]),
  #             length(cpue2[,4][!is.na(cpue2[,4])])) # number of observations for each strata
  I1_pos <- as.vector(unlist(I1_pos))
  
  #Extra variance to consider?
  cpue2; I1_pos
  i1_tau_vec<-rep(0,length(I1_pos))  #same length as I1_obs and set up as time blocks for extra variance term; see S_Bobs... 
  #if you want different tau's by strata you'll need to coordinate that in this vector
  # only interested in one extra variance term across strata so that is how this is set up now
  i1_tau_no<-length(unique(i1_tau_vec[i1_tau_vec != 0])) #number of extra variance terms to estimate... 
}

# 4) Survey index 2:
# either remove random values, or specify what values you want to get rid of
# There is a section to declare extra variance if you want to consider it...
{
  surv2 <- apply(surv,2,function(x) {
    x[sample(c(1:years), floor(years/3))] <- NA
    x
  })
  
  surv2[runif(3,1,30),1]<-NA
  surv2[runif(3,1,30),strata]<-NA
  
  # or mimick the iphc cpue and lop off the last and the first ~12 years of samples
  surv2<-surv
  surv2[c(1:14,years),]<-NA
  
  I2_cv2 <- ifelse(is.na(surv2),NA,min(s_err))
  
  N_I2obs <- sum (!is.na(surv2))
  #N_I2obs <- c(length(surv2[,1][!is.na(surv2[,1])]),
  #             length(surv2[,2][!is.na(surv2[,2])]),
  #             length(surv2[,3][!is.na(surv2[,3])]),
  #             length(surv2[,4][!is.na(surv2[,4])]))
  
  I2_pos <- list()
  for (i in 1:strata) {
    I2_pos[[i]] <- which(!is.na(surv2[,i]))
  }
  #I2_pos <- list(which(!is.na(surv2[,1])),
  #               which(!is.na(surv2[,2])),
  #               which(!is.na(surv2[,3])),
  #               which(!is.na(surv2[,4])))
  
  #I guess we vectorize this to deal with the ragged arrays...
  I2_obs <- as.vector(surv2[!is.na(surv2)])
  I2_cv <- as.vector(I2_cv2[!is.na(I2_cv2)])
  N_I2obs <- N_I2obs #sum(N_I2obs) #total number of observations
  S_I2obs <- vector()
  for (i in 1:strata) {
    S_I2obs[i] <- length(surv2[,i][!is.na(surv2[,i])])
  }
  #S_I2obs <- c(length(surv2[,1][!is.na(surv2[,1])]),
  #             length(surv2[,2][!is.na(surv2[,2])]),
  #             length(surv2[,3][!is.na(surv2[,3])]),
  #             length(surv2[,4][!is.na(surv2[,4])])) # number of observations for each strata
  I2_pos <- as.vector(unlist(I2_pos))
  
  #Extra variance to consider?
  surv2; I2_pos
  i2_tau_vec<-rep(0,length(I2_pos))  #same length as B_obs and set up as time blocks for extra variance term; see S_Bobs... 
  #if you want different tau's by strata you'll need to coordinate that in this vector
  # only interested in one extra variance term across strata so that is how this is set up now
  i2_tau_no<-length(unique(i2_tau_vec[i2_tau_vec != 0])) #number of extra variance terms to estimate...
}

# 5) Take a LOOK at your fake data set
{
  par(mfrow=c(3,2))
  plot(C[,1],ylim=c(0,max(C)), type="l", cex=1, col = pal[1])
  for (s in 2:strata){ #s<-4
    lines(C[,s],type="l", cex=0.5, col=pal[s])
  }
  
  plot(D[,1],ylim=c(0,max(C)), type="l", cex=1, col = pal[1])
  for (s in 2:strata){ #s<-4
    lines(D[,s],type="l", cex=0.5, col=pal[s])
  }
  
  plot(cpue2[,1],ylim=c(0,max(cpue2, na.rm=T)), type="b", cex=0.5, col = pal[1])
  for (s in 2:strata){
    lines(cpue2[,s],type="b", cex=0.5, col=pal[s])
  }
  
  plot(surv2[,1],ylim=c(0,max(surv2, na.rm=T)), type="b", cex=0.5, col = pal[1])
  for (s in 2:strata){
    lines(surv2[,s],type="b", cex=0.5, col=pal[s])
  }
  
  plot(bio[,1],ylim=c(0,max(bio)*1.2), type="l", cex=1, col = pal[1])
  points(bio_est2[,1], type="p", cex=1, col = pal[1], pch=18)
  for (s in 2:strata){
    lines(bio[,s],type="l", cex=0.5, col=pal[s])
    points(bio_est2[,s],type="p", cex=1, col=pal[s], pch=18)
  }
  
  plot(pe); abline(h=0, col="red")
}

#-------------------------------------------------------------------------------
# SWITCHES
# What indices do we want to fit to? Create dummy data if we aren't fitting to it

# do we have biomass data?
bio_switch <- 1 #0 means no biomass data, 1 means there is biomass data
if (bio_switch == 0){
  #fake data to be ignored
  B_obs <- rep(1,strata)
  B_cv <- rep(1,strata)
  N_Bobs <- strata
  B_pos <- rep(1,strata)
  S_Bobs <- rep(1,strata)
  bio_tau_vec<-rep(0,length(B_pos))  
  bio_tau_no<-length(unique(bio_tau_vec[bio_tau_vec != 0]))
}

# do we have an index 1?
ind1_switch <- 1 #0 means no index data, 1 means there is index data
if (ind1_switch == 0){
  #fake data to be ignored
  I1_obs <- rep(1,strata)
  I1_cv <- rep(1,strata)
  N_I1obs <- strata
  I1_pos <- rep(1,strata)
  S_I1obs <- rep(1,strata)
  i1_tau_vec<-rep(0,length(I1_pos))  
  i1_tau_no<-length(unique(i1_tau_vec[i1_tau_vec != 0]))
}

# do we have an index 2?
ind2_switch <- 0 #0 means no second index data, 1 means there is secondary index data
if (ind2_switch == 0){
  #fake data to be ignored
  I2_obs <- rep(1,strata)
  I2_cv <- rep(1,strata)
  N_I2obs <- strata
  I2_pos <- rep(1,strata)
  S_I2obs <- rep(1,strata)
  i2_tau_vec<-rep(0,length(I2_pos))  
  i2_tau_no<-length(unique(i2_tau_vec[i2_tau_vec != 0]))
}

# which parameterization do we want to use
paramaterization = "NC" #NC #C = centered, NC = non-centered
if (paramaterization == "C") {
  pte <- pte_p
  msy_frac <- get_pmsy_p(pte_p)
} else if (paramaterization == "NC") {
  pte <- pte_m
  msy_frac <- get_pmsy_m(pte_m)
}

#-------------------------------------------------------------------------------
# Bundle the data! OK if you have extra data relative to the model you want to run
# stan will ignore the extraneous stuff... 

# FLAG you can change PRIOR values in this list!!
# FLAG you can decide if yu want to estimate PROCESS ERROR in this list! 
C_obs <- C+0.000001 #can't have catch of 0...
C_cv <- cbind(rep(0.05,years),rep(0.05,years),
              rep(0.05,years),rep(0.05,years)); C_cv <- C_cv[,c(1:strata)]

D_obs <- D+0.000001 #can't have catch of 0...
D_cv <- cbind(rep(0.4,years),rep(0.4,years),
              rep(0.4,years),rep(0.4,years)); D_cv <- D_cv[,c(1:strata)]

fish_dat2 <- list(N = years,
                  S = strata,
                  #catch data
                  C_obs = C_obs,  
                  C_cv = C_cv,
                  #discard data
                  D_obs = D_obs,  #can't have catch of 0...
                  D_cv = D_cv,
                  #observed biomass
                  B_obs = B_obs, 
                  B_cv = B_cv,
                  N_Bobs = N_Bobs,
                  B_pos = B_pos,
                  S_Bobs = S_Bobs, 
                  bio_tau_vec = bio_tau_vec,
                  bio_tau_no = bio_tau_no,
                  # Index 1
                  I1_obs = I1_obs,
                  I1_cv = I1_cv, 
                  N_I1obs = N_I1obs,
                  I1_pos = I1_pos,
                  S_I1obs = S_I1obs,
                  i1_tau_vec = i1_tau_vec,
                  i1_tau_no = i1_tau_no,
                  #index 2
                  I2_obs = I2_obs,
                  I2_cv = I2_cv, 
                  N_I2obs = N_I2obs,
                  I2_pos = I2_pos,
                  S_I2obs = S_I2obs,
                  i2_tau_vec = i2_tau_vec,
                  i2_tau_no = i2_tau_no,
                  #pella tomlinson settings
                  pte = pte,  #pella tomlinson exponent 0.18815 to set Bmsy at B40; 0.626099 to set Bmsy at B46
                  msy_frac = msy_frac,
                  #priors
                  rp1 = 1, #1st parameter for r prior
                  rp2 = 1,# 2nd parameter for r prior
                  k1 = 10, #mean for lognormal K parameter
                  k2 = 2, #sd for lognormal K parameter
                  alpha = rep(1,strata), #c(1,1,1,1),
                  pe_bound = -4,
                  P1_mu = log(1),  #make a switch for this... 
                  P1_sig = 0.1,
                  isig2_mu = 5,  #prior mu for fit of P
                  isig2_sig = 0.102,  #prior sig for fit of P
                  #model switches
                  r_prior_switch = 0, #0=beta, 1=gamma, 2=lognormal
                  pe_switch = 0, #0 = no process error, 1 means random pe, 2 means autocorrelated pe (in development)
                  bio_switch = bio_switch,
                  ind1_switch = ind1_switch,
                  ind2_switch = ind2_switch
)     #0 means no process error, 1 means estimate pe

#-------------------------------------------------------------------------------
# Establish initial values: 
#one chain for starters... 
true_inits <- list(r=r*1.1, K=K*0.9, 
                   iq1=c(1/q1[1],1/q1[2],1/q1[3],1/q1[4]), 
                   iq2=c(1/q2[1],1/q2[2],1/q2[3],1/q2[4]), 
                   isigma2=50, logvar=-5, eps = rep(0,years),
                   P=P1, C = C+0.000001, pi = c(0.25,0.25,0.25,0.25))

#4 chains for awesome
P1guess <- exp(fish_dat2$P1_mu)
Pinit <-matrix(nrow=years,ncol=strata)
Pinit[1,] <- rep(P1guess,strata)
if (bio_switch == 1) {
  for (i in 1:strata) { #i <- 2
    j <- 1
    bioests <- bio_est2[,i][!is.na(bio_est2[,i])]
    samps <- length(bioests)
    pests <- bioests / (max(bioests) * 1/P1guess) # - (1 - P1guess)
    pests[c((length(pests) + 1):(length(pests) + 5))] <- pests[length(pests)]
    for (t in 2:years) { #t <- 3
      if (is.na(bio_est2[t,i])) {
        if (pests[j] == P1guess) { #if you haven't gotten to the first biomass estimate... 
          Pinit[t,i] <- Pinit[t-1,i] - (pests[j] - pests[j+1]) * 1 / samps
        } else { #have passed your first biomass est
          Pinit[t,i] <- Pinit[t-1,i] - (pests[j] - pests[j+1]) * 1 / (length(bioests))
        }
        Pinit[t,i] 
      } else {
        #Pinit[t,i] <- P1guess * bio_est2[t,i] / max(bio_est2[,i], na.rm=T)
        Pinit[t,i] <- pests[j]
        j <- j + 1
      }
    }
  }
} else {
  for (s in 1:strata) {
    low <- min(cpue2[,s], na.rm=T) / max(cpue2[,s] - P1guess, na.rm=T)
    by <- (P1guess - low) / years
    Pinit[,s] <- seq(P1guess,
                     low,
                     -by)
  }
}

{plot(Pinit[,1], col="red", ylim=c(0,1)); points(Pinit[,2], col="blue"); points(Pinit[,3], col="forestgreen")
  lines(prop[,1], col="red"); lines(prop[,2], col="blue");lines(prop[,3], col="forestgreen")
  points(bio_est2[,1]/max(bio_est2[,1], na.rm=T)*P1guess, col="red", pch=18) 
  points(bio_est2[,2]/max(bio_est2[,2], na.rm=T)*P1guess, col="blue", pch=18) 
  points(bio_est2[,3]/max(bio_est2[,3], na.rm=T)*P1guess, col="forestgreen", pch=18) }

P1<-Pinit
P2<-Pinit * cbind(rnorm(years,1,0.1),rnorm(years,1,0.1),rnorm(years,1,0.1))
P3<-Pinit * cbind(rnorm(years,1,0.1),rnorm(years,1,0.1),rnorm(years,1,0.1))
P4<-Pinit * cbind(rnorm(years,1,0.1),rnorm(years,1,0.1),rnorm(years,1,0.1))

muq <- fish_dat2$P1_mu
mu1<-cbind(rep(muq,years),rep(muq,years),rep(muq,years),rep(muq,years)); mu1 <- mu1[,c(1:strata)]
mu2<-cbind(rep(muq,years),rep(muq,years),rep(muq,years),rep(muq,years))+rnorm(years*4,0,0.05); mu2 <- mu2[,c(1:strata)] 
mu3<-cbind(rep(muq,years),rep(muq,years),rep(muq,years),rep(muq,years))+rnorm(years*4,0,0.05); mu3 <- mu3[,c(1:strata)]
mu4<-cbind(rep(muq,years),rep(muq,years),rep(muq,years),rep(muq,years))+rnorm(years*4,0,0.05); mu4 <- mu4[,c(1:strata)]

#mu1<-log(P1)
#mu2<-log(P1)
#mu3<-log(P1)
#mu4<-log(P1)
maxbio <- 0; pis <- vector()
for (s in 1:strata) {
  maxbio <- maxbio + max(bio_est2[,s], na.rm=T)
}
for (s in 1:strata) {
  if(s == strata) {
    pis[s] <- 1 - sum(pis)
  } else {
    pis[s] <- max(bio_est2[,s], na.rm=T)/maxbio
  }
}

pi1 <- pis; sum(pi1)
pi2 <- c(pis[c(1:(strata-1))] * runif(strata-1, 1, 1.1)); pi2 <- c(pi2, 1-sum(pi2)); sum(pi2)
pi3 <- c(pis[c(1:(strata-1))] * runif(strata-1, 1, 1.1)); pi3 <- c(pi3, 1-sum(pi3)); sum(pi3)
pi4 <- c(pis[c(1:(strata-1))] * runif(strata-1, 1, 1.1)); pi4 <- c(pi4, 1-sum(pi4)); sum(pi4)
# or better initial P values from last run... jitter them a little... 

bioguess <- c() #guess a rough estimate of biomass to get started.. 
iq1 <- vector()
iq2 <- vector()
if (bio_switch == 1) { 
  for (s in 1:strata) {
    iq1[s] <- 1 / (mean(cpue2[,s], na.rm=T) / mean(bio_est2[,s], na.rm=T))
    iq2[s] <- 1 / (mean(surv2[,s], na.rm=T) / mean(bio_est2[,s], na.rm=T))
  }
} else { 
  for (s in 1:strata) {
    iq1[s] <- 1 / (mean(cpue2[,s], na.rm=T) / mean(bioguess[,s], na.rm=T))
    iq2[s] <- 1 / (mean(surv2[,s], na.rm=T) / mean(bioguess[,s], na.rm=T))
  }
  
}

inits1 <- list(r=0.05, K=5000, 
               iq1=iq1,#c(50,50,50,50), 
               iq2=iq2,#c(100,600,1000,100),
               i1tau2 = rep(20,strata), #c(20,40,40,40),
               i2tau2 = rep(20,strata), #c(40,20,20,20),
               isigma2=100,logvar=-6, eps = rep(0,years),  
               P=P1, mu = mu1,
               medP1 = c(0.75,0.65,0.7),
               C = C_obs, D = D_obs, 
               pi = pi1) #c(0.25,0.25,0.25,0.25))
inits2 <- list(r=0.03, K=3000, 
               iq1=iq1 * runif(strata,1,1.1),#c(50,50,50,50), 
               iq2=iq2 * runif(strata,1,1.1),
               i1tau2 = c(30,50,30,15),
               i2tau2 = c(2,2,2,2),
               isigma2=140, logvar=-5, eps = rep(0,years),
               P=P2, mu = mu2,
               medP1 = c(0.8,0.6,0.77),
               C = C_obs, D = D_obs,  pi = pi2)
inits3 <- list(r=0.02, K=4500, 
               iq1=iq1 * runif(strata,1,1.1),#c(50,50,50,50), 
               iq2=iq2 * runif(strata,1,1.1),
               i1tau2 = rep(30,strata),
               i2tau2 = rep(25,strata),
               isigma2=80, logvar=-7, eps = rep(0,years),
               P=P3, mu = mu3,
               medP1 = c(0.65,0.8,0.72),
               C = C_obs, D = D_obs, pi = pi3)
inits4 <- list(r=0.06, K=2000, 
               iq1=iq1 * runif(strata,1,1.1),#c(50,50,50,50), 
               iq2=iq2 * runif(strata,1,1.1),
               i1tau2 = rep(30,strata),
               i2tau2 = rep(30,strata),
               isigma2=100, logvar=-8, eps = rep(0,years),
               P=P4, mu = mu4,
               medP1 = c(0.75,0.65,0.7),
               C = C_obs, D = D_obs,  pi = pi4)

################################################################################
# OK, lets run some short chain models over various parameterization and model structure
# and se if we can't do something about those divergent transisitions... 

mod_list <- c('PT_NC_margQ_xV_D.stan',
              'PT_NC_xV_D.stan')

#metrics
metrics<-c("diag_e","unit_e","dense_e")

#isigma2
sig_mus <- c(3,5,7)
sig_sigs <- c(0.1,0.05,0.01)

#r
r_opts <- 4
r_dists <- c(0,0,2,2)
rp1s <- c(1,1.483,log(0.04),log(0.03))
rp2s <- c(1,22.908,0.1,0.1)

#K none for now

# hmc settings:
adapt_delta = 0.999
stepsize = 0.001
max_treedepth = 15
stepsize_jitter = 0


results <- data.frame()
i <- 1
for (m in mod_list) {   #m <- mod_list[1]
  for (r in 1:r_opts) {   # r <- 1
    for (s in 1:length(sig_mus)) { # s <- 1
      for (ss in 1:length(sig_sigs)) { #  ss <- 1
        for (mm in metrics) { #mm <-metrics[1]
          fish_dat2$isig2_mu <- sig_mus[s]
          fish_dat2$isig2_sig <- sig_sigs[ss]
          fish_dat2$r_prior_switch <- r_dists[r]
          fish_dat2$rp1 <- rp1s[r]
          fish_dat2$rp2 <- rp2s[r]
          
          results[i,"model"] <- m
          results[i,"iters"] <- iters <- 10000
          results[i,"metric"] <- mm
          results[i,"adapt_delta"] <- adapt_delta
          results[i,"stepsize"] <- stepsize
          results[i,"max_treedepth"] <- max_treedepth
          results[i,"stepsize_jitter"] <- stepsize_jitter
          results[i,"isig2mu"] <- sig_mus[s]
          results[i,"isig2_sig"] <- sig_sigs[ss]
          results[i,"r_dist"] <-r_dists[r]
          results[i,"rp1"] <- rp1s[r]
          results[i,"rp2"] <- rp2s[r]
          
          tstart <- Sys.time()
          fit <- stan(file = m, data = fish_dat2, 
                      iter = iters, chains = 1, cores=1,
                      init=list(inits1), warmup=0.8 * iters, verbose=F, thin=20,
                      control = list(adapt_delta = adapt_delta, stepsize = stepsize, 
                                     max_treedepth = max_treedepth, 
                                     stepsize_jitter = stepsize_jitter,
                                     metric = mm))
          runtime <- Sys.time() - tstart
          
          results[i,"energy"] <- get_bfmi(fit)
          results[i,"treedepth"] <- get_num_max_treedepth(fit)
          results[i,"prop_divergent"] <- get_num_divergent(fit)/ ((0.2 * iters ) / 20)
          results[i,"run_time"] <- runtime
          
          i <- i+1
          
          print(results)
          write.csv(results, paste0("divergency_exam.csv"))
        }
      }
    }
  }
}

#
gs <-rgamma(10000,0.4, 2); hist(gs[gs < 1.001], breaks = 100)


#priors
rp1 = 1, #1st parameter for r prior
rp2 = 1,# 2nd parameter for r prior
k1 = 10, #mean for lognormal K parameter
k2 = 2, #sd for lognormal K parameter
alpha = rep(1,strata), #c(1,1,1,1),
pe_bound = -4,
P1_mu = log(1),  #make a switch for this... 
P1_sig = 0.1,
isig2_mu = 5,  #prior mu for fit of P
isig2_sig = 0.102,  #prior sig for fit of P
#model switches
r_prior_switch = 0, #0=beta, 1=gamma, 2=lognormal
pe_switch = 0, #0 = no process error, 1 means random pe, 2 means autocorrelated pe (in development)
bio_switch = bio_switch,
ind1_switch = ind1_switch,
ind2_switch = ind2_switch

gamma(1.709, 0.00861)


