rm(list=ls())
gc()

library(MASS) ## needed for the multivariate normal distribution
library(dplyr)

Niter = 1000   # number of simulation iterations
Narea = 3      # number of sub-area


Catch_type = "F"    # "MSY", or "F"

#for (i in 1:Niter){
i = 123

### Setting up the OM 

set.seed(i)
rs = runif(Narea, 0.045, 0.055)    # variability in rs, but will estimate on r in model
Ks = runif(Narea, 100000, 200000)
Nyear = 50
Year_start = 11   # when we have survey/CPUE info from 
p = 0.18815      # Pella Tom productivity curve exponent
var_ratio = 1  # ration proc / obs error variance; PJ: Swith this off and have separate PE and OE?
sigE = 0.05
varE = rep(sigE^2, Narea)
varO = var_ratio*varE ## OE will be estimated... 
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
#    if (Catch_type == "MSY") val = sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-10), seq(Hmax, 0.5, length.out=10))*MSYs[x], 0.05*MSYs[x]))
    
  ## Catch is dependent of stock size: HMSy * scaler describing fishing pressure over time... 0.5HMSY to 3*HMSYs then down again.
    # Kotaro's eg
#    if (Catch_type == "F")
C_br1 <- 25
C_br2 <- 40
# Alt 1 Kotaro's; increasing fishing pressure than gradual decline in fishing pressure 
H1 = sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-20), seq(Hmax, 0.5, length.out=20))*HMSYs[x], 0.1*HMSYs[x]))
matplot(1:Nyear, H1, type="l", ylim = c(0,max(H1)))
    # Alt 2: increasing fishing pressure, then rapid decline and low levels
H2 = sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-C_br1), seq(Hmax,0.1*Hmax,length.out=(C_br2-C_br1)), seq(Hmax*0.1, 0.25, length.out=Nyear-C_br2))*HMSYs[x], 0.1*HMSYs[x]))
matplot(1:Nyear, H2, type="l", ylim = c(0,max(H2)))
    # Alt 3
H3 = sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, 0.75* Hmax, length.out=Nyear-C_br1), seq(0.5* Hmax,0.1*Hmax,length.out=(C_br2-C_br1)), seq(Hmax*0.1, 0.25, length.out=Nyear-C_br2))*HMSYs[x], 0.1*HMSYs[x]))
#H3 = sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-15), seq(Hmax*0.5, 0.5, length.out=15))*HMSYs[x], 0.1*HMSYs[x]))
matplot(1:Nyear, H3, type="l", ylim = c(0,max(H3)))
    # Alt 4: increasing fishing pressure, then rapid decline and minimal fishing
H4 = sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-C_br1), seq(1,0.01*Hmax,length.out=(C_br2-C_br1)), seq(Hmax*0.01, 0.25, length.out=Nyear-C_br2))*HMSYs[x], 0.1*HMSYs[x]))
matplot(1:Nyear, H4, type="l", ylim = c(0,max(H4)))
# Alt 5: increasing fishing pressure, then slightly tappering, then rapid decline
H5 = sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-C_br1), seq(Hmax,0.75*Hmax,length.out=(C_br2-C_br1)), seq(Hmax*0.75, 0.5*Hmax, length.out=Nyear-C_br2))*HMSYs[x], 0.1*HMSYs[x]))
matplot(1:Nyear, H5, type="l", ylim = c(0,max(H5)))
# Alt 6.. gradual increase in fishing pressure over entire time
H6 = sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.1, 0.5*Hmax, length.out=Nyear-C_br1), seq(0.5*Hmax,0.75*Hmax,length.out=(C_br2-C_br1)), seq(Hmax*0.75, 0.7*Hmax, length.out=Nyear-C_br2))*HMSYs[x], 0.1*HMSYs[x]))
matplot(1:Nyear, H6, type="l", ylim = c(0,max(H6)))

Hlist <- list(H1,H2,H3,H4,H5, H6)

Hmix <- data.frame(matrix(nrow=Nyear))
for (i in 1:Narea){
  Hx <- Hlist[[sample(length(Hlist), 1)]]
  Hmix[,i] <- Hx[,sample(1:Narea,1)]
  if(i > 1){
    if (identical(Hmix[,i],Hmix[,i-1])) {
      Hx <- Hlist[[sample(length(Hlist), 1)]]
      Hmix[,i] <- Hx[,sample(1:Narea,1)]
    }
  }
}
matplot(1:Nyear, Hmix, type="l", ylim = c(0,max(Hmix)))

par(mfrow=c(3,2), mar=c(2,4,1,1))
matplot(1:Nyear, H6, type="l", ylim = c(0,0.1))
matplot(1:Nyear, H5, type="l", ylim = c(0,0.1))
matplot(1:Nyear, H1, type="l", ylim = c(0,0.1))
matplot(1:Nyear, H2, type="l", ylim = c(0,0.1))
matplot(1:Nyear, H3, type="l", ylim = c(0,0.1))
matplot(1:Nyear, H4, type="l", ylim = c(0,0.1))

#------------------------------------------------------------------------------  
# Lets generate a gazzilion simulations to get a sense of how these are doing
sim_dat <- list()
sim_stats <-data.frame()
nsims <- 50
#h_hist <- "opt6" #opt2, opt3, opt4, opt5, opt6, opt_mix
#harv_opts <-c("opt1","opt2","opt3","opt4","opt5","opt6") = c("IncMax_Taper","IncMax_Taper_plat","IncMod_Drop_taper",
#                                                             "IncMax_shutdown","IncMax_modTaper","IncMod_notaper")
harv_opts <- c("Harv1","Harv2","Harv3","Harv4","Harv5",'Harv6')

fpres <- c(3,4,5,6)

iter <- 1

for (opts in harv_opts) {
  h_hist <- opts
  for (fs in fpres){
    Hmax <- fs
    for (i in 1:nsims){ #i <- 1
      if (h_hist == "Harv3") { #Kotaro; increasing and then declining fishing pressure, Opt 1 above
        H = sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-20), seq(Hmax, 0.5, length.out=20))*HMSYs[x], 0.1*HMSYs[x]))
      }
      
      if (h_hist == "Harv4") { #increasing, then quicker decline, then low level: opt 2 above
        H = sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-C_br1), seq(Hmax,0.1*Hmax,length.out=(C_br2-C_br1)), seq(Hmax*0.1, 0.25, length.out=Nyear-C_br2))*HMSYs[x], 0.1*HMSYs[x]))
      }
      
      if (h_hist == "Harv6") { #increasing gradual, quick decline then declining to end; opt3 above
        H = sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, 0.75* Hmax, length.out=Nyear-C_br1), seq(0.5* Hmax,0.1*Hmax,length.out=(C_br2-C_br1)), seq(Hmax*0.1, 0.25, length.out=Nyear-C_br2))*HMSYs[x], 0.1*HMSYs[x]))
      }
      
      if (h_hist == "Harv5") { #increasing, then near total shut down and low level fishing; opt 4 above
        H = sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-C_br1), seq(1,0.01*Hmax,length.out=(C_br2-C_br1)), seq(Hmax*0.01, 0.25, length.out=Nyear-C_br2))*HMSYs[x], 0.1*HMSYs[x]))
      }
      
      if (h_hist == "Harv1") { # Increasing then slow tapr; opt 5 above
        H = sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-C_br1), seq(Hmax,0.75*Hmax,length.out=(C_br2-C_br1)), seq(Hmax*0.75, 0.5*Hmax, length.out=Nyear-C_br2))*HMSYs[x], 0.1*HMSYs[x]))
      }
      
      if (h_hist == "Harv2") { #gradual increase through all with slight tapering at end; opt 6 above
        H = sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.1, 0.5*Hmax, length.out=Nyear-C_br1), seq(0.5*Hmax,0.75*Hmax,length.out=(C_br2-C_br1)), seq(Hmax*0.75, 0.7*Hmax, length.out=Nyear-C_br2))*HMSYs[x], 0.1*HMSYs[x]))
      }
      
      if (h_hist == "Harv_mix") {
        H1 <- sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-20), seq(Hmax, 0.5, length.out=20))*HMSYs[x], 0.1*HMSYs[x]))
        H2 <- sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-C_br1), seq(Hmax,0.1*Hmax,length.out=(C_br2-C_br1)), seq(Hmax*0.1, 0.5, length.out=Nyear-C_br2))*HMSYs[x], 0.1*HMSYs[x]))
        H3 <- sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-15), seq(Hmax*0.5, 0.5, length.out=15))*HMSYs[x], 0.1*HMSYs[x]))
        H4 <- sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-C_br1), seq(1,0.1*Hmax,length.out=(C_br2-C_br1)), seq(Hmax*0.5, 0.5, length.out=Nyear-C_br2))*HMSYs[x], 0.1*HMSYs[x]))
        H5 <- sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-C_br1), seq(Hmax,0.75*Hmax,length.out=(C_br2-C_br1)), seq(Hmax*0.75, Hmax*0.5, length.out=Nyear-C_br2))*HMSYs[x], 0.1*HMSYs[x]))
        H6 <- sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.1, 0.5*Hmax, length.out=Nyear-C_br1), seq(0.5*Hmax,0.75*Hmax,length.out=(C_br2-C_br1)), seq(Hmax*0.75, 0.7*Hmax, length.out=Nyear-C_br2))*HMSYs[x], 0.1*HMSYs[x]))
        Hlist <- list(H1,H2,H3,H4,H5,H6)
        
        Hmix <- data.frame(matrix(nrow=Nyear))
        for (na in 1:Narea){
          Hx <- Hlist[[sample(length(Hlist), 1)]]
          H[,na] <- Hx[,sample(1:Narea,1)]
          if(na > 1){
            if (identical(Hmix[,na],Hmix[,na-1])) {
              Hx <- Hlist[[sample(length(Hlist), 1)]]
              Hmix[,na] <- Hx[,sample(1:Narea,1)]
            }
          }
        }
      }
      
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
      
      # Record True contrasts:
      #total contrast
      sim_stats[iter,"iter"] <- i
      sim_stats[iter,"harv_hist"] <- h_hist
      sim_stats[iter,"Hmax"] <- Hmax
      for(j in 1:Narea){ #j <-1
        #total contrast
        sim_stats[iter,paste0("true_tot_cont",j)] <- (max(B[,j], na.rm=T) - min(B[,j],na.rm=T))/min(B[,j],na.rm=T)
        #rebound of the population from the low point
        sim_stats[iter,paste0("true_reb_cont",j)] <- (tail(na.omit(B[,j]), 1) - min(B[,j], na.rm=T))/min(B[,j], na.rm=T)
      } 
      
      ## abundance indices
      qs <- runif(Narea, 0.00001, 0.0001)
      
      epsilonO = mvrnorm(Nyear, -0.5*varO, diag(Narea)*varO)
      
      IAs = B
      for (yr in 1: (Nyear)){
        for (sp in 1:Narea){
          IAs[yr,sp] <- qs[sp]*B[yr,sp]*exp(epsilonO[yr,sp])
        }
      }
      
      # Biomass estimates: 
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
      
      B_cv <- ifelse(is.na(bio_est),NA,sigB)
      BEs <- bio_est
      
      ### Now making the abundance indices not available for the first years & calculating total catch
      Catch_all <- apply(Catch, 1, sum)
      B_all <- apply(B, 1, sum)
      K_all <- sum(Ks)
      
      IA_obs <- data.frame()
      for (k in 1:(Year_start-1)){
        IA_obs[k,c(1:Narea)] <- NA
      }
      IA_obs[Year_start:Nyear,] <- IAs[Year_start:Nyear,]
      Btrue <- B[1:Nyear,]
      
      # Record Observed contrasts:
      #total contrast
      for(j in 1:Narea){ #j <-1
        sim_stats[iter,paste0("obs_tot_bio_cont",j)] <- (max(BEs[,j], na.rm=T) - min(BEs[,j],na.rm=T))/min(BEs[,j],na.rm=T)
        sim_stats[iter,paste0("obs_reb_bio_cont",j)] <- (tail(na.omit(BEs[,j]), 1) - min(BEs[,j], na.rm=T))/min(BEs[,j], na.rm=T)
        sim_stats[iter,paste0("obs_tot_ind_cont",j)] <- (max(IA_obs[,j], na.rm=T) - min(IA_obs[,j],na.rm=T))/min(IA_obs[,j],na.rm=T)
        sim_stats[iter,paste0("obs_reb_ind_cont",j)] <- (tail(na.omit(IA_obs[,j]), 1) - min(IA_obs[,j], na.rm=T))/min(IA_obs[,j], na.rm=T)[1]
        sim_stats[iter,paste0("Depl_",j)] <- Btrue[Nyear,j]/Btrue[1,j]
      } 
      
      sim_dat[[iter]] <- list("harv_hist" = h_hist,"Hmax" = Hmax,
                              "Btrue" = Btrue,"B_ests" = BEs,"Index" = IA_obs, "Catch" = Catch)
      iter <- iter+1
    }
  }
}

str(sim_stats)
sd(sim_stats$true_tot_cont1)

sim_stats %>% group_by(Hmax, harv_hist) %>%
  dplyr::summarise(true_cont1 = mean(true_tot_cont1),
                   true_cont1_sd = sd(as.numeric(true_tot_cont1),na.rm=T),
                   true_cont2 = mean(true_tot_cont2),
                   true_cont2_sd = sd(true_tot_cont2,na.rm=T),
                   true_cont3 = mean(true_tot_cont3),
                   true_cont3_sd = sd(true_tot_cont3,na.rm=T),
                   true_reb1 = mean(true_reb_cont1),
                   true_reb_sd1 = sd(true_reb_cont1,na.rm=T),
                   true_reb2 = mean(true_reb_cont2),
                   true_reb_sd2 = sd(true_reb_cont2,na.rm=T),
                   true_reb3 = mean(true_reb_cont3),
                   true_reb_sd3 = sd(true_reb_cont3,na.rm=T),
                   obs_bio_cont1 = mean(obs_tot_bio_cont1),
                   obs_bio_cont1_sd = sd(as.numeric(obs_tot_bio_cont1),na.rm=T),
                   obs_bio_cont2 = mean(obs_tot_bio_cont2),
                   obs_bio_cont2_sd = sd(obs_tot_bio_cont2,na.rm=T),
                   obs_bio_cont3 = mean(obs_tot_bio_cont3),
                   obs_bio_cont3_sd = sd(obs_tot_bio_cont3,na.rm=T),
                   obs_reb_bio1 = mean(obs_reb_bio_cont1),
                   obs_reb_bio_sd1 = sd(obs_reb_bio_cont1,na.rm=T),
                   obs_reb_bio2 = mean(obs_reb_bio_cont2),
                   obs_reb_bio_sd2 = sd(obs_reb_bio_cont2,na.rm=T),
                   obs_reb_bio3 = mean(obs_reb_bio_cont3),
                   obs_reb_bio_sd3 = sd(obs_reb_bio_cont3,na.rm=T),
                   depletion_1 = mean(Depl_1),
                   depletion_2 = mean(Depl_2),
                   depletion_3 = mean(Depl_3)
                   ) -> Summary

View(Summary)

library(ggplot2)
library(ggpubr)
library(viridis)
library(ggnewscale)


{ggplot(sim_stats, aes(x = factor(Hmax), y = true_tot_cont1, fill = harv_hist)) + #, fill = factor(cyl))) +
  geom_boxplot(notch = TRUE) +
  labs(title = "Contrast in Area 1",
       x = "Hmax",
       y = "Biomass contrast") +
  theme_minimal() -> cont1

ggplot(sim_stats, aes(x = factor(Hmax), y = true_tot_cont2, fill = harv_hist)) + #, fill = factor(cyl))) +
  geom_boxplot(notch = TRUE) +
  labs(title = "Contrast in Area 2",
       x = "Hmax",
       y = "Biomass contrast") +
  theme_minimal() -> cont2

ggplot(sim_stats, aes(x = factor(Hmax), y = true_tot_cont3, fill = harv_hist)) + #, fill = factor(cyl))) +
  geom_boxplot(notch = TRUE) +
  labs(title = "Contrast in Area 3",
       x = "Hmax",
       y = "Biomass contrast") +
  theme_minimal() -> cont3

ggplot(sim_stats, aes(x = factor(Hmax), y = obs_tot_bio_cont1, fill = harv_hist)) + #, fill = factor(cyl))) +
  geom_boxplot(notch = TRUE) +
  labs(title = "Obs contrast in Area 1",
       x = "Hmax",
       y = "Biomass contrast") +
  theme_minimal() -> ocont1

ggplot(sim_stats, aes(x = factor(Hmax), y = obs_tot_bio_cont2, fill = harv_hist)) + #, fill = factor(cyl))) +
  geom_boxplot(notch = TRUE) +
  labs(title = "Obs contrast in Area 2",
       x = "Hmax",
       y = "Biomass contrast") +
  theme_minimal() -> ocont2

ggplot(sim_stats, aes(x = factor(Hmax), y = obs_tot_bio_cont3, fill = harv_hist)) + #, fill = factor(cyl))) +
  geom_boxplot(notch = TRUE) +
  labs(title = "Obs contrast in Area 3",
       x = "Hmax",
       y = "Biomass contrast") +
  theme_minimal() -> ocont3

ggarrange(cont1,cont2,cont3,ocont1,ocont2,ocont3)}

str(sim_stats)

{ggplot(sim_stats, aes(x = factor(Hmax), y = true_reb_cont1, fill = harv_hist)) + #, fill = factor(cyl))) +
    geom_boxplot(notch = TRUE) +
    labs(title = "Rebound in Area 1",
         x = "Hmax",
         y = "Biomass contrast") +
    theme_minimal() -> reb1
  
  ggplot(sim_stats, aes(x = factor(Hmax), y = true_reb_cont2, fill = harv_hist)) + #, fill = factor(cyl))) +
    geom_boxplot(notch = TRUE) +
    labs(title = "Rebound in Area 2",
         x = "Hmax",
         y = "Biomass contrast") +
    theme_minimal() -> reb2
  
  ggplot(sim_stats, aes(x = factor(Hmax), y = true_reb_cont3, fill = harv_hist)) + #, fill = factor(cyl))) +
    geom_boxplot(notch = TRUE) +
    labs(title = "Rebound in Area 3",
         x = "Hmax",
         y = "Biomass contrast") +
    theme_minimal() -> reb3
  
  ggplot(sim_stats, aes(x = factor(Hmax), y = obs_reb_bio_cont1, fill = harv_hist)) + #, fill = factor(cyl))) +
    geom_boxplot(notch = TRUE) +
    labs(title = "Obs Rebound in Area 1",
         x = "Hmax",
         y = "Biomass contrast") +
    theme_minimal() -> oreb1
  
  ggplot(sim_stats, aes(x = factor(Hmax), y = obs_reb_bio_cont2, fill = harv_hist)) + #, fill = factor(cyl))) +
    geom_boxplot(notch = TRUE) +
    labs(title = "Obs Rebound in Area 2",
         x = "Hmax",
         y = "Biomass contrast") +
    theme_minimal() -> oreb2
  
  ggplot(sim_stats, aes(x = factor(Hmax), y = obs_reb_bio_cont3, fill = harv_hist)) + #, fill = factor(cyl))) +
    geom_boxplot(notch = TRUE) +
    labs(title = "Obs Rebound in Area 3",
         x = "Hmax",
         y = "Biomass contrast") +
    theme_minimal() -> oreb3
  
ggarrange(reb1,reb2,reb3,oreb1,oreb2,oreb3)}

{ggplot(sim_stats, aes(x = factor(Hmax), y = Depl_1, fill = harv_hist)) + #, fill = factor(cyl))) +
    geom_boxplot(notch = TRUE) +
    labs(title = "Stock Status in Area 1",
         x = "Hmax",
         y = "Stock Status (B/B0)") +
    theme_minimal() -> depl1
  
  ggplot(sim_stats, aes(x = factor(Hmax), y = Depl_2, fill = harv_hist)) + #, fill = factor(cyl))) +
    geom_boxplot(notch = TRUE) +
    labs(title = "Stock Status in Area 2",
         x = "Hmax",
         y = "Stock Status (B/B0)") +
    theme_minimal() -> depl2
  
  ggplot(sim_stats, aes(x = factor(Hmax), y = Depl_3, fill = harv_hist)) + #, fill = factor(cyl))) +
    geom_boxplot(notch = TRUE) +
    labs(title = "Stock Status in Area 3",
         x = "Hmax",
         y = "Stock Status (B/B0)") +
    theme_minimal() -> depl3
  
  ggarrange(depl1,depl2,depl3)}


# Make data frames from the master list

length(sim_dat)
for (i in 1:length(sim_dat)) { #i <- 1
  dat <- sim_dat[[i]]
  if (i == 1) {
    bio_true <- as.data.frame(dat$Btrue)
    bio_true$iter <- i
    bio_true$Hmax <- dat$Hmax
    bio_true$harv_hist <- dat$harv_hist
    bio_true$year <- c(1:Nyear)
    
    bio_est <- as.data.frame(dat$B_est)
    bio_est$iter <- i
    bio_est$Hmax <- dat$Hmax
    bio_est$harv_hist <- dat$harv_hist
    bio_est$year <- c(1:Nyear)
    
    catch <- as.data.frame(dat$Catch)
    catch$iter <- i
    catch$Hmax <- dat$Hmax
    catch$harv_hist <- dat$harv_hist
    catch$year <- c(1:(Nyear-1))
  } else {
    bio_true2 <- as.data.frame(dat$Btrue)
    bio_true2$iter <- i
    bio_true2$Hmax <- dat$Hmax
    bio_true2$harv_hist <- dat$harv_hist
    bio_true2$year <- c(1:Nyear)
    bio_true <- rbind(bio_true,bio_true2)
    
    bio_est2 <- as.data.frame(dat$B_est)
    bio_est2$iter <- i
    bio_est2$Hmax <- dat$Hmax
    bio_est2$harv_hist <- dat$harv_hist
    bio_est2$year <- c(1:Nyear)
    bio_est <- rbind(bio_est,bio_est2)
    
    catch2 <- as.data.frame(dat$Catch)
    catch2$iter <- i
    catch2$Hmax <- dat$Hmax
    catch2$harv_hist <- dat$harv_hist
    catch2$year <- c(1:(Nyear-1))
    catch <- rbind(catch,catch2)
  }
}

names(bio_true)<-c("Area1","Area2","Area3","iter","Hmax","harv_hist", "year")
names(bio_est)<-c("Area1","Area2","Area3","iter","Hmax","harv_hist","year")
names(catch)<-c("Area1","Area2","Area3","iter","Hmax","harv_hist","year")

colors1<- viridis(length(sim_dat), begin = 0.7, end = 0.9, option = "B")
colors2<- viridis(length(sim_dat), begin = 0.3, end = 0.5, option = "D")
colors3<- viridis(length(sim_dat), begin = 0.7, end = 0.9, option = "A")

h_pat <- "Harv6"

ggplot(bio_true %>% filter(harv_hist == h_pat)) +
  geom_line(aes(x=year ,y=Area1, col=as.factor(iter)), alpha = 0.2) +
  geom_point(data = bio_est %>% filter(harv_hist == h_pat),
             aes(x=year,y=Area1,col=as.factor(iter),alpha=0.2)) +
  scale_color_manual(values = colors1) +
  new_scale_color()+
  geom_line(aes(x=year ,y=Area2, col=as.factor(iter)), alpha = 0.2) +
  geom_point(data = bio_est %>% filter(harv_hist == h_pat),
             aes(x=year,y=Area2,col=as.factor(iter),alpha=0.2)) +
  scale_color_manual(values = colors2) +
  new_scale_color()+
  geom_line(aes(x=year ,y=Area3, col=as.factor(iter)), alpha = 0.2) +
  geom_point(data = bio_est %>% filter(harv_hist == h_pat),
             aes(x=year,y=Area3,col=as.factor(iter),alpha=0.2)) +
  scale_color_manual(values = colors3) +
  facet_wrap(~Hmax, labeller = "label_both") +
  guides(alpha = "none") +
  theme(legend.position = "none") +
  labs(title = paste0("Biomass; Harvest = ",h_pat),
         x = "Year",
         y = "Biomass")

# plot catch: 
ggplot(catch %>% filter(harv_hist == h_pat)) +
  geom_line(aes(x=year ,y=Area1, col=as.factor(iter)), alpha = 0.2) +
  scale_color_manual(values = colors1) +
  new_scale_color()+
  geom_line(aes(x=year ,y=Area2, col=as.factor(iter)), alpha = 0.2) +
  scale_color_manual(values = colors2) +
  new_scale_color()+
  geom_line(aes(x=year ,y=Area3, col=as.factor(iter)), alpha = 0.2) +
  scale_color_manual(values = colors3) +
  facet_wrap(~Hmax, labeller = "label_both") +
  guides(alpha = "none") +
  theme(legend.position = "none") +
  labs(title = paste0("Catch; Harvest = ",h_pat),
       x = "Year",
       y = "Catch")


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{color1 <- rgb(195/255,33/255,72/255,alpha=0.1)
color2<- rgb(33/255,171/255,205/255,alpha=0.1)
color3<- rgb(135/255,169/255,107/255,alpha=0.1)

par(mfrow=c(1,1), mar=c(4,4,1,1))}
for (i in 1:nsims){
  if (i == 1 ) {
    matplot(subsim[[1]]$Btrue, type="l",ylim = c(0,max(subsim[[1]]$Btrue)*1.1),col=c(color1,color2,color3),
            main = paste0("Hmax: ", Hmax, ", Harvest Pattern = ", H_pat), 
            ylab="Biomass", xlab="Year")
    matpoints(1:Nyear, subsim[[1]]$B_ests, type="p", pch=16, col=c(color1,color2,color3))
    for (i in 1:ncol(subsim[[1]]$B_ests)) {
      arrows(1:Nyear, subsim[[1]]$B_ests[,i] - sigB*subsim[[1]]$B_ests[,i]*1.96, 
             1:Nyear, subsim[[1]]$B_ests[,i] + sigB*subsim[[1]]$B_ests[,i]*1.96, 
             angle=90, code=3, length=0.05, col=c(color1,color2,color3)[i])
    } 
  } else {
    matlines(subsim[[i]]$Btrue, type="l", col=c(color1,color2,color3))
    matpoints(1:Nyear, subsim[[i]]$B_ests, type="p", pch=16, col=c(color1,color2,color3))
  }
  if (i == nsims) {
    mtext(paste0("Depletion = ",round(mean(sim_stats$Depl_1),2),", ",
                 round(mean(sim_stats$Depl_2),2),", ",
                 round(mean(sim_stats$Depl_3),2)),
          side=1, line=0.5, cex=1, adj = 0, padj=-14)
    mtext(paste0("True Contrast = ",round(mean(sim_stats$true_tot_cont1),2),", ",
                 round(mean(sim_stats$true_tot_cont2),2),", ",
                 round(mean(sim_stats$true_tot_cont3),2)),
          side=1, line=0.5, cex=1, adj = 0, padj=-12)
    mtext(paste0("Obs Contrast = ",round(mean(sim_stats$obs_tot_bio_cont1),2),", ",
                 round(mean(sim_stats$obs_tot_bio_cont2),2),", ",
                 round(mean(sim_stats$obs_tot_bio_cont3),2)),
          side=1, line=0.5, cex=1, adj = 0, padj=-10)
    mtext(paste0("True rebound = ",round(mean(sim_stats$true_reb_cont1),2),", ",
                 round(mean(sim_stats$true_reb_cont2),2),", ",
                 round(mean(sim_stats$true_reb_cont3),2)),
          side=1, line=0.5, cex=1, adj = 0, padj=-8)
    mtext(paste0("Obs rebound = ",round(mean(sim_stats$obs_reb_bio_cont1),2),", ",
                 round(mean(sim_stats$obs_reb_bio_cont2),2),", ",
                 round(mean(sim_stats$obs_reb_bio_cont3),2)),
          side=1, line=0.5, cex=1, adj = 0, padj=-6)
  }
}

#Plot index:
for (i in 1:nsims){
  if (i == 1 ) {
    matplot(sim_dat[[1]]$Index, type="l",ylim = c(0,max(sim_dat[[1]]$Index,na.rm=T)*2),
            col=c(color1,color2,color3),
            main = paste0("Hmax: ", Hmax, ", Harvest Pattern = ", h_hist), 
            ylab="Index", xlab="Year")
    }  else {
    matlines(sim_dat[[i]]$Index, type="l", col=c(color1,color2,color3))
  }
}

#Plot catch:
for (i in 1:nsims){
  if (i == 1 ) {
    matplot(sim_dat[[1]]$Catch, type="l",ylim = c(0,max(sim_dat[[1]]$Catch,na.rm=T)*1.5),col=c(color1,color2,color3),
            main = paste0("Hmax: ", Hmax, ", Harvest Pattern = ", h_hist), 
            ylab="Catch", xlab="Year")
  }  else {
    matlines(sim_dat[[i]]$Catch, type="l", col=c(color1,color2,color3))
  }
}
  
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Prepping Data for STAN: 
# Vectorize the index and biomass estimates to deal with incomplete data:
  B_obs <- as.vector(BEs[!is.na(BEs)])
  B_cv <- as.vector(B_cv[!is.na(B_cv)])
  N_Bobs <- sum (!is.na(BEs))
  B_pos <- list()
  for (i in 1:Narea) {
    B_pos[[i]] <- which(!is.na(BEs[,i]))
  }
  B_pos <- as.vector(unlist(B_pos))
  S_Bobs <- vector()
  for (i in 1:Narea) {
    S_Bobs[i] <- length(BEs[,i][!is.na(BEs[,i])])
  }

  I_obs <- as.vector(IA_obs[!is.na(IA_obs)])
  
  I_cv <- ifelse(is.na(IA_obs),NA,sigO)
  I_cv <- as.vector(I_cv[!is.na(I_cv)])
  N_Iobs <- sum (!is.na(IA_obs))
  I_pos <- list()
  for (i in 1:Narea) {
    I_pos[[i]] <- which(!is.na(IA_obs[,i]))
  }
  I_pos <- as.vector(unlist(I_pos))
  S_Iobs <- vector()
  for (i in 1:Narea) {
    S_Iobs[i] <- length(IA_obs[,i][!is.na(IA_obs[,i])])
  }

### If using STan

library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

iters <- 2000
chains <- 3
burnin <- 0.6 #proportion of chain length used as warmup/burnin

adapt_delta <- 0.99
stepsize <- 0.01
max_treedepth <- 15
stepsize_jitter <- 0

## Play around with informative priors, start at true, add more noisiness
if (Narea == 1){
  data = list(Npre = length(Year_start:Nyear), 
              S=Narea, 
              p = 0.18815,
              C_obs = as.vector(Catch[Year_start:(Nyear-1),]),
              I_obs=as.vector(IA_obs), 
              ratio = 1,
              #observed biomass
              B_obs = B_obs, 
              B_cv = B_cv,
              N_Bobs = N_Bobs,
              B_pos = B_pos,
              S_Bobs = S_Bobs, 
              # Index 1
              I_obs = I_obs,
              I_cv = I_cv, 
              N_1obs = N_Iobs,
              I_pos = I_pos,
              S_Iobs = S_Iobs
              )
  
  init_create <- function(chain_id=1){
    set.seed(chain_id+123)
    inits <- list(
      # MSY = 15000*apply( data$C_obs, 1, mean)/apply( data$C_obs, 1, mean)[1],
      # MSY = MSYs,
      r = runif(data$S, 0.01, 0.1),
      # PP_init = runif(4, 0.8, 1),
      PP_init = runif(data$S, 0.8, 1),
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
              I_obs=t(IA_obs), ratio = 1)
  
  init_create <- function(chain_id=1){
    set.seed(chain_id+123)
    inits <- list(
    # MSY = 15000*apply( data$C_obs, 1, mean)/apply( data$C_obs, 1, mean)[1],
    # MSY = MSYs,
    r = runif(data$S, 0.01, 0.1),
    # PP_init = runif(4, 0.8, 1),
    PP_init = runif(data$S, 0.8, 1),
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

tstart <- Sys.time()
fit <- stan(file = paste0("Production_models_STAN/Kotaro/kotaro.stan"), 
            data = data, init = init_ll, #inits, inits),
            iter = iters, chains = chains, cores=chains, seed=123,
            warmup=burnin*iters, verbose=F, thin=1,
            control = list(adapt_delta = adapt_delta, stepsize = stepsize, 
                           max_treedepth = max_treedepth, stepsize_jitter = stepsize_jitter)) #stepsize_jitter default(0), values between 0 and 1
#metric (string, one of "unit_e", "diag_e", "dense_e", defaults to "diag_e")
runtime <- Sys.time() - tstart; runtime

fit
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
launch_shinystan(fit)


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



