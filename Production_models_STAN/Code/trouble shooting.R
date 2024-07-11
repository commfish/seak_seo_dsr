sim_stats <-data.frame()
nsims <- 5
data_iter <- 1
for (i in 1:nsims){ #i <- 1
  #rs = runif(Narea, 0.045, 0.055)    # variability in rs, but will estimate on r in model data_iter <- data_iter+1
  # set.seed(123)
  rs = rlnorm(Narea,log(0.05),0.025)
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
  
  C_br1 <- 25
  C_br2 <- 40
  
  for(j in 1:Narea){ #j <-1
    sim_stats[data_iter,"iter"] <- data_iter
    sim_stats[data_iter,"harv_hist"] <- h_hist
    sim_stats[data_iter,"Hmax"] <- Hmax
    sim_stats[data_iter,paste0("r_",j)] <- rs[j]
    sim_stats[data_iter,paste0("K_",j)] <- Ks[j]
    sim_stats[data_iter,paste0("MSY_t_",j)] <- MSYs[j]
    sim_stats[data_iter,paste0("HMSYs_t_",j)] <- HMSYs[j]
    sim_stats[data_iter,paste0("Bmsy_t_",j)] <- 0.4 * Ks[j]
    sim_stats[data_iter,paste0("Fmsy_t_",j)] <- MSYs[j] / sim_stats[data_iter,paste0("Bmsy_t_",j)]
  } 
  
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
    #H6 <- sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.1, 0.5*Hmax, length.out=Nyear-C_br1), seq(0.5*Hmax,0.75*Hmax,length.out=(C_br2-C_br1)), seq(Hmax*0.75, 0.7*Hmax, length.out=Nyear-C_br2))*HMSYs[x], 0.1*HMSYs[x]))
    Hlist <- list(H1,H2,H3,H4,H5)
    
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
      B[yr+1,sp] = (max( B[yr,sp] + rs[sp]*B[yr,sp]/p*(1 - (B[yr,sp]/Ks[sp])^p) - B[yr,sp]*H[yr,sp], 1))*exp(epsilon[yr,sp])
      Catch[yr,sp] = B[yr,sp]*H[yr,sp]
    }}
  
  ## abundance indices
  qs <- runif(Narea, 0.00001, 0.0001)
  
  for(j in 1:Narea){ #j <-1
    sim_stats[data_iter,paste0("true_term_bio",j)] <- B[Nyear,j]
    sim_stats[data_iter,paste0("true_ss",j)] <- B[Nyear,j] / Ks[j]
    #total contrast
    sim_stats[data_iter,paste0("true_tot_cont",j)] <- (max(B[,j], na.rm=T) - min(B[,j],na.rm=T))/min(B[,j],na.rm=T)
    #rebound of the population from the low point
    sim_stats[data_iter,paste0("true_reb_cont",j)] <- (tail(na.omit(B[,j]), 1) - min(B[,j], na.rm=T))/min(B[,j], na.rm=T)
    # number of years between low point and terminal year
    sim_stats[data_iter,paste0("true_reb_len",j)] <- Nyear - which(B[,j] == min(B[,j]))
    sim_stats[data_iter,paste0("q",j)] <- qs[j]
  } 
  
  epsilonO = mvrnorm(Nyear, -0.5*varO, diag(Narea)*varO)
  
  IAs = sim_IAs = B
  for (yr in 1: (Nyear)){
    for (sp in 1:Narea){
      IAs[yr,sp] <- qs[sp]*B[yr,sp]*exp(epsilonO[yr,sp])
      sim_IAs[yr,sp] <- qs[sp]*B[yr,sp]
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
  
  B_cvs <- ifelse(is.na(bio_est),NA,sigB)
  B_ests <- bio_est
  
  #need to vectorize the biomass estimates:B_obs <- as.vector(data$B_ests[!is.na(data$B_ests)])
  B_obs <- as.vector(B_ests[!is.na(B_ests)])
  B_cv <- as.vector(B_cvs[!is.na(B_cvs)])
  N_Bobs <- sum (!is.na(B_ests))
  B_pos <- list()
  for (as in 1:Narea) {
    B_pos[[as]] <- which(!is.na(B_ests[,as]))
  }
  B_pos <- as.vector(unlist(B_pos))
  S_Bobs <- vector()
  for (as in 1:Narea) {
    S_Bobs[as] <- length(B_ests[,as][!is.na(B_ests[,as])])
  }
  
  ### Now making the abundance indices not available for the first years & calculating total catch
  Catch_all <- apply(Catch, 1, sum)
  B_all <- apply(B, 1, sum)
  K_all <- sum(Ks)
  
  #Year_start <- 1
  IA_obs <- I_obs_plot <- IAs[1:Nyear,] 
  #matplot(IA_obs, type="l")
  Btrue <- B[1:Nyear,]
  #matplot(Btrue, type="l")
  # matplot(Catch[Year_start:Nyear,], type="l")
  #Btrue[1,]/Ks
  
  #par(mfrow=c(4,1), mar=c(2,4,1,1))
  #matplot(Btrue, type="l")
  #matplot(epsilon[1:Nyear,], type="l")
  #matplot(IA_obs, type="l")
  #matplot(Catch[Year_start:(Nyear-1),], type="l")
  
  #Eventually vectorize the index...
  #try getting rid of first five years of index data:
  cpue_start <- cpue_start
  
  cpue2 <- data.frame()
  for (ys in 1:Nyear){
    cpue2[ys,c(1:Narea)] <- NA
  }
  for (ys in 1:Narea){
    cpue2[seq(cpue_start,Nyear,1),ys] <- 
      IA_obs[seq(cpue_start,Nyear,1),ys]
  }
  
  #I_obs <- as.vector(IA_obs[!is.na(IA_obs)])
  I_obs_mat <- IAs[cpue_start:Nyear,] 
  I_obs <- as.vector(cpue2[!is.na(cpue2)])
  I_cv <- I_cv_plot <- ifelse(is.na(cpue2),NA,sigO) 
  I_cv <- as.vector(I_cv[!is.na(I_cv)])
  N_Iobs <- sum (!is.na(I_obs))
  I_pos <- list()
  for (as in 1:Narea) {
    I_pos[[as]] <- which(!is.na(cpue2[,as]))
  }
  I_pos <- as.vector(unlist(I_pos))
  S_Iobs <- vector()
  for (as in 1:Narea) {
    S_Iobs[as] <- length(cpue2[,as][!is.na(cpue2[,as])])
  }
  
  # Record Observed contrasts:
  # total contrast
  for(j in 1:Narea){ #j <-1
    sim_stats[data_iter,paste0("obs_tot_bio_cont",j)] <- (max(BEs[,j], na.rm=T) - min(BEs[,j],na.rm=T))/min(BEs[,j],na.rm=T)
    sim_stats[data_iter,paste0("obs_reb_bio_cont",j)] <- (tail(na.omit(BEs[,j]), 1) - min(BEs[,j], na.rm=T))/min(BEs[,j], na.rm=T)
    sim_stats[data_iter,paste0("obs_reb_bio_len",j)] <- Nyear - which(BEs[,j] == min(BEs[,j], na.rm=T))
    sim_stats[data_iter,paste0("obs_tot_ind_cont",j)] <- (max(IA_obs[,j], na.rm=T) - min(IA_obs[,j],na.rm=T))/min(IA_obs[,j],na.rm=T)
    sim_stats[data_iter,paste0("obs_reb_ind_cont",j)] <- (tail(na.omit(IA_obs[,j]), 1) - min(IA_obs[,j], na.rm=T))/min(IA_obs[,j], na.rm=T)[1]
    sim_stats[data_iter,paste0("obs_reb_ind_len",j)] <- Nyear - which(IA_obs[,j] == min(IA_obs[,j], na.rm=T))
    sim_stats[data_iter,paste0("Depl_",j)] <- Btrue[Nyear,j]/Btrue[1,j]
  } 
  
  sim_dat[[data_iter]] <- list("harv_hist" = h_hist,"Hmax" = Hmax,
                               "Btrue" = Btrue,"B_ests" = BEs,"Index" = IA_obs, "Catch" = Catch)
  
  saveRDS(sim_dat, file = "Production_models_STAN/Output/sim_res/simulated_data.Rds")
  saveRDS(sim_dat, file = "H://Documents/SEO_DSR/stan_development/sim_backup/simulated_data.Rds")
  write.csv(sim_stats,'Production_models_STAN/Output/sim_res/sim_stats.csv')
  write.csv(sim_stats,'H://Documents/SEO_DSR/stan_development/sim_backup/sim_stats.csv')
  
  stan_iters <- 2000
  chains <- 3
  burnin <- 0.5 #proportion of chain length used as warmup/burnin
  
  adapt_delta <- 0.99
  stepsize <- 0.01
  max_treedepth <- 15
  stepsize_jitter <- 0
  
  data1 = list(Npre = length(cpue_start:Nyear), 
               N = Nyear, 
               S=Narea, 
               p = 0.18815,
               C_obs = t(Catch[1:(Nyear-1),]),
               I_obs=t(I_obs_mat), #t(IA_obs), 
               ratio = 1,
               N_Bobs = N_Bobs,
               B_obs = B_obs,
               B_cv = B_cv,
               S_Bobs = S_Bobs,
               B_pos = B_pos)
  
  data2 = list(Npre = length(cpue_start:Nyear), 
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
  
  #init_ll_3r <- lapply(1:chains, function(id) init_create_3r(chain_id = id))
  init_ll_1r <- lapply(1:chains, function(id) init_create_1r(chain_id = id))
  # Run Index only, assumed PE = OE, complete index
  
  
  data_iter <- data_iter + 1
}

sim_stats
