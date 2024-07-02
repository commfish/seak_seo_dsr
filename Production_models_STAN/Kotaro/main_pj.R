rm(list=ls())
gc()

library(MASS) ## needed for the multivariate normal distribution

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
# Alt 1 Kotaro's 
H1 = sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-20), seq(Hmax, 0.5, length.out=20))*HMSYs[x], 0.1*HMSYs[x]))
matplot(1:Nyear, H1, type="l", ylim = c(0,max(H1)))
    # Alt 1: 
H2 = sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-C_br1), seq(Hmax,0.1*Hmax,length.out=(C_br2-C_br1)), seq(Hmax*0.1, 0.5, length.out=Nyear-C_br2))*HMSYs[x], 0.1*HMSYs[x]))
matplot(1:Nyear, H2, type="l", ylim = c(0,max(H2)))
    # Alt 2
H3 = sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-15), seq(Hmax*0.5, 0.5, length.out=15))*HMSYs[x], 0.1*HMSYs[x]))
matplot(1:Nyear, H3, type="l", ylim = c(0,max(H3)))
    # Alt 3
H4 = sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-C_br1), seq(1,0.1*Hmax,length.out=(C_br2-C_br1)), seq(Hmax*0.1, 0.5, length.out=Nyear-C_br2))*HMSYs[x], 0.1*HMSYs[x]))
matplot(1:Nyear, H4, type="l", ylim = c(0,max(H4)))

Hlist <- list(H1,H2,H3,H4)

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

# Pick a harvest strategy
H<-H3

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

  matplot(1:Nyear, B, type="l", ylim = c(0,max(B)))

  #True contrasts:
  #total contrast
  sapply(1:Narea, function(x) (max(B[,x]) - min(B[,x]))/min(B[,x]))
  
  #rebound contrast; difference between low point and end point
  sapply(1:Narea, function(x) (B[Nyear,x] - min(B[,x]))/min(B[,x]))
  
  #diff between start and end
  sapply(1:Narea, function(x) (B[1,x] - B[Nyear,x])/B[Nyear,x])
    

### Now trying to create some abundance indices
  qs <- runif(Narea, 0.00001, 0.0001)
 
  epsilonO = mvrnorm(Nyear, -0.5*varO, diag(Narea)*varO)

  IAs = B
  for (yr in 1: (Nyear)){
    for (sp in 1:Narea){
      IAs[yr,sp] <- qs[sp]*B[yr,sp]*exp(epsilonO[yr,sp])
    }
  }

  matplot(1:Nyear, IAs, type="l")
  
# Biomass estimates:   
  
  sigB = 0.15
  varB = rep(sigB^2, Narea)
  epsilonB = mvrnorm(Nyear, varB, diag(Narea)*varB)
  
  BEs = B
  
  for (yr in 1: (Nyear)){
    for (sp in 1:Narea){
      BEs[yr,sp] <- B[yr,sp]*exp(epsilonB[yr,sp])
    }
  }

  matplot(1:Nyear, B, type="l", ylim = c(0,max(B)*2))
  matpoints(1:Nyear, BEs, type="p", pch=16)
  for (i in 1:ncol(BEs)) {
    arrows(1:Nyear, BEs[,i] - sigB*BEs[,i]*1.96, 1:Nyear, BEs[,i] + sigB*BEs[,i]*1.96, angle=90, code=3, length=0.05, col=c("black", "red", "green")[i])
    #polygon(c(BEs[,i], rev(BEs[,i])), c(BEs[,i] + sigB*BEs[,i]*1.96, rev(BEs[,i] - sigB*BEs[,i]*1.96)), col = c("black", "red", "green")[i], border = NA)
  }
  
# But each area is only visited every third year starting in the 11th year of the fishery:
  bio_srv_frq <- 3
  bio_srv_start <- 11
  
  bio_est <- data.frame()
  for (i in 1:Nyear){
    bio_est[i,c(1:Narea)] <- NA
  }
  for (i in 1:Narea){
    bio_est[seq(bio_srv_start+i,Nyear,bio_srv_frq),i] <- 
      BEs[seq(bio_srv_start+i,Nyear,bio_srv_frq),i]
  }
  
  B_cv <- ifelse(is.na(bio_est),NA,sigB)
  BEs <- bio_est
  
  matplot(1:Nyear, B, type="l", ylim = c(0,max(B)*1.1))
  matpoints(1:Nyear, BEs, type="p", pch=16)
  for (i in 1:ncol(BEs)) {
    arrows(1:Nyear, BEs[,i] - sigB*BEs[,i]*1.96, 1:Nyear, BEs[,i] + sigB*BEs[,i]*1.96, angle=90, code=3, length=0.05, col=c("black", "red", "green")[i])
  }
  
### Now making the abundance indices not available for the first years & calculating total catch
  Catch_all <- apply(Catch, 1, sum)
  B_all <- apply(B, 1, sum)
  K_all <- sum(Ks)
  
  IA_obs <- IAs[Year_start:Nyear,]
  matplot(IA_obs, type="l")
  Btrue <- B[Year_start:Nyear,]
  matplot(Btrue, type="l")
  # matplot(Catch[Year_start:Nyear,], type="l")
  Btrue[1,]/Ks
  
  par(mfrow=c(2,2), mar=c(2,4,1,1))
  matplot(Btrue, type="l",ylim = c(0,max(Btrue)*1.1))
  matpoints(1:Nyear, BEs, type="p", pch=16)
  for (i in 1:ncol(BEs)) {
    arrows(1:Nyear, BEs[,i] - sigB*BEs[,i]*1.96, 1:Nyear, BEs[,i] + sigB*BEs[,i]*1.96, angle=90, code=3, length=0.05, col=c("black", "red", "green")[i])
  }
  matplot(epsilon[Year_start:Nyear,], type="l")
  matplot(IA_obs, type="l")
  matplot(Catch[Year_start:(Nyear-1),], type="l")

#------------------------------------------------------------------------------  
# Lets generate a gazzilion simulations to get a sense of how these are doing
sim_dat <- list()
sim_stats <-data.frame()
nsims <- 10
h_hist <- "opt1" #opt2, opt3, opt4, opt_mix
Hmax <- 3

for (i in 1:nsims){ #i <- 1
  if (h_hist == "opt1") {
    H = sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-20), seq(Hmax, 0.5, length.out=20))*HMSYs[x], 0.1*HMSYs[x]))
  }
  
  if (h_hist == "opt2") {
    H = sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-C_br1), seq(Hmax,0.1*Hmax,length.out=(C_br2-C_br1)), seq(Hmax*0.1, 0.5, length.out=Nyear-C_br2))*HMSYs[x], 0.1*HMSYs[x]))
  }
  
  if (h_hist == "opt3") {
    H = sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-15), seq(Hmax*0.5, 0.5, length.out=15))*HMSYs[x], 0.1*HMSYs[x]))
  }
  
  if (h_hist == "opt4") {
    H = sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-C_br1), seq(1,0.1*Hmax,length.out=(C_br2-C_br1)), seq(Hmax*0.5, 0.5, length.out=Nyear-C_br2))*HMSYs[x], 0.1*HMSYs[x]))
  }
  
  if (h_hist == "opt_mix") {
    H1 <- sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-20), seq(Hmax, 0.5, length.out=20))*HMSYs[x], 0.1*HMSYs[x]))
    H2 <- sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-C_br1), seq(Hmax,0.1*Hmax,length.out=(C_br2-C_br1)), seq(Hmax*0.1, 0.5, length.out=Nyear-C_br2))*HMSYs[x], 0.1*HMSYs[x]))
    H3 <- sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-15), seq(Hmax*0.5, 0.5, length.out=15))*HMSYs[x], 0.1*HMSYs[x]))
    H4 <- sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-C_br1), seq(1,0.1*Hmax,length.out=(C_br2-C_br1)), seq(Hmax*0.5, 0.5, length.out=Nyear-C_br2))*HMSYs[x], 0.1*HMSYs[x]))
    
    Hlist <- list(H1,H2,H3,H4)
    
    Hmix <- data.frame(matrix(nrow=Nyear))
    for (i in 1:Narea){
      Hx <- Hlist[[sample(length(Hlist), 1)]]
      H[,i] <- Hx[,sample(1:Narea,1)]
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
  sim_stats[i,"iter"] <- i
  sim_stats[i,"harv_hist"] <- h_hist
  sim_stats[i,"Hmax"] <- Hmax
  for(j in 1:Narea){ #j <-1
    #total contrast
    sim_stats[i,paste0("true_tot_cont",j)] <- (max(B[,j], na.rm=T) - min(B[,j],na.rm=T))/min(B[,j],na.rm=T)
    #rebound of the population from the low point
    sim_stats[i,paste0("true_reb_cont",j)] <- (tail(na.omit(B[,j]), 1) - min(B[,j], na.rm=T))/min(B[,j], na.rm=T)
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
  epsilonB = mvrnorm(Nyear, -0.5*varB, diag(Narea)*varB)
  
  logB <- log(B)
  
  logBEs = B
  BEs <- B
  BElo <- B
  BEhi <- B
  
  for (yr in 1: (Nyear)){
    for (sp in 1:Narea){
      logBEs[yr,sp] <- rlnorm(1,logB[yr,sp],sigB * exp(logB[yr,sp]))
      BEs[yr,sp] <- B[yr,sp]*logBEs[yr,sp]
      BElo[yr,sp] <- exp(logB[yr,sp] - 1.96 * sigB * logB[yr,sp])
      BEhi[yr,sp] <- exp(logB[yr,sp] + 1.96 * sigB * logB[yr,sp])
      #BEs[yr,sp] <- B[yr,sp]*exp(epsilonB[yr,sp])
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
  #matplot(IA_obs, type="l")
  Btrue <- B[Year_start:Nyear,]
  #matplot(Btrue, type="l")
  # matplot(Catch[Year_start:Nyear,], type="l")
  #Btrue[1,]/Ks
  
  # Record Observed contrasts:
  #total contrast
  for(j in 1:Narea){ #j <-1
    sim_stats[i,paste0("obs_tot_bio_cont",j)] <- (max(BEs[,j], na.rm=T) - min(BEs[,j],na.rm=T))/min(BEs[,j],na.rm=T)
    sim_stats[i,paste0("obs_reb_bio_cont",j)] <- (tail(na.omit(BEs[,j]), 1) - min(BEs[,j], na.rm=T))/min(BEs[,j], na.rm=T)
    sim_stats[i,paste0("obs_tot_ind_cont",j)] <- (max(IA_obs[,j], na.rm=T) - min(IA_obs[,j],na.rm=T))/min(IA_obs[,j],na.rm=T)
    sim_stats[i,paste0("obs_reb_ind_cont",j)] <- (tail(na.omit(IA_obs[,j]), 1) - min(IA_obs[,j], na.rm=T))/min(IA_obs[,j], na.rm=T)[1]
  } 
  
  if (i == 1) {
    par(mfrow=c(1,1), mar=c(2,4,1,1))
    matplot(Btrue, type="l",ylim = c(0,max(Btrue)*1.1))
    matpoints(1:Nyear, BEs, type="p", pch=16, col=c("black", "red", "forestgreen", alpha=0.05))
    for (i in 1:ncol(BEs)) {
      arrows(1:Nyear, BEs[,i] - sigB*BEs[,i]*1.96, 1:Nyear, BEs[,i] + sigB*BEs[,i]*1.96, angle=90, code=3, length=0.05, col=c("black", "red", "forestgreen", alpha=0.05)[i])
    }
  } else {
    matlines(Btrue, type="l", col=c("black", "red", "forestgreen", alpha=0.05))
    matpoints(1:Nyear, BEs, type="p", pch=16, col=c("black", "red", "forestgreen", alpha=0.05))
  }
  
  sim_dat[["Btrue"]] <- Btrue
  sim_dat[["Bobs"]] <- BEs
  sim_dat[["Iobs"]] <- IA_obs
  sim_dat[["Catch"]] <- Catch
 
  
  
  par(mfrow=c(2,2), mar=c(2,4,1,1))
  matplot(Btrue, type="l",ylim = c(0,max(Btrue)*1.1))
  matpoints(1:Nyear, BEs, type="p", pch=16)
  for (i in 1:ncol(BEs)) {
    arrows(1:Nyear, BEs[,i] - sigB*BEs[,i]*1.96, 1:Nyear, BEs[,i] + sigB*BEs[,i]*1.96, angle=90, code=3, length=0.05, col=c("black", "red", "forestgreen")[i])
  }
  matplot(epsilon[Year_start:Nyear,], type="l")
  matplot(IA_obs, type="l")
  matplot(Catch[Year_start:(Nyear-1),], type="l")
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



