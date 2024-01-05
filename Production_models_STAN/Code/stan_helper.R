################################################################################
## Canned functions for SS-SPM

#-------------------------------------------------------------------------------
# get pella-tomlinson exponent for a desired Bmsy point:

pella_toml <- function(pmsy = 0.46) {
  pt_p<-function(p){
    return(exp((1/p)*log(1/(p+1)))-pmsy)
  }
  pte_p <- uniroot(pt_p,c(0.01,5))$root
  
  pt_m<-function(m){
    return(m ^ (-1/(m-1))-pmsy)
  }
  
  pte_m <- uniroot(pt_m,c(0.5,5))$root
  
  res<-list(pte_m, pte_p)
  names(res) <- c("pte_m","pte_p")
  return(res)
}

#pella_toml(0.46)


#--------------------------------------------------------------------------------
# Create simulated data for simulation testing of models

simulate <- function(years = years, strata=strata, r=r, K=K, rho=rho,
                     tau_b1 = tau_b1, tau_b2 = tau_b2, 
                     tau_f1 = tau_f1, tau_f2 = tau_f2,
                     depletion = depletion,
                     by_to_C = by_to_C #bycatch to catch ratio: is the fishery a bycatchfishery, directed fishery or some mix
) 
{ 
  pi <- c(runif(strata-1,0.1,0.5)); pi <- c(pi,1-sum(pi))#c(0.3,0.1,0.4,0.2)
  q1 <- runif(strata,0.015, 0.025) #c(0.025,0.02,0.025,0.02)  #fishing q
  q2 <- runif(strata, 0.001,0.004)#c(0.006,0.0020,0.001,0.008) #survey q
  q3 <- runif(strata,0.0001,0.001) #discard q
  
  # observation error terms for simulating... 
  q_err <- c(rep(0.15 + tau_f1,years-15),rep(0.15 + tau_f2, 15)) #fishery catchability deviation
  d_err <- c(rep(0.5,(years - 20)), rep(0.05,(20))) #discard measurement error, bad estimates then well documented
  b_err <- c(rep(0.15 + tau_b1,years-15),rep(0.15 + tau_b2,15)); length(b_err) #biomass estimate uncertainty
  by_err <- c(rep(0.15, years)) # bycatch uncertainty
  s_err <- c(rep(0.15,years-20),rep(0.15,20)); length(s_err) #survey uncertainty
  
  #P1: where the population is in year 1 relative to virgin biomass
  #depletion <- depletion #population is 30% depleted; ie at B70
  
  #fishery F start low, go up, decline, close it and then low fishing pressure
  f_shape <- c(rep(0.25 * r, 0.1 * years), # early phase of low fishing pressure)
               seq(0.3 * r, 2 * r, (2 * r - 0.25 * r) / (0.2 * years)), # building fishery 
               seq(2 * r, 0.5 * r, (0.5 * r - 2 * r) / (0.3 * years)), # decreasing fishery
               rep(0.5 * r, 5), #slow down fishing for a few years
               rep(0.1 * r, 5), #slow it down even more
               rep(0.35 * r, 5))
  f_shape <- c(f_shape ,rep(0.25 * r, years - length(f_shape))); length(f_shape)
  
  # lets offset that curve so that the strata have slightly different histories
  f_shape2 <- c(rep(0.5 * r, 0.1 * years),
                f_shape[1:(0.9 * years)] * 1.3)
  f_shape2 <- c(f_shape2 ,rep(0.25 * r, years - length(f_shape2))); length(f_shape2)
  
  f_shape3 <- c(rep(0.05 * r, 0.1 * years),
                f_shape[1:(0.9 * years)] * 1.2)
  f_shape3 <- c(f_shape3 ,rep(0.25 * r, years - length(f_shape3))); length(f_shape3)
  
  f <- cbind(f_shape, f_shape2, f_shape3)
  
  f <- f * by_to_C[2] 
  
  effort <- f/ q1

  #d_effort <- d / q3
  
  #bycatch fishing F method 2
  # modelled such that constant bycatch rate that goes up and down over the course of the time series... 
  by_shape <- c(seq(0.1 * r, 0.15 * r, (0.15 * r - 0.1 * r) / (0.3 * years)), # building fishery 
                seq(0.15 * r, 0.10 * r, (0.1 * r - 0.15 * r) / (0.5 * years)), # decreasing fishery
                rep(0.1 * r, 2), #slow down fishing for a few years
                rep(0.15 * r, 3))
  by_shape <- c(by_shape ,rep(0.06 * r, years - length(by_shape))); length(by_shape)
  
  # lets offset that curve so that the strata have slightly different histories
  by_shape2 <- by_shape * runif(length(by_shape),0.8,1.2); length(by_shape2)
  by_shape3 <- by_shape * runif(length(by_shape),0.8,1.2); length(by_shape3)
  
  by <- cbind(by_shape, by_shape2, by_shape3)
  
  by <- by * by_to_C[1] * 10
  
  by_effort <- by / q3
  
  # retention and landing of bycatch.. for simple 0 retention in early years and then 95% landed
  retention <- c(rep(0, 0.6 * years), # early phase of low fishing pressure)
                 rep(0.95, 0.4 * years)  ); length(retention)
  #retention <- c(f_shape ,rep(0.25 * r, years - length(f_shape))); length(f_shape)
  
  #process error with autocorrelation:  
  pesd<- 0.5 * r #0.02
  rho <- 0.75
  ar<-scale(arima.sim(n=years-1, list(ar=rho), innov=rnorm(years-1))) * pesd
  pe<-c(0,ar)
  
  C<-matrix(nrow=years, ncol=strata)
  D<-matrix(nrow=years, ncol=strata)
  d_est<-matrix(nrow=years, ncol=strata)
  By<-matrix(nrow=years, ncol=strata)
  LBy<-matrix(nrow=years, ncol=strata)
  Byc<-matrix(nrow=years, ncol=strata)
  C_tot<-matrix(nrow=years, ncol=strata)
  
  bio<-matrix(nrow=years, ncol=strata)
  bio[1,]<-pi*K*(1) * (1 - depletion)
  bio_est<-matrix(nrow=years, ncol=strata)
  bio_est[1,]<-bio[1,]*rnorm(strata,1,b_err)
  
  prop<-matrix(nrow=years, ncol=strata)
  prop[1,]<-1 * (1 - depletion)
  prop_est<-matrix(nrow=years, ncol=strata)
  prop_est[1,]<-1 
  
  prop2<-matrix(nrow=years, ncol=strata)
  prop2[1,]<-1 * (1 - depletion)
  
  C[1,]<- f[1,] * bio[1,] #- (exp(f_shape[1]) + 1) #bio[1,]*q1*rnorm(strata,1,q_err)*f[1]
  #D[1,]<- d[1,] * bio[1,] #bio[1,]*q3*rnorm(strata,1,d_err)*df[1]
  By[1,]<-by[1,] * bio[1,]
  LBy[1,] <- By[1,] * retention[1]
  D[1,]<- By[1,] -  LBy[1,]
  d_est[1,] <- D[1,] * rnorm(strata,1,d_err)
  C_tot[1,] <- C[1,] + LBy[1,]
  
  cpue<-matrix(nrow=years, ncol=strata)
  cpue[1,]<-C[1,] / effort[1,]
  
  surv<-matrix(nrow=years,ncol=strata)
  surv[1,]<-q2*bio[1,]
  
  for (s in 1:strata){  #s<-1
    for (i in 2:years){ #i<-2
      bio[i,s] <- (bio[i-1,s] + (r/pte_p)*bio[i-1,s]*(1-(bio[i-1,s]/(K*pi[s]))^pte_p) - C[i-1,s] - By[i-1,s])*(1+pe[i])#bio[i-1]*q*rnorm(1,1,0.1)
      prop[i,s] <- (prop[i-1,s] + r/(pte_m-1) * prop[i-1,s] * (1 - prop[i-1,s]^(pte_m-1)) - (C[i-1, s] + By[i-1,s])/(K*pi[s]))*(1+pe[i])
      prop2[i,s] <- bio[i,s] / (K * pi[s])
      #q1[s]*rnorm(1,1,q_err[i])*bio[i,s]*f[i,s]
      By[i,s] <- max(0, by[i,s] * bio[i,s] * rnorm(1,1,by_err[i]))
      LBy[i,s] <- retention[i] * By[i,s]
      C[i,s] <- max(0, f[i,s] * bio[i,s] * rnorm(1,1,q_err[i])) 
      C_tot[i,s] <- C[i,s] + LBy[i,s]
      D[i,s] <- By[i,s] - LBy[i,s] #max(0, d[i,s] * bio[i,s] * rnorm(1,1,d_err[i]))#max(q3[s]*rnorm(1,1,d_err[i])*bio[i,s]*df[i,s],0)
      # estimates
      bio_est[i,s] <- bio[i,s]*rlnorm(1,0,b_err[i])
      d_est[i,s] <- D[i,s] * rlnorm(1,0,d_err[i])
      
      surv[i,s] <- bio[i,s] * q2[s] * rlnorm(1,0,s_err[i])
      if (f[i,s] == 0) {
        cpue[i,s]<- NA
      } else {
        cpue[i,s]<- C[i,s] / effort[i,s] #C[i,s]/f[i,s]
      }                               
    }
  }
  sim <- list(years, strata, r, K, pi, q1, q2, q3, q_err, d_err, b_err, s_err, 
              f, df, pesd, rho, pe, depletion, retention, 
              bio, bio_est, prop, cpue, surv, 
              C, By, LBy, D, d_est, pte_m, pte_p)
  names(sim)<-c("years", "strata", "r", "K", "pi", "q1", "q2", "q3",
                "q_err", "d_err","b_err", "s_err", 
                "f", "df", "pesd", "rho", "pe", "depletion","retention",
                "bio","bio_est","prop", "cpue", "surv", 
                "C","By","LBy","D","d_est", "pte_m","pte_p")
  return(sim)
}

# plot your sims...

plot_sim <- function(sim) {
  library("wesanderson")
  names(wes_palettes)
  pal <- wes_palette("Darjeeling1", strata, type = "discrete")
  
  par(mfrow=c(3,3))
  plot(C[,1],ylim=c(0,max(C)), type="l", cex=1, col = pal[1])
  for (s in 2:strata){ #s<-4
    lines(C[,s],type="l", cex=0.5, col=pal[s])
  }
  
  plot(D[,1],ylim=c(0,max(By)), type="l", cex=1, col = pal[1])
  lines(By[,1], type="l", cex=1, col = pal[1], lty = 2)
  for (s in 2:strata){ #s<-4
    lines(D[,s],type="l", cex=0.5, col=pal[s])
    lines(By[,s],type="l", cex=0.5, col=pal[s], lty=2)
  }
  
  plot(cpue[,1],ylim=c(0,max(cpue, na.rm=T)), type="l", cex=1, col = pal[1])
  for (s in 2:strata){
    lines(cpue[,s],type="l", cex=0.5, col=pal[s])
  }
  
  plot(bio[,1],ylim=c(0,max(bio)), type="l", cex=1, col = pal[1])
  points(bio_est[,1], type="p", cex=1, col = pal[1], pch=18)
  for (s in 2:strata){
    lines(bio[,s],type="l", cex=0.5, col=pal[s])
    points(bio_est[,s],type="p", cex=1, col=pal[s], pch=18)
  }
  
  plot(surv[,1],ylim=c(0,max(surv)), type="l", cex=1, col = pal[1])
  points(surv[,1],pch=18, cex=1, col = pal[1])
  for (s in 2:strata){
    lines(surv[,s],type="l", cex=0.5, col=pal[s])
    points(surv[,s],type="p", cex=1, col=pal[s], pch=18)
  }
  
  plot(prop[,1],ylim=c(0,max(prop)), type="l", cex=1, col = pal[1])
  for (s in 2:strata){
    lines(prop[,s],type="l", cex=0.5, col=pal[s])
    #points(surv[,s],type="p", cex=1, col=pal[s], pch=18)
  }
  
  plot(pe); abline(h=0, col="red") 
}

#-------------------------------------------------------------------------------
# 1st step of simulation includes fake data for all years... which isn't likely
# So lets simulate what data is actually available...

simulate2 <- function(sim = sim, 
                      bio_srv_frq = 4, bio_srv_start = 10, #staggered surveys... only one strata sampled per year... 
                      cpue_frq = 1, cpue_start = 14, #not staggered... cpue index all strata for each year that is happens... 
                      srv_frq = 1, srv_start = 20,
                      c_cv = 0.05) {
  # simulate biomass estimates:
  bio_est2 <- data.frame()
  for (i in 1:sim$years){
    bio_est2[i,c(1:sim$strata)] <- NA
  }
  for (i in 1:sim$strata){
    bio_est2[seq(bio_srv_start+i,sim$years,bio_srv_frq),i] <- 
      sim$bio_est[seq(bio_srv_start+i,sim$years,bio_srv_frq),i]
  }
  
  B_cv2 <- ifelse(is.na(bio_est2),NA,min(sim$b_err)) 
  
  cpue2 <- data.frame()
  for (i in 1:sim$years){
    cpue2[i,c(1:sim$strata)] <- NA
  }
  for (i in 1:sim$strata){
    cpue2[seq(cpue_start,sim$years,cpue_frq),i] <- 
      sim$cpue[seq(cpue_start,sim$years,cpue_frq),i]
  }
  
  I1_cv2 <- ifelse(is.na(cpue2),NA,min(sim$q_err))
  
  srv2 <- data.frame()
  for (i in 1:sim$years){
    srv2[i,c(1:sim$strata)] <- NA
  }
  for (i in 1:sim$strata){
    srv2[seq(srv_start,sim$years,srv_frq),i] <- 
      sim$surv[seq(srv_start,sim$years,srv_frq),i]
  }
  
  I2_cv2 <- ifelse(is.na(srv2),NA,min(sim$s_err))
  
  C_obs <- sim$C + 0.000001
  C_cv <- matrix(nrow=sim$years, ncol = sim$strata)
  for (i in 1:sim$strata) {
    C_cv[,i] <- rep(c_cv,sim$years)
   }
  
  D_obs <- sim$D + 0.000001
  D_cv <- matrix(nrow=sim$years, ncol = sim$strata)
  for (i in 1:sim$strata) {
    D_cv[,i] <- sim$d_err
  }
  
  pal <- wes_palette("Darjeeling1", strata, type = "discrete")
  
  par(mfrow=c(3,3))
  plot(sim$C[,1],ylim=c(0,max(C)), type="l", cex=1, col = pal[1])
  for (s in 2:sim$strata){ #s<-4
    lines(sim$C[,s],type="l", cex=0.5, col=pal[s])
  }
  
  plot(sim$D[,1],ylim=c(0,max(C)), type="l", cex=1, col = pal[1])
  for (s in 2:sim$strata){ #s<-4
    lines(sim$D[,s],type="l", cex=0.5, col=pal[s])
  }
  
  plot(cpue2[,1],ylim=c(0,max(cpue, na.rm=T)), type="l", cex=1, col = pal[1])
  for (s in 2:sim$strata){
    lines(cpue2[,s],type="l", cex=0.5, col=pal[s])
  }
  
  plot(sim$bio[,1],ylim=c(0,max(bio)), type="l", cex=1, col = pal[1])
  points(bio_est2[,1], type="p", cex=1, col = pal[1], pch=18)
  for (s in 2:sim$strata){
    lines(sim$bio[,s],type="l", cex=0.5, col=pal[s])
    points(bio_est2[,s],type="p", cex=1, col=pal[s], pch=18)
  }
  
  plot(sim$surv[,1],ylim=c(0,max(surv)), type="l", cex=1, col = pal[1])
  points(srv2[,1],pch=18, cex=1, col = pal[1])
  for (s in 2:sim$strata){
    lines(sim$surv[,s],type="l", cex=0.5, col=pal[s])
    points(srv2[,s],type="p", cex=1, col=pal[s], pch=18)
  }
  
  plot(sim$prop[,1],ylim=c(0,max(prop)), type="l", cex=1, col = pal[1])
  for (s in 2:sim$strata){
    lines(sim$prop[,s],type="l", cex=0.5, col=pal[s])
    #points(surv[,s],type="p", cex=1, col=pal[s], pch=18)
  }
  
  plot(sim$pe); abline(h=0, col="red") 
  
  sim_dat <- list(years = sim$years,
                  strata = sim$strata,
                  C_obs=C_obs,C_cv=C_cv, 
                  D_obs=D_obs,  D_cv=D_cv,
                  #observed biomass
                  B_ests = bio_est2, B_cv = B_cv2,
                  # Index 1
                  I1_ests = cpue2, I1_cv = I1_cv2,
                  #index 2
                  I2_ests=srv2 ,I2_cv = I2_cv2 #,
                  )
  return(sim_dat)
}

#------------------------------------------------------------------------------
# Function to bundle data, model switches and changable priors:
prep_data <- function(data = sim_dat, # set up to accept data from simulation set up
                      bio_switch = 1, #1 for biomass data, 0 for no biomass data
                      ind1_switch = 1, #1 for cpue index data, 0 for no cpue data
                      ind2_switch = 0, #1 for 2nd cpue index, 0 for none
                      bio_xV_lastY = 1, 
                      cpue_xV_lastY = 1,
                      srv_xV_lastY = 1,
                      parameterization = "NC", #"C" for centered parameterization, "NC" for non-centered model
                      pmsy = 0.46,
                      pe_switch = 0, #0 = no process error, 1 = est. process error
                      r_prior_switch = 0, #0=beta, 1=gamma, 2=lognormal
                      rp1 = 1, #1st parameter for r prior
                      rp2 = 1,# 2nd parameter for r prior
                      k1 = 10, #mean for lognormal K parameter
                      k2 = 2, #sd for lognormal K parameter
                      pe_bound = -4,
                      P1_mu = log(1),  #make a switch for this... 
                      P1_sig = 0.01,
                      isig2_mu = 5,  #prior mu for fit of P
                      isig2_sig = 0.102  #prior sig for fit of P
) {
  N = data$years
  S = data$strata
  
  if (bio_switch == 0){
    #fake data to be ignored
    B_obs <- rep(1,S)
    B_cv <- rep(1,S)
    N_Bobs <- S
    B_pos <- rep(1,S)
    S_Bobs <- rep(1,S)
    bio_tau_vec<-rep(0,length(B_pos))  
    bio_tau_no<-length(unique(bio_tau_vec[bio_tau_vec != 0]))
  } else {
    B_obs <- as.vector(data$B_ests[!is.na(data$B_ests)])
    B_cv <- as.vector(data$B_cv[!is.na(data$B_cv)])
    N_Bobs <- sum (!is.na(data$B_ests))
    B_pos <- list()
    for (i in 1:S) {
      B_pos[[i]] <- which(!is.na(data$B_ests[,i]))
    }
    B_pos <- as.vector(unlist(B_pos))
    S_Bobs <- vector()
    for (i in 1:S) {
      S_Bobs[i] <- length(data$B_ests[,i][!is.na(data$B_ests[,i])])
    }
    bio_tau_vec<-ifelse(B_pos < bio_xV_lastY,0,0)  
    bio_tau_no<-length(unique(bio_tau_vec[bio_tau_vec != 0]))
  }
  
  if (ind1_switch == 0){
    #fake data to be ignored
    I1_obs <- rep(1,S)
    I1_cv <- rep(1,S)
    N_I1obs <- S
    I1_pos <- rep(1,S)
    S_I1obs <- rep(1,S)
    i1_tau_vec<-rep(0,length(I1_pos))  
    i1_tau_no<-length(unique(i1_tau_vec[i1_tau_vec != 0]))
  } else {
    I1_obs <- as.vector(data$I1_ests[!is.na(data$I1_ests)])
    I1_cv <- as.vector(data$I1_cv[!is.na(data$I1_cv)])
    N_I1obs <- sum (!is.na(data$I1_ests))
    I1_pos <- list()
    for (i in 1:strata) {
      I1_pos[[i]] <- which(!is.na(data$I1_ests[,i]))
    }
    I1_pos <- as.vector(unlist(I1_pos))
    S_I1obs <- vector()
    for (i in 1:strata) {
      S_I1obs[i] <- length(data$I1_ests[,i][!is.na(data$I1_ests[,i])])
    }
    i1_tau_vec<-ifelse(I1_pos < cpue_xV_lastY,0,0)  
    i1_tau_no<-length(unique(i1_tau_vec[i1_tau_vec != 0]))
  }
  
  # do we have an index 2?
  if (ind2_switch == 0){
    #fake data to be ignored
    I2_obs <- rep(1,S)
    I2_cv <- rep(1,S)
    N_I2obs <- S
    I2_pos <- rep(1,S)
    S_I2obs <- rep(1,S)
    i2_tau_vec<-rep(0,length(I2_pos))  
    i2_tau_no<-length(unique(i2_tau_vec[i2_tau_vec != 0]))
  } else {
    I2_obs <- as.vector(data$I2_ests[!is.na(data$I2_ests)])
    I2_cv <- as.vector(data$I2_cv[!is.na(data$I2_cv)])
    N_I2obs <- sum (!is.na(data$I2_ests))
    I2_pos <- list()
    for (i in 1:strata) {
      I2_pos[[i]] <- which(!is.na(data$I2_ests[,i]))
    }
    I2_pos <- as.vector(unlist(I2_pos))
    S_I2obs <- vector()
    for (i in 1:strata) {
      S_I2obs[i] <- length(data$I2_ests[,i][!is.na(data$I2_ests[,i])])
    }
    i2_tau_vec<-ifelse(I2_pos < srv_xV_lastY,0,0)  
    i2_tau_no<-length(unique(i2_tau_vec[i2_tau_vec != 0]))
  }
  
  # which parameterization do we want to use
  pt <- pella_toml(pmsy)
  msy_frac <- pmsy
  
  pte_p <- pt$pte_p
  pte_m <- pt$pte_m
  
  if (parameterization == "C") {
    pte <- pte_p
  } else if (parameterization == "NC") {
    pte <- pte_m
  }
  
  stan_dat <- list(N = data$years,
                   S = data$strata,
                   #catch data
                   C_obs = data$C_obs,  
                   C_cv = data$C_cv,
                   #discard data
                   D_obs = data$D_obs,  #can't have catch of 0...
                   D_cv = data$D_cv,
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
                   rp1 = rp1, #1st parameter for r prior
                   rp2 = rp2,# 2nd parameter for r prior
                   k1 = k1, #mean for lognormal K parameter
                   k2 = k2, #sd for lognormal K parameter
                   alpha = rep(1,S), #c(1,1,1,1),
                   pe_bound = pe_bound,
                   P1_mu = P1_mu,  #make a switch for this... 
                   P1_sig = P1_sig,
                   isig2_mu = isig2_mu,  #prior mu for fit of P
                   isig2_sig = isig2_sig,  #prior sig for fit of P
                   #model switches
                   r_prior_switch = r_prior_switch, #0=beta, 1=gamma, 2=lognormal
                   pe_switch = pe_switch, #0 = no process error, 1 means random pe, 2 means autocorrelated pe (in development)
                   bio_switch = bio_switch,
                   ind1_switch = ind1_switch,
                   ind2_switch = ind2_switch
  )
  
  return(stan_dat)
}

#------------------------------------------------------------------------------
# get some initial P values for model running:

Pinit_guess <- function(bio_dat = bio_est2, strata, P1guess, bio_switch = 1, cpue2 = cpue2) {
  P1guess <- P1guess
  Pinit <-matrix(nrow=years,ncol=strata)
  Pinit[1,] <- rep(P1guess,strata)
  if (bio_switch == 1) {
    for (i in 1:strata) { #i <- 1
      
      bioests <- bio_dat[,i] #[!is.na(bio_est2[,i])]
      samps <- length(bioests)
      pests <- bioests / (max(bioests, na.rm=T) * 1/P1guess) # - (1 - P1guess)
      pests_nonna <- pests[!is.na(pests)]
    
      step1 <- (min(pests, na.rm=T)-max(pests, na.rm=T))/
        (which(pests == min(pests, na.rm=T))-1)
      step2 <- (pests_nonna[length(pests_nonna)]-min(pests, na.rm=T))/
        ( years-1-which(pests == min(pests, na.rm=T)))
      
      Pinit[,i] <- c(seq(max(pests, na.rm=T),min(pests, na.rm=T),
                         by=step1),
                     seq(min(pests, na.rm=T),pests_nonna[length(pests_nonna)],
                         by=step2))
      
    }
  } else {
    for (s in 1:strata) { #s <- 1
      low <- min(cpue2[,s], na.rm=T) / max(cpue2[,s] - P1guess, na.rm=T)
      by <- (P1guess - low) / (years - 1)
      Pinit[,s] <- seq(P1guess,
                       low,
                       -by)
    }
  }
    
  library("wesanderson")
  #names(wes_palettes)
  pal <- wes_palette("Darjeeling1", strata, type = "discrete")
    plot(Pinit[,1], col=pal[1], ylim=c(0,1))
    
    for (i in 2:strata) {
      points(Pinit[,i], col=pal[i])
    }
    
    for (i in 1:strata) {
      lines(prop[,i], col=pal[i])
      points(bio_dat[,i]/max(bio_dat[,i], na.rm=T)*P1guess, col=pal[i], pch=18)
    }
    
    return(Pinit)
}

#Pinit <- Pinit_guess(bio_dat = bio_est2, strata, P1guess, bio_switch = 1, cpue2 = cpue2)
#--------------------------------------------------------------------------------
# get initial pi estimates for running model: 

pi_init_guess <- function(maxbio = 0, stan_dat = stan_dat, raw_dat = raw_dat) {
  pis <- vector()
  if (stan_dat$bio_switch == 1) {
    for (s in 1:strata) { #s<-1
      maxbio <- maxbio + max(raw_dat$B_ests[,s], na.rm=T)
    }
    for (s in 1:strata) {
      if(s == strata) {
        pis[s] <- 1 - sum(pis)
      } else {
        pis[s] <- max(raw_dat$B_ests[,s], na.rm=T)/maxbio
      }
    }
  } else {
    for (s in 1:strata) { #s<-1
      maxbio <- maxbio + max(raw_dat$I1_ests[,s], na.rm=T)
    }
    for (s in 1:strata) {
      if(s == strata) {
        pis[s] <- 1 - sum(pis)
      } else {
        pis[s] <- max(raw_dat$I1_ests[,s], na.rm=T)/maxbio
      }
    }
  }
  return(pis)
}

#pis<-pi_init_guess(maxbio=0, bio_dat = bio_est2)

#-------------------------------------------------------------------------------
# get initial iq estimatess:

iq_init_guess<-function(stan_dat = stan_dat, raw_dat = raw_dat, bioguess=bioguess){
  bioguess <- bioguess #guess a rough estimate of biomass to get started.. 
  iq1 <- vector()
  iq2 <- vector()
  if (stan_dat$bio_switch == 1 & stan_dat$ind1_switch == 1 & stan_dat$ind2_switch == 1) { 
    for (s in 1:stan_dat$S) {
      iq1[s] <- 1 / (mean(raw_dat$I1_ests[,s], na.rm=T) / mean(raw_dat$B_ests[,s], na.rm=T))
      iq2[s] <- 1 / (mean(raw_dat$I2_ests[,s], na.rm=T) / mean(raw_dat$B_ests[,s], na.rm=T))
    }
  } 
  if (stan_dat$bio_switch == 0 & stan_dat$ind1_switch == 1 & stan_dat$ind2_switch == 1) { 
    for (s in 1:stan_dat$S) {
      iq1[s] <- 1 / (mean(raw_dat$I1_ests[,s], na.rm=T) / bioguess[s])
      iq2[s] <- 1 / (mean(raw_dat$I2_ests[,s], na.rm=T) / bioguess[s])
    }
  }
  if (stan_dat$ind1_switch == 1 & stan_dat$ind2_switch == 0) {
    for (s in 1:stan_dat$S) {
      iq1[s] <- 1 / (mean(raw_dat$I1_ests[,s], na.rm=T) / bioguess[s])
      iq2[s] <- 200
    }
  }
  if (stan_dat$ind1_switch == 0) {
    iq1[s] <- 200
    iq2[s] <- 200
  }
  return(list(iq1,iq2))
}

#iqs <- iq_init_guess(bio_dat = bio_est2,cpue = cpue2,surv = surv2)

#-------------------------------------------------------------------------------
# Get F initial values...
# only available with biomass estimates right now... 
Finit_guess <- function(stan_dat = stan_dat, raw_dat = raw_dat, bioguess=bioguess, iq1 = iq1) #,
                        #bio_dat = bio_est2, C_dat = C_obs, strata, bio_switch = 1, cpue2 = cpue2) 
{
  Finit <-matrix(nrow=stan_dat$N,ncol=stan_dat$S)
  if (bio_switch == 1) {
    for (i in 1:stan_dat$S) { #i <- 1
      bioests <- raw_dat$B_ests[,i] #[!is.na(bio_est2[,i])]
      Cs <- raw_dat$C_obs[,i]
      Fs <- Cs / bioests
      Fs[1] <- mean(Fs, na.rm=T)
      Fs[length(Fs)] <- median(Fs, na.rm=T)
      Fs <- na.approx(Fs)
      Finit[,i] <- Fs
    }
  } else {
    for (i in 1:stan_dat$S) {
      bioests <- raw_dat$I1_ests[,i] * iq1[i]#[!is.na(bio_est2[,i])]
      Cs <- raw_dat$C_obs[,i]
      Fs <- Cs / bioests
      Fs[1] <- mean(Fs, na.rm=T)
      Fs[length(Fs)] <- median(Fs, na.rm=T)
      Fs <- na.approx(Fs)
      Finit[,i] <- Fs
    }
  }
  
  #names(wes_palettes)
  pal <- wes_palette("Darjeeling1", strata, type = "discrete")
  plot(Finit[,1], col=pal[1], ylim=c(0,max(Finit)*1.1))
  
  for (i in 2:strata) {
    points(Finit[,i], col=pal[i])
  }
  
  for (i in 1:strata) {
    points(C_obs[,i] / bio_dat[,i], col=pal[i], pch=18)
 #   points(bio_dat[,i]/max(bio_dat[,i], na.rm=T)*P1guess, col=pal[i], pch=18)
  }
  
  return(Finit)
}

#------------------------------------------------------------------------------
# inital values for 1st time model run...
new_inits <- function(stan_dat = stan_dat,
                      raw_dat = data,
                      chains = 3,
                      r_guess = 0.04,
                      K_guess = 5000,
                      i1tau2 = i1tau2,
                      i2tau2 = i2tau2,
                      isigma2 = isigma2,
                      logvar = logvar,
                      medP1 = medP1,
                      fsigma2 = fsigma2) {
  P1guess <- exp(stan_dat$P1_mu)
  Pinit <- Pinit_guess(bio_dat = raw_dat$B_ests, strata=stan_dat$S, 
                       P1guess=P1guess, bio_switch = stan_dat$bio_switch, 
                       cpue2 = raw_dat$I1_ests)
  
  Ps <- list()
  for (i in 1:chains) {
    if (i == 1) {
      P1 = Pinit
      Ps[[i]] <- P1
    } else {
      P = Pinit * rnorm(stan_dat$N*stan_dat$S,1,0.1)
      assign(paste0("P",i),P)
      Ps[[i]] <- P
    }
    
  }
  
  muq <- stan_dat$P1_mu
  mus <- list()
  for (i in 1:chains) { # i <- 1
    for (s in 1:stan_dat$S) {
      if (s == 1) {
        mu_s <- rep(muq,stan_dat$N)
        mu <- mu_s
      } else {
        mu_s <- rep(muq,stan_dat$N) + rnorm(stan_dat$N,0,0.05)
        mu <- cbind(mu,mu_s)
      }
    }
    mu <- as.data.frame(mu)
    names(mu) <- NULL
    assign(paste0("mu",i),mu)
    mus[[i]] <- assign(paste0("mu",i),mu)
  }
  
  maxbio <- 0; 
  
  pi_s<-pi_init_guess(maxbio=0, stan_dat = stan_dat, raw_dat = raw_dat)
  pis <- list()
  for (i in 1:chains) {
    if (i == 1) {
      pig = pi_s
    } else {
      pig <- c(pi_s[c(1:(stan_dat$S-1))] * runif(stan_dat$S-1, 1, 1.1))
      pig <- c(pig, 1-sum(pig))
      #assign(paste0("pi",i),pig)
    }
    pis[[i]] <- assign(paste0("pi",i),pig)
  }
 
  bioguess <- vector()
  for (i in 1:stan_dat$S) {
    bioguess[i] <- median(raw_dat$C_obs[,i], na.rm=T) / (0.5 * r_guess)
  }
  
  iqs <- iq_init_guess(stan_dat = stan_dat, raw_dat = raw_dat, bioguess=bioguess)
  
  #Finit function only available with biomass data at the moment: 
  Finit <- Finit_guess(stan_dat = stan_dat, raw_dat = raw_dat, bioguess=bioguess, iq1 = iqs[[1]])
  
  for (i in 1:chains) {
    if (i == 1) {
      F1<-Finit
    } else {
      Fg <- Finit * rnorm(stan_dat$N * stan_dat$S,1,0.1) 
      assign(paste0("F",i),Fg)
    }
  }
  
  inits <- list()
  for (i in 1:chains) { #i <- 2
    one_init <- list(r = r_guess * runif(1,1,1.2),
                     K = K_guess * runif(1,1,1.2),
                     iq1=iqs[[1]] * runif(stan_dat$S,1,1.1),
                     iq2=iqs[[2]] * runif(stan_dat$S,1,1.1),
                     i1tau2 = i1tau2[[i]], #c(20,40,40,40),
                     i2tau2 = i2tau2[[i]], 
                     isigma2=isigma2[i],
                     logvar=logvar[i], 
                     eps = rep(0,years), 
                     fsigma2=fsigma2[i],
                     P=Ps[[i]], 
                     mu = mus[[i]],
                     medP1 = medP1[i],
                     C = stan_dat$C_obs, D = stan_dat$D_obs, #F = F1,
                     pi = pis[[i]])
    inits[[i]] <- assign(paste0("inits",i),one_init)
  }
  
  return(inits)
  
}

#-------------------------------------------------------------------------------
# Save model run values for initial values in future model runs:

save_inits<-function(fit = fit, parameterization = "NC", margQ = TRUE){
  if (parameterization == "C") {
    Ps<-data.frame()
    for (i in 1:years){
      for (j in 1:strata) {
        Ps[i,j] <- median(rstan::extract(fit)$P[,i,j])
      }
    }
    
    mus <- 1
  } else {
    mus<-data.frame()
    for (i in 1:years){
      for (j in 1:strata) {
        mus[i,j] <- median(rstan::extract(fit)$mu[,i,j])
      }
    }
    Ps <- 1
  }
  
  tune_inits <- list()
  for (j in 1:chains){ #j<-1
    pis<-vector()
    for (s in 1:strata) {
      if (s < strata) {
        pis[s]<-median(rstan::extract(fit)$pi[,s])*rnorm(1,1,0.1)
      } else {
        pis[s]<-1-sum(pis)
      }
    }
    
    #taus are for margQ models
    if (margQ == TRUE) {
      taus1<-vector()
      taus2<-vector()
      for (s in 1:strata) {
        taus1[s]<-median(rstan::extract(fit)$i1tau2[,s])*rnorm(1,1,0.1)
        taus2[s]<-median(rstan::extract(fit)$i1tau2[,s])*rnorm(1,1,0.1)
      }
      iq1s <- rep(1,strata) #c(1,1,1,1)
      iq2s <- rep(1,strata) #c(1,1,1,1)
    } else {
      taus1 <- rep(1,strata) #c(1,1,1,1)
      taus2 <- rep(1,strata) #c(1,1,1,1)
      iq1s <- vector()
      iq2s <- vector()
      for (s in 1:strata) {
        iq1s[s] <- median(rstan::extract(fit)$iq1[,s])*rnorm(1,1,0.1)
        iq2s[s] <- median(rstan::extract(fit)$iq2[,s])*rnorm(1,1,0.1)
      }
    }
    
    new_inits <- list(r=median(rstan::extract(fit)$r)*rnorm(1,1,0.2), 
                      K=median(rstan::extract(fit)$K)*rnorm(1,1,0.1), 
                      ## turn off iq1 and iq2 if using a margQ model
                      iq1=iq1s,#c(80,90,40,60), 
                      iq2=iq2s,
                      isigma2=median(rstan::extract(fit)$isigma2)*rnorm(1,1,0.1), 
                      #taus are for margQ models
                      i1tau2 = taus1,
                      i2tau2 = taus2,
                      logvar=median(rstan::extract(fit)$logvar)*rnorm(1,1,0.1), 
                      eps = rep(0,years),
                      #extra variance terms:
                      #tau_b = median(rstan::extract(tune)$tau_b)*rnorm(1,1,0.1), 
                      #tau_i1 = median(rstan::extract(tune)$tau_i1)*rnorm(1,1,0.1), 
                      #tau_i2 = median(rstan::extract(tune)$tau_i2)*rnorm(1,1,0.1), 
                      P = as.matrix(Ps*rnorm(years,1,0.1)),
                      mu = as.matrix(mus*rnorm(years,1,0.1)), 
                      C = C+0.000001, 
                      D = D_obs+0.000001,
                      pi = pis
    )
    if (parameterization == "C") {
      new_inits$P[1,] <- rep(0.99,strata) #c(0.99,0.99,0.99,0.99)
    }
    assign(paste0("tune_inits",j),new_inits)
    tune_inits[[j]] <- new_inits
    #new_inits$r
  }
  return(tune_inits)
}

#try <- save_inits(fit = fit, parameterization = "NC", margQ = TRUE)

#saveRDS(try,"inits.rds")

#------------------------------------------------------------------------------
# plot posterior functions:
plot_biomass <- function(fit=fit,years=years,strata=strata,
                         save=TRUE, 
                         bio_dat = bio_est2, bio_cv = B_cv2,
                         sim=TRUE, 
                         sim_bio = bio, mod_name="Testing",
                         units = "mt",
                         year1 = 1980) {
  samples <- rstan::extract(fit, permuted = TRUE)
  N<-years
  true_years = seq(year1,year1+years-1,1)
  
  quants4<-array(dim=c(3,years,strata))
  for (i in 1:years){
    for (s in 1:strata){
      quants4[,i,s] <- quantile(samples$Biomass[,i,s],probs=,c(0.05,0.5,0.95))
    }
  }
  quants44 <- data.frame()
  for (s in 1:strata) {
    ssamp <- t(quants4[,,s]) %>% data.frame() %>% mutate(strata = s,
                                                         year = seq(1,years,1))
    colnames(ssamp) <- c("lo90","median","hi90","strata","year")
    if (s == 1) {
      quants44 <- ssamp
    } else {
      quants44 <- rbind(quants44,ssamp)
    }
  }
  
  quants44$true_year <- rep(true_years,strata)
  #names(wes_palettes)
  pal <- wes_palette("Darjeeling1", 4, type = "discrete")
  
  str_ref <- strata
  bio_dat %>% data.frame() %>% 
    pivot_longer(cols = seq(1,str_ref,1), names_to = "strata", values_to = "est" ) %>%
    arrange(strata) %>% 
    mutate(strata = gsub("[^0-9.-]", "", strata),
           year = rep(seq(1,years,1),str_ref)) %>%
    full_join(bio_cv %>% data.frame() %>%
                pivot_longer(cols = seq(1,str_ref,1), names_to = "strata", values_to = "cv") %>%
                arrange(strata) %>% 
                mutate(strata = gsub("[^0-9.-]", "", strata),
                       year = rep(seq(1,years,1),str_ref)),
              by=c("strata","year")) %>% 
    mutate(lo90 = exp(log(est) - 1.96*(log(est + est * cv) - log(est))),
           hi90 = exp(log(est) + 1.96*(log(est + est * cv) - log(est)))) -> plot_data
  plot_data$true_year <- rep(true_years,strata)
  
  if (sim == TRUE) {
    sim_bio %>% data.frame() %>% 
    pivot_longer(cols = seq(1,str_ref,1), names_to = "strata", values_to = "est" ) %>%
    arrange(strata) %>% 
    mutate(strata = gsub("[^0-9.-]", "", strata),
           year = rep(seq(1,years,1),str_ref)) -> true_bio
    true_bio$true_year <- rep(true_years,strata)
    
    ggplot(quants44) +
      scale_fill_manual(values = pal) +
      geom_line(data = true_bio, aes(true_year,est), col=pal[3], size=1.25, linetype=1) +
      geom_ribbon(aes(ymin=lo90,ymax=hi90,x=true_year),col=NA,fill=pal[2], alpha=0.2) +
      geom_line(aes(true_year,median),col=pal[2]) +
      geom_point(aes(true_year,median),col=pal[2]) +
      facet_grid(rows = vars(strata)) +  #facet_wrap?
      #facet_wrap(~ strata) +
      geom_point(data=plot_data, aes(true_year,est), col=pal[4]) +
      geom_errorbar(data=plot_data,aes(x = true_year, ymin = lo90, ymax = hi90), col=pal[4])+ 
      labs(title = "Biomass") +
      ylab(paste0("Biomass (",units,")")) +
      xlab("Year") + 
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust = 0.5
        )) + 
      scale_x_continuous(breaks=seq(year1,years+year1-1,5))-> plot
    
  } else {
    ggplot(quants44) +
      scale_fill_manual(values = pal) +
      geom_ribbon(aes(ymin=lo90,ymax=hi90,x=true_year),col=NA,fill=pal[2], alpha=0.2) +
      geom_line(aes(true_year,median),col=pal[2]) +
      geom_point(aes(true_year,median),col=pal[2]) +
      facet_grid(rows = vars(strata)) +
      #facet_wrap(~ strata) +
      geom_point(data=plot_data, aes(true_year,est), col=pal[4]) +
      geom_errorbar(data=plot_data,aes(x = true_year, ymin = lo90, ymax = hi90), col=pal[4])+ 
      labs(title = "Biomass") +
      ylab(paste0("Biomass (",units,")")) +
      xlab("Year") + 
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust = 0.5
        )) + 
      scale_x_continuous(breaks=seq(year1,years+year1-1,5))-> plot
  }
  
  if (save == TRUE) {
    ggsave(paste0("Production_models_STAN/Figures/biomass_", mod_name, ".png"),plot=plot,
           dpi=300,  height=6, width=7, units="in")
  }
  print(plot)

}

#-------------------------------------------------------------------------------
plot_p <- function(fit=fit,years=years,strata=strata,
                         save=TRUE, 
                         bio_dat = bio_est2, bio_cv = B_cv2,
                         sim=TRUE, 
                         sim_bio = bio, mod_name="Testing",
                         units = "mt",
                         year1 = 1980) {
  samples <- rstan::extract(fit, permuted = TRUE)
  N<-years
  true_years = seq(year1,year1+years-1,1)
  
  quants4<-array(dim=c(3,years,strata))
  for (i in 1:years){
    for (s in 1:strata){
      quants4[,i,s] <- quantile(samples$P[,i,s],probs=,c(0.05,0.5,0.95))
    }
  }
  quants44 <- data.frame()
  for (s in 1:strata) {
    ssamp <- t(quants4[,,s]) %>% data.frame() %>% mutate(strata = s,
                                                         year = seq(1,years,1))
    colnames(ssamp) <- c("lo90","median","hi90","strata","year")
    if (s == 1) {
      quants44 <- ssamp
    } else {
      quants44 <- rbind(quants44,ssamp)
    }
  }
  
  quants44$true_year <- rep(true_years,strata)
  #names(wes_palettes)
  pal <- wes_palette("Darjeeling1", 4, type = "discrete")
  
  prop <- vector(); j <- 1
  for (i in 1:strata) {
    for (y in 1:years) {
      prop[j] <- pi[i]
      j <- j + 1
    }
  }
  
  str_ref <- strata
  bio_dat %>% data.frame() %>% 
    pivot_longer(cols = seq(1,str_ref,1), names_to = "strata", values_to = "est" ) %>%
    arrange(strata) %>% 
    mutate(strata = gsub("[^0-9.-]", "", strata),
           year = rep(seq(1,years,1),str_ref),
           prop = prop) %>%
    full_join(bio_cv %>% data.frame() %>%
                pivot_longer(cols = seq(1,str_ref,1), names_to = "strata", values_to = "cv") %>%
                arrange(strata) %>% 
                mutate(strata = gsub("[^0-9.-]", "", strata),
                       year = rep(seq(1,years,1),str_ref)),
              by=c("strata","year")) %>% 
    mutate(prop_K = est/(K * prop),
           lo90 = exp(log(prop_K) - 1.96*(log(prop_K + prop_K * cv) - log(prop_K))),
           hi90 = exp(log(prop_K) + 1.96*(log(prop_K + prop_K * cv) - log(prop_K)))) -> plot_data
  plot_data$true_year <- rep(true_years,strata)
  
  if (sim == TRUE) {
    sim_bio %>% data.frame() %>% 
      pivot_longer(cols = seq(1,str_ref,1), names_to = "strata", values_to = "est" ) %>%
      arrange(strata) %>% 
      mutate(strata = gsub("[^0-9.-]", "", strata),
             year = rep(seq(1,years,1),str_ref),
             prop = prop,
             prop_K = est/(K * prop)) -> true_p
    
    true_p$true_year <- rep(true_years,strata)
    
    ggplot(quants44) +
      scale_fill_manual(values = pal) +
      geom_line(data = true_p, aes(true_year,prop_K), col=pal[3], size=1.25, linetype=1) +
      geom_ribbon(aes(ymin=lo90,ymax=hi90,x=true_year),col=NA,fill=pal[2], alpha=0.2) +
      geom_line(aes(true_year,median),col=pal[2]) +
      geom_point(aes(true_year,median),col=pal[2]) +
      facet_grid(rows = vars(strata)) +  #facet_wrap?
      #facet_wrap(~ strata) +
      geom_point(data=plot_data, aes(true_year,prop_K), col=pal[4]) +
      geom_errorbar(data=plot_data,aes(x = true_year, ymin = lo90, ymax = hi90), col=pal[4])+ 
      labs(title = "Proportion of K") +
      ylab(paste0("Proportion of K")) +
      xlab("Year") + 
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust = 0.5
        )) + 
      scale_x_continuous(breaks=seq(year1,years+year1-1,5))-> plot
    
  } else {
    ggplot(quants44) +
      scale_fill_manual(values = pal) +
      geom_ribbon(aes(ymin=lo90,ymax=hi90,x=true_year),col=NA,fill=pal[2], alpha=0.2) +
      geom_line(aes(true_year,median),col=pal[2]) +
      geom_point(aes(true_year,median),col=pal[2]) +
      facet_grid(rows = vars(strata)) +
      #facet_wrap(~ strata) +
      geom_point(data=plot_data, aes(true_year,est), col=pal[4]) +
      geom_errorbar(data=plot_data,aes(x = true_year, ymin = lo90, ymax = hi90), col=pal[4])+ 
      labs(title = "Biomass") +
      ylab(paste0("Biomass (",units,")")) +
      xlab("Year") + 
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust = 0.5
        )) + 
      scale_x_continuous(breaks=seq(year1,years+year1-1,5))-> plot
  }
  
  if (save == TRUE) {
    ggsave(paste0("Production_models_STAN/Figures/propK_", mod_name, ".png"),plot=plot,
           dpi=300,  height=6, width=7, units="in")
  }
  print(plot)
  
}

#-------------------------------------------------------------------------------
plot_index1<- function(fit=fit,years=years,strata=strata,
                      save=TRUE, 
                      cpue_dat = cpue2, cpue_cv = I1_cv2,
                      sim=TRUE, 
                      sim_cpue = cpue, mod_name="Testing",
                      units = "mt",
                      year1 = 1980,
                      cpue_name = "Fishery CPUE",
                      cpue_type = "kgs/box") {
  samples <- rstan::extract(fit, permuted = TRUE)
  N<-years
  true_years = seq(year1,year1+years-1,1)
  
  quants4<-array(dim=c(3,years,strata))
  for (i in 1:years){
    for (s in 1:strata){
      quants4[,i,s] <- quantile(samples$I1est[,i,s],probs=,c(0.05,0.5,0.95))
    }
  }
  quants44 <- data.frame()
  for (s in 1:strata) {
    ssamp <- t(quants4[,,s]) %>% data.frame() %>% mutate(strata = s,
                                                         year = seq(1,years,1))
    colnames(ssamp) <- c("lo90","median","hi90","strata","year")
    if (s == 1) {
      quants44 <- ssamp
    } else {
      quants44 <- rbind(quants44,ssamp)
    }
  }
  
  quants44$true_year <- rep(true_years,strata)
  #names(wes_palettes)
  pal <- wes_palette("GrandBudapest1", 4, type = "discrete")
  
  str_ref <- strata
  cpue_dat %>% data.frame() %>% 
    pivot_longer(cols = seq(1,str_ref,1), names_to = "strata", values_to = "est" ) %>%
    arrange(strata) %>% 
    mutate(strata = gsub("[^0-9.-]", "", strata),
           year = rep(seq(1,years,1),str_ref)) %>%
    full_join(cpue_cv %>% data.frame() %>%
                pivot_longer(cols = seq(1,str_ref,1), names_to = "strata", values_to = "cv") %>%
                arrange(strata) %>% 
                mutate(strata = gsub("[^0-9.-]", "", strata),
                       year = rep(seq(1,years,1),str_ref)),
              by=c("strata","year")) %>% 
    mutate(lo90 = exp(log(est) - 1.96*(log(est + est * cv) - log(est))),
           hi90 = exp(log(est) + 1.96*(log(est + est * cv) - log(est)))) -> plot_data
  plot_data$true_year <- rep(true_years,strata)
  
  if (sim == TRUE) {
    sim_cpue %>% data.frame() %>% 
      pivot_longer(cols = seq(1,str_ref,1), names_to = "strata", values_to = "est" ) %>%
      arrange(strata) %>% 
      mutate(strata = gsub("[^0-9.-]", "", strata),
             year = rep(seq(1,years,1),str_ref)) -> true_cpue
    true_cpue$true_year <- rep(true_years,strata)
    
    ggplot(quants44) +
      scale_fill_manual(values = pal) +
      geom_line(data = true_cpue, aes(true_year,est), col=pal[3], size=0.75, linetype=2) +
      geom_ribbon(aes(ymin=lo90,ymax=hi90,x=true_year),col=NA,fill=pal[2], alpha=0.2) +
      geom_line(aes(true_year,median),col=pal[2]) +
      geom_point(aes(true_year,median),col=pal[2]) +
      facet_grid(rows = vars(strata)) +  #facet_wrap?
      #facet_wrap(~ strata) +
      geom_point(data=plot_data, aes(true_year,est), col=pal[3]) +
      geom_errorbar(data=plot_data,aes(x = true_year, ymin = lo90, ymax = hi90), col=pal[3])+ 
      labs(title = cpue_name) +
      ylab(paste0(cpue_name," (",cpue_type,")" )) +
      xlab("Year") + 
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust = 0.5
        )) + 
      scale_x_continuous(breaks=seq(year1,years+year1-1,5))-> plot
    
  } else {
    ggplot(quants44) +
      scale_fill_manual(values = pal) +
      geom_ribbon(aes(ymin=lo90,ymax=hi90,x=true_year),col=NA,fill=pal[2], alpha=0.2) +
      geom_line(aes(true_year,median),col=pal[2]) +
      geom_point(aes(true_year,median),col=pal[2]) +
      facet_grid(rows = vars(strata)) +  #facet_wrap?
      #facet_wrap(~ strata) +
      geom_point(data=plot_data, aes(true_year,est), col=pal[3]) +
      geom_errorbar(data=plot_data,aes(x = true_year, ymin = lo90, ymax = hi90), col=pal[3])+ 
      labs(title = cpue_name) +
      ylab(paste0(cpue_name," (",cpue_type,")" )) +
      xlab("Year") + 
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust = 0.5
        )) + 
      scale_x_continuous(breaks=seq(year1,years+year1-1,5))-> plot
  }
  
  if (save == TRUE) {
    ggsave(paste0("Production_models_STAN/Figures/index1_", mod_name, ".png"),plot=plot,
           dpi=300,  height=6, width=7, units="in")
  }
  print(plot)
  
}

#-------------------------------------------------------------------------------
plot_index2<- function(fit=fit,years=years,strata=strata,
                       save=TRUE, 
                       cpue_dat = surv2, cpue_cv = I2_cv2,
                       sim=TRUE, 
                       sim_cpue = cpue, mod_name="Testing",
                       units = "mt",
                       year1 = 1980,
                       cpue_name = "Fishery CPUE",
                       cpue_type = "kgs/box") {
  samples <- rstan::extract(fit, permuted = TRUE)
  N<-years
  true_years = seq(year1,year1+years-1,1)
  
  quants4<-array(dim=c(3,years,strata))
  for (i in 1:years){
    for (s in 1:strata){
      quants4[,i,s] <- quantile(samples$I2est[,i,s],probs=,c(0.05,0.5,0.95))
    }
  }
  quants44 <- data.frame()
  for (s in 1:strata) {
    ssamp <- t(quants4[,,s]) %>% data.frame() %>% mutate(strata = s,
                                                         year = seq(1,years,1))
    colnames(ssamp) <- c("lo90","median","hi90","strata","year")
    if (s == 1) {
      quants44 <- ssamp
    } else {
      quants44 <- rbind(quants44,ssamp)
    }
  }
  
  quants44$true_year <- rep(true_years,strata)
  #names(wes_palettes)
  pal <- wes_palette("GrandBudapest1", 4, type = "discrete")
  
  str_ref <- strata
  cpue_dat %>% data.frame() %>% 
    pivot_longer(cols = seq(1,str_ref,1), names_to = "strata", values_to = "est" ) %>%
    arrange(strata) %>% 
    mutate(strata = gsub("[^0-9.-]", "", strata),
           year = rep(seq(1,years,1),str_ref)) %>%
    full_join(cpue_cv %>% data.frame() %>%
                pivot_longer(cols = seq(1,str_ref,1), names_to = "strata", values_to = "cv") %>%
                arrange(strata) %>% 
                mutate(strata = gsub("[^0-9.-]", "", strata),
                       year = rep(seq(1,years,1),str_ref)),
              by=c("strata","year")) %>% 
    mutate(lo90 = exp(log(est) - 1.96*(log(est + est * cv) - log(est))),
           hi90 = exp(log(est) + 1.96*(log(est + est * cv) - log(est)))) -> plot_data
  plot_data$true_year <- rep(true_years,strata)
  
  if (sim == TRUE) {
    sim_cpue %>% data.frame() %>% 
      pivot_longer(cols = seq(1,str_ref,1), names_to = "strata", values_to = "est" ) %>%
      arrange(strata) %>% 
      mutate(strata = gsub("[^0-9.-]", "", strata),
             year = rep(seq(1,years,1),str_ref)) -> true_cpue
    true_cpue$true_year <- rep(true_years,strata)
    
    ggplot(quants44) +
      scale_fill_manual(values = pal) +
      geom_line(data = true_cpue, aes(true_year,est), col=pal[3], size=0.75, linetype=2) +
      geom_ribbon(aes(ymin=lo90,ymax=hi90,x=true_year),col=NA,fill=pal[2], alpha=0.2) +
      geom_line(aes(true_year,median),col=pal[2]) +
      geom_point(aes(true_year,median),col=pal[2]) +
      facet_grid(rows = vars(strata)) +  #facet_wrap?
      #facet_wrap(~ strata) +
      geom_point(data=plot_data, aes(true_year,est), col=pal[3]) +
      geom_errorbar(data=plot_data,aes(x = true_year, ymin = lo90, ymax = hi90), col=pal[3])+ 
      labs(title = cpue_name) +
      ylab(paste0(cpue_name," (",cpue_type,")" )) +
      xlab("Year") + 
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust = 0.5
        )) + 
      scale_x_continuous(breaks=seq(year1,years+year1-1,5))-> plot
    
  } else {
    ggplot(quants44) +
      scale_fill_manual(values = pal) +
      geom_ribbon(aes(ymin=lo90,ymax=hi90,x=true_year),col=NA,fill=pal[2], alpha=0.2) +
      geom_line(aes(true_year,median),col=pal[2]) +
      geom_point(aes(true_year,median),col=pal[2]) +
      facet_grid(rows = vars(strata)) +  #facet_wrap?
      #facet_wrap(~ strata) +
      geom_point(data=plot_data, aes(true_year,est), col=pal[3]) +
      geom_errorbar(data=plot_data,aes(x = true_year, ymin = lo90, ymax = hi90), col=pal[3])+ 
      labs(title = cpue_name) +
      ylab(paste0(cpue_name," (",cpue_type,")" )) +
      xlab("Year") + 
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust = 0.5
        )) + 
      scale_x_continuous(breaks=seq(year1,years+year1-1,5))-> plot
  }
  
  if (save == TRUE) {
    ggsave(paste0("Production_models_STAN/Figures/index1_", mod_name, ".png"),plot=plot,
           dpi=300,  height=6, width=7, units="in")
  }
  print(plot)
  
}

#-------------------------------------------------------------------------------
# plot catch estimates
plot_catch <- function(fit=fit,years=years,strata=strata,
                         save=TRUE, 
                         dat = fish_dat2, #bio_cv = B_cv2,
                         sim=TRUE, 
                         mod_name="Testing",
                         units = "mt",
                         year1 = 1980) {
  samples <- rstan::extract(fit, permuted = TRUE)
  N<-years
  true_years = seq(year1,year1+years-1,1)
  
  quants4<-array(dim=c(3,years,strata))
  for (i in 1:years){
    for (s in 1:strata){
      quants4[,i,s] <- quantile(samples$C[,i,s],probs=,c(0.05,0.5,0.95))
    }
  }
  quants44 <- data.frame()
  for (s in 1:strata) {
    ssamp <- t(quants4[,,s]) %>% data.frame() %>% mutate(strata = s,
                                                         year = seq(1,years,1))
    colnames(ssamp) <- c("lo90","median","hi90","strata","year")
    if (s == 1) {
      quants44 <- ssamp
    } else {
      quants44 <- rbind(quants44,ssamp)
    }
  }
  
  quants44$true_year <- rep(true_years,strata)
  #names(wes_palettes)
  pal <- wes_palette("FantasticFox1", 4, type = "discrete")
  
  str_ref <- strata
  dat$C_obs %>% data.frame() %>% 
    pivot_longer(cols = seq(1,str_ref,1), names_to = "strata", values_to = "est" ) %>%
    arrange(strata) %>% 
    mutate(strata = gsub("[^0-9.-]", "", strata),
           year = rep(seq(1,years,1),str_ref)) %>%
    full_join(dat$C_cv %>% data.frame() %>%
                pivot_longer(cols = seq(1,str_ref,1), names_to = "strata", values_to = "cv") %>%
                arrange(strata) %>% 
                mutate(strata = gsub("[^0-9.-]", "", strata),
                       year = rep(seq(1,years,1),str_ref)),
              by=c("strata","year")) %>% 
    mutate(lo90 = exp(log(est) - 1.96*(log(est + est * cv) - log(est))),
           hi90 = exp(log(est) + 1.96*(log(est + est * cv) - log(est)))) -> plot_data
  plot_data$true_year <- rep(true_years,strata)
  
  ggplot(quants44) +
      scale_fill_manual(values = pal) +
      #geom_line(data = true_bio, aes(true_year,est), col=pal[3], size=1.25, linetype=1) +
      geom_ribbon(aes(ymin=lo90,ymax=hi90,x=true_year),col=NA,fill=pal[2], alpha=0.2) +
      geom_line(aes(true_year,median),col=pal[2]) +
      geom_point(aes(true_year,median),col=pal[2]) +
      facet_grid(rows = vars(strata)) +  #facet_wrap?
      #facet_wrap(~ strata) +
      geom_point(data=plot_data, aes(true_year,est), col=pal[4]) +
      geom_errorbar(data=plot_data,aes(x = true_year, ymin = lo90, ymax = hi90), col=pal[4])+ 
      labs(title = "Catch") +
      ylab(paste0("Catch (",units,")")) +
      xlab("Year") + 
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust = 0.5
        )) + 
      scale_x_continuous(breaks=seq(year1,years+year1-1,5))-> plot
    
  if (save == TRUE) {
    ggsave(paste0("Production_models_STAN/Figures/catch_", mod_name, ".png"),plot=plot,
           dpi=300,  height=6, width=7, units="in")
  }
  print(plot)
  
}

#-------------------------------------------------------------------------------
plot_discards <- function(fit=fit,years=years,strata=strata,
                       save=TRUE, 
                       dat = fish_dat2, #bio_cv = B_cv2,
                       sim=TRUE, 
                       mod_name="Testing",
                       units = "mt",
                       year1 = 1980) {
  samples <- rstan::extract(fit, permuted = TRUE)
  N<-years
  true_years = seq(year1,year1+years-1,1)
  
  quants4<-array(dim=c(3,years,strata))
  for (i in 1:years){
    for (s in 1:strata){
      quants4[,i,s] <- quantile(samples$D[,i,s],probs=,c(0.05,0.5,0.95))
    }
  }
  quants44 <- data.frame()
  for (s in 1:strata) {
    ssamp <- t(quants4[,,s]) %>% data.frame() %>% mutate(strata = s,
                                                         year = seq(1,years,1))
    colnames(ssamp) <- c("lo90","median","hi90","strata","year")
    if (s == 1) {
      quants44 <- ssamp
    } else {
      quants44 <- rbind(quants44,ssamp)
    }
  }
  
  quants44$true_year <- rep(true_years,strata)
  #names(wes_palettes)
  pal <- wes_palette("FantasticFox1", 4, type = "discrete")
  
  str_ref <- strata
  dat$D_obs %>% data.frame() %>% 
    pivot_longer(cols = seq(1,str_ref,1), names_to = "strata", values_to = "est" ) %>%
    arrange(strata) %>% 
    mutate(strata = gsub("[^0-9.-]", "", strata),
           year = rep(seq(1,years,1),str_ref)) %>%
    full_join(dat$D_cv %>% data.frame() %>%
                pivot_longer(cols = seq(1,str_ref,1), names_to = "strata", values_to = "cv") %>%
                arrange(strata) %>% 
                mutate(strata = gsub("[^0-9.-]", "", strata),
                       year = rep(seq(1,years,1),str_ref)),
              by=c("strata","year")) %>% 
    mutate(lo90 = exp(log(est) - 1.96*(log(est + est * cv) - log(est))),
           hi90 = exp(log(est) + 1.96*(log(est + est * cv) - log(est)))) -> plot_data
  plot_data$true_year <- rep(true_years,strata)
  
  ggplot(quants44) +
    scale_fill_manual(values = pal) +
    #geom_line(data = true_bio, aes(true_year,est), col=pal[3], size=1.25, linetype=1) +
    geom_ribbon(aes(ymin=lo90,ymax=hi90,x=true_year),col=NA,fill=pal[2], alpha=0.2) +
    geom_line(aes(true_year,median),col=pal[2]) +
    geom_point(aes(true_year,median),col=pal[2]) +
    facet_grid(rows = vars(strata)) +  #facet_wrap?
    #facet_wrap(~ strata) +
    geom_point(data=plot_data, aes(true_year,est), col=pal[4]) +
    geom_errorbar(data=plot_data,aes(x = true_year, ymin = lo90, ymax = hi90), col=pal[4])+ 
    labs(title = "Discards") +
    ylab(paste0("Discards (",units,")")) +
    xlab("Year") + 
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 0.5
      )) + 
    scale_x_continuous(breaks=seq(year1,years+year1-1,5))-> plot
  
  if (save == TRUE) {
    ggsave(paste0("Production_models_STAN/Figures/discards_", mod_name, ".png"),plot=plot,
           dpi=300,  height=6, width=7, units="in")
  }
  print(plot)
  
}

#------------------------------------------------------------------------------
# process error
plot_pe <- function(fit = fit, sim_pe = pe, sim=TRUE, save=TRUE, 
                    mod_name="Testing", year1 = 1980) {
  samples <- rstan::extract(fit, permuted = TRUE)
  N<-years
  true_years = seq(year1,year1+years-1,1)
  
  quants4<-array(dim=c(3,years))
  for (i in 1:years){
    #for (s in 1:strata){
    quants4[,i] <- quantile(samples$eps[,i],probs=,c(0.05,0.5,0.95))
    #}
  }
  quants4<-t(quants4)
  quants4<-cbind(quants4,seq(1,years))
  colnames(quants4) <- c("lo90","median","hi90","year")
  
  pal <- wes_palette("FantasticFox1", 4, type = "discrete")
  
  if (sim == TRUE){
    true_pe<-cbind(seq(1,years),as.data.frame(sim_pe))
    colnames(true_pe)<-c("year","true_pe")
    ggplot(quants4 %>% data.frame()) +
      scale_fill_manual(values = pal) +
      geom_ribbon(aes(ymin=lo90,ymax=hi90,x=year),col=NA,fill=pal[1], alpha=0.2) +
      geom_line(aes(year,median),col=pal[1]) +
      geom_point(aes(year,median),col=pal[1]) +
      geom_point(data=true_pe, aes(year,true_pe), col=pal[3]) +
      labs(title = "process error") -> plot
  } else {
    ggplot(quants4 %>% data.frame()) +
      scale_fill_manual(values = pal) +
      geom_ribbon(aes(ymin=lo90,ymax=hi90,x=year),col=NA,fill=pal[1], alpha=0.2) +
      geom_line(aes(year,median),col=pal[1]) +
      geom_point(aes(year,median),col=pal[1]) +
      labs(title = "process error") -> plot
  }
  if (save == TRUE) {
    ggsave(paste0("Production_models_STAN/Figures/proc_err_", mod_name, ".png"),plot=plot,
           dpi=300,  height=6, width=7, units="in")
  }
  print(plot)

}

#-------------------------------------------------------------------------------
# Function to bundle and record 1 simulation run. Assumes 

