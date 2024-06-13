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
rs = runif(Narea, 0.04, 0.06)    # variability in rs
Ks = runif(Narea, 100000, 200000)
Nyear = 40
Year_start = 1   # when we have survey/CPUE info from 
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
    
  ## Catch is dependent of stock size
    if (Catch_type == "F")  H = sapply(1:Narea, function(x) rnorm(Nyear, c(seq(0.5, Hmax, length.out=Nyear-10), seq(Hmax, 0.5, length.out=10))*HMSYs, 0.1*HMSYs))

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

  IAs = B
  for (yr in 1: (Nyear)){
    for (sp in 1:Narea){
      IAs[yr,sp] <- qs[sp]*B[yr,sp]*exp(epsilonO[yr,sp])
    }
  }

  matplot(1:Nyear, IAs, type="l")

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
  
  par(mfrow=c(4,1), mar=c(2,4,1,1))
  matplot(Btrue, type="l")
  matplot(epsilon[Year_start:Nyear,], type="l")
  matplot(IA_obs, type="l")
  matplot(Catch[Year_start:(Nyear-1),], type="l")

  
  ## Does it exist an r_all that fits the above?
  # 
  # func <- function(x) {
  #   SS = sum((log(B_all[-1]) - log(pmax(B_all[-length(B_all)] + x*B_all[-length(B_all)]/p*(1 - (B_all[-length(B_all)]/K_all)^p) - Catch_all[-length(B_all)], 1)))^2)
  #   return(SS)
  # }
  # 
  # out <- optimize(func, lower = 0, upper = 5, maximum = FALSE)
  # 
  # r_est$r_est[i] <- out$minimum
  # r_est$pass[i] <- !any(B < 100)
  # #}
  # 
  # table(r_est$pass)
  # 
  # fake <- data.frame(true_mean = rowMeans(r_OM), true_median = apply(r_OM, 1, median), est =  r_est$r_est, col=apply(r_OM, 1, sd))
  # ggplot(fake, aes(x= true_mean, y = est, col=col)) + geom_point() + scale_color_viridis_c() + theme_bw() + geom_abline(slope=1, intercept=0)
  # ggplot(fake, aes(x= true_median, y = est, col=col)) + geom_point() + scale_color_viridis_c() + theme_bw() + geom_abline(slope=1, intercept=0)
  # 



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
              I_obs=as.vector(IA_obs), ratio = 1)
  
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
fit <- stan(file = paste0("src/kotaro.stan"), 
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



