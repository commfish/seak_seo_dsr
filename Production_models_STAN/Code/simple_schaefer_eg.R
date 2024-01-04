#rstan not working, trying things from the web
#remove.packages(c("StanHeaders","rstan")); install.packages(c("StanHeaders","rstan"),type="source")

{library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
#source('stan_utility.R')
#lsf.str()
library("bayesplot")
library("rstanarm")
library("ggplot2")
library("shinystan")
library("hexbin")}

#getwd()
#util <- new.env()
#source('C:/Users/pjjoy/Documents/Groundfish Biometrics/SEO_DSR/stan_development/code/stan_utility.R', local=util)

setwd("code/")


C<-c(15.9,25.7,28.5,23.7,25.0,33.3,28.2,19.7,17.5,19.3,21.6,23.1,22.5,22.5,23.6,29.1,14.4,13.2,28.4,34.6,37.5,25.9,25.3) 
I<-c(61.89,78.98,55.59,44.61,56.89,38.27,33.84,36.13,41.95,36.63,36.33,38.82,34.32,37.64,34.01,32.16,26.88,36.61,30.07,30.75,23.36,22.36,21.91)
N<-23

C_cv<-rep(0.01,N)
I_cv<-rep(0.01,N)

I_cv*I

hist(I)

length(I)

tuna_data <- list(C=C,I=I,N=N, C_cv = C_cv, I_cv=I_cv)

set.seed(1901)

#initial values
inits1 <- list(r=0.2, K=500, iq=10, isigma2=100, itau2=100, 
               P=c(0.99,0.98,0.96,0.94,0.92,0.90,0.88,0.86,0.84,0.82,0.80,0.78,0.76,0.74,0.72,0.70,0.68,0.66,0.64,0.62,0.60,0.58,0.56))
inits2 <- list(r=0.3, K=300, iq=8, isigma2=140, itau2=200, #isigma2 is too high...
               P=c(0.99,0.99,0.99,0.99,0.90,0.90,0.90,0.90,0.90,0.85,0.85,0.85,0.75,0.75,0.7,0.7,0.7,0.7,0.65,0.65,0.6,0.5,0.4))
inits3 <- list(r=0.1, K=600, iq=4, isigma2=80, itau2=300, 
               P=c(0.99,0.99,0.99,0.99,0.90,0.90,0.90,0.90,0.90,0.85,0.85,0.85,0.75,0.75,0.7,0.7,0.7,0.7,0.65,0.65,0.6,0.5,0.4))

inits <- list(inits1,inits2, inits3)

fit <- stan(file = 'code/Schaefer_margq2.stan', data = tuna_data, 
              iter = 20000, chains = 3, cores=3,
            init=inits, warmup=5000,verbose=F, thin=100, #test_grad = T,#thin=X,
            control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 30))

#note defaults: adapt_delta=0.8, stepsize = 1, max_treedepth=10
#note values to deal with issues: adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 15

#changing sigma prior to gamma(3.785, 0.102) divergences down to ~20%
#notes: changing r prior to uninformative beta reduced divergences from ~20% to ~11&
#changing sigma to gamma(5,0.102) got down to ~8%
#changing sigma to gamma(6,0.102) and gamma(7,0.102) got down to ~6%
# went down to sigma gamma (9,0.1) but divergences back up to ~9%

launch_shinystan(fit)

check_energy(fit)
check_treedepth(fit)
get_divergent_iterations(fit)
check_divergences(fit)

get_num_divergent(fit)
get_num_max_treedepth(fit) 
get_bfmi(fit) 
get_low_bfmi_chains(fit)


samples <- rstan::extract(fit)
#graphics.off()
#par("mar")
#par(mar=c(0.1,0.1,0.1,0.1))

{par(mfrow=c(2,2))
quants <- matrix(0,nrow=3,ncol=N)
for (II in 1:N)
  quants[,II] <- quantile(samples$Biomass[,II],probs=,c(0.05,0.5,0.95))
ymax <- max(quants)*1.05
plot(1:N,quants[2,],xlab="Year",ylab="Biomass",type="b",lty=1,pch=16,ylim=c(0,ymax))
lines(1:N,quants[1,],lwd=1,lty=2)
lines(1:N,quants[3,],lwd=1,lty=2)

quants <- matrix(0,nrow=3,ncol=N)
for (II in 1:N)
  quants[,II] <- quantile(samples$Inew[,II],probs=,c(0.05,0.5,0.95))
ymax <- max(quants)*1.05
plot(1:N,I,xlab="Year",ylab="Index",type="b",lty=1,pch=16,ylim=c(0,ymax))
lines(1:N,quants[2,],lwd=2,lty=1)
lines(1:N,quants[1,],lwd=1,lty=2)
lines(1:N,quants[3,],lwd=1,lty=2)}

pairs(fit)

util$check_all_diagnostics(fit)
partition <- partition_div(fit)

posterior <- as.matrix(fit)
print(fit)

traceplot(fit, pars = c("r","K","iq","isigma2","itau2"))

 #, permuted = TRUE) # return a list of arrays 

par(mfrow=c(1, 1))

c_light_trans <- c("#DCBCBC80")
c_dark_trans <- c("#8F272780")
c_green_trans <- c("#00FF0080")

{par(mfrow=c(3, 4))

plot(samples$r, samples$K,
     col=c_dark_trans, 
     pch=16, cex=0.8,
     xlab="r", ylab="K")

plot(samples$r, samples$iq,
     col=c_dark_trans, 
     pch=16, cex=0.8,
     xlab="r", ylab="iq")

plot(samples$r, samples$isigma2,
     col=c_dark_trans, 
     pch=16, cex=0.8,
     xlab="r", ylab="isigma2")

plot(samples$r, samples$itau2,
     col=c_dark_trans, 
     pch=16, cex=0.8,
     xlab="r", ylab="itau2")

plot(samples$K, samples$iq,
     col=c_dark_trans, 
     pch=16, cex=0.8,
     xlab="K", ylab="iq")

plot(samples$K, samples$isigma2,
     col=c_dark_trans, 
     pch=16, cex=0.8,
     xlab="K", ylab="isigma2")

plot(samples$K, samples$itau2,
     col=c_dark_trans, 
     pch=16, cex=0.8,
     xlab="K", ylab="itau2")

plot(samples$iq, samples$isigma2,
     col=c_dark_trans, 
     pch=16, cex=0.8,
     xlab="iq", ylab="isigma2")

plot(samples$iq, samples$itau2,
     col=c_dark_trans, 
     pch=16, cex=0.8,
     xlab="iq", ylab="itau2")

plot(samples$isigma2, samples$itau2,
     col=c_dark_trans, 
     pch=16, cex=0.8,
     xlab="isigma2", ylab="itau2")}

unpermuted_samples <- extract(fit, permute=FALSE)

par(mfrow=c(2, 2))

plot(samples$r, samples$K,
     main="Chain 1", col=c_light_trans, pch=16, cex=0.8,
     xlab="r", ylab="K")
points(unpermuted_samples[,1,1], unpermuted_samples[,1,2],
       col=c_dark_trans, pch=16, cex=0.8)

plot(samples$r, samples$K,
     main="Chain 2", col=c_light_trans, pch=16, cex=0.8,
     xlab="r", ylab="K")
points(unpermuted_samples[,2,1], unpermuted_samples[,2,2],
       col=c_dark_trans, pch=16, cex=0.8)

plot(samples$r, samples$K,
     main="Chain 3", col=c_light_trans, pch=16, cex=0.8,
     xlab="r", ylab="K")
points(unpermuted_samples[,3,1], unpermuted_samples[,3,2],
       col=c_dark_trans, pch=16, cex=0.8)

mcmc_areas(posterior,
           pars = c("r"),
           prob = 0.8)
mcmc_areas(posterior,
           pars = c("K"),
           prob = 0.8)
mcmc_areas(posterior,
           pars = c("iq"),
           prob = 0.8)
mcmc_areas(posterior,
           pars = c("isigma2"),
           prob = 0.8)
mcmc_areas(posterior,
           pars = c("itau2"),
           prob = 0.8)

color_scheme_set("darkgray")
mcmc_scatter(
  as.matrix(fit),
  pars = c("r", "K"), 
  np = nuts_params(fit), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
)

mcmc_hex(
  as.matrix(fit),
  pars = c("r", "K"), 
  np = nuts_params(fit), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
)

posterior <- as.array(fit)
np <- nuts_params(fit)


mcmc_pairs(posterior, pars = c("r", "K","iq","isigma2","itau2"))

mcmc_pairs(
  posterior,
  pars = c("r", "K","iq","isigma2","itau2"),
  #transformations = list(sigma = "log"), # show log(sigma) instead of sigma
  #off_diag_fun = "hex" # use hexagonal heatmaps instead of scatterplots
  condition = pairs_condition(nuts = "accept_stat__"),
  np = np
)

bayesplot::color_scheme_set("brightblue")
mcmc_pairs(
  fit,
  pars = c("r", "K","iq","isigma2","itau2"),
  #transformations = list(sigma = "log"),
  condition = pairs_condition(nuts = "accept_stat__"),
  off_diag_args = list(size = 3/4, alpha = 1/3), # size and transparency of scatterplot points
  np_style = pairs_style_np(div_color = "green", div_shape = 2, div_size=0.5), # color and shape of the divergences
  np = np
)

pairs(
  fit,
  pars = c("r", "K","iq","isigma2","itau2"),
  condition = pairs_condition(nuts = "divergent__")
)

mcmc_pairs(
  fit,
  pars = c("r", "K","iq","isigma2","itau2"),
  diag_fun = "dens",
  off_diag_fun = "hex",
  np_style = pairs_style_np(div_color = "green", div_shape = 2)
)

pairs(
  fit,
  pars = c("r", "K","iq","isigma2","itau2"),
  #diag_fun = "dens",
  #off_diag_fun = "hex",
  #np_style = pairs_style_np(div_color = "green", div_shape = 2), # color and shape of the divergences
  condition = pairs_condition(nuts = "divergent__")
)

dev.off()
mcmc_parcoord(fit)

mcmc_pairs(fit)

mcmc_rank_overlay(fit)

color_scheme_set("red")
ppc_dens_overlay(y = fit$r,
                 yrep = posterior_predict(fit, draws = 50))

library("dplyr")
color_scheme_set("brightblue")
fit %>%
  posterior_predict(draws = 500) %>%
  ppc_stat_grouped(y = I,
                   group = C,
                   stat = "median")

color_scheme_set("mix-blue-pink")
p <- mcmc_trace(la,  pars = c("r", "K"), n_warmup = 300,
                facet_args = list(nrow = 2, labeller = label_parsed))
p + facet_text(size = 15)

par(mfrow=c(2,2))

samples <- extract(fit, permuted = TRUE)

quants <- matrix(0,nrow=3,ncol=N)
for (II in 1:N)
  quants[,II] <- quantile(samples$Biomass[,II],probs=,c(0.05,0.5,0.95))
ymax <- max(quants)*1.05
plot(1:N,quants[2,],xlab="Year",ylab="Biomass",type="b",lty=1,pch=16,ylim=c(0,ymax))
lines(1:N,quants[1,],lwd=1,lty=2)
lines(1:N,quants[3,],lwd=1,lty=2)

quants <- matrix(0,nrow=3,ncol=N)
for (II in 1:N)
  quants[,II] <- quantile(samples$Inew[,II],probs=,c(0.05,0.5,0.95))
ymax <- max(quants)*1.05
plot(1:N,I,xlab="Year",ylab="Index",type="b",lty=1,pch=16,ylim=c(0,ymax))
lines(1:N,quants[2,],lwd=2,lty=1)
lines(1:N,quants[1,],lwd=1,lty=2)
lines(1:N,quants[3,],lwd=1,lty=2)

#---------------------------------------------------------------------------------
# Lots of deviations with this model!!! I think it has a hard time with q,K and r being correlated
# Try reparameterizing according to Best and Punt 2019:

set.seed(1901)

#initial values
inits1 <- list(r=0.2, K=500, #Z=rep(0.1,23), #qhat=exp(0.1),#Z=rep(0.1,23), 
               isigma2=100, itau2=100, 
               P=c(0.99,0.98,0.96,0.94,0.92,0.90,0.88,0.86,0.84,0.82,0.80,0.78,0.76,0.74,0.72,0.70,0.68,0.66,0.64,0.62,0.60,0.58,0.56))
inits2 <- list(r=0.3, K=300, #qhat=exp(0.1),#Z=rep(0.1,23), 
               isigma2=140, itau2=200, #isigma2 is too high...
               P=c(0.99,0.99,0.99,0.99,0.90,0.90,0.90,0.90,0.90,0.85,0.85,0.85,0.75,0.75,0.7,0.7,0.7,0.7,0.65,0.65,0.6,0.5,0.4))
inits3 <- list(r=0.1, K=600, #qhat=exp(0.1),#Z=rep(0.1,23), 
               isigma2=80, itau2=300, 
               P=c(0.99,0.99,0.99,0.99,0.90,0.90,0.90,0.90,0.90,0.85,0.85,0.85,0.75,0.75,0.7,0.7,0.7,0.7,0.65,0.65,0.6,0.5,0.4))

inits <- list(inits1,inits2, inits3)

fit <- stan(file = 'Schaefer_margq.stan', data = tuna_data, 
            iter = 50000, chains = 3, cores=3,
            init=inits, warmup=5000,verbose=F, thin=100, #test_grad = T,#thin=X,
            control = list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth = 30))

launch_shinystan(fit)

par(mfrow=c(2,1))
samples <- extract(fit, permuted = TRUE)
quants <- matrix(0,nrow=3,ncol=N)
for (II in 1:N)
  quants[,II] <- quantile(samples$Biomass[,II],probs=,c(0.05,0.5,0.95))
ymax <- max(quants)*1.05
plot(1:N,quants[2,],xlab="Year",ylab="Biomass",type="b",lty=1,pch=16,ylim=c(0,ymax))
lines(1:N,quants[1,],lwd=1,lty=2)
lines(1:N,quants[3,],lwd=1,lty=2)

quants <- matrix(0,nrow=3,ncol=N)
for (II in 1:N)
  quants[,II] <- quantile(samples$Inew[,II],probs=,c(0.05,0.5,0.95))
ymax <- max(quants)*1.05
plot(1:N,I,xlab="Year",ylab="Index",type="b",lty=1,pch=16,ylim=c(0,ymax))
lines(1:N,quants[2,],lwd=2,lty=1)
lines(1:N,quants[1,],lwd=1,lty=2)
lines(1:N,quants[3,],lwd=1,lty=2)

posterior <- as.array(fit)
np <- nuts_params(fit)

bayesplot::color_scheme_set("brightblue")
mcmc_pairs(
  fit,
  pars = c("r", "K","isigma2","itau2"),
  #transformations = list(sigma = "log"),
  condition = pairs_condition(nuts = "accept_stat__"),
  off_diag_args = list(size = 3/4, alpha = 1/3), # size and transparency of scatterplot points
  np_style = pairs_style_np(div_color = "green", div_shape = 2, div_size=0.5), # color and shape of the divergences
  np = np
)

pairs(
  fit,
  pars = c("r", "K","isigma2","itau2"),
  condition = pairs_condition(nuts = "divergent__")
)

mcmc_pairs(
  fit,
  pars = c("r", "K","isigma2","itau2"),
  diag_fun = "dens",
  off_diag_fun = "hex",
  np_style = pairs_style_np(div_color = "green", div_shape = 2)
)

pairs(
  fit,
  pars = c("r", "K","iq","isigma2","itau2"),
  #diag_fun = "dens",
  #off_diag_fun = "hex",
  #np_style = pairs_style_np(div_color = "green", div_shape = 2), # color and shape of the divergences
  condition = pairs_condition(nuts = "divergent__")
)

#-------------------------------------------------------------------------------
# another example with rstanarm
color_scheme_set("purple")

fit <- stan_glmer(mpg ~ wt + (1|cyl), data = mtcars)
ppc_intervals(
  y = mtcars$mpg,
  yrep = posterior_predict(fit),
  x = mtcars$wt,
  prob = 0.5
) +
  labs(
    x = "Weight (1000 lbs)",
    y = "MPG",
    title = "50% posterior predictive intervals \nvs observed miles per gallon",
    subtitle = "by vehicle weight"
  ) +
  panel_bg(fill = "gray95", color = NA) +
  grid_lines(color = "white")


