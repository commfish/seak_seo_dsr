#Extract B1980 for phi prior in 1980 model:
{library(tidyverse)
  library(ggridges)
  library(ggplot2)
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
  library(scales)}

Year <-2023
source("Production_models/Code/SPM_helper.R")

#source("Code/Posterior_Plotting/YE_SPM_posterior_exams_Func.R")

{phase2.3d0<-"PT2i_fullcatch_B1-1_B2-1_upv-3_derb_0_1600k"
  phase2.3dundest<-"PT2i_fullcatch_B1-1_B2-1_upv-3_derb_0.3_1600k"
  phase2.3dovest<-"PT2i_fullcatch_B1-1_B2-1_upv-3_derb_-0.3_1600k"
  
  phase2.5d0<-"PT2i_fullcatch_B1-1_B2-1_upv-5_derb_0_1600k"
  phase2.5dundest<-"PT2i_fullcatch_B1-1_B2-1_upv-5_derb_0.3_1600k"
  phase2.5dovest<-"PT2i_fullcatch_B1-1_B2-1_upv-5_derb_-0.3_1600k"}

{load(file=paste("Production_models/Output/",phase2.3d0,"/post.Rdata", sep=""))
lv3_derb0<-post

load(file=paste("Production_models/Output/",phase2.3dundest,"/post.Rdata", sep=""))
lv3_derb3p<-post

load(file=paste("Production_models/Output/",phase2.3dovest,"/post.Rdata", sep=""))
lv3_derb3o<-post

load(file=paste("Production_models/Output/",phase2.5d0,"/post.Rdata", sep=""))
lv5_derb0<-post

load(file=paste("Production_models/Output/",phase2.5dundest,"/post.Rdata", sep=""))
lv5_derb3p<-post

load(file=paste("Production_models/Output/",phase2.5dovest,"/post.Rdata", sep=""))
lv5_derb3o<-post}

ModList.st2<-list(
  lv3_derb0=lv3_derb0,
  lv3_derb3p=lv3_derb3p,
  lv3_derb3o=lv3_derb3o,
  lv5_derb0=lv5_derb0,
  lv5_derb3p=lv5_derb3p,
  lv5_derb3o=lv5_derb3o)

ModList.st2.sig<-list(
  lv3_derb0=lv3_derb0,
  lv5_derb0=lv5_derb0)

ModList.st2.db<-list(
  lv3_derb0=lv3_derb0,
  lv3_derb3p=lv3_derb3p,
  lv3_derb3o=lv3_derb3o)

B1980.sig3<-MCMCchains(  object=lv3_derb0,  params = "B1980",  excl = NULL,
                    ISB = TRUE,  mcmc.list = FALSE,chain_num = NULL)
B1980.sig5<-MCMCchains(  object=lv5_derb0,  params = "B1980",  excl = NULL,
                         ISB = TRUE,  mcmc.list = FALSE,chain_num = NULL)
K.sig3<-MCMCchains(  object=lv3_derb0,  params = "K",  excl = NULL,
                   ISB = TRUE,  mcmc.list = FALSE,chain_num = NULL)
logK.sig3<-log(K.sig3)

K.sig5<-MCMCchains(  object=lv5_derb0,  params = "K",  excl = NULL,
                     ISB = TRUE,  mcmc.list = FALSE,chain_num = NULL)
logK.sig5<-log(K.sig5)

B1980.d0<-MCMCchains(  object=lv3_derb0,  params = "B1980",  excl = NULL,
                         ISB = TRUE,  mcmc.list = FALSE,chain_num = NULL)
B1980.d3p<-MCMCchains(  object=lv3_derb3p,  params = "B1980",  excl = NULL,
                       ISB = TRUE,  mcmc.list = FALSE,chain_num = NULL)
B1980.d3o<-MCMCchains(  object=lv3_derb3o,  params = "B1980",  excl = NULL,
                       ISB = TRUE,  mcmc.list = FALSE,chain_num = NULL)

K.d0<-MCMCchains(  object=lv3_derb0,  params = "K",  excl = NULL,
                     ISB = TRUE,  mcmc.list = FALSE,chain_num = NULL)
logK.d0<-log(K.d0)

K.d3p<-MCMCchains(  object=lv3_derb3p,  params = "K",  excl = NULL,
                   ISB = TRUE,  mcmc.list = FALSE,chain_num = NULL)
logK.d3p<-log(K.d3p)

K.d3o<-MCMCchains(  object=lv3_derb3o,  params = "K",  excl = NULL,
                   ISB = TRUE,  mcmc.list = FALSE,chain_num = NULL)
logK.d3o<-log(K.d3o)
#--------------------------------------------------------------------
# plot and save posteriors for phi and logK for 4x4
png(paste("Production_models/Figures/B1980_K_prior.png",sep=""),
    width=7,height=6,#width=9.5,height=8.5,
    units="in",res=1200)

par(mfrow=c(2,2), mar=c(4,4,2,1),omi=c(0,0,0,0))
cexc=1.3
# 1) plot K vs PE
{
cols<-c("darkblue","darkred")
plot(density(logK.sig3), main="log(K)",ylab="density",xlab="",col=cols[1],
     ylim=c(0,1.4),cex.axis=cexc, cex.lab=cexc, lwd=1.5)
lines(density(logK.sig5), col=cols[2], lwd=1.5)
legend(x="topright", c("max sigma = 0.22","max sigma = 0.08"), 
  col=c(cols[1],cols[2]),
  text.col=c(cols[1],cols[2]),border=F, bty="n", pch=c(NA),cex=cexc)

#3) plot phi vs PE
cols<-c("darkblue","darkred")
plot(density(B1980.sig3), main="phi (B1980/K)",ylab="",xlab="",col=cols[1],
     ylim=c(0,2.4),cex.axis=cexc, cex.lab=cexc, lwd=1.5)
lines(density(B1980.sig5), col=cols[2], lwd=1.5)
legend(x="topright", c("max sigma = 0.22","max sigma = 0.08"), 
       col=c(cols[1],cols[2]),
       text.col=c(cols[1],cols[2]),border=F, bty="n", pch=c(NA),cex=cexc)

#2) plot K versus halibut
cols<-c("black","forestgreen","darkcyan")
plot(density(logK.d0), main="",ylab="density",xlab="log (K)",col=cols[1],
     ylim=c(0,1.4),cex.axis=cexc, cex.lab=cexc, lwd=1.5)
lines(density(logK.d3p), col=cols[2], lwd=1.5)
lines(density(logK.d3o), col=cols[3], lwd=1.5)
legend(x="topright", c("wcpue ~ pre-IFQ","wcpue < pre-IFQ","wcpue > pre-IFQ"), 
       col=c(cols[1],cols[2],cols[3]),
       text.col=c(cols[1],cols[2],cols[3]),border=F, bty="n", pch=c(NA),cex=cexc)

#4) plot phi vs halibut
cols<-c("black","forestgreen","darkcyan")
plot(density(B1980.d0), main="",ylab="",xlab="phi (B1980/K)",col=cols[1],
     ylim=c(0,2.4),cex.axis=cexc, cex.lab=cexc, lwd=1.5)
lines(density(B1980.d3p), col=cols[2], lwd=1.5)
lines(density(B1980.d3o), col=cols[3], lwd=1.5)
legend(x="topright", c("wcpue ~ pre-IFQ","wcpue < pre-IFQ","wcpue > pre-IFQ"), 
       col=c(cols[1],cols[2],cols[3]),
       text.col=c(cols[1],cols[2],cols[3]),border=F, bty="n", pch=c(NA),cex=cexc)
}
dev.off()
#---------------------------------------------------------------------
# phi prior; normal distribution the best in 2022
# need six: sig3-d0, sig3-d3p, sig3-dPo
#           sig5-d0, sig5-d5p, sig5-dPo
#do one by one and save to Stage-3 multimodel run

ModList.st2<-list(
  lv3_derb0=lv3_derb0,
  lv3_derb3p=lv3_derb3p,
  lv3_derb3o=lv3_derb3o,
  lv5_derb0=lv5_derb0,
  lv5_derb3p=lv5_derb3p,
  lv5_derb3o=lv5_derb3o)

prior.post<-lv3_derb0

phi.prior.post<-MCMCchains(  object=prior.post,  params = "B1980",  excl = NULL,
                             ISB = TRUE,  mcmc.list = FALSE,chain_num = NULL)
logK.prior.post<-log(MCMCchains(  object=prior.post,  params = "K",  excl = NULL,
                             ISB = TRUE,  mcmc.list = FALSE,chain_num = NULL))

K.prior.post<-MCMCchains(object=prior.post,  params = "K",  excl = NULL,
                           ISB = TRUE,  mcmc.list = FALSE,chain_num = NULL)

phi.prior<-rnorm(10000,mean(phi.prior.post),sd(phi.prior.post))
logK.prior<-rnorm(10000,mean(logK.prior.post),sd(logK.prior.post))
K.prior<-rlnorm(10000,log(mean(K.prior.post)),log(sd(K.prior.post)))

par(mfrow=c(1,1))
plot(density(phi.prior.post)); lines(density(phi.prior),col="red")

#trim?
trim<-0.999
{phi.prior.post<-MCMCchains(  object=prior.post,  params = "B1980",  excl = NULL,
                             ISB = TRUE,  mcmc.list = FALSE,chain_num = NULL)
phi.prior.post<-phi.prior.post[phi.prior.post<quantile(phi.prior.post,c(trim))]
phi.prior<-rnorm(1000,mean(phi.prior.post),sd(phi.prior.post))
lines(density(phi.prior),col="blue")}

mean(phi.prior.post)
sd(phi.prior.post)
sd(phi.prior.post)/mean(phi.prior.post)

#------------------------------------------------------------------------------
par(mfrow=c(1,1))
plot(density(logK.prior.post)); lines(density(logK.prior),col="red")

#trim?
trim<-0.95
exp<-1.0
{ logK.prior.post2<-logK.prior.post[logK.prior.post<quantile(logK.prior.post,c(trim))]
  logK.prior2<-rnorm(10000,mean(logK.prior.post2),sd(logK.prior.post)*exp)
  lines(density(logK.prior2),col="blue")}

mean(logK.prior.post)
sd(logK.prior.post)
sd(logK.prior.post)/mean(logK.prior.post)

mean(logK.prior.post2)
sd(logK.prior.post2)
sd(phi.prior.post)/mean(phi.prior.post)

#------------------------------------------------------------------------------
par(mfrow=c(1,1))
K.prior<-rlnorm(10000,log(mean(K.prior.post)),sd(logK.prior.post))
plot(density(K.prior.post)); lines(density(K.prior),col="red")

#trim?
trim<-0.93
exp<-1.0
{ K.prior.post2<-K.prior.post[K.prior.post<quantile(K.prior.post,c(trim))]
  K.prior2<-rlnorm(10000,log(mean(K.prior.post2)),sd(logK.prior.post)*exp)
  lines(density(K.prior2),col="blue",lty=2,lwd=2)}

log(mean(K.prior.post))
log(sd(K.prior.post))
log(sd(K.prior.post))/log(mean(K.prior.post))

log(mean(K.prior.post2))
log(sd(K.prior.post2))
log(sd(K.prior.post2))/log(mean(K.prior.post2))

exp(log(sd(K.prior.post2)))/exp(log(mean(K.prior.post2)))

(sd(K.prior.post2))/(mean(K.prior.post2))






















