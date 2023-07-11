############################################################################################################
# Pella Tomlinson Surplus Production Model (Bayesian State-Space)
# Use with Data from "SEOsubd_YE_Production_Discard_Model_DATA_LOAD.r"
# Abundance indices from submersible surveys and IPHC longline survey CPUE data
# Harvests from ADF&G-CFD
# Discards estimated from IPHC longline survey data and Halibut harvests (see IPHC_Survey_Discard_Est.r)
# with future projections
###########################################################################################################
{
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
  source("Code/2022_DSR_SAFE_models/Phase1/DATALOAD_SEO_YE_SPM_Func_1980.R")
}
###########################################################################################################

#=========================================================================================
#name the model
Mod.title<-"PHASE_3_sim"
#=========================================================================================
cat('data { 

#Bycatch Model and Discard Calculations =========================================
for (i in 1:Subd){
for (t in 1:15) {   #Derby years pre-1995
  logBy[i,t] <-log(realBy[i,t]) #~ dunif(-10,10)
  tau.log.By[i,t] <- 1 / log(cv.ExpByc[i,t]*cv.ExpByc[i,t] + 1)
  ExpByc[i,t] ~ dlnorm(logBy[i,t],tau.log.By[i,t]) #T(0.01,)  #(50,)
  By[i,t]<-exp(logBy[i,t])
  logD[i,t]<-log(max(1,(By[i,t]*exp(derb.err[t]))-Lnd.By[i,t]))
  D[i,t]<-exp(logD[i,t])  
}
for (t in 16:N) {   #IFQ years 1995 on
  logBy[i,t] <-log(realBy[i,t]) #~ dunif(-10,10) 
  tau.log.By[i,t] <- 1 / log(cv.ExpByc[i,t]*cv.ExpByc[i,t] + 1)
  ExpByc[i,t] ~ dlnorm(logBy[i,t],tau.log.By[i,t])T(0.01,)  #(50,)
  By[i,t]<-exp(logBy[i,t])
  logD[i,t]<-log(max(1,By[i,t]-Lnd.By[i,t]))
  D[i,t]<-exp(logD[i,t]) 
}
}

tau.log.derb<-1/(DEsd*DEsd+1)
for (t in 1:N){
  derb.err[t] <- depsilon[t]-(DEsd*DEsd/2)
  depsilon[t] ~ dnorm(Derby.Eff,tau.log.derb)T(Derby.Eff-0.3,Derby.Eff+0.3)
}

#Known Catches
for (i in 1:Subd){
for (t in 1:N){
  logKnC[i,t] <- log(realC[i,t]) #~ dunif(-10,10)
  KnC.obs[i,t] ~ dlnorm(logKnC[i,t], tau.log.KnC[i,t])  
  tau.log.KnC[i,t]<-1/log(cv.KnC[i,t]*cv.KnC[i,t]+1)  
  KnC[i,t]<-exp(logKnC[i,t])

#Total catches
  logC[i,t]<-log(D[i,t]+KnC[i,t])
  C[i,t]<-exp(logC[i,t])
}
}

## POPULATION MODEL...=====================================
#Process error 
### !!! SIMULATED!!! logvar ~ dunif(-10,upvar) #dunif(-10,0.5)   #log variance prior = very restrictive to keep it small
sigma <- sqrt(exp(logvar)) 
tau.log.pe <- 1/log(sigma*sigma+1)
for (t in 1:(N+Fu)) {
  PE[t]<-epsilon[t]-(sigma*sigma/2)
  epsilon[t] ~ dnorm(0,tau.log.pe)T(-0.1,0.1) #epsilon truncation = +/- 10%  
  #epsilon[t] <- 0
}

#Year 1
for (i in 1:Subd){
  #K[i]<-exp(logK[i])
  logB[i,1] <-log(phi[i]*K[i]) 
  B[i,1]<-exp(logB[i,1])
  S[i,1]<-0
}

#Year 2 through N
for (i in 1:Subd){
for( t in 2:N )   #PT with Biomass starting in year 2
  { 
    logB[i,t]<-log(max(B[i,t-1]+(r[i]/p)*B[i,t-1]*(1-(B[i,t-1]/K[i])^p)-C[i,t-1],1)*exp(PE[t-1]))
    B[i,t]<-exp(logB[i,t])
    S[i,t]<-B[i,t-1]-B[i,t]+C[i,t]
  }
}

#surplus production 
for (i in 1:Subd){
for (t in 1:N){
    Surplus[i,t]<-(r[i]/p)*B[i,t]*(1-(B[i,t]/K[i])^p) 
}
}

#Biomass likelihoods...============================================
for (i in 1:Subd){
for (t in 1:30) {          #Submersible surveys
    use.cv1[i,t]<-sqrt(cv.B[i,t]*cv.B[i,t]+Tau1*Tau1)
    tau.log.B1[i,t] <- 1 / log(use.cv1[i,t]*use.cv1[i,t] + 1)
    B.obs[i,t] ~ dlnorm(logB[i,t],tau.log.B1[i,t])
}

for (t in 31:N) {        #ROV surveys
    use.cv2[i,t]<-sqrt(cv.B[i,t]*cv.B[i,t]+Tau2*Tau2)
    tau.log.B2[i,t] <- 1 / log(use.cv2[i,t]*use.cv2[i,t] + 1)
    B.obs[i,t] ~ dlnorm(logB[i,t],tau.log.B2[i,t])
}

for (t in 1:N){         #IPHC CPUE 
    use.cv3[i,t]<-sqrt(cv.sCPUE[i,t]*cv.sCPUE[i,t]+Tau3*Tau3)
    log.qsCPUEsubmean[i,t] <- log(qsCPUE[i] * B[i,t])
    tau.log.sCPUE[i,t] <- 1 / log(use.cv3[i,t]*use.cv3[i,t] + 1)
    sCPUE[i,t] ~ dlnorm(log.qsCPUEsubmean[i,t],tau.log.sCPUE[i,t])
}
}


#PROJECTIONS:================================================================
#future catches and biomass...

#1) future error based on miscalculating biomass and inaccuracy of fishing
logvarF ~ dunif(-10,-4) #dunif(-10,0.5)   #log variance prior = very restrictive to keep it small
sigmaF <- sqrt(exp(logvar)) 
tau.log.fe<-1/log(sigmaF*sigmaF+1)
for (t in (N+1):(N+Fu)) {
  FE[t]<-F.eps[t]-(sigmaF*sigmaF/2)
  F.eps[t] ~ dnorm(0,tau.log.fe)T(-0.1,0.1) #epsilon truncation = +/- 10%  
}

  for (i in 1:Subd){
  for (t in (N+1):(N+Fu)){                              #this might not work because need the cv for B?
    logFC[i,t]<-log(B[i,t-1]*exp(FE[t])*OFLABC)
    FC[i,t]<-exp(logFC[i,t])
  }
  }
  
  #gap year to deal with offset catches and biomass...
  for (i in 1:Subd){
    logB[i,N+1]<-log(max(B[i,N]+(r[i]/p)*B[i,N]*(1-(B[i,N]/K[i])^p)-C[i,N],1)*exp(PE[N]))
    B[i,N+1]<-exp(logB[i,N+1])
  }
  
  #rest of projection
  for (i in 1:Subd){
  for( t in (N+2):(N+Fu)){ 
    logB[i,t]<-log(max(B[i,t-1]+(r[i]/p)*B[i,t-1]*(1-(B[i,t-1]/K[i])^p)-FC[i,t-1],1)*exp(PE[t-1]))
    B[i,t]<-exp(logB[i,t])
  }
  }


#Subdistrict metrics==============================================================
for (i in 1:Subd){
    CBtoK[i]<-B[i,N]/K[i]
    FBtoK[i]<-B[i,N+Fu]/K[i]
    FBtoCB[i] <-  B[i,N+Fu]/B[i,N]
    MSY[i]<-r[i]*K[i]/((p+1)^((p+1)/p)) #r*K/4 for Schaefer
    Bmsy[i]<-0.4*K[i]  #0.5 for Schaefer
    Fmsy[i]<-MSY[i]/Bmsy[i]
    Hmsy[i]<-r[i]/(1+p) 
    Stock.Status[i]<-B[i,N]/(0.4*K[i])
    
    for (t in 1:N){
    CtoB[i,t]<-C[i,t]/B[i,t]
    }
}

#SEO calculations and summations
for (t in 1:(N+Fu)){
   Bseo[t]<-sum(B[,t])
}

for (t in 1:N){
  Cseo[t]<-sum(C[,t])
  KnCseo[t]<-sum(KnC[,t])
  Byseo[t]<-sum(By[,t])
  Dseo[t]<-sum(D[,t])
  S.seo[t]<-sum(S[,t])
  Surplus.seo[t]<-sum(Surplus[,t])
}

Stock.Status.SEO<-Bseo[N]/(0.4*Kseo)  
#Kseo<-sum(K)
Bmsyseo<-sum(Bmsy)
Fmsyseo<-sum(MSY)/Bmsyseo
Fmsyseo2<-(Fmsy[1]*Bmsy[1]+Fmsy[2]*Bmsy[2]+Fmsy[3]*Bmsy[3]+Fmsy[4]*Bmsy[4])/Bmsyseo
FBtoKseo<-Bseo[N+Fu]/Kseo
CBtoKseo<-Bseo[N]/Kseo

#Proj1_biomass<-Bseo[N+1]

#subdistrict level priors =================================================
for (i in 1:Subd){
  ## !!! SIMULATED !!! r[i] ~ dbeta(rB1,rB2)T(0.0001,0.2) #dlnorm(log(R.hyp),0.5)T(0.0001,0.1) #dunif(0.01,0.4)
  ## !!! SIMULATED !!!  logqs[i] ~ dunif(-10,20) #dunif(2,20)
  qsCPUE[i]<-exp(logqs[i])
  
  #new priors from 1888 model output
  K[i] <- Kseo*pi[i] #exp(logK[i])
  
  ## !!! SIMULATED !!!  phi[i] ~ dnorm(bigphi,invtau2)T(0.01,1)
  }

#priors/hyper priors ======================================================
p <- 0.18815 #0.18815 = Bmsy =0.4K; 1=Schaefer; 1e-08 ~ modfied fox
## !!! SIMULATED !!!  Tau1 ~ dunif(0.01,1)
Tau2 <- 0 #~ dunif(0.01,1)  #ROV years with no extra variance
## !!! SIMULATED !!!  Tau3 ~ dunif(0.01,1)

#Rhyper prior... 
## !!! SIMULATED !!!  R.hyp ~ dbeta(B1,B2)T(0.0001,0.2) #set to M=0.02,Z=0.05,AaM=17
## !!! SIMULATED !!!  eta<-exp(logeta)
## !!! SIMULATED !!!  logeta ~ dlogis(logn,1)T(-5,7.55) #cap of 10 gives max eta of 20K; 9.5 ~ 13K, 9.25~10K
## !!! SIMULATED !!!  logn<-log(100)  #indicates shrinkage is uniformly distributed between 0 and 1
## !!! SIMULATED !!!  rB1<-R.hyp*eta
## !!! SIMULATED !!!  rB2<-(1-R.hyp)*eta

#new priors from 1888 model
#bigKtau <- 0 #1/(log(bigKsigma*bigKsigma+1))
## !!! SIMULATED !!!  logKseo ~ dnorm(bigKmu,bigKtau)
Kseo <- exp(logKseo)
alpha <- c(1,1,1,1)
## !!! SIMULATED !!!  pi[1:Subd] ~ ddirch(alpha)

## !!! SIMULATED !!!  phitau <- 1/(log(phisig*phisig+1))
## !!! SIMULATED !!!  bigphi ~ dnorm(phimu,phitau) #T(0.01,1)
## !!! SIMULATED !!!  invtau2 ~ dgamma(a_t, b_t) #with ones is uninformative... 
## !!! SIMULATED !!!  a_t<-1
## !!! SIMULATED !!!  b_t<-1
## !!! SIMULATED !!!  phiTau<-1/invtau2
}
model{
fake<-0
}', file=Mod.title) 

#beta formulation from Albert, J. and J. Hu.  2019.  Probability and Bayesian Modeling.  Chapman and Hall/CRC New York.  552. pg. 

#===========================================================================================================
# prep model setting with data, initial values
{Derby.list<-c(0) #,0.3,-0.3)
DEsdlist<-c(0.1,0.1,0.1)
B1list<-c(1) #c(1,1.483,1.247,1.241,1.374)
B2list<-c(1)}

Data<-load.data(YEAR=2022,
                Derby.Eff = Derby.list[1],
                DEsd=0.1,  #this is CV for derby
                B1=B1list[1],
                B2=B2list[1])
list2env(Data,.GlobalEnv)
#Set future management strategy
OFLABC<-0.00  #0.02 (ABC) 0.032 (OFL). 0=nofishing
HarvStrat<-0 #2.58 #1.96 #1.68 #0, 1.15
Fu<-0

res.to.sim<-"PHASE3_B1-1_B2-1_upv-5_phmu-0.7_phsig-1.2_Kmu-10.6_Ksig-9.2_derb_0_1500k"

ests<-read.csv(paste0("Model Output/",res.to.sim,"/results_sum.csv"))
ests<-as.data.frame(ests)

r<-ests$X50[grep(c("r\\["),ests$parameter)]
logqs<-ests$X50[grepl("qsCPUE",ests$parameter)==T]
Tau1<-ests$mean[grepl("Tau1",ests$parameter)==T]
Tau3<-ests$mean[grepl("Tau3",ests$parameter)==T]
pi<-ests$mean[grepl("pi",ests$parameter)==T]
phi<-ests$mean[grep(c("phi\\["),ests$parameter)]
logKseo<-log(ests$X50[ests$parameter=="Kseo"])
logvar<-log(ests$mean[ests$parameter=="sigma"])
bigKsigma<-c(0.269)   #c(0.269)
 
{

data<-list(N=N,Subd=Subd,
           realC=KnC.obs, 
           cv.KnC=cv.KnC,
           realBy = ExpByc, 
           cv.ExpByc = cv.ExpByc,
           Lnd.By=Lnd.By, 
           #B.obs=B.obs,
           cv.B=cv.B, #CKP=CKP,
           #sCPUE=sCPUE, 
           cv.sCPUE=cv.sCPUE, #,
           #fCPUE=sCPUE, cv.fCPUE=cv.fCPUE, 
           Fu=Fu, 
           OFLABC=OFLABC, #B.err=B.err,
           Derby.Eff = Derby.Eff, DEsd=DEsd,
           #B1=B1,B2=B2
           #r = c(0.02,0.029,0.018,0.018),  #should I be using mean or median?
           r = r,
           logqs=logqs,
           Tau1=Tau1,
           Tau3=Tau3,
           pi=pi,
           phi=phi,
           logKseo=logKseo,
           bigKsigma=bigKsigma,
           logvar=logvar
           )


parameters=c(
  "B.obs","sCPUE" ,"KnC.obs","ExpByc"
)               
}

#=========================================================================================================
#### run JAGS ####


{
  set.seed(2054)
  sim <- jagsUI::jags(model.file=Mod.title, data=data,
                               parameters.to.save=parameters, #inits=inits,
                               n.chains=1, parallel=F, n.iter=1,
                               n.burnin=0 )#, n.thin=niter/1000)  # can mess with n.adapt=
  #or
  #sim <- run.jags(model.file=Mod.title, data = data, monitor=parameters, sample=1, n.chains=1, summarise=FALSE)
  
  #dir.create(paste("Sim.Dat/",res.to.sim,"_",niter/1000000,"m",sep=""))
  #save(sim,file=paste("Model Output/", Mod.title,"_",niter/1000000,"m/post.Rdata", sep=""))
}

#***********************************************************************************
#* Now need to extract fake data for years where we have real data...
Simulated <- coda::as.mcmc(sim)
Simulated
## Markov Chain Monte Carlo (MCMC) output:
## Start = 5001 
## End = 5001 
## Thinning interval = 1 
str(Simulated)
## [1]  1 30
dat = as.data.frame(Simulated$sims.list)

sim.B<-Simulated$sims.list$B.obs
sim.KnC<-Simulated$sims.list$KnC.obs
sim.sCPUE<-Simulated$sims.list$sCPUE
sim.eby<-Simulated$sims.list$ExpByc

sim.B<-matrix(sim.B,nrow=4)
sim.B[which(is.na(B.obs),arr.ind=T)]<-NA

sim.sCPUE<-matrix(sim.sCPUE,nrow=4)
sim.sCPUE[which(is.na(sCPUE),arr.ind=T)]<-NA

sim.KnC<-matrix(sim.KnC,nrow=4)
sim.eby<-matrix(sim.eby,nrow=4)
#export some fake data:


#======================================================================================================
#Now run the fake data in the real model
#go to check_sims.R
B1list<-c(1,1.483,1.247,1.241,1.374)
B2list<-c(1,22.908,31.478,53.481,106.057)
B1list<-c(1)  
B2list<-c(1)  

phiB1<-c(2.339) #from Phase2
phiB2<-c(1.22)

phimu.list<-c(0.7420787,0.6726209,0.7702393) #c(0.6802432) #c(0.7420787) #c(0.6802432)
phisig.list<-c(0.1796734,0.1690635,0.1795257) #c(0.2054276) #c(0.1796734) #c(0.2054276)

bigKmu.list<-c(10.62973,10.75576,10.58545) #c(10.65388) #c(10.6366) 
bigKsigma.list<-c(9.237758,9.390445,9.140554) #c(0.2812395) #c(0.2940348)

upvarlist<- c(-5)

HarvStrat<-0 #2.58 #1.96 #1.68 #0, 1.15
Fu<-1
#Derby.Eff<-Derby.list[d]

Data<-load.data(YEAR=2022,
                Derby.Eff = Derby.list[d],
                DEsd=0.1,  #this is CV for derby
                B1=B1list[b],
                B2=B2list[b])
list2env(Data,.GlobalEnv)

niter<-1000000 #300000 #350000
burnin<-300000
#5K, 2.5m.
#600K - 4.8 hours; really close to converged for both norm- and beta-phi models
#700K ~ 6.3 hours; not converged but close.  small r only bounded to 0.3.  to 0.2 will help  Probably 1m like phase 1 gets us there...
set.seed(1234)
#=======================================================================================================
#Preliminary check of model run
summary(post)

tbl<-jags.View(post, title="", digits=3)
write.csv(tbl,
          file=paste("Model Output/", Mod.title,".csv", sep=""))

MCMCtrace(post, params=c("R.hyp" ,"r","Kseo","K","phi","phi.hyp",#"qfCPUE","qsCPUE",
                         "Tau1","Tau2","sigma","epsilon"),
          ISB=TRUE, pdf=FALSE, Rhat=TRUE, file="Model Output/Param_trace.pdf")
MCMCtrace(post, params=c("MSY","Bmsy","Fmsy","Hmsy","CBtoK","FBtoK","FBtoCB"), ISB=TRUE, pdf=FALSE, Rhat=TRUE)
MCMCtrace(post, params=c("Bseo"),
          ISB=TRUE, pdf=FALSE, Rhat=TRUE, file="Model Output/Param_trace.pdf")
MCMCtrace(post, params=c("FC"),
          ISB=TRUE, pdf=FALSE, Rhat=TRUE, file="Model Output/Param_trace.pdf")
MCMCtrace(post, params=c("epsilon"),
          ISB=TRUE, pdf=FALSE, Rhat=TRUE, file="Model Output/Param_trace.pdf")

#===========================================================================================
#PLOTTING  moved to separate files:
# 1) Code/Posterior_Plotting/YE_SPM_posterior_exams_Func.R
# 2) Code/YE_SPM_posterior_exams.R

##########################################################################################
##########################################################################################
### SCRAP ################################################################################
##########################################################################################
##########################################################################################

load(file="Model Output/SEOsubd_PTmod_FULL4_yes94_unif_hB_rchange_ABC005.Rdata")  # will be Testing or TrueB

densityplot(post, parameters=c('r','K'), layout=NULL, ask=NULL)
tbl<-jags.View(post, title="", digits=3)
tbl<-as.data.frame(tbl)
write.csv(tbl,
     file=paste(wd,"/Model Output/", Mod.title,".csv", sep=""))
pp.check(post, observed="B[", simulated="B[", xlab='Observed data', ylab='Simulated data',
         main='Posterior Predictive Check', ...)
traceplot(post, parameters=c('r','K'), Rhat_min=1.5, layout=c(2,2), ask=FALSE)
whiskerplot(post,parameters=c('r','K'))

MCMCchains(  object=post,  params = "r",  excl = NULL,
  ISB = TRUE,  #ignore square brackets
  mcmc.list = FALSE,chain_num = NULL
)

MCMCplot(post,params=c('B')) 

MCMCpstr(post, func = function(x) quantile(x, probs = c(0.01, 0.99)))

MCMCsummary(post, params = 'r', round = 5)

MCMCtrace(mdl, params = 'r', pdf = FALSE)

#dev.off()
MCMCtrace(post, params=c("R.hyp" ,"r","Kseo","K","phi","phi.hyp","qfCPUE","qsCPUE",
                         "Tau1","Tau2","Tau3","Tau4","sigma","epsilon","hB1","hB2","pB1","pB2"),
          ISB=TRUE, pdf=FALSE, Rhat=TRUE, file="Model Output/Param_trace.pdf")
MCMCtrace(post, params=c("MSY","Bmsy","Fmsy","Hmsy","CBtoK","FBtoK","FBtoCB"), ISB=TRUE, pdf=FALSE, Rhat=TRUE)
MCMCtrace(post, params=c("Bseo"),
          ISB=TRUE, pdf=FALSE, Rhat=TRUE, file="Model Output/Param_trace.pdf")
MCMCtrace(post, params=c("CBtoK"),
          ISB=TRUE, pdf=FALSE, Rhat=TRUE, file="Model Output/Param_trace.pdf")
MCMCtrace(post, params=c("B"), ISB=TRUE, pdf=FALSE, Rhat=TRUE, ind=TRUE,file="SimpT1.pdf")
MCMCtrace(diag, params=c("C"), ISB=TRUE, pdf=FALSE, Rhat=TRUE, ind=TRUE,file="SimpT1.pdf")
MCMCtrace(diag, params=c("D"), ISB=TRUE, pdf=FALSE, Rhat=TRUE, ind=TRUE,file="SimpT1.pdf")
MCMCtrace(post, params=c("K"), ISB=TRUE, pdf=FALSE, Rhat=TRUE, ind=TRUE,file="SimpT1.pdf")
MCMCtrace(post, params=c("r","R.hyp"), ISB=TRUE, pdf=FALSE, Rhat=TRUE, ind=TRUE,file="SimpT1.pdf")
MCMCtrace(diag, params=c("sigma"), ISB=TRUE, pdf=FALSE, Rhat=TRUE, ind=TRUE,file="SimpT1.pdf")
MCMCtrace(diag, params=c("CtoB"), ISB=TRUE, pdf=FALSE, Rhat=TRUE, ind=FALSE,file="Model Output/Param_trace.pdf")

#============================================================================================================
# FUNCTIONS FOR POSTERIOR PLOTTING
{
envplot <- function(df, xaxis, n=0, cols=c(4,2), lwd50=3, lwd90=1, ylim=NULL,...) {   #### redo this as envelope
  means <- colMeans(df)
  qmat <- apply(df, 2, quantile, p=c(0.05,0.25,0.5,0.75,0.95))
  if(is.null(ylim)) ylim <- c(0,max(qmat,na.rm=T))
  if(is.null(xaxis)) xaxis <- 1:ncol(df)
  plot(NA,xlim=range(xaxis),ylim=ylim,...=...)
  # segments(1:ncol(df),qmat[1,],y1=qmat[5,],lwd=lwd90,col=col,lend=1)
  # segments(1:ncol(df),qmat[2,],y1=qmat[4,],lwd=lwd50,col=col,lend=1)
  x1 <- xaxis[1:(ncol(df)-n)]
  polygon(c(x1,rev(x1)),c(qmat[1,1:(ncol(df)-n)],qmat[5,rev(1:(ncol(df)-n))]),border=NA,col=adjustcolor(cols[1], alpha.f=.15))
  polygon(c(x1,rev(x1)),c(qmat[2,1:(ncol(df)-n)],qmat[4,rev(1:(ncol(df)-n))]),border=NA,col=adjustcolor(cols[1], alpha.f=.2))
  x2 <- xaxis[(ncol(df)-n):ncol(df)]
  polygon(c(x2,rev(x2)),c(qmat[1,(ncol(df)-n):ncol(df)],qmat[5,rev((ncol(df)-n):ncol(df))]),border=NA,col=adjustcolor(cols[2], alpha.f=.15))
  polygon(c(x2,rev(x2)),c(qmat[2,(ncol(df)-n):ncol(df)],qmat[4,rev((ncol(df)-n):ncol(df))]),border=NA,col=adjustcolor(cols[2], alpha.f=.2))
  lines(x2,qmat[3,(ncol(df)-n):ncol(df)],col=adjustcolor(cols[2],alpha.f=.4),lwd=3)
  lines(x1,qmat[3,1:(ncol(df)-n)],col=adjustcolor(cols[1],alpha.f=.4),lwd=3)
  #points(x2,means[(ncol(df)-n):ncol(df)],col=adjustcolor(cols[2],alpha.f=.4),pch=16)
  #points(x1,means[1:(ncol(df)-n)],col=adjustcolor(cols[1],alpha.f=.4),pch=16)
  #abline(v=xaxis[(ncol(df)-n)],lty=3)
  # points(means)
}
addenv<-function(df, cols=c("blue","blue")){ #df<-B.ey
  means<-colMeans(df)
  qmat <- apply(df, 2, quantile, p=c(0.05,0.25,0.5,0.75,0.95))
  x1<-Y #length(Y)
  polygon(c(x1,rev(x1)),c(qmat[1,1:(ncol(df)-0)],qmat[5,rev(1:(ncol(df)-0))]),border=NA,col=adjustcolor(cols[1], alpha.f=.15))
  polygon(c(x1,rev(x1)),c(qmat[2,1:(ncol(df)-0)],qmat[4,rev(1:(ncol(df)-0))]),border=NA,col=adjustcolor(cols[1], alpha.f=.2))
  
  lines(x1,qmat[3,1:(ncol(df)-0)],col=adjustcolor(cols[1],alpha.f=.4),lwd=3)
}
par.ext<-function(par,years,areai){
  for (i in 1:years){  #i<-1  A<-1
    cl<-paste(par,"[",areai,",",i,"]",sep="")
    if (i == 1){
      cl<-paste(par,"[",areai,",",i,"]",sep="")
      parl<-c(cl)
    } else {
      parl<-c(parl,cl)
    }
  }
  out<-MCMCchains(post,params=parl,ISB=FALSE,exact=TRUE)
  return(out)
}
midpoint<-function(post,par){
  mid<-quantile(MCMCchains(post,par,ISB=FALSE,exact=TRUE),c(0.5))
  return(mid)
}
plot4<-function(parp="B",years=62,datap=datap,datap2=datap2,scaler=scaler,ymax=0.1){  #
  EY<-par.ext(par=parp,years=years,areai=1)
  NS<-par.ext(parp,years,2)
  CS<-par.ext(parp,years,3)
  SS<-par.ext(parp,years,4)
  par(mfrow=c(2,2))
  vars<-list(EY,NS,CS,SS)
  sdlist<-c("EYKT", "NSEO", "CSEO", "SSEO")
  for (i in 1:4){  #i<-1
    envplot(vars[[i]],1980:(1980+years),cols=c("blue","blue"),n=0,ylab=parp,xlab="Year",
            ylim=c(0,ymax), main=sdlist[i])
    points(Y[2:30],datap[i,2:30], col="blue", pch=18)
    points(Y[31:42],datap[i,31:42], col="blue", pch=18)
    if (scaler == 0){
      qsc<-1
    } else {
      parpl<-paste(scaler,"[",i,"]",sep="")
      qsc<-midpoint(post=post,par=parpl)
    }
    points(Y[1:42],datap2[i,]/qsc, col="orange", pch=18)
    abline(v=2022,  col="purple", lty=2)
  }
}
}
#===============================================================
# MCMCchains example
Bs <-MCMCchains(object=post,  params = c("B[1,1]"),  excl = NULL,
                ISB = FALSE,  #ignore square brackets
                mcmc.list = FALSE,chain_num = NULL, exact=TRUE
); str(Bs); summary(Bs)
#===============================================================
{N+Fu
B.ey<-par.ext(par="B",years=N+Fu,areai=1)
B.ns<-par.ext("B",N+Fu,2)
B.cs<-par.ext("B",N+Fu,3)
B.ss<-par.ext("B",N+Fu,4)
B.all<-MCMCchains(post,"Bseo")
#Bs<-Bs[,1:42]
par(mfrow=c(1,1))
envplot(B.all,1980:2046,cols=c("goldenrod4","goldenrod4"),n=0,ylab="Biomass (t)",xlab="Year",ylim=c(0,40000))
Y<-seq(1980,2046,1)
addenv(B.ey,cols=c("blue","blue"))
addenv(B.ns, cols=c("red","red"))
addenv(B.cs, cols=c("forestgreen","forestgreen"))
addenv(B.ss, cols=c("purple","purple"))
}
#Y<-seq(1980,2021,1)
par(mfrow=c(2,2))
vars<-list(B.ey,B.ns,B.cs,B.ss)
for (i in 1:4){  #i<-1
  envplot(vars[[i]],1980:2046,cols=c("blue","blue"),n=0,ylab="Biomass (t)",xlab="Year",ylim=c(0,12000))
  points(Y[2:30],B.obs[i,2:30], col="blue", pch=18)
  points(Y[31:42],B.obs[i,31:42], col="darkblue", pch=18)
  #parpl<-paste("qsCPUE[",i,"]",sep="")
  #qCPUEm<-midpoint(post=post,par=parpl)
  #points(Y[1:42],sCPUE[i,]/qCPUEm, col="orange", pch=18)#qsCPUEm, col="orange", pch=18)
  abline(v=2022,  col="purple", lty=2)
}
Y<-seq(1980,2046,1)
plot4(parp="B", years=N+Fu,datap=B.obs,datap2=sCPUE,scaler="qsCPUE", ymax=15000)
plot4(parp="B", years=N+Fu,datap=B.obs,datap2=fCPUE,scaler="qfCPUE", ymax=15000)
#plot4(parp="YEHA.ratio", years=42,datap=wcpue,datap2=By.obs/Hal,scaler=0, ymax=0.12)
plot4(parp="C", years=42,datap=KnC.obs,datap2=KnC.obs,scaler=0, ymax=400)
plot4(parp="By", years=42,datap=ExpByc,datap2=Lnd.By,scaler=0, ymax=200)
plot4(parp="D", years=42,datap=ExpByc-Lnd.By,datap2=ExpByc-Lnd.By,scaler=0, ymax=150)
plot4(parp="KnC", years=42,datap=KnC.obs,datap2=KnC.obs,scaler=0, ymax=400)
plot4(parp="Surplus", years=42,datap=KnC.obs,datap2=KnC.obs,scaler=0, ymax=250)
plot4(parp="S", years=42,datap=KnC.obs,datap2=KnC.obs,scaler=0, ymax=800)

Surplus.all<-MCMCchains(post,"Surplus.seo")
par(mfrow=c(1,1))
envplot(Surplus.all,1980:2022,cols=c("goldenrod4","goldenrod4"),n=0,ylab="Surplus (t)",xlab="Year",ylim=c(0,500))
Cg<-MCMCchains(post,"Cseo")
addenv(Cg, cols=c("red","red"))

{
par(mfrow=c(2,2))
sdlist<-c("EYKT", "NSEO", "CSEO", "SSEO")
for (i in 1:4){
Cg<-par.ext(par="C",years=42,areai=i)  
Byg<-par.ext("By",42,i)
Dg<-par.ext("D",42,i)
KnCg<-par.ext("KnC",42,i)
envplot(Cg,1980:2021,cols=c("blue","blue"),n=0,ylab="Catches (t)",xlab="Year",
        ylim=c(0,400), main=sdlist[i])
Y<-seq(1980,2021,1)
addenv(Byg,cols=c("red","red"))
lines(Lnd.By[i,]~Y, col="darkred", lwd=2, lty=2)
addenv(Dg, cols=c("black","black"))
addenv(KnCg, cols=c("forestgreen","forestgreen"))
}
}

{
  par(mfrow=c(1,1))
    Cg<-MCMCchains(post,"Cseo")
    Byg<-MCMCchains(post,"Byseo")
    Dg<-MCMCchains(post,"Dseo")
    KnCg<-MCMCchains(post,"KnCseo")
    envplot(Cg,1980:2021,cols=c("blue","blue"),n=0,ylab="Catches (t)",xlab="Year",ylim=c(0,1000))
    Y<-seq(1980,2021,1)
    addenv(Byg,cols=c("red","red"))
    Lndby.tot<-Lnd.ey+Lnd.ns+Lnd.cs+Lnd.ss
    lines(Lndby.tot~Y, col="darkred", lwd=2, lty=2)
    addenv(Dg, cols=c("black","black"))
    addenv(KnCg, cols=c("forestgreen","forestgreen"))
 
}
#==============================================================================
#Plots from SEO as a whole for reference...
envplot(C_try,1980:2022,n=1,ylab="Total Catch (t)",xlab="Year", cols=c("grey","grey"))
addenv(By_try, cols=c("blue","blue"))
lines(Lnd.By~Y, col="royalblue", lwd=2, lty=2)
addenv(D_try, cols=c("black","black"))
addenv(KnC_try, cols=c("orange","orange"))

C_try <-MCMCchains(post,"C")# pullpost(Simple,"C[")
envplot(C_try,1980:2022,n=1,ylab="Total Catch (t)",xlab="Year", cols=c("red","red"))
points(Y,KnC.obs+D.obs,col="purple", pch=18)
points(Y,KnC.obs+wcpue*Hal-Lnd.By, col="seagreen", pch=12)
points(Y,KnC.obs+Hal*(mean(wcpue, na.rm=T))-Lnd.By, col="forestgreen", pch=20)

By_try <- MCMCchains(post,"By")
envplot(By_try,1980:2022,n=1,ylab="Bycatch (t)",xlab="Year")
#qDm<-quantile(pullpost(Simple,"q.D"),c(0.5))
points(Y,By.obs, col="red", pch=18)
points(Y,Hal*wcpue,col="purple", pch=18)
points(Y,Hal*(mean(wcpue, na.rm=T)),col="seagreen", pch=20)
#points(Y,KnC.obs+D.surv/qDm, col="seagreen", pch=12)

D_try <- MCMCchains(post,"D")
envplot(D_try,1980:2022,n=1,ylab="Discards (t)",xlab="Year", cols=c("black","black"))
points(Y,D.obs,col="red", pch=18)
points(Y,Hal*wcpue-Lnd.By, col="violet", pch=16)
points(Y,Hal*(mean(wcpue, na.rm=T))-Lnd.By, col="orange", pch=17)

KnC_try <- MCMCchains(post,"KnC")
envplot(KnC_try,1980:2022,cols=c("orange","orange"), n=1,ylab="Knownw Catch (t)",xlab="Year")
points(Y,KnC.obs,col="darkcyan", pch=18)

envplot(C_try,1980:2022,n=1,ylab="Total Catch (t)",xlab="Year", cols=c("grey","grey"))
addenv(By_try, cols=c("blue","blue"))
lines(Lnd.By~Y, col="royalblue", lwd=2, lty=2)
addenv(D_try, cols=c("black","black"))
addenv(KnC_try, cols=c("orange","orange"))

CtB_try <- MCMCchains(post,"CtoB")
envplot(CtB_try,1980:2022,n=1,ylab="C to B",xlab="Year",cols=c("red","red"))

S_try<-MCMCchains(post,"S")
envplot(S_try,1980:2022,n=1,ylab="Surplus",xlab="Year")

Surplus_try<-MCMCchains(post,"Surplus")
envplot(Surplus_try,1980:2022,n=1,ylab="Surplus",xlab="Year",cols=c("darkorange","darkorange"))

B1_try<-MCMCchains(post,"B1")
envplot(B1_try,1980:2022,n=1,ylab="Surplus",xlab="Year",cols=c("darkorange","darkorange"))

B2_try<-MCMCchains(post,"B2")
envplot(B2_try,1980:2022,n=1,ylab="Surplus",xlab="Year",cols=c("darkorange","darkorange"))

#Check Tau's
t1<-pullpost(Simple,"Tau1")
quantile(t1,c(0.05,0.1,0.5,0.90,0.95))
t2<-pullpost(Simple,"Tau2")
quantile(t2,c(0.05,0.1,0.5,0.90,0.95))


SB<-lm(S_try[2:41] ~ Bhist[2:41])
length(S_try)
length(B_try)
Bhist<-B_try[2:41]

str(S_try)
str(Bhist)

S_try[1]   #S estimates from 1st year
