###############################################################################
## development of r prior for yelloweye using McAllister et al. 2001 methods
## and life history from Love et al.
###############################################################################
{library(MCMCvis)
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
library(fishmethods)
library(FSA)}

#================================================================================
# LWA data for yelloweye
Port1<-read.csv("Data/port_sampling_bio_data.csv")
Port2<-read.csv("Data/port_sampling_bio_data2.csv")
Port<-rbind(Port1,Port2)
colnames(Port)  #Port <- subset(Port,!is.na(Weight.Kilograms) & !is.na(Length.Millimeters))

Fem<-Port[Port$Sex == "Female",]
L<-mean(Fem$Length.Millimeters, na.rm=T)

#=================================================================================
# natural mortality prior dev.
getmode <- function(v,r) {
  rnds<-round(v,r)
  uniqv <- unique(rnds)
  uniqv[which.max(tabulate(match(rnds, uniqv)))]
}
lognorms<-function(mean,sd){
  rphypln<-rlnorm(5000000,log(mean),sd)
  hist(rphypln[rphypln<0.2],
       breaks=100,col="orange")
  abline(v=mean(rphypln), col="blue", lwd=1.5)
  abline(v=getmode(rphypln), col="green")
  abline(v=median(rphypln), col="red")
  print("mean");print(mean(rphypln))
  print("mode");print(getmode(rphypln,3))
  print("median");print(median(rphypln))
  print("variance"); print(var(rphypln))
}
par(mfrow=c(3,1))
lognorms(0.02,0.3)

gammas<-function(shape,rate){
  rphypgam<-rgamma(5000000,shape=shape,rate=rate)
  hist(rphypgam,breaks=100,col="orange")
  abline(v=mean(rphypgam), col="blue", lwd=1.5)
  abline(v=getmode(rphypgam), col="green")
  abline(v=median(rphypgam), col="red")
  print("mean");print(mean(rphypgam))
  print("mode");print(getmode(rphypgam))
  print("median");print(median(rphypgam))
  print("variance"); print(var(rphypgam))
}
#================================================================================
# fecundity at length from Dick et al. 2017 (lengths in mm)
# should be between 1.2 -2.7 million based on Love et al.
#Sebastes, unobserved subgenus... 
feca<-6.538e-06  #(1.63 e- 11, 1.03 e- 04)
fecb<-4.043   #(3.43, 4.71)

adist<-rlnorm(100000,log(6.538e-06),1.4); quantile(adist,c(0.025,0.5,0.975)); hist(adist,breaks=100)
mean(adist)

adist<-rlnorm(100000,-11.938,2); quantile(adist,c(0.025,0.5,0.975)); hist(adist,breaks=100)

bdist<-rlnorm(100000,log(4.043),0.081); quantile(bdist,c(0.025,0.5,0.975)); hist(bdist,breaks=100)
mean(bdist)
#subgenus Sebastodes(unobserved)
#feca<-1.566e-05  #-11.065
#fecb<-3.778
#exp(log(feca)+fecb*log(L))

#subgenus Sebastes
#feca<-3.256e-07  #-11.065
#fecb<-4.423
#exp(log(feca)+fecb*log(L))
F<-feca*L^fecb 

F<-exp(log(feca)+fecb*log(L))
#================================================================================
## get some LVB parameters  #Port <- subset(Port,!is.na(Weight.Kilograms) & !is.na(Length.Millimeters))
#get females  colnames(Fem)
nrow(Fem)
flage<-subset(Fem,!is.na(Age) & !is.na(Length.Millimeters))
#fwage<-subset(Fem,!is.na(Weight.Kilograms) & !is.na(Age))

laa_st<-vbStarts(data=flage, Length.Millimeters~Age)

lvb_f <- fishmethods::growth(unit = 1, # length (2 = weight)
                                 size = flage$Length.Millimeters, 
                                 age = flage$Age,
                                 error = 1, # additive (2 = multiplicative)
                                 # starting values from Love et al.)
 #                                Sinf = 65.93, K = 0.04, t0 = -11.65)
                                 Sinf = laa_st$Linf, K = laa_st$K, t0 = laa_st$t0)
lvb_f
summary(lvb_f$vout)$parameters #[,1:3]

Linf<-summary(lvb_f$vout)$parameters[1,1]  #65.93
K<-summary(lvb_f$vout)$parameters[2,1] #0.04
t0<-summary(lvb_f$vout)$parameters[3,1] #-11.65

Linf.se<-summary(lvb_f$vout)$parameters[1,2]  #65.93
K.se<-summary(lvb_f$vout)$parameters[2,2] #0.04
t0.se<-summary(lvb_f$vout)$parameters[3,2]

a<-20
ltest<-Linf*(1-exp(-K*(a-t0)))

feca*(ltest^fecb)

par(mfrow=c(3,1))
lognorms(Linf,Linf.se)
ble<-rlnorm(10000,log(Linf),Linf.se); hist(ble, breaks=100)
Linfdis<-rnorm(10000,Linf,Linf.se); hist(Linfdis, breaks=100)
Kdis<-rnorm(10000,K,K.se); hist(Kdis, breaks=100)
t0dis<-rnorm(10000,t0,t0.se); hist(t0dis, breaks=100)

### weight ###================

waa_st_f<-vbStarts(data=fwage, Weight.Kilograms~Age)

wvb_f <- fishmethods::growth(unit = 2, # 1 = length, 2 = weight
                             size = fwage$Weight.Kilograms, 
                             age = fwage$Age,
                             error = 1, # 1 = additive, 2 = multiplicative log(w_i) = log(w_inf) + beta * log(1 - exp * (-k * (age_i - t0))) + error
                             # starting values from Hanselman et al. 2007 (Appendix C, Table 5)
                             Sinf = waa_st_f$Linf, K = waa_st_f$K, t0 = waa_st_f$t0) #,
                             #B = allom_pars %>% filter(Sex == "Female" & Parameter == "b") %>% pull(Estimate))
wvb_f # gompertz failed but that wasn't our target so we're good

Winf<-summary(wvb_f$vout)$parameters[1,1]  #65.93
wK<-summary(wvb_f$vout)$parameters[2,1] #0.04
wt0<-summary(wvb_f$vout)$parameters[3,1] 

Winf.se<-summary(wvb_f$vout)$parameters[1,2]  #65.93
wK.se<-summary(wvb_f$vout)$parameters[2,2] #0.04
wt0.se<-summary(wvb_f$vout)$parameters[3,2]

a<-3
Winf*(1-exp(-wK*(a-wt0)))

Winfdis<-rnorm(10000,Winf,Winf.se); hist(Winfdis, breaks=100)
wKdis<-rnorm(10000,wK,wK.se); hist(wKdis, breaks=100)
wt0dis<-rnorm(10000,wt0,wt0.se); hist(wt0dis, breaks=100)
#==================================================================================
Mod.title<-"r_prior_dev_leslie"
#
cat('model { 
## survival to age 1
  egg.surv<-exp(-Z*365)   #0.95^365
  egg.start<-egg.surv

for (a in 1:A) {
  l[a]<-S^(a)
  l.start[a]<-0.98^(a)
}

#fecundity at age
for (a in 1:A){
  length[a]<-max(Linfp*(1-exp(-Kp*(a-t0p))),1)
  eggs[a]<-ifelse(a < AaM,0,feca*(length[a]^fecb))
  m[a]<-eggs[a]*egg.surv
}

#YEAR 1
for (a in 1:A){
  N[a,1]<-dround(1000*l.start[a],0)
}

for (t in 1:Years){
for (a in 1:A) {
   logbirths[a,t]<-log(max(m[a]*N[a,t]*0.5,1))
   births[a,t]<-exp(logbirths[a,t])   
}
}

for (t in 2:(Years)){
  N[1,t]<-max(sum(births[,t-1]),2)  
}

for (t in 2:(Years)){
  for (a in 2:A){
    N[a,t]<-max(S*N[a-1,t-1],1)
  }
}

for (t in 1:Years){
  P[t]<-max(sum(N[,t]),1)
}

for (a in 1:A){
  dnum[a]<-(N[a,Years]/P[Years])-(N[a,Years-1]/P[Years-1])
  dden[a]<-N[a,Years-1]/P[Years-1]
  d[a]<-dnum[a]/dden[a]
}

delta<-100*(1/A)*sum(d)

r.check<-log(P[(Years-10)]/P[(Years-11)])

r<-log(P[Years]/P[Years-1])

#priors
S <- exp(-M)  #1-M     #based on 0.98 or S=exp(-M)
#M ~ dlnorm(log(mean.M),sigma.M)T(0.001,0.08) #dunif()
tau.M<-1/log(sigma.M*sigma.M+1)
M ~ dlnorm(log(mean.M),tau.M)T(0.001,0.08) #dunif()

#AaM ~ dunif(14,22)
tau.AaM<-1/log(sigma.AaM*sigma.AaM+1)
AaM ~ dnorm(mean.AaM,tau.AaM)T(12,22)

Z <-mean.Z #~ dnorm(mean.Z,sigma.Z)T(0.005,0.1)

feca ~ dlnorm(log(6.538e-06),1/log(1.4*1.4+1))T(1.63e-11,1.03e-04)
fecb ~ dlnorm(log(4.043),1/log(0.081*0.081+1))T(3.43,4.71) #<- 4.043 #

Linfp ~ dnorm(Linf,1/log(Linf.se*Linf.se+1))T(645,653)
Kp ~ dnorm(K,1/log(K.se*K.se+1))T(0.05,0.55)
t0p ~ dnorm(t0,1/(t0.se*t0.se+1))T(-6,-4)

}', file=Mod.title) 

#zs<-rnorm(100000,0.05,0.02); hist(zs)
#================================================================================
# load and run... 
Years<-500
A<-120

hist(rlnorm(10000,log(0.02),0.3))
mean.M<-0.02
sigma.M<-0.3

hist(rlnorm(10000,log(0.05),0.05))
mean.Z<-0.05
sigma.Z<-0.05

hist(rnorm(10000,17,2))
mean.AaM<-17
sigma.AaM<-2
  
data<-list(Years=Years,A=A,
           Linf=Linf,K=K,t0=t0,
           Linf.se=Linf.se,K.se=K.se,t0.se=t0.se,
           mean.M=mean.M,sigma.M=sigma.M,
           mean.Z=mean.Z, #sigma.Z=sigma.Z,
           mean.AaM=mean.AaM,sigma.AaM=sigma.AaM) #,
           
parameters=c("r","r.check",
             "delta","AaM",
             "M","S","Z",#,
             "P" #,"N"#"births",
             )

{
  tstart <- Sys.time()
  print(tstart)
  niter <- 1000 #1000000
  set.seed(2054)
  post <- jagsUI::jags(model.file=Mod.title, data=data,
                       parameters.to.save=parameters, #inits=inits,
                       n.chains=1, parallel=F, n.iter=niter,
                       n.burnin=niter/4, n.thin=niter/500)  # can mess with n.adapt=
  print(Sys.time() - tstart)
  dir.create(paste("Model Output/",Mod.title,"_",niter/1000,"k",sep=""))
  save(post,file=paste("Model Output/", Mod.title,"_",niter/1000,"k/post.Rdata", sep=""))
  
  #print(Sys.time() - tstart)
}

#500 years and 1000 iters ~ 7 minutes
#2000 years and 2000 iters ~ blows up - can't calculate r; do P's get fucked up bc pop gets too big? 

#================================================================================
# look at results

tbl<-jags.View(post, title="", digits=3)

MCMCtrace(post, params=c("M","S","r","AaM","delta"),
          ISB=TRUE, pdf=FALSE, Rhat=FALSE, file="Model Output/Param_trace.pdf")
MCMCtrace(post, params=c("P[10]"),
          ISB=FALSE, pdf=FALSE, Rhat=TRUE, file="Model Output/Param_trace.pdf")
MCMCtrace(post, params=c("P[100]"),
          ISB=FALSE, pdf=FALSE, Rhat=TRUE, file="Model Output/Param_trace.pdf")
MCMCtrace(post, params=c("P[490]"),
          ISB=FALSE, pdf=FALSE, Rhat=TRUE, file="Model Output/Param_trace.pdf")

ps100<-MCMCchains(  object=post,  params = "P[200]",  excl = NULL,
                 ISB = FALSE,  #ignore square brackets
                 mcmc.list = FALSE,chain_num = NULL
); hist(log(ps100),breaks=100)

quantile(ps100,c(0.5,0.95,0.99,0.999))

rs<-MCMCchains(  object=post,  params = "r",  excl = NULL,
                    ISB = TRUE,  #ignore square brackets
                    mcmc.list = FALSE,chain_num = NULL
)
rs<-rs[rs>0]
hist(rs,breaks=100); 

{mean(rs)
median(rs)
getmode(rs,2)}

{plot(density(rs))
abline(v=mean(rs), col="red")
abline(v=median(rs), col="blue")
abline(v=getmode(rs,2), col="violet")}

#MCMCplot(post,params=c('r'))
#===============================================================================
# Find distribution that matches the posterior for rB.stable.... 
plotr<-function(){plot(density(rs),lwd=2, 
                       main=paste0(" M = ",m," Z = ",z," AaM = ", as))
abline(v=mean(rs), col="darkcyan",lwd=2)
abline(v=median(rs), col="blue",lwd=2)
abline(v=getmode(rs,2), col="violet",lwd=2)
mtext(round(mean(rs),3), col="darkcyan",side=3,padj=2,adj=0.75)
mtext(round(median(rs),3), col="blue",side=3,padj=3,adj=0.75)
mtext(round(getmode(rs,3),2), col="violet",side=3,padj=4,adj=0.75)
}

prior.ln<-rlnorm(length(rs),log(0.025),1.8); lines(density(prior.ln),col="red")
for (i in 1:100) {
  prior.ln<-rlnorm(length(rs),log(0.03),2) 
  lines(density(prior.ln),
        col=adjustcolor("red", alpha.f=.15))}

shape=myshape
rate=1/myscale
{ plotr()
  for (i in 1:100) {
  prior.gam<-rgamma(length(rs),shape,rate) 
  lines(density(prior.gam),
        col=adjustcolor("green", alpha.f=.15))
  abline(v=mean(prior.gam),col=adjustcolor("darkcyan", alpha.f=.15),lty=2)
  abline(v=median(prior.gam),col=adjustcolor("blue", alpha.f=.15),lty=2)
  abline(v=getmode(prior.gam,3),col=adjustcolor("violet", alpha.f=.15),lty=2)}
mtext(round(mean(rgamma(10000,shape,rate)),3), col="darkcyan",side=3,padj=2,adj=0.95)
mtext(round(median(rgamma(10000,shape,rate)),3), col="blue",side=3,padj=3,adj=0.95)
mtext(round(getmode(rgamma(10000,shape,rate),3),3), col="violet",side=3,padj=4,adj=0.95)}

B1=betas$alpha
B2=betas$beta
{ plotr()
  for (i in 1:100) {
    prior.beta<-rbeta(length(rs),B1,B2) 
    lines(density(prior.beta),
          col=adjustcolor("green", alpha.f=.15))
    abline(v=mean(prior.beta),col=adjustcolor("darkcyan", alpha.f=.15),lty=2)
    abline(v=median(prior.beta),col=adjustcolor("blue", alpha.f=.15),lty=2)
    abline(v=getmode(prior.beta,3),col=adjustcolor("violet", alpha.f=.15),lty=2)}
  mtext(round(mean(rbeta(10000,shape,rate)),3), col="darkcyan",side=3,padj=2,adj=0.95)
  mtext(round(median(rbeta(10000,shape,rate)),3), col="blue",side=3,padj=3,adj=0.95)
  mtext(round(getmode(rbeta(10000,shape,rate),3),3), col="violet",side=3,padj=4,adj=0.95)}

#********************************************************************************
#*Run a bunch of models to see how they effect our prior...
#priortry[priortry<0.]
Ms<-c(0.02) #c(0.02,0.03,0.04,0.05)
Zs<-c(0.04)#c(0.045,0.05,0.055)#c(0.04,0.05,0.06)
AaMs<-c(15,17,20)
Years<-500
A<-120

R.out<-data.frame()
niter<-100000
dir.create(paste("Model Output/",Mod.title,"_",niter/1000,"k_r_search4",sep=""))

i<-1
for (m in Ms) {  #m<-Ms[1]
  for (z in Zs){ #z<-Zs[1]
    for (as in AaMs){  #as<-AaM[1]
      
      #if(i< 4) {i<-i+1} else {  #as<-AaMs[3]
      
      mean.M<-m   #mean.M<-0.03
      sigma.M<-0.3
      mean.Z<-z
      mean.AaM<-as
      sigma.AaM<-2
      
      data<-list(Years=Years,A=A,
                 Linf=Linf,K=K,t0=t0,
                 Linf.se=Linf.se,K.se=K.se,t0.se=t0.se,
                 mean.M=mean.M,sigma.M=sigma.M,
                 mean.Z=mean.Z,#sigma.Z=sigma.Z,
                 mean.AaM=mean.AaM,sigma.AaM=sigma.AaM) #,
      
      parameters=c("r","r.check",
                   "delta","AaM",
                   "M","S","Z",#,
                   "P")
        tstart <- Sys.time()
        print(tstart)
        #niter <- 5000 #1000000
        set.seed(2054)
        post <- jagsUI::jags(model.file=Mod.title, data=data,
                             parameters.to.save=parameters, #inits=inits,
                             n.chains=2, parallel=T, n.iter=niter,
                             n.burnin=niter/4, n.thin=niter/500)  # can mess with n.adapt=
        print(Sys.time() - tstart)
        #
        save(post,file=paste("Model Output/", Mod.title,"_",niter/1000,
                             "k/post_",
                             "M_",m ,"_Z_",z,"_AaM_",as,
                             ".Rdata", sep=""))
        
        rs<-MCMCchains(  object=post,  params = "r",  excl = NULL,
                         ISB = TRUE,  mcmc.list = FALSE,chain_num = NULL)
        rs<-rs[rs>0]
        
        R.out[i,"mean.M"]<-m
        R.out[i,"sd.M"]<-sigma.M
        R.out[i,"mean.Z"]<-z
        R.out[i,"sd.Z"]<-NA
        R.out[i,"mean.AaM"]<-as
        R.out[i,"sd.AaM"]<-sigma.AaM
        R.out[i,"mean.r"]<-mean(rs)
        R.out[i,"var.r"]<-var(rs)
        R.out[i,"median.r"]<-median(rs)
        R.out[i,"mode.r"]<-getmode(rs,2)
        R.out[i,"lo95.r"]<-quantile(rs,c(0.025))
        R.out[i,"hi95.r"]<-quantile(rs,c(0.975))
        #search for gamma distribution
        rs95<-quantile(rs,c(0.95))
        gammafind <- function(shape) {
          scale <- mean(rs)/shape
          pgamma(rs95, shape, scale=scale) - 0.95
        }
        
        tmp <- uniroot( gammafind, lower=0.1, upper=100 )
        
        R.out[i,"gam.shape"]<- tmp$root
        R.out[i,"gam.scale"]<-1/(mean(rs)/tmp$root)
        #myscale <- mean(rs)/tmp$root
        
        estBetaParams <- function(mu, var) {
          alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
          beta <- alpha * (1 / mu - 1)
          return(params = list(alpha = alpha, beta = beta))
        }
        
        betas<-estBetaParams(mean(rs),var(rs))
        
        R.out[i,"beta.B1"]<- betas$alpha
        R.out[i,"beta.B2"]<-betas$beta
        i<-i+1
        write.csv(R.out,file=paste("Model Output/", Mod.title,"_",niter/1000,"k_r_search4.csv", sep=""))
     # }#skipping one AaM because I fucked up
    }
  }
}

oldr<-read.csv("Model Output/r_prior_dev_leslie_100k/r_search1.csv")
oldr<-subset(oldr,select=-1)
ncol(R.out)
R.out<-rbind(oldr,R.out)
R.out<-unique(Rss)
write.csv(R.out,file=paste("Model Output/", Mod.title,"_",niter/1000,"k/r_search.csv", sep=""))

#++++++++++++++++****************************+++++++++++++++++
#+Quick check... 
R.out<-read.csv("Model Output/r_prior_dev_leslie_100k_r_search_07302022.csv")
R.out<-R.out[order(R.out$mean.Z),]

hist(R.out$mean.r, breaks=25)
hist(R.out$median.r, breaks=5)

par(mfrow=c(2,2))
plot(R.out$median.r~R.out$mean.M)
plot(R.out$median.r~R.out$mean.Z)
plot(R.out$median.r~R.out$mean.AaM)

par(mfrow=c(2,2))
plot(R.out$var.r~R.out$mean.M)
plot(R.out$var.r~R.out$mean.Z)
plot(R.out$var.r~R.out$mean.AaM)

plot(density(R.out$median.r))
medbetas<-estBetaParams(mean(R.out$median.r),var(R.out$median.r))
B1=medbetas$alpha
B2=medbetas$beta
{ #plotr()
  for (i in 1:100) {
    prior.beta<-rbeta(length(R.out$median.r),B1,B2) 
    lines(density(prior.beta),
          col=adjustcolor("purple", alpha.f=.15))
    abline(v=mean(prior.beta),col=adjustcolor("darkcyan", alpha.f=.15),lty=2)
    abline(v=median(prior.beta),col=adjustcolor("blue", alpha.f=.15),lty=2)
    abline(v=getmode(prior.beta,3),col=adjustcolor("violet", alpha.f=.15),lty=2)}
  mtext(round(mean(rbeta(10000,shape,rate)),3), col="darkcyan",side=3,padj=2,adj=0.95)
  mtext(round(median(rbeta(10000,shape,rate)),3), col="blue",side=3,padj=3,adj=0.95)
  mtext(round(getmode(rbeta(10000,shape,rate),3),3), col="violet",side=3,padj=4,adj=0.95)}

#***********************************************************************************
#*Examine models and compare predicted prior with posterior one by one...  
m<-0.04
z<-0.05
as<-17

load(file=paste("Model Output/", Mod.title,"_",niter/1000,
                     "k/post_",
                     "M_",m ,"_Z_",z,"_AaM_",as,
                     ".Rdata", sep=""))
rs<-MCMCchains(  object=post,  params = "r",  excl = NULL,
                 ISB = TRUE,  #ignore square brackets
                 mcmc.list = FALSE,chain_num = NULL
)
rs<-rs[rs>0]


tbl<-jags.View(post, title="", digits=3)

MCMCtrace(post, params=c("M","S","r","AaM","delta"),
          ISB=TRUE, pdf=FALSE, Rhat=TRUE, file="Model Output/Param_trace.pdf")

par(mfrow=c(2,2))

shape=R.out$gam.shape[R.out$mean.M==m & R.out$mean.Z==z & R.out$mean.AaM==as]
rate=R.out$gam.scale[R.out$mean.M==m & R.out$mean.Z==z & R.out$mean.AaM==as]
{ plotr()
  for (i in 1:100) {
    prior.gam<-rgamma(length(rs),shape,rate) 
    lines(density(prior.gam),
          col=adjustcolor("green", alpha.f=.15))
    abline(v=mean(prior.gam),col=adjustcolor("darkcyan", alpha.f=.15),lty=2)
    abline(v=median(prior.gam),col=adjustcolor("blue", alpha.f=.15),lty=2)
    abline(v=getmode(prior.gam,3),col=adjustcolor("violet", alpha.f=.15),lty=2)}
  mtext(round(mean(rgamma(10000,shape,rate)),3), col="darkcyan",side=3,padj=2,adj=0.95)
  mtext(round(median(rgamma(10000,shape,rate)),3), col="blue",side=3,padj=3,adj=0.95)
  mtext(round(getmode(rgamma(10000,shape,rate),3),3), col="violet",side=3,padj=4,adj=0.95)}

B1=R.out$beta.B1[R.out$mean.M==m & R.out$mean.Z==z & R.out$mean.AaM==AaM]
B2=R.out$beta.B2[R.out$mean.M==m & R.out$mean.Z==z & R.out$mean.AaM==AaM]
{ plotr()
  for (i in 1:100) {
    prior.beta<-rbeta(length(rs),B1,B2) 
    lines(density(prior.beta),
          col=adjustcolor("purple", alpha.f=.15))
    abline(v=mean(prior.beta),col=adjustcolor("darkcyan", alpha.f=.15),lty=2)
    abline(v=median(prior.beta),col=adjustcolor("blue", alpha.f=.15),lty=2)
    abline(v=getmode(prior.beta,3),col=adjustcolor("violet", alpha.f=.15),lty=2)}
  mtext(round(mean(rbeta(10000,shape,rate)),3), col="darkcyan",side=3,padj=2,adj=0.95)
  mtext(round(median(rbeta(10000,shape,rate)),3), col="blue",side=3,padj=3,adj=0.95)
  mtext(round(getmode(rbeta(10000,shape,rate),3),3), col="violet",side=3,padj=4,adj=0.95)}

#********************************************************************************

library(RColorBrewer)
dev.off()
# set the colour palette
r.select<-R.out
r.select<-R.out[R.out$mean.M > 0.015,]
#r.select<-R.out[R.out$hi95.r < 0.2,]
nrow(r.select)
Ms<-unique(r.select$mean.M)
Zs<-unique(r.select$mean.Z)
Zs<-Zs[Zs<0.06]
AaMs<-unique(r.select$mean.AaM)



par(mfrow=c(1,1))
cols <- brewer.pal(nrow(r.select),'Spectral')
fun_color_range <- colorRampPalette(c("#1b98e0", "red"))
cols <- fun_color_range(26) # nrow(r.select)) 
tran<-0.5

i<-1
for (z in Zs) {  #m<-Ms[1]
  for (m in Ms){ #z<-Zs[1]
    for (as in AaMs){  #as<-AaMs[1]
      check<-r.select[r.select$mean.M==m & r.select$mean.Z==z & r.select$mean.AaM==as,]
      if (nrow(check) > 0) {
        B1=r.select$beta.B1[r.select$mean.M==m & r.select$mean.Z==z & r.select$mean.AaM==as]
        B2=r.select$beta.B2[r.select$mean.M==m & r.select$mean.Z==z & r.select$mean.AaM==as]
        prior<-rbeta(10000,B1,B2)
        if (i == 1) {
          plot(density(prior), xlim=c(-0.01,0.3),col=cols[i], ylim=c(0,20))
          mtext(paste("M=",m," Z=",z," A=",as),adj=1,padj=2,side=3,cex=0.5,col=cols[i],lwd=2)
          i<-i+1
        } else {
          lines(density(prior), col=cols[i])
          mtext(paste("M=",m," Z=",z," A=",as),adj=1,padj=1+i,side=3,cex=0.5,col=cols[i],lwd=2)
          i<-i+1
        }
        
      } else{}
    }
  }
}

fun_color_range <- colorRampPalette(c("black", "goldenrod3"))
cols <- fun_color_range(26) # nrow(r.select)) 
tran<-0.5
lwdc<-2
bgc<-"lightgrey"
avgcol<-"darkred"

png(paste("Figures//r_prior_fits.png", sep=""),
    width=8,height=4,#width=9.5,height=8.5,
    units="in",res=1200)

par(mfrow=c(1,4), bg="white")

Zs<-Zs[order(Zs)]
for (z in Zs) {  #m<-Ms[1]
  zplot<-r.select[r.select$mean.Z==z,]
  cols <- fun_color_range(nrow(zplot))
  i<-1
  for (m in Ms){ #z<-Zs[1]
    for (as in AaMs){  #as<-AaMs[1]
      check<-zplot[zplot$mean.M==m & zplot$mean.AaM==as,]
      if (nrow(check) > 0) {
        B1=zplot$beta.B1[zplot$mean.M==m & zplot$mean.Z==z & zplot$mean.AaM==as]
        B2=zplot$beta.B2[zplot$mean.M==m & zplot$mean.Z==z & zplot$mean.AaM==as]
        prior<-rbeta(10000,B1,B2)
        #prior<-seq(0,1,length=100)
        if (i == 1) {
          #plot(prior,dbeta(prior,B1,B2),type="l")
          plot(density(prior), xlim=c(-0.01,0.2),
               col=alpha(cols[i],tran), 
               ylim=c(0,35),
               main=paste("Z = ",z), xlab="r prior")
          
          rect(par("usr")[1], par("usr")[3],
               par("usr")[2], par("usr")[4],
               col = bgc)
          
          # adding the new plot
          par(new = TRUE)
          
          plot(density(prior), xlim=c(-0.01,0.2),
               col=cols[i], 
               ylim=c(0,35), lwd=lwdc,
               main=paste("Z = ",z), xlab="r prior")
          mtext(paste("M=",m," A=",as),adj=1,padj=2,side=3,
                cex=0.75,col=cols[i],lwd=2)
          i<-i+1
        } else {
          #lines(prior,dbeta(prior,B1,B2),type="l")
          lines(density(prior), col=alpha(cols[i],tran+0.015*i),lwd=lwdc)
          mtext(paste("M=",m," A=",as),adj=1,padj=1+i*1.1,side=3,
                cex=0.75,col=cols[i],lwd=2)
          i<-i+1
        }
        
      } else{}
    }
  }
  #B1=mean(zplot$beta.B1)
  B1=quantile(zplot$beta.B1, c(0.5))
  B2=quantile(zplot$beta.B2, c(0.5))
  avgprior<-rbeta(100000,B1,B2)
  lines(density(avgprior),col=avgcol, lwd=2, lty=2)
  mtext(paste("med.beta(",round(B1,3),", ",round(B2,3),")",sep=""),
        adj=0,padj=-0.2,side=3,
        cex=0.75,col=avgcol,lwd=2)
}
dev.off()

#********************************************************************************
#*If you want to plot all the variations separately
#*Takes a while to run this
#*******************************************************************************
#*Loop through and plot them all! 
Ms<-unique(R.out$mean.M)
Zs<-unique(R.out$mean.Z)
AaMs<-unique(R.out$mean.AaM)

pc<-1
for (m in Ms) {  #m<-Ms[1]
  for (z in Zs){ #z<-Zs[1]
    for (as in AaMs){  #as<-AaMs[2]
      mean.M<-m
      sigma.M<-0.3
      mean.Z<-z
      mean.AaM<-as
      sigma.AaM<-2
      
      check<-R.out[R.out$mean.M==m & R.out$mean.Z==z & R.out$mean.AaM==as,]
      if (nrow(check) > 0) {
        
        load(file=paste("Model Output/", Mod.title,"_",niter/1000,
                        "k/post_",
                        "M_",m ,"_Z_",z,"_AaM_",as,
                        ".Rdata", sep=""))
        rs<-MCMCchains(  object=post,  params = "r",  excl = NULL,
                         ISB = TRUE,  #ignore square brackets
                         mcmc.list = FALSE,chain_num = NULL
        )
        rs<-rs[rs>0]
        
        plotno<-seq(from=1,to=nrow(R.out),by=4)
        if(any(pc==plotno)) {
          png(paste("Model Output/", Mod.title,"_",niter/1000,
                    "k/r_prior_fit_",pc,
                    ".png", sep=""),
              width=7,height=6,#width=9.5,height=8.5,
              units="in",res=1200)
          par(mfrow=c(2,2), mar=c(4,5,3,1))
          par(mfrow=c(2,2))
        }
        
        shape=R.out$gam.shape[R.out$mean.M==m & R.out$mean.Z==z & R.out$mean.AaM==as]
        rate=R.out$gam.scale[R.out$mean.M==m & R.out$mean.Z==z & R.out$mean.AaM==as]
        plotr()
        for (i in 1:100) {
          prior.gam<-rgamma(length(rs),shape,rate) 
          lines(density(prior.gam),
                col=adjustcolor("green", alpha.f=.15))
          abline(v=mean(prior.gam),col=adjustcolor("darkcyan", alpha.f=.1),lty=2)
          abline(v=median(prior.gam),col=adjustcolor("blue", alpha.f=.1),lty=2)
          abline(v=getmode(prior.gam,3),col=adjustcolor("violet", alpha.f=.1),lty=2)}
        #mtext(round(mean(rgamma(10000,shape,rate)),3), col="darkcyan",side=3,padj=2,adj=0.95)
        #mtext(round(median(rgamma(10000,shape,rate)),3), col="blue",side=3,padj=3,adj=0.95)
        #mtext(round(getmode(rgamma(10000,shape,rate),3),3), col="violet",side=3,padj=4,adj=0.95)}
        
        B1=R.out$beta.B1[R.out$mean.M==m & R.out$mean.Z==z & R.out$mean.AaM==as]
        B2=R.out$beta.B2[R.out$mean.M==m & R.out$mean.Z==z & R.out$mean.AaM==as]
        { #plotr()
          for (i in 1:100) {
            prior.beta<-rbeta(length(rs),B1,B2) 
            lines(density(prior.beta),
                  col=adjustcolor("purple", alpha.f=.15))
            abline(v=mean(prior.beta),col=adjustcolor("darkcyan", alpha.f=.1),lty=2)
            abline(v=median(prior.beta),col=adjustcolor("blue", alpha.f=.1),lty=2)
            abline(v=getmode(prior.beta,3),col=adjustcolor("violet", alpha.f=.1),lty=2)}
          mtext(round(mean(rbeta(10000,shape,rate)),3), col="darkcyan",side=3,padj=2,adj=0.95)
          mtext(round(median(rbeta(10000,shape,rate)),3), col="blue",side=3,padj=3,adj=0.95)
          mtext(round(getmode(rbeta(10000,shape,rate),3),3), col="violet",side=3,padj=4,adj=0.95)}
        
        if(any(pc==(plotno-1))){
          dev.off()
        }
        pc<-pc+1} else {}
    }
  }
}
#***********************************************************************************
meanrs<- mean(rs)
rs95<-quantile(rs,c(0.95))

myfun <- function(shape) {
  scale <- meanrs/shape
  pgamma(rs95, shape, scale=scale) - 0.95
}




#===============================================================================
Mod.title<-"r_prior_dev_gentime"
#
cat('model { 
## survival to age x
  #Z<-0.05  #daily mortality of eggs/larva
  egg.surv<-exp(-Z*365)   #0.95^365
  egg.start<-egg.surv

for (a in 1:A) {
  l[a]<-S^(a)
}

for (a in 1:A){
  length[a]<-max(Linfp*(1-exp(-Kp*(a-t0p))),1)
  eggs[a]<-ifelse(a < AaM,0,feca*(length[a]^fecb))
  m[a]<-eggs[a]*egg.surv*0.5
  
  G.num[a]<-a*l[a]*m[a]
  G.den[a]<-l[a]*m[a]
}

#Generation time
G<-sum(G.num)/sum(G.den)

r<-log(sum(G.den))/G

#priors
S <- exp(-M)  #1-M     #based on 0.98 or S=exp(-M)
M ~ dlnorm(log(mean.M),sigma.M)T(0.001,0.08) #dunif()

#AaM ~ dunif(14,22)
AaM ~ dnorm(mean.AaM,sigma.AaM)T(12,22)

Z <-0.05 #~ dlnorm(log(mean.Z),sigma.Z)T(0.005,0.1)

feca ~ dlnorm(log(6.538e-06),1.4)T(1.63e-11,1.03e-04)
fecb ~ dlnorm(log(4.043),0.081)T(3.43,4.71) #<- 4.043 #

Linfp ~ dnorm(Linf,Linf.se)T(645,653)
Kp ~ dnorm(K,K.se)T(0.05,0.55)
t0p ~ dnorm(t0,t0.se)T(-6,-4)

}', file=Mod.title) 

hist(rlnorm(10000,log(0.02),0.3))
mean.M<-0.02
sigma.M<-0.3

hist(rlnorm(10000,log(0.05),0.05))
mean.Z<-0.05
sigma.Z<-0.05

hist(rnorm(10000,17,2))
mean.AaM<-17
sigma.AaM<-2

data<-list(Years=Years,A=A,
           Linf=Linf,K=K,t0=t0,
           Linf.se=Linf.se,K.se=K.se,t0.se=t0.se,
           mean.M=mean.M,sigma.M=sigma.M,
           #mean.Z=mean.Z,sigma.Z=sigma.Z,
           mean.AaM=mean.AaM,sigma.AaM=sigma.AaM) #,
#Winf=Winf,wK=wK,wt0=wt0,
#Winf.se=Winf.se,wK.se=wK.se,wt0.se=wt0.se)
parameters=c("r",#"r.check",
  "G",
  "M","S"#,
  #"P" #,"N"#"births",
)

{
  tstart <- Sys.time()
  print(tstart)
  niter <- 2000 #1000000
  set.seed(2054)
  post <- jagsUI::jags(model.file=Mod.title, data=data,
                       parameters.to.save=parameters, #inits=inits,
                       n.chains=2, parallel=T, n.iter=niter,
                       n.burnin=niter/4, n.thin=niter/1000)  # can mess with n.adapt=
  print(Sys.time() - tstart)
  #dir.create(paste("Model Output/",Mod.title,"_",niter/1000,"k",sep=""))
  #save(post,file=paste("Model Output/", Mod.title,"_",niter/1000,"k/post.Rdata", sep=""))
  
  #print(Sys.time() - tstart)
}






#===============================================================================
#===============================================================================
#===============================================================================
# Population Model Scratch
S<-0.98
## survival to age x
#l[1]<-0.8
l<-vector()
#for (a in 1:A) {
#  l[a]<-S^a  
#}

## survival to age x
l[1]<-0.7
l.start[1]<-0.7

for (a in 1:A) {
  #l[a]<-S^a  
  #or if a in 2:A and juvy survival sucky... 
  l[a]<-l[1]*S^(a-1)
  l.start[a]<-0.98^a
}

m<-vector()
length<-vector()
weight<-vector()

#fecundity at age
for (a in 1:19){
  m[a]<-0
  length[a]<-max(Linf*(1-exp(-K*(a-t0))),1)
  weight[a]<-max(Winf*(1-exp(-wK*(a-wt0))),0.01) #max(0.0287*(length[a]^2.822),0.01)
}

for (a in 20:A){
  length[a]<-Linf*(1-exp(-K*(a-t0)))
  weight[a]<-Winf*(1-exp(-wK*(a-wt0))) #0.0287*(length[a]^2.822)
  #rel.w[a]<-weight[a]/weight[A]
  m[a]<-feca*(length[a]^fecb) #fec*rel.w[a]
}

#YEAR 1
N<-data.frame()
B<-data.frame()

for (a in 1:A){
  N[a,1]<-10000*l[a]
  B[a,1]<-N[a,1]*weight[a]
}

for (t in 2:(Years)){
  for (a in 2:A){
    N[a,t]<-max(S*N[a-1,t-1],1)
    B[a,t]<-N[a,t]*weight[a]
  }
}

births<-data.frame()
for (t in 1:Years){
  for (a in 1:A) {
    births[a,t]<-m[a]*N[a,t]
  }
}

births[20:40,]

for (t in 2:(Years)){
  N[1,t]<-max(sum(births[,t-1]),1000)
  B[1,t]<-N[1,t]*weight[1]
}























