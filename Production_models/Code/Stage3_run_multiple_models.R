###############################################################################
## This is code for setting up and running multiple SS-SPM models
## Use PT_models.R or Fox_models.R to generate list of model files to pull from
###############################################################################
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
  source("Code/2022_DSR_SAFE_models/Phase3/DATAPREP_SPM_1980_PHASE3.R")
  source("Code/2022_DSR_SAFE_models/Phase1/PLOT_SPM80.R")
  source("Code/Posterior_Plotting/YE_SPM_posterior_exams_Func.R")
}

#============================================================
#Organize which models to run, set derby if needed
#ModPHASE3.1<-"PT2i_base_PHASE3_beta-phi"
ModPHASE3.2<-"PT2i_base_PHASE3_norm-phi"

Mod.list<-c(ModPHASE3.2)

Derby.list<-c(0.3, -0.3)#c(0,0.3,-0.3)
DEsdlist<-c(0.1)

#Blist from r prior developed such that

#B1list<-c(1,1.483,1.247,1.241)
#B2list<-c(1,22.908,31.478,53.481)
B1list<-c(1.483)  
B2list<-c(22.908)  

phiB1<-c(2.339) #from Phase2
phiB2<-c(1.22)

#model order of priors is phi(mu,sd); logK(mu,sig)
#1)   lv3_derb0: phi(0.6802432,0.2054276) logK(10.74306,9.652545)
#2)  lv3_derb3p: phi(0.6333498,0.1911159) logK(10.78994,9.600169)
#3)  lv3_derb3o: phi(0.7340861,0.207902) logK(10.61891,9.37281)

#4)  lv5_derb0: phi(0.7420787,0.1796734) logK(10.62973,9.237758)
#5)  lv5_derb3p: phi(0.6726209,0.1690635) logK(10.75576,9.390445)
#6)  lv5_derb3o: phi(0.7702393,0.1795257) logK(10.58545,9.140554)

phimu.list<-c(0.6333498,0.7340861)#c(0.6802432,0.6333498,0.7340861) #c(0.6802432) #c(0.7420787) #c(0.6802432)
phisig.list<-c(0.1911159,0.207902)#c(0.2054276,0.1911159,0.207902) #c(0.2054276) #c(0.1796734) #c(0.2054276)

bigKmu.list<-c(10.78994,10.61891)#c(10.74306,10.78994,10.61891) #c(10.65388) #c(10.6366) 
bigKsigma.list<-c(9.600169,9.37281)#c(0.2940348)

upvarlist<- c(-3) #c(-5)

#===================================================================
#FLAG: good idea to run list on tiny chains (2K) to check for bugs in model code! 
#beta500k not fully converged, 3.5 hours
niter<-1500000 #1000000 #300000 #350000
burnin<-500000 #300000
#5K, 2.5m.
#600K - 4.8 hours; really close to converged for both norm- and beta-phi models
#700K ~ 6.3 hours; not converged but close.  small r only bounded to 0.3.  to 0.2 will help  Probably 1m like phase 1 gets us there...
set.seed(1234)

#!!!: FLAG: adjust initial values and data in loop for specific models!!! 
#m<-2; d<-1; b<-1
for (v in 1:length(upvarlist)){ #v<-1
for (m in 1:length(Mod.list)) {  #m<-1
  for (d in 1:length(Derby.list)){ #d<-1
     for (b in 1:length(B1list)) { #b<-1
  OFLABC<-c(0.0,0.0,0.00,0.0)  #0.02 (ABC) 0.032 (OFL). 0=nofishing
  HarvStrat<-0 #2.58 #1.96 #1.68 #0, 1.15
  Fu<-1
  #Derby.Eff<-Derby.list[d]
  
  Data<-load.data(YEAR=2022,
                  Derby.Eff = Derby.list[d],
                  DEsd=0.1,  #this is CV for derby
                  B1=B1list[b],
                  B2=B2list[b])
  list2env(Data,.GlobalEnv)
  
  phiB1<-phiB1
  phiB2<-phiB2
  
  phimu<-phimu.list[d]
  phisig<-phisig.list[d]
  
  bigKmu<-bigKmu.list[d]
  bigKsigma<-bigKsigma.list[d]
  
  upvar<-upvarlist[v]
  # adjust data as needed for specific model
  #if (m == 1) {
  #B.obs[2:4,15]<-NA
  #cv.B[2:4,15]<-1
  #}
  
  prep<-data.prep.PHASE3()
  list2env(prep,.GlobalEnv)
  
  #adjust initial values for particular models as needed...
  #if (m == 1) {
  #  inits[[1]]$Tau1<-0; inits[[2]]$Tau1<-0; inits[[3]]$Tau1<-0
  #}
  
  tstart <- Sys.time()
  print(tstart)
  post <- jagsUI::jags(model.file=Mod.list[m], data=data,
                       parameters.to.save=parameters, inits=inits,
                       n.chains=3, parallel=T, n.iter=niter,
                       n.burnin=burnin, n.thin=(niter-burnin)/1000)  # can mess with n.adapt=
  print(Sys.time() - tstart)
  
  if (grepl("beta-phi",Mod.list[m],fixed=TRUE) == TRUE) {
    filename<-paste("PHASE3_B1-",round(B1list[b],2),
                    "_B2-",round(B2list[b]),"_",
                    "upv",(upvarlist[v]),"_",
                    "phB1-",round(phiB1,2),"_",
                    "phB2-",round(phiB2,2),"_",
                    "Kmu-",round(bigKmu,1),"_",
                    "Ksig-",round(bigKsigma,1),"_",
                    "derb_",Derby.list[d],"_",
                    niter/1000,"k",sep="")
  } else {
  filename<-paste("PHASE3_B1-",round(B1list[b],2),
                  "_B2-",round(B2list[b]),"_",
                  "upv",(upvarlist[v]),"_",
                  "phmu-",round(phimu,1),"_",
                  "phsig-",round(phiB2,1),"_",
                  "Kmu-",round(bigKmu,1),"_",
                  "Ksig-",round(bigKsigma,1),"_",
                  "derb_",Derby.list[d],"_",
                  niter/1000,"k",sep="")
  }
  
  dir.create(paste("Model Output/",filename,sep=""))
  save(post,file=paste("Model Output/", filename,"/post.Rdata", sep=""))
  #save(post,file=paste("Model Output/", Mod.list[m],"_",niter/1000000,"m/post.Rdata", sep=""))
  save(post,file=paste("D:/Groundfish Biometrics/MCMC_backups/", filename,".Rdata", sep=""))
  dir.create(paste("Figures/",filename,sep=""))
  
  if (grepl("beta-phi",Mod.list[m],fixed=TRUE) == TRUE){
    MCMCtrace(post, params=c("R.hyp" ,"r","Kseo","K","phi",#"qfCPUE","qsCPUE",
                             "Tau1","Tau2","Tau3", "sigma","eta", "rB1","rB2",
                             "bigphi","phi.eta","littlephiB1","littlephiB2","pi"),
              ISB=TRUE, pdf=TRUE, Rhat=TRUE, 
              file=paste("Figures/",filename,"/Param_trace.pdf",sep=""))
  } else {
    MCMCtrace(post, params=c("R.hyp" ,"r","Kseo","K","phi",#"qfCPUE","qsCPUE",
                             "Tau1","Tau2","Tau3", "sigma","eta", "rB1","rB2",
                             "bigphi","phiTau","pi"),
              ISB=TRUE, pdf=TRUE, Rhat=TRUE, 
              file=paste("Figures/",filename,"/Param_trace.pdf",sep=""))
  }
  
  ggmcmc(ggs(post$samples),family="r",file=paste("Model Output/",filename,"/DIAG_r.pdf", sep=""))
  ggmcmc(ggs(post$samples),family="K",file=paste("Model Output/",filename,"/DIAG_K.pdf", sep=""))
  ggmcmc(ggs(post$samples),family="phi",file=paste("Model Output/",filename,"/DIAG_phi.pdf", sep=""))
  ggmcmc(ggs(post$samples),family="sigma",file=paste("Model Output/",filename,"/DIAG_sigma.pdf", sep=""))
 # ggmcmc(ggs(post$samples),family="Tau1",file=paste("Model Output/",filename,"/DIAG_tau1.pdf", sep=""))
#  ggmcmc(ggs(post$samples),family="Tau3",file=paste("Model Output/",filename,"/DIAG_tau3.pdf", sep=""))
  print(Sys.time() - tstart)
  
    #load(file=paste("Model Output/",Mod1,"_",niter/1000000,"m/post.Rdata", sep=""))
  
  #Function for making plots...
    #fix broken fucking plots in function that I cant figure out why they are fucking up the scale!!
  #1) Biomass observed versus modeled 
  B.ey<-par.ext(par="B",years=N+1,areai=1)
  B.ns<-par.ext("B",N+1,2)
  B.cs<-par.ext("B",N+1,3)
  B.ss<-par.ext("B",N+1,4)
  
  All.years<-seq(min(Years),max(Years)+1,by=1)
  FuYear<-Years+1
  sdlist<-c("EYKT", "NSEO", "CSEO", "SSEO")
  
  colch<-c("blue","red","purple","orange","forestgreen")
  vars<-list(B.ey,B.ns,B.cs,B.ss)
  
  png(paste("Figures/",filename,"/Biomass_fit_FIX.png",sep=""),
      width=7,height=6,#width=9.5,height=8.5,
      units="in",res=1200)
  par(mfrow=c(2,2), mar=c(4,5,3,1))
  for (i in 1:4){  #i<-4
    envplot(vars[[i]],All.years,cols=c(colch[1],colch[1]),n=0,ylab="Biomass (t)",xlab="Year",
            ylim=c(0,17500),main=sdlist[i], cex.axis=0.8, cex.lab=1.25, cex.main=1.5)
    addpoints(post,Years=Years+0.1 ,Points=sCPUE[i,], errordat=cv.sCPUE[i,], scaler="qsCPUE", #scaler="qfCPUE" for example...
              bar=1.95 ,col=colch[4], pch=17, cex=0.75)
    addpoints(post,Years=Years ,Points=B.obs[i,], errordat=cv.B[i,], scaler=0, #scaler="qfCPUE" for example...
              bar=1.95 ,col=colch[2], pch=18, cex=1.2)
    #abline(v=2022,  col=colch[3], lty=2)
    if (i == 2){
      legend(x="topright", c("posterior est.","observed ROV/sub","IPHC CPUE"), col=c(colch[1],colch[2],colch[4]),
             text.col=c(colch[1],colch[2],colch[4]),border=F, bty="y", pch=c(NA,18,17), cex=1)
    } else {}
  }
  dev.off()
  
  N.subd<-N; Years.cont<-Years
  plot_spm(post, Mod.title=Mod.list[m], filename=filename)

}}}}

  dev.off()

inits[[1]]$Tau1





