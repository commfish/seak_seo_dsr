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
  library(IDPmisc)
  library(scales)
  library(tidyverse)
  library(ggridges)
  
  Year <-2023
  source("Production_models/Code/SPM_helper.R")
}

#============================================================
#Organize which models to run
#ModPHASE3.1<-"PT2i_base_PHASE3_beta-phi"
ModPHASE3.2<-"PT2i_base_PHASE3_norm-phi"

phase3<-ModPHASE3.2
simMod<-"PHASE_3_sim"

#Mod.list<-c(ModPHASE3.1,ModPHASE3.2)

#pick which model run you want to simulate:
res.to.sim<-"PHASE3_B1-1_B2-1_upv-5_phmu-0.7_phsig-1.2_Kmu-10.6_Ksig-9.2_derb_0_1500k"

#set the appropriate settings for hyper priors
Derby.list<-c(0)#c(0,0.3,-0.3)
DEsdlist<-c(0.1)

#B1list<-c(1,1.483,1.247,1.241,1.374)
#B2list<-c(1,22.908,31.478,53.481,106.057)
B1list<-c(1)  
B2list<-c(1)  

phiB1<-c(2.339) #from Phase2
phiB2<-c(1.22)

phimu.list<-c(0.7420787,0.6726209,0.7702393) #c(0.6802432) #c(0.7420787) #c(0.6802432)
phisig.list<-c(0.1796734,0.1690635,0.1795257) #c(0.2054276) #c(0.1796734) #c(0.2054276)

bigKmu.list<-c(10.62973,10.75576,10.58545) #c(10.65388) #c(10.6366) 
bigKsigma.list<-c(9.237758,9.390445,9.140554) #c(0.2812395) #c(0.2940348)

upvarlist<- c(-5)

#===================================================================
#FLAG: good idea to run list on tiny chains (2K) to check for bugs in model code! 
#beta500k not fully converged, 3.5 hours
nsims<-50
niter<-1500000#300000 #350000
burnin<-500000
#5K, 2.5m.
#600K - 4.8 hours; really close to converged for both norm- and beta-phi models
#700K ~ 6.3 hours; not converged but close.  small r only bounded to 0.3.  to 0.2 will help  Probably 1m like phase 1 gets us there...
set.seed(1234)

#!!!: FLAG: adjust initial values and data in loop for specific models!!! 
for (i in 1:nsims){  #i<-1
  ests<-read.csv(paste0("Model Output/",res.to.sim,"/results_sum.csv"))
  ests<-as.data.frame(ests)
  
  r<-ests$X50[grep(c("r\\["),ests$parameter)]
  logqs<-log(ests$X50[grepl("qsCPUE",ests$parameter)==T])
  Tau1<-ests$mean[grepl("Tau1",ests$parameter)==T]
  Tau3<-ests$mean[grepl("Tau3",ests$parameter)==T]
  pi<-ests$mean[grepl("pi",ests$parameter)==T]
  phi<-ests$mean[grep(c("phi\\["),ests$parameter)]
  logKseo<-log(ests$X50[ests$parameter=="Kseo"])
  logvar<-log(ests$mean[ests$parameter=="sigma"])
  #bigKsigma<-c(0.269)   #c(0.269)
  Data<-load.data(YEAR=2022,
                  Derby.Eff = Derby.list[1],
                  DEsd=0.1,  #this is CV for derby
                  B1=B1list[1],
                  B2=B2list[1])
  list2env(Data,.GlobalEnv)
  #Set future management strategy
  OFLABC<-c(0,0,0,0)  #0.02 (ABC) 0.032 (OFL). 0=nofishing
  HarvStrat<-0 #2.58 #1.96 #1.68 #0, 1.15
  Fu<-0
  
  data<-list(N=N,Subd=Subd,
             realC=KnC.obs, 
             cv.KnC=cv.KnC,
             realBy = ExpByc, 
             cv.ExpByc = cv.ExpByc,
             Lnd.By=Lnd.By, 
             cv.B=cv.B, cv.sCPUE=cv.sCPUE, #,
             Fu=Fu, OFLABC=OFLABC, #B.err=B.err,
             Derby.Eff = Derby.Eff, DEsd=DEsd,
             r = r,
             logqs=logqs,
             Tau1=Tau1,
             Tau3=Tau3,
             pi=pi,
             phi=phi,
             logKseo=logKseo,
             logvar=logvar)
  realC <- KnC.obs
  realBy <- ExpByc
  
  parameters=c("B.obs","sCPUE" ,"KnC.obs","ExpByc") 
  
  sim <- jagsUI::jags(model.file=simMod, data=data,
                      parameters.to.save=parameters, #inits=inits,
                      n.chains=1, parallel=F, n.iter=1,
                      n.burnin=0)
  save(sim,file=paste("Model Output/", res.to.sim,"/simdat.",i,"post.Rdata", sep=""))

  #load(file=paste("Model Output/",res.to.sim,"/simdat.",i,"post.Rdata", sep=""))
  
  Simulated <- coda::as.mcmc(sim)
  
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
  
  #-----------------------------------------------------------------------------
  #run model with simulated data
  Fu<-1
  Data<-list(Subd=Subd,N=N,Years=Years,
             YEAR=YEAR,
             KnC.obs=sim.KnC,
             cv.KnC=cv.KnC,
             Lnd.By=Lnd.By,
             ExpByc=sim.eby,
             cv.ExpByc=cv.ExpByc,
             B.obs=sim.B,
             cv.B=cv.B,
             sCPUE=sim.sCPUE,
             cv.sCPUE=cv.sCPUE,#fCPUE=fCPUE,cv.fCPUE=cv.fCPUE,
             Derby.Eff = Derby.list[1], DEsd=DEsd,
             B1=B1list[1],B2=B2list[1], Fu=Fu)
  
  KnC.obs=sim.KnC
  ExpByc=sim.eby
  B.obs=sim.B
  sCPUE=sim.sCPUE
  
  phiB1<-phiB1
  phiB2<-phiB2
  
  phimu<-phimu.list[1]
  phisig<-phisig.list[1]
  
  bigKmu<-bigKmu.list[1]
  bigKsigma<-bigKsigma.list[1]
  
  upvar<-upvarlist[1]
  Fu<-1
  
  prep<-data.prep.PHASE3()
  list2env(prep,.GlobalEnv)
  
  tstart <- Sys.time()
  print(tstart)
  post <- jagsUI::jags(model.file=phase3, data=data,
                       parameters.to.save=parameters, inits=inits,
                       n.chains=3, parallel=T, n.iter=niter,
                       n.burnin=burnin, n.thin=(niter-burnin)/1000)  # can mess with n.adapt=
  print(Sys.time() - tstart)
  
  save(post,file=paste("Model Output/", res.to.sim,"/postsim",i,".Rdata", sep=""))
  
  #save(post,file=paste("Model Output/", Mod.list[m],"_",niter/1000000,"m/post.Rdata", sep=""))
  #save(post,file=paste("D:/Groundfish Biometrics/MCMC_backups/", filename,".Rdata", sep=""))
  load(file=paste("Model Output/",res.to.sim,"/postsim",i,".Rdata", sep=""))
  
  dir.create(paste("Figures/",res.to.sim,"/sim",i,sep=""))
  dir.create(paste("Model Output/",res.to.sim,"/sim",i,sep=""))
  filename<-paste(res.to.sim,"/sim",i,sep="")
  
  if (grepl("beta-phi",res.to.sim,fixed=TRUE) == TRUE){
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
  
  B.ey<-par.ext(par="B",years=N+Fu,areai=1)
  B.ns<-par.ext("B",N+Fu,2)
  B.cs<-par.ext("B",N+Fu,3)
  B.ss<-par.ext("B",N+Fu,4)
  
  All.years<-seq(min(Years),max(Years)+Fu,by=1)
  FuYear<-Years+Fu
  sdlist<-c("EYKT", "NSEO", "CSEO", "SSEO")
  
  colch<-c("blue","red","purple","orange","forestgreen")
  vars<-list(B.ey,B.ns,B.cs,B.ss)
  
  #png(paste("Figures/",filename,"/Biomass_fit_FIX.png",sep=""),
  #    width=7,height=6,#width=9.5,height=8.5,
  #    units="in",res=1200)
  par(mfrow=c(2,2), mar=c(4,5,3,1))
  for (j in 1:4){  #j<-4
    envplot(vars[[j]],All.years,cols=c(colch[1],colch[1]),n=0,ylab="Biomass (t)",xlab="Year",
            ylim=c(0,17500),main=sdlist[j], cex.axis=0.8, cex.lab=1.25, cex.main=1.5)
    addpoints(post,Years=Years+0.1 ,Points=sCPUE[j,], errordat=cv.sCPUE[j,], scaler="qsCPUE", #scaler="qfCPUE" for example...
              bar=1.95 ,col=colch[4], pch=17, cex=0.75)
    addpoints(post,Years=Years ,Points=B.obs[j,], errordat=cv.B[j,], scaler=0, #scaler="qfCPUE" for example...
              bar=1.95 ,col=colch[2], pch=18, cex=1.2)
    abline(v=2022,  col=colch[3], lty=2)
    if (j == 2){
      legend(x="topright", c("posterior est.","observed ROV/sub","IPHC CPUE"), col=c(colch[1],colch[2],colch[4]),
             text.col=c(colch[1],colch[2],colch[4]),border=F, bty="y", pch=c(NA,18,17), cex=1)
    } else {}
  }
  #dev.off()
  
  N.subd<-N; Years.cont<-Years
  plot_spm(post, Mod.title=phase3, filename=filename)
  
} 
  
#--------------------------------------------------------------------------------
  if (i == max(nsims)) {
    load(file=paste("Model Output/",res.to.sim,"/post.Rdata", sep=""))
    orig<-post
    #Mods<-list()
    #Mods[[1]]<-orig # str(Mods)
    for (i in 1:nsims){  #i<-2
      load(file=paste("Model Output/", res.to.sim,"/postsim",i,".Rdata", sep=""))
      assign(paste0("postsim",i),post)
      #postsim1
      #Mods[[i+1]]<-c(paste0("postsim",i) = paste0("postsim",i))
    }
  }

#----------------------------------------------------------------------------
# compare simulations to original model!!! 

    Mods<-list(orig=orig, 
               postsim1=postsim1, 
               postsim2=postsim2)  #Mods[1]
    Modnames<-c("True Data","sim. data 1","sim. data 2")
    
    paramcomp(SEOq="FALSE",Mods=Mods,Modnames=Modnames,param="K", 100,75000,10000,
              pscale=2,xlabel="K (t)")
    ggsave(paste0("Figures/",res.to.sim,"/simK_comp_subs.png"), 
           dpi=300, height=8, width=8, units="in")
    
    paramcomp(SEOq="TRUE",Mods=Mods,Modnames=Modnames,param="Kseo", 100,250000,10000,
              pscale=1.25,xlabel="K (t)")
    ggsave(paste0("Figures/",res.to.sim,"/simK_comp.png"), 
           dpi=300, height=8, width=8, units="in")
    
    paramcomp(SEOq="FALSE",Mods=Mods,Modnames=Modnames,param="r", -0.005,0.2,10000,
              pscale=2,xlabel="r")
    ggsave(paste0("Figures/",res.to.sim,"/sim_r_comp_subs.png"), 
           dpi=300, height=8, width=8, units="in")
    
    paramcomp(SEOq="TRUE",Mods=Mods,Modnames=Modnames,param="R.hyp", -0.005,0.2,10000,
              pscale=1.25,xlabel="hyper r")
    ggsave(paste0("Figures/",res.to.sim,"/sim_r_comp.png"), 
           dpi=300, height=8, width=8, units="in")
    
    paramcomp(SEOq="FALSE",Mods=Mods,Modnames=Modnames,param="CBtoK", 0,1,10000,
              pscale=2,xlabel="Current biomass:virgin biomass")
    ggsave(paste0("Figures/",res.to.sim,"/sim_CBtoK_subs.png"), 
           dpi=300, height=8, width=8, units="in")
    
    paramcomp(SEOq="TRUE",Mods=Mods,Modnames=Modnames,param="CBtoKseo", 0,0.9,10000,
              pscale=1.3,xlabel="Current biomass:virgin biomass")
    ggsave(paste0("Figures/",res.to.sim,"/sim_CBtoK.png"), 
           dpi=300, height=8, width=8, units="in")
    
    paramcomp(SEOq="FALSE",Mods=Mods,Modnames=Modnames,param="Fmsy", -0.005,0.1,10000,
              pscale=2,xlabel="Fmsy")
    ggsave(paste0("Figures/",res.to.sim,"/sim_Fmsy_subs.png"), 
           dpi=300, height=8, width=8, units="in")
    
    paramcomp(SEOq="TRUE",Mods=Mods,Modnames=Modnames,param="Fmsyseo", -0.005,0.075,10000,
              pscale=1.3,xlabel="Fmsy")
    ggsave(paste0("Figures/",res.to.sim,"/sim_Fmsy.png"), 
           dpi=300, height=8, width=8, units="in")
    
    paramcomp(SEOq="FALSE",Mods=Mods,Modnames=Modnames,param="Bmsy", 0,20000,10000,
              pscale=2,xlabel="Bmsy (mt)")
    ggsave(paste0("Figures/",res.to.sim,"/sim_Bmsy_subs.png"), 
           dpi=300, height=8, width=8, units="in")
    
    paramcomp(SEOq="TRUE",Mods=Mods,Modnames=Modnames,param="Bmsyseo", 0,80000,10000,
              pscale=1.5,xlabel="Bmsy (mt)")
    ggsave(paste0("Figures/",res.to.sim,"/sim_Bmsy.png"), 
           dpi=300, height=8, width=8, units="in")
    
    paramcomp(SEOq="FALSE",Mods=Mods,Modnames=Modnames,param="Stock.Status", 0,2.2,10000,
              pscale=2,xlabel="Stock Status (B2022:B40)")
    ggsave(paste0("Figures/",res.to.sim,"/sim_StockStatus_subs.png"), 
           dpi=300, height=8, width=8, units="in")
    
    paramcomp(SEOq="TRUE",Mods=Mods,Modnames=Modnames,param="Stock.Status.SEO", 0,2.2,10000,
              pscale=1.25,xlabel="Stock Status (B2022:B40)")
    ggsave(paste0("Figures/",res.to.sim,"/sim_StockStatus.png"), 
           dpi=300, height=8, width=8, units="in")
    
  ##############################################################################
    nsims<-15
    for (i in 3:nsims){  #i<-2

      load(file=paste("Model Output/",res.to.sim,"/simdat.",i,"post.Rdata", sep=""))
      
      Simulated <- coda::as.mcmc(sim)
      
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
      
      #save(post,file=paste("Model Output/", Mod.list[m],"_",niter/1000000,"m/post.Rdata", sep=""))
      #save(post,file=paste("D:/Groundfish Biometrics/MCMC_backups/", filename,".Rdata", sep=""))
      load(file=paste("Model Output/",res.to.sim,"/postsim",i,".Rdata", sep=""))
      
      #dir.create(paste("Figures/",res.to.sim,"/sim",i,sep=""))
      #dir.create(paste("Model Output/",res.to.sim,"/sim",i,sep=""))
      filename<-paste(res.to.sim,"/sim",i,sep="")
      
      B.ey<-par.ext(par="B",years=N+Fu,areai=1)
      B.ns<-par.ext("B",N+Fu,2)
      B.cs<-par.ext("B",N+Fu,3)
      B.ss<-par.ext("B",N+Fu,4)
      
      All.years<-seq(min(Years),max(Years)+Fu,by=1)
      FuYear<-Years+Fu
      sdlist<-c("EYKT", "NSEO", "CSEO", "SSEO")
      
      colch<-c("blue","red","purple","orange","forestgreen")
      vars<-list(B.ey,B.ns,B.cs,B.ss)
      
      #png(paste("Figures/",filename,"/Biomass_fit_FIX.png",sep=""),
      #    width=7,height=6,#width=9.5,height=8.5,
      #    units="in",res=1200)
      par(mfrow=c(2,2), mar=c(4,5,3,1))
      for (j in 1:4){  #j<-4
        envplot(vars[[j]],All.years,cols=c(colch[1],colch[1]),n=0,ylab="Biomass (t)",xlab="Year",
                ylim=c(0,17500),main=sdlist[j], cex.axis=0.8, cex.lab=1.25, cex.main=1.5)
        addpoints(post,Years=Years+0.1 ,Points=sCPUE[j,], errordat=cv.sCPUE[j,], scaler="qsCPUE", #scaler="qfCPUE" for example...
                  bar=1.95 ,col=colch[4], pch=17, cex=0.75)
        addpoints(post,Years=Years ,Points=B.obs[j,], errordat=cv.B[j,], scaler=0, #scaler="qfCPUE" for example...
                  bar=1.95 ,col=colch[2], pch=18, cex=1.2)
        abline(v=2022,  col=colch[3], lty=2)
        if (j == 2){
          legend(x="topright", c("posterior est.","observed ROV/sub","IPHC CPUE"), col=c(colch[1],colch[2],colch[4]),
                 text.col=c(colch[1],colch[2],colch[4]),border=F, bty="y", pch=c(NA,18,17), cex=1)
        } else {}
      }
      #dev.off()
      
      N.subd<-N; Years.cont<-Years
      plot_spm(post, Mod.title=phase3, filename=filename)
      
    } 
  



###############################################################################
###############################################################################
### SCRAP SCRAP SCRAP 
###############################################################################
###############################################################################

m<-2; d<-1; b<-1
for (v in 1:length(upvarlist)){
for (m in 1:length(Mod.list)) {  #m<-1
  for (d in 1:length(Derby.list)){ #d<-1
     for (b in 1:length(B1list)) { #b<-1
  OFLABC<-0.00  #0.02 (ABC) 0.032 (OFL). 0=nofishing
  HarvStrat<-0 #2.58 #1.96 #1.68 #0, 1.15
  Fu<-50
  #Derby.Eff<-Derby.list[d]
  
  Data<-load.data(YEAR=2022,
                  Derby.Eff = Derby.list[d],
                  DEsd=0.1,  #this is CV for derby
                  B1=B1list[b],
                  B2=B2list[b])
  list2env(Data,.GlobalEnv)
  
  phiB1<-phiB1
  phiB2<-phiB2
  
  phimu<-phimu
  phisig<-phisig
  
  bigKmu<-bigKmu
  bigKsigma<-bigKsigma
  
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
    filename<-paste(Mod.list[m],"_B1-",round(B1list[b],2),
                    "_B2-",round(B2list[b]),"_",
                    "upv",(upvarlist[v]),"_",
                    "phB1-",round(phiB1,2),"_",
                    "phB2-",round(phiB2,2),"_",
                    "Kmu-",round(bigKmu,1),"_",
                    "Ksig-",round(bigKsigma,1),"_",
                    "derb_",Derby.list[d],"_",
                    niter/1000,"k",sep="")
  } else {
  filename<-paste(Mod.list[m],"_B1-",round(B1list[b],2),
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
  B.ey<-par.ext(par="B",years=N+Fu,areai=1)
  B.ns<-par.ext("B",N+Fu,2)
  B.cs<-par.ext("B",N+Fu,3)
  B.ss<-par.ext("B",N+Fu,4)
  
  All.years<-seq(min(Years),max(Years)+Fu,by=1)
  FuYear<-Years+Fu
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





