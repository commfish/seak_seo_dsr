###############################################################################
## This is code for simulating data and biomass based on the results of SS-SPM 
## fits to the true data.  This code creates simulated data sets and plots results 
## relevant to the original model run
##
## Code order to conduct simulations
## identify which model you are working with and make sure consistent between scripts!
## 1) Simmulate data to look at and save values and plots: Data_simulator.R
## 2) Run the model on the simulated data (looooong time): Model_simulated_data.R
## 3) Examine simulation results: Sim_exam.R
##
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
  
  source("Code/2022_DSR_SAFE_models/Phase1/DATALOAD_SEO_YE_SPM_Func_1980.R")
  source("Code/2022_DSR_SAFE_models/Phase3/DATAPREP_SPM_1980_PHASE3.R")
  source("Code/2022_DSR_SAFE_models/Phase1/PLOT_SPM80.R")
  source("Code/Posterior_Plotting/YE_SPM_posterior_exams_Func.R")
}

#============================================================
#Organize which models to run
#ModPHASE3.1<-"PT2i_base_PHASE3_beta-phi"
ModPHASE3.2<-"PT2i_base_PHASE3_norm-phi"

phase3<-ModPHASE3.2
simMod<-"PHASE_3_sim"
#simMod<-"PHASE_3_sim_rBvers"

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

#phimu.list<-c(0.739)#c(0.7420787,0.6726209,0.7702393) #c(0.6802432) #c(0.7420787) #c(0.6802432)
#phisig.list<-c(0.213)#c(0.1796734,0.1690635,0.1795257) #c(0.2054276) #c(0.1796734) #c(0.2054276)

#bigKmu.list<-c(log(41497)) #c(10.62973,10.75576,10.58545) #c(10.65388) #c(10.6366) 
#bigKsigma.list<-c(log(8747)) #c(9.237758,9.390445,9.140554) #c(0.2812395) #c(0.2940348)

#upvarlist<- c(-5)

#===================================================================
#FLAG: good idea to run list on tiny chains (2K) to check for bugs in model code! 
#beta500k not fully converged, 3.5 hours

set.seed(1234)

#!!!: FLAG: adjust initial values and data in loop for specific models!!! 

  ests<-read.csv(paste0("Model Output/",res.to.sim,"/results_sum.csv"))
  ests<-as.data.frame(ests)
  
  view(ests)
  
  r<-ests$X50[grep(c("r\\["),ests$parameter)]
  #r<-ests$mean[grep(c("r\\["),ests$parameter)]
  #r<-c(0.03,0.03,0.03,0.03)
  R.hyp<-rep(ests$mean[grepl("R.hyp",ests$parameter)==T],4)
  
  rB1<-ests$X50[grep(c("rB1"),ests$parameter)]
  rB2<-ests$X50[grep(c("rB2"),ests$parameter)]
  
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
             logvar=logvar,
             rB1=rB1,rB2=rB2)
  realC <- KnC.obs
  realBy <- ExpByc
  
  str(post)
  
  parameters=c("B.obs","sCPUE" ,"KnC.obs","ExpByc","B", "epsilon") 
  
#--------------------------------------------------------------------------------

  #plot simulated data...
  
  ests<-read.csv(paste0("Model Output/",res.to.sim,"/results_sum.csv"))
  load(file=paste("Model Output/",res.to.sim,"/post.Rdata", sep=""))
  
  All.years<-seq(min(Years),max(Years)+1,by=1)
  FuYear<-Years+1
  sdlist<-c("EYKT", "NSEO", "CSEO", "SSEO")
  
  colch<-c("blue","red","purple","orange","forestgreen")
  
    #png(paste("Figures/",filename,"/Biomass_fit_FIX.png",sep=""),
  #    width=7,height=6,#width=9.5,height=8.5,
  #    units="in",res=1200)
  nsims<-100
  
  #par(mfrow=c(2,2), mar=c(4,5,3,1))
  {biomass.simvals<-as.data.frame(matrix(nrow=nsims, ncol=ncol(B.obs)))
  colnames(biomass.simvals)<-seq(1980,(1980+ncol(B.obs)-1),1)
  sCPUE.simvals<-biomass.simvals
  KnC.simvals<-biomass.simvals
  eby.simvals<-biomass.simvals
  
  EYKT.biomass.simvals<-biomass.simvals
  NSEO.biomass.simvals<-biomass.simvals
  CSEO.biomass.simvals<-biomass.simvals
  SSEO.biomass.simvals<-biomass.simvals
  
  EYKT.sCPUE.simvals<-sCPUE.simvals
  NSEO.sCPUE.simvals<-sCPUE.simvals
  CSEO.sCPUE.simvals<-sCPUE.simvals
  SSEO.sCPUE.simvals<-sCPUE.simvals
  
  EYKT.KnC.simvals<-sCPUE.simvals
  NSEO.KnC.simvals<-sCPUE.simvals
  CSEO.KnC.simvals<-sCPUE.simvals
  SSEO.KnC.simvals<-sCPUE.simvals
  
  EYKT.eby.simvals<-sCPUE.simvals
  NSEO.eby.simvals<-sCPUE.simvals
  CSEO.eby.simvals<-sCPUE.simvals
  SSEO.eby.simvals<-sCPUE.simvals
  
  modB.simvals<-as.data.frame(matrix(nrow=nsims, ncol=ncol(B.obs)+1))
  colnames(modB.simvals)<-seq(1980,(1980+ncol(B.obs)),1)
  EYKT.modB.simvals<-sCPUE.simvals
  NSEO.modB.simvals<-sCPUE.simvals
  CSEO.modB.simvals<-sCPUE.simvals
  SSEO.modB.simvals<-sCPUE.simvals
  
  eps.simvals<-as.data.frame(matrix(nrow=nsims, ncol=ncol(B.obs)+1))
  colnames(eps.simvals)<-seq(1980,(1980+ncol(B.obs)),1)
  EYKT.eps.simvals<-sCPUE.simvals
  NSEO.eps.simvals<-sCPUE.simvals
  CSEO.eps.simvals<-sCPUE.simvals
  SSEO.eps.simvals<-sCPUE.simvals}

  #save values from simulations... 
  getwd()
  write.csv(data, file=paste("Model Output/", res.to.sim,"/simulations/simulated_values.csv", sep=""))
  
  for (j in 1:nsims){  #j<-1
      sim <- jagsUI::jags(model.file=simMod, data=data,
                          parameters.to.save=parameters, #inits=inits,
                          n.chains=1, parallel=F, n.iter=1,
                          n.burnin=0)
      #save(sim,file=paste("Model Output/", res.to.sim,"/simdat.",i,"post.Rdata", sep=""))
      #load(file=paste("Model Output/",res.to.sim,"/simdat.",i,"post.Rdata", sep=""))
      
      Simulated <- coda::as.mcmc(sim)
      
      sim.B<-Simulated$sims.list$B.obs
      sim.KnC<-Simulated$sims.list$KnC.obs
      sim.sCPUE<-Simulated$sims.list$sCPUE
      sim.eby<-Simulated$sims.list$ExpByc
      sim.modB<-Simulated$sims.list$B
      sim.eps<-Simulated$sims.list$epsilon
      
      sim.B<-matrix(sim.B,nrow=4)
      sim.B[which(is.na(B.obs),arr.ind=T)]<-NA
      
      sim.sCPUE<-matrix(sim.sCPUE,nrow=4)
      sim.sCPUE[which(is.na(sCPUE),arr.ind=T)]<-NA
      
      sim.KnC<-matrix(sim.KnC,nrow=4)
      sim.eby<-matrix(sim.eby,nrow=4)
      sim.modB<-matrix(sim.modB,nrow=4)
      sim.eps<-matrix(sim.eps,nrow=4)
      
      sim.cv.c<-as.vector(rep(0,ncol(cv.sCPUE)))
      sim.sv.b<-as.vector(rep(0,ncol(cv.B)))
      
      EYKT.biomass.simvals[j,]<-sim.B[1,]
      EYKT.sCPUE.simvals[j,]<-sim.sCPUE[1,]
      EYKT.KnC.simvals[j,]<-sim.KnC[1,]
      EYKT.eby.simvals[j,]<-sim.eby[1,]
      EYKT.modB.simvals[j,]<-sim.modB[1,]
      
      NSEO.biomass.simvals[j,]<-sim.B[2,]
      NSEO.sCPUE.simvals[j,]<-sim.sCPUE[2,]
      NSEO.KnC.simvals[j,]<-sim.KnC[2,]
      NSEO.eby.simvals[j,]<-sim.eby[2,]
      NSEO.modB.simvals[j,]<-sim.modB[2,]
      
      CSEO.biomass.simvals[j,]<-sim.B[3,]
      CSEO.sCPUE.simvals[j,]<-sim.sCPUE[3,]
      CSEO.KnC.simvals[j,]<-sim.KnC[3,]
      CSEO.eby.simvals[j,]<-sim.eby[3,]
      CSEO.modB.simvals[j,]<-sim.modB[3,]
      
      SSEO.biomass.simvals[j,]<-sim.B[4,]
      SSEO.sCPUE.simvals[j,]<-sim.sCPUE[4,]
      SSEO.KnC.simvals[j,]<-sim.KnC[4,]
      SSEO.eby.simvals[j,]<-sim.eby[4,]
      SSEO.modB.simvals[j,]<-sim.modB[4,]
      
      eps.simvals[j,]<-sim.eps
    }
    
    biosims<-list(EYKT.biomass.simvals,NSEO.biomass.simvals,CSEO.biomass.simvals,SSEO.biomass.simvals)
    cpuesims<-list(EYKT.sCPUE.simvals,NSEO.sCPUE.simvals,CSEO.sCPUE.simvals,SSEO.sCPUE.simvals)
    KnCsims<-list(EYKT.KnC.simvals,NSEO.KnC.simvals,CSEO.KnC.simvals,SSEO.KnC.simvals)
    ebysims<-list(EYKT.eby.simvals,NSEO.eby.simvals,CSEO.eby.simvals,SSEO.eby.simvals)
    modBsims<-list(EYKT.modB.simvals,NSEO.modB.simvals,CSEO.modB.simvals,SSEO.modB.simvals)
    epssims<-eps.simvals
#-------------------------------------------------------------------------------
#plot simulated biomass versus posterior biomass... 
    colch<-c("blue","red","purple","orange","forestgreen")
    
    B.ey<-par.ext(par="B",years=N+1,areai=1)
    B.ns<-par.ext("B",N+1,2)
    B.cs<-par.ext("B",N+1,3)
    B.ss<-par.ext("B",N+1,4)
    vars<-list(B.ey,B.ns,B.cs,B.ss)
    
    png(file=paste("Model Output/", res.to.sim,"/simulations/simulated_vs_posterior_biomass.png", sep=""),
        width=8,height=8,#width=9.5,height=8.5,
        units="in",res=1200)
    
    par(mfrow=c(2,2), mar=c(4,5,3,1))
    for (i in 1:4){  #i<-1
      envplot(vars[[i]],All.years,cols=c(colch[1],colch[1]),n=0,ylab="Biomass (t)",xlab="Year",
              ylim=c(0,17500),main=sdlist[i], cex.axis=0.8, cex.lab=1.25, cex.main=1.5)
      
      mean.simbio<-data.frame()
      
      for (k in 1:ncol(B.obs)) { #k<-1
        mean.simbio["mean.simbio",k]<-mean(as.data.frame(modBsims[[i]])[,k], na.rm=T)
        mean.simbio["sd.simbio",k]<-sqrt(var(as.data.frame(modBsims[[i]])[,k], na.rm=T))
        mean.simbio["cv.simbio",k]<-mean.simbio["sd.simbio",k]/mean.simbio["mean.simbio",k]
        mean.simbio["lo99.simbio",k]<-quantile(as.data.frame(modBsims[[i]])[,k],0.01, na.rm=T)
        mean.simbio["hi99.simbio",k]<-quantile(as.data.frame(modBsims[[i]])[,k],0.99, na.rm=T)
        mean.simbio["lo95.simbio",k]<-quantile(as.data.frame(modBsims[[i]])[,k],0.025, na.rm=T)
        mean.simbio["hi95.simbio",k]<-quantile(as.data.frame(modBsims[[i]])[,k],0.975, na.rm=T)
        mean.simbio["lo50.simbio",k]<-quantile(as.data.frame(modBsims[[i]])[,k],0.25, na.rm=T)
        mean.simbio["hi50.simbio",k]<-quantile(as.data.frame(modBsims[[i]])[,k],0.75, na.rm=T)
      }
      x1<-seq(1980,(1980+ncol(B.obs)-1),1)
      y1<-c(t(mean.simbio["lo99.simbio",]))
      y2<-rev(c(t(mean.simbio["hi99.simbio",])))
      polygon(c(x1,rev(x1)),c(y1,y2), col=alpha("goldenrod", alpha=0.2), border = NA)   #,border=NA,col=adjustcolor(cols[1], alpha.f=.15))
      
      x1<-seq(1980,(1980+ncol(B.obs)-1),1)
      y1<-c(t(mean.simbio["lo95.simbio",]))
      y2<-rev(c(t(mean.simbio["hi95.simbio",])))
      polygon(c(x1,rev(x1)),c(y1,y2), col=alpha("goldenrod", alpha=0.2), border = NA)   #,border=NA,col=adjustcolor(cols[1], alpha.f=.15))
      
      y1<-c(t(mean.simbio["lo50.simbio",]))
      y2<-rev(c(t(mean.simbio["hi50.simbio",])))
      polygon(c(x1,rev(x1)),c(y1,y2), col=alpha("goldenrod", alpha=0.3), border = NA) 
      lines(c(t(mean.simbio["mean.simbio",]))~seq(1980,(1980+ncol(B.obs)-1),1),
             col="goldenrod", pch=7, cex=1.5, lwd=2)

      if (i == 2){
        legend(x="topright", c("posterior est. biomass","simulated biomass"),
               col=c(colch[1],"goldenrod"),
               text.col=c(colch[1],"goldenrod"),border=F, bty="y", pch=c(15,15,45,0,7), cex=1.25)
      } else {}
    }  
  dev.off()
#-------------------------------------------------------------------------------- 
#PLOT simulated biomass
    colch<-c("blue","red","purple","orange","forestgreen")
    
    B.ey<-par.ext(par="B",years=N+1,areai=1)
    B.ns<-par.ext("B",N+1,2)
    B.cs<-par.ext("B",N+1,3)
    B.ss<-par.ext("B",N+1,4)
    vars<-list(B.ey,B.ns,B.cs,B.ss)
    
    png(file=paste("Model Output/", res.to.sim,"/simulations/simulated_vs_posterior_biomass_horsehair.png", sep=""),
        width=8,height=8,#width=9.5,height=8.5,
        units="in",res=1200)
    
  par(mfrow=c(2,2), mar=c(4,5,3,1))
  for (i in 1:4){  #i<-1
    envplot(vars[[i]],All.years,cols=c(colch[1],colch[1]),n=0,ylab="Biomass (t)",xlab="Year",
            ylim=c(0,17500),main=sdlist[i], cex.axis=0.8, cex.lab=1.25, cex.main=1.5)
    for (j in 1:nsims) {  #j<-1
      pts<-as.vector(as.numeric(biosims[[i]][j,]))
      points(pts~seq(1980,(1980+ncol(B.obs)-1),1), pch=0, col=alpha("darkred",0.2))
      
      modB<-as.vector(as.numeric(modBsims[[i]][j,]))
      lines(modB~seq(1980,(1980+ncol(B.obs)-1),1), pch=0, col=alpha("darkcyan",0.2))
    }
    
    mean.simbio<-data.frame()
    
    for (k in 1:ncol(B.obs)) { #k<-1
      mean.simbio["mean.simbio",k]<-mean(as.data.frame(biosims[[i]])[,k], na.rm=T)
      mean.simbio["sd.simbio",k]<-sqrt(var(as.data.frame(biosims[[i]])[,k], na.rm=T))
      mean.simbio["cv.simbio",k]<-mean.simbio["sd.simbio",k]/mean.simbio["mean.simbio",k]
      mean.simbio["lo90.simbio",k]<-quantile(as.data.frame(biosims[[i]])[,k],0.05, na.rm=T)
      mean.simbio["hi90.simbio",k]<-quantile(as.data.frame(biosims[[i]])[,k],0.95, na.rm=T)
      
    }
    
    points(c(t(mean.simbio["mean.simbio",]))~seq(1980,(1980+ncol(B.obs)-1),1),
           col="darkred", pch=7, cex=1.5)
    arrows(x0=seq(1980,(1980+ncol(B.obs)-1),1), 
           y0=c(t(mean.simbio["lo90.simbio",])), 
           x1=seq(1980,(1980+ncol(B.obs)-1),1), 
           y1=c(t(mean.simbio["hi90.simbio",])), 
           code=3, 
           col="darkred", lwd=1, angle=90, length=0.05)
    
    addpoints(post,Years=Years ,Points=B.obs[i,], errordat=cv.B[i,], scaler=0, #scaler="qfCPUE" for example...
              bar=1.95 ,col=colch[2], pch=18, cex=1.2)
    
    if (i == 2){
      legend(x="topright", c("posterior est.","observed ROV/sub",
                             "simulated biomass","simulated ROV",
                             "mean sim. ROV"),
             col=c(colch[1],colch[2],"darkcyan",alpha("darkred",0.2),"darkred"),
             text.col=c(colch[1],colch[2],"darkcyan","darkred","darkred"),border=F, bty="y", pch=c(NA,18,45,0,7), cex=1.25)
    } else {}
    
  }
  dev.off()
  
  #-------------------------------------------------------------------------------
  #PLOT simulated CPUE
  
  B.ey<-par.ext(par="B",years=N+1,areai=1)
  B.ns<-par.ext("B",N+1,2)
  B.cs<-par.ext("B",N+1,3)
  B.ss<-par.ext("B",N+1,4)
  vars<-list(B.ey,B.ns,B.cs,B.ss)
  
  png(file=paste("Model Output/", res.to.sim,"/simulations/simulated_vs_posterior_biomass_CPUE_ht.png", sep=""),
      width=8,height=8,#width=9.5,height=8.5,
      units="in",res=1200)
  
  par(mfrow=c(2,2), mar=c(4,5,3,1))
  for (i in 1:4){  #i<-1
    envplot(vars[[i]],All.years,cols=c(colch[1],colch[1]),n=0,ylab="Biomass (t)",xlab="Year",
            ylim=c(0,17500),main=sdlist[i], cex.axis=0.8, cex.lab=1.25, cex.main=1.5)
    for (j in 1:nsims) {  #j<-1
      pts2<-as.vector(as.numeric(cpuesims[[i]][j,]))/exp(logqs[i])
      points(pts2~seq(1980,(1980+ncol(B.obs)-1),1), pch=0, col=alpha("yellow4",0.1))
      
      modB<-as.vector(as.numeric(modBsims[[i]][j,]))
      lines(modB~seq(1980,(1980+ncol(B.obs)-1),1), pch=0, col=alpha("darkcyan",0.2))
    }
    
    mean.simbio<-data.frame()
    mean.simcpue<-data.frame()
    
    for (k in 1:ncol(B.obs)) { #k<-1
      mean.simcpue["mean.simcpue",k]<-mean(as.data.frame(cpuesims[[i]])[,k], na.rm=T)
      mean.simcpue["sd.simcpue",k]<-sqrt(var(as.data.frame(cpuesims[[i]])[,k], na.rm=T))
      mean.simcpue["cv.simcpue",k]<-mean.simcpue["sd.simcpue",k]/mean.simcpue["mean.simcpue",k]
      mean.simcpue["lo90.simcpue",k]<-quantile(as.data.frame(cpuesims[[i]])[,k],0.05, na.rm=T)
      mean.simcpue["hi90.simcpue",k]<-quantile(as.data.frame(cpuesims[[i]])[,k],0.95, na.rm=T)
    }
    logqs
    points(c(t(mean.simcpue["mean.simcpue",]))/exp(logqs[i])~seq(1980,(1980+ncol(B.obs)-1),1),
           col="yellow4", pch=7, cex=1.5)
    lines(c(t(mean.simcpue["mean.simcpue",]))/exp(logqs[i])~seq(1980,(1980+ncol(B.obs)-1),1),
           col="yellow4", pch=7, cex=1.5)
    arrows(x0=seq(1980,(1980+ncol(B.obs)-1),1), 
           y0=c(t(mean.simcpue["lo90.simcpue",]))/exp(logqs[i]), 
           x1=seq(1980,(1980+ncol(B.obs)-1),1), 
           y1=c(t(mean.simcpue["hi90.simcpue",]))/exp(logqs[i]), 
           code=3, 
           col="yellow4", lwd=1, angle=90, length=0.05)
    
    addpoints(post,Years=Years+0.1 ,Points=sCPUE[i,], errordat=cv.sCPUE[i,], scaler="qsCPUE", #scaler="qfCPUE" for example...
              bar=1.95 ,col=colch[4], pch=17, cex=1)
    
    if (i == 2){
      legend(x="topright", c("posterior est.","observed IPHC CPUE","simulated biomass",
                             "simulated IPHC CPUE","mean sim. IPHC CPUE"),
             col=c(colch[1],colch[4],"darkcyan",alpha("yellow4",0.2),"yellow4"),
             text.col=c(colch[1],colch[4],"darkcyan","yellow4","yellow4"),
             border=F, bty="y", pch=c(NA,17,45,0,7), cex=1.25)
    } else {}
  }
  dev.off()
  
  #--------------------------------------------------------------------------
  #PLOT simulated catch and bycatch... 
  
  KnC.ey<-par.ext(par="KnC",years=N,areai=1)
  KnC.ns<-par.ext("KnC",N,2)
  KnC.cs<-par.ext("KnC",N,3)
  KnC.ss<-par.ext("KnC",N,4)
  
  colch<-c("blue","forestgreen","purple")
  vars<-list(KnC.ey,KnC.ns,KnC.cs,KnC.ss)
  
  png(file=paste("Model Output/", res.to.sim,"/simulations/simulated_vs_posterior_catch.png", sep=""),
      width=8,height=8,#width=9.5,height=8.5,
      units="in",res=1200)
  
  par(mfrow=c(2,2), mar=c(4,5,3,1))
  for (i in 1:4){  #i<-1
    envplot(vars[[i]],Years,cols=c(colch[1],colch[1]),n=0,ylab="Known catches (t)",xlab="Year",
            ylim=c(0,max(KnC.obs*1.25)),main=sdlist[i],cex.axis=0.8, cex.lab=1.25, cex.main=1.5)
    for (j in 1:nsims) {  #j<-1
      pts<-as.vector(as.numeric(KnCsims[[i]][j,]))
      points(pts~seq(1980,(1980+ncol(B.obs)-1),1), pch=0, col=alpha("limegreen",0.1))
    }
    
    mean.simKnC<-data.frame()
    #mean.simcpue<-data.frame()
    
    for (k in 1:ncol(B.obs)) { #k<-1
      mean.simKnC["mean.simKnC",k]<-mean(as.data.frame(KnCsims[[i]])[,k], na.rm=T)
      mean.simKnC["sd.simKnC",k]<-sqrt(var(as.data.frame(KnCsims[[i]])[,k], na.rm=T))
      mean.simKnC["cv.simKnC",k]<-mean.simKnC["sd.simKnC",k]/mean.simKnC["mean.simKnC",k]
    }
    
    addpoints(post,Years=Years+0.1 ,Points=c(t(mean.simKnC["mean.simKnC",])), 
              errordat= c(t(mean.simKnC["cv.simKnC",])),#mean.simbio["cv.simbio",], 
              scaler=0, #scaler="qfCPUE" for example...
              bar=1.95 ,col="limegreen", pch=16, cex=1.5)
    
    addpoints(post,Years=Years ,Points=KnC.obs[i,] ,errordat=cv.KnC[i,] ,scaler=0 ,
              bar=1.95 ,col=colch[2], pch=18,cex=1.2)
    
    if (i == 2){
      legend(x="topright", c("posterior est.","observed known catch",
                             "simulated knownn catch","mean sim. known catch"),
             col=c(colch[1],colch[2],alpha("limegreen",0.2),"limegreen"),
             text.col=c(colch[1],colch[2],"limegreen","limegreen"),border=F, bty="y",
             pch=c(NA,16,0,7), cex=1.25)
    } else {}
  }
  dev.off()
#------------------------------------------------------------------------------
# Expected Bycatch
  
  By.ey<-par.ext(par="By",years=N,areai=1)
  By.ns<-par.ext("By",N,2)
  By.cs<-par.ext("By",N,3)
  By.ss<-par.ext("By",N,4)
  
  colch<-c("blue","orange","purple")
  vars<-list(By.ey,By.ns,By.cs,By.ss)
  
  png(file=paste("Model Output/", res.to.sim,"/simulations/simulated_vs_posterior_bycatch.png", sep=""),
      width=8,height=8,#width=9.5,height=8.5,
      units="in",res=1200)
  
  par(mfrow=c(2,2), mar=c(4,5,3,1))
  for (i in 1:4){  #i<-1
    envplot(vars[[i]],Years,cols=c(colch[1],colch[1]),n=0,ylab="Bycatch (t)",xlab="Year",
            ylim=c(0,max(ExpByc*1.25)),main=sdlist[i],cex.axis=0.8, cex.lab=1.25, cex.main=1.5)
    for (j in 1:nsims) {  #j<-1
      pts<-as.vector(as.numeric(ebysims[[i]][j,]))
      points(pts~seq(1980,(1980+ncol(B.obs)-1),1), pch=0, col=alpha("darkslategray",0.1))
    }
    
    mean.simeby<-data.frame()
    #mean.simcpue<-data.frame()
    
    for (k in 1:ncol(B.obs)) { #k<-1
      mean.simeby["mean.simeby",k]<-mean(as.data.frame(ebysims[[i]])[,k], na.rm=T)
      mean.simeby["sd.simeby",k]<-sqrt(var(as.data.frame(ebysims[[i]])[,k], na.rm=T))
      mean.simeby["cv.simeby",k]<-mean.simeby["sd.simeby",k]/mean.simeby["mean.simeby",k]
      mean.simeby["lo90.simeby",k]<-quantile(as.data.frame(ebysims[[i]])[,k],0.05, na.rm=T)
      mean.simeby["hi90.simeby",k]<-quantile(as.data.frame(ebysims[[i]])[,k],0.95, na.rm=T)
    }
    
    points(c(t(mean.simeby["mean.simeby",]))~seq(1980,(1980+ncol(B.obs)-1),1),
           col="darkslategray", pch=7, cex=1.5)
    
    arrows(x0=seq(1980,(1980+ncol(B.obs)-1),1), 
           y0=c(t(mean.simeby["lo90.simeby",])), 
           x1=seq(1980,(1980+ncol(B.obs)-1),1), 
           y1=c(t(mean.simeby["hi90.simeby",])), 
           code=3, 
           col="darkslategray", lwd=1, angle=90, length=0.05)
    
    #addpoints(post,Years=Years+0.1 ,Points=c(t(mean.simeby["mean.simeby",])), 
    #          errordat= c(t(mean.simeby["cv.simeby",])),#mean.simbio["cv.simbio",], 
    #          scaler=0, #scaler="qfCPUE" for example...
    #          bar=1.95 ,col="darkslategray", pch=7, cex=1.5)
    
    addpoints(post,Years=Years ,Points=ExpByc[i,] ,errordat=cv.ExpByc[i,] ,scaler=0 ,
              bar=1.95 ,col="violet", pch=18, cex=1.2)
    
    if (i == 2){
      legend(x="topright", c("posterior est.","observed expected bycatch",
                             "simulated expected bycatch","mean sim. expected bycatch"),
             col=c(colch[1],"violet",alpha("darkslategray",0.2),"darkslategray"),
             text.col=c(colch[1],"violet","darkslategray","darkslategray"),
             border=F, bty="y", pch=c(NA,16,0,7), cex=1.25)
    } else {}
  }
  dev.off()
  #------------------------------------------------------------------------------
  #Epsilon
  
  pe<-MCMCchains(post,"epsilon")
  pe2<-MCMCchains(post,"PE")
  par(mfrow=c(1,1))
  #envplot(pe,1980:(1980+N),cols=c("black","black"),n=0,ylab="epsilon",xlab="Year",ylim=c(-0.2,0.2))
  #abline(h=0,col="red",lty=2)
  Years2<-seq(1980,2023,1)
  
  png(file=paste("Model Output/", res.to.sim,"/simulations/simulated_vs_posterior_PE.png", sep=""),
      width=5,height=5,#width=9.5,height=8.5,
      units="in",res=1200)
  
  envplot(pe,Years2, #1980:(1980+N),
          cols=c("black","black"),n=0,ylab="process error",xlab="Year",ylim=c(-0.12,0.12))
  abline(h=0,col="red",lty=2)
  for (j in 1:nsims) {  #j<-1
    pts<-as.vector(as.numeric(epssims[j,]))
    points(pts~seq(1980,(1980+ncol(B.obs)),1), pch=0, col=alpha("coral3",0.1))
  }
  
  mean.simeps<-data.frame()
  for (k in 1:(ncol(B.obs)+1)) { #k<-1
    mean.simeps["mean.simeps",k]<-mean(as.data.frame(epssims)[,k], na.rm=T)
    mean.simeps["sd.simeps",k]<-sqrt(var(as.data.frame(epssims)[,k], na.rm=T))
    mean.simeps["cv.simeps",k]<-mean.simeps["sd.simeps",k]/mean.simeps["mean.simeps",k]
    mean.simeps["lo90.simeps",k]<-quantile(as.data.frame(epssims)[,k],0.05, na.rm=T)
    mean.simeps["hi90.simeps",k]<-quantile(as.data.frame(epssims)[,k],0.95, na.rm=T)
  }
  
  points(c(t(mean.simeps["mean.simeps",]))~seq(1980,(1980+ncol(B.obs)),1),
         col="coral3", pch=7, cex=1.5)
  arrows(x0=seq(1980,(1980+ncol(B.obs)),1), 
         y0=c(t(mean.simeps["lo90.simeps",])), 
         x1=seq(1980,(1980+ncol(B.obs)),1), 
         y1=c(t(mean.simeps["hi90.simeps",])), 
         code=3, 
         col="coral3", lwd=1, angle=90, length=0.05)
  
  legend(x="top", c("posterior est.",
                             "simulated epsilon","mean sim. epsilon"),
             col=c("darkgrey",alpha("coral3",0.2),"coral3"),
             text.col=c("darkgrey","coral3","coral3"),
             border=F, bty="y", pch=c(15,0,7), cex=1.25)

  dev.off()
  