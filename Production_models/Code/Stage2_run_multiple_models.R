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
  
  source("Code/2022_DSR_SAFE_models/Phase2/DATALOAD_SEO_YE_SPM_Func_1888.R")
  source("Code/2022_DSR_SAFE_models/Phase2/DATAPREP_SPM_SEO1888.R")
  source("Code/2022_DSR_SAFE_models/Phase2/PLOT_SPM1888.R")
  source("Code/Posterior_Plotting/YE_SPM_posterior_exams_Func.R")
}

#Mod0<-"PT_2i_pe02_T1e_r0405_d1_nofishing"
{
Mod1<-"PT2i_fullcatch"
biomassmod1 <-"PT_2i_pe01_T1e_rbeta_B1-1_B2-1_upv-3_derb_0_1400k"
biomassmod2 <-"PT_2i_pe01_T1e_rbeta_B1-1_B2-1_upv-3_derb_0.3_1400k"
biomassmod3 <-"PT_2i_pe01_T1e_rbeta_B1-1_B2-1_upv-3_derb_-0.3_1400k"

biomassmod4 <-"PT_2i_pe01_T1e_rbeta_B1-1_B2-1_upv-5_derb_0_1400k"
biomassmod5 <-"PT_2i_pe01_T1e_rbeta_B1-1_B2-1_upv-5_derb_0.3_1400k"
biomassmod6 <-"PT_2i_pe01_T1e_rbeta_B1-1_B2-1_upv-5_derb_-0.3_1400k"
#biomassmod2 <- "PT_2i_pe01_T1e_rbeta_B1-1.43_B2-19__derb_0_1000k"
#biomassmod1 <-"PT_2i_pe01_T1e_rbeta_B1-1_B2-1__derb_0_1200k"
#biomassmod3 <-"PT_2i_pe01_T1e_rbeta_B1.52-1.43_B2-35__derb_0_1000k"
#biomassmod4 <-"PT_2i_pe01_T1e_rbeta_B1.61-1.43_B2-67__derb_0_1000k"
}
#============================================================
#Organize which models to run, set derby if needed

Mod.list<-c(Mod1)

biomass.list<-c(biomassmod6,biomassmod5,biomassmod4) #,biomassmod4, biomassmod5, biomassmod6)

Derby.list<-c(-0.3,0.3,0)
DEsdlist<-c(0.1)

#Blist from r prior developed such that
B1list<-c(1) #c(1,1.483,1.247,1.241,1.374)
B2list<-c(1) #c(1,22.908,31.478,53.481,106.057)

upvarlist<-c(-5)
#===================================================================
#FLAG: good idea to run list on tiny chains (2K) to check for bugs in model code! 
#beta500k not fully converged, 3.5 hours
niter<-1600000 #300000 #350000
burnin<-0.4*niter
#50K = 56 minits ; not converged
#150K = 2.65 hours; not converged
#500K = 9 hours, not converged.  loosened K bounds and tightened upper r from 0.3 to 0.2
#800K = 15 hours still not converged.  A little wierdness on B1980 phi with values above 1?
#1200K = 24 hours, still not converged with K limited to 250. looks like truncating at 150-200K would be fair...  
#1400K = 1.07 days; 175 ceiling K; r not converged.  Looking at trace plots, need longer run
#1600 converged nicely; 1.26 days(!) 
set.seed(1234)

#!!!: FLAG: adjust initial values and data in loop for specific models!!! 
#m<-1; d<-1; b<-1
for (v in 1:length(upvarlist)){
for (m in 1:length(biomass.list)) {  #m<-1
  for (d in 1:length(Derby.list)){ #d<-1
     for (b in 1:length(B1list)) { #b<-1
  OFLABC<-0.00  #0.02 (ABC) 0.032 (OFL). 0=nofishing
  HarvStrat<-0 #2.58 #1.96 #1.68 #0, 1.15
  Fu<-1
  #Derby.Eff<-Derby.list[d]
  
  Data<-load.data.1888(YEAR=2022,
                  biomassmod = biomass.list[m],
                       bioY1 = 1994, bioYe = 2021,
                  Derby.Eff = Derby.list[d],
                  DEsd=0.1,
                  B1=B1list[b],
                  B2=B2list[b])
  list2env(Data,.GlobalEnv)
  # adjust data as needed for specific model
  #if (m == 1) {
  #B.obs[2:4,15]<-NA
  #cv.B[2:4,15]<-1
  #}
  upvar<- upvarlist[v]
  
  prep<-data.prep.1888()
  list2env(prep,.GlobalEnv)
  
    #adjust initial values for particular models as needed...
  #if (m == 1) {
  #  inits[[1]]$Tau1<-0; inits[[2]]$Tau1<-0; inits[[3]]$Tau1<-0
  #}
  
  tstart <- Sys.time()
  print(tstart)
  post <- jagsUI::jags(model.file=Mod.list[1], data=data,
                       parameters.to.save=parameters, inits=inits,
                       n.chains=3, parallel=T, n.iter=niter,
                       n.burnin=burnin, n.thin=(niter-burnin)/1000)  # can mess with n.adapt=
  print(Sys.time() - tstart)
  #5K = 5.6 minutes
  
  filename<-paste(Mod.list[1],"_B1-",round(B1list[b],2),
                  "_B2-",round(B2list[b]),"_",
                  "upv",upvarlist[v],"_",
                  "derb_",Derby.list[d],"_",
                  niter/1000,"k",sep="")
  
  dir.create(paste("Model Output/",filename,sep=""))
  save(post,file=paste("Model Output/", Mod.list[m],"_B1-",round(B1list[b],2),
                       "_B2-",round(B2list[b]),"_",
                       "upv",upvarlist[v],"_",
                       "derb_",Derby.list[d],"_",
                       
                       niter/1000,"k/post.Rdata", sep=""))
  #save(post,file=paste("Model Output/", Mod.list[m],"_",niter/1000000,"m/post.Rdata", sep=""))
  save(post,file=paste("D:/Groundfish Biometrics/MCMC_backups/", Mod.list[m],"_B1-",
                       round(B1list[b],2),
                       "_B2-",round(B2list[b]),"_",
                       "upv",upvarlist[v],"_",
                       "derb_",Derby.list[d],"_",
                       niter/1000,"k.Rdata", sep=""))
  
  #ggmcmc(ggs(post$samples),family="r",file=paste("Model Output/",Mod.list[i],"_",niter/1000000,"m/DIAG_r.pdf", sep=""))
  #ggmcmc(ggs(post$samples),family="K",file=paste("Model Output/",Mod.list[i],"_",niter/1000000,"m/DIAG_K.pdf", sep=""))
  #ggmcmc(ggs(post$samples),family="phi",file=paste("Model Output/",Mod.list[i],"_",niter/1000000,"m/DIAG_phi.pdf", sep=""))
  #ggmcmc(ggs(post$samples),family="sigma",file=paste("Model Output/",Mod.list[i],"_",niter/1000000,"m/DIAG_sigma.pdf", sep=""))
  #ggmcmc(ggs(post$samples),family="Tau1",file=paste("Model Output/",Mod.list[i],"_",niter/1000000,"m/DIAG_tau1.pdf", sep=""))
  #ggmcmc(ggs(post$samples),family="Tau3",file=paste("Model Output/",Mod.list[i],"_",niter/1000000,"m/DIAG_tau3.pdf", sep=""))
  print(Sys.time() - tstart)
  
    #load(file=paste("Model Output/",Mod1,"_",niter/1000000,"m/post.Rdata", sep=""))
  dir.create(paste("Figures/",filename,sep=""))
  MCMCtrace(post, params=c("r","K",#"qfCPUE","qsCPUE",
                           "Tau1","sigma"),
            ISB=TRUE, pdf=TRUE, Rhat=TRUE, 
            file=paste("Figures/",filename,"/Param_trace.pdf",sep=""))
  #Function for making plots...
    #fix broken fucking plots in function that I cant figure out why they are fucking up the scale!!
  #1) Biomass observed versus modeled 
  
  All.years<-seq(min(Years),max(Years)+Fu,by=1)
  FuYear<-Years+Fu
  
  plot_spm1888(post, Mod.title=Mod.list[m], filename=filename)
  Y<-seq(1888,(1888+N-1),1)
  colch<-c("forestgreen","black")
  Surplus<-MCMCchains(post,"Surplus")
  Tot.Catch<-MCMCchains(post,"C")
  png(paste("Figures/", filename,"/Surplus3.png",sep=""),
      width=7,height=6,#width=9.5,height=8.5,
      units="in",res=1200)
  par(mfrow=c(1,1), mar=c(4,5,3,1))
  #i<-1
  envplot(Tot.Catch,Years,cols=c(colch[1],colch[1]),n=0,
          ylab="Total Catch & Surplus Production (t)",xlab="Year",
          ylim=c(0,max(quantile(Tot.Catch,c(0.991)),quantile(Tot.Catch,c(0.99)))),
          main="",cex.axis=0.8, cex.lab=1.25, cex.main=1.5)
  addenv(Surplus,cols=c(colch[2],colch[2]), Y=Y)
  addenv(Tot.Catch,cols=c(colch[1],colch[1]), Y=Y)
  
  legend(x="topleft", c("Total Removals","Surplus"), col=c(colch[1],colch[2],colch[3]),
         text.col=c(colch[1],colch[2],colch[3]),border=F, bty="n", pch=c(NA),cex=1.2)
  
  dev.off()
}}}}

dev.off()

inits[[1]]$Tau1





