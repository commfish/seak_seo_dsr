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
  source("Code/2022_DSR_SAFE_models/Phase1/DATAPREP_SPM_1980.R")
  source("Code/2022_DSR_SAFE_models/Phase1/PLOT_SPM80.R")
  source("Code/Posterior_Plotting/YE_SPM_posterior_exams_Func.R")
}

#Mod0<-"PT_2i_pe02_T1e_r0405_d1_nofishing"
{
Mod1<-"PT_2i_pe015_T1e_r0405_d1_nofishing"
Mod2<-"PT_2i_pe015_T1n_r0405_d1_nofishing"  #need to change Tau1 initial values...  
Mod3<-"PT_2i_pe015_T12e_r0405_d1_nofishing"
Mod4<-"PT_2i_pe01_T1e_r0405_d1_nofishing"

#drop '94 from best model
Mod5<-"PT_2i_pe015_T1n_r0405_no94_d1_nofishing"
Mod6<-"PT_2i_pe01_T1n_r0405_no94_d1_nofishing"

#r sensitivity
Mod7<-"PT_2i_pe015_T1e_r0410_d1_nofishing"
Mod8<-"PT_2i_pe015_T1e_r0510_d1_nofishing"
Mod9<-"PT_2i_pe015_T1e_r0610_d1_nofishing"
Mod10<-"PT_2i_pe015_T1e_runi_d1_nofishing"

Mod11<-"PT_2i_pe015_T1e_rg2r50_d1_nofishing"
Mod12<-"PT_2i_pe015_T1e_rg5r100_d1_nofishing"
Mod13<-"PT_2i_pe015_T1e_rg9r200_d1_nofishing"

Mod14<-"PT_2i_pe015_T1e_rgamv035_d1_nofishing"
Mod15<-"PT_2i_pe015_T1e_rgamv009_d1_nofishing"
Mod16<-"PT_2i_pe015_T1e_rgamv0034_d1_nofishing"
Mod17<-"PT_2i_pe015_T1e_rgamv0012_d1_nofishing"
Mod18<-"PT_2i_pe015_T1e_runif_d1_nofishing"

Mod19<-"PT_2i_pe015_T1e_rln0405_d1_nofishing"
Mod20<-"PT_2i_pe015_T1e_rln0410_d1_nofishing"
Mod21<-"PT_2i_pe015_T1e_rln0505_d1_nofishing"
Mod22<-"PT_2i_pe015_T1e_rln1010_d1_nofishing"

Mod23<-"PT_2i_pe01_T1e_rbeta"
Mod24<-"PT_2i_pe01_T1e_rbeta_uni"
}
#============================================================
#Organize which models to run, set derby if needed
Mod23<-"PT_2i_pe01_T1e_rbeta"

Mod.list<-c(#Mod14,Mod17, Mod18,
            Mod23)


Derby.list<-c(0,0.3,-0.3)
DEsdlist<-c(0.1,0.1,0.1)

#Blist from r prior developed such that
#1) uninformative beta(1,1)
#2) broadest prior Z=0.045, M=0.02, AaM=17
#3) best guess in my opinion, middle Z=0.05, M= 0.025,AaM=17
#4) narrowest, Z=0.055, M=0.02, AaM=17
B1list<-c(1) #c(1,1.483,1.247,1.241,1.374)
B2list<-c(1) #c(1,22.908,31.478,53.481,106.057)

upvarlist<- c(-3,-5)

#===================================================================
#FLAG: good idea to run list on tiny chains (2K) to check for bugs in model code! 
#beta500k not fully converged, 3.5 hours
niter<-1400000 #300000 #350000
burnin<-500000
#500K, 4 hours; mostly converged.  Further trimming of etas...
#800K 6 hours almost there.  etas converged with trunc, some r and k not quite there
#1m ~ 8.4 hours; pretty well converged.  A tad off on NSEO K - run 1.20m for final
#    other diagnostics not so great.  need more than 25% burnin and longer chains

set.seed(1234)

#!!!: FLAG: adjust initial values and data in loop for specific models!!! 
#m<-1; d<-1; b<-1
for (v in 1:length(upvarlist)){
for (m in 1:length(Mod.list)) {  #m<-1
  for (d in 1:length(Derby.list)){ #d<-1
     for (b in 1:length(B1list)) { #b<-1
  OFLABC<-0.00  #0.02 (ABC) 0.032 (OFL). 0=nofishing
  HarvStrat<-0 #2.58 #1.96 #1.68 #0, 1.15
  Fu<-1
  #Derby.Eff<-Derby.list[d]
  
  Data<-load.data(YEAR=2022,
                  Derby.Eff = Derby.list[d],
                  DEsd=0.1,  #this is CV for derby
                  B1=B1list[b],
                  B2=B2list[b])
  list2env(Data,.GlobalEnv)
  # adjust data as needed for specific model
  #if (m == 1) {
  #B.obs[2:4,15]<-NA
  #cv.B[2:4,15]<-1
  #}
  
  upvar<-upvarlist[v]
  
  prep<-data.prep()
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
  
  filename<-paste(Mod.list[m],"_B1-",round(B1list[b],2),
                  "_B2-",round(B2list[b]),"_",
                  "upv",round(upvarlist[v]),"_",
                  "derb_",Derby.list[d],"_",
                  niter/1000,"k",sep="")
  
  #load(file=paste("D:/Groundfish Biometrics/MCMC_backups/",Mod23,".Rdata", sep=""))
  #b<-1; d<-1; m<-1; niter<-1000000
  
  dir.create(paste("Model Output/",filename,sep=""))
  save(post,file=paste("Model Output/", filename,"/post.Rdata", sep=""))
  #save(post,file=paste("Model Output/", Mod.list[m],"_",niter/1000000,"m/post.Rdata", sep=""))
  save(post,file=paste("D:/Groundfish Biometrics/MCMC_backups/", filename,"k.Rdata", sep=""))
  dir.create(paste("Figures/",filename,sep=""))
  MCMCtrace(post, params=c("R.hyp" ,"r","Kseo","K","phi",#"qfCPUE","qsCPUE",
                           "Tau1","Tau2","Tau3", "sigma","eta", "rB1","rB2"),
            ISB=TRUE, pdf=TRUE, Rhat=TRUE, 
            file=paste("Figures/",filename,"/Param_trace.pdf",sep=""))
  ggmcmc(ggs(post$samples),family="r",file=paste("Model Output/",filename,"/DIAG_r.pdf", sep=""))
  ggmcmc(ggs(post$samples),family="K",file=paste("Model Output/",filename,"/DIAG_K.pdf", sep=""))
  ggmcmc(ggs(post$samples),family="phi",file=paste("Model Output/",filename,"/DIAG_phi.pdf", sep=""))
  ggmcmc(ggs(post$samples),family="sigma",file=paste("Model Output/",filename,"/DIAG_sigma.pdf", sep=""))
  #ggmcmc(ggs(post$samples),family="Tau1",file=paste("Model Output/",Mod.list[i],"_",niter/1000000,"m/DIAG_tau1.pdf", sep=""))
  #ggmcmc(ggs(post$samples),family="Tau3",file=paste("Model Output/",Mod.list[i],"_",niter/1000000,"m/DIAG_tau3.pdf", sep=""))
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





