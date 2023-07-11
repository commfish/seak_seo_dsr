###############################################################################
## Code for projecting YE populations into the future under alternative harvests
## Projections involve RISK ANALYSIS to determine probability population
## remains above reference points in XX years into future
## 9/1/22 is prob above B40 in 50 years
## For use with posterior data from jags 
## Phil Joy
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
  
  YEAR<-2022
}

#============================================================
#OFL<-c(0.019,0.026,0.016,0.015)  #means
BestBio<-c(4604,1409,5994,6019)  #min PE, unif R
BestBio<-c(4619,1416,6070,6046) #mod PE, broad R

#SQ Bio
SQ<-read.csv("Data/SEO_YE_Biomass_subd_082322.csv") %>% filter(Year == YEAR)
SQ %>% mutate(harv.EY = Biomass_mt_lo90*Bio.EYKT.mt/Biomass_mt,
              harv.NS = Biomass_mt_lo90*Bio.NSEO.mt/Biomass_mt,
              harv.CS = Biomass_mt_lo90*Bio.CSEO.mt/Biomass_mt,
              harv.SS = Biomass_mt_lo90*Bio.SSEO.mt/Biomass_mt) -> SQ

SQBio.lo<-as.vector(c(SQ$harv.EY, SQ$harv.NS, SQ$harv.CS, SQ$harv.SS) )
sum(SQBio.lo); 0.032*sum(SQBio.lo); 0.026*sum(SQBio.lo); 0.02*sum(SQBio.lo)
SQBio.pt<-as.vector(c(SQ$Bio.EYKT.mt, SQ$Bio.NSEO.mt,SQ$Bio.CSEO.mt,SQ$Bio.SSEO.mt))
sum(SQBio.pt); 0.032*sum(SQBio.pt); 0.026*sum(SQBio.pt); 0.02*sum(SQBio.pt)

OFL<-c(0.013,0.019,0.012,0.011)  #medians  #min PE, unif R
OFL<-c(0.019,0.024,0.018,0.016)  #mod PE, broad R
ABC<-OFL*0.75

predux<-0.0
recABCp<-ABC*(1-predux)

recABC<-recABCp*BestBio
#===================================================================

set.seed(1234)

OUT<-data.frame()
SEOcalc<-data.frame()

FU<-50

FY<-2022
p <- 0.18815
sdlist<-c("EYKT", "NSEO", "CSEO", "SSEO")

ABClist<-list(BestBio*ABC,BestBio*(ABC*0.9),BestBio*(ABC*0.75) ,SQBio.lo*0.02, SQBio.pt*0.02)
harvpol<-c("SPM_maxABC","SPM_redABC10","SPM_redABC25","SQ_lo90","SQ_pt")

D.mods<-list("RISK3_B1-1_B2-1_upv-5_derb_0_recABC_-0perc_1500k/post.Rdata",   #min PE, unif R
             "RISK3_B1-1_B2-1_upv-5_derb_0.3_recABC_-0perc_1500k/post.Rdata",
             "RISK3_B1-1_B2-1_upv-5_derb_-0.3_recABC_-0perc_1500k/post.Rdata")
D.mods<-list("PHASE3_B1-1.48_B2-23_upv-3_phmu-0.7_phsig-1.2_Kmu-10.7_Ksig-9.7_derb_0_1500k/post.Rdata",   #mod PE, broad R R
             "PHASE3_B1-1.48_B2-23_upv-3_phmu-0.6_phsig-1.2_Kmu-10.8_Ksig-9.6_derb_0.3_1500k/post.Rdata",
             "PHASE3_B1-1.48_B2-23_upv-3_phmu-0.7_phsig-1.2_Kmu-10.6_Ksig-9.4_derb_-0.3_1500k/post.Rdata")
derb.lbl<-c("unbiased","biased low","biased hi")

j<-1

for (h in 1:length(ABClist)){  #h<-1
  for (d in 1:length(D.mods)){  #d<-1
    load(file=paste("Model Output/",D.mods[d],sep="")); post<-post
    
    #get vector of posteriors for R, K, epsilon, B22
    sigs<-MCMCchains(object=post, params=c("sigma"),ISB=FALSE,exact=TRUE)
    
    Biolist<-data.frame()
    
    for (s in 1:4){  #s<-1
      rs<-MCMCchains(object=post, params=c(paste("r[",s,"]",sep="")),ISB=FALSE,exact=TRUE)
      Ks<-MCMCchains(object=post, params=c(paste("K[",s,"]",sep="")),ISB=FALSE,exact=TRUE)
      B22<-MCMCchains(object=post, params=c(paste("B[",s,",",FY-1980,"]",sep="")),ISB=FALSE,exact=TRUE)
      Bmsy<-MCMCchains(object=post, params=c(paste("Bmsy[",s,"]",sep="")),ISB=FALSE,exact=TRUE)
      
      FC<-ABClist[[h]][s]
      
      sim<-data.frame()
      
      #Biolist<-vector()
      #for each posterior sample project the pop forward FU years
      for (i in 1:length(sigs)) {  #i<-1
        
        logB<-vector()
        logB[1] <- log(B22[i])
        
        B<-vector()
        B[1] <- B22[i]
        
        for (t in 1:FU) {   #t<-1
          #tau.log.pe <- 1/log(sigs[i]*sigs[i]+1)
          epsilon <- rnorm(1,0,sigs[i]*sigs[i])
          PE <- epsilon-(sigs[i]*sigs[i]/2)
          
          logvarF <- runif(1,-10,-4) #dunif(-10,0.5)   #log variance prior = very restrictive to keep it small
          sigmaF <- sqrt(exp(logvarF)) 
          F.eps <- rnorm(1,0,sigmaF*sigmaF)
          FE <- F.eps-(sigmaF*sigmaF/2)
          logFC <- log(FC*exp(FE))
          FCe <- exp(logFC)
          
          logB[t+1]<-log(max(B[t]+(rs[i]/p)*B[t]*(1-(B[t]/Ks[i])^p)-FCe,1)*exp(PE))
          B[t+1]<-exp(logB[t+1])
          
          #logB[i,t+1]<-log(max(B[i,t]+(r[i]/p)*B[i,t]*(1-(B[i,t]/K[i])^p)-FCe,1)*exp(PE))
          #B[i,t+1]<-exp(logB[i,t+1])
        }
        #after projection save the appropriate stuff
        sim[i,"MU"]<-sdlist[s]
        
        if (B22[i]-Bmsy[i] >= 0) {
          sim[i,"abv.msy.now"]<-1
        } else {
          sim[i,"abv.msy.now"]<-0
        }
        
        if (B[FU+1]-Bmsy[i] >= 0) {
          sim[i,"abv.msy.fu"]<-1
        } else {
          sim[i,"abv.msy.fu"]<-0
        }
        
        if (B[FU+1]-B22[i] >= 0) {
          sim[i,"pop.inc"]<-1
        } else {
          sim[i,"pop.inc"]<-0
        }
        
        Biolist[i,s]<-B[FU+1]
      }
      #calculate mu risks
      OUT[j,"harv_policy"]<-harvpol[h]#"SPM"#harvpol[hp]
      OUT[j,"Derb.eff"]<-derb.lbl[d]
      OUT[j,"MU"]<-sdlist[s]
      OUT[j,"Above_Bmsy_now"]<-sum(sim$abv.msy.now)/length(sigs)
      OUT[j,"Above_Bmsy_future"]<-sum(sim$abv.msy.fu)/length(sigs)
      OUT[j,"Pop_inc_future"]<-sum(sim$pop.inc)/length(sigs)
      j<-j+1
      
    }
    #claculate SEO risk
    OUT[j,"harv_policy"]<-harvpol[h]#"SPM"
    OUT[j,"Derb.eff"]<-derb.lbl[d]
    OUT[j,"MU"]<-"SEO"
    
    B22seo<-MCMCchains(object=post, params=c(paste("Bseo[",FY-1980,"]",sep="")),ISB=FALSE,exact=TRUE)
    Bmsyseo<-MCMCchains(object=post, params=c("Bmsyseo"),ISB=FALSE,exact=TRUE)
    
    SEO_Fu<-Biolist$V1+Biolist$V2+Biolist$V3+Biolist$V4
    
    seosim<-data.frame()
    for (i in 1:length(sigs)){
      if (B22seo[i]-Bmsyseo[i] >= 0) {
        seosim[i,"abv.msy.now"]<-1
      } else {
        seosim[i,"abv.msy.now"]<-0
      }
      
      if (SEO_Fu[i]-Bmsyseo[i] >= 0) {
        seosim[i,"abv.msy.fu"]<-1
      } else {
        seosim[i,"abv.msy.fu"]<-0
      }
      
      if (SEO_Fu[i]-B22seo[i] >= 0) {
        seosim[i,"pop.inc"]<-1
      } else {
        seosim[i,"pop.inc"]<-0
      }
    }
    
    OUT[j,"Above_Bmsy_now"]<-sum(seosim$abv.msy.now)/length(sigs)
    OUT[j,"Above_Bmsy_future"]<-sum(seosim$abv.msy.fu)/length(sigs)
    OUT[j,"Pop_inc_future"]<-sum(seosim$pop.inc)/length(sigs)
    
    j<-j+1
  }
}

write.csv(OUT,"Model Output/Risk_Analysis_minPE_unifR.csv")  
write.csv(OUT,"Model Output/Risk_Analysis_modPE_broadR.csv")
###################################################################################
###################################################################################

for (d in 1:length(D.mods)){  #d<-1
  load(file=paste("Model Output/",D.mods[d],sep="")); post<-post
  
  #get vector of posteriors for R, K, epsilon, B22
  sigs<-MCMCchains(object=post, params=c("sigma"),ISB=FALSE,exact=TRUE)
  
  Biolist<-data.frame()
  
  for (s in 1:4){  #s<-1
    rs<-MCMCchains(object=post, params=c(paste("r[",s,"]",sep="")),ISB=FALSE,exact=TRUE)
    Ks<-MCMCchains(object=post, params=c(paste("K[",s,"]",sep="")),ISB=FALSE,exact=TRUE)
    B22<-MCMCchains(object=post, params=c(paste("B[",s,",",FY-1980,"]",sep="")),ISB=FALSE,exact=TRUE)
    Bmsy<-MCMCchains(object=post, params=c(paste("Bmsy[",s,"]",sep="")),ISB=FALSE,exact=TRUE)
   
    FC<-ABClist[[d]][s]
    
    sim<-data.frame()
    
    #Biolist<-vector()
    #for each posterior sample project the pop forward FU years
    for (i in 1:length(sigs)) {  #i<-1
      
      logB<-vector()
      logB[1] <- log(B22[i])
      
      B<-vector()
      B[1] <- B22[i]
      
      for (t in 1:FU) {   #t<-1
        #tau.log.pe <- 1/log(sigs[i]*sigs[i]+1)
        epsilon <- rnorm(1,0,sigs[i]*sigs[i])
        PE <- epsilon-(sigs[i]*sigs[i]/2)
        
        logvarF <- runif(1,-10,-4) #dunif(-10,0.5)   #log variance prior = very restrictive to keep it small
        sigmaF <- sqrt(exp(logvarF)) 
        F.eps <- rnorm(1,0,sigmaF*sigmaF)
        FE <- F.eps-(sigmaF*sigmaF/2)
        logFC <- log(FC*exp(FE))
        FCe <- exp(logFC)
        
        logB[t+1]<-log(max(B[t]+(rs[i]/p)*B[t]*(1-(B[t]/Ks[i])^p)-FCe,1)*exp(PE))
        B[t+1]<-exp(logB[t+1])
        
        #logB[i,t+1]<-log(max(B[i,t]+(r[i]/p)*B[i,t]*(1-(B[i,t]/K[i])^p)-FCe,1)*exp(PE))
        #B[i,t+1]<-exp(logB[i,t+1])
      }
      #after projection save the appropriate stuff
      sim[i,"MU"]<-sdlist[s]
      
      if (B22[i]-Bmsy[i] >= 0) {
        sim[i,"abv.msy.now"]<-1
      } else {
        sim[i,"abv.msy.now"]<-0
      }
      
      if (B[FU+1]-Bmsy[i] >= 0) {
        sim[i,"abv.msy.fu"]<-1
      } else {
        sim[i,"abv.msy.fu"]<-0
      }
      
      if (B[FU+1]-B22[i] >= 0) {
        sim[i,"pop.inc"]<-1
      } else {
        sim[i,"pop.inc"]<-0
      }
      
      Biolist[i,s]<-B[FU+1]
    }
    #calculate mu risks
    OUT[j,"harv_policy"]<-harvpol[h]#"SPM"#harvpol[hp]
    OUT[j,"Derb.eff"]<-derb.lbl[d]
    OUT[j,"MU"]<-sdlist[s]
    OUT[j,"Above_Bmsy_now"]<-sum(sim$abv.msy.now)/length(sigs)
    OUT[j,"Above_Bmsy_future"]<-sum(sim$abv.msy.fu)/length(sigs)
    OUT[j,"Pop_inc_future"]<-sum(sim$pop.inc)/length(sigs)
    j<-j+1
    
  }
  #claculate SEO risk
  OUT[j,"harv_policy"]<-harvpol[h]#"SPM"
  OUT[j,"Derb.eff"]<-derb.lbl[d]
  OUT[j,"MU"]<-"SEO"
  
  B22seo<-MCMCchains(object=post, params=c(paste("Bseo[",FY-1980,"]",sep="")),ISB=FALSE,exact=TRUE)
  Bmsyseo<-MCMCchains(object=post, params=c("Bmsyseo"),ISB=FALSE,exact=TRUE)
  
  SEO_Fu<-Biolist$V1+Biolist$V2+Biolist$V3+Biolist$V4
  
  seosim<-data.frame()
  for (i in 1:length(sigs)){
    if (B22seo[i]-Bmsyseo[i] >= 0) {
      seosim[i,"abv.msy.now"]<-1
    } else {
      seosim[i,"abv.msy.now"]<-0
    }
    
    if (SEO_Fu[i]-Bmsyseo[i] >= 0) {
      seosim[i,"abv.msy.fu"]<-1
    } else {
      seosim[i,"abv.msy.fu"]<-0
    }
    
    if (SEO_Fu[i]-B22seo[i] >= 0) {
      seosim[i,"pop.inc"]<-1
    } else {
      seosim[i,"pop.inc"]<-0
    }
  }
  
  OUT[j,"Above_Bmsy_now"]<-sum(seosim$abv.msy.now)/length(sigs)
  OUT[j,"Above_Bmsy_future"]<-sum(seosim$abv.msy.fu)/length(sigs)
  OUT[j,"Pop_inc_future"]<-sum(seosim$pop.inc)/length(sigs)
  
  j<-j+1
}





#####################################################################################
#####################################################################################
#####################################################################################

#SCRAP
for (v in 1:length(upvarlist)){ #v<-1
  for (m in 1:length(Mod.list)) {  #m<-1
    for (d in 1:length(Derby.list)){ #d<-1
      for (b in 1:length(B1list)) { #b<-1
        OFLABC<-recABC  #0.02 (ABC) 0.032 (OFL). 0=nofishing
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
          filename<-paste("RISK3_B1-",round(B1list[b],2),
                          "_B2-",round(B2list[b]),"_",
                          "upv",(upvarlist[v]),"_",
                          "derb_",Derby.list[d],"_",
                          "recABC_-",predux*100,"perc_",
                          niter/1000,"k",sep="")
        } else {
          filename<-paste("RISK3_B1-",round(B1list[b],2),
                          "_B2-",round(B2list[b]),"_",
                          "upv",(upvarlist[v]),"_",
                          "derb_",Derby.list[d],"_",
                          "recABC_-",predux*100,"perc_",
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
        
        #ggmcmc(ggs(post$samples),family="r",file=paste("Model Output/",filename,"/DIAG_r.pdf", sep=""))
        #ggmcmc(ggs(post$samples),family="K",file=paste("Model Output/",filename,"/DIAG_K.pdf", sep=""))
        #ggmcmc(ggs(post$samples),family="phi",file=paste("Model Output/",filename,"/DIAG_phi.pdf", sep=""))
        #ggmcmc(ggs(post$samples),family="sigma",file=paste("Model Output/",filename,"/DIAG_sigma.pdf", sep=""))
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