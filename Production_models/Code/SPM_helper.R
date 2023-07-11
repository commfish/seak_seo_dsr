#################################################################################
## SPM_helper 
##
## Canned functions for running surplus production models
##
################################################################################
library(dplyr)

## Load Data for Stage 1 of 3 stage approach

load.data<-function (YEAR = 2022, Derby.Eff = Derby.Eff, DEsd=DEsd,
                     B1=B1,B2=B2){
  #=========================================================
  #number of management districts
  Subd<-4
  
  ## 1) Halibut bycatch - longest time series; use this to establish number of years
  # Load expected bycatch generated in Code/Data_Generation/IPHC_Survey_Expected_Bycatch.R
  
  #=========================
  ## Known Catches fo KnC.obs... 
  #1) contemporary removals from 1980 to present... 
  #Removals.subd<-read.csv("Data/SEsubdistrict_YE_removals.csv")  #YE reconstruction from Rhea
  Removals.subd<-read.csv(paste0("Data_processing/Data/SE_YE_known_removals_1980-",YEAR,".csv"))
  #str(Removals.subd); unique(Removals.subd$Year)
  
  
  #*************************************
  #* Length of Data set; catch data will always be the longest data series
  Subd<-4
  #YEAR <-2022
  Year1<-min(Removals.subd$Year)
  Yearmax<-max(Removals.subd$Year)
  
  Years<-seq(min(Removals.subd$Year),YEAR,by=1)
  N<-length(Years)  
  
  #***************************************
  
  ey.rem<-c(Removals.subd$tot.rem.mt[Removals.subd$Mgt.Area == "EYKT"]); length(ey.rem)
  ns.rem<-c(Removals.subd$tot.rem.mt[Removals.subd$Mgt.Area == "NSEO"]); length(ns.rem)
  cs.rem<-c(Removals.subd$tot.rem.mt[Removals.subd$Mgt.Area == "CSEO"]); length(cs.rem)
  ss.rem<-c(Removals.subd$tot.rem.mt[Removals.subd$Mgt.Area == "SSEO"]); length(ss.rem)
  #seo.rem<-c(Removals.pre$tot.rem.mt,colSums(rbind(ey.rem,ns.rem,cs.rem,ss.rem))); length(seo.rem)
  
  cv.subd<-0.1
  cv.eyrem<-c(rep(cv.subd,N)); length(cv.eyrem)
  cv.nsrem<-c(rep(cv.subd,N)); length(cv.nsrem)
  cv.csrem<-c(rep(cv.subd,N)); length(cv.csrem)
  cv.ssrem<-c(rep(cv.subd,N)); length(cv.ssrem)
  #cv.seorem<-c(Removals.pre$Remove.cv,rep(cv.subd,N.subd)); length(cv.seorem)
  
  #Data for subdistrict models:
  KnC.obs<-rbind(ey.rem,ns.rem,cs.rem,ss.rem); length(KnC.obs)
  cv.KnC<-rbind(cv.eyrem, cv.nsrem, cv.csrem, cv.ssrem); length(cv.KnC)
  
  #Data for SEO model
  #KnCseo.obs<-seo.rem
  #cv.KnCseo<-cv.seorem
  
  #YE landings in Halibut fishery
  Lnd.ey<-c(Removals.subd$tot.hal.by[Removals.subd$Mgt.Area == "EYKT"])
  Lnd.ns<-c(Removals.subd$tot.hal.by[Removals.subd$Mgt.Area == "NSEO"])
  Lnd.cs<-c(Removals.subd$tot.hal.by[Removals.subd$Mgt.Area == "CSEO"])
  Lnd.ss<-c(Removals.subd$tot.hal.by[Removals.subd$Mgt.Area == "SSEO"])
  #Lnd.seo<-c(rep(0,N.pre), colSums(rbind(Lnd.ey,Lnd.ns,Lnd.cs,Lnd.ss)))
  
  Lnd.By<-rbind(Lnd.ey,Lnd.ns,Lnd.cs,Lnd.ss); length(Lnd.By)
  #Lndseo.By<-Lnd.seo
  
  #=========================================================================================
  # Load expected bycatch generated in Code/Data_Generation/IPHC_Survey_Expected_Bycatch.R
  ExpBy<-read.csv(paste0("Data_processing/Data/SEO_expBy_",YEAR,".csv"))
  N.exp<-nrow(ExpBy)/4
  expY1<-min(ExpBy$Year)
  expYe<-max(ExpBy$Year)
  
  ExpBy.ey<-c(rep(NA,N-N.exp-(YEAR-expYe)),ExpBy$expBy32_mt.mean[ExpBy$mngmt.area=="EYKT"],
              rep(mean(ExpBy$expBy32_mt.mean[ExpBy$mngmt.area=="EYKT"][(N.exp-4):N.exp]),(YEAR-expYe))); length(ExpBy.ey)
  ExpBy.ns<-c(rep(NA,N-N.exp-(YEAR-expYe)),ExpBy$expBy32_mt.mean[ExpBy$mngmt.area=="NSEO"],
              rep(mean(ExpBy$expBy32_mt.mean[ExpBy$mngmt.area=="NSEO"][(N.exp-4):N.exp]),(YEAR-expYe))); length(ExpBy.ey)
  ExpBy.cs<-c(rep(NA,N-N.exp-(YEAR-expYe)),ExpBy$expBy32_mt.mean[ExpBy$mngmt.area=="CSEO"],
              rep(mean(ExpBy$expBy32_mt.mean[ExpBy$mngmt.area=="CSEO"][(N.exp-4):N.exp]),(YEAR-expYe))); length(ExpBy.ey)
  ExpBy.ss<-c(rep(NA,N-N.exp-(YEAR-expYe)),ExpBy$expBy32_mt.mean[ExpBy$mngmt.area=="SSEO"],
              rep(mean(ExpBy$expBy32_mt.mean[ExpBy$mngmt.area=="SSEO"][(N.exp-4):N.exp]),(YEAR-expYe))); length(ExpBy.ey)
  
  cv.ExpBy.ey<-c(rep(1,N-N.exp-(YEAR-expYe)),ExpBy$expBy32.cv2[ExpBy$mngmt.area=="EYKT"],rep(1,(YEAR-expYe))); length(cv.ExpBy.ey)
  cv.ExpBy.ns<-c(rep(1,N-N.exp-(YEAR-expYe)),ExpBy$expBy32.cv2[ExpBy$mngmt.area=="NSEO"],rep(1,(YEAR-expYe))); length(cv.ExpBy.ns)
  cv.ExpBy.cs<-c(rep(1,N-N.exp-(YEAR-expYe)),ExpBy$expBy32.cv2[ExpBy$mngmt.area=="CSEO"],rep(1,(YEAR-expYe))); length(cv.ExpBy.cs)
  cv.ExpBy.ss<-c(rep(1,N-N.exp-(YEAR-expYe)),ExpBy$expBy32.cv2[ExpBy$mngmt.area=="SSEO"],rep(1,(YEAR-expYe))); length(cv.ExpBy.ss)
  
  ExpByc<-rbind(ExpBy.ey,ExpBy.ns,ExpBy.cs,ExpBy.ss); length(ExpByc)
  cv.ExpByc<-rbind(cv.ExpBy.ey,cv.ExpBy.ns,cv.ExpBy.cs,cv.ExpBy.ss); length(cv.ExpByc)
  
  #====================================================================================================
  #ROV Biomass Data
  #set up right now that biomass is not estimated in terminal year... 
  # thus the -1 in first list of NA,s and NA at the end of the data string
  Biomass<-read.csv(paste0("Data_processing/Data/SEO_YE_Biomass_subdistrict_",YEAR,".csv")); min(Biomass$Year); max(Biomass$Year)   #biomass estimates in metric tons
  Bio.first<-min(Biomass$Year)
  N.bio<-length(unique(Biomass$Year))
  bioY1<-min(Biomass$Year)
  bioYe<-max(Biomass$Year)
  
  B.obs.ey<-c(rep(NA,N-N.bio-(YEAR-bioYe)),Biomass$Biomass.mt[Biomass$Subdistrict == "EYKT"],rep(NA,(YEAR-bioYe))); length(B.obs.ey) 
  cv.B.obs.ey<-c(rep(NA,N-N.bio-(YEAR-bioYe)),Biomass$Biomass.cv[Biomass$Subdistrict == "EYKT"],rep(NA,(YEAR-bioYe))); length(cv.B.obs.ey) 
  B.obs.ns<-c(rep(NA,N-N.bio-(YEAR-bioYe)),Biomass$Biomass.mt[Biomass$Subdistrict == "NSEO"],rep(NA,(YEAR-bioYe))); length(B.obs.ns) 
  cv.B.obs.ns<-c(rep(NA,N-N.bio-(YEAR-bioYe)),Biomass$Biomass.cv[Biomass$Subdistrict == "NSEO"],rep(NA,(YEAR-bioYe))); length(cv.B.obs.ns) 
  B.obs.cs<-c(rep(NA,N-N.bio-(YEAR-bioYe)),Biomass$Biomass.mt[Biomass$Subdistrict == "CSEO"],rep(NA,(YEAR-bioYe))); length(B.obs.cs) 
  cv.B.obs.cs<-c(rep(NA,N-N.bio-(YEAR-bioYe)),Biomass$Biomass.cv[Biomass$Subdistrict == "CSEO"],rep(NA,(YEAR-bioYe))); length(cv.B.obs.cs) 
  B.obs.ss<-c(rep(NA,N-N.bio-(YEAR-bioYe)),Biomass$Biomass.mt[Biomass$Subdistrict == "SSEO"],rep(NA,(YEAR-bioYe))); length(B.obs.ss) 
  cv.B.obs.ss<-c(rep(NA,N-N.bio-(YEAR-bioYe)),Biomass$Biomass.cv[Biomass$Subdistrict == "SSEO"],rep(NA,(YEAR-bioYe))); length(cv.B.obs.ss) 
  
  #trying without blank for whole BS... 
  B.obs<-rbind(B.obs.ey, B.obs.ns, B.obs.cs, B.obs.ss); length(B.obs)
  cv.B<-rbind(cv.B.obs.ey, cv.B.obs.ns, cv.B.obs.cs, cv.B.obs.ss); length(cv.B)
  #cv.B[,(N.pre+1):N][is.na(cv.B[,(N.pre+1):N])]<-1
  cv.B[is.na(cv.B)]<-1
  #===============================================================
  # Get rid of '94???
  #B.obs[2:4,15]<-NA
  #cv.B[2:4,15]<-1
  #==========================================================
  #IPHC Survey CPUE RPN; 1998-2020
  IPHC.cpue<-read.csv(paste0("Data_processing/Data/IPHC.cpue.SEO_non0_",YEAR,".csv")); min(IPHC.cpue$Year); max(IPHC.cpue$Year)
  N.iphc<-length(unique(IPHC.cpue$Year))
  ipY1<-min(IPHC.cpue$Year)
  ipYe<-max(IPHC.cpue$Year)
  #IPHC.cpue$cv.cpue<-sqrt(IPHC.cpue$YE_CPUE.var)/IPHC.cpue$YE_CPUE
  
  ey.scpue<-c(rep(NA,N-N.iphc-(YEAR-ipYe)),IPHC.cpue$CPUE.mean[IPHC.cpue$mngmt.area == "EYKT"],rep(NA,(YEAR-ipYe))); length(ey.scpue)
  ns.scpue<-c(rep(NA,N-N.iphc-(YEAR-ipYe)),IPHC.cpue$CPUE.mean[IPHC.cpue$mngmt.area == "NSEO"],rep(NA,(YEAR-ipYe))); length(ns.scpue)
  cs.scpue<-c(rep(NA,N-N.iphc-(YEAR-ipYe)),IPHC.cpue$CPUE.mean[IPHC.cpue$mngmt.area == "CSEO"],rep(NA,(YEAR-ipYe))); length(cs.scpue)
  ss.scpue<-c(rep(NA,N-N.iphc-(YEAR-ipYe)),IPHC.cpue$CPUE.mean[IPHC.cpue$mngmt.area == "SSEO"],rep(NA,(YEAR-ipYe))); length(ss.scpue)
  
  cv.ey.scpue<-c(rep(1,N-N.iphc-(YEAR-ipYe)),IPHC.cpue$CPUE.cv[IPHC.cpue$mngmt.area == "EYKT"],rep(1,(YEAR-ipYe))); length(cv.ey.scpue)
  cv.ns.scpue<-c(rep(1,N-N.iphc-(YEAR-ipYe)),IPHC.cpue$CPUE.cv[IPHC.cpue$mngmt.area == "NSEO"],rep(1,(YEAR-ipYe))); length(cv.ns.scpue)
  cv.cs.scpue<-c(rep(1,N-N.iphc-(YEAR-ipYe)),IPHC.cpue$CPUE.cv[IPHC.cpue$mngmt.area == "CSEO"],rep(1,(YEAR-ipYe))); length(cv.cs.scpue)
  cv.ss.scpue<-c(rep(1,N-N.iphc-(YEAR-ipYe)),IPHC.cpue$CPUE.cv[IPHC.cpue$mngmt.area == "SSEO"],rep(1,(YEAR-ipYe))); length(cv.ss.scpue)
  
  sCPUE<-rbind(ey.scpue,ns.scpue,cs.scpue,ss.scpue)
  sCPUE<-sCPUE*1000000; length(sCPUE)             #scale CPUE to biomass level numbers
  cv.sCPUE<-rbind(cv.ey.scpue,cv.ns.scpue, cv.cs.scpue, cv.ss.scpue); length(cv.sCPUE)
  
  #===============================================================================
  #Directed DSR Fishery CPUE data
  #Fish.cpue<-read.csv("Data/YE_dirfishery_nom_CPUE.csv")
  #min(Fish.cpue$Year); max(Fish.cpue$Year)
  #N.fish<-length(unique(Fish.cpue$Year))
  #N.fish.closed<-YEAR-max(Fish.cpue$Year)
  #str(Fish.cpue)
  #View(Fish.cpue)
  #directed fishery has been closed and many years where there was none... 
  #Years
  #ey.fcpue<-c(rep(NA,N.subd-N.fish-N.fish.closed),Fish.cpue$fsh_cpue[Fish.cpue$Subd == "EYKT"],
  #            rep(NA,N.fish.closed)); length(ey.fcpue)
  #ns.fcpue<-c(rep(NA,N.subd-N.fish-N.fish.closed),Fish.cpue$fsh_cpue[Fish.cpue$Subd == "NSEO"],
  #            rep(NA,N.fish.closed)); length(ns.fcpue)
  #cs.fcpue<-c(rep(NA,N.subd-N.fish-N.fish.closed),Fish.cpue$fsh_cpue[Fish.cpue$Subd == "CSEO"],
  #            rep(NA,N.fish.closed)); length(cs.fcpue)
  #ss.fcpue<-c(rep(NA,N.subd-N.fish-N.fish.closed),Fish.cpue$fsh_cpue[Fish.cpue$Subd == "SSEO"],
  #            rep(NA,N.fish.closed)); length(ss.fcpue)
  #ss.fcpue[32]; Years.cont[32]
  
  #need cv's of 1 during 80+ years for model... 
  #cv.ey.fcpue<-c(rep(1,N-N.fish-N.fish.closed-N.pre),
  #            Fish.cpue$cv[Fish.cpue$Subd == "EYKT"],
  #            rep(1,N.fish.closed)); length(cv.ey.fcpue)
  #ey.fcpue[33]; Years.cont[33]; cv.ey.fcpue[33]
  #cv.ns.fcpue<-c(rep(1,N-N.fish-N.fish.closed-N.pre),
  #               Fish.cpue$cv[Fish.cpue$Subd == "NSEO"],
  #               rep(1,N.fish.closed)); length(cv.ns.fcpue)
  #cv.cs.fcpue<-c(rep(1,N-N.fish-N.fish.closed-N.pre),
  #               Fish.cpue$cv[Fish.cpue$Subd == "CSEO"],
  #               rep(1,N.fish.closed)); length(cv.cs.fcpue)
  #cv.ss.fcpue<-c(rep(1,N-N.fish-N.fish.closed-N.pre),
  #               Fish.cpue$cv[Fish.cpue$Subd == "SSEO"],
  #               rep(1,N.fish.closed)); length(cv.ss.fcpue)
  
  #fCPUE<-rbind(ey.fcpue,ns.fcpue,cs.fcpue,ss.fcpue)
  #fCPUE<-fCPUE*100000; length(fCPUE)             #scale CPUE to biomass level numbers
  #cv.fCPUE<-rbind(cv.ey.fcpue,cv.ns.fcpue, cv.cs.fcpue, cv.ss.fcpue); length(cv.fCPUE)
  
  #*******************************************************************************
  #* SPECIAL DATA :)
  #* Beta parameters and Derby effect
  #* # can be changed in MULTIMODEL RUN for ease, but need to load base stuff here...
  #================================================================================
  # END DATA LOAD for 2022 SAFE
  #Below are other data sources for different model structures.
  Derby.Eff <- Derby.Eff
  DEsd<-DEsd
  B1<-B1
  B2<-B2
  
  Data<-list(Subd=Subd,N=N,Years=Years,
             YEAR=YEAR,
             KnC.obs=KnC.obs,cv.KnC=cv.KnC,
             Lnd.By=Lnd.By,
             ExpByc=ExpByc,cv.ExpByc=cv.ExpByc,
             B.obs=B.obs,cv.B=cv.B,
             sCPUE=sCPUE,cv.sCPUE=cv.sCPUE,#fCPUE=fCPUE,cv.fCPUE=cv.fCPUE,
             Derby.Eff = Derby.Eff, DEsd=DEsd,
             B1=B1,B2=B2)
  return(Data)
}

################################################################################
################################################################################
# Function to package the data for runnning the model, including initial values

data.prep<-function(){
  
  data<-list(N=N,Subd=Subd,
             KnC.obs=KnC.obs, cv.KnC=cv.KnC,
             ExpByc = ExpByc, cv.ExpByc = cv.ExpByc,
             Lnd.By=Lnd.By, 
             B.obs=B.obs,cv.B=cv.B, #CKP=CKP,
             sCPUE=sCPUE, cv.sCPUE=cv.sCPUE, #,
             #fCPUE=sCPUE, cv.fCPUE=cv.fCPUE, 
             Fu=Fu, 
             OFLABC=OFLABC, #B.err=B.err,
             HarvStrat=HarvStrat,
             Derby.Eff = Derby.Eff, DEsd=DEsd,
             B1=B1,B2=B2,upvar=upvar
  )
  
  KnCinits1<-log(KnC.obs+runif(N*Subd,-0.5*KnC.obs,0.5*KnC.obs))
  KnCinits2<-log(KnC.obs+runif(N*Subd,-0.5*KnC.obs,0.5*KnC.obs))
  KnCinits3<-log(KnC.obs+runif(N*Subd,-0.5*KnC.obs,0.5*KnC.obs))
  KnCinits4<-log(KnC.obs+runif(N*Subd,-0.5*KnC.obs,0.5*KnC.obs))
  KnCinits5<-log(KnC.obs+runif(N*Subd,-0.5*KnC.obs,0.5*KnC.obs))
  
  Kinits1<-c(log(20000),log(20000),log(20000),log(20000))
  Kinits2<-c(log(30000),log(10000),log(8000),log(25000))
  Kinits3<-c(log(15000),log(6000),log(25000),log(5000))
  
  phiinits1<-log(c(0.8,0.7,0.45,0.8))
  phiinits2<-log(c(0.6,0.5,0.7,0.6))
  phiinits3<-log(c(0.5,0.8,0.6,0.43))
  
  qsinits1<-log(runif(4,0.2,10))#c(3,4,0.5,0.3)
  qsinits2<-log(runif(4,0.2,10))#c(3,10,18,15)
  qsinits3<-log(runif(4,0.2,10))#c(10,10,12,12)
  
  #qfinits1<-log(runif(4,0.2,10))#c(3,4,5,10)
  #qfinits2<-log(runif(4,0.2,10))#c(3,10,18,15)
  #qfinits3<-log(runif(4,0.2,10))#c(10,10,12,12)
  
  inits1<-list(R.hyp=0.04, logK=Kinits1, #phi=phiinits1,
               logphi=phiinits1, logvar=-9,
               logqs=qsinits1, #logqf=qfinits1, 
               Tau1=0.4, #Tau3=0.03, Tau4 = 0.3, 
               Tau3=0.01, 
               logKnC=KnCinits1)
  inits2<-list(R.hyp=0.07, logK=Kinits2, #phi=phiinits2, 
               logphi=phiinits2, logvar = -8,
               logqs=qsinits2, #logqf=qfinits2,  
               Tau1=0.12, #Tau3=0.1, Tau4 = 0.2,
               Tau3=0.3, 
               logKnC=KnCinits2)
  inits3<-list(R.hyp=0.02, logK=Kinits3, #phi=phiinits3, 
               logphi=phiinits3, logvar = -7,
               logqs=qsinits3, #logqf=qfinits3, 
               Tau1=0.2, #Tau3=0.13, Tau4=0.02,
               Tau3=0.5, 
               logKnC=KnCinits3)
  
  #inits <- list(inits1,inits2,inits3,inits4,inits5)
  inits<-list(inits1,inits2,inits3)
  
  parameters=c(
    "r","K","phi", "R.hyp",#"R.sig",#"phi.hyp", 
    "Kseo", 
    "sigma","eta","rB1","rB2",
    "epsilon","PE",
    "B", "D","C","FC","KnC", "By",
    "Bseo","Cseo","Byseo","Dseo","KnCseo",
    "S", "Surplus", "S.seo","Surplus.seo",
    "qsCPUE", #"qfCPUE",
    "Tau1" ,"Tau2", "Tau3", #"Tau4",
    "CtoB","CBtoK","CBtoKseo","FBtoK","FBtoKseo","FBtoCB",
    "MSY","Bmsy","Fmsy","Hmsy","Stock.Status",
    "Bmsyseo","Fmsyseo","Stock.Status.SEO",
    "fit","subfit.new","iphcfit.new"
  )
  
  prep<-list(data=data,inits=inits,parameters=parameters)
  return(prep)
}

################################################################################
################################################################################
# Plotting function for quick plotting of SPM results:

plot_spm<-function(post=post, Mod.title="test",filename){
  #tbl<-jags.View(post, title="", digits=3)
  sdlist<-c("EYKT", "NSEO", "CSEO", "SSEO")
  tbl<-as.data.frame(jags.View(post, title="", digits=3))
  write.csv(tbl,
            file=paste("Production_models/Output/", filename,"/results_sum.csv", sep=""))
  
  Bseo.sum<-tbl[grepl("Bseo",tbl$parameter),]
  Bseo.sum$cv<-Bseo.sum$sd/Bseo.sum$mean
  write.csv(Bseo.sum,file=paste("Production_models/Output/",filename,"/Bseo_post_", Mod.title,".csv", sep=""))
  
  png(paste("Production_models/Figures/", filename,"/PP_check_submersibles.png",sep=""),
      width=6,height=6,#width=9.5,height=8.5,
      units="in",res=1200)
  par(mfrow=c(2,2), omi=c(0,0,0.1,0.01))
  
  for (i in 1:4) {
    obs<-paste("fit[",i,"]", sep="")
    sim<-paste("subfit.new[",i,"]", sep="")
    pp.check(post, observed=obs, simulated=sim, 
             xlab='log Observed biomass', ylab='log Simulated biomass',
             main=paste(sdlist[i],""), col="grey")
  }
  mtext("Posterior predictive check for submersible/ROV biomass", side=3,  
        padj=1,outer=TRUE, cex=1.4)
  dev.off()
  
  png(paste("Production_models/Figures/", filename,"/PP_check_iphc.png",sep=""),
      width=6,height=6,#width=9.5,height=8.5,
      units="in",res=1200)
  par(mfrow=c(2,2), omi=c(0,0,0.1,0.01))
  for (i in 1:4) {
    obs<-paste("fit[",i,"]", sep="")
    sim<-paste("iphcfit.new[",i,"]", sep="")
    pp.check(post, observed=obs, simulated=sim, 
             xlab='log Observed biomass', ylab='log Simulated biomass',
             main=paste(sdlist[i],""), col="grey")
  }
  mtext("Posterior predictive check for IPHC CPUE", side=3, cex=1.4, 
        padj=1,outer=TRUE)
  dev.off()
  #===============================================================
  #Posterior of all biomass estimates
  {
    png(paste("Production_models/Figures/", filename,"/Biomass_overview.png",sep=""),
        width=7,height=6,#width=9.5,height=8.5,
        units="in",res=1200)
    B.ey<-par.ext(par="B",years=N.subd+1,areai=1)
    B.ns<-par.ext("B",N.subd+1,2)
    B.cs<-par.ext("B",N.subd+1,3)
    B.ss<-par.ext("B",N.subd+1,4)
    B.all<-MCMCchains(post,"Bseo")
    B.all2<-B.all[,1:N.subd]
    #Bs<-Bs[,1:42]
    par(mfrow=c(1,1))
    #envplot(B.all,1980:(1980+N.subd-1),cols=c("goldenrod3","goldenrod3"),n=0,ylab="Biomass (t)",xlab="Year",ylim=c(0,40000))
    envplot(B.all2,1980:(1980+N.subd),cols=c("goldenrod3","goldenrod3"),n=0,ylab="Biomass (t)",xlab="Year",ylim=c(0,40000))
    Y<-seq(1980,(1980+N.subd+1-1),1)
    addenv(B.ey,cols=c("blue","blue"), Y=Y)
    addenv(B.ns, cols=c("red","red"), Y=Y)
    addenv(B.cs, cols=c("forestgreen","forestgreen"), Y=Y)
    addenv(B.ss, cols=c("purple","purple"), Y=Y)
    legend(x="topright", c("Entire SEO","EYKT","NSEO","CSEO","SSEO"), 
           col=c("goldenrod","blue","red","forestgreen","purple"),
           text.col=c("goldenrod","blue","red","forestgreen","purple"),
           border=F, bty="n")
    dev.off()
  }
  
  #i) Plotting setup
  #Years.cont
  All.years<-seq(min(Years.cont),max(Years.cont)+1,by=1)#seq(min(Years.cont),max(Years.cont)+Fu,by=1)
  FuYear<-Years.cont+1 #Years.cont+Fu
  sdlist<-c("EYKT", "NSEO", "CSEO", "SSEO")
  
  #1) Biomass observed versus modeled 
  colch<-c("blue","red","purple","orange","forestgreen")
  vars<-list(B.ey,B.ns,B.cs,B.ss)
  
  png(paste("Production_models/Figures/", filename,"/Biomass_fitb.png",sep=""),
      width=7,height=6,#width=9.5,height=8.5,
      units="in",res=1200)
  par(mfrow=c(2,2), mar=c(4,5,3,1))
  for (i in 1:4){  #i<-4
    envplot(vars[[i]],All.years,cols=c(colch[1],colch[1]),n=0,ylab="Biomass (t)",xlab="Year",
            ylim=c(0,17500),main=sdlist[i], cex.axis=0.8, cex.lab=1.25, cex.main=1.5)
    #addpoints(post,Years=Years.cont+0.2 ,Points=fCPUE[i,] ,errordat=cv.fCPUE[i,] ,scaler="qfCPUE" , #scaler="qfCPUE" for example...
    #                      bar=1.95 ,col=colch[5], pch=18, cex=0.75)
    #points(fCPUE[i,]/3.808~Years.cont, col="green")
    #points(sCPUE[i,]/3.821~Years.cont,col="darkorange")
    addpoints(post,Years=Years.cont+0.1 ,Points=sCPUE[i,], errordat=cv.sCPUE[i,], scaler="qsCPUE", #scaler="qfCPUE" for example...
              bar=1.95 ,col=colch[4], pch=17, cex=0.75)
    addpoints(post,Years=Years.cont ,Points=B.obs[i,], errordat=cv.B[i,], scaler=0, #scaler="qfCPUE" for example...
              bar=1.95 ,col=colch[2], pch=18, cex=1.2)
    #3rd index  addpoints(post,Years=Years.cont+0.2 ,Points=B.obs[i,] ,errordat=cv.B[i,] ,scaler=0 , #scaler="qfCPUE" for example...
    #            bar=1.95 ,col=colch[5], pch=18)
    #addenv(vars[[i]], cols=c(colch[1],colch[1]),alpha.f=0.05)
    abline(v=2022,  col=colch[3], lty=2)
    if (i == 2){
      legend(x="topright", c("posterior est.","observed ROV/sub","IPHC CPUE"), col=c(colch[1],colch[2],colch[4]),
             text.col=c(colch[1],colch[2],colch[4]),border=F, bty="y", pch=c(NA,18,17), cex=1)
    } else {}
  }
  dev.off()
  
  #1.B) Biomass observed versus modeled with Bmsy lines
  {colch<-c("blue","red","purple","orange","forestgreen")
    vars<-list(B.ey,B.ns,B.cs,B.ss)
    
    Bmsy.ey.p<-as.vector(MCMCchains(post,"Bmsy[1]",ISB = FALSE,  #ignore square brackets
                                    mcmc.list = FALSE,chain_num = NULL, exact=TRUE))
    #str(Bmsy.ey.p); quantile(Bmsy.ey.p,p=c(0.25))
    Bmsy.ns.p<-as.vector(MCMCchains(post,"Bmsy[2]",ISB = FALSE,  #ignore square brackets
                                    mcmc.list = FALSE,chain_num = NULL, exact=TRUE))
    Bmsy.cs.p<-as.vector(MCMCchains(post,"Bmsy[3]",ISB = FALSE,  #ignore square brackets
                                    mcmc.list = FALSE,chain_num = NULL, exact=TRUE))
    Bmsy.ss.p<-as.vector(MCMCchains(post,"Bmsy[4]",ISB = FALSE,  #ignore square brackets
                                    mcmc.list = FALSE,chain_num = NULL, exact=TRUE))
    
    Bmsy.ey<-tbl$mean[tbl$parameter=="Bmsy[1]"]
    Bmsy.ey_lo<-tbl$`2.5%`[tbl$parameter=="Bmsy[1]"]
    Bmsy.ey_hi<-tbl$`97.5%`[tbl$parameter=="Bmsy[1]"]
    
    Bmsy.ns<-tbl$mean[tbl$parameter=="Bmsy[2]"]
    Bmsy.ns_lo<-tbl$`2.5%`[tbl$parameter=="Bmsy[2]"]
    Bmsy.ns_hi<-tbl$`97.5%`[tbl$parameter=="Bmsy[2]"]
    
    Bmsy.cs<-tbl$mean[tbl$parameter=="Bmsy[3]"]
    Bmsy.cs_lo<-tbl$`2.5%`[tbl$parameter=="Bmsy[3]"]
    Bmsy.cs_hi<-tbl$`97.5%`[tbl$parameter=="Bmsy[3]"]
    
    Bmsy.ss<-tbl$mean[tbl$parameter=="Bmsy[4]"]
    Bmsy.ss_lo<-tbl$`2.5%`[tbl$parameter=="Bmsy[4]"]
    Bmsy.ss_hi<-tbl$`97.5%`[tbl$parameter=="Bmsy[4]"]
    
    var1<-list(Bmsy.ey,Bmsy.ns,Bmsy.cs,Bmsy.ss)
    var2<-list(Bmsy.ey_lo,Bmsy.ns_lo,Bmsy.cs_lo,Bmsy.ss_lo)
    var3<-list(Bmsy.ey_hi,Bmsy.ns_hi,Bmsy.cs_hi,Bmsy.ss_hi)
    var4<-list(Bmsy.ey.p,Bmsy.ns.p,Bmsy.cs.p,Bmsy.ss.p)
  }
  png(paste("Production_models/Figures/", filename,"/Biomass_wBmsy.png",sep=""),
      width=7,height=6,#width=9.5,height=8.5,
      units="in",res=1200)
  par(mfrow=c(2,2), mar=c(4,5,3,1))
  
  for (i in 1:4){  #i<-1
    envplot(vars[[i]],All.years,cols=c("blue","blue"),n=0,ylab="Biomass (t)",xlab="Year",
            ylim=c(0,17500),main=sdlist[i], cex.axis=0.8, cex.lab=1.25, cex.main=1.5)
    polygon(x=c(All.years,rev(All.years)),
            y=c(rep(quantile(as.numeric(var4[[i]]),p=c(0.25)),length(All.years)),
                rep(quantile(as.numeric(var4[[i]]),p=c(0.75)),length(All.years))),
            border=NA,
            col=adjustcolor("grey", alpha.f=.75))
    polygon(x=c(All.years,rev(All.years)),
            y=c(rep(quantile(as.numeric(var4[[i]]),p=c(0.025)),length(All.years)),
                rep(quantile(as.numeric(var4[[i]]),p=c(0.975)),length(All.years))),
            border=NA,
            col=adjustcolor("grey", alpha.f=.50))
    abline(h=quantile(as.numeric(var4[[i]]),p=c(0.5)),col="black", lwd=1.8, lty=2)
    addenv(vars[[i]],cols=c("blue","blue"), Y=All.years)
    
    #abline(h=var1[i],col=colch[5], lty=1, lwd=1)
    
    if (i == 2){
      legend(x="topright", c("posterior est.","Bmsy"), col=c("blue","black"),
             text.col=c("blue","black"), border=F, bty="y", pch=c(NA,NA,NA,NA,NA), cex=1)
    } else {}
  }
  dev.off()
  
  #2) Known catch
  {KnC.ey<-par.ext(par="KnC",years=N.subd,areai=1)
    KnC.ns<-par.ext("KnC",N.subd,2)
    KnC.cs<-par.ext("KnC",N.subd,3)
    KnC.ss<-par.ext("KnC",N.subd,4)
    
    colch<-c("blue","forestgreen","purple")
    vars<-list(KnC.ey,KnC.ns,KnC.cs,KnC.ss)}
  
  png(paste("Figures/", filename,"/KnownCatch_fit.png",sep=""),
      width=7,height=6,#width=9.5,height=8.5,
      units="in",res=1200)
  par(mfrow=c(2,2), mar=c(4,5,3,1))
  for (i in 1:4){  #i<-1
    envplot(vars[[i]],Years.cont,cols=c(colch[1],colch[1]),n=0,ylab="Known catches (t)",xlab="Year",
            ylim=c(0,max(KnC.obs*1.25)),main=sdlist[i],cex.axis=0.8, cex.lab=1.25, cex.main=1.5)
    addpoints(post,Years=Years.cont ,Points=KnC.obs[i,] ,errordat=cv.KnC[i,] ,scaler=0 ,
              bar=1.95 ,col=colch[2], pch=18,cex=1.2)
    abline(v=2022,  col=colch[3], lty=2)
    if (i == 2){
      legend(x="topright", c("posterior est.","observed"), col=c(colch[1],colch[2]),
             text.col=c(colch[1],colch[2]),border=F, bty="n", pch=c(NA,18))
    } else {}
  }
  dev.off()
  
  #3) Expected, landed and modeled bycatch... 
  {By.ey<-par.ext(par="By",years=N.subd,areai=1)
    By.ns<-par.ext("By",N.subd,2)
    By.cs<-par.ext("By",N.subd,3)
    By.ss<-par.ext("By",N.subd,4)
    
    colch<-c("blue","orange","purple")
    vars<-list(By.ey,By.ns,By.cs,By.ss)}
  
  png(paste("Production_models/Figures/", filename,"/Bycatch_fit.png",sep=""),
      width=7,height=6,#width=9.5,height=8.5,
      units="in",res=1200)
  par(mfrow=c(2,2), mar=c(4,5,3,1))
  for (i in 1:4){  #i<-1
    envplot(vars[[i]],Years.cont,cols=c(colch[1],colch[1]),n=0,ylab="Bycatch (t)",xlab="Year",
            ylim=c(0,max(ExpByc*1.25)),main=sdlist[i],cex.axis=0.8, cex.lab=1.25, cex.main=1.5)
    addpoints(post,Years=Years.cont ,Points=ExpByc[i,] ,errordat=cv.ExpByc[i,] ,scaler=0 ,
              bar=1.95 ,col=colch[2], pch=18, cex=1.2)
    lines(Years.cont,Lnd.By[i,],type="b",col=colch[3], lwd=1.5, lty=2, pch=15, cex=0.5)
    #abline(v=2022,  col=colch[3], lty=2)
    if (i == 2){
      legend(x="topright", c("posterior est.","observed 'expected'","landed bycatch"), col=c(colch[1],colch[2],colch[3]),
             text.col=c(colch[1],colch[2],colch[3]),border=F, bty="n", pch=c(NA,18,15),cex=1)
    } else {}
  }
  dev.off()
  
  #4) Discard plots
  {D.ey<-par.ext(par="D",years=N.subd,areai=1)
    D.ns<-par.ext("D",N.subd,2)
    D.cs<-par.ext("D",N.subd,3)
    D.ss<-par.ext("D",N.subd,4)
    
    colch<-c("blue","orange","purple")
    vars<-list(D.ey,D.ns,D.cs,D.ss)}
  
  png(paste("Production_models/Figures/", filename,"/Discards.png",sep=""),
      width=7,height=6,#width=9.5,height=8.5,
      units="in",res=1200)
  par(mfrow=c(2,2), mar=c(4,5,3,1))
  for (i in 1:4){  #i<-1
    envplot(vars[[i]],Years.cont,cols=c(colch[1],colch[1]),n=0,ylab="Estimated Discards (t)",xlab="Year",
            ylim=c(0,max(ExpByc*1.25)),main=sdlist[i],cex.axis=0.8, cex.lab=1.25, cex.main=1.5)
    #  addpoints(post,Years=Years.cont ,Points=ExpByc[i,] ,errordat=cv.ExpByc[i,] ,scaler=0 ,
    #            bar=1.95 ,col=colch[2], pch=18)
  }
  dev.off()
  
  #5) Catch plots of total, known, bycatch, landed bycatch, and discards
  {
    png(paste("Production_models/Figures/", filename,"/Catch_breakdown_subd.png",sep=""),
        width=7,height=6,#width=9.5,height=8.5,
        units="in",res=1200)
    par(mfrow=c(2,2), mar=c(4,5,3,1))
    for (i in 1:4){
      Cg<-par.ext(par="C",years=N.subd,areai=i)  
      Byg<-par.ext("By",N.subd,i)
      Dg<-par.ext("D",N.subd,i)
      KnCg<-par.ext("KnC",N.subd,i)
      envplot(Cg,1980:(1980+N.subd-1),cols=c("blue","blue"),n=0,ylab="Catches (t)",xlab="Year",
              ylim=c(0,400), main=sdlist[i],cex.axis=0.8, cex.lab=1.25, cex.main=1.5)
      Y<-seq(1980,(1980+N.subd-1),1)
      addenv(Byg,cols=c("red","red"), Y=Y)
      lines(Lnd.By[i,]~Y, col="darkred", lwd=2, lty=2)
      addenv(Dg, cols=c("black","black"), Y=Y)
      addenv(KnCg, cols=c("forestgreen","forestgreen"), Y=Y)
      if (i == 2) {
        legend(x="topright", c("Total removals","Know catches","Expected bycatch (IPHC)","Landed Bycatch","Discards"),
               col=c("blue","forestgreen","red","darkred","darkgrey"),
               text.col=c("blue","forestgreen","red","darkred","black"),
               border=F, bty="n", pch=c(NA),cex=1)
      } else {}
    }
    dev.off()
  }
  
  
  {
    par(mfrow=c(1,1))
    Cg<-MCMCchains(post,"Cseo")
    Byg<-MCMCchains(post,"Byseo")
    Dg<-MCMCchains(post,"Dseo")
    KnCg<-MCMCchains(post,"KnCseo")
    png(paste("Production_models/Figures/", filename,"/Catch_breakdown.png",sep=""),
        width=7,height=6,#width=9.5,height=8.5,
        units="in",res=1200)
    par(mfrow=c(1,1))
    envplot(Cg,1980:(1980+N.subd-1),cols=c("blue","blue"),n=0,ylab="Catches (t)",xlab="Year",ylim=c(0,1000),
            cex.axis=1, cex.lab=1.25)
    Y<-seq(1980,(1980+N.subd-1),1)
    addenv(Byg,cols=c("red","red"), Y=Y)
    Lndby.tot<-Lnd.By[1,]+ Lnd.By[2,]+Lnd.By[3,]+Lnd.By[4,]#Lnd.ey+Lnd.ns+Lnd.cs+Lnd.ss
    lines(Lndby.tot~Y, col="darkred", lwd=2, lty=2)
    addenv(Dg, cols=c("black","black"), Y=Y)
    addenv(KnCg, cols=c("forestgreen","forestgreen"), Y=Y)
    legend(x="topright", c("Total removals","Know catches","Expected bycatch (IPHC)","Landed Bycatch","Discards"),
           col=c("blue","forestgreen","red","darkred","darkgrey"),
           text.col=c("blue","forestgreen","red","darkred","black"),
           border=F, bty="n", pch=c(NA),cex=1)
    dev.off()
    
  }
  
  #6) Surplus production estimates
  {S.ey<-par.ext(par="Surplus",years=N.subd,areai=1)
    S.ns<-par.ext("Surplus",N.subd,2)
    S.cs<-par.ext("Surplus",N.subd,3)
    S.ss<-par.ext("Surplus",N.subd,4)
    
    colch<-c("forestgreen","black","purple")
    vars<-list(S.ey,S.ns,S.cs,S.ss)}
  
  png(paste("Production_models/Figures/", filename,"/Surplus.png",sep=""),
      width=7,height=6,#width=9.5,height=8.5,
      units="in",res=1200)
  par(mfrow=c(2,2), mar=c(4,5,3,1))
  for (i in 1:4){  #i<-1
    envplot(vars[[i]],Years.cont,cols=c(colch[1],colch[1]),n=0,ylab="Surplus Production (t)",xlab="Year",
            ylim=c(0,300),main=sdlist[i],cex.axis=0.8, cex.lab=1.25, cex.main=1.5)
    addpoints(post,Years=Years.cont ,Points=KnC.obs[i,] ,errordat=cv.KnC[i,] ,scaler=0 ,
              bar=1.95 ,col=colch[2], pch=18,cex=1)
    if (i == 2){
      legend(x="topright", c("Surplus","Known Catches"), col=c(colch[1],colch[2],colch[3]),
             text.col=c(colch[1],colch[2],colch[3]),border=F, bty="n", pch=c(NA,18,15),cex=1)
    } else {}
  }
  dev.off()
  
  #7) CtoB
  {CtoB.ey<-par.ext(par="CtoB",years=N.subd,areai=1)
    CtoB.ns<-par.ext("CtoB",N.subd,2)
    CtoB.cs<-par.ext("CtoB",N.subd,3)
    CtoB.ss<-par.ext("CtoB",N.subd,4)
    
    Fmsy.ey<-tbl$mean[tbl$parameter=="Fmsy[1]"]#posteriors$Fmsy[1,4]#par.ext("Fmsy")
    Fmsy.ns<-tbl$mean[tbl$parameter=="Fmsy[2]"]#posteriors$Fmsy[2,4]
    Fmsy.cs<-tbl$mean[tbl$parameter=="Fmsy[3]"]#posteriors$Fmsy[3,4]
    Fmsy.ss<-tbl$mean[tbl$parameter=="Fmsy[4]"]#posteriors$Fmsy[4,4]
    
    Flo.ey<-tbl$`2.5%`[tbl$parameter=="Fmsy[1]"]; Fhi.ey<-tbl$`97.5%`[tbl$parameter=="Fmsy[1]"]
    Flo.ns<-tbl$`2.5%`[tbl$parameter=="Fmsy[2]"]; Fhi.ns<-tbl$`97.5%`[tbl$parameter=="Fmsy[2]"]
    Flo.cs<-tbl$`2.5%`[tbl$parameter=="Fmsy[3]"]; Fhi.cs<-tbl$`97.5%`[tbl$parameter=="Fmsy[3]"]
    Flo.ss<-tbl$`2.5%`[tbl$parameter=="Fmsy[4]"]; Fhi.ss<-tbl$`97.5%`[tbl$parameter=="Fmsy[4]"]
    
    colch<-c("black","yellow2","orange")
    vars<-list(CtoB.ey,CtoB.ns,CtoB.cs,CtoB.ss); vars2<-list(Fmsy.ey,Fmsy.ns,Fmsy.cs,Fmsy.ss)
    varlo<-list(Flo.ey,Flo.ns,Flo.cs,Flo.ss); varup<-list(Fhi.ey,Fhi.ns,Fhi.cs,Fhi.ss)}
  
  png(paste("Production_models/Figures/", filename,"/CtoB.png",sep=""),
      width=7,height=6,#width=9.5,height=8.5,
      units="in",res=1200)
  par(mfrow=c(2,2), mar=c(4,5,3,1))
  for (i in 1:4){  #i<-1
    envplot(vars[[i]],Years.cont,cols=c(colch[1],colch[1]),n=0,ylab="Catch:Biomass",xlab="Year",
            ylim=c(0,0.08),main=sdlist[i],cex.axis=0.8, cex.lab=1.25, cex.main=1.5)
    abline(h=vars2[i],col=colch[3], lty=1, lwd=1.5)
    abline(h=varlo[i],col=colch[3], lty=2, lwd=1)
    abline(h=varup[i],col=colch[3], lty=2, lwd=1)
    polygon(x=c(Years.cont,rev(Years.cont)),
            y=c(rep(varlo[i],length(Years.cont)),rep(varup[i],length(Years.cont))),
            border=NA,
            col=adjustcolor(colch[2], alpha.f=.10))
    # addenv(vars[[i]], cols=c(colch[1],colch[1]))
    legend(x="topright", c("Catch:Biomass",paste0("Fmsy = ",vars2[i])), col=c(colch[1],colch[3],colch[3]),
           text.col=c(colch[1],colch[3],colch[3]),border=F, bty="n", pch=c(NA,NA,15),cex=1)
  }
  dev.off()
  
  #8) Epsilons/ Process Error
  png(paste("Production_models/Figures/", filename,"/epsilon.png",sep=""),
      width=5,height=4,#width=9.5,height=8.5,
      units="in",res=1200)
  #pe<-MCMCchains(post,"epsilon")
  pe2<-MCMCchains(post,"PE")
  pe<-pe2[,1:N.subd]
  par(mfrow=c(1,1))
  #envplot(pe,1980:(1980+N.subd),cols=c("black","black"),n=0,ylab="epsilon",xlab="Year",ylim=c(-0.2,0.2))
  #abline(h=0,col="red",lty=2)
  envplot(pe,1980:(1980+N.subd),cols=c("black","black"),n=0,ylab="process error",xlab="Year",ylim=c(-0.12,0.12))
  abline(h=0,col="red",lty=2)
  dev.off()
  
  #==============================================================================
  # RISK ANALYSIS: the probability that Fu years in the future population is at or
  # above Bmsy; 
  
  {
    Bmsy.ey <-MCMCchains(object=post,  params = c("Bmsyseo"),  excl = NULL,
                         ISB = FALSE,  #ignore square brackets
                         mcmc.list = FALSE,chain_num = NULL, exact=TRUE
    )#; nrow(Bmsy.ey)
    
    BCurrent<-MCMCchains(object=post,  params = c(paste("Bseo[",N.subd,"]",sep="")),  excl = NULL,
                         ISB = FALSE,  #ignore square brackets
                         mcmc.list = FALSE,chain_num = NULL, exact=TRUE
    )#; nrow(BCurrent)
    BFuture<-MCMCchains(object=post,  params = c(paste("Bseo[",N.subd+Fu,"]",sep="")),  excl = NULL,
                        ISB = FALSE,  #ignore square brackets
                        mcmc.list = FALSE,chain_num = NULL, exact=TRUE
    )#; nrow(BFuture)
    Risk<-cbind(Bmsy.ey,BCurrent,BFuture)
    
    colnames(Risk)<-c("Bmsy","Bcur","Bfut")
    Risk<-as.data.frame(Risk)
    Risk$current.stat<-Risk$Bcur - Risk$Bmsy
    Risk$future.stat<-Risk$Bfut - Risk$Bmsy
    Risk$increaseQ<-Risk$Bfut - Risk$Bcur
    
    CurrentRisk<-nrow(Risk[Risk$current.stat > 1,])/nrow(Risk)
    FutureRisk<-nrow(Risk[Risk$future.stat > 1,])/nrow(Risk)
    IncreaseProb<-nrow(Risk[Risk$increaseQ > 1,])/nrow(Risk)
    
    Risk.sum<-data.frame()
    i<-1
    for (i in 1:4){
      Bmsyi <-MCMCchains(object=post,  params = c(paste("Bmsy[",i,"]",sep="")),  excl = NULL,
                         ISB = FALSE,  #ignore square brackets
                         mcmc.list = FALSE,chain_num = NULL, exact=TRUE
      )
      BCur<-MCMCchains(object=post,  params = c(paste("B[",i,",",N.subd,"]",sep="")),  excl = NULL,
                       ISB = FALSE,  #ignore square brackets
                       mcmc.list = FALSE,chain_num = NULL, exact=TRUE
      )
      BFut<-MCMCchains(object=post,  params = c(paste("B[",i,",",N.subd+Fu,"]",sep="")),  excl = NULL,
                       ISB = FALSE,  #ignore square brackets
                       mcmc.list = FALSE,chain_num = NULL, exact=TRUE
      )
      Risk<-cbind(Bmsyi,BCur,BFut) #head(Risk)
      colnames(Risk)<-c("Bmsy","Bcur","Bfut")
      Risk<-as.data.frame(Risk)
      Risk$current.stat<-Risk$Bcur - Risk$Bmsy
      Risk$future.stat<-Risk$Bfut - Risk$Bmsy
      Risk$increaseQ<-Risk$Bfut - Risk$Bcur
      
      CR<-nrow(Risk[Risk$current.stat > 1,])/nrow(Risk)
      FR<-nrow(Risk[Risk$future.stat > 1,])/nrow(Risk)
      IP<-nrow(Risk[Risk$increaseQ > 1,])/nrow(Risk)
      
      Risk.sum[i,"Subd"]<-sdlist[i]
      Risk.sum[i,"Above_Bmsy_now"]<-CR
      Risk.sum[i,"Above_Bmsy_future"]<-FR
      Risk.sum[i,"Pop_inc_future"]<-IP
    }
    Risk.sum[5,"Subd"]<-"SEO"
    Risk.sum[5,"Above_Bmsy_now"]<-CurrentRisk
    Risk.sum[5,"Above_Bmsy_future"]<-FutureRisk
    Risk.sum[5,"Pop_inc_future"]<-IncreaseProb
  }
  write.csv(Risk.sum,file=paste("Production_models/Output/", filename,"/Risk_Summary_", Mod.title,".csv", sep=""))
  
}

################################################################################
################################################################################
## Stage 2 Dataloading function: 
load.data.1888<-function (YEAR=2022,biomassmod = "PT_2i_pe01_T1e_rbeta_B1-1.61_B2-67__derb_0_900k",
                          bioY1 = 1994, bioYe = 2021,
                          Derby.Eff = Derby.Eff, DEsd=DEsd, B1=B1,B2=B2){
  #=========================================================
  YEAR<-YEAR  
  ## 1) Halibut bycatch - longest time series; use this to establish number of years
  # Load expected bycatch generated in Code/Data_Generation/IPHC_Survey_Expected_Bycatch.R
  # Load expected bycatch generated in Code/Data_Generation/IPHC_Survey_Expected_Bycatch.R
  ExpBy<-read.csv(paste0("Data_processing/Data/SE_expBy_1888_",YEAR,".csv"))
  Year1<-min(ExpBy$Year)
  Yearmax<-max(ExpBy$Year)
  
  Years<-seq(min(ExpBy$Year),YEAR,by=1)
  N<-length(Years)  
  
  recavg<-mean(ExpBy$expBy32_mt.mean[(N-(YEAR-Yearmax)-4):(N-(YEAR-Yearmax))])
  
  ExpByc<-c(ExpBy$expBy32_mt.mean,rep(recavg,YEAR-Yearmax)); length(ExpByc)
  cv.ExpByc<-c(ExpBy$expBy32.cv,rep(1,YEAR-Yearmax)); length(cv.ExpByc)  
  
  #=========================
  ## Known Catches fo KnC.obs... 
  #1) contemporary removals from 1980 to present... 
  Removals<-read.csv(paste0("Data_processing/Data/SE_YE_known_removals_1980-",YEAR,".csv"))  #YE reconstruction from Rhea
  #2) pre-1980 removals from SEO as a whole... 
  # Removals.pre<-read.csv("Data/XXX.csv") #Waiting on Donnie; 
  # generate fake data here for model development...
  Removals %>% filter (Mgt.Area != "SEO" & Mgt.Area != "SSEI" & Mgt.Area != "NSEI") %>%
    group_by(Year) %>%
    summarise(tot.rem = sum(tot.rem.mt),
              tot.hal.by = sum(tot.hal.by))-> rem
  N.rem<-nrow(rem)
  remY1<-min(rem$Year)
  remYe<-max(rem$Year)
  
  #recavg2<-mean(rem$expBy32_mt.mean[(N-(YEAR-Yearmax)-4):(N-(YEAR-Yearmax))])
  
  KnC.obs<-c(rep(0.0001,N-N.rem-(YEAR-remYe)),rem$tot.rem,rep(NA,(YEAR-remYe))); length(KnC.obs)
  
  cv<-0.1
  cv.KnC<-c(rep(0.01,N-N.rem-(YEAR-remYe)),rep(cv,N.rem),rep(NA,(YEAR-remYe))); length(cv.KnC)
  
  #YE landings in Halibut fishery
  Lnd.By<-c(rep(0,N-N.rem-(YEAR-remYe)),rem$tot.hal.by,rep(NA,(YEAR-remYe))); length(Lnd.By)
  
  
  #=========================================================================================
  # Load expected bycatch generated in Code/Data_Generation/IPHC_Survey_Expected_Bycatch.R
  Foreign<-read.csv("Data_processing/Data/Harvests/Foreign_YERF_SEAK.csv",skip=1)
  N.for<-nrow(Foreign)
  forY1<-min(Foreign$Year)
  forYe<-max(Foreign$Year)
  
  For<-c(rep(0.0001,N-N.for-(YEAR-forYe)),Foreign$Estimate,rep(0.0001,(YEAR-forYe))); length(For)
  For[For == 0]<-0.0001
  forcv<-0.75
  cv.For<-c(rep(0.01,N-N.for-(YEAR-forYe)),rep(forcv,N.for),rep(0.01,(YEAR-forYe))); length(cv.For)
  
  #====================================================================================================
  # Load posterior biomass stuff
  load(file=paste("Production_models/Output/",biomassmod,"/post.Rdata", sep=""))
  bio<-MCMCchains(object=post, params="Bseo",ISB=TRUE,exact=TRUE)
  #select years we want to use - limit to years when biomass measured by ROV/SUB: 1994-2021
  
  subdys<-seq(1980,YEAR,1)
  
  bioY1<-bioY1; Y1ref<-which(subdys==bioY1)
  bioYe<-bioYe; Yeref<-which(subdys==bioYe)
  vecref1<-which(Years==bioY1)
  vecref2<-which(Years==bioYe)
  
  bioys<-seq(Y1ref,Yeref,1)
  vecref<-seq(vecref1,vecref2,1)
  
  N.bio<-length(bioys)
  
  biotab<-data.frame()
  biom<-c(rep(NA,N))
  cv.biom<-c(rep(1,N))
  
  j<-1
  for (i in bioys){ #i<-bioys[1]
    biopost<-bio[,i]
    biom[vecref[j]]<-median(biopost)
    cv.biom[vecref[j]]<-sd(biopost)/mean(biopost)
    j<-j+1
  }
  
  #*******************************************************************************
  #* SPECIAL DATA :)
  #* Beta parameters and Derby effect
  #* # can be changed in MULTIMODEL RUN for ease, but need to load base stuff here...
  #================================================================================
  # END DATA LOAD for 2022 SAFE
  #Below are other data sources for different model structures.
  Derby.Eff <- Derby.Eff
  DEsd<-DEsd
  B1<-B1
  B2<-B2
  Subd<- 4
  
  Data<-list(N=N,Years=Years,Subd=Subd,
             KnC.obs=KnC.obs,cv.KnC=cv.KnC,
             Lnd.By=Lnd.By,
             ExpByc=ExpByc,cv.ExpByc=cv.ExpByc,
             For=For,cv.For=cv.For,
             biom=biom, cv.biom=cv.biom,
             #B.obs=B.obs,cv.B=cv.B,
             #sCPUE=sCPUE,cv.sCPUE=cv.sCPUE,
             Derby.Eff = Derby.Eff, DEsd=DEsd,
             B1=B1,B2=B2)
  return(Data)
}
################################################################################
################################################################################
# Data preparation function for Stage 2 SPM model
data.prep.1888<-function(){
  
  data<-list(N=N,Subd=Subd,
             KnC.obs=KnC.obs, cv.KnC=cv.KnC,
             ExpByc = ExpByc, cv.ExpByc = cv.ExpByc,
             Lnd.By=Lnd.By, 
             For=For,cv.For=cv.For,
             biom=biom, cv.biom=cv.biom,
             #B.obs=B.obs,cv.B=cv.B,
             #sCPUE=sCPUE,cv.sCPUE=cv.sCPUE, 
             Fu=Fu, 
             OFLABC=OFLABC, #B.err=B.err,
             HarvStrat=HarvStrat,
             Derby.Eff = Derby.Eff, DEsd=DEsd,
             B1=B1,B2=B2,upvar=upvar
  )
  
  # KnCinits1<-log(KnC.obs+runif(N,-0.5*KnC.obs,0.5*KnC.obs))
  #  KnCinits2<-log(KnC.obs+runif(N,-0.5*KnC.obs,0.5*KnC.obs))
  # KnCinits3<-log(KnC.obs+runif(N,-0.5*KnC.obs,0.5*KnC.obs))
  #  KnCinits4<-log(KnC.obs+runif(N,-0.5*KnC.obs,0.5*KnC.obs))
  # KnCinits5<-log(KnC.obs+runif(N,-0.5*KnC.obs,0.5*KnC.obs))
  
  #Kinits1<-c(log(20000),log(20000),log(20000),log(20000))
  #Kinits2<-c(log(30000),log(10000),log(8000),log(25000))
  #Kinits3<-c(log(15000),log(6000),log(25000),log(5000))
  
  qsinits1<-log(runif(4,0.2,10))#c(3,4,0.5,0.3)
  qsinits2<-log(runif(4,0.2,10))#c(3,10,18,15)
  qsinits3<-log(runif(4,0.2,10))#c(10,10,12,12)
  
  #qfinits1<-log(runif(4,0.2,10))#c(3,4,5,10)
  #qfinits2<-log(runif(4,0.2,10))#c(3,10,18,15)
  #qfinits3<-log(runif(4,0.2,10))#c(10,10,12,12)
  
  inits1<-list(r=0.04, logK=log(100000), #phi=phiinits1,
               logvar=-9,
               logqs=qsinits1, #logqf=qfinits1, 
               Tau1=0.4 #Tau3=0.03, Tau4 = 0.3, 
  ) #, 
  #logKnC=KnCinits1)
  inits2<-list(r=0.07, logK=log(130000), #phi=phiinits2, 
               logvar = -8,
               logqs=qsinits2, #logqf=qfinits2,  
               Tau1=0.12)#, 
  #                logKnC=KnCinits2)
  inits3<-list(r=0.02, logK=log(80000), #phi=phiinits3, 
               logvar = -7,
               logqs=qsinits3, #logqf=qfinits3, 
               Tau1=0.2) #, 
  #logKnC=KnCinits3)
  
  #inits <- list(inits1,inits2,inits3,inits4,inits5)
  inits<-list(inits1,inits2,inits3)
  
  parameters=c(
    "r","K","phi", #"R.hyp",#"R.sig",#"phi.hyp", 
    "sigma","epsilon","PE",
    "B", "D","C","FC","KnC", "By","For.catch",#"Bsub",
    #"Bseo","Cseo","Byseo","Dseo","KnCseo",
    "Surplus", 
    #"qsCPUE", #"qfCPUE",
    "Tau1" ,#"Tau2", "Tau3", #"Tau4",
    "CtoB","CBtoK","B1980","Stock.Status",# CBtoKseo","FBtoK","FBtoKseo","FBtoCB",
    "MSY","Bmsy","Fmsy","Hmsy",
    #"Bmsyseo","Fmsyseo",
    "fit","subfit.new","iphcfit.new"
  )
  
  prep<-list(data=data,inits=inits,parameters=parameters)
  return(prep)
}
################################################################################
################################################################################
## Plotting function for Stage 2 data:
plot_spm1888<-function(post=post, Mod.title="test",filename){
  #tbl<-jags.View(post, title="", digits=3)
  
  tbl<-as.data.frame(jags.View(post, title="", digits=3))
  write.csv(tbl,
            file=paste("Production_models/Output/", filename,"/results_sum.csv", sep=""))
  
  Bseo.sum<-tbl[grepl("B",tbl$parameter),]
  #Bseo.sum$cv<-Bseo.sum$sd/Bseo.sum$mean
  #write.csv(Bseo.sum,file=paste("Model Output/",filename,"/Bseo_post_", Mod.title,".csv", sep=""))
  
  
  
  png(paste("Production_models/Figures/", filename,"/PP_check_submersibles.png",sep=""),
      width=6,height=6,#width=9.5,height=8.5,
      units="in",res=1200)
  par(mfrow=c(1,1), omi=c(0,0,0.1,0.01))
  
  
  obs<-paste("fit", sep="")
  sim<-paste("subfit.new", sep="")
  pp.check(post, observed=obs, simulated=sim, 
           xlab='log Observed biomass', ylab='log Simulated biomass',
           main="pp.check", col="grey")
  
  mtext("Posterior predictive check for 1888 model posteriors", side=3,  
        padj=1,outer=TRUE, cex=1.4)
  dev.off()
  
  #===============================================================
  #Posterior of all biomass estimates
  {
    png(paste("Production_models/Figures/", filename,"/Biomass_overview.png",sep=""),
        width=7,height=6,#width=9.5,height=8.5,
        units="in",res=1200)
    
    B.all<-MCMCchains(post,"B")
    #Bs<-Bs[,1:42]
    par(mfrow=c(1,1))
    envplot(B.all,1888:(1888+N+Fu),cols=c("black","black"),n=0,
            ylab="Biomass (t)",xlab="Year",ylim=c(0,quantile(B.all,c(0.99))))
    addpoints(post,Years=Years ,Points=biom, errordat=cv.biom, scaler=0, 
              #scaler="qfCPUE" for example...
              bar=1.95 ,col="blue", pch=17, cex=0.75)
    Y<-seq(1888,(1888+N),1)
    addenv(B.all, cols=c("black","black"),alpha.f=0.7, Y=Y)
    legend(x="topright", c("yelloweye biomass", "posterior biomass est. from stage 1 model"), 
           col=c("black","blue"),
           text.col=c("black","blue"),
           border=F, bty="n")
    dev.off()
  }
  
  #i) Plotting setup
  #Years.cont
  
  
  #1) Biomass observed versus modeled 
  
  
  #1.B) Biomass observed versus modeled with Bmsy lines
  
  
  #2) Known catch
  {KnC<-MCMCchains(post,"KnC")
    ExpBy<-MCMCchains(post,"By")
    Disc<-MCMCchains(post,"D")
    Foreign<-MCMCchains(post,"For.catch")
    Tot.Catch<-MCMCchains(post,"C")
    
    colch<-c("black","blue","black","black","red")
  }
  
  {png(paste("Production_models/Figures/", filename,"/KnownCatch_fit.png",sep=""),
       width=7,height=6,#width=9.5,height=8.5,
       units="in",res=1200)
    par(mfrow=c(3,1), mar=c(4,5,3,1))
    #i<-1
    envplot(KnC,Years,cols=c(colch[1],colch[1]),n=0,ylab="Known catches (t)",xlab="Year",
            ylim=c(0,max(KnC.obs*1.25)),
            main="Catch Data Fit",cex.axis=1.2, cex.lab=1.25, cex.main=1.5)
    addpoints(post,Years=Years ,Points=KnC.obs ,errordat=cv.KnC ,scaler=0 ,
              bar=1.95 ,col=colch[1], pch=18,cex=0.8)
    legend(x="topleft", c("Known Catches"), cex=1.5,
           col=c(colch[1],colch[1],colch[1]),
           text.col=c(colch[1],colch[1],colch[1]),border=F, bty="n", pch=c(NA,NA,NA))
    Y<-seq(1888,(1888+N-1),1)
    envplot(ExpBy,Years,cols=c(colch[2],colch[2]),n=0,ylab="Expected Bycatch (t)",xlab="Year",
            ylim=c(0,quantile(ExpBy,c(0.999))), #max(KnC.obs*1.25)),
            main="",cex.axis=1.2, cex.lab=1.25, cex.main=1.5)
    
    addpoints(post,Years=Years ,Points=ExpByc ,errordat=cv.ExpByc ,scaler=0 ,
              bar=1.95 ,col=colch[2], pch=18,cex=0.8)
    addenv(Disc-5, cols=c(colch[4],colch[4]),alpha.f=0.35, Y=Y)
    lines(Years,Lnd.By,type="l",col=colch[5], lwd=1.5, lty=1, pch=15, cex=0.5)
    legend(x="topleft", c("Expected bycatch","Discards","Landed Bycatch"), 
           col=c(colch[2],colch[4],colch[5]),cex=1.25,
           text.col=c(colch[2],colch[4],colch[5]),border=F, bty="n", pch=c(NA,NA,NA))
    
    envplot(Foreign,Years,cols=c(colch[3],colch[3]),n=0,ylab="Foreign removals (t)",xlab="Year",
            ylim=c(0,quantile(Foreign,c(0.999))), #max(KnC.obs*1.25)),
            main="",cex.axis=1.2, cex.lab=1.25, cex.main=1.5)
    #addenv(Foreign, cols=c(colch[3],colch[3]),alpha.f=0.7, Y=Y)
    addpoints(post,Years=Years ,Points=For ,errordat=cv.For ,scaler=0 ,
              bar=1.95 ,col=colch[3], pch=18,cex=0.8)
    legend(x="topleft", c("Foreign fleet harvest"), 
           col=c(colch[3],colch[3],colch[3]),cex=1.5,
           text.col=c(colch[3],colch[3],colch[3]),border=F, bty="n", pch=c(NA,NA,NA))
    #abline(v=2022,  col=colch[3], lty=2)
    #if (i == 2){
    #  legend(x="topright", c("posterior est.","observed"), col=c(colch[1],colch[2]),
    #         text.col=c(colch[1],colch[2]),border=F, bty="n", pch=c(NA,18))
    #} else {}
    
    dev.off()}
  
  #3) Expected, landed and modeled bycatch... 
  
  #4) Discard plots
  
  #5) Catch plots of total, known, bycatch, landed bycatch, and discards
  
  {
    png(paste("Production_models/Figures/", filename,"/Catch_breakdown.png",sep=""),
        width=6,height=8,#width=9.5,height=8.5,
        units="in",res=1200)
    par(mfrow=c(2,1))
    envplot(Tot.Catch,Years,cols=c("black","black"),n=0,ylab="Catches (t)",xlab="Year",
            ylim=c(0,quantile(Tot.Catch,c(0.999))),
            cex.axis=1, cex.lab=1.25, border=TRUE)
    legend(x="topleft", c("Est. Total removals"),
           col=c("black","blue","darkorange","purple","darkgrey"),
           text.col=c("black","blue","darkorange","purple","darkgrey"), 
           border=F, bty="n", pch=c(NA),cex=1.2)
    envplot(Foreign,Years,cols=c("red","red"),n=0,ylab="Catches (t)",xlab="Year",
            ylim=c(0,quantile(Tot.Catch,c(0.99))),
            cex.axis=1, cex.lab=1.25)
    #Y<-seq(1980,(1980+N.subd-1),1)
    addenv(KnC,cols=c("blue","blue"), Y=Y,alpha.f=0.6)
    addenv(Disc-10,cols=c("black","black"), Y=Y,alpha.f=0.6)
    addenv(Foreign,cols=c("red","red"), Y=Y,alpha.f=0.1)
    #addenv(Tot.Catch,cols=c("black","black"),Y=Y, alpha.f=0.01)
    #lines(colMedians(Tot.Catch)~Years, col="black", lwd=2)
    
    legend(x="topleft", c("Know catches","Est. Discards","Est. Foreign Fleet catches"),
           col=c("red","blue","black","purple","darkgrey"),
           text.col=c("blue","black","red","purple","darkgrey"), 
           border=F, bty="n", pch=c(NA),cex=1.2)
    dev.off()
    
  }
  
  #6) Surplus production estimates
  {Surplus<-MCMCchains(post,"Surplus")
    
    colch<-c("forestgreen","black","purple")
    
    png(paste("Production_models/Figures/", filename,"/Surplus.png",sep=""),
        width=7,height=6,#width=9.5,height=8.5,
        units="in",res=1200)
    par(mfrow=c(1,1), mar=c(4,5,3,1))
    #i<-1
    envplot(Surplus,Years,cols=c(colch[1],colch[1]),n=0,ylab="Surplus Production (t)",xlab="Year",
            ylim=c(0,max(quantile(Tot.Catch,c(0.999)),quantile(Tot.Catch,c(0.99)))),
            main="",cex.axis=0.8, cex.lab=1.25, cex.main=1.5)
    addenv(Tot.Catch,cols=c("black","black"), Y=Y)
    
    legend(x="topleft", c("Surplus","Removals"), col=c(colch[1],colch[2],colch[3]),
           text.col=c(colch[1],colch[2],colch[3]),border=F, bty="n", pch=c(NA),cex=1.2)
    
    dev.off()
  }
  
  colch<-c("black","forestgreen")
  png(paste("Production_models/Figures/", filename,"/Surplus2.png",sep=""),
      width=7,height=6,#width=9.5,height=8.5,
      units="in",res=1200)
  par(mfrow=c(1,1), mar=c(4,5,3,1))
  #i<-1
  envplot(Tot.Catch,Years,cols=c(colch[1],colch[1]),n=0,
          ylab="Total Catch & Surplus Production (t)",xlab="Year",
          ylim=c(0,max(quantile(Tot.Catch,c(0.999)),quantile(Tot.Catch,c(0.99)))),
          main="",cex.axis=0.8, cex.lab=1.25, cex.main=1.5)
  addenv(Surplus,cols=c(colch[2],colch[2]), Y=Y)
  addenv(Tot.Catch,cols=c(colch[1],colch[1]), Y=Y)
  
  legend(x="topleft", c("Total Removals","Surplus"), col=c(colch[1],colch[2],colch[3]),
         text.col=c(colch[1],colch[2],colch[3]),border=F, bty="n", pch=c(NA),cex=1.2)
  
  dev.off()
  
  
  
  #7) CtoB
  {CtoB<-MCMCchains(post,"CtoB")
    
    Fmsy<-tbl$mean[tbl$parameter=="Fmsy"]#posteriors$Fmsy[1,4]#par.ext("Fmsy")
    
    Flo<-tbl$`2.5%`[tbl$parameter=="Fmsy"]; Fhi<-tbl$`97.5%`[tbl$parameter=="Fmsy"]
    
    colch<-c("black","yellow2","orange")
    #vars<-list(CtoB.ey,CtoB.ns,CtoB.cs,CtoB.ss); vars2<-list(Fmsy.ey,Fmsy.ns,Fmsy.cs,Fmsy.ss)
    #varlo<-list(Flo.ey,Flo.ns,Flo.cs,Flo.ss); varup<-list(Fhi.ey,Fhi.ns,Fhi.cs,Fhi.ss)}
    
    png(paste("Production_models/Figures/", filename,"/CtoB.png",sep=""),
        width=7,height=6,#width=9.5,height=8.5,
        units="in",res=1200)
    par(mfrow=c(1,1), mar=c(4,5,3,1))
    envplot(CtoB,Years,cols=c(colch[1],colch[1]),n=0,ylab="Catch:Biomass",xlab="Year",
            ylim=c(0,max(quantile(CtoB,c(0.999)))),
            main="",cex.axis=0.8, cex.lab=1.25, cex.main=1.5)
    abline(h=Fmsy,col=colch[3], lty=1, lwd=1.5)
    abline(h=Flo,col=colch[3], lty=2, lwd=1)
    abline(h=Fhi,col=colch[3], lty=2, lwd=1)
    polygon(x=c(Years,rev(Years)),
            y=c(rep(Flo,length(Years)),rep(Fhi,length(Years))),
            border=NA,
            col=adjustcolor(colch[2], alpha.f=.10))
    # addenv(vars[[i]], cols=c(colch[1],colch[1]))
    legend(x="topright", c("Catch:Biomass",paste0("Fmsy = ",Fmsy)), col=c(colch[1],colch[3],colch[3]),
           text.col=c(colch[1],colch[3],colch[3]),border=F, bty="n", pch=c(NA,NA,15),cex=1)
    
    dev.off()
    
    #8) Epsilons/ Process Error
    png(paste("Production_models/Figures/", filename,"/epsilon.png",sep=""),
        width=7,height=6,#width=9.5,height=8.5,
        units="in",res=1200)
    pe<-MCMCchains(post,"epsilon")
    pe2<-MCMCchains(post,"PE")
    par(mfrow=c(1,1))
    envplot(pe,1888:(1888+N+Fu),cols=c("black","black"),n=0,ylab="epsilon",xlab="Year",ylim=c(-0.1,0.1))
    abline(h=0,col="red",lty=2)
    #envplot(pe2,1888:(1888+N+Fu),cols=c("black","black"),n=0,ylab="process error",xlab="Year",ylim=c(-0.1,0.1))
    #abline(h=0,col="red",lty=2)
    dev.off()
    
    #==============================================================================
    # RISK ANALYSIS: the probability that Fu years in the future population is at or
    # above Bmsy; 
    
    {
      Bmsy <-MCMCchains(object=post,  params = c("Bmsy"),  excl = NULL,
                        ISB = FALSE,  #ignore square brackets
                        mcmc.list = FALSE,chain_num = NULL, exact=TRUE
      )#; nrow(Bmsy.ey)
      
      BCurrent<-MCMCchains(object=post,  params = c(paste("B[",N,"]",sep="")),  excl = NULL,
                           ISB = FALSE,  #ignore square brackets
                           mcmc.list = FALSE,chain_num = NULL, exact=TRUE
      )#; nrow(BCurrent)
      BFuture<-MCMCchains(object=post,  params = c(paste("B[",N+Fu,"]",sep="")),  excl = NULL,
                          ISB = FALSE,  #ignore square brackets
                          mcmc.list = FALSE,chain_num = NULL, exact=TRUE
      )#; nrow(BFuture)
      Risk<-cbind(Bmsy,BCurrent,BFuture)
      
      colnames(Risk)<-c("Bmsy","Bcur","Bfut")
      Risk<-as.data.frame(Risk)
      Risk$current.stat<-Risk$Bcur - Risk$Bmsy
      Risk$future.stat<-Risk$Bfut - Risk$Bmsy
      Risk$increaseQ<-Risk$Bfut - Risk$Bcur
      
      CurrentRisk<-nrow(Risk[Risk$current.stat > 1,])/nrow(Risk)
      FutureRisk<-nrow(Risk[Risk$future.stat > 1,])/nrow(Risk)
      IncreaseProb<-nrow(Risk[Risk$increaseQ > 1,])/nrow(Risk)
      
      Risk.sum<-data.frame()
      i<-1
      
      Bmsyi <- Bmsy #MCMCchains(object=post,  params = c(paste("Bmsy[",i,"]",sep="")),  excl = NULL,
      #    ISB = FALSE,  #ignore square brackets
      #   mcmc.list = FALSE,chain_num = NULL, exact=TRUE
      
      BCur<-BCurrent #MCMCchains(object=post,  params = c(paste("B[",i,",",N.subd,"]",sep="")),  excl = NULL,
      #ISB = FALSE,  #ignore square brackets
      #mcmc.list = FALSE,chain_num = NULL, exact=TRUE
      
      BFut<- BFuture #MCMCchains(object=post,  params = c(paste("B[",i,",",N.subd+Fu,"]",sep="")),  excl = NULL,
      # ISB = FALSE,  #ignore square brackets
      #mcmc.list = FALSE,chain_num = NULL, exact=TRUE
      
      Risk<-cbind(Bmsy,BCur,BFut) #head(Risk)
      colnames(Risk)<-c("Bmsy","Bcur","Bfut")
      Risk<-as.data.frame(Risk)
      Risk$current.stat<-Risk$Bcur - Risk$Bmsy
      Risk$future.stat<-Risk$Bfut - Risk$Bmsy
      Risk$increaseQ<-Risk$Bfut - Risk$Bcur
      
      CR<-nrow(Risk[Risk$current.stat > 1,])/nrow(Risk)
      FR<-nrow(Risk[Risk$future.stat > 1,])/nrow(Risk)
      IP<-nrow(Risk[Risk$increaseQ > 1,])/nrow(Risk)
      
      Risk.sum[i,"Subd"]<-"SEO"
      Risk.sum[i,"Above_Bmsy_now"]<-CR
      Risk.sum[i,"Above_Bmsy_future"]<-FR
      Risk.sum[i,"Pop_inc_future"]<-IP
      
      Risk.sum[5,"Subd"]<-"SEO"
      Risk.sum[5,"Above_Bmsy_now"]<-CurrentRisk
      Risk.sum[5,"Above_Bmsy_future"]<-FutureRisk
      Risk.sum[5,"Pop_inc_future"]<-IncreaseProb
    }
    write.csv(Risk.sum,file=paste("Production_models/Output/", filename,"/Risk_Summary_", Mod.title,".csv", sep=""))
    
  }
  
  #-----------------------------------------------------------------------------
  #Extract B1980 for phi prior in 1980 model:
  B1980<-MCMCchains(  object=post,  params = "B1980",  excl = NULL,
                      ISB = TRUE,  mcmc.list = FALSE,chain_num = NULL)
  logB1980<-log(B1980)
  
  betas1<-estBetaParams(mean(B1980),var(B1980))
  
  Q95<-quantile(B1980,c(0.95))
  
  gammafind <- function(shape) {
    scale <- mean(B1980)/shape
    pgamma(Q95, shape, scale=scale) - 0.95
    
  }
  
  tmp <- uniroot( gammafind, lower=0.1, upper=100 )
  gam.shape<-tmp$root
  gam.scale<-1/(mean(B1980)/tmp$root)
  
  png(paste("Production_models/Figures/", filename,"/B1980_prior.png",sep=""),
      width=7,height=6,#width=9.5,height=8.5,
      units="in",res=1200)
  plot(density(B1980), main="B1980 (phi)")
  for (i in 1:100){ #i<-1
    betad<-rbeta(length(B1980),betas1$alpha,betas1$beta)
    lines(density(betad),col="purple")
    norm<-rnorm(length(B1980),mean(B1980),sd(B1980))
    lines(density(norm),col="orange")
    lnorm<-rlnorm(length(B1980),log(mean(B1980)),sd(B1980))
    lines(density(lnorm),col="red")
    gams<-rgamma(length(B1980),gam.shape,gam.scale)
    lines(density(gams),col="green")
  }
  lines(density(B1980),col="black",lwd=2)
  legend(x="topright", c(paste0("beta: ",round(betas1$alpha,3),", ",round(betas1$beta,3),sep=""),
                         paste0("norm: ",round(mean(B1980),3),", ",round(sd(B1980),3),sep=""),
                         paste0("lnorm: ",round(log(mean(B1980)),3),", ",round(sd(B1980),3),sep=""),
                         paste0("gamma: ",round(gam.shape,3),", ",round(gam.scale,3),sep="")), 
         col=c("purple","orange","red","green"),
         text.col=c("purple","orange","red","forestgreen"),border=F, bty="n", pch=c(NA),cex=0.75)
  dev.off()
  
  #logphi---------------------------------------------------------------------------
  betas2<-estBetaParams(mean(logB1980),var(logB1980))
  
  png(paste("Production_models/Figures/", filename,"/logB1980_prior.png",sep=""),
      width=7,height=6,#width=9.5,height=8.5,
      units="in",res=1200)
  plot(density(logB1980), main="log phi")
  for (i in 1:100){ #i<-1
    norm<-rnorm(length(logB1980),mean(logB1980),sd(logB1980))
    lines(density(norm),col="orange")
  }
  lines(density(logB1980),col="black",lwd=2)
  legend(x="topright", c(paste0("norm; ",round(mean(logB1980),3),", ",round(sd(logB1980),3))), col=c("orange"),
         text.col=c("orange","red","green"),border=F, bty="n", pch=c(NA),cex=1)
  dev.off()
  
  ##extract K posterior as well!!---------------------------------------------------
  bigK<-MCMCchains(  object=post,  params = "K",  excl = NULL,
                     ISB = TRUE,  mcmc.list = FALSE,chain_num = NULL)
  logK<-log(bigK)
  
  betas3<-estBetaParams(mean(bigK),var(bigK))
  
  Q95<-quantile(bigK,c(0.95))
  gammafind <- function(shape) {
    scale <- mean(bigK)/shape
    pgamma(Q95, shape, scale=scale) - 0.95
    
  }
  tmp <- uniroot( gammafind, lower=0.1, upper=100 )
  gam.shape<-tmp$root
  gam.scale<-1/(mean(bigK)/tmp$root)
  
  png(paste("Production_models/Figures/", filename,"/bigK_prior.png",sep=""),
      width=7,height=6,#width=9.5,height=8.5,
      units="in",res=1200)
  plot(density(bigK), main="K")
  for (i in 1:100){ #i<-1
    norm<-rnorm(length(bigK),mean(bigK),sd(bigK))
    lines(density(norm),col="orange")
    gams<-rgamma(length(bigK),gam.shape,gam.scale)
    lines(density(gams),col="green")
  }
  lines(density(bigK),col="black",lwd=2)
  legend(x="topright", c(paste0("norm; ",round(mean(bigK),0),", ",round(sd(bigK),0)),
                         paste0("gamma; ",round(gam.shape,3),", ",round(gam.scale,7))),
         col=c("orange","green"),
         text.col=c("orange","forestgreen"),border=F, bty="n", pch=c(NA),cex=1)
  dev.off()
  
  #logK-----------------------------------------------------------------------------
  betas4<-estBetaParams(mean(logK),var(logK))
  
  Q95<-quantile(logK,c(0.95))
  #gammafind <- function(shape) {
  #  scale <- mean(logK)/shape
  #  pgamma(Q95, shape, scale=scale) - 0.95
  #  
  #}
  #tmp <- uniroot( gammafind, lower=0.1, upper=1000 )
  #gam.shape<-tmp$root
  #gam.scale<-1/(mean(logK)/tmp$root)
  #logK<-log(bigK)
  #logK<-logK[logK<quantile(logK,c(0.90))]
  png(paste("Production_models/Figures/", filename,"/logbigK_prior3.png",sep=""),
      width=7,height=6,#width=9.5,height=8.5,
      units="in",res=1200)
  plot(density(logK),main="log K")
  for (i in 1:100){ #i<-1
    norm<-rnorm(length(logK),mean(logK),sd(logK))
    lines(density(norm),col="orange")
    norm2<-rnorm(length(logK),median(logK),sd(logK))
    lines(density(norm2),col="blue")
  }
  lines(density(logK),col="black",lwd=2)
  legend(x="topright",c(paste0("norm-mean,",round(mean(logK),3),", ",round(sd(logK),3)),
                        paste0("norm-med, ",round(median(logK),3),", ",round(sd(logK),3))), 
         col=c("orange","blue"),
         text.col=c("orange","blue"),border=F, bty="n", pch=c(NA),cex=1)
  dev.off()
}
################################################################################
################################################################################
## Prepare data for Stage 3 model:
data.prep.PHASE3<-function(){
  
  data<-list(N=N,Subd=Subd,
             KnC.obs=KnC.obs, cv.KnC=cv.KnC,
             ExpByc = ExpByc, cv.ExpByc = cv.ExpByc,
             Lnd.By=Lnd.By, 
             B.obs=B.obs,cv.B=cv.B, #CKP=CKP,
             sCPUE=sCPUE, cv.sCPUE=cv.sCPUE, #,
             #fCPUE=sCPUE, cv.fCPUE=cv.fCPUE, 
             Fu=Fu, 
             OFLABC=OFLABC, #B.err=B.err,
             HarvStrat=HarvStrat,
             Derby.Eff = Derby.Eff, DEsd=DEsd,
             B1=B1,B2=B2,
             phiB1=phiB1,phiB2=phiB2,
             bigKmu=bigKmu,bigKsigma=bigKsigma,
             phimu=phimu,phisig=phisig,upvar=upvar
  )
  
  KnCinits1<-log(KnC.obs+runif(N*Subd,-0.5*KnC.obs,0.5*KnC.obs))
  KnCinits2<-log(KnC.obs+runif(N*Subd,-0.5*KnC.obs,0.5*KnC.obs))
  KnCinits3<-log(KnC.obs+runif(N*Subd,-0.5*KnC.obs,0.5*KnC.obs))
  KnCinits4<-log(KnC.obs+runif(N*Subd,-0.5*KnC.obs,0.5*KnC.obs))
  KnCinits5<-log(KnC.obs+runif(N*Subd,-0.5*KnC.obs,0.5*KnC.obs))
  
  #Kinits1<-c(log(20000),log(20000),log(20000),log(20000))
  #Kinits2<-c(log(30000),log(10000),log(8000),log(25000))
  #Kinits3<-c(log(15000),log(6000),log(25000),log(5000))
  
  phiinits1<-log(c(0.8,0.7,0.6,0.8))
  phiinits2<-log(c(0.6,0.5,0.7,0.6))
  phiinits3<-log(c(0.5,0.8,0.6,0.7))
  
  qsinits1<-log(runif(4,0.2,10))#c(3,4,0.5,0.3)
  qsinits2<-log(runif(4,0.2,10))#c(3,10,18,15)
  qsinits3<-log(runif(4,0.2,10))#c(10,10,12,12)
  
  #qfinits1<-log(runif(4,0.2,10))#c(3,4,5,10)
  #qfinits2<-log(runif(4,0.2,10))#c(3,10,18,15)
  #qfinits3<-log(runif(4,0.2,10))#c(10,10,12,12)
  
  inits1<-list(R.hyp=0.04, #logK=Kinits1, #phi=phiinits1,
               logphi=phiinits1, logvar=-9,
               logqs=qsinits1, #logqf=qfinits1, 
               Tau1=0.4, #Tau3=0.03, Tau4 = 0.3, 
               Tau3=0.01, 
               logKnC=KnCinits1)
  inits2<-list(R.hyp=0.07, #logK=Kinits2, #phi=phiinits2, 
               logphi=phiinits2, logvar = -8,
               logqs=qsinits2, #logqf=qfinits2,  
               Tau1=0.12, #Tau3=0.1, Tau4 = 0.2,
               Tau3=0.3, 
               logKnC=KnCinits2)
  inits3<-list(R.hyp=0.02, #logK=Kinits3, #phi=phiinits3, 
               logphi=phiinits3, logvar = -7,
               logqs=qsinits3, #logqf=qfinits3, 
               Tau1=0.2, #Tau3=0.13, Tau4=0.02,
               Tau3=0.5, 
               logKnC=KnCinits3)
  
  #inits <- list(inits1,inits2,inits3,inits4,inits5)
  inits<-list(inits1,inits2,inits3)
  
  parameters=c(
    "r","K","phi","bigphi", "R.hyp",#"R.sig",#"phi.hyp", 
    "Kseo", 
    "sigma","eta","rB1","rB2","pi",
    "phi.eta","littlephiB1","littlephiB2",
    "phiTau",
    "epsilon","PE",
    "B", "D","C","FC","KnC", "By",
    "Bseo","Cseo","Byseo","Dseo","KnCseo",
    "S", "Surplus", "S.seo","Surplus.seo",
    "qsCPUE", #"qfCPUE",
    "Tau1" ,"Tau2", "Tau3", #"Tau4",
    "CtoB","CBtoK","CBtoKseo","FBtoK","FBtoKseo","FBtoCB",
    "MSY","Bmsy","Fmsy","Hmsy","Stock.Status","Proj1_biomass",
    "Bmsyseo","Fmsyseo","Fmsyseo2","Stock.Status.SEO",
    "fit","subfit.new","iphcfit.new"
  )
  
  prep<-list(data=data,inits=inits,parameters=parameters)
  return(prep)
}

################################################################################
################################################################################
## Posterior plotting functions: 

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

#Add second envelope plot to existing plot
addenv<-function(df, cols=c("blue","blue"),alpha.f=0.4, Y=Y){ #df<-B.all
  means<-colMeans(df)
  qmat <- apply(df, 2, quantile, p=c(0.05,0.25,0.5,0.75,0.95))
  x1<-Y #length(Y)
  polygon(c(x1,rev(x1)),c(qmat[1,1:(ncol(df)-0)],qmat[5,rev(1:(ncol(df)-0))]),border=NA,col=adjustcolor(cols[1], alpha.f=.15))
  polygon(c(x1,rev(x1)),c(qmat[2,1:(ncol(df)-0)],qmat[4,rev(1:(ncol(df)-0))]),border=NA,col=adjustcolor(cols[1], alpha.f=.2))
  
  lines(x1,qmat[3,1:(ncol(df)-0)],col=adjustcolor(cols[1],alpha.f=alpha.f),lwd=3)
}

#Add data points with error bars to existing envelope plot
#scaler is the posterior q for index... will extract midpoint; eg: scaler="qsCPUE"
# scaler=0 if same scale
addpoints<-function(Post=post ,Years=Years.cont ,Points=B.obs[1,] ,errordat=cv.B[1,] ,scaler=0 ,
                    bar=1.68 ,col="red" ,pch, cex, ...){  #pch=18  #
  cv<-errordat
  if (scaler == 0){
    qsc<-1
    arrows(x0=Years, y0=Points-bar*cv*Points, x1=Years, y1=Points+bar*cv*Points, code=3, 
           col=col, lwd=1, angle=90, length=0.01)
  } else {
    parpl<-paste(scaler,"[",i,"]",sep="")
    qsc<-midpoint(post=Post,par=parpl)
    arrows(x0=Years, y0=Points/qsc-bar*cv*Points/qsc, x1=Years, y1=Points/qsc+bar*cv*Points/qsc, code=3, 
           col=col, lwd=1, angle=90, length=0.01)
  }
  points(Years,Points/qsc,col=col,pch=pch, cex=cex)
  #cv<-errordat
  #arrows(x0=Years, y0=Points-bar*cv*Points, x1=Years, y1=Points+bar*cv*Points, code=3, 
  #       col=col, lwd=1, angle=90, length=0.01)
}

##################################################################################'
##################################################################################
# Function Development: 


# get MCMC posterior values for parameter, year and area/subdistrict
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

#extract the midpoint of a posterior
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

#=================================================================================
# Parameter plotting a
paramcomp<-function(SEOq="FALSE",Mods,Modnames,param, xlow,xhi,xbreak,pscale=1,xlabel){
  
  Mod.num<-length(Mods)
  for (m in 1:Mod.num){   #m<-2
    sdlist<-c("EYKT", "NSEO", "CSEO", "SSEO")
    M<-Mods[[m]]
    
    if (SEOq=="TRUE") { ##!! for SEO level 
      P<-MCMCchains(object=M, params=c(param),ISB=FALSE,exact=TRUE)
      posterior<-data.frame(matrix(nrow=length(P)))
      posterior$Mod<-Modnames[m]
      #Param$subd<-sdlist[1]
      posterior$posterior<-P[1:length(P)]
      comp.mod<-posterior
      
    } else { ##!! for subdistrict parameters
      
      for (i in 1:4){  #i<-1
        cl<-paste(param,"[",i,"]",sep="")
        P<-MCMCchains(object=M, params=c(cl),ISB=FALSE,exact=TRUE)
        posterior<-data.frame(matrix(nrow=length(P)))
        posterior$Mod<-Modnames[m]
        posterior$subd<-sdlist[i]
        posterior$posterior<-P[1:length(P)]
        if (i == 1) {
          comp.mod<-posterior
        } else {
          comp.mod<-rbind(comp.mod,posterior)
        }
      }
    }
    
    if (m == 1) {
      comp<-comp.mod
    } else {
      comp<-rbind(comp,comp.mod)
    }
  }
  
  #plot
  if (SEOq=="TRUE") {
    comp %>% 
      ggplot(aes(posterior, Mod, group = Mod, fill = 0.5 - abs(0.5-stat(ecdf)))) +#Mod)) + 
      #   geom_density_ridges(aes(point_fill = Mod, point_color = Mod),
      #                        alpha = 0.3, scale=pscale, jittered_points = FALSE,
      #                        rel_min_height = 0.001) +
      #  stat_density_ridges(aes(point_fill = Mod, point_color = Mod, color=Mod),
      stat_density_ridges(aes(point_fill = Mod, point_color = Mod),
                          geom = "density_ridges_gradient", calc_ecdf =  TRUE,
                          quantile_lines = FALSE, # quantiles = c(0.025,0.5,0.975), 
                          alpha=0.3,rel_min_height = 0.001, scale=pscale) +
      scale_fill_viridis_c(direction=1, option="cividis", #alpha=0.5,
                           begin=0.5,end=0.9) + #theme_ridges()
      xlim(xlow, xhi)+
      #geom_vline(aes(xintercept=mean(posterior), group = Mod, color=Mod),size=1) +
      xlab(xlabel) + 
      ylab("Model") +
      theme(legend.position = "none") + 
      theme(text=element_text(size=16),
            axis.text.x = element_text(size=14,angle=45, hjust=1),
            axis.text.y = element_text(size=14)) +
      scale_x_continuous(label=comma)#+
    #scale_x_continuous(limits=c(xlow,xhi),breaks = c(seq(from=xlow, to=xhi, by=xbreak))) 
  } else {
    comp %>% 
      ggplot(aes(posterior, Mod, group = Mod, fill = 0.5 - abs(0.5-stat(ecdf)))) +#Mod)) + 
      stat_density_ridges(aes(point_fill = Mod, point_color = Mod),
                          geom = "density_ridges_gradient", calc_ecdf =  TRUE,
                          quantile_lines = FALSE, # quantiles = c(0.025,0.5,0.975), 
                          alpha=0.3,rel_min_height = 0.001, scale=pscale) +
      scale_fill_viridis_c(direction=1, option="cividis",#option="inferno", #alpha=0.5,
                           begin=0.5,end=0.9) +
      #ggplot(aes(posterior, Mod, group = Mod, fill = Mod)) + 
      #geom_density_ridges(aes(point_fill = Mod, point_color = Mod),
      #                    alpha = 0.3, scale=pscale, jittered_points = FALSE,
      #                    rel_min_height = 0.001) +
      #geom_vline(xintercept=mean(posterior),linetype="solid",color=Mod,size=1) +
      facet_wrap(~ subd, ncol=1) +
      xlim(xlow, xhi)+
      xlab(xlabel) + 
      ylab("Model") +
      theme(legend.position = "none") + 
      theme(text=element_text(size=16),
            axis.text.x = element_text(size=11,angle=45, hjust=1),
            axis.text.y = element_text(size=11)) +
      scale_x_continuous(label=comma)#limits=c(xlow,xhi),breaks = c(seq(from=xlow, to=xhi, by=xbreak))) 
  }
}

#================================================================================
paramcomp2<-function(SEOq="FALSE",Mods,Modnames,param, xlow,xhi,xbreak,pscale=1,xlabel,vvalue){
  
  Mod.num<-length(Mods)
  for (m in 1:Mod.num){   #m<-2
    sdlist<-c("EYKT", "NSEO", "CSEO", "SSEO")
    M<-Mods[[m]]
    
    if (SEOq=="TRUE") { ##!! for SEO level 
      P<-MCMCchains(object=M, params=c(param),ISB=FALSE,exact=TRUE)
      posterior<-data.frame(matrix(nrow=length(P)))
      posterior$Mod<-Modnames[m]
      #Param$subd<-sdlist[1]
      posterior$posterior<-P[1:length(P)]
      comp.mod<-posterior
      
    } else { ##!! for subdistrict parameters
      
      for (i in 1:4){  #i<-1
        cl<-paste(param,"[",i,"]",sep="")
        P<-MCMCchains(object=M, params=c(cl),ISB=FALSE,exact=TRUE)
        posterior<-data.frame(matrix(nrow=length(P)))
        posterior$Mod<-Modnames[m]
        posterior$subd<-sdlist[i]
        posterior$posterior<-P[1:length(P)]
        posterior$vvalue2<-vvalue[i]
        if (i == 1) {
          comp.mod<-posterior
        } else {
          comp.mod<-rbind(comp.mod,posterior)
        }
      }
    }
    
    if (m == 1) {
      comp<-comp.mod
    } else {
      comp<-rbind(comp,comp.mod)
    }
  }
  
  #plot
  if (SEOq=="TRUE") {
    comp %>% 
      ggplot(aes(posterior, Mod, group = Mod, fill = 0.5 - abs(0.5-stat(ecdf)))) +#Mod)) + 
      #   geom_density_ridges(aes(point_fill = Mod, point_color = Mod),
      #                        alpha = 0.3, scale=pscale, jittered_points = FALSE,
      #                        rel_min_height = 0.001) +
      #  stat_density_ridges(aes(point_fill = Mod, point_color = Mod, color=Mod),
      stat_density_ridges(aes(point_fill = Mod, point_color = Mod),
                          geom = "density_ridges_gradient", calc_ecdf =  TRUE,
                          quantile_lines = FALSE, # quantiles = c(0.025,0.5,0.975), 
                          alpha=0.3,rel_min_height = 0.001, scale=pscale) +
      scale_fill_viridis_c(direction=1, option="cividis", #alpha=0.5,
                           begin=0.5,end=0.9) + #theme_ridges()
      geom_vline(xintercept=vvalue, color="red")+
      xlim(xlow, xhi)+
      #geom_vline(aes(xintercept=mean(posterior), group = Mod, color=Mod),size=1) +
      xlab(xlabel) + 
      ylab("Model") +
      theme(legend.position = "none") + 
      theme(text=element_text(size=16),
            axis.text.x = element_text(size=14,angle=45, hjust=1),
            axis.text.y = element_text(size=14)) +
      scale_x_continuous(label=comma)#+
    #scale_x_continuous(limits=c(xlow,xhi),breaks = c(seq(from=xlow, to=xhi, by=xbreak))) 
  } else {
    comp %>% 
      ggplot(aes(posterior, Mod, group = Mod, fill = 0.5 - abs(0.5-stat(ecdf)))) +#Mod)) + 
      stat_density_ridges(aes(point_fill = Mod, point_color = Mod),
                          geom = "density_ridges_gradient", calc_ecdf =  TRUE,
                          quantile_lines = FALSE, # quantiles = c(0.025,0.5,0.975), 
                          alpha=0.3,rel_min_height = 0.001, scale=pscale) +
      scale_fill_viridis_c(direction=1, option="cividis",#option="inferno", #alpha=0.5,
                           begin=0.5,end=0.9) +
      #ggplot(aes(posterior, Mod, group = Mod, fill = Mod)) + 
      #geom_density_ridges(aes(point_fill = Mod, point_color = Mod),
      #                    alpha = 0.3, scale=pscale, jittered_points = FALSE,
      #                    rel_min_height = 0.001) +
      #geom_vline(xintercept=mean(posterior),linetype="solid",color=Mod,size=1) +
      facet_wrap(~ subd, ncol=1) +
      geom_vline(mapping=aes(xintercept=vvalue2), color="red")+
      xlim(xlow, xhi)+
      xlab(xlabel) + 
      ylab("Model") +
      theme(legend.position = "none") + 
      theme(text=element_text(size=16),
            axis.text.x = element_text(size=11,angle=45, hjust=1),
            axis.text.y = element_text(size=11)) +
      scale_x_continuous(label=comma)#limits=c(xlow,xhi),breaks = c(seq(from=xlow, to=xhi, by=xbreak))) 
  }
}

#=================================================================================
# old pull post stuff
pullpost <- function(x, varname) {
  x[,substr(names(x), 1, nchar(varname))==varname] 
}

get.post=function(post.samp, var, do.post=F, do.plot=F, windows=T){
  require(coda)
  parmfrow <- par("mfrow")
  
  #coerce to matrix if mcmc.list
  if (!is.mcmc.list(post.samp) & !is.matrix(post.samp)) stop("post.samp is not of class mcmc.list or matrix")
  if (is.mcmc.list(post.samp)) post.samp=as.matrix(post.samp)
  
  #pull out posteriors for requested variable
  if (substr(var,nchar(var), nchar(var))=="[") post=post.samp[,substr(colnames(post.samp), 1, nchar(var))==var]
  else post=post.samp[,var]
  
  #if it has subscripts, apply
  if (is.matrix(post)) {
    post.est=apply(post, 2, function(x) c(mean=mean(x), sd=sd(x), quantile(x, c(0.5, 0.025, 0.975))))
    if (do.plot==T ){
      n.iter=dim(post)[1]/2
      if(windows) windows(record=T)
      par(mfrow=c(4,2), mar=c(2, 2, 1.5, 1.5), oma=c(1,1,1,1))
      
      chain1=post[1:n.iter,]
      chain2=post[(n.iter+1):dim(post)[1],]
      
      for (i in 1:dim(post)[2]){
        plot(density(chain1[,i]), type="l", col="Blue", xlim=c(min(chain1[,i], chain2[,i]), max(chain1[,i], chain2[,i])),
             xlab=" ", ylab="Density", main=paste("Posterior of ", colnames(post)[i], sep=""))
        
        lines(density(chain2[,i]), col="Red")
        
        plot(chain1[,i], type="l", col="Blue", ylim=c(min(chain1[,i], chain2[,i]), max(chain1[,i], chain2[,i])),
             xlab="Iteration", ylab=" ", main=paste("Trace of ", colnames(post)[i], sep=""))
        lines(chain2[,i], col="Red")
      }
    }
  }
  
  if (is.vector(post)) {
    post.est=c(mean=mean(post), sd=sd(post), quantile(post, c(0.5, 0.025, 0.975)))
    
    if (do.plot==T){
      n.iter=length(post)/2
      chain1=post[1:n.iter]
      chain2=post[(n.iter+1):length(post)]
      
      if(windows) windows(width=12, height=8)
      par(mfrow=c(1,2))
      
      plot(density(chain1), col="Blue", xlab="", xlim=c(min(chain1, chain2), max(chain1, chain2)), main=paste("Density of ", var, sep=""))
      lines(density(chain2), col="Red")
      
      plot(chain1, type="l", col="Blue",xlab="Iteration", ylab="", ylim=c(min(chain1, chain2), max(chain1, chain2)), main=paste("Trace of ", var, sep=""))
      lines(chain2, col="Red")
    }
  }
  #par(mfrow=parmfrow)
  
  if (do.post==T) list(posterior=post, summary=post.est)
  else post.est
}

#########################################################################
# Crosscorrelation plotting exam functions: 
# from: https://www.r-bloggers.com/2011/12/mcmc-chain-analysis-and-convergence-diagnostics-with-coda-in-r/#:~:text=To%20check%20for%20pairwise%20correlations%20is%20quite%20easy,pair%20correlation%20plot%20that%20you%20can%20see%20below.

library(IDPmisc)

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="blue4", ...)
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, method = "spearman"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex * r)
}

betterPairs <- function(YourData){
  return(pairs(YourData, lower.panel=function(...) {par(new=TRUE);ipanel.smooth(...)}, diag.panel=panel.hist, upper.panel=panel.cor))
}

rKpairs<-function(postlist,i){
  post<-postlist[[i]]
  rchains<-MCMCchains(  object=post,  params = "r",  excl = NULL,
                        ISB = TRUE,  #ignore square brackets
                        mcmc.list = FALSE,chain_num = NULL
  )
  
  Kchains<-MCMCchains(  object=post,  params = "K",  excl = NULL,
                        ISB = TRUE,  #ignore square brackets
                        mcmc.list = FALSE,chain_num = NULL
  )
  
  rK<-cbind(rchains,Kchains)
  
  betterPairs(data.frame(rK))
  mtext(names(postlist)[i], side=3, padj=-4)
}


estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

################################################################################
################################################################################

