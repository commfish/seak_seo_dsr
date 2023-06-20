###########################################################################
## Bootstrap function for dealing with ROIV lengths
## bootstrap resamples lengths using sd of estimate
## resamples those new lengths with replacement
## estimaes weight from those lengths using prescribe lm of logW from logL from port (or other) data set
##########################################################################
library(truncnorm)

w.from.l<-function (LM, length){
	syx<-summary(LM)$sigma
	cf<-exp((syx^2)/2)
	pred.log<-predict(LM,data.frame(logL=log(length)), interval = "c")
	bias.pred.orig <- exp(pred.log) 
	pred.orig<-cf*bias.pred.orig
	pred.orig[1]
}

################################################################################
par(mfrow = c(4,1))

plot(density(Port.SSEO.last3$Length.Millimeters ), xlim=c(0,1000),
	xlab="Length", 
	col="goldenrod", main="Length Dist - basic boot")

for (i in 1:nboot){	#i<-1
	ROV.rand<-sample(ROV.good$Length, size=length(ROV.good$randL),replace=TRUE)
	lines(density(ROV.rand),col=adjustcolor("blue",alpha.f=0.1), lty=1)
	abline(v=mean(ROV.rand), col=adjustcolor("blue",alpha.f=0.1))
}

polygon(density(Port.SSEO.last3$Length.Millimeters ), col=adjustcolor("red",alpha.f=0.25))
polygon(density(ROV.L$Length), col=adjustcolor("forestgreen",alpha.f=0.25))
abline(v=mean(Port.SSEO.last3$Length.Millimeters), col="red", lwd=3)
abline(v=mean(ROV.L$Length), col="forestgreen", lwd=3)

mtext("SSEO Port Lengths", 
	side=3, col="red",adj=1,padj=1.5)
mtext("ROV Lengths - point est", 
	side=3, col="forestgreen",adj=1,padj=3)
mtext("ROV bootstrap2 lengths", 
	side=3, col="blue",adj=1,padj=4.5)

###############################################################################
nboot<-500
plot(density(Port.SSEO.last3$Length.Millimeters ), xlim=c(0,1000),
	xlab="Length",col="goldenrod", main="Length Dist - boot ind. lengths on mean and sd")

for (i in 1:nboot){	#i<-1
	ROV.good$randL<-rtruncnorm(n=nrow(ROV.good),a=0,b=Inf,mean=ROV.good$Length, sd=ROV.good$L.prec)
	lines(density(ROV.good$randL),col=adjustcolor("blue",alpha.f=0.1), lty=1)
	abline(v=mean(ROV.good$randL), col=adjustcolor("blue",alpha.f=0.1))
}

polygon(density(Port.SSEO.last3$Length.Millimeters ), col=adjustcolor("red",alpha.f=0.25))
polygon(density(ROV.L$Length), col=adjustcolor("forestgreen",alpha.f=0.5))
abline(v=mean(Port.SSEO.last3$Length.Millimeters), col="red", lwd=3)
abline(v=mean(ROV.L$Length), col="forestgreen", lwd=3)

mtext("SSEO Port Lengths", 
	side=3, col="red",adj=1,padj=1.5)
mtext("ROV Lengths - point est", 
	side=3, col="forestgreen",adj=1,padj=3)
mtext("ROV bootstrap1 lengths", 
	side=3, col="blue",adj=1,padj=4.5)

##################################################################################
plot(density(Port.SSEO.last3$Length.Millimeters ), xlim=c(0,1000),
	xlab="Length", 
	col="goldenrod", main="Double Boot Lengths")

for (i in 1:nboot){	#i<-1
	ROV.good$randL<-rtruncnorm(n=nrow(ROV.good),a=0,b=Inf,mean=ROV.good$Length, sd=ROV.good$L.prec)
	ROV.rand<-sample(ROV.good$randL, size=length(ROV.good$randL),replace=TRUE)
	lines(density(ROV.rand),col=adjustcolor("blue",alpha.f=0.1), lty=1)
	abline(v=mean(ROV.rand), col=adjustcolor("blue",alpha.f=0.1))
}

polygon(density(Port.SSEO.last3$Length.Millimeters ), col=adjustcolor("red",alpha.f=0.25))
polygon(density(ROV.L$Length), col=adjustcolor("forestgreen",alpha.f=0.25))
abline(v=mean(Port.SSEO.last3$Length.Millimeters), col="red", lwd=3)
abline(v=mean(ROV.L$Length), col="forestgreen", lwd=3)

mtext("SSEO Port Lengths", 
	side=3, col="red",adj=1,padj=1.5)
mtext("ROV Lengths - point est", 
	side=3, col="forestgreen",adj=1,padj=3)
mtext("ROV bootstrap2 lengths", 
	side=3, col="blue",adj=1,padj=4.5)

##################################################################################
## Look at proportion > 270 in bootstraps
CO<-270

for (i in 1:nboot){	#i<-4
	ROV.rand<-sample(ROV.good$Length, size=length(ROV.good$randL),replace=TRUE)
	prop.gt<-length(ROV.rand[ROV.rand > CO])/length(ROV.rand)

	if (i == 1){
		Props.bas<-prop.gt
	} else {
		Props.bas<-rbind(Props.bas, prop.gt)
	}
}

for (i in 1:nboot){	#i<-1
	ROV.good$randL<-rtruncnorm(n=nrow(ROV.good),a=0,b=Inf,mean=ROV.good$Length, sd=ROV.good$L.prec)
	ROV.rand<-ROV.good$randL
	prop.gt<-length(ROV.rand[ROV.rand > CO])/length(ROV.rand)

	if (i == 1){
		Props.lengthboot<-prop.gt
	} else {
		Props.lengthboot<-rbind(Props.lengthboot, prop.gt)
	}
}

for (i in 1:nboot){	#i<-1
	ROV.good$randL<-rtruncnorm(n=nrow(ROV.good),a=0,b=Inf,mean=ROV.good$Length, sd=ROV.good$L.prec)
	ROV.rand<-sample(ROV.good$randL, size=length(ROV.good$randL),replace=TRUE)

	prop.gt<-length(ROV.rand[ROV.rand > CO])/length(ROV.rand)

	if (i == 1){
		Props.double<-prop.gt
	} else {
		Props.double<-rbind(Props.lengthboot, prop.gt)
	}
}

plot(density(Props.bas), xlim=c(0.9,1),
	xlab="Proportion", 
	col="goldenrod", main="Proportion over cutoff (270)")

polygon(density(Props.bas), col=adjustcolor("red",alpha.f=0.25))
polygon(density(Props.lengthboot), col=adjustcolor("forestgreen",alpha.f=0.25))
polygon(density(Props.double), col=adjustcolor("blue",alpha.f=0.25))
abline(v=mean(Props.bas), col="red", lwd=1)
abline(v=mean(Props.lengthboot), col="forestgreen", lwd=3)
abline(v=mean(Props.double), col="blue", lwd=1, lty=2)

mtext("Basic boot", 
	side=3, col="red",adj=0,padj=1.5)
mtext("Length mean + sd boot", 
	side=3, col="forestgreen",adj=0,padj=3)
mtext("Double boot", 
	side=3, col="blue",adj=0,padj=4.5)

################################################################################
### Check out how weights get changed in bootstrap...
################################################################################
Port.SSEO.last3$logL <- log(Port.SSEO.last3$Length.Millimeters)
Port.SSEO.last3$logW <- log(Port.SSEO.last3$Weight.Kilograms)
LWreg<-lm(logW~logL,data=Port.SSEO.last3)
maxW<-w.from.l(lm.P,length=max(ROV.good$Length)+2*max(ROV.good$L.prec))
par(mfrow = c(3,1))

#############################################################################
plot(density(Port.SSEO.last3$Weight.Kilograms), col="darkorange", 
	xlab="Weight (kg)", main="Weights - basic boot", xlim=c(0,maxW))

for (i in 1:nboot){	#i<-1
	ROV.rand<-sample(ROV.good$Length, size=length(ROV.good$Length),replace=TRUE)

	Ws<-sapply(ROV.rand, function(x)w.from.l(LM=LWreg, length=x))	
	
	lines(density(Ws),col=adjustcolor("blue",alpha.f=0.1), lty=1)
	abline(v=mean(Ws), col=adjustcolor("blue",alpha.f=0.1))
}

polygon(density(Port.SSEO.last3$Weight.Kilograms ), col=adjustcolor("red",alpha.f=0.25))
abline(v=mean(Port.SSEO.last3$Weight.Kilograms), col="red", lwd=3)

ROV.pw<-sapply(ROV.good$Length, function(x)w.from.l(LM=lm.P, length=x))
polygon(density(ROV.pw), col=adjustcolor("forestgreen",alpha.f=0.25))
abline(v=mean(ROV.pw), col="forestgreen", lwd=2, lty=2)

mtext("SSEO Port Weights", 
	side=3, col="red",adj=1,padj=1.5)
mtext("ROV point weights", 
	side=3, col="forestgreen",adj=1,padj=3)
mtext("ROV Bootstrapped weights", 
	side=3, col="blue",adj=1,padj=4.5)

##################################################################################

plot(density(Port.SSEO.last3$Weight.Kilograms), col="darkorange", 
	xlab="Weight (kg)", main="Weights - boot individual lengths", xlim=c(0,maxW))

for (i in 1:nboot){	#i<-1
	ROV.good$randL<-rtruncnorm(n=nrow(ROV.good),a=0,b=900,mean=ROV.good$Length, sd=ROV.good$L.prec)
	ROV.rand<-ROV.good$randL

	Ws<-sapply(ROV.rand, function(x)w.from.l(LM=LWreg, length=x))	
	
	lines(density(Ws),col=adjustcolor("blue",alpha.f=0.1), lty=1)
	abline(v=mean(Ws), col=adjustcolor("blue",alpha.f=0.1))
}

polygon(density(Port.SSEO.last3$Weight.Kilograms ), col=adjustcolor("red",alpha.f=0.25))
abline(v=mean(Port.SSEO.last3$Weight.Kilograms), col="red", lwd=3)

ROV.pw<-sapply(ROV.good$Length, function(x)w.from.l(LM=lm.P, length=x))
polygon(density(ROV.pw), col=adjustcolor("forestgreen",alpha.f=0.25))
abline(v=mean(ROV.pw), col="forestgreen", lwd=2, lty=2)

mtext("SSEO Port Weights", 
	side=3, col="red",adj=1,padj=1.5)
mtext("ROV point weights", 
	side=3, col="forestgreen",adj=1,padj=3)
mtext("ROV Bootstrapped weights", 
	side=3, col="blue",adj=1,padj=4.5)

##################################################################################
plot(density(Port.SSEO.last3$Weight.Kilograms), col="darkorange", 
	xlab="Weight (kg)", main="Weights - double boot", xlim=c(0,maxW))

for (i in 1:nboot){	#i<-1
	ROV.good$randL<-rtruncnorm(n=nrow(ROV.good),a=0,b=Inf,mean=ROV.good$Length, sd=ROV.good$L.prec)
	ROV.rand<-sample(ROV.good$randL, size=length(ROV.good$randL),replace=TRUE)

	Ws<-sapply(ROV.rand, function(x)w.from.l(LM=LWreg, length=x))	
	
	lines(density(Ws),col=adjustcolor("blue",alpha.f=0.1), lty=1)
	abline(v=mean(Ws), col=adjustcolor("blue",alpha.f=0.1))
}

polygon(density(Port.SSEO.last3$Weight.Kilograms ), col=adjustcolor("red",alpha.f=0.25))
abline(v=mean(Port.SSEO.last3$Weight.Kilograms), col="red", lwd=3)

ROV.pw<-sapply(ROV.good$Length, function(x)w.from.l(LM=lm.P, length=x))
polygon(density(ROV.pw), col=adjustcolor("forestgreen",alpha.f=0.25))
abline(v=mean(ROV.pw), col="forestgreen", lwd=2, lty=2)

mtext("SSEO Port Weights", 
	side=3, col="red",adj=1,padj=1.5)
mtext("ROV point weights", 
	side=3, col="forestgreen",adj=1,padj=3)
mtext("ROV Bootstrapped weights", 
	side=3, col="blue",adj=1,padj=4.5)

##################################################################################
## just looking at l:w regressions....

plot(Port.SSEO.last3$Weight.Kilograms~Port.SSEO.last3$Length.Millimeters, col="red",
	xlim=c(200,1200), ylim=c(0,30), xlab="length", ylab="weight", main="W:L")
for (i in 1:nboot){	#i<-1
	ROV.good$randL<-rtruncnorm(n=nrow(ROV.good),a=0,b=Inf,mean=ROV.good$Length, sd=ROV.good$L.prec)
	ROV.rand<-sample(ROV.good$randL, size=length(ROV.good$randL),replace=TRUE)
	
	Ref.rand<-round(runif(n=nrow(Port.SSEO.last3), min=1, max=nrow(Port.SSEO.last3)),0)
		LW.rand<-Port.SSEO.last3[Ref.rand,]		#hist(Ref.rand)	#str(LW.rand)	#nrow(LW.rand)
		LWreg<-lm(logW~logL,data=LW.rand)

	Ws<-sapply(ROV.rand, function(x)w.from.l(LM=LWreg, length=x))	
	points(Ws~ROV.rand, col=adjustcolor("blue",alpha.f=0.1), pch=20)
}
points(Port.SSEO.last3$Weight.Kilograms~Port.SSEO.last3$Length.Millimeters, col="red", pch=20)
mtext("SSEO Port L:W", 
	side=3, col="red",adj=1,padj=6)
mtext("ROV backtransformed L:W", 
	side=3, col="blue",adj=1,padj=7.5)


###############################################################################
plot(density(Port.SSEO.last3$Length.Millimeters ), xlim=c(0,1000),
	xlab="Length",col="goldenrod", main="Length Dist - boot ind. lengths on mean and sd")

for (i in 1:nboot){	#i<-1
	ROV.good$randL<-rtruncnorm(n=nrow(ROV.good),a=0,b=Inf,mean=ROV.good$Length, sd=ROV.good$L.prec)
	lines(density(ROV.good$randL),col=adjustcolor("blue",alpha.f=0.1), lty=1)
	abline(v=mean(ROV.good$randL), col=adjustcolor("blue",alpha.f=0.1))
}

polygon(density(Port.SSEO.last3$Length.Millimeters ), col=adjustcolor("red",alpha.f=0.25))
polygon(density(ROV.L$Length), col=adjustcolor("forestgreen",alpha.f=0.5))
abline(v=mean(Port.SSEO.last3$Length.Millimeters), col="red", lwd=3)
abline(v=mean(ROV.L$Length), col="forestgreen", lwd=3)

mtext("SSEO Port Lengths", 
	side=3, col="red",adj=1,padj=1.5)
mtext("ROV Lengths - point est", 
	side=3, col="forestgreen",adj=1,padj=3)
mtext("ROV bootstrap1 lengths", 
	side=3, col="blue",adj=1,padj=4.5)

##################################################################################
plot(density(Port.SSEO.last3$Length.Millimeters ), xlim=c(0,1000),
	xlab="Length", 
	col="goldenrod", main="Double Boot Lengths")

for (i in 1:nboot){	#i<-1
	ROV.good$randL<-rtruncnorm(n=nrow(ROV.good),a=0,b=Inf,mean=ROV.good$Length, sd=ROV.good$L.prec)
	ROV.rand<-sample(ROV.good$randL, size=length(ROV.good$randL),replace=TRUE)
	lines(density(ROV.rand),col=adjustcolor("blue",alpha.f=0.1), lty=1)
	abline(v=mean(ROV.rand), col=adjustcolor("blue",alpha.f=0.1))
}

polygon(density(Port.SSEO.last3$Length.Millimeters ), col=adjustcolor("red",alpha.f=0.25))
polygon(density(ROV.L$Length), col=adjustcolor("forestgreen",alpha.f=0.25))
abline(v=mean(Port.SSEO.last3$Length.Millimeters), col="red", lwd=3)
abline(v=mean(ROV.L$Length), col="forestgreen", lwd=3)

mtext("SSEO Port Lengths", 
	side=3, col="red",adj=1,padj=1.5)
mtext("ROV Lengths - point est", 
	side=3, col="forestgreen",adj=1,padj=3)
mtext("ROV bootstrap2 lengths", 
	side=3, col="blue",adj=1,padj=4.5)
################################################################################


ROV.L$badL<-ROV.L$Length-3*ROV.L$L.prec
ROV.goodL<-ROV.L[ROV.L$badL > 0,]
nrow(ROV.L)
nrow(ROV.goodL)
hist(ROV.goodL$L.cv) 

nboot<-1000
O<-data.frame(NA)
pb = txtProgressBar(min = 0, max = length(nboot), initial = 0, style=3) 
for (i in 1:nboot){	#i<-1
	#get lengths by randomly sampling from mean and sd of measured ROV fish...
	ROV.goodL$randL<-rtruncnorm(n=nrow(ROV.goodL),a=0,b=Inf,mean=ROV.goodL$Length, sd=ROV.goodL$L.prec)		#head(ROV.goodL,10)
	#Next take random sample of those lengths... so now we've accounted for inaccuracy of measurements and randomness of sampling!
	ROV.rand<-sample(ROV.goodL$randL, size=length(ROV.goodL$randL),replace=TRUE)
	#get the mean length from that sample and record mean and var
	O$mean.l.all<-mean(ROV.rand)
	O$var.l.all<-var(ROV.rand)

	#get and record the proportion of lengths greater than CUTOFF value... 
	Prop.gt<-length(ROV.rand[ROV.rand>CUTOFF])/length(ROV.rand)
	O$Prop.gt<-Prop.gt
	
	#get mean length of fish greater than CUTOFF value for biomass calculations... 
	ROV.gt.co<-ROV.rand[ROV.rand>CUTOFF]
	O$mean.l.gtco<-mean(ROV.gt.co)
	O$var.l.gtco<-var(ROV.gt.co)

	#now lets get 2 things : the weight of the average fish and the average or the converted weight of each fish... (very different!)
	#first lets get the weight of the average fish... 
	O$W.mean.L<-w.from.l(LM=LW.lm , length=mean(ROV.rand))
	O$W.mean.L.gtco<-w.from.l(LM=LW.lm , length=mean(ROV.gt.co))
	#now weight each ROV fish gt cutoff...
	Ws<-sapply(ROV.rand[ROV.rand > 0], function(x)w.from.l(LM=LW.lm, length=x))	#plot(Ws~ROV.rand[ROV.rand>CUTOFF])
	mean(Ws)
	O$mean.W<-mean(Ws)					
	O$mean.W.var<-var(Ws)
	Ws.gtco<-sapply(ROV.gt.co, function(x)w.from.l(LM=LW.lm, length=x))	#plot(Ws~ROV.rand[ROV.rand>CUTOFF])
	O$mean.W.gtco<-mean(Ws.gtco)					
	O$mean.W.gtco.var<-var(Ws.gtco)

	if (i == 1) {
	ROV.boot<-O
	} else {
	ROV.boot<-rbind(ROV.boot,O)
	}
setTxtProgressBar(pb,i/nboot)
}
