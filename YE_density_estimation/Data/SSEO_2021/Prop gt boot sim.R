Length.true<-ROV.good$Length
Length.cv<-ROV.good$L.prec
Length.samp<-rtruncnorm(n=nrow(ROV.good),a=0,b=Inf,mean=ROV.good$Length, sd=ROV.good$L.prec)

True.prop<-length(Length.true[Length.true > 270])/length(Length.true)
True.se<-sqrt(True.prop*(1-True.prop)/length(Length.true))
False.prop<-length(Length.samp[Length.samp> 270])/length(Length.samp)
False.se<-sqrt(False.prop*(1-False.prop)/length(Length.samp))

CO<-270
nboot=1000
for (i in 1:nboot){	#i<-4
	ROV.rand<-sample(Length.samp, size=length(Length.samp),replace=TRUE)
	prop.gt<-length(ROV.rand[ROV.rand > CO])/length(ROV.rand)

	if (i == 1){
		Props.bas<-prop.gt
	} else {
		Props.bas<-rbind(Props.bas, prop.gt)
	}
}

for (i in 1:nboot){	#i<-1
	R2<-rtruncnorm(n=length(Length.samp),a=0,b=Inf,mean=Length.samp, sd=Length.cv)
	ROV.rand<-R2
	prop.gt<-length(ROV.rand[ROV.rand > CO])/length(ROV.rand)

	if (i == 1){
		Props.lengthboot<-prop.gt
	} else {
		Props.lengthboot<-rbind(Props.lengthboot, prop.gt)
	}
}

for (i in 1:nboot){	#i<-1
	R2<-rtruncnorm(n=length(Length.samp),a=0,b=Inf,mean=Length.samp, sd=Length.cv)
	ROV.rand<-sample(R2, size=length(R2),replace=TRUE)

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
abline(v=True.prop, col="orange")
mtext("Basic boot", 
	side=3, col="red",adj=0,padj=1.5)
mtext("Length mean + sd boot", 
	side=3, col="forestgreen",adj=0,padj=3)
mtext("Double boot", 
	side=3, col="blue",adj=0,padj=4.5)

