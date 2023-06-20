ROV.L[ROV.L$L.cv > 0.4,]
ROV.goodL[ROV.goodL$L.cv > 0.3,]

mean(ROV.goodL$Length)
mean(ROV.L$Length)

plot(density(ROV.L$Length),  col="blue")
lines(density(ROV.goodL$Length), col="red")

MEAN<-326.845
SD<-105.815

Trunc<-rtruncnorm(10000,a=0,b=Inf,mean=MEAN , sd=SD)
DTrunc<-rtruncnorm(10000,a=0,b=2.5*SD+MEAN,mean=MEAN , sd=SD)

Norm<-rnorm(10000,mean=MEAN, sd= SD)

plot(density(Trunc), col="red")
lines(density(DTrunc), col="blue")
lines(density(Norm), col="green")

mean(Trunc)
mean(DTrunc)
mean(Norm)

min(Trunc)
min(DTrunc)
min(Norm)

((mean(Trunc)-MEAN)/MEAN)*100
((mean(DTrunc)-MEAN)/MEAN)*100
((mean(Norm)-MEAN)/MEAN)*100

Prior.boot<-ROV.boot

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

Prior.boot
New.boot<-ROV.boot

str(New.boot)

plot(density(Prior.boot$mean.l.all), col="red")
lines(density(New.boot$mean.l.all), col="blue")
abline(v=mean(Prior.boot$mean.l.all), col="red")
abline(v=mean(New.boot$mean.l.all), col="blue")

plot(density(Prior.boot$var.l.all), col="red")
lines(density(New.boot$var.l.all), col="blue")
abline(v=mean(Prior.boot$var.l.all), col="red")
abline(v=mean(New.boot$var.l.all), col="blue")

plot(density(Prior.boot$Prop.gt), col="red")
lines(density(New.boot$Prop.gt), col="blue")
abline(v=mean(Prior.boot$Prop.gt), col="red")
abline(v=mean(New.boot$Prop.gt), col="blue")

plot(density(Prior.boot$mean.l.gtco), col="red")
lines(density(New.boot$mean.l.gtco), col="blue")
abline(v=mean(Prior.boot$mean.l.gtco), col="red")
abline(v=mean(New.boot$mean.l.gtco), col="blue")

plot(density(Prior.boot$var.l.gtco), col="red")
lines(density(New.boot$var.l.gtco), col="blue")
abline(v=mean(Prior.boot$var.l.gtco), col="red")
abline(v=mean(New.boot$var.l.gtco), col="blue")

plot(density(Prior.boot$mean.W ), col="red")
lines(density(New.boot$mean.W ), col="blue")
abline(v=mean(Prior.boot$mean.W ), col="red")
abline(v=mean(New.boot$mean.W ), col="blue")

plot(density(Prior.boot$mean.W.var ), col="red")
lines(density(New.boot$mean.W.var ), col="blue")
abline(v=mean(Prior.boot$mean.W.var ), col="red")
abline(v=mean(New.boot$mean.W.var ), col="blue")

plot(density(Prior.boot$mean.W.gtco ), col="red")
lines(density(New.boot$mean.W.gtco ), col="blue")
abline(v=mean(Prior.boot$mean.W.gtco ), col="red")
abline(v=mean(New.boot$mean.W.gtco ), col="blue")

plot(density(Prior.boot$mean.W.gtco.var ), col="red")
lines(density(New.boot$mean.W.gtco.var ), col="blue")
abline(v=mean(Prior.boot$mean.W.gtco.var ), col="red")
abline(v=mean(New.boot$mean.W.gtco.var ), col="blue")











