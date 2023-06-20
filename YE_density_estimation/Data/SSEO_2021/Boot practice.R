library(boot)

meanfun <- function(x, d) {
  return(mean(x[d]))
  return(var(x[d]))
}

DATA<-ROV.L$randL

set.seed(1)
bo <- boot(DATA, statistic=meanfun, R=1000)
summary(bo)
bo$t0		#original mean
mean(bo$t)	#bootstrap mean
sd(bo$t)	#bootstrap sd
mean(bo$t)-bo$t0
sd(bo$t)

print(bo)
plot(bo)

bo.ci<-boot.ci(bo, conf=0.95, type="bca")
bo.ci[2]	#original mean
bo.ci$bca[1,4]	#lower CI
bo.ci$bca[1,5] 	#upper CI
bo.ci$bca[1,2]	#SD?

mean(bo$t)+1.96*sd(bo$t)



str(bo.ci)
summary(bo.ci)

CUTOFF<-270	#cutoff length

nboot<-1000
O<-data.frame(matrix(NA, ncol=5))
pb = txtProgressBar(min = 0, max = length(nboot), initial = 0, style=3) 
for (i in 1:nboot){	#i<-1
	ROV.L$randL<-rnorm(n=nrow(ROV.L),mean=ROV.L$Length, sd=ROV.L$L.prec)
	#ROV.rand1<-rnorm(n=nrow(ROV.L), mean=ROV.L[,which(colnames(ROV.L) == "Length")], sd=ROV.L[,which(colnames(ROV.L) == "L.prec")])
	ROV.rand<-sample(ROV.L$randL, size=length(ROV.L$randL),replace=TRUE)
	ROV.mean.l<-mean(ROV.rand)
	Prop.gt<-length(ROV.rand[ROV.rand>CUTOFF])/length(ROV.rand)
	O$Prop.gt<-Prop.gt
	ROV.gt.cu<-ROV.rand[ROV.rand>CUTOFF]
	ROV.mean.L.gt<-mean(ROV.gt.cu)
	O$ROV.mean.length.gt<-ROV.mean.L.gt

	O$ROV.W.mean<-w.from.l(LM=LW.lm , length=ROV.mean.L.gt)
	Ws<-sapply(ROV.rand[ROV.rand>CUTOFF], function(x)w.from.l(LM=LW.lm, length=x))	#plot(Ws~ROV.rand[ROV.rand>CUTOFF])
	O$ROV.mean.W<-mean(Ws)					
	O$ROV.mean.W.var<-var(Ws)
	if (i == 1) {
	OUT<-O
	} else {
	OUT<-rbind(OUT,O)
	}
setTxtProgressBar(pb,i/nboot)
}

Avg.ROV.Length<-mean(OUT$ROV.mean.length.gt)
Prop.gt.CU<-mean(OUT$Prop.gt)
Avg.ROV.Weight<-mean(OUT$ROV.W.mean)
Avg.ROV.W2<-mean(OUT$ROV.mean.W)
Avg.W2.var<-mean(OUT$ROV.mean.W.var)
Avg.ROV.Length
Prop.gt.CU
Avg.ROV.Weight
Avg.ROV.W2
sqrt(Avg.W2.var)


WLfun <- function(x, d) {
  return(mean(x[d]))
  return(var(x[d]))
}

set.seed(1)
bo <- boot(DATA, statistic=meanfun, R=1000)

N<-Best$Nhat
N.se<-Best$Nhat.se
N.var<-N.se^2

D<-Best$Dhat
D.se<-Best$Dhat.se
D.var<-D.se^2

AA<-1056

varN.q<-D.var+AA^2
seN.q<-sqrt(varN.q)

N.var-D.var
AA^2























