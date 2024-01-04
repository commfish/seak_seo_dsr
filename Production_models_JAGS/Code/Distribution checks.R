YY<-rlnorm(10000,log(1000),0.5); hist(YY, breaks=100); max(YY); min(YY); mean(YY)

rr <- rlnorm(10000,-1.38,3.845); rr<-rr[rr>0.01 & rr<1.2]; hist(rr, breaks=100); max(rr); min(rr); mean(rr)
rr <- rlnorm(10000,-1.38,1.5); rr<-rr[rr>0.01 & rr<0.5]; hist(rr, breaks=100); max(rr); min(rr); mean(rr)
rrr<-runif(10000,0.01,0.3); hist(rrr)

KK <- rlnorm(10000,log(200000),1); KK<-KK[KK>40000 & KK<1000000];hist(KK, breaks=100); max(KK); min(KK); mean(KK)

rr2 <- rlnorm(10000,log(0.2),0.5);rr2<-rr2[rr2>0.01 & rr2<0.4]; hist(rr2, breaks=100); max(rr2); min(rr2); mean(rr2)

logr<-rnorm(10000,0,1.0E-6); hist(logr, breaks=100)
r<-exp(logr); hist(r, breaks=100)

C<- rlnorm(10000,log(1.2*CK[t]),0.5); hist(C, breaks=100) 

CV<-rlnorm(100000,log(0.17),0.1); hist(CV, breaks=100); max(CV)

phi<-rnorm(100000,0.6,0.4); phi<-phi[phi < 0.9 & phi > 0]; hist(phi, breaks=50)

qB<-rnorm(10000,1,0.2); hist(qB)  #qB<-qB[qB>0 & qB < 2]; 
qB3<-rnorm(10000,)

qq<-14
exp(qq); exp(qq)*.03



qBl<-rlnorm(10000,0,0.15); hist(qBl)
qB2<-rlnorm(10000,0,0.1); hist(qB2, xlab="q prior", ylab="",main="")
qB3<-rlnorm(10000,log(1000000),0.1); hist(qB3)


p<-rlnorm(10000,log(0.18815), 0.05); hist(p)

plot(qB); points(qBl,col="red")
log(1)

X<-rnorm(20,20,1)
X[1:10]

DD<-vector()
C<-vector()
D<-vector()
logD<-vector()
for (t in 1:N){   t<-N[1]
  DD[t]<-Dobs[t]
  C[t] <-rlnorm(1,log(DD[t]+CK[t]),0.5)  #Total catches
  D[t]<-C[t]-CK[t]                       #Discards=total minus known (CK)
  logD[t]<-log(D[t])
  Dobs[t] ~ dlnorm(logD[t], tau.log.D[t])  #Observed discard estimate
  tau.log.D[t]<-1/log(cv.D[t]*cv.D[t]+1)   #estimated cvs of discards
}

#==========================================================================
#Process error
sss<--3
logvar <- runif(100000,-10,sss)
sigma <- sqrt(exp(logvar)) ; mean(sigma)  #as per Ono et al. 
hist(sigma, breaks=100, main=sss); hist(sigma^2, breaks=100, col="blue")
max(sigma^2)

eps<-rnorm(10000,0,max(sigma^2)); hist(eps, breaks=100,col="orange")
hist(exp(eps),breaks=100)

prec<-1/(0.3*0.3+1)

logvar<-runif(10000,-10,3)
sigma<-sqrt(exp(logvar)) # ; plot(density(sigma))

epsilon<-rnorm(10000,0,sigma^2); epsilon<-epsilon[epsilon > -0.1 & epsilon < 0.1]
lines(density(epsilon)) #plot(density(epsilon)); max(epsilon);min(epsilon)

taulog<-1/(log(sigma*sigma+1)); plot(density(taulog))
taulog2<-(1/(sigma))^2; plot(density(taulog2))

logq <- runif(10000,-10,10)
q <- (exp(logq))   #as per Ono et al. 
hist(q, breaks=100)
min(q)
max(q)

#==========
#Epsilon log check
Bt<-23; e<-0.1; sig<-0.3
log(Bt*exp(e-sig))
log(Bt)+e-sig

#==================================================================
# r priors
getmode <- function(v) {
  rnds<-round(v,3)
  uniqv <- unique(rnds)
  uniqv[which.max(tabulate(match(rnds, uniqv)))]
}

gammas<-function(shape,rate){
  rphypgam<-rgamma(5000000,shape=shape,rate=rate)
  hist(rphypgam,breaks=100,col="orange")
  abline(v=mean(rphypgam), col="blue", lwd=1.5)
  abline(v=getmode(rphypgam), col="green")
  abline(v=median(rphypgam), col="red")
  print("mean");print(mean(rphypgam))
  print("mode");print(getmode(rphypgam))
  print("median");print(median(rphypgam))
  print("variance"); print(var(rphypgam))
}

lognorms<-function(mean,sd){
  rphypln<-rlnorm(5000000,log(mean),sd)
  hist(rphypln[rphypln<0.2],
       breaks=100,col="orange")
  abline(v=mean(rphypln), col="blue", lwd=1.5)
  abline(v=getmode(rphypln), col="green")
  abline(v=median(rphypln), col="red")
  print("mean");print(mean(rphypln))
  print("mode");print(getmode(rphypln))
  print("median");print(median(rphypln))
  print("variance"); print(var(rphypln))
}

gammas(2,25)
gammas(2,24.3) #<-mode ~0.04, var=0.0034
gammas(2.5,38) #mode ~0.04, var=0.0017

gammas(3,50) #<-mode ~0.04, var = 0.0012

gammas(5,100) #mode ~0.04, var=0.0005

gammas(1.5,12.75) #mode ~0.04, var=0.009

gammas(1.25,6) #mode ~0.04, var=0.035

par(mfrow=c(3,1))
lognorms(0.04,0.5)
lognorms(0.04,1)
lognorms(0.10,1)
lognorms(0.05,0.5)

#sub r prior...
rs<-rlnorm(100000,log(0.04),1); rs<-rs[rs<0.1]
hist(rs,breaks=100,col="green")

rs<-rnorm(100000,0.04,0.1); mean(rs) #; rs<-rs[rs<0.1]
hist(rs,breaks=100,col="green")

#K priors
logKp<-runif(100000,7.5,11.5); hist(exp(logKp), breaks=100)
exp(7)
exp(mean(logKp))  
  
rsigs<-rnorm(100000,0.5,0.2); hist(rsigs)
rsigs<-rlnorm(100000,log(0.5),1); hist(rsigs)
  
#-----------------------------------------------
logn<-log(100)
logetas<-rlogis(100000,logn,1); hist(logetas); max(logetas); min(logetas)
logetas<-logetas[logetas > -5 & logetas < 7.55]

etas<-exp(logetas); hist(etas, breaks=100); max(etas); min(etas)

shrinkage<-etas/(etas+100)
hist(shrinkage); quantile(shrinkage,c(0.01,0.025,0.5,0.0975,0.99,1))

shrinkeg<-409/(409+100)

R.hyp<-rbeta(length(etas),1,1)
rB1<-R.hyp*etas; plot(density(rB1)); max(rB1); min(rB1)
rB2<-(1-R.hyp)*etas; hist(rB2); max(rB2); min(rB2)

#-------------------------------------------------------------------
#Prior for foreign fleet
meanFor<-mean(For[For>1])
hist(For[For>1],breaks=25)

plot(density(For[For>1]),lwd=2, xlim=c(0,1300))
ff<-rlnorm(100000,log(meanFor),log(2.5)); ff<-ff[ff < 1300]; lines(density(ff),col="blue")  #hist(ff,breaks=100)
ff<-rlnorm(100000,log(150),log(2.5)); ff<-ff[ff < 1300]; lines(density(ff),col="red")  #hist(ff,breaks=100)

#-----------------------------------------------------------------
#hyper sig for normal 
sigs<-rgamma(100000,1,1); plot(density(sigs))

invtau2 <- rgamma(100000, 1, 1) #with ones is uninformative... 
X<-6
trunc<- invtau2[invtau2 < X]
plot(density(invtau2)); lines(density(trunc), col="red")
  
phiTau<-sqrt(invtau2^-1)
ttau<-sqrt(trunc^-1)

plot(density(phiTau)); lines(density(ttau), col="red")
#over 1K very rare... truncate

Y<-phiTau[phiTau < 1000]
lines(density(Y), col="blue")

max(phiTau)
max(ttau)
max(Y)


## priors
for (j in 1:J){
  mu_j[j] ~ dnorm(mu, invtau2)
}
invsigma2 ~ dgamma(a_s, b_s)
sigma <- sqrt(pow(invsigma2, -1))
## hyperpriors
mu ~ dnorm(mu0, g0)
invtau2 ~ dgamma(a_t, b_t)
tau <- sqrt(pow(invtau2, -1))
}
"










