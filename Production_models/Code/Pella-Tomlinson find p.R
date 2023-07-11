
#relationship of z (or p) to Bmsy:K ratio
G<-function(z){
  sol<-exp((1/z)*log(1/(z+1)))
  return(sol)
}

G(2.39) #whale example in class where ratio is 0.6
G(2)
G(1.6)
G(1)    #sets at Bmsy:K at 0.5 which is same as base Schaefer model
G(0.5)
G(0.43)
G(0.4)
G(0.3)
G(0.2)
G(0.18815) #this sets Bmsy:K at 0.4
G(0)
G(0.00001)
G(1/10000)
G(1/100000)
G(1/1000000)
G(1/1000000000000000)
G(1/1000000000000000)
G(-0.09338958436888) #this sets Bmsy:K at 0.35; see www.wolframalpha.com



G(-0.09338958436888)

r<- 0.032
K<- 41497
p<--0.09338958436888
MSY<-r*K/((p+1)^((p+1)/p))

pop<-seq(0,42000,250)
Surplus<-(r/p)*pop*(1-(pop/K)^p) 
pop.pro<-as.data.frame(cbind(pop,Surplus))
pop.pro[pop.pro$Surplus == max(Surplus, na.rm=T),]
pop.pro$pop[pop.pro$Surplus == max(Surplus, na.rm=T)]/K

plot(Surplus ~ pop)

Bmsy<-0.4*K

Bmsy[i]<-0.4*K[i]  #0.5 for Schaefer
Fmsy[i]<-MSY[i]/Bmsy[i]
Hmsy[i]<-r[i]/(1+p) 
Stock.Status[i]<-B[i,N]/(0.4*K[i])