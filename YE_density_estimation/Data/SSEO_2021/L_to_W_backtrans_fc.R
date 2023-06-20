###########################################################################
## Function for backtransforming predicted weight from LW regression
## YE SSEO ROV project
## Phil Joy July 2021
###########################################################################
##LM is the linear model of logW versus logL
## length is the length you want to back transform

w.from.l<-function (LM, length){
	syx<-summary(LM)$sigma
	cf<-exp((syx^2)/2)
	pred.log<-predict(LM,data.frame(logL=log(length)), interval = "c")
	bias.pred.orig <- exp(pred.log) 
	pred.orig<-cf*bias.pred.orig
	pred.orig[1]
}

#w.from.l(LM=lm.Port, length=600)
