quants4<-array(dim=c(3,(Nyear-1),Narea))
for (i in 1:(Nyear-1)){
  for (s in 1:Narea){
    quants4[,i,s] <- quantile(samples$PE[,s,i],probs=,c(0.05,0.5,0.95))
  }
}


mean(quants4[2,,1])
mean(quants4[2,,2])
mean(quants4[2,,3])
median(quants4[2,,])

mean(epsilon)
mean(epsilon[,1])
mean(epsilon[,2])
mean(epsilon[,3])

mean(epsilon[,1])
mean(epsilon[,2])
mean(epsilon[,3])

pe_bias_1 <- (mean(quants4[2,,1]) - mean(epsilon[,1]))/mean(epsilon[,1])
pe_bias_2 <- (mean(quants4[2,,2]) - mean(epsilon[,2]))/mean(epsilon[,2])
pe_bias_3 <- (mean(quants4[2,,3]) - mean(epsilon[,3]))/mean(epsilon[,3])

pe_bias_1 <- (mean(abs(quants4[2,,1])) - mean(abs(epsilon[,1])))/mean(abs(epsilon[,1]))
pe_bias_2 <- (mean(abs(quants4[2,,2])) - mean(abs(epsilon[,2])))/mean(abs(epsilon[,2]))
pe_bias_3 <- (mean(abs(quants4[2,,3])) - mean(abs(epsilon[,3])))/mean(abs(epsilon[,3]))


quants44 <- data.frame()
for (s in 1:Narea) {
  ssamp <- t(quants4[,,s]) %>% data.frame() %>% mutate(strata = s,
                                                       year = seq(1,(Narea-1),1))
  colnames(ssamp) <- c("lo90","median","hi90","strata","year")
  if (s == 1) {
    quants44 <- ssamp
  } else {
    quants44 <- rbind(quants44,ssamp)
  }
}

quants44$true_year <- rep(true_years,strata)

init_ll_3r[[2]]$PE <- matrix(runif(3*49,-0.01,0.01), nrow=data1$S, ncol=Nyear - 1)

init_ll_3r[[2]]$iqs <- init_ll_3r[[1]]$iqs

init_ll_3r[[1]]$r
init_ll_3r[[2]]$r
init_ll_3r[[3]]$r

mean(init_ll_3r[[1]]$K)
mean(init_ll_3r[[2]]$K)
mean(init_ll_3r[[3]]$K)

init_ll_3r[[2]]$r[3] <- init_ll_3r[[3]]$r[3]
init_ll_3r[[2]]$K[3] <- init_ll_3r[[3]]$K[3]
init_ll_3r[[2]]$iqs[3] <- init_ll_3r[[3]]$iqs[3]

?runif

crap_inits <- init_ll_3r
init_ll_3r <- crap_inits

hist(exp(rnorm(10000,12,0.5)))
