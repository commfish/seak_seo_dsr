data {
  int<lower=0> N; // number of years 
  real C[N]; // estimated treatment effects
  real I[N]; // estimated treatment effects
}
parameters {
  real<lower=0> r;
  real<lower=0> K;
  real<lower=0> isigma2;
  real<lower=0> itau2;
  vector[N] P;
  //vector[N] Z;
}

transformed parameters {
 real sigma; real tau; 
 sigma = 1.0/sqrt(isigma2);
 tau = 1.0/sqrt(itau2);
}

model {
  vector[N] Pmed;
  vector[N] Zmed;
  real qhat; 
  
  //priors
  r ~  beta(1,1) T[0.01,0.99];;// r ~ beta(1,1) T[0.01,0.99]; //lognormal(-1.38,3.845) T[0.01,1.2]
  isigma2 ~ gamma(6, 0.115); //gamma(3.785, 0.0102);
  itau2 ~ gamma(1.709, 0.00861); //gamma(1.709, 0.00861);
  K ~ lognormal(5.0429, 3.7603) T[10,1500];
  
  // Set initial state 
  Pmed[1] = 0;
  P[1] ~ lognormal(Pmed[1],sigma) T[0.05,1.6]; 
  
  // time steps of the model
  for (t in 2:N)
   {
    Pmed[t] = log(fmax(P[t - 1] + (r * P[t - 1]) * (1 - P[t - 1]) - C[t - 1] / K, 0.001) );
    P[t] ~ lognormal(Pmed[t], sigma) T[0.05, 1.6];
   }
  
  // Likelihood
  for (t in 1:N)
   {
    Zmed[t] = log(I[t] / (P[t]  * K));
   }    

   qhat = sum(Zmed) / N;  //qhat = sum(Z)/N;

  for (t in 1:N)
   {
    Zmed[t] ~ normal(qhat,tau);
   }  

}  // end model

generated quantities {
  vector[N] Imed; 
  vector[N] Inew; //for posterior predictive checks... 
  vector[N] Biomass;
  vector[N] Zmed;
  real qhat;
  real q;
  real MSY;
  real EMSY;

  // Other ouputs
   for (t in 1:N)
   {
    Zmed[t] = log(I[t] / (P[t]  * K));
   }  
   qhat = sum(Zmed) / N;
   q = exp(qhat);
   MSY = r*K/4;
   EMSY = r/(2*q);

  //posterior predictions (hint, the parameterization of dlnorm is not the same as in R)
  for (t in 1:N)
   {
    Imed[t] = log((q * K) * P[t]);
    Inew[t] = lognormal_rng(Imed[t], tau);
    Biomass[t] = K*P[t];
   }
  }
  
  



