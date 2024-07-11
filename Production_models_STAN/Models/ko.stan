// Multiarea version

data {
  int Npre; // number of years with cpue
  int N; // number of years of catch data
  int S; // number of strata 
  real p; // exponent of PT
  real C_obs[S, N-1]; // catch
  real I_obs[S, Npre]; // indices
  real ratio; 
}

transformed data {
  real n = p + 1;
  real bb = pow(n,(n/(n-1)));
}

parameters {
  // real<lower=0> MSY[S];
  real<lower=0, upper=2> r[S];
  real<lower=0, upper=1> PP_init[S];
  real K[S];
  real iqs[S];
  real<lower=0> isigma2;
  real PE[S, N-1];  
}

transformed parameters {
  matrix[S, N] logB;  
  matrix[S, N] B;  
  real logqs[S];
  real sigma_proc = 1/sqrt(isigma2);
  real sigma_obs = ratio * sigma_proc;
	for (i in 1:S){
		logqs[i] = log(1/iqs[i]);
		logB[i,1] =  log(K[i]); //logB[i,1] =  log(PP_init[i]) + log(K[i]);
		B[i,1] = exp(logB[i,1]);
	}

	// when we have the disaggregated data
	for (t in 1:(N-1))
	{    
		for (i in 1:S)
		{
      logB[i,t+1] = log( fmax(B[i,t]+(r[i]/p)*B[i,t]*(1-pow(B[i,t]/K[i],p))-C_obs[i,t],1) *exp(PE[i,t]) );
			B[i,t+1] = exp(logB[i,t+1]);
		}
	}
}

model {
	
// Priors
	isigma2 ~ gamma(100, 0.25);

	for (i in 1:S){
		K[i] ~ lognormal(12, 0.5);
		r[i] ~ lognormal(log(0.1), 0.2);
		iqs[i] ~ uniform(1000,100000);
		PP_init[i] ~ normal(0.99,0.01);
  }

// Likelihoods process error
	for (t in 1:(N-1)) 
	{   
		for (i in 1:S)
		{
		  PE[i,t] ~ normal(-0.5*sigma_proc*sigma_proc, sigma_proc); 
	  }
	}
	
		
// Likelihoods observation error
	for (i in 1:S)
	{
		for (t in 1:Npre) 
		{       
			I_obs[i,t] ~ lognormal((logqs[i]+logB[i,t+(N-Npre)]-0.5*sigma_obs*sigma_obs), sigma_obs); 
		}
	}
} // end model

generated quantities {
	matrix[S, N] Iest;
 matrix[S, N] Imed;
 real qs[S];
	for (s in 1:S) {
		qs[s] = exp(logqs[s]);
    for (t in 1:N) {
      Imed[s,t] = log(qs[s] * B[s,t]);
      Iest[s,t] = exp(Imed[s,t]);
   }
  }

  real MSY[S];
  real Bmsy[S];
  real Fmsy[S];
  real Stock_status[S];

  for (s in 1:S){
   MSY[s] = r[s] * K[s] / (pow(p + 1,(p + 1) / p));
   Bmsy[s] = 0.4 * K[s];
   Fmsy[s] = MSY[s] / Bmsy[s];
   Stock_status[s] = B[s,N] / K[s];
  }
}
