// Multiarea version

data {
  int Npre; // number of years 
  int S; // number of strata 
  real p; // exponent of PT
  real C_obs[S, Npre-1]; // catch
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
  // real<lower=0> sigma_proc;
  real PE[S, Npre-1];  
}

transformed parameters {
  matrix[S, Npre] logB;  
  matrix[S, Npre] B;  
  real logqs[S];
  real sigma_proc = 1/sqrt(isigma2);
  // matrix[S, Npre] PP;
  real sigma_obs = ratio * sigma_proc;
	for (i in 1:S){
		logqs[i] = log(1/iqs[i]);
		// PP[i,1] = PP_init[i];
		// logB[i,1] =  log(PP[i,1]) + log(K[i]);
		logB[i,1] =  log(PP_init[i]) + log(K[i]);
		B[i,1] = exp(logB[i,1]);
	}

	// when we have the disaggregated data
	for (t in 1:(Npre-1))
	{    
		for (i in 1:S)
		{
      logB[i,t+1] = log( fmax(B[i,t]+(r[i]/p)*B[i,t]*(1-pow(B[i,t]/K[i],p))-C_obs[i,t],1) *exp(PE[i,t]) );
      // B[i,t] = exp(logB[i,t])
			// PP[i,t+1] = fmax(PP[i,t] + MSY[i]/K[i]*bb*(PP[i,t]-pow(PP[i,t],n))-C_obs[i,t]/K[i],0.000001)*exp(PE[i,t]);
			// logB[i,t+1] = log(PP[i,t+1]) + log(K[i]);
			B[i,t+1] = exp(logB[i,t+1]);
		}
	}
}

model {
	
// Priors
	isigma2 ~ gamma(100, 0.25);

	for (i in 1:S){
		K[i] ~ lognormal(12, 0.5);
		// MSY[i] ~ lognormal(10, 0.5);
		r[i] ~ lognormal(log(0.05), 0.25);
		iqs[i] ~ uniform(1000,100000);
		PP_init[i] ~ normal(0.9,0.1);
  }

// Likelihoods process error
	for (t in 1:(Npre-1))
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
			I_obs[i,t] ~ lognormal((logqs[i]+logB[i,t]-0.5*sigma_obs*sigma_obs), sigma_obs); 
		}
	}
	
}
