// Multiarea version

data {
  int Npre; // number of years with cpue
  int N; // number of years of catch data 
  int S; // number of strata 
  real p; // exponent of PT
  real C_obs[S, N-1]; // catch
  real ratio; 

  //add in biomass
  int<lower=0> N_Bobs; //  total number of non-NA values
  real B_obs[N_Bobs];
  real B_cv[N_Bobs];
  int<lower=0> S_Bobs[S]; //number of observations per strata, corresponds to s on 133-134 Stan manual
  int<lower=0> B_pos[N_Bobs]; // index positions as one vector but relevant to each strata

  //add index with cvs
  int<lower=0> N_Iobs; //  total number of non-NA values
  real I_obs[N_Iobs];
  real I_cv[N_Iobs];
  int<lower=0> S_Iobs[S]; //number of observations per strata, corresponds to s on 133-134 Stan manual
  int<lower=0> I_pos[N_Iobs]; // index positions as one vector but relevant to each strata
}

transformed data {
  real n = p + 1;
  real bb = pow(n,(n/(n-1)));
}

parameters {
  // real<lower=0> MSY[S];
  real<lower=0, upper=2> r;
  real<lower=0, upper=1> PP_init[S];
  real K[S];
  real<lower=0.1> iqs[S];
  real<lower=0> isigma2;
  // real<lower=0> sigma_proc;
  real PE[S, N-1];  
}

transformed parameters {
  matrix[S, N] logB;  
  matrix[S, N] B;  
  real logqs[S];
  real sigma_proc = 1/sqrt(isigma2);
	for (i in 1:S){
		logqs[i] = log(1/iqs[i]);
		logB[i,1] =  log(K[i]); //log(PP_init[i]) + log(K[i]);
		B[i,1] = exp(logB[i,1]);
	}

	// when we have the disaggregated data
	for (t in 1:(N-1))
	{    
		for (i in 1:S)
		{
      logB[i,t+1] = log( fmax(B[i,t]+(r/p)*B[i,t]*(1-pow(B[i,t]/K[i],p))-C_obs[i,t],1) *exp(PE[i,t]) );
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
		r ~ lognormal(log(0.1), 0.2);
		iqs[i] ~ uniform(1000,100000);
		PP_init[i] ~ normal(0.9,0.1);
  }

// Likelihoods process error
	for (t in 1:(N-1))
	{   
		for (i in 1:S)
		{
		  PE[i,t] ~ normal(-0.5*sigma_proc*sigma_proc, sigma_proc); 
	  }
	}
	
// Index likelihoods: 
	int pos_i;
  pos_i = 1;
  real I_sig[N_Iobs];
  matrix[S, N] logI;
  real qs[S];

  for (s in 1:S) {
  	qs[s] = exp(logqs[s]);
    for (t in 1:N) {
      logI[s,t] = log(qs[s] * B[s,t]);
    }
  }
  for (i in 1:N_Iobs) {
    I_sig[i] = log(I_obs[i] + I_obs[i] * I_cv[i]) - log(I_obs[i]);
  }
  for (s in 1:S){
    segment(I_obs, pos_i, S_Iobs[s]) ~ lognormal(logI[s,segment(I_pos, pos_i, S_Iobs[s])], segment(I_sig, pos_i, S_Iobs[s]));
    pos_i = pos_i + S_Iobs[s];
  }	

// Biomass likelihoods: 
	int pos_b;
  pos_b = 1;
  real B_sig[N_Bobs];

  for (i in 1:N_Bobs) {
    B_sig[i] = log(B_obs[i] + B_obs[i] * B_cv[i]) - log(B_obs[i]);
  }
  for (s in 1:S){
    segment(B_obs, pos_b, S_Bobs[s]) ~ lognormal(logB[s,segment(B_pos, pos_b, S_Bobs[s])], segment(B_sig, pos_b, S_Bobs[s]));
    pos_b = pos_b + S_Bobs[s];
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
   MSY[s] = r * K[s] / (pow(p + 1,(p + 1) / p));
   Bmsy[s] = 0.4 * K[s];
   Fmsy[s] = MSY[s] / Bmsy[s];
   Stock_status[s] = B[s,N] / K[s];
  }
}
