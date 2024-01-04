data {
  int<lower=0> N; // number of years 
  int<lower=0> S; // number of strata 
  real C_obs[N,S]; // estimated treatment effects
  real C_cv[N,S]; // estimated treatment effects
  int<lower=0> N_I1obs; //  total number of non-NA values
  real I1_obs[N_I1obs]; // estimated treatment effects
  real I1_cv[N_I1obs]; // estimated treatment effects
  int<lower=0> S_I1obs[S]; //number of observations per strata, corresponds to s on 133-134 Stan manual
  int<lower=0> I1_pos[N_I1obs]; // index positions as one vector but relevant to each strata
  int<lower=0> N_I2obs; //  total number of non-NA values
  real I2_obs[N_I2obs]; // estimated treatment effects
  real I2_cv[N_I2obs]; // estimated treatment effects
  int<lower=0> S_I2obs[S]; //number of observations per strata, corresponds to s on 133-134 Stan manual
  int<lower=0> I2_pos[N_I2obs]; // index positions as one vector but relevant to each strata
  int<lower=0> N_Bobs; //  total number of non-NA values
  real B_obs[N_Bobs];
  real B_cv[N_Bobs];
  int<lower=0> S_Bobs[S]; //number of observations per strata, corresponds to s on 133-134 Stan manual
  int<lower=0> B_pos[N_Bobs]; // index positions as one vector but relevant to each strata
  real pte; //pella tomlinson exponent
  real msy_frac; //where Bmsy occurs
  vector<lower=0>[S] alpha; // dirchlet pi prior
  real pe_bound;
  real isig2_mu; //mean for prior on isigma2
  real isig2_sig; //sigma for prior on isigma2
  real<lower=0,upper=1.1> P1_mu;  //mean for prior on P1
  real P1_sig; //sigma for prior on P1
  int pe_switch; //switch to turn process error estimation on and off; 0 = off, 1 = on
  int bio_switch; // switch for fitting to a biomass estimate
  int ind1_switch; // switch for fitting to an index
  int ind2_switch; // switch for fitting to an index
}
parameters {
  real<lower=0> r;
  real<lower=0> K;
  real<lower=0> isigma2;
  real logvar; //process error
//  real rho; // autocorrelation term
  real<lower=0> i1tau2[S];
  real<lower=0> i2tau2[S];
  simplex[S] pi;
  real eps[N];
  matrix[N,S] mu;
  matrix<lower=0>[N,S] C;  
}

transformed parameters {
 real sigma;
 real tau1[S]; 
 real tau2[S];
 sigma = 1.0/sqrt(isigma2);
 for (s in 1:S){
    tau1[s] = 1.0/sqrt(i1tau2[s]);
    tau2[s] = 1.0/sqrt(i2tau2[s]);
 } 
}

model {
  matrix[N,S] Pmed;
  matrix[N,S] P;
  matrix[N,S] I1med;
  matrix[N,S] I2med;
  matrix[N,S] C_sig;
  real B_sig[N_Bobs];
  real I1_sig[N_I1obs];
  real I2_sig[N_I2obs];
  matrix[N,S] Biomass;
  matrix[N,S] log_bio;
  matrix[N,S] logC;
  vector[N_I1obs] Z1;
  vector[N_I2obs] Z2;
  real q1hat[S];
  real q2hat[S];

// Priors:  
  r ~ beta(1,1) T[0.01,0.99]; //r ~ lognormal(-1.38,3.845) T[0.01,1.2];
  isigma2 ~ gamma(isig2_mu, isig2_sig) T[0.01,1000]; //isigma2 ~ gamma(3.75, 0.1); // T[0.01,1000];
  K ~ lognormal(10, 2) T[10,10000];; //lognormal(10, 2) T[10,10000]; // uniform(100, 10000);
  pi ~ dirichlet(alpha);
  for (s in 1:S){
    i1tau2[s] ~ gamma(1.709, 0.00861);
    i2tau2[s] ~ gamma(1.709, 0.00861);
    for (t in 1:N){
      mu[t,s] ~ normal(0,1);
    }
  } 

  //PROCESS ERROR
  real PE[N];
  real pe_sigma;
  real tau_log_pe;
  if (pe_switch == 0) {  // no process error
    pe_sigma = 0; //sqrt(exp(-12));
    tau_log_pe = 2000000000; //1 / log(pe_sigma * pe_sigma + 1);
    for (t in 1:N){
      eps[t] ~ normal(0,tau_log_pe) T[-0.01,0.01];
      PE[t] = 0; 
    }
  }
  else if (pe_switch == 1) {  // estimate process error
    logvar ~ uniform(-12,pe_bound);
    pe_sigma = sqrt(exp(logvar));
    tau_log_pe = 1 / log(pe_sigma * pe_sigma + 1);
    for (t in 1:N){
      eps[t] ~ normal(0,tau_log_pe) T[-0.1,0.1];
      PE[t] = eps[t] - (pe_sigma * pe_sigma / 2);
    }
  }

// Set initial state
  for (s in 1:S) {
//    P[1,s] ~ lognormal(P1_mu,P1_sig) T[0.05,1.01]; 
    P[1,s] = exp(sigma*mu[1,s]);
//    Pmed[1,s] = log(P[1,s]);
    Biomass[1,s] = P[1,s] * K * pi[s];
  }

// time steps of the model
for (s in 1:S){
  for (t in 2:N) {
    P[t,s] = fmax((P[t - 1,s] + (r / (pte - 1)) * P[t - 1,s] * (1 - pow(P[t - 1,s] , (pte - 1))) - C[t - 1,s] / (K * pi[s])) * exp(sigma*mu[t,s]), 0.001) * exp(PE[t-1]); //P[t,s] ~ lognormal(Pmed[t,s], sigma) T[0.05, 1.6];
    Biomass[t,s] = K * P[t,s] * pi[s];
  }
}

// Likelihoods
for (s in 1:S) {
  for (t in 1:N){
    logC[t,s] = log(C[t,s]);
    C_sig[t,s] = log(C_obs[t,s] + C_obs[t,s] * C_cv[t,s]) - log(C_obs[t,s]);   
    C_obs[t,s] ~ lognormal(logC[t,s],C_sig[t,s]);
    log_bio[t,s] = log(Biomass[t,s]);
  }
}
// biomass likelihoods
if (bio_switch == 0){  //if no biomass dont try to fit to the data
}
else if (bio_switch == 1){
  int pos;
  pos = 1;
  for (i in 1:N_Bobs) {
    B_sig[i] = log(B_obs[i] + B_obs[i] * B_cv[i]) - log(B_obs[i]);
  }
  for (s in 1:S){
    segment(B_obs, pos, S_Bobs[s]) ~ lognormal(log_bio[segment(B_pos, pos, S_Bobs[s]),s], segment(B_sig, pos, S_Bobs[s]));
    pos = pos + S_Bobs[s];
  }
}

// index1 likelihoods
if (ind1_switch == 0) { //if no index dont try to fit to the data
}
else if (ind1_switch == 1) {
  int pos2;
  pos2 = 1;
  int pos3;
  pos3 = 1;
  real Pmed1;
  real Kpi;

  int geez;
  geez = 0;
  int s_start;

  for (s in 1:S) {  //vectorize K * pi and P to align with observed values...
    s_start = pos3;
    for (j in 1:S_I1obs[s]){
      Kpi = K * pi[s];
      Pmed1 = P[I1_pos[j + geez] ,s];
      Z1[pos3] = log(I1_obs[pos3] / (Pmed1 * Kpi));
      pos3 = pos3 + 1;
    }
    geez = geez + S_I1obs[s];
    vector[S_I1obs[s]] Zs;
    Zs = segment(Z1, pos2, S_I1obs[s]);
    q1hat[s] = sum(Zs) / S_I1obs[s]; //sum(segment(Z1med, pos2, S_I1obs[s]));
    for (i in s_start:(pos3 - 1)){
      Z1[i] ~ normal(q1hat[s], tau1[s]);
    }
    pos2 = pos2 + S_I1obs[s];
  }
}

// index2 likelihoods
if (ind2_switch == 0) { //if no index dont try to fit to the data
}
else if (ind2_switch == 1) {
  int pos4;
  pos4 = 1;
  int pos5;
  pos5 = 1;
  real Pmed2;
  real Kpi;

  int geez;
  geez = 0;
  int s_start;

for (s in 1:S) {  //vectorize K * pi and P to align with observed values...
    s_start = pos4;
    for (j in 1:S_I2obs[s]){
      Kpi = K * pi[s];
      Pmed2 = P[I2_pos[j + geez] ,s];
      Z2[pos4] = log(I2_obs[pos4] / (Pmed2 * Kpi));
      pos4 = pos4 + 1;
    }
    geez = geez + S_I2obs[s];
    vector[S_I2obs[s]] Zs;
    Zs = segment(Z2, pos5, S_I2obs[s]);
    q2hat[s] = sum(Zs) / S_I2obs[s]; 
    for (i in s_start:(pos4 - 1)){
      Z2[i] ~ normal(q2hat[s], tau2[s]);
    }
    pos5 = pos5 + S_I2obs[s];
  }
}

}  //end model
   
generated quantities {
  matrix[N,S] P;
  matrix[N,S] I1med;
  matrix[N,S] I1est;
  matrix[N,S] I2med;
  matrix[N,S] I2est;
  matrix[N,S] log_bio; 
  matrix[N,S] Biomass;
  real log_bio_pp[N_Bobs];
  real Bio_new[N_Bobs];
  real B_sig[N_Bobs];
  real I1med_pp[N_I1obs];
  real I1_new[N_I1obs];
  real I1_sig[N_I1obs];
  vector[N_I1obs] Z1;
  real I2med_pp[N_I2obs];
  real I2_new[N_I2obs];
  real I2_sig[N_I2obs];
  vector[N_I2obs] Z2;
  real MSY;
//  real EMSY;
  real Bmsy;
  real Fmsy;
  real Stock_status; 
  real Tot_biomass[N];
  int j; 
  int pos;
  int pos1;
  int pos2;
  int pos3;
  real q1[S];
  real q2[S];

  //posterior predictions (hint, the parameterization of dlnorm is not the same as in R)
  for (s in 1:S) {
    P[1,s] = exp(sigma*mu[1,s]);
    for (t in 2:N) {
      P[t,s] = fmax((P[t - 1,s] + (r / (pte - 1)) * P[t - 1,s] * (1 - pow(P[t - 1,s] , (pte - 1))) - C[t - 1,s] / (K * pi[s])) * exp(sigma*mu[t,s]), 0.001);
    }
    for (t in 1:N) {
      Biomass[t,s] = K * P[t,s] * pi[s];
      log_bio[t,s] = log(Biomass[t,s]);
   }
  }
  //pp for biomass
  if (bio_switch == 0) {
    for (i in 1:N_Bobs) {
      B_sig[i] = 0;
      Bio_new[i] = 0;
    }
  }
  else if (bio_switch == 1) {
    j = 1;
    pos = 1;
    for (s in 1:S){    
      vector[S_Bobs[s]] ugh;
      ugh = log_bio[segment(B_pos, pos, S_Bobs[s]),s];
      for (ii in 1:rows(ugh)) {
        log_bio_pp[j] = ugh[ii];
        j = j + 1;
      }
      pos = pos + S_Bobs[s];
    }
    for (i in 1:N_Bobs) {
      B_sig[i] = log(B_obs[i] + B_obs[i] * B_cv[i]) - log(B_obs[i]);
      Bio_new[i] = lognormal_rng(log_bio_pp[i], B_sig[i]);
    }  
  }
 
  //pp for index 1
  if (ind1_switch == 0) {
    for (i in 1:N_I1obs) {
      I1_sig[i] = 0;
      I1_new[i] = 0;
      I1med_pp[i] = 0;
      Z1[i] = 0;
    }
    for (s in 1:S){
      for (t in 1:N) {
        I1med[t,s] = 0;
        I1est[t,s] = exp(I1med[t,s]);
      }
    }
  }
  else if (ind1_switch == 1) {
    pos = 1;
    pos1 = 1;
    real q1hat[S];
  //  real q1[S];
    real Pmed1;
    real Kpi;

    int geez;
    geez = 0;
    int s_start;
    int k;
    k = 1;

    for (s in 1:S) {  //vectorize K * pi and P to align with observed values...
      s_start = pos1;
      for (a in 1:S_I1obs[s]){
        Kpi = K * pi[s];
        Pmed1 = P[I1_pos[a + geez] ,s];
        Z1[pos1] = log(I1_obs[pos1] / (Pmed1 * Kpi));
        pos1 = pos1 + 1;
      }
      geez = geez + S_I1obs[s];
      vector[S_I1obs[s]] Zs;
      Zs = segment(Z1, pos, S_I1obs[s]);
      q1hat[s] = sum(Zs) / S_I1obs[s]; //sum(segment(Z1med, pos2, S_I1obs[s]));
      q1[s] = exp(q1hat[s]);
      for (t in 1:N) {
        I1med[t,s] = log(q1[s] * K * pi[s] * P[t,s]);
        I1est[t,s] = exp(I1med[t,s]);
      }
      vector[S_I1obs[s]] ugh;
      ugh = I1med[segment(I1_pos, pos, S_I1obs[s]),s];
      for (ii in 1:rows(ugh)) {
        I1med_pp[k] = ugh[ii];
        k = k + 1;
      }
      pos = pos + S_I1obs[s];
    }

    for (i in 1:N_I1obs) {
      I1_sig[i] = log(I1_obs[i] + I1_obs[i] * I1_cv[i]) - log(I1_obs[i]);
      I1_new[i] = lognormal_rng(I1med_pp[i], I1_sig[i]);
    }
  }

  //pp for index 2
  if (ind2_switch == 0) {
    for (i in 1:N_I2obs) {
      I2_sig[i] = 0;
      I2_new[i] = 0;
      I2med_pp[i] = 0;
      Z2[i] = 0;
    }
    for (s in 1:S) {
      for (t in 1:N) {
        I2med[t,s] = 0;
        I2est[t,s] = exp(I2med[t,s]);
      }
    }
    
  }
  else if (ind1_switch == 2) {
    pos = 1;
    pos1 = 1;
    real q2hat[S];
  //  real q2[S];
    real Pmed2;
    real Kpi;

    int geez;
    geez = 0;
    int s_start;
    int k;
    k = 1;

    for (s in 1:S) {  //vectorize K * pi and P to align with observed values...
      s_start = pos1;
      for (a in 1:S_I1obs[s]){
        Kpi = K * pi[s];
        Pmed2 = P[I2_pos[a + geez] ,s];
        Z2[pos1] = log(I2_obs[pos1] / (Pmed2 * Kpi));
        pos1 = pos1 + 1;
      }
      geez = geez + S_I1obs[s];
      vector[S_I2obs[s]] Zs;
      Zs = segment(Z2, pos, S_I2obs[s]);
      q2hat[s] = sum(Zs) / S_I2obs[s]; //sum(segment(Z2med, pos2, S_I2obs[s]));
      q2[s] = exp(q2hat[s]);
      for (t in 1:N) {
        I2med[t,s] = log(q2[s] * K * pi[s] * P[t,s]);
        I2est[t,s] = exp(I2med[t,s]);
      }
      vector[S_I2obs[s]] ugh;
      ugh = I2med[segment(I2_pos, pos, S_I2obs[s]),s];
      for (ii in 1:rows(ugh)) {
        I2med_pp[k] = ugh[ii];
        k = k + 1;
      }
      pos = pos + S_I2obs[s];
    }

    for (i in 1:N_I1obs) {
      I2_sig[i] = log(I2_obs[i] + I2_obs[i] * I2_cv[i]) - log(I2_obs[i]);
      I2_new[i] = lognormal_rng(I2med_pp[i], I2_sig[i]);
    }
  }
  
  // Other ouputs
   MSY = r * K / (pow(pte + 1,(pte + 1) / pte));
  // EMSY = r / (q[s] * (pte + 1)); //r/(2*q);
   Bmsy = msy_frac * K;
   Fmsy = MSY / Bmsy;
   Stock_status = sum(Biomass[N,]) / K;
   for (i in 1:N) {
    Tot_biomass[i] = sum(Biomass[i,]);
   }
  }
  
  



