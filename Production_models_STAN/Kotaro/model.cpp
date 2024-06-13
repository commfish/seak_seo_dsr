// model need to have some hierarchy to help constrain the model 


// State-space SP
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
// data:
DATA_INTEGER(Npre);
DATA_INTEGER(S);
DATA_SCALAR(p);
DATA_MATRIX(C_obs);
DATA_MATRIX(I_obs);
DATA_SCALAR(ratio);

// parameters:
PARAMETER_VECTOR(logr); // population growth rate parameter
PARAMETER_VECTOR(PP_init_logit); // density dependence parameter
PARAMETER_VECTOR(logK); // density dependence parameter
PARAMETER_VECTOR(log_invq); // density dependence parameter
PARAMETER(log_sigma_proc); // log(process SD)
PARAMETER_MATRIX(PE); // unobserved state vector

// procedures: (transformed parameters)
Type sigma_proc = exp(log_sigma_proc);
Type sigma_obs = sigma_proc*sqrt(ratio);
vector<Type> r = exp(logr);
vector<Type> PP_init = invlogit(PP_init_logit);
vector<Type> K = exp(logK);
vector<Type> q = 1/exp(log_invq);

// data transformation
Type n = p + 1; 
Type bb = pow(n, n/(n-1));

// reports on transformed parameters:
Type nll = 0.0; // initialize negative log likelihood

// process model:
matrix<Type> B_temp(S,Npre); 
matrix<Type> logB(S,Npre); 
matrix<Type> B(S,Npre); 

	for (int i=0; i<S; i++ ){
		logB(i,0) =  log(PP_init(i)) + logK(i);
		B(i,0) = exp(logB(i,0));
	}
		
	for (int t=0; t < (Npre-1); t++)
	{    
		for (int i=0; i<S; i++ )
		{
      B_temp(i,t+1) = B(i,t) + (r(i)/p)*B(i,t)*(1-pow(B(i,t)/K(i),p)) - C_obs(i,t);
			if (B_temp(i,t+1) < 0 ) B_temp(i,t+1) = 1; 
			logB(i,t+1) = log(B_temp(i,t+1)) *exp(PE(i,t)) ;
			B(i,t+1) = exp(logB(i,t+1));
		}
	}
	
	for (int t=0; t<(Npre-1); t++)
	{   
		for (int i=0; i<S; i++)
		{
			nll -= dnorm(PE(i,t), -0.5*sigma_proc*sigma_proc, sigma_proc, true);
	  }
	}

// observation model:
	for (int t=0; t<Npre; t++)
	{   
		for (int i=0; i<S; i++)
		{
			nll -= dnorm(log(I_obs(i,t)), log(q(i)) + logB(i,t) - 0.5*sigma_obs*sigma_obs, sigma_obs, true) - log(I_obs(i,t));
	  }
	}
	
// penalty
  for (int i=0; i<S; i++)
  {
    nll -= dnorm(logr(i), Type(-2), Type(0.5), true); 
    nll -= dnorm(logK(i), Type(12), Type(1), true); 
    // vector<Type> vals = I_obs.row(i);
    // vector<Type> logvals = log(vals);
    // nll -= dnorm(log(q(i)), (sum(logvals)-sum(vector<Type>(logB.row(i))))/Type(Npre), Type(0.1), true); 
  }

ADREPORT(r);
ADREPORT(K);
ADREPORT(q);
REPORT(logB);
ADREPORT(PP_init);
ADREPORT(sigma_proc);

return nll;
}
