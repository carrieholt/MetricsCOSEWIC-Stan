// Linear regression. Missing values estimated as parameters
// Adapted from Stan code (manual section 3.3)
// https://mc-stan.org/docs/2_21/stan-users-guide/sliced-missing-data.html

data {
int<lower=0> N; // number of data items
real logAbd[N]; // Y value
//vector[N] logAbd_miss; // Y value
int<lower=0, upper=1> logAbd_obs[N];  // missingness indicator for Y 
vector[N] Year;   // X value
real intercept_mean;
real intercept_sig;
//real Sig_Gam_Dist;
real exp_rate;
real slope_mean;
real slope_sig;
}

transformed data {
  int n_obs = sum(logAbd_obs);  // number of observed cases
  int ns[n_obs];  // indices of observed cases
  int nmiss = N - n_obs; // number of missing cases
  int ns_miss[nmiss];  // indices of missing cases
  int ny = 1;
  int ny_miss = 1;
  for (n in 1:N) {
    if (logAbd_obs[n]) {
      ns[ny] = n;
      ny += 1;
    }
    else {
      ns_miss[ny_miss] = n;
      ny_miss += 1;
    }
  }
}




parameters {
real intercept; // intercept
real slope; // coefficients for predictors
real<lower=0> sigma; // error scale
real logAbd_miss[nmiss]; //estimate missing values
}

transformed parameters {
  real logAbd_sim[N];  // simulated data based on model
  logAbd_sim[ns] = logAbd[ns];  // where data are observed
  logAbd_sim[ns_miss] = logAbd_miss; // where data are missing
}
model {
//logAbd_sim ~ normal(intercept + slope * Year, sigma); // likelihood
for(i in 1:N){
  target += normal_lpdf(logAbd_sim[i] | intercept + slope * Year[i], sigma);
}

// PRIORS
// normal prior on intercept
//intercept ~ normal( intercept_mean, intercept_sig);
target += normal_lpdf(intercept | intercept_mean, intercept_sig);
// exponential on variance (sigma^2), need jacobian adjustment
//sigma ~ exponential(exp_rate);
target += exponential_lpdf( sigma | exp_rate);

// inverse gamma on variance (sigma^2), need jacobian adjustment
// pow(sigma, 2) ~ inv_gamma(Sig_Gam_Dist, Sig_Gam_Dist);
// Jacobian adjustment needed here
//target += log(2*sigma); // log|x^2 d/dx|

// normal prior on slope
//slope ~ normal( slope_mean, slope_sig);
target += normal_lpdf( slope | slope_mean, slope_sig);

}



generated quantities {
  vector[N] logAbd_Fit_out = intercept + Year * slope;  //Fitted values
  real logAbd_Pred[N];  // predictions based on model with variability
  for (i in 1:N) {
    logAbd_Pred[i] = normal_rng(logAbd_Fit_out[i], sigma);
  }
  
  real intercept_pp = normal_rng(intercept_mean, intercept_sig);
  real slope_pp = normal_rng(slope_mean, slope_sig);
  real sigma_pp = exponential_rng(exp_rate);
  array[N] real logAbd_PriorPred = normal_rng(intercept_pp + slope_pp * Year, sigma_pp);
  //vector[N] logAbd_PriorPred;
  // for (i in 1:N) {
  //  logAbd_PriorPred[i] = normal_rng(intercept_pp + slope_pp * Year);
  //}

}





