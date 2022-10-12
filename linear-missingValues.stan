// Linear regression that allows missing values, 
// Estimated in "generated quantities" section
// adapted from https://psyc-bayes-notes.netlify.app/missing-data.html#bayesian-approaches-for-missing-data

data {
int<lower=0> N; // number of data items
vector[N] logAbd; // Y value
int<lower=0, upper=1> logAbd_obs[N];  // missingness indicator for Y 
vector[N] Year;   // X value
real intercept_mean;
real intercept_sig;
real Sig_Gam_Dist;
real slope_mean;
real slope_sig;
}

transformed data {
  int n_obs = sum(logAbd_obs);  // number of observed cases
  int ns[n_obs];  // indices of observed cases
  int ny = 1; 
  for (n in 1:N) {
    if (logAbd_obs[n]) {
      ns[ny] = n;
      ny += 1;
    }
  }
}


parameters {
real intercept; // intercept
real slope; // coefficients for predictors
real<lower=0> sigma; // error scale
}

model {
logAbd[ns] ~ normal(intercept + slope * Year[ns], sigma); // likelihood

//Priors
// normal prior on intercept
intercept ~ normal( intercept_mean, intercept_sig);
// inverse gamma on variance (sigma^2), need jacobian adjustment
pow(sigma, 2) ~ inv_gamma(Sig_Gam_Dist, Sig_Gam_Dist);
// Jacobian adjustment
target += log(2*sigma); // log|x^2 d/dx|
// normal prior on slope
slope ~ normal( slope_mean, slope_sig);
}



generated quantities {
  real logAbd_Pred[N];  // simulated data based on model
  vector[N] logAbd_Fit_out = intercept + Year * slope;  
  for (i in 1:N) {
    logAbd_Pred[i] = normal_rng(logAbd_Fit_out[i], sigma);
  }
}




