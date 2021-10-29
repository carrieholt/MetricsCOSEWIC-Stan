//Linear regression

data {
int<lower=0> N; // number of data items
vector[N] logAbd; // Y value
vector[N] Year;   // X value
real intercept_mean;
real intercept_sig;
real Sig_Gam_Dist;
real slope_mean;
real slope_sig;
}

parameters {
real intercept; // intercept
real slope; // coefficients for predictors
real<lower=0> sigma; // error scale
}



model {
vector[N] logAbd_Fit;
logAbd_Fit = intercept + slope * Year;
logAbd ~ normal(logAbd_Fit, sigma); // likelihood

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
vector[N] logAbd_Pred;
vector[N] logAbd_Fit_out;
for (i in 1:N){
 logAbd_Pred[i] = normal_rng(intercept + slope * Year[i], sigma);
 logAbd_Fit_out[i] = intercept + slope * Year[i];
}
}

