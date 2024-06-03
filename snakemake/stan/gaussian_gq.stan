data {

  int<lower=1> T; // number of time points
  int<lower=1> L; // last time point for training
  vector[T] Y; // observations
  int<lower=0> oos; // out-of-sample data index
  int<lower=0> p; // AR order

  // sigma prior sd
  real<lower=0> sigma_sd; // sd of sigma prior

  // phi prior sd
  real<lower=0> phi_sd; // sd of phi prior

  // intercept prior
  real alpha_mean;
  real<lower=0> alpha_sd;

  int<lower=0, upper=1> alpha_cent;
  int<lower=0, upper=1> phi_cent;
  int<lower=0, upper=1> sigma_cent;
}

transformed data {
  matrix[(T-p), p] Y_matrix;
  for(t in 1:(T-p)) {
    for(i in 1:p) {
      Y_matrix[t, p-i+1] = Y[t + (i-1)];
    }
  }
}

parameters {

  vector<multiplier=(phi_cent < 1 ? phi_sd : 1.0)>[p] phi;
  real<offset=(alpha_cent < 1 ? alpha_mean : 0.0),
        multiplier=(alpha_cent < 1 ? alpha_sd : 1.0)> alpha; // intercept
  real<lower=0> sigma;

}

generated quantities {
  // log-likelihood
  real log_lik;
  if (oos != 0) {
    log_lik = normal_id_glm_lpdf(Y[oos] | to_matrix(Y_matrix[oos-p]), alpha, phi, sigma);
  } else {
    log_lik = normal_id_glm_lpdf(Y[(L+1):T] | Y_matrix[(L+1-p):(T-p),], alpha, phi, sigma);
  }
  // posterior predictive
  real pred;
  if (oos != 0) {
    pred = normal_rng(alpha + Y_matrix[oos-p] * phi, sigma);
  } else {
    pred = 0;
  }
}
