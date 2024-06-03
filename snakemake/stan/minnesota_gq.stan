data {

  int<lower=1> T; // number of time points
  int<lower=1> L; // last time point for training
  vector[T] Y; // observations
  int<lower=0> oos; // out-of-sample data index
  int<lower=0> p; // AR order

  // sigma prior sd
  real<lower=0> sigma_sd; // sd of sigma prior
  
  // intercept prior
  real alpha_mean;
  real<lower=0> alpha_sd;
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

  real alpha;
  real<lower=0> kappa2;
  vector[p] phi;
  real<lower=0> sigma;

}

transformed parameters {
  vector[p] phi_sd;
  for (i in 1:p) {
    phi_sd[i] = sqrt(kappa2 / i^2);
  }
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
