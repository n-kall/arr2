data {

  int<lower=1> T; // number of time points
  int<lower=1> L; // last time point for training
  vector[T] Y; // observations
  int<lower=0> oos; // out-of-sample data index
  int<lower=0> p; // AR order

  // data for the R2D2 prior
  real<lower=0> mean_R2;  // mean of the R2 prior
  real<lower=0> prec_R2;  // precision of the R2 prior

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

  // Variance estimate of y
  real<lower=0> var_y;
  var_y = variance(Y);

  // concentration vector of the Dirichlet prior
  vector<lower=0>[p] cons;
  for (i in 1:p) {
    cons[i] = 0.1;
  }
}

parameters {
  real alpha; // intercept
  real<lower=0, upper=1> R2;
  real<lower=0> sigma;
  vector[p] phi_z;
  vector<lower=0>[p] gamma;
}

transformed parameters {
  real<lower=0> tau2 = R2 / (1 - R2);
  vector<lower=0,upper=1>[p] psi = gamma / sum(gamma);
  vector[p] phi_sd = sqrt(sigma^2 * tau2 * psi);
  vector[p] phi = phi_z .* sqrt(sigma^2/var_y * tau2 * psi);
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
