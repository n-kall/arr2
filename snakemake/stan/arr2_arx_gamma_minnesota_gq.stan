data {

  int<lower=1> T; // number of time points
  int<lower=1> L; // last time point for training
  vector[T] Y; // observations
  int<lower=0> oos; // out-of-sample data index
  int<lower=0> p; // AR order

  int<lower = 0> K; // number of covariates
  matrix[T,K] X; //covariate matrix

  // data for the R2D2 prior
  real<lower=0> mean_R2;  // mean of the R2 prior
  real<lower=0> prec_R2;  // precision of the R2 prior

  // sigma prior sd
  real<lower=0> sigma_sd; // sd of sigma prior

  // intercept prior
  real alpha_mean;
  real<lower=0> alpha_sd;

  // Variance estimates of y and X
  real<lower=0> var_y;
  vector<lower=0>[K] var_x;

}

transformed data {
  matrix[(T-p), p] Y_matrix;
  matrix[(T-p), K] X_matrix = X[(p+1):T,];
  matrix[(T-p), (K + p)] XY_matrix;
  for(t in 1:(T-p)) {
    for(i in 1:p) {
      Y_matrix[t, p-i+1] = Y[t + (i-1)];
    }
    XY_matrix[t] = append_col(X_matrix[t],Y_matrix[t]);
  }
  // concentration vector of the Dirichlet prior
  vector<lower=0>[(p+K)] cons;
  for (i in 1:p) {
    cons[i] = (p^2)/10 * 1/(i^2);
  }
  for (i in (p + 1):(p + K)) {
    cons[i] = 0.1;
  }

}

parameters {
  real alpha; // intercept
  vector<lower=0>[(p + K)] gamma;
  real<lower=0, upper=1> R2;
  real<lower=0> sigma;
  vector[p] phi_z;
  vector[K] beta_z;
}

transformed parameters {
  real<lower=0> tau2 = R2 / (1 - R2);
  vector<lower=0,upper=1>[(p+K)] psi = gamma / sum(gamma);
  vector[p] phi = phi_z .* sqrt(sigma^2/var_y * tau2 * psi[1:p]);
  vector[K] beta = beta_z .* sqrt(sigma^2 * tau2 * psi[(p+1):(p+K)] ./ var_x);
}

generated quantities {
  // log-likelihood
  real log_lik;
  // posterior predictive
  real pred;

  if (oos != 0) {
    log_lik = normal_id_glm_lpdf(Y[oos] | to_matrix(XY_matrix[oos-p]), alpha, append_row(beta, phi), sigma);
  } else {
    log_lik = normal_id_glm_lpdf(Y[(L+1):T] | to_matrix(XY_matrix[(L+1-p):(T-p)]), alpha, append_row(beta, phi), sigma);
  }

  if (oos != 0) {
    pred = normal_rng(alpha + XY_matrix[oos-p] * append_row(beta, phi), sigma);
  } else {
    pred = 0;
  }
}
