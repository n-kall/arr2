data {

  int<lower=1> T; // number of time points
  int<lower=1> L; // last time point for training
  int<lower=0> oos; // out-of-sample data index
  vector[T] Y; // observations
  int<lower=0> p; // AR order

  int<lower = 0> K; // number of covariates
  matrix[T,K] X; //covariate matrix

  // sigma prior sd
  real<lower=0> sigma_sd; // sd of sigma prior

  // intercept prior
  real alpha_mean;
  real<lower=0> alpha_sd;

  // AR(4) residual variances of y and X
  real<lower=0> var_y_minn;
  vector<lower=0>[K] var_x_minn;

}

transformed data {
  matrix[(T-p), p] Y_matrix;
  matrix[(T-p),K] X_matrix = X[(p+1):T,];
  matrix[(T-p),K + p] XY_matrix;
  for(t in 1:(T-p)) {
    for(i in 1:p) {
      Y_matrix[t, p-i+1] = Y[t + (i-1)];
    }
    XY_matrix[t] = append_col(X_matrix[t],Y_matrix[t]);
  }
}

parameters {

  real alpha;
  real<lower=0> kappa2_lags;
  vector[p] phi;
  real<lower=0> sigma;
  real<lower=0> kappa2_k;
  vector[K] beta;

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
