functions {
  /* Efficient computation of the horseshoe prior
   * see Appendix C.1 in https://projecteuclid.org/euclid.ejs/1513306866
   * Args:
   *   z: standardized population-level coefficients
   *   lambda: local shrinkage parameters
   *   tau: global shrinkage parameter
   *   c2: slab regularization parameter
   * Returns:
   *   population-level coefficients following the horseshoe prior
   */
  vector horseshoe(vector z, vector lambda, real tau, real c2) {
    int K = rows(z);
    vector[K] lambda2 = square(lambda);
    vector[K] lambda_tilde = sqrt(c2 * lambda2 ./ (c2 + tau ^ 2 * lambda2));
    return z .* lambda_tilde * tau;
  }
}

data {
  int<lower=1> T; // number of time points
  int<lower=1> L; // last time point for training
  vector[T] Y; // observations
  int<lower=0> oos; // out-of-sample data index
  int<lower=0> p; // AR order
  int<lower = 0> K; // number of covariates
  matrix[T,K] X; //covariate matrix
  
  // sigma prior sd
  real<lower=0> sigma_sd; // sd of sigma prior

  // data for the horseshoe prior
  real<lower=0> hs_df; // local degrees of freedom
  real<lower=0> hs_df_global; // global degrees of freedom
  real<lower=0> hs_df_slab; // slab degrees of freedom
  real<lower=0> hs_scale_slab; // slab prior scale
  real<lower=0> p0;
  
  // intercept prior
  real alpha_mean;
  real<lower=0> alpha_sd;
}

transformed data {
  int<lower=1> par_tot = p+K;
  matrix[(T-p), p] Y_matrix;
  matrix[(T-p), K] X_matrix = X[(p+1):T,];
  matrix[(T-p), K + p] XY_matrix;
  for(t in 1:(T-p)) {
    for(i in 1:p) {
      Y_matrix[t, p-i+1] = Y[t + (i-1)];
    }
    XY_matrix[t] = append_col(X_matrix[t], Y_matrix[t]);
  }
}

parameters {
  real alpha; // intercept
  real<lower=0> sigma; // noise
  
  // local parameters for the horseshoe prior
  vector[(par_tot)] phi_z;
  vector<lower=0>[(par_tot)] hs_local;
  // horseshoe shrinkage parameters
  real<lower=0> hs_global; // global shrinkage parameter
  real<lower=0> hs_slab; // slab regularization parameter
}

transformed parameters {
  vector[(par_tot)] phi_beta = horseshoe(phi_z, hs_local, hs_global, hs_scale_slab ^ 2 * hs_slab);
  vector[p] phi = phi_beta[1:p];
  vector[K] beta = phi_beta[(p+1):(p+K)];
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
