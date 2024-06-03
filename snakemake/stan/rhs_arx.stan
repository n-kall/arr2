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
  vector[T] Y; // observations
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

model {

  // priors
  alpha ~ normal(alpha_mean, alpha_sd);
  phi_z ~ std_normal();
  hs_global ~ student_t(hs_df_global, 0, p0 / (p + K - p0) * sigma / sqrt(T));
  hs_slab ~ inv_gamma(0.5 * hs_df_slab, 0.5 * hs_df_slab);
  sigma ~ normal(0, sigma_sd);
  hs_local ~ student_t(hs_df, 0, 1);
  
  // likelihood
  Y[(p+1):T] ~ normal_id_glm(XY_matrix, alpha, append_row(beta, phi), sigma);
}


generated quantities{
  vector[(T-p)] mu = alpha + XY_matrix * append_row(beta, phi);
  real<lower=0, upper=1> postR2 = variance(mu) / (variance(mu) + sigma^2 );
  //  real<lower=0, upper=1> postR2_lags = (phi' * Y_matrix' * Y_matrix * phi) /
  // (append_row(beta, phi)' * XY_matrix' * XY_matrix * append_row(beta, phi));
  //  real<lower=0, upper=1> postR2_reg = (beta' * X_matrix' * X_matrix * beta) /
  //    (append_row(beta, phi)' * XY_matrix' * XY_matrix * append_row(beta, phi));
}
  
