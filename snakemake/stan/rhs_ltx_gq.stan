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

  int<lower=0> K; // number of covariates
  matrix[T, K] X; //covariate matrix

  // AR(1) variance parameter for the states
  real<lower=0> phi_sd;

  // Variance of state standard deviations
  real<lower=0> sigma_tau_sd;
  real<lower=0> sigma_tau1_sd; //for first state

  // intercept prior
  real<lower=0> alpha_sd;

  // data for the horseshoe prior
  real<lower=0> hs_df; // local degrees of freedom
  real<lower=0> hs_df_global; // global degrees of freedom
  real<lower=0> hs_df_slab; // slab degrees of freedom
  real<lower=0> hs_scale_slab; // slab prior scale
  real<lower=0> p0;
  int<lower=1> L; // last time point for training
  int<lower=1> Tc; // train data length
  int<lower=0> oos; // out-of-sample data index

}

transformed data{
  real Ymean;
  vector[Tc] Yc;
  matrix[K, K] xtx;
  real Yc_sd;
  xtx = X' * X;

  Ymean = mean(Y[1:Tc]);
  for (i in 1:Tc){
    Yc[i] = Y[i] - Ymean;
  }

  Yc_sd = sd(Yc);
}

parameters {
  // States
  vector[Tc] z_tau;
  real z_tau1;
  // Observation variance
  real<lower=0> sigma;
  // State Variances
  real<lower=0> sigma_tau;
  real<lower=0> sigma_tau1;
  // State AR Weights
  real<lower=-1, upper=1> phi;
  // local parameters for the horseshoe prior
  vector[(K)] phi_z;
  vector<lower=0>[K] hs_local;
  // horseshoe shrinkage parameters
  real<lower=0> hs_global; // global shrinkage parameter
  real<lower=0> hs_slab; // slab regularization parameter

  real alpha;
}

transformed parameters {
  vector[K] beta = horseshoe(phi_z, hs_local, hs_global, hs_scale_slab^2 * hs_slab);
  vector[Tc] tau;

  tau[1] = z_tau[1] * sigma_tau1;
  for ( t in 2:Tc){
    tau[t] = phi * tau[t-1] + z_tau[t] * sigma_tau;
  }

}

generated quantities{
  real log_lik;
  real new_tau;

  new_tau = phi * tau[Tc] + z_tau[Tc] * sigma_tau;

  if (oos != 0) {
    log_lik = normal_id_glm_lpdf(Y[oos] | to_matrix(X[oos,]), Ymean + alpha + new_tau, beta, sigma);
  } else {
    log_lik = normal_id_glm_lpdf(Y[(L+1):T] | to_matrix(X[(L+1):T,]), Ymean + alpha + tau[(L+1):T], beta, sigma);
  }
}
