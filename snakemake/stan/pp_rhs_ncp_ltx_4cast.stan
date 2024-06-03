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
  vector[T] y; // observations
  int<lower=0> p; // AR order of covariates

  int<lower = 0> K; // number of covariates
  matrix[T,K] X; //covariate matrix
  
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

  int prior_only;
}

transformed data{
  real Ymean;
  vector[T] yc;
  matrix[K,K] xtx = X'*X;
  
  Ymean= mean(y);  
  for (i in 1:T){
    yc[i] = y[i] - Ymean;
  }
  
}

parameters {
  // States
  vector[T] z_tau;
  real z_tau1;
  // Observation variance
  real<lower=0> sigma_y;
  // State Variances
  real<lower=0> sigma_tau;
  real<lower=0> sigma_tau1;
  // State AR Weights
  real<lower = -1, upper = 1> phi;
  // local parameters for the horseshoe prior
  vector[(K)] phi_z;
  vector<lower=0>[(K)] hs_local;
  // horseshoe shrinkage parameters
  real<lower=0> hs_global; // global shrinkage parameter
  real<lower=0> hs_slab; // slab regularization parameter
  
  real alpha;
}

transformed parameters {
  vector[(K)] beta = horseshoe(phi_z, hs_local, hs_global, hs_scale_slab ^ 2 * hs_slab);
  vector[T] tau;
  
  tau[1] = z_tau[1] * sigma_tau1;
  for ( t in 2:T){
    tau[t] = phi * tau[(t-1)] + z_tau[t] * sigma_tau; 
  }
  
}

model {

  // priors
  sigma_y ~ student_t(3,0,sd(yc));
  sigma_tau ~ normal(0,sigma_tau_sd);
  sigma_tau1 ~ normal(0,sigma_tau1_sd);
  phi_z ~ std_normal();
  hs_global ~ student_t(hs_df_global, 0, p0 / (p + K - p0) * sigma_y / sqrt(T));
  hs_slab ~ inv_gamma(0.5 * hs_df_slab, 0.5 * hs_df_slab);
  hs_local ~ student_t(hs_df, 0, 1);

  // likelihood
  z_tau ~ std_normal();
  alpha ~ normal(0, alpha_sd); // Intercept

  if (!prior_only) {
    yc ~ normal_id_glm(X,tau+alpha,beta, sigma_y);
  }
}

generated quantities{
  vector[T] y_pred;
  vector[T] mu_pred;
  real Intercept;
  real<lower=0, upper=1> R2_pp;

  Intercept = Ymean + alpha;
  mu_pred = Intercept + X*beta + tau;
  R2_pp = variance(mu_pred) / (variance(mu_pred) + sigma_y^2 );
}


