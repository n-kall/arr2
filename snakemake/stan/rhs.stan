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
  matrix[(T-p), p] Y_matrix;
  matrix[p,p] YY;
  for(t in 1:(T-p)) {
    for(i in 1:p) {
      Y_matrix[t, p-i+1] = Y[t + (i-1)];
    }
  }
  
  YY = Y_matrix'*Y_matrix;
}


parameters {
  real alpha; // intercept
  real<lower=0> sigma; // noise
  // local parameters for the horseshoe prior
  vector[p] phi_z;
  vector<lower=0>[p] hs_local;
  // horseshoe shrinkage parameters
  real<lower=0> hs_global; // global shrinkage parameter
  real<lower=0> hs_slab; // slab regularization parameter
}

transformed parameters {
  vector[p] phi = horseshoe(phi_z, hs_local, hs_global, hs_scale_slab ^ 2 * hs_slab);
}

model {

  // priors
  alpha ~ normal(alpha_mean, alpha_sd);
  phi_z ~ std_normal();
  hs_global ~ student_t(hs_df_global, 0, p0 / (p - p0) * sigma / sqrt(T-p));
  hs_slab ~ inv_gamma(0.5 * hs_df_slab, 0.5 * hs_df_slab);
  sigma ~ normal(0, sigma_sd);
  hs_local ~ student_t(hs_df, 0, 1);
  
  // likelihood
  Y[(p+1):T] ~ normal_id_glm(Y_matrix, alpha, phi, sigma);
}


generated quantities{
  vector[(T-p)] mu = alpha +  Y_matrix*phi;
  real<lower=0, upper=1> postR2 = variance(mu) / (variance(mu) + sigma^2 );
  vector[p] rel_R2;
  
  for (i in 1:p){
    rel_R2[i] =  phi[i]^2 * YY[i,i] / (phi'*YY*phi + (T-p) * sigma^2);
  }
  }
