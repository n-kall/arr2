data {

  int<lower=1> T; // number of time points
  vector[T] Y; // observations
  int<lower=0> p; // AR order of covariates
  int<lower=1> L; // last time point for training
  int<lower=1> Tc; // train data length
  int<lower=0> oos; // out-of-sample data index
  int<lower = 0> K; // number of covariates
  matrix[T,K] X; //covariate matrix

  // AR(1) variance parameter for the states
  real<lower=0> phi_sd;
  
  // Variance of state standard deviations
  real<lower=0> sigma_tau_sd;
  real<lower=0> sigma_tau1_sd; //for first state

  // intercept prior
  real<lower=0> alpha_sd;

  // AR(12) residual variances of y and X
  real<lower=0> var_y_minn;
  vector<lower=0>[K/p] var_x_minn;

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
  real<lower = -1, upper = 1> phi;
  // Regression Weights
  vector[K] z_beta; 
  real alpha;
  // Hyperparameters for regression weights
  real<lower=0> kappa2_k;
}

transformed parameters{
  vector[Tc] tau;
  vector<lower=0>[(K)] beta_sd;
  vector[K] beta;
  
  for (i in 1:K/p){
    for (j in 1:p){
      beta_sd[(j + (i-1)*p)] = sqrt(var_y_minn/var_x_minn[i]*kappa2_k/(j^2));
    }
  }
  
  beta = z_beta .* beta_sd;
  
  tau[1] = z_tau[1] * sigma_tau1;
  
  for ( t in 2:Tc){
    tau[t] = phi * tau[(t-1)] + z_tau[t] * sigma_tau; 
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
