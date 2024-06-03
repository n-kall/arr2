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

  // AR(12) residual variances of y and X
  real<lower=0> var_y_minn;
  vector<lower=0>[K/p] var_x_minn;

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
  // Regression Weights
  vector[K] z_beta; 
  real alpha;
  // Hyperparameters for regression weights
  real<lower=0> kappa2_k;
}

transformed parameters{
  vector[T] tau;
  vector<lower=0>[(K)] beta_sd;
  vector[K] beta;
  
  for (i in 1:K/p){
  for (j in 1:p){
    beta_sd[(j + (i-1)*p)] = sqrt(var_y_minn/var_x_minn[i]*kappa2_k/(j^2));
  }
  }
  
  beta = z_beta .* beta_sd;
  
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
  kappa2_k ~ gamma(1,1/(0.04^2));
  z_beta ~ std_normal();
  // likelihood
  z_tau ~ std_normal();
  alpha ~ normal(0, alpha_sd); // Intercept

  if (!prior_only) {
    yc ~ normal_id_glm(X, tau + alpha, beta, sigma_y);
  }

}

generated quantities{
 vector[T] y_pred;
 vector[T] mu_pred;
 real Intercept;
 real<lower=0, upper=1> R2_pp;

  Intercept = Ymean+alpha;
  mu_pred = Intercept +  X*beta + tau;
  R2_pp = variance(mu_pred) / (variance(mu_pred) + sigma_y^2);
    
}
