data {

  int<lower=1> T; // number of time points
  vector[T] Y; // observations
  int<lower=0> p; // AR order of covariates

  int<lower=0> K; // number of covariates
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
  vector<lower=0>[K / p] var_x_minn;

}

transformed data{
  real Ymean;
  vector[T] Yc;
  matrix[K, K] xtx = X'*X;
  real Yc_sd;

  
  Ymean = mean(Y);  
  for (i in 1:T){
    Yc[i] = Y[i] - Ymean;
  }

  Yc_sd = sd(Yc);
}


parameters {
  // States
  vector[T] z_tau;
  real z_tau1;
  // Observation variance
  real<lower=0> sigma;
  // State Variances
  real<lower=0> sigma_tau;
  real<lower=0> sigma_tau1;
  // State AR Weights
  real<lower=-1, upper=1> phi;
  // Regression Weights
  vector[K] z_beta; 
  real alpha;
  // Hyperparameters for regression weights
  real<lower=0> kappa2_k;
}

transformed parameters{
  vector[T] tau;
  vector<lower=0>[K] beta_sd;
  vector[K] beta;
  
  for (i in 1:(K / p)) {
    for (j in 1:p){
      beta_sd[(j + (i - 1) * p)] = sqrt(var_y_minn / var_x_minn[i] * kappa2_k / (j^2));
    }
  }
  
  beta = z_beta .* beta_sd;
  
  tau[1] = z_tau[1] * sigma_tau1;
  
  for ( t in 2:T){
    tau[t] = phi * tau[(t - 1)] + z_tau[t] * sigma_tau; 
  }
}

model {

  // priors
  sigma ~ student_t(3, 0, Yc_sd);
  sigma_tau ~ normal(0, sigma_tau_sd);
  sigma_tau1 ~ normal(0, sigma_tau1_sd);
  kappa2_k ~ gamma(1, 1 / (0.04^2));
  z_beta ~ std_normal();
  phi ~ normal(0, phi_sd);
  z_tau ~ std_normal();
  alpha ~ normal(0, alpha_sd); // Intercept
  
  // likelihood
  Yc ~ normal_id_glm(X, tau + alpha, beta, sigma);

}

generated quantities{
  vector[T] mu_pred;
  real Intercept;
  real postR2;
  real postR2_trend;
  real postR2_reg;
 
  Intercept = Ymean+alpha;
  mu_pred = Intercept +  X*beta + tau;
  postR2 = variance(mu_pred) / (variance(mu_pred) + sigma^2 ); 
 
  postR2_trend = (tau' * tau) / (Intercept^2 + tau' * tau + beta' * xtx * beta + T * sigma^2);
  postR2_reg = (beta' * xtx*beta) / (Intercept^2 + tau' * tau + beta' * xtx * beta + T * sigma^2);
 

    
}

