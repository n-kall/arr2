data {
  int<lower=1> T;
  vector[T] Y;
  int<lower=0> K;
  matrix[T, K] X;
  // data for the R2D2 prior
  real<lower=0> mean_R2;  // mean of the R2 prior
  real<lower=0> prec_R2;  // precision of the R2 prior

  int<lower=0> p;
  
  // AR(1) variance parameter for the states
  real<lower=0> phi_sd;
  
  // Variance for intercept
  real<lower=0> alpha_sd;
  
  vector<lower=0>[K] var_x;
}

transformed data{
  real Ymean;
  vector[T] Yc;
  matrix[K, K] xtx;
  real Yc_sd;
  xtx = X' * X;

  Ymean = mean(Y);  
  for (i in 1:T){
    Yc[i] = Y[i] - Ymean;
  }

  Yc_sd = sd(Yc);

  // concentration vector of the Dirichlet prior
  vector<lower=0>[K+1] cons;
  for (i in 1:K) {
    cons[i] = (i^2)/10 * 1/(i^2);
  }
  cons[K+1] = (p^2) / 10;
}

parameters {
  // States
  vector[T] z_tau;
  real z_tau1;
  // Variance Parameters
  real<lower=0> sigma;
  // State AR Weights
  real<lower=-1, upper=1> phi;
  // Regression Weights
  vector[K] z_beta; 
  // R2 Parameters
  vector<lower=0>[1 + K] gamma; // when normalised, yield Dirichlet(concentrations)
  real<lower=0, upper=1> R2;
  // Intercept parameters
  real alpha;
}


transformed parameters{
  real<lower=0> tau2 = R2 / (1 - R2);
  vector<lower=0,upper=1>[1 + K] psi = gamma / sum(gamma);
  vector<lower=0>[K] beta_sd = sqrt(sigma^2 * tau2 * psi[1:K] ./ var_x);
  vector[K] beta = z_beta .* beta_sd;
  real<lower=0> sigma_tau = sqrt(sigma^2 * (1-phi^2) * tau2 * psi[K + 1]);
  vector[T] tau;
  tau[1] = z_tau[1] * sigma_tau;
  
  for (t in 2:T){
    tau[t] = phi * tau[(t-1)] + z_tau[t] * sigma_tau; 
  }
}

model {
  // Priors 
  sigma ~ student_t(3, 0, Yc_sd);

  phi ~ normal(0,phi_sd);
  
  z_beta ~ std_normal();
  z_tau ~ std_normal();
  
  gamma ~ gamma(cons, 1); // prior over unnormalised simplex;
  R2 ~ beta(mean_R2 * prec_R2, (1 - mean_R2) * prec_R2);

  alpha ~ normal(0, alpha_sd); // Intercept
  
  // Observation Likelihood
  Yc ~ normal_id_glm(X, tau + alpha, beta, sigma);
}

generated quantities{
  vector[T] y_pred;
  vector[T] mu_pred;
  vector[T] log_lik;
  real lprior;
  real Intercept;
  real postR2;
  real postR2_trend;
  real postR2_reg;
 
  Intercept = Ymean + alpha;
  mu_pred = Intercept + X * beta + tau;
  postR2 = variance(mu_pred) / (variance(mu_pred) + sigma^2); 
 
  postR2_trend = (tau' * tau) / (Intercept^2 + tau' * tau + beta' * xtx * beta + T * sigma^2);
  postR2_reg = (beta' * xtx * beta) / (Intercept^2 + tau'*tau + beta' * xtx * beta + T * sigma^2);
 

    
}

