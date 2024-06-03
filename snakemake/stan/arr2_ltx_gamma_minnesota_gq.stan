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

  int<lower=1> L; // last time point for training
  int<lower=1> Tc; // train data length
  int<lower=0> oos; // out-of-sample data index

}

transformed data {
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


  // concentration vector of the Dirichlet prior
  vector<lower=0>[K+1] cons;
  for (i in 1:K+1) {
    cons[i] = (i^2)/ 10 * 1/(i^2);
  }
  
}

parameters {
  // States
  vector[Tc] z_tau;
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
  vector[Tc] tau;
  tau[1] = z_tau[1] * sigma_tau;

  for (t in 2:Tc){
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
