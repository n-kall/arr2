functions {
  vector repeat_vector(vector input, int K) {
    int N = rows(input);
    vector[N*K] repvec; // stack N-vector K times

    for (k in 1:K) {
      for (i in 1:N) {
	repvec[i+(k-1)*N] = input[i]; // assign i-th value of input to i+(k-1)*N -th value of repvec
      }
    }
    return repvec;
  }
}

data {
  int<lower=1> T;
  vector[T] Y;
  int<lower=0> K;
  matrix[T,K] X;
  int<lower = 1> p; // number of lags for each of the covariates

  // data for the R2D2 prior
  real<lower=0> mean_R2;  // mean of the R2 prior
  real<lower=0> prec_R2;  // precision of the R2 prior

  // AR(1) variance parameter for the states
  real<lower=0> phi_sd;

  // Variance for intercept
  real<lower=0> alpha_sd;

  vector<lower=0>[K] var_x;

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
  int<lower=1> K_dc = K/p;
  vector<lower=0>[(K_dc+1)] cons; // concentration vector of the Dirichlet prior
  vector[p] Ldc; // decomposition of the weight for each covariate
  vector[p] Ldc_partial;

  xtx = X' * X;

  Ymean = mean(Y[1:Tc]);
  for (i in 1:Tc){
    Yc[i] = Y[i] - Ymean;
  }

  Yc_sd = sd(Yc);


  for (i in 1:(K_dc+1)) {
    cons[i] = 1;
  }

  for (i in 1:p) {
    Ldc_partial[i] = 1/(i^2);
  }

    for (i in 1:(p)) {
    Ldc[i] = (1/(i^2)) / sum(Ldc_partial) ;
  }
}

parameters {
  // States
  vector[Tc] z_tau;
  real z_tau1;
  // Variance Parameters
  real<lower=0> sigma;
  // State AR Weights
  real<lower=0,upper=1> phi;
  // Regression Weights
  vector[K] z_beta;
  // R2 Parameters
  simplex[(K_dc+1)] psi;
  real<lower=0, upper=1> R2;
  // Intercept parameters
  real alpha;
}


transformed parameters{
  real<lower=0> tau2 = R2 / (1 - R2);
  real<lower=0> sigma_tau = sqrt(sigma^2 * (1-phi^2) * tau2 * psi[(K_dc+1)]);
  vector[K] beta;
  vector<lower=0>[K] beta_sd;
  vector[Tc] tau;

  for (j in 1:(K_dc)){
    beta_sd[(p*j-p+1):(p*j)] = sqrt(sigma^2 * tau2 * psi[j] * Ldc  ./ var_x[(p*j-p+1):(p*j)]);
  }
  beta  = z_beta .* beta_sd;

  tau[1] = z_tau[1] * sigma_tau;

  for ( t in 2:Tc) {
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
