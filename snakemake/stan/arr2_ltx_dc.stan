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
  int<lower = 1> p; // numer of lags for each of the covariates

  // data for the R2D2 prior
  real<lower=0> mean_R2;  // mean of the R2 prior
  real<lower=0> prec_R2;  // precision of the R2 prior

  // AR(1) variance parameter for the states
  real<lower=0> phi_sd;

  // Variance for intercept
  real<lower=0> alpha_sd;

  vector<lower=0>[K] var_x;
}

transformed data{
  real Ymean;
  real sd_Yc;
  vector[T] Yc;
  matrix[K,K] xtx;
  int<lower= 1> K_dc = K/p;
  vector<lower=0>[(K_dc+1)] cons; // concentration vector of the Dirichlet prior
  vector[p] Ldc; // decomposition of the weight for each covariate
  vector[p] Ldc_partial;

  xtx = X' * X;

  Ymean = mean(Y);
  for (i in 1:T){
    Yc[i] = Y[i] - Ymean;
  }

  sd_Yc = sd(Yc);


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
  vector[T] z_tau;
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
  vector[T] tau;

  for (j in 1:(K_dc)){
    beta_sd[(p*j-p+1):(p*j)] = sqrt(sigma^2 * tau2 * psi[j] * Ldc  ./ var_x[(p*j-p+1):(p*j)]);
  }
  beta  = z_beta .* beta_sd;

  tau[1] = z_tau[1] * sigma_tau;

  for ( t in 2:T){
    tau[t] = phi * tau[(t-1)] + z_tau[t] * sigma_tau;
  }
}

model {
  // Priors
  sigma ~ student_t(3,0, sd_Yc);

  z_beta ~ std_normal();
  z_tau ~ std_normal();

  phi ~ normal(0,phi_sd);

  psi ~ dirichlet(cons);
  R2 ~ beta(mean_R2 * prec_R2, (1 - mean_R2) * prec_R2);

  // Observation Likelihood
  target += normal_lpdf(alpha | 0, alpha_sd); // Intercept
  Yc ~ normal_id_glm(X,tau+alpha,beta, sigma);
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
  vector<lower=0>[p] Ldc_view;

  Ldc_view = Ldc;

  lprior = beta_lpdf(R2|mean_R2 * prec_R2, (1 - mean_R2) * prec_R2);
  Intercept = Ymean+alpha;
  mu_pred = Intercept +  X*beta + tau;
  postR2 = variance(mu_pred) / (variance(mu_pred) + sigma^2 );

  postR2_trend = (tau'*tau)/(Intercept^2 + tau'*tau + beta'*xtx*beta + T*sigma^2);
  postR2_reg = (beta'*xtx*beta)/(Intercept^2 + tau'*tau + beta'*xtx*beta + T*sigma^2);
}
