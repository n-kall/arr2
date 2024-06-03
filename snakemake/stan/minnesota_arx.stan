data {

  int<lower=1> T; // number of time points
  vector[T] Y; // observations
  int<lower=0> p; // AR order

  int<lower = 0> K; // number of covariates
  matrix[T,K] X; //covariate matrix

  // sigma prior sd
  real<lower=0> sigma_sd; // sd of sigma prior

  // intercept prior
  real alpha_mean;
  real<lower=0> alpha_sd;

  // AR(4) residual variances of y and X
  real<lower=0> var_y_minn;
  vector<lower=0>[K] var_x_minn;

}

transformed data {
  matrix[(T-p), p] Y_matrix;
  matrix[(T-p),K] X_matrix = X[(p+1):T,];
  matrix[(T-p),K + p] XY_matrix;
  for(t in 1:(T-p)) {
    for(i in 1:p) {
      Y_matrix[t, p-i+1] = Y[t + (i-1)];
    }
    XY_matrix[t] = append_col(X_matrix[t],Y_matrix[t]);
  }
}

parameters {

  real alpha;
  real<lower=0> kappa2_lags;
  vector[p] phi;
  real<lower=0> sigma;
  real<lower=0> kappa2_k;
  vector[K] beta;

}

model {
  // priors
  alpha ~ normal(alpha_mean, alpha_sd);
  sigma ~ normal(0, sigma_sd);
  kappa2_lags ~ gamma(1, 1/0.04);
  kappa2_k ~ gamma(1,1/(0.04^2));
  for (i in 1:p) {
    phi[i] ~ normal(0, sqrt(kappa2_lags / i^2));
  }
  for (j in 1:K){
    beta[j] ~ normal(0, sqrt(var_y_minn / var_x_minn[j] * kappa2_k));
  }
  // likelihood
  Y[(p+1):T] ~ normal_id_glm(XY_matrix, alpha, append_row(beta, phi), sigma);
}


generated quantities{
  vector[(T-p)] mu = alpha +  XY_matrix*append_row(beta, phi);
  real<lower=0, upper=1> postR2 = variance(mu) / (variance(mu) + sigma^2 );
  //  real<lower=0, upper=1> postR2_lags = (phi'*Y_matrix'*Y_matrix*phi)/(append_row(beta, phi)'*XY_matrix'*XY_matrix*append_row(beta, phi));
  //  real<lower=0, upper=1> postR2_reg = (beta'*X_matrix'*X_matrix*beta)/(append_row(beta, phi)'*XY_matrix'*XY_matrix*append_row(beta, phi));
  }
