functions {
  /** ARR2 prior
   *
   * @param phi Real
   * @param psi Simplex
   * @param R2 Real
   * @param sigma Real
   * @param sigma_sd Real
   * @param mean_R2 Real
   * @param prec_R2 Real
   * @param cons Vector
   * @param var_y Real
   */
  real arr2_lpdf(vector phi, vector psi, real R2, real sigma,
		 data real sigma_sd, data real mean_R2,
		 data real prec_R2, data vector cons, data real var_y) {
    return normal_lpdf(phi | 0, sqrt(sigma^2/var_y * (R2 / (1 - R2)) * psi)) +
      beta_lpdf(R2 | mean_R2 * prec_R2, (1 - mean_R2) * prec_R2) +
      normal_lpdf(sigma | 0, sigma_sd) +
      dirichlet_lpdf(psi | cons);
  }
}

data {
  int<lower=1> T; // number of time points
  vector[T] Y; // observations
  int<lower=0> p; // AR order
  // concentration vector of the Dirichlet prior
  vector<lower=0>[p] cons;
  // data for the R2D2 prior
  real<lower=0> mean_R2; // mean of the R2 prior
  real<lower=0> prec_R2; // precision of the R2 prior
  real<lower=0> sigma_sd; // sd of sigma prior
}

transformed data {
  // Variance estimate of y
  real<lower=0> var_y;
  var_y = variance(Y);
}

parameters {
  vector[p] phi; // AR coefficients
  simplex[p] psi; // decomposition simplex
  real<lower=0, upper=1> R2; // coefficient of determination
  real<lower=0> sigma; // observation model sd
  real alpha; // intercept
}
transformed parameters {
  vector[T] mu = rep_vector(0.0, T);
  for (t in (p+1):T) {
    for (i in 1:p) {
      mu[t] += alpha + phi[i] * Y[t-i];
    }
  }
}
model {
  // priors
  target += arr2_lpdf(phi | psi, R2, sigma, sigma_sd, mean_R2, prec_R2, cons, var_y);
  target += normal_lpdf(alpha | 0, 1);
  // likelihood
  target += normal_lpdf(Y | mu, sigma);
}
