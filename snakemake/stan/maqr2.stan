data {
  int<lower=0> Q;       // num previous noise terms
  int<lower=3> T;       // num observations
  vector[T] y;          // observation at time t
  
  // data for the R2D2 prior
  real<lower=0> mean_R2;  // mean of the R2 prior
  real<lower=0> prec_R2;  // precision of the R2 prior
  vector<lower=0>[Q] cons; // concentrations on the R2 prior

}
parameters {
  real mu;              // mean
  real<lower=0> sigma;  // error scale
  vector<lower=0>[Q] gamma;
  real<lower=0, upper=1> R2;
  vector[Q] phi_z;
}
transformed parameters {
  real<lower=0> tau2 = R2 / (1 - R2);
  vector<lower=0,upper=1>[Q] psi = gamma / sum(gamma);
  vector[Q] theta = phi_z .* sqrt(tau2 * psi);
  vector[T] epsilon;    // error term at time t
  for (t in 1:T) {
    epsilon[t] = y[t] - mu;
    for (q in 1:min(t - 1, Q)) {
      epsilon[t] = epsilon[t] - theta[q] * epsilon[t - q];
    }
  }
}

model {
  vector[T] eta;
  mu ~ normal(0, 1);
  gamma ~ gamma(cons, 1); // prior over unnormalised simplex
  R2 ~ beta_proportion(mean_R2, prec_R2);
  phi_z ~ std_normal();
  sigma ~ normal(0, 2.5);
  for (t in 1:T) {
    eta[t] = mu;
    for (q in 1:min(t - 1, Q)) {
      eta[t] = eta[t] + theta[q] * epsilon[t - q];
    }
  }
  y ~ normal(eta, sigma);
}

