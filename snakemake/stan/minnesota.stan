data {

  int<lower=1> T; // number of time points
  vector[T] Y; // observations
  int<lower=0> p; // AR order

  // sigma prior sd
  real<lower=0> sigma_sd; // sd of sigma prior
  
  // intercept prior
  real alpha_mean;
  real<lower=0> alpha_sd;
}


transformed data {
  matrix[(T-p), p] Y_matrix;
  matrix[p,p] YY;
  for(t in 1:(T-p)) {
    for(i in 1:p) {
      Y_matrix[t, p-i+1] = Y[t + (i-1)];
    }
  }
  
  YY = Y_matrix'*Y_matrix;
}

parameters {

  real alpha;
  real<lower=0> kappa2;
  vector[p] phi;
  real<lower=0> sigma;

}

transformed parameters {

  vector[p] phi_sd;
  for (i in 1:p) {
    phi_sd[i] = sqrt(kappa2 / i^2);
  }
}


model {

  // priors
  alpha ~ normal(alpha_mean, alpha_sd);
  sigma ~ normal(0, sigma_sd);
  kappa2 ~ gamma(1, 1/0.04);
  for (i in 1:p) {
    phi[i] ~ normal(0, sqrt(kappa2 / i^2));
  }

  // likelihood
  Y[(p+1):T] ~ normal_id_glm(Y_matrix, alpha, phi, sigma);

}

generated quantities{
  vector[(T-p)] mu = alpha +  Y_matrix*phi;
  real<lower=0, upper=1> postR2 = variance(mu) / (variance(mu) + sigma^2 );
  vector[p] rel_R2;
  
  for (i in 1:p){
    rel_R2[i] =  phi[i]^2 * YY[i,i] / (phi'*YY*phi + (T-p) * sigma^2);
  }
  }

