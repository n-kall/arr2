data {

  int<lower=1> T; // number of time points
  vector[T] Y; // observations
  int<lower=0> p; // AR order

  // sigma prior sd
  real<lower=0> sigma_sd; // sd of sigma prior

  // phi prior sd
  real<lower=0> phi_sd; // sd of phi prior

  // intercept prior
  real alpha_mean;
  real<lower=0> alpha_sd;

  int<lower=0, upper=1> alpha_cent;
  int<lower=0, upper=1> phi_cent;
  int<lower=0, upper=1> sigma_cent;
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

  vector<multiplier=(phi_cent < 1 ? phi_sd : 1.0)>[p] phi;
  real<offset=(alpha_cent < 1 ? alpha_mean : 0.0),
        multiplier=(alpha_cent < 1 ? alpha_sd : 1.0)> alpha; // intercept
  real<lower=0> sigma;

}

model {

  // priors
  alpha ~ normal(alpha_mean, alpha_sd);
  sigma ~ normal(0, sigma_sd);
  phi ~ normal(0, phi_sd);

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
