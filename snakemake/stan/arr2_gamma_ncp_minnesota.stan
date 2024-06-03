data {

  int<lower=1> T; // number of time points
  vector[T] Y; // observations
  int<lower=0> p; // AR order

  // data for the R2D2 prior
  real<lower=0> mean_R2;  // mean of the R2 prior
  real<lower=0> prec_R2;  // precision of the R2 prior

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

  // Variance estimate of y
  real<lower=0> var_y;
  var_y = variance(Y);

  // concentration vector of the Dirichlet prior
  vector<lower=0>[p] cons;
  for (i in 1:p) {
    cons[i] = (p^2)/10 * 1/(i^2);
  }
  YY = Y_matrix'*Y_matrix;
}

parameters {
  real alpha; // intercept
  real<lower=0, upper=1> R2;
  real<lower=0> sigma;
  vector[p] phi_z;
  vector<lower=0>[p] gamma;
}

transformed parameters {
  real<lower=0> tau2 = R2 / (1 - R2);
  vector<lower=0,upper=1>[p] psi = gamma / sum(gamma);
  vector[p] phi_sd = sqrt(sigma^2 * tau2 * psi);
  vector[p] phi = phi_z .* sqrt(sigma^2/var_y * tau2 * psi);
}


model {

  // priors
  alpha ~ normal(alpha_mean, alpha_sd);
  sigma ~ normal(0, sigma_sd);
  gamma ~ gamma(cons, 1); // prior over unnormalised simplex
  R2 ~ beta_proportion(mean_R2, prec_R2);
  phi_z ~ std_normal();

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
  
  

