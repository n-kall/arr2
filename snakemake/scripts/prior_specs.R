minn_sig_create <- function(X, lags) {
  X <- as.matrix(X)
  T <- dim(X)[1]
  K <- dim(X)[2]
  # Outputs
  var_x <- array(0,c(K,1))
  for (j in 1:K){
    # Create Lags
    Y_matrix <- array(0,c((T - lags), lags))
    for (t in 1:(T - lags)) {
      for (i in 1:lags) {
        Y_matrix[t, (lags-i+1)] <- X[(t + (i-1)),j]
      }
    }
    # Create AR(lags) Residual Variances
    Y <- X[(lags + 1):T, j]
    beta <- solve(t(Y_matrix) %*% Y_matrix) %*% t(Y_matrix) %*% Y
    var_x[j] <- var(Y-Y_matrix %*% beta)
  }
  return(var_x)
}

minn_sig_create_ltx <- function(X,lags){
  X <- as.matrix(X)
  T <- dim(X)[1]
  K <- dim(X)[2]
  
  # Outputs
  var_x <- array(0,c(K,1))
  
  for (j in 1:K){
    # Create Lags
    Y_matrix <- array(0,c((T-lags), lags))
    for(t in 1:(T-lags)) {
      for(i in 1:lags) {
        Y_matrix[t, (lags-i+1)] = X[(t + (i-1)),j];
      }
    }
    # Create AR(p) Residual Variances
    Y <- X[(lags+1):T,j]
    beta <- solve(t(Y_matrix)%*%Y_matrix)%*%t(Y_matrix)%*%Y
    var_x[j] <- var(Y-Y_matrix%*%beta)
  }
  
  return(var_x)
}

gaussian_prior <- list(
  phi_sd = 1,
  sigma_tau_sd = 5,
  sigma_tau1_sd = 5
)

all_prior <- list(
  sigma_sd = 5,
  alpha_mean = 0,
  alpha_sd = 10,
  beta_sd = 1
)

arr2_prior <- list(
  mean_R2 = 0.33,
  prec_R2 = 3,
  alpha_cent = 1,
  phi_cent = 0,
  sigma_cent = 1
)

rhs_pars <- list(
  hs_df = 1,
  hs_df_global = 1,
  hs_df_slab = 4,
  hs_scale_slab = 2
)

prior_settings <- c(
  all_prior, gaussian_prior, arr2_prior, rhs_pars
)
