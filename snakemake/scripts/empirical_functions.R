
# more stable than log(sum(exp(x)))
log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

# more stable than log(mean(exp(x)))
log_mean_exp <- function(x) {
  log_sum_exp(x) - log(length(x))
}


Xlagsim <- function(p,K,T) {
  f <- rnorm(T)
  X_all <- array(0,c((T-p),K*p))
  for (j in 1:K){
    rho <- 0.3 + 0.1*rnorm(1)
    x_temp <- sqrt(rho)*f + sqrt(1-rho)*rnorm(T)
    lagtemp <- lagmatrix(x_temp,1:p)
    X_all[,(1 + (j-1)*p ): (j*p)] <- lagtemp[(p+1):T,]
  }
  return(X_all)
}

minn_sig_create <- function(X,lags) {
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
        Y_matrix[t, (lags-i+1)] <- X[(t + (i-1)),j]
      }
    }
    # Create AR(p) Residual Variances
    Y <- X[(lags+1):T, j]
    beta <- solve(t(Y_matrix)%*%Y_matrix)%*%t(Y_matrix)%*%Y
    var_x[j] <- var(Y-Y_matrix%*%beta)
  }
  return(var_x)
}

# get matrix of lags
lagmatrix <- function(x, lag) {
  
  # Construct matrix with lead and lags
  n <- length(x)
  k <- length(lag)
  
  # How much to expand for leads and lags
  mlg <- max(c(0,lag[lag>0]))
  mld <- max(abs(c(0,lag[lag<0])))
  
  # Assign values
  lmat <- array(NA,c(n+mlg+mld,k))
  for (i in 1:k){
    lmat[(1+lag[i]+mld):(n+lag[i]+mld),i] <- x
  }
  
  # Trim lmat for expansion
  lmat <- lmat[(mld+1):(mld+n),,drop=FALSE]
  
  return(lmat)
  
}
