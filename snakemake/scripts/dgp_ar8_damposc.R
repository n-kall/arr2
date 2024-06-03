dgp_dampened_oscillations <- function(rep_id, n_obs, phi, sigma, ...) {
  
  # Number of lags
  p <- length(phi)
  # save parameters
  dgp_args <- list(n_obs = n_obs, phi = phi, sigma = sigma)
  
  true_pars <- list(phi = phi, sigma = sigma)
  
  # generate data
  simy <- arima.sim(list(ar = phi), sd = sigma, n = (n_obs+p))
  
  lags <- lagmatrix(simy,1:p)
  lags <- lags[(p+1):(n_obs+p),]
  
  dgp_r2 <- t(phi)%*% t(lags) %*% lags %*% phi / (t(phi)%*% t(lags) %*% lags %*% phi + n_obs * sigma^2) 
  
  true_pars$postR2 <- dgp_r2
  
  return(
    list(
      data = list(Y = simy[(p+1):(n_obs+p)]),
      dgp_args = dgp_args,
      true_pars = true_pars
    )
  )
}

dampened_oscillator <- function(order) {
  coeff <- array(0,c(order,1))
  A0 <- 1
  b <- 1
  m <- 2
  omega <- 2 # was at 2
  phi <- 2
  for (i in 1:order){
    coeff[i] <- A0 * exp(-b / (2 * m) * i) * cos(omega * i + phi)
  }
  return(coeff)
}

dgp <- function(rep_id, ...) {
  set.seed(SEED + rep_id)

  dgp_dampened_oscillations(
    rep_id = rep_id,
    phi = dampened_oscillator(8),
    sigma = 1,
    ...
  )
}
