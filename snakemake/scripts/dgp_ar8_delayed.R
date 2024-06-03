dgp_ar8 <- function(rep_id, n_obs, phi, sigma, ...) {
  
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
      rep_id = rep_id,
      data = list(Y = simy[(p+1):(n_obs+p)]),
      dgp_args = dgp_args,
      true_pars = true_pars
    )
  )
}

dgp <- function(rep_id, ...) {
  set.seed(SEED + rep_id)

  dgp_ar8(
    rep_id = rep_id,
    phi = c(0, 0, 0, 0, 0.7, 0.2, 0.05, 0.025),
    u = 1,
    sigma = 1,
    ...
  )
}
