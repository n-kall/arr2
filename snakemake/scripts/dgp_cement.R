dgp_cement <- function(rep_id, n_obs, phi, m, rho, ...) {

  # save parameters
  dgp_args <- list(n_obs = n_obs, phi, m = m, rho = rho)

  # Generate data according to Piironen and Vehtari 2017
  # define xi given rho
  if (rho == 0) {
    xi <- 0.59
  } else if (rho == 0.5) {
    xi <- 0.34
  } else if (rho == 0.9) {
    xi <- 0.28
  } else {
    warning("Invalid rho provided.")
  }

  # build correlation matrix
  block_size <- 5
  num_matrices <- m / block_size
  listOfMatrices <- vector("list", num_matrices)
  for (i in 1:num_matrices) {
    listOfMatrices[[i]] <- matrix(
      rep(rho, block_size * block_size), nrow = block_size
    )
  }
  R <- matrix(Matrix::bdiag(listOfMatrices), nrow = m)
  diag(R) <- 1

  # build associated covariate weights
  w1 <- rep(xi, block_size)
  w2 <- rep(xi * 0.5, block_size)
  w3 <- rep(xi * 0.25, block_size)
  w4 <- rep(0, (num_matrices - 3) * block_size)
  w <- c(w1, w2, w3, w4)

  true_pars <- list(beta = w, phi = phi)

  # define zero mean vector
  sigma <- 1
  mu <- rep(0, m)

  # generate X and y
  X <- MASS::mvrnorm(n = (n_obs+length(phi)), mu = mu, Sigma = R)
  y <- array(0,c((n_obs + length(phi)), 1))

  for (i in (length(phi) + 1):(n_obs + length(phi))){
    lagtemp <- lagmatrix(y, 1:length(phi))
    y[i] <- lagtemp[i,] %*% phi + X[i,] %*% w + rnorm(1) * sigma
  }

  lags <- lagtemp[(length(phi) + 1):(n_obs + length(phi)), ]

  X <- X[(length(phi) + 1):(n_obs + length(phi)),]
  y <- y[(length(phi) + 1):(n_obs + length(phi))]


  dgp_r2 <- (t(w) %*% t(X) %*% X %*% w + t(phi) %*% t(lags) %*% lags %*% phi) /
    (t(w) %*% t(X) %*% X %*% w + t(phi) %*% t(lags) %*% lags %*% phi + n_obs * sigma^2)

  true_pars$postR2 <- dgp_r2

  return(
    list(
      rep_id = rep_id,
      dgp_args = dgp_args,
      true_pars = true_pars,
      data = list(X = X, Y = y)
    )
  )
}

dgp <- function(rep_id, ...) {
  set.seed(SEED + rep_id)
  dgp_cement(
    rep_id = rep_id,
    phi = c(0.6, 0.15, 0.067, 0.038, 0.024, 0.017, 0.012, 0.009),
    m = 400,
    ...
  )
}
