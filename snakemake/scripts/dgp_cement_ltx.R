# Function for generating the lags

Xlagsim <- function(X, p, K, T) {
  X_all <- array(0, c((T - p), K * p))
  for (j in 1:K){
    lagtemp <- lagmatrix(X[, j], 1:p)
    X_all[, (1 + (j - 1) * p):(j * p)] <- lagtemp[(p + 1):T, ]
  }
  return(X_all)
}


dgp_cement <- function(rep_id, n_obs, p, m, rho, state_phi, state_snr, ...) {
  
  # save parameters
  dgp_args <- list(n_obs = n_obs, p = p, m = m, rho = rho, state_phi = state_phi, state_snr = state_snr)

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
  w1 <- rep(xi/(1:p)^2, block_size)
  w2 <- rep(xi * 0.5/(1:p)^2, block_size)
  w3 <- rep(xi * 0.25/(1:p)^2, block_size)
  w4 <- rep(0, (num_matrices - 3) * block_size * p)
  w <- c(w1, w2, w3, w4)/3

  true_pars <- list(beta = w, p = p, state_phi = state_phi, state_snr = state_snr)

  # define zero mean vector
  sigma <- 1
  mu <- rep(0, m)

  # generate X and  y
  X <- MASS::mvrnorm(n = (n_obs + p), mu = mu, Sigma = R)
  y <- array(0, c((n_obs + p), 1))

  # Placeholder for tau
  tau <- array(0, c((n_obs + p), 1))

  # Generate the lags
  X <- Xlagsim(X, p, m, n_obs + p)


  for (i in 2:n_obs){
    tau[i] <- tau[(i - 1)] * state_phi + rnorm(1) * state_snr * 1/sigma
    y[i] <- tau[i] + X[i, ] %*% w + sigma * rnorm(1)
  }

  y <- y[1:n_obs]
  tau <- tau[1:n_obs]

  dgp_r2 <- (t(w) %*% t(X) %*% X %*% w + t(tau) %*% tau) /
    (t(w) %*% t(X) %*% X %*% w + t(tau) %*% tau + n_obs * sigma^2)

  true_pars$postR2 <- dgp_r2

  true_pars$tau <- tau

  true_pars$sigma_tau <- state_snr * 1 / sigma

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
    p = 12,
    state_phi = 0.95,
    ...
  )
}

