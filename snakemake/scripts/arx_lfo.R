library(cmdstanr)
library(posterior)
library(tidyverse)

dataset <- readRDS(snakemake@input[["dataset"]])
functions <- snakemake@input[["functions"]]
prior_specs <- snakemake@input[["prior_specs"]]
exe_file <- snakemake@input[["exe_file"]]
output_file <- snakemake@output[["fit_summary"]]
SEED <- as.numeric(snakemake@config[["seed"]])
npars_x <- as.numeric(snakemake@wildcards[["m"]])

source(functions)
source(prior_specs)

# refit a model on new data
update_model <- function(model, newx, newy, standata) {
  standata_copy <- standata
  standata_copy$X <- newx
  standata_copy$Y <- newy
  standata_copy$T <- length(newy)

  fit <- model$sample(data = standata_copy, parallel_chains = 4, iter_warmup = 1000, iter_sampling = 4000, refresh = 0, seed = SEED)

  return (fit)
}

# extract a model's log-likelihood
extract_model_gq_draws <- function(fit, model_gq, standata, oos, L) {
  gq_standata <- c(standata, list(oos = oos, L = L))
  fit_gq <- model_gq$generate_quantities(
    fit,
    data = gq_standata,
    parallel_chains = 4
  )
  return ( posterior::merge_chains(fit_gq$draws()) )
}

# primary function of the experiment: compute approx. 1SAP LFO-CV
compute_lfo <- function(stanfile, stanfile_gq, standata, L) {
  # retrieve model
  model <- cmdstan_model(exe_file = stanfile)

  # fit generated quantities model
  model_gq <- cmdstan_model(exe_file = stanfile_gq)
  fit <- model$sample(
    data = standata,
    refresh = 0,
    iter_warmup = 1000,
    iter_sampling = 4000,
    parallel_chains = 4,
    seed = SEED
  )

  # extract dimensions of the problem
  N <- standata$T

  gq_full_draws <- extract_model_gq_draws(
    fit = fit,
    model_gq = model_gq,
    standata = standata,
    oos = 0,
    L = L
  )

  loglik_exact_full <- posterior::extract_variable(gq_full_draws, "log_lik")

  # compute exact full lpd
  exact_lpd_full <- log_mean_exp(loglik_exact_full)

  # initialise log-likelihood matrix
  num_draws <- ncol(fit$draws()) * nrow(fit$draws())
  loglik_exact <- matrix(nrow = num_draws, ncol = N)
  pred_exact <- list()


  rhats <- list()
  
  # iterate over
  for (i in L:(N - 1)) {
    past <- 1:i
    oos <- i + 1
    Y_past <- standata$Y[past, drop = FALSE]
    X_past <- standata$X[past,]
    fit_i <- update_model(model = model, newx = X_past, newy = Y_past, standata = standata)
    gq_draws <- extract_model_gq_draws(
      fit = fit_i,
      model_gq = model_gq,
      standata = standata,
      oos = oos,
      L = L
    )

    max_rhat <- max(summarise_draws(fit_i$draws(), rhat)$rhat)
    rhats <- c(rhats, list(max_rhat))

    loglik_exact[, i + 1] <- posterior::extract_variable_matrix(gq_draws, "log_lik")

    current_pred <- posterior::subset_draws(gq_draws, "pred")
    variables(current_pred) <- paste0("pred[", i + 1 - L, "]")
    pred_exact <- c(pred_exact, list(current_pred))
  }

  pred_exact <- bind_draws(pred_exact)

  # compute exact LFO-CV elpd
  exact_elpds_1sap <- apply(loglik_exact, 2, log_mean_exp)
  
  exact_elpd_1sap <- sum(exact_elpds_1sap[-(1:L)])
  folds <- nrow(exact_elpds_1sap)

  out <- tibble(
    elpds = exact_elpds_1sap[-(1:L)],
    fold = folds[-(1:L)],
    max_rhat = max_rhat
  )

  return(out)
}

standata <- c(prior_settings, dataset$data)
standata$T <- length(standata$Y)
standata$X <- standata$X[,1:npars_x]
standata$p <- 12
standata$K <- ncol(standata$X)

standata$var_y_minn <- as.numeric(minn_sig_create(standata$Y, 4))
standata$var_x_minn <- as.numeric(minn_sig_create(standata$X, 4))

standata$var_y <- as.numeric(var(standata$Y))

standata$var_x <- matrixStats::colVars(standata$X)


standata$p0 <- floor(standata$p / 2)

# run experiment
results <- compute_lfo(
  stanfile = exe_file,
  stanfile_gq = paste0(exe_file, "_gq"),
  standata = standata,
  L = floor(120 / 2) + 1 # initial n_obs for LFO
)

# save results
write_csv(as_tibble(results), file = output_file)
