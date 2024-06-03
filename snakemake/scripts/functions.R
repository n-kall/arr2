bayes_R2 <- function(y_pred,var_res) {
  var_fit <- apply(y_pred, 1, var)
  r2 <- var_fit / (var_fit + var_res)
  return(r2)
}

vroom_chunked <- function(files, altrep = FALSE, show_col_types = FALSE, ...) {
  file_limit <- as.numeric(system("ulimit -n", intern = TRUE)) - 100
  nfiles <- length(files)
  nconnections <- 10
  out <- tibble()
  first_file <- 1
  count <- 0
  while(count < nfiles) {
    chunk_size <- min(file_limit - nconnections, nfiles - first_file + 1)
    last_file <- first_file + chunk_size - 1
    current <- vroom(files[first_file:last_file], altrep = altrep, show_col_types = FALSE, ...)

    probs <- vroom::problems(current)
    if (nrow(probs) > 0) {
      print(probs)
    }

    first_file <- first_file + chunk_size
    count <- count + chunk_size
    out <- bind_rows(out, current)
  }

  out
}

all_metrics <- function(fit, true_pars, zero_padding = TRUE) {

  draws <- fit$draws()

  summary <- posterior::summarise_draws(
    draws,
    mean,
    sd,
    quantile2,
    rhat,
    ess_bulk,
    ess_tail
  )

  matched <- match_var_length(
    true_pars,
    draws,
    zero_padding
  )

  true_vars <- Reduce(c, matched$true_vars)
  true_var_names <- names(true_vars)

  true_par_vals <- tibble(
    variable = true_var_names,
    true_par = true_vars
  )

  diagnostics <- fit$sampler_diagnostics() |>
    summarise_draws(mean, sd)

  # lp <- logprob_score(true_pars, fit)
  rmse <- overall_rmse(true_pars, draws, zero_padding)
  energy <- energy_score(true_pars, draws, zero_padding)
  mae <- overall_mae(true_pars, draws, zero_padding)

  out <- full_join(
    summary,
    diagnostics
  ) |>
    full_join(
      true_par_vals
    ) |>
    mutate(
      rmse = rmse,
      energy_score = energy,
      mae = mae,
      zero_padding = zero_padding
    )

  return(out)
}

##' Energy score calculation
##'
##' @param true_pars named list of true parameter values
##' @param draws draws object of posterior draws
##' @param ... additional arguments passed to functions
##' @return numeric value of energy score
energy_score <- function(true_pars, draws, zero_padding, ...) {

  matched <- match_var_length(
    true_pars,
    draws,
    zero_padding
  )

  true_vars <- Reduce(c, matched$true_vars)
  draws <- matched$draws
  var_names <- names(true_vars)

  # convert to matrix
  draws <- posterior::as_draws_matrix(draws, ...) |>
    # ensure order is correct
    posterior::subset_draws(variable = var_names, ...) |>
    # transpose
    t()
  # calculate score
  scoringRules::es_sample(y = true_vars, dat = draws, w = weights(draws))
}


##' Overall RMSE
##'
##' @param true_pars named list of true parameter values
##' @param draws draws object of posterior draws
##' @param ... additional arguments passed to functions
##' @return numeric value of root mean square error
overall_rmse <- function(true_pars, draws, zero_padding = TRUE, ...) {

  matched <- match_var_length(
    true_pars,
    draws,
    zero_padding
  )

  true_vars <- Reduce(c, matched$true_vars)
  draws <- matched$draws
  var_names <- names(true_vars)

  # convert to matrix
  means <- posterior::as_draws_matrix(draws, ...) |>
    # ensure order is correct
    posterior::subset_draws(variable = var_names, ...) |>
    # column means
    colMeans()

  rmse <- sqrt(sum((true_vars - means)^2 / length(means)))

  return(rmse)

}


##' Overall MAE
##'
##' @param true_pars named list of true parameter values
##' @param draws draws object of posterior draws
##' @param ... additional arguments passed to functions
##' @return numeric value of mean absolute error
overall_mae <- function(true_pars, draws, zero_padding = TRUE, ...) {

  matched <- match_var_length(
    true_pars,
    draws,
    zero_padding
  )

  true_vars <- Reduce(c, matched$true_vars)
  draws <- matched$draws
  var_names <- names(true_vars)

  # convert to matrix
  means <- posterior::as_draws_matrix(draws, ...) |>
    # ensure order is correct
    posterior::subset_draws(variable = var_names, ...) |>
    # column means
    colMeans()

  mae <- sum(abs(true_vars - means) / length(means))

  return(mae)

}

match_var_length <- function(true_variable_values,
                             draws,
                             zero_padding = TRUE) {

  all_draws <- list()
  for (var in names(true_variable_values)) {
    var_draws <- posterior::subset_draws(draws, variable = var)

    n_true_vars <- length(true_variable_values[[var]])
    n_model_vars <- posterior::nvariables(var_draws)

    if (n_true_vars == 0) {
      print("white noise has no true vars")
      true_variable_values[[var]] <- c(
        true_variable_values[[var]],
        rep(0, n_model_vars)
      )
      n_true_vars <- length(true_variable_values[[var]])
    }
    n_extra_model_vars <- n_model_vars - n_true_vars
    n_total_vars <- max(n_true_vars, n_model_vars)

    if (n_extra_model_vars < 0) {
      # more true variables than model variables
      if (zero_padding) {
        # create additional var = zero variables in the model
        zero_draws <- array(
          rep(0, times = posterior::ndraws(var_draws) * abs(n_extra_model_vars)),
          dim = c(
            posterior::ndraws(var_draws) / posterior::nchains(var_draws),
            posterior::nchains(var_draws),
            abs(n_extra_model_vars)
          )
        )
        zero_draws <- posterior::as_draws_array(zero_draws)
        variables(zero_draws) <- paste0(var, "[", (n_model_vars+1):n_total_vars, "]")
        all_draws[[var]] <- posterior::bind_draws(var_draws, zero_draws)
      } else {
        # no zero padding, then just subset the true variables
        true_variable_values[[var]] <- true_variable_values[[var]][1:n_model_vars]
        all_draws[[var]] <- posterior::subset_draws(var_draws, variable = var)
      }
    } else if (n_extra_model_vars > 0) {
      # more model variables than true variables
      if (zero_padding) {
        # if zero padding is on, then add more true variables with
        # zeros
        true_variable_values[[var]] <- c(
          true_variable_values[[var]],
          rep(0, n_total_vars - n_true_vars)
        )
        all_draws[[var]] <- posterior::subset_draws(var_draws, variable = var)
      } else {
        if (n_true_vars > 1) {
          var_subset <- paste0(var, "[", 1:n_true_vars, "]")
        } else {
          var_subset <- var
        }
        all_draws[[var]] <- posterior::subset_draws(var_draws, variable = var_subset)
      }
    } else {
      # number of true vars and model vars equal
      # then just extract the corresponding draws
      if (n_true_vars > 1) {
        var_subset <- paste0(var, "[", 1:n_true_vars, "]")
      } else {
        var_subset <- var
      }
      all_draws[[var]] <- posterior::subset_draws(var_draws, variable = var_subset)
    }

    if (n_true_vars > 1) {
      names(true_variable_values[[var]]) <- paste0(var, "[", 1:length(true_variable_values[[var]]), "]")
    } else {
      names(true_variable_values[[var]]) <- var
    }
  }

  return(list(draws = Reduce(posterior::bind_draws, all_draws), true_vars = true_variable_values))
}

# more stable than log(sum(exp(x)))
log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

# more stable than log(mean(exp(x)))
log_mean_exp <- function(x) {
  log_sum_exp(x) - log(length(x))
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
