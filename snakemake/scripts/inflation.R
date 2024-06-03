
## Iterative Forecast Script: ltx ##
library(cmdstanr)
library(resample)
library(posterior)
library(scoringRules)
# Load Data

dataset <- readRDS(snakemake@input[["dataset"]])
functions <- snakemake@input[["functions"]]
prior_specs <- snakemake@input[["prior_specs"]]
exe_file <- snakemake@input[["exe_file"]]
output_file <- snakemake@output[["fit_object"]]
fold <- as.numeric(snakemake@wildcards[["fold"]])
SEED <- as.numeric(snakemake@config[["seed"]])
reparam <- 1
source(functions)

# Forecast Preliminaries
rwindow <- 240
T_all <- length(dataset$y)
iter <- T_all - rwindow
chains <- 4
mcmc <- 1000
K <- dim(dataset$X)[2]
p <- 12
K_ind <- 1:K


# Create Data List
t_start <- fold
t_end <- fold + rwindow - 1
dat <- list(
  T = length(dataset$y[t_start:t_end]),
  y = as.vector(dataset$y[t_start:t_end]),
  Y = as.vector(dataset$y[t_start:t_end]),
  K = K,
  K_ltx = K,
  X = scale(dataset$X[t_start:t_end,K_ind]),
  X_ltx = scale(dataset$X[t_start:t_end,K_ind]),
  Xf = as.numeric(scale(dataset$X[1:(t_end + 1),])[(t_end+1),K_ind]),
  Xf_ltx = as.numeric(scale(dataset$X[1:(t_end + 1),])[(t_end+1),K_ind]),
  yf = dataset$y[(t_end + 1)],
  mean_R2 = 0.33,
  prec_R2 = 3,
  cons = c(rep(1,p)),
  phi_sd = 2,
  alpha_sd = sd(as.vector(dataset$y[t_start:t_end])),
  p = p,
  beta_sd = 2,
  hs_df = 3,
  hs_df_global = 1,
  hs_df_slab = 4,
  hs_scale_global = 1,
  hs_scale_slab = 2,
  p0 = round((K+p)/2),
  sigma_tau_sd = 1,
  sigma_tau1_sd = 1,
  var_x_minn = as.numeric(minn_sig_create(scale(dataset$X[t_start:t_end,seq(1, K, by = p)]), p)),
  var_y_minn = as.numeric(minn_sig_create(dataset$y[t_start:t_end],p)),
  var_x = as.numeric(colVars(scale(dataset$X[t_start:t_end,K_ind]))),
  var_x_ltx = as.numeric(colVars(scale(dataset$X[t_start:t_end,K_ind]))),
  var_y = as.numeric(var(as.vector(dataset$y[t_start:t_end])))
)

dat$centeredness <- rep(0, times = dat$T + dat$K)

# Run Model
mod <- cmdstan_model(exe_file = exe_file)
fit <- mod$sample(
  data = dat,
  seed = SEED,
  chains = chains,
  iter_warmup = mcmc,
  iter_sampling = mcmc,
  adapt_delta = 0.99,
  max_treedepth = 15,
  parallel_chains = 4,
  )

if (reparam == 1) {
  print("reparametersing on")
  num_div_before <- sum(fit$diagnostic_summary()$num_divergent)
  if (num_div_before > -5) {
    print("trying reparemeterisation")
    d <- fit$draws()
    d <- merge_chains(d)

    log_sigma_tau1 <- t(as_draws_matrix(subset_draws(merge_chains(d), "log_sigma_tau1")))
    class(log_sigma_tau1) <- "matrix"
    log_sigma_tau <- t(as_draws_matrix(subset_draws(merge_chains(d), "log_sigma_tau")))
    class(log_sigma_tau) <- "matrix"
    log_beta_sd <- t(as_draws_matrix(subset_draws(merge_chains(d), "log_beta_sd")))
    class(log_beta_sd) <- "matrix"
    tau_post <- t(as_draws_matrix(subset_draws(merge_chains(d), "tau")))
    class(tau_post) <- "matrix"
    beta_post <- t(as_draws_matrix(subset_draws(merge_chains(d), "beta")))
    class(beta_post) <- "matrix"
    
    log_sigma_tau_all <- rbind(log_sigma_tau1, do.call(rbind, replicate(dat$T, log_sigma_tau, simplify = FALSE)))
    
    log_scales <- rbind(log_sigma_tau_all, log_beta_sd)
    parameters <- rbind(tau_post, beta_post)
    
    dat$centeredness <- fix_centeredness(log_scales, parameters)
    fit <- mod$sample(data = dat, chains = 4, parallel_chains = 4, seed = SEED,   adapt_delta = 0.99,
  max_treedepth = 15, iter_sampling = mcmc, iter_warmup = mcmc
)
  }
}

# Save output
fit$save_object(file = output_file)
