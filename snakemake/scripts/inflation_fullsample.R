
## Full sample fit
library(cmdstanr)
library(resample)
library(posterior)
library(scoringRules)
library(readr)
# Load Data

dataset <- readRDS(snakemake@input[["dataset"]])
functions <- snakemake@input[["functions"]]
prior_specs <- snakemake@input[["prior_specs"]]
exe_file <- snakemake@input[["exe_file"]]
output_file <- snakemake@output[["fit_object"]]
SEED <- as.numeric(snakemake@config[["seed"]])
prior_only <- as.numeric(snakemake@wildcards[["priorpost"]] == "prior")
source(functions)

# Forecast Preliminaries
chains <- 4
mcmc <- 1000
K <- dim(dataset$X)[2]
p <- 12
K_ind <- 1:K


# Create Data List
dat <- list(
  T = length(dataset$y),
  y = as.vector(dataset$y),
  Y = as.vector(dataset$y),
  K = K,
  K_ltx = K,
  X = scale(dataset$X[,K_ind]),
  X_ltx = scale(dataset$X[,K_ind]),
  mean_R2 = 0.33,
  prec_R2 = 3,
  cons = c(rep(1,p)),
  phi_sd = 2,
  alpha_sd = sd(as.vector(dataset$y)),
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
  var_x_minn = as.numeric(minn_sig_create(scale(dataset$X[,seq(1, K, by = p)]), p)),
  var_y_minn = as.numeric(minn_sig_create(dataset$y,p)),
  var_x = as.numeric(colVars(scale(dataset$X[,K_ind]))),
  var_x_ltx = as.numeric(colVars(scale(dataset$X[,K_ind]))),
  var_y = as.numeric(var(as.vector(dataset$y))),
  prior_only = prior_only
)

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

R2_draws <- fit$draws(variable = "R2_pp") |>
  merge_chains() |>
  mutate_variables(priorpost = prior_only) |>
  as_draws_df()

write_csv(R2_draws, file = output_file)
