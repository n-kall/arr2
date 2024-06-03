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

prior_settings <- c(
  all_prior, gaussian_prior, arr2_prior, rhs_pars
)

dat <- c(prior_settings, dataset$data)
dat$X <- dat$X[,1:npars_x]
dat$T <- length(dat$Y)

dat$var_y_minn <- as.numeric(minn_sig_create(dat$Y, 4))
dat$var_x_minn <- as.numeric(minn_sig_create(dat$X, 4))

dat$var_y <- as.numeric(var(dat$Y))

dat$var_x <- matrixStats::colVars(dat$X)

dat$p <- 12

dat$K <- npars_x


print("assigned concentration")

dat$p0 <- floor(dat$p / 2)

# load model
model <- cmdstan_model(exe_file = exe_file)

# fit model
fit <- model$sample(data = dat, chains = 4, parallel_chains = 4, seed = SEED, adapt_delta = 0.99)

fit$save_object(paste0(output_file, ".RDS"))

#conv_diagnostics <- summarise_draws(fit$draws(), default_convergence_measures())

#max_rhat <- max(conv_diagnostics$rhat, na.rm = TRUE)
#min_ess_bulk <- min(conv_diagnostics$ess_bulk, na.rm = TRUE)
#min_ess_tail <- min(conv_diagnostics$ess_tail, na.rm = TRUE)

#if (min_ess_bulk < 400 || min_ess_tail < 400 || max_rhat > 1.1) {
#    print("low ess or high rhat, refitting with adapt_delta = 0.99, iter_sampling = 4000, max_treedepth = 12")

#    fit <- model$sample(data = dat, adapt_delta = 0.99, iter_sampling = 4000, chains = 4, parallel_chains = 4, seed = SEED)

#}

true_pars_phi <- list(phi = as.numeric(dataset$true_pars$phi))

true_pars_beta <- list(beta = as.numeric(dataset$true_pars$beta))

true_pars_beta_nonzero <- list(beta = true_pars_beta$beta[1:15])

true_pars_r2 <- list(postR2 = as.numeric(dataset$true_pars$postR2))

metrics_phi <- all_metrics(fit, true_pars = true_pars_phi, zero_padding = TRUE) |>
  bind_rows(
    all_metrics(fit, true_pars = true_pars_phi, zero_padding = FALSE)
  ) |>
  mutate(pars = "phi")

metrics_beta <- all_metrics(fit, true_pars = true_pars_beta, zero_padding = TRUE) |>
  bind_rows(
    all_metrics(fit, true_pars = true_pars_beta_nonzero, zero_padding = FALSE)
  ) |>
  mutate(pars = "beta")

metrics_phibeta <- all_metrics(fit, true_pars = c(true_pars_phi, true_pars_beta), zero_padding = TRUE) |>
  bind_rows(
    all_metrics(fit, true_pars = c(true_pars_phi, true_pars_beta_nonzero), zero_padding = FALSE)
  ) |>
  mutate(pars = "phi+beta")

metrics_r2 <- all_metrics(fit, true_pars = true_pars_r2, zero_padding = FALSE) |>
  mutate(pars = "postR2")

out <- bind_rows(metrics_phibeta, metrics_phi, metrics_r2, metrics_beta)  |>
  bind_cols(prior_settings)

write_csv(out, output_file)
