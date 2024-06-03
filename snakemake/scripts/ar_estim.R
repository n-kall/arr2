library(cmdstanr)
library(posterior)
library(tidyverse)

dataset <- readRDS(snakemake@input[["dataset"]])
functions <- snakemake@input[["functions"]]
npars <- as.numeric(snakemake@wildcards[["npars"]])
prior_specs <- snakemake@input[["prior_specs"]]
exe_file <- snakemake@input[["exe_file"]]
output_file <- snakemake@output[["fit_summary"]]
SEED <- as.numeric(snakemake@config[["seed"]])

source(functions)
source(prior_specs)

prior_settings <- c(
  all_prior, gaussian_prior, arr2_prior, rhs_pars
)

dat <- c(prior_settings, dataset$data)
dat$T <- length(dat$Y)

# scale npars depending on observation size
npars <- (dat$T / 24) * npars

dat$p <- npars

# set RHS p0
dat$p0 <- floor(dat$p / 2)

# load model
model <- cmdstan_model(exe_file = exe_file)

# fit model
fit <- model$sample(data = dat, seed = SEED, parallel_chains = 4, chains = 4)

# diagnostics
conv_diagnostics <- summarise_draws(fit$draws(), default_convergence_measures())

max_rhat <- max(conv_diagnostics$rhat, na.rm = TRUE)
min_ess_bulk <- min(conv_diagnostics$ess_bulk, na.rm = TRUE)
min_ess_tail <- min(conv_diagnostics$ess_tail, na.rm = TRUE)

if (min_ess_bulk < 400 || min_ess_tail < 400 || max_rhat > 1.1) {
  print("low ess or high rhat, refitting with adapt_delta = 0.99, iter_sampling = 4000")

  fit <- model$sample(data = dat, adapt_delta = 0.99, iter_sampling = 4000, parallel_chains = 4, seed = SEED)

}

# check metrics
true_pars_phi <- list(phi = as.numeric(dataset$true_pars$phi))
true_pars_r2 <- list(postR2 = as.numeric(dataset$true_pars$postR2))

metrics_phi <- all_metrics(fit, true_pars = true_pars_phi, zero_padding = TRUE)

if (!any(is.null(true_pars_phi$phi))) {

  metrics_phi <- metrics_phi |>
    bind_rows(
      all_metrics(fit, true_pars = true_pars_phi, zero_padding = FALSE)
    )
}

metrics_phi <- metrics_phi |>
  mutate(pars = "phi")

metrics_r2 <- all_metrics(fit, true_pars = true_pars_r2, zero_padding = FALSE)

metrics_r2 <- metrics_r2 |>
  mutate(pars = "postR2")

out <- bind_rows(metrics_phi, metrics_r2) |>
  bind_cols(prior_settings)


write_csv(out, output_file)
