library(cmdstanr)
library(posterior)
library(tidyverse)

dataset <- readRDS(snakemake@input[["dataset"]])
functions <- snakemake@input[["functions"]]
npars <- as.numeric(snakemake@wildcards[["npars"]])
prior_specs <- snakemake@input[["prior_specs"]]
gq_draws <- snakemake@input[["gq_draws"]]
output_file <- snakemake@output[["fold_summary"]]
SEED <- as.numeric(snakemake@config[["seed"]])

source(functions)
source(prior_specs)



Y <- dataset$data$Y
T <- length(Y)
L <- floor(T / 2) + 1 # initial n_obs for LFO

pred_exact <- posterior::subset_draws(gq_draws, "pred")
variables(pred_exact) <- paste0("pred[", fold, "]")

loglik_exact <- posterior::extract_variable_matrix(gq_draws, "log_lik")


# compute exact LFO-CV elpd
exact_elpd_1sap <- log_mean_exp(loglik_exact)

out <- tibble(
    elpd = exact_elpd_1sap
)

write_csv(out, output_file)
