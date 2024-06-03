library(cmdstanr)
library(resample)
library(posterior)
library(scoringRules)
library(tidyverse)

fit <- readRDS(snakemake@input[["fit_object"]])
dataset <- readRDS(snakemake@input[["dataset"]])
output_file <- snakemake@output[["fit_summary"]]
functions <- snakemake@input[["functions"]]
fold <- as.numeric(snakemake@wildcards[["fold"]])

source(functions)

# forecast metrics
log_lik <- fit$draws(variable = c("yf_lpdf"), format = "matrix")
elpd <- log_mean_exp(log_lik)

R2_summ <- summarise_draws(fit$draws(variable = c("postr22")))
R2_mean <- R2_summ$mean
R2_sd <- R2_summ$sd

R2_reg_summ <- summarise_draws(fit$draws(variable = c("R2_reg")))
R2_reg_mean <- R2_reg_summ$mean
R2_reg_sd <- R2_reg_summ$sd

R2_trend_summ <- summarise_draws(fit$draws(variable = c("R2_trend")))
R2_trend_mean <- R2_trend_summ$mean
R2_trend_sd <- R2_trend_summ$sd


drw <- subset_draws(fit$draws(), setdiff(variables(fit$draws()), "z_tau1"))
s <- summarise_draws(drw, rhat, ess_tail, ess_bulk)

#fit_sum <- fit$summary()
min_ess_bulk <- min(s[["ess_bulk"]], na.rm = TRUE)
min_ess_tail <- min(s[["ess_tail"]], na.rm = TRUE)
max_rhat <- max(s[["rhat"]], na.rm = TRUE)
divergences <- sum(fit$sampler_diagnostics()[,,2])
res <- data.frame(
  elpd = elpd,
  R2_mean = R2_mean,
  R2_reg_mean = R2_reg_mean,
  R2_trend_mean = R2_trend_mean,
  min_ess_bulk,
  min_ess_tail,
  max_rhat,
  divergences
)

write_csv(res, output_file)
