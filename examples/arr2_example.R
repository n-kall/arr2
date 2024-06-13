library(cmdstanr)
library(bayesplot)

arr2_noncentered <- cmdstan_model("arr2_noncentered.stan")


Y <- arima.sim(list(ar = c(0, 0, 0,  0, 0, 0, 0.8), sd = 1), n = 36)

sdata1 <- list(
  T = length(Y),
  Y = Y,
  p = 12,
  mean_R2 = 1/3,
  prec_R2 = 3,
  sigma_sd = 2.5,
  alpha_mean = 0,
  alpha_sd = 10
)

arr2_rhs_noncentered_fit <- arr2_noncentered$sample(
  data = c(sdata1, list(cons = rep(0.1, times = 12))),
  adapt_delta = 0.99
)

arr2_minn_noncentered_fit <- arr2_noncentered$sample(
  data = c(sdata1, list(cons = (12^2 * 0.1) / (1:12)^2)),
  adapt_delta = 0.99
)

mcmc_intervals(arr2_rhs_noncentered_fit$draws("phi"))

mcmc_intervals(arr2_minn_noncentered_fit$draws("phi"))


Y2 <- arima.sim(list(ar = c(0.6, 0.2, 0.1, 0.05), sd = 1), n = 36)

sdata2 <- list(
  T = length(Y2),
  Y = Y2,
  p = 12,
  mean_R2 = 1/3,
  prec_R2 = 3,
  sigma_sd = 2.5,
  alpha_mean = 0,
  alpha_sd = 10
)

arr2_rhs_noncentered_fit2 <- arr2_noncentered$sample(data = c(sdata2, list(cons = rep(0.1, times = 12))))

arr2_minn_noncentered_fit2 <- arr2_noncentered$sample(data = c(sdata2,   list(cons = (12^2 * 0.1) / (1:12)^2)), adapt_delta = 0.99)

mcmc_intervals(arr2_rhs_noncentered_fit2$draws("phi"))

mcmc_intervals(arr2_minn_noncentered_fit2$draws("phi"))
