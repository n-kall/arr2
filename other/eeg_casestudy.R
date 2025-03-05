library(tidyverse)
library(cmdstanr)
library(posterior)
library(ggdist)
library(lemon)

source("revision/induced_priors_functions.R")

eeg_data <- read_csv("revision/eeg400.csv")

standata <- list(
  Y = eeg_data$EEG,
  T = length(eeg_data$EEG),
  p = 30,
  mean_R2 = 1/3,
  prec_R2 = 3,
  sigma_sd = 2.5 * sd(eeg_data$EEG),
  alpha_mean = mean(eeg_data$EEG),
  alpha_sd = 2.5 * sd(eeg_data$EEG),
  hs_df = 1,
  hs_df_global = 1,
  hs_df_slab = 4,
  hs_scale_slab = 2,
  p0 = 15,
  alpha_cent = 1,
  phi_cent = 1,
  sigma_cent = 1,
  phi_sd = 2.5
)

ndraws <- 1000

arr2_flat_model <- cmdstan_model("snakemake/stan/arr2_gamma_ncp_flat.stan")
arr2_flat_fit <- arr2_flat_model$sample(data = standata, iter_sampling = ndraws, seed = 1994)

arr2_flat_draws <- add_induced_priors(arr2_flat_fit$draws()) |>
  as_draws_df() |>
  as_tibble() |>
    bind_cols(arr2_flat_fit$draws("postR2") |> as_draws_df() |> as_tibble() |> select("postR2")) |>
  mutate(prior = "ARR2 (flat)")


arr2_minn_model <- cmdstan_model("snakemake/stan/arr2_gamma_ncp_minnesota.stan")
arr2_minn_fit <- arr2_minn_model$sample(data = standata, iter_sampling = ndraws, seed = 1994)

arr2_minn_draws <- add_induced_priors(arr2_minn_fit$draws()) |>
  as_draws_df() |>
  as_tibble() |>
  bind_cols(arr2_minn_fit$draws("postR2") |> as_draws_df() |> as_tibble() |> select("postR2")) |>
  mutate(prior = "ARR2 (Minn.)")


minn_model <- cmdstan_model("snakemake/stan/minnesota.stan")
minn_fit <- minn_model$sample(data = standata, iter_sampling = ndraws, seed = 1994)

minn_draws <- add_induced_priors(minn_fit$draws()) |>
  as_draws_df() |>
  as_tibble() |>
  bind_cols(minn_fit$draws("postR2") |> as_draws_df() |> as_tibble() |> select("postR2")) |>
  mutate(prior = "Minnesota")


rhs_model <- cmdstan_model("snakemake/stan/rhs.stan")
rhs_fit <- rhs_model$sample(data = standata, adapt_delta = 0.999, iter_sampling = ndraws, seed = 1994)

rhs_draws <- add_induced_priors(rhs_fit$draws()) |>
  as_draws_df() |>
  as_tibble() |>
  bind_cols(rhs_fit$draws("postR2") |> as_draws_df() |> as_tibble() |> select("postR2")) |>
  mutate(prior = "RHS")


gaussian_model <- cmdstan_model("snakemake/stan/gaussian.stan")
gaussian_fit <- gaussian_model$sample(data = standata, iter_sampling = ndraws, seed = 1994)

gaussian_draws <- add_induced_priors(gaussian_fit$draws()) |>
  as_draws_df() |>
  as_tibble() |>
  bind_cols(gaussian_fit$draws("postR2") |> as_draws_df() |> as_tibble() |> select("postR2")) |>
  mutate(prior = "Gaussian")


all_draws <- bind_rows(
  gaussian_draws,
  rhs_draws,
  minn_draws,
  arr2_minn_draws,
  arr2_flat_draws
) |>
  mutate(prior = factor(prior, levels = rev(unique(prior)))) 


papertheme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
    axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black")
  )

colours <- c("ARR2 (flat)" = "#4477AA", "ARR2 (Minn.)" = "#EE6677", "Gaussian" = "#228833", "Minnesota" = "#CCBB44", "RHS" = "#AA3377", ARR2_bump = "#66CCEE")


r2_draws <- arr2_minn_fit$draws("postR2") |>
  as_draws_df() |>
  as_tibble()


# plot coefficients

p3 <- all_draws |>
  select(prior, starts_with("phi")) |>
  filter(prior != "Gaussian") |>
  pivot_longer(cols = starts_with("phi")) |>
  mutate(lag = as.numeric(str_extract(name, "\\d+"))) |>
  ggplot(aes(x = lag)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  stat_slabinterval(aes(y = value, fill = prior, colour = prior), slab_alpha = 0, linewidth = 4, point_size = 1.5, .width = c(0.95)) +
  guides(fill = "none") +
  ylab("$\\phi$") +
  xlab("Lag") +
  facet_wrap(~prior, ncol = 1, labeller = as_labeller(c(Minnesota = "Minn.", RHS = "RHS", "ARR2 (Minn.)" = "ARR2 (Minn.)", "ARR2 (flat)" = "ARR2 (flat)"))) +
  guides(colour = "none") +
  papertheme


# plot max period

p1_period <- all_draws |>
  filter(prior != "Gaussian") |>
  ggplot(aes(x = period_max_modulus, fill = prior)) +
#  stat_slabinterval(.width = c(0.95)) +
  stat_histinterval(.width = c(0.9)) +
  guides(fill = "none") +
  ylab("AR coefficient") +
  xlab("Corresponding\nperiod") +
  facet_rep_wrap(~prior, ncol = 1) +
  scale_x_continuous(breaks = c(13, 14)) +
  coord_cartesian(xlim = c(12.5, 14.5)) +
  papertheme


# plot max modulus

p1_mod <- all_draws |>
  filter(prior != "Gaussian") |>
  ggplot(aes(x = max_complex_modulus, fill = prior)) +
#  stat_slabinterval(.width = c(0.95)) +
  stat_histinterval(.width = c(0.9)) +
  guides(fill = "none") +
  ylab("Prior") +
  xlab("Max\nmodulus") +
  facet_rep_wrap(~prior, ncol = 1) +
  papertheme +
scale_x_continuous(breaks = c(0.95, 1), limits = c(0.94, 1.01))


# r2 plot
p1_r2 <- all_draws |>
  filter(prior != "Gaussian") |>
  ggplot(aes(x = postR2, fill = prior)) +
  #  stat_slabinterval(.width = c(0.95)) +
  stat_histinterval(.width = c(0.9)) +
  guides(fill = "none") +
  ylab("Prior") +
  xlab("$R^2$") +
  facet_rep_wrap(~prior, ncol = 1) +
  papertheme +
  scale_x_continuous(breaks = c(0.4, 0.5, 0.6))


library(patchwork)

p_all <- p3 + (p1_r2 + p1_mod + p1_period & theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), strip.text = element_blank())) &   scale_colour_manual(values = colours) &   scale_fill_manual(values = colours)


bayesflow::save_tikz_plot(p_all, "eeg_plot.tex", width = 6, height = 4)



# raw data plot

eeg_raw <- eeg_data |>
  mutate(time = 1:400) |>
  filter(time <= 100) |>
  ggplot(aes(x = time, y = EEG)) +
  geom_line() +
  papertheme +
  xlab("Time (ms)") +
  ylab("EEG signal")


bayesflow::save_tikz_plot(eeg_raw, "eeg_raw_plot.tex", width = 5, height = 1.8)


# plot coefficients

p3_gauss <- all_draws |>
  select(prior, starts_with("phi")) |>
  filter(prior == "Gaussian") |>
  pivot_longer(cols = starts_with("phi")) |>
  mutate(lag = as.numeric(str_extract(name, "\\d+"))) |>
  ggplot(aes(x = lag)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_slabinterval(aes(y = value, fill = prior, colour = prior), size = 0.5, slab_alpha = 0, .width = c(0.95)) +
  guides(fill = "none") +
  ylab("$\\phi$") +
  xlab("Lag") +
  facet_wrap(~prior, ncol = 1) +
  guides(colour = "none") +
  papertheme


# plot max period

p1_period_gauss <- all_draws |>
  filter(prior == "Gaussian") |>
  ggplot(aes(x = period_max_modulus, fill = prior)) +
  stat_halfeye() +
  guides(fill = "none") +
  ylab("$\\phi$") +
  xlab("Corresponding period") +
  facet_rep_wrap(~prior, ncol = 1) +
  coord_cartesian(xlim = c(0, 35)) +
  papertheme


# plot max modulus

p1_mod_gauss <- all_draws |>
  filter(prior == "Gaussian") |>
  ggplot(aes(x = max_complex_modulus, fill = prior)) +
  stat_halfeye() +
  guides(fill = "none") +
  ylab("Prior") +
  xlab("Largest complex modulus") +
  facet_rep_wrap(~prior, ncol = 1) +
  papertheme


# r2 plot
p1_r2_gauss <- all_draws |>
  filter(prior == "Gaussian") |>
  ggplot(aes(x = postR2, fill = prior)) +
  stat_slabinterval(.width = c(0.95)) +
  guides(fill = "none") +
  ylab("Prior") +
  xlab("$R^2$") +
  facet_rep_wrap(~prior, ncol = 1) +
  papertheme 


p_gauss <- p3_gauss + (p1_r2_gauss + p1_mod_gauss + p1_period_gauss & theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), strip.text = element_blank()))  & scale_colour_manual(values = colours) & scale_fill_manual(values = colours)



layout <- "A
A
A
A
B"
