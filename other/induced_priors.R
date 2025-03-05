library(tidyverse)
library(extraDistr)
library(posterior)
library(furrr)
library(ggdist)
library(directlabels)

papertheme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
    axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black")
  )


arr2_prior <- function(p, mean_R2 = 1/3, prec_R2 = 3, sigma_sd = 1, alpha = rep(0.1, times = p)) {

  sigma <- 1
  
  shapes <- simstudy::betaGetShapes(mean = mean_R2, precision = prec_R2)
  r2 <- rbeta(1, shape1 = shapes[["shape1"]], shape2 = shapes[["shape2"]])

  psi <- extraDistr::rdirichlet(
    1, alpha = alpha
  )

  tau2 <- r2 / (1 - r2)
    
  phi <- rnorm(p, mean = rep(0, times = p), sd = sqrt(tau2 * sigma^2 * psi))
  names(phi) <- paste0("phi[", 1:p, "]")
  
  return(phi)
}

gauss_prior <- function(p, phi_sd = 1) {

  phi <- rnorm(p, sd = phi_sd)
  names(phi) <- paste0("phi[", 1:p, "]")
  
  return(phi)

}

rhs_prior <- function(p, p0 = p / 2, T = 24) {

   tau2 <- abs(rt(1,1)*(p0 / (p - p0) * 1 / sqrt(T)))^2
   lambda2 <- abs(rt(p,1))^2
  lambda <- tau2 * diag(lambda2)
  phi <- MASS::mvrnorm(n = 1, mu = rep(0, p), Sigma = lambda)
  names(phi) <- paste0("phi[", 1:p, "]")

  return(phi)
}

minn_prior <- function(p) {
  n <- 1 / (1:p)^2
    
  tau2 <- rgamma(n = 1, shape = 1, scale = 0.04)
      
  lambda <- tau2 * 2 * diag(n)
  phi <- MASS::mvrnorm(n = 1, mu = rep(0, p), Sigma = lambda)

  names(phi) <- paste0("phi[", 1:p, "]")
  return(phi)
}
  
phi_to_R2 <- function(phi) {
  sigma <- 1
  return(c(R2 = t(phi) %*% phi / (t(phi) %*% phi + sigma^2)))
}

phi_to_nonstationary <- function(phi) {
  return(nonstationary = any(Mod(pracma::roots(c(-1, -phi))) >= 1))
}

phi_to_relR2 <- function(phi) {
  out <- phi^2 / c((t(phi) %*% phi))
  names(out) <- paste0("relR2[", 1:length(phi), "]")
  return(out)
}

phi_to_pacf <- function(phi) {
  out <- INLA::inla.ar.phi2pacf(phi = phi)
  names(out) <- paste0("pacf[", 1:length(phi), "]")
  return(out)
}

phi_to_roots <- function(phi) {
  roots <- pracma::roots(c(1, -phi))
  complex_roots <- roots[Im(roots) != 0]

  moduli <- unique(Mod(roots))
  complex_moduli <- unique(Mod(complex_roots))

  periods <- 2 * pi / unique(abs(Arg(roots)))
  
  out <- c(
    max_modulus = max(moduli),
    max_complex_modulus = max(complex_moduli),
    period_max_modulus = periods[which(moduli == max(complex_moduli))],
    max_period = max(periods),
    modulus_max_period = complex_moduli[which(periods == max(periods))]
  )
  
  return(out)
}


add_induced_priors <- function(draws) {

  
  r2_draws <- draws |>
    as_draws_matrix() |>
    apply(1, phi_to_R2) |>
    as.matrix() |>
    as_draws() |>
    set_variables("R2")

  pacf_draws <- draws |>
    as_draws_matrix() |>
    apply(1, phi_to_pacf) |>
    t() |>
    as_draws()

  rel_r2_draws <- draws |>
    as_draws_matrix() |>
    apply(1, phi_to_relR2) |>
    t() |>
    as_draws()

  roots_draws <- draws |>
    as_draws_matrix() |>
    apply(1, phi_to_roots) |>
    as.matrix() |>
    t() |>
    as_draws()

  nonstationary_draws <- draws |>
    as_draws_matrix() |>
    apply(1, phi_to_nonstationary) |>
    as.matrix() |>
    as_draws() |>
    set_variables("nonstationary")

  bind_draws(draws, r2_draws, rel_r2_draws, pacf_draws, roots_draws, nonstationary_draws)
  
}


set.seed(1994)
niters <- 1e5

minn_prior_draws <- future_map_dfr(1:niters, ~minn_prior(p = 12), .options = furrr_options(seed = 1994)) |>
  as_draws() |>
  add_induced_priors() |>
  as_tibble() |>
  mutate(prior = "Minn.")


rhs_prior_draws <- future_map_dfr(
  1:niters, ~rhs_prior(p = 12),
  .options = furrr_options(seed = 123)) |>
  as_draws() |>
  add_induced_priors() |>
  as_tibble() |>
  mutate(prior = "RHS")


gauss_prior_draws <- future_map_dfr(
  1:niters,
  ~gauss_prior(p = 12),
  .options = furrr_options(seed = 123)
) |>
  as_draws() |>
  add_induced_priors() |>
  as_tibble() |>
  mutate(prior = "Gaussian")


arr2_flat_prior_draws <- future_map_dfr(
  1:niters,
  ~arr2_prior(p = 12),
  .options = furrr_options(seed = 123)
) |>
  as_draws() |>
  add_induced_priors() |>
  as_tibble() |>
  mutate(prior = "ARR2 (flat)")



arr2_bump_prior_draws <- future_map_dfr(
  1:niters,
  ~arr2_prior(p = 12, alpha = dnbinom(1:12, prob = 0.4, size = 4)),
  .options = furrr_options(seed = 123)
) |>
  as_draws() |>
  add_induced_priors() |>
  as_tibble() |>
  mutate(prior = "ARR2 (bump)")

minnesota_cons <- numeric(12)

for (i in 1:length(minnesota_cons)) {
  minnesota_cons[i] <- 12^2 / 10 * 1 / (i^2)
}

arr2_minn_prior_draws <- map_dfr(
  1:niters,
  ~arr2_prior(p = 12, alpha = minnesota_cons),
  .options = furrr_options(seed = 123)
) |>
  as_draws() |>
  add_induced_priors() |>
  as_tibble() |>
  mutate(prior = "ARR2 (Minn.)")

arr2_minn_flatR2_prior_draws <- future_map_dfr(
  1:niters,
  ~arr2_prior(mean_R2 = 0.5, prec_R2 = 2, p = 12, alpha = minnesota_cons),
  .options = furrr_options(seed = 123)
) |>
  as_draws() |>
  add_induced_priors() |>
  as_tibble() |>
  mutate(prior = "ARR2 (Minn., flat R2)")

arr2_minn_bathtubR2_prior_draws <- future_map_dfr(
  1:niters,
  ~arr2_prior(p = 12, mean_R2 = 0.5, prec_R2 = 1, alpha = minnesota_cons),
  .options = furrr_options(seed = 123)
) |>
  as_draws() |>
  add_induced_priors() |>
  as_tibble() |>
  mutate(prior = "ARR2 (Minn., bathtub R2)")


plot_data <- bind_rows(
  minn_prior_draws,
  gauss_prior_draws,
  rhs_prior_draws,
  arr2_minn_prior_draws,
  arr2_flat_prior_draws,
  arr2_bump_prior_draws,
  arr2_minn_flatR2_prior_draws,
  arr2_minn_bathtubR2_prior_draws
)

colours <- c(
  "ARR2 (flat)" = "#4477AA",
  "ARR2 (Minn.)" = "#EE6677",
  "Gauss." = "#228833",
  "Gaussian" = "#228833",
  "Minn." = "#CCBB44",
  "RHS" = "#AA3377",
  "ARR2 (bump)" = "#66CCEE",
  "ARR2 (Minn., flat R2)" = "#FF8C42",
  "ARR2 (Minn., bathtub R2)" = "#FFB3D9"
)


# R2 plot
r2_plot <- plot_data |>
  filter(!(prior %in% c("ARR2 (Minn.)", "ARR2 (bump)", "ARR2 (Minn., flat R2)", "ARR2 (Minn., bathtub R2)"))) |>
  mutate(prior = ifelse(prior == "ARR2 (flat)", "ARR2", prior)) |>
  mutate(prior = ifelse(prior == "Gaussian", "Gauss.", prior)) |>
  ggplot() +
  geom_density(aes(x = R2, colour = prior), bounds = c(0, 1)) +
    papertheme +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_colour_manual(values = colours, guide = "none") +
  xlab("$R^2$")

r2_plot <- direct.label(r2_plot)




library(patchwork)


# max modulus for stationary check
stationarity_plot <- plot_data |>
   ggplot() +
   geom_vline(xintercept = 1, linetype = "dotted") +
  geom_density(aes(x = max_modulus, colour = prior)) +
  xlim(0, 2.2) +
  xlab("Max modulus") +
  scale_colour_manual(values = colours, labels = c(Gaussian = "Gauss."), guide = "none") +
  papertheme +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank()
  )

stationarity_plot <- direct.label(stationarity_plot, method = list(dl.trans(y = y + 0.1, x = x + 0.8), "top.points", "bumpup", cex = 0.8))

bayesflow::save_tikz_plot(stationarity_plot, "stationarity_plot.tex", width = 4.5, height = 3)



  

# phi
phi_plot <- plot_data |>
  select(starts_with("phi"), "prior") |>
  pivot_longer(cols = starts_with("phi")) |>
  mutate(lag = as.numeric(str_extract(name, "\\d+")))

# pacf

pacf_plot <- plot_data |>
  filter(nonstationary == 0) |>
  select(starts_with("pacf"), "prior") |>
  pivot_longer(cols = starts_with("pacf")) |>
  mutate(lag = as.numeric(str_extract(name, "\\d+"))) |>
  filter(prior != "Gaussian", abs(value) <= 1)


# relative R2
relr2_plot <- plot_data |>
  select(starts_with("relR2"), "prior") |>
  pivot_longer(cols = starts_with("relR2")) |>
  mutate(lag = as.numeric(str_extract(name, "(?<=\\[)\\d+(?=\\])")))



relr2_plot_out <- relr2_plot |>
  filter(!(prior %in% c("ARR2 (Minn., flat R2)", "ARR2 (Minn., bathtub R2)"))) |>
  mutate(prior = ifelse(prior == "Gaussian", "Gauss.", prior)) |>
  group_by(prior, lag) |>
  summarise(mean = mean(value), median = median(value), q5 = quantile2(value, 0.05), q95 = quantile2(value, 0.95), sd = sd(value)) |>
  ggplot(aes(x = lag, colour = prior)) +
  #  geom_pointinterval(aes(y = mean, ymin = q5, ymax = q95)) +
  geom_point(aes(y = mean)) +
  geom_line(aes(y = mean)) + 
#  stat_interval(aes(y = value), .width = c(0.25, 0.5, 0.95)) +
  scale_colour_manual(values = colours) +
  guides(colour = "none") +
  scale_x_continuous(breaks = seq(1, 12, by = 1)) +
  xlab("Lag") +
  papertheme +
theme(
  strip.placement = "outside",
  axis.title.y = element_blank()
)

relr2_plot_out <- direct.label(relr2_plot_out, method = list(dl.trans(y = y + 0.1, x = x + 0.8), "top.points", "bumpup"))

r2_paper_plot <- r2_plot / relr2_plot_out

bayesflow::save_tikz_plot(r2_paper_plot, "r2_plot_priors.tex", width = 5.5, height = 5.5)


all_plot_data <- bind_rows(
  relr2_plot,
  pacf_plot,
  phi_plot
) |>
  mutate(quantity = str_remove_all(name, "[^a-zA-Z]"))

colours <- c(
  "ARR2 (flat)" = "#4477AA",
  "ARR2 (Minn.)" = "#EE6677",
  "Gauss." = "#228833",
  "Gaussian" = "#228833",
  "Minn." = "#CCBB44",
  "RHS" = "#AA3377",
  "ARR2 (bump)" = "#66CCEE",
  "ARR2 (Minn., flat R2)" = "#FF8C42",
  "ARR2 (Minn., bathtub R2)" = "#FFB3D9"
)


p_quants <- all_plot_data |>
  filter(!(prior %in% c("ARR2 (bump)", "ARR2 (Minn., flat R2)", "ARR2 (Minn., bathtub R2)"))) |>
  filter(abs(value) < 100) |>
  group_by(prior, lag, quantity) |>
  summarise(mean = mean(value), median = median(value), q5 = quantile2(value, 0.05), q95 = quantile2(value, 0.95), sd = sd(value)) |>
  ggplot(aes(x = lag, colour = prior)) +
  geom_pointinterval(aes(y = mean, ymin = q5, ymax = q95)) +
#  stat_interval(aes(y = value), .width = c(0.25, 0.5, 0.95)) +
  facet_grid(factor(quantity, levels = c("phi", "pacf", "relR")) ~ prior, scales = "free", labeller = as_labeller(c(
     phi = "AR coefficient", pacf = "Partial autocorrelation", relR = "Relative $R^2$",
    Gaussian = "Gauss.", RHS = "RHS", "Minn." = "Minn.", "ARR2 (flat)" = "ARR2 (flat)", "ARR2 (Minn.)" = "ARR2 (Minn.)")
),     switch = "y") +
  scale_colour_manual(values = colours) +
  guides(colour = "none") +
  scale_x_continuous(breaks = seq(1, 12, by = 2)) +
  xlab("Lag") +
  papertheme +
theme(
  strip.placement = "outside",
  axis.title.y = element_blank()
)


bayesflow::save_tikz_plot(p_quants, "induced_priors.tex", width = 5.5, height = 5.5)




