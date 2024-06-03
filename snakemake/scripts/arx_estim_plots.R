library(tidyverse)
library(patchwork)
library(vroom)
library(directlabels)
library(lemon)

if (exists("snakemake")) {
  summary_path <- snakemake@config[["summary_path"]]
  results_file <- snakemake@input[["results"]]
} else {
  summary_path <- "snakemake/summary/"
  results_file <- "snakemake/summary/arx_estim_results.csv"
}

colours <- c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#AA3377")
colours <- c("#4477AA", "#EE6677", "#CCBB44", "#AA3377")

dataset_labeller <- c(
  "0" = "Uncorrelated $X$\n($\\rho = 0$)",
  "0.5" = "Moderately correlated $X$\n($\\rho = 0.5$)",
  "0.9" = "Highly correlated $X$\n($\\rho = 0.9$)"
)

model_labels <- c(
  arxr2_flat = "ARR2 (simplex, flat)",
  arxr2_gamma_flat = "ARR2 (gamma, flat)",
  arr2_arx_gamma_minnesota = "ARR2 (Minn.)",
  arr2_arx_gamma_flat = "ARR2 (flat)",
  arr2_arx_gamma = "ARR2 (gamma)",
  arxr2_minnesota = "ARR2 (simplex, minnesota)",
  arxr2_gamma_minnesota = "ARXR2 (gamma, minnesota)",
  arxr2_dc_flat = "ARR2 (dc, flat)",
  arxr2_dc_minnesota = "ARR2 (dc, minnesota)",
  gaussian_arx = "Gauss.",
  minnesota_arx =  "Minn.",
  rhs_arx = "RHS"
)

papertheme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
    axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black")
  )

r2_plot <- function(dat, dodge_width = 0) {
  dat |>
    filter(
      variable == "postR2",
      pars == "postR2",
      dataset == "cement",
      zero_padding == FALSE,
      !(is.na(mean))
    ) |>
    group_by(npars_x, model, dataset, rho, n_obs) |>
    summarise(
      mean_post_R2 = mean(mean),
      q5_post_R2 = mean(q5),
      q95_post_R2 = mean(q95),
      mean_true_R2 = mean(true_par)
    ) |>
    mutate(model_label = model_labels[model]) |>
    ggplot(aes(x = npars_x, color = model, group = model)) +
    geom_dl(aes(y = mean_post_R2, label = model_label), method = list(dl.trans(x = x + 0.2), "last.points", "bumpup", cex = 0.8)) +
    geom_pointrange(
      aes(y = mean_post_R2,
          ymin = q5_post_R2,
          ymax = q95_post_R2,
          group = model,
          colour = model
          ), position = position_dodge2(width = dodge_width)) +
    geom_line(aes(x = npars_x, y = mean_true_R2), inherit.aes = FALSE, colour = "black") +
    scale_color_manual(
      name = "Prior",
      labels = model_labels,
      guide = "none",
      values = colours
    ) +
    facet_rep_wrap(rho~dataset, ncol = 2, scales = "free_y") +
    scale_x_continuous(expand = expansion(mult = c(0.06, 0.4)), breaks = c(20, 200, 400)) +
    xlab("AR order")+
    ylab("Posterior $R^2$")

}

metric_plot <- function(dat, nobs, rho, zero_padding, metric, metric_name, pars, direct_labels = FALSE, xbreaks = waiver(), ybreaks = waiver(), ylim = waiver(), xlim = waiver()) {

  p <- dat |>
    filter(
      variable == "phi[1]",
      n_obs == !!nobs,
      zero_padding == !!zero_padding,
      pars == !!pars,
      rho == !!rho
    ) |>
    group_by(npars_x, model, dataset, rho) |>
    summarise(mean_metric = mean(.data[[metric]]), stdev = sd(.data[[metric]]), nreps = n())  |>
    mutate(model_label = model_labels[[model]])  |>
    ggplot(aes(x = npars_x, color = model, group = model)) +
    geom_pointrange(aes(
      y = mean_metric,
      ymin = mean_metric - stdev / sqrt(nreps),
      ymax = mean_metric + stdev / sqrt(nreps)
    ),
    position = position_dodge(width = 15),
    size = 0.001
    ) +
    geom_line(aes(
      y = mean_metric),
      position = position_dodge(width = 15)
      ) +
    facet_rep_wrap(~rho, ncol = 3, labeller = as_labeller(dataset_labeller), repeat.tick.labels = TRUE) +
    scale_color_manual(
      name = "Prior",
      labels = model_labels,
      guide = "none",
      values = colours
    ) +
    ylab(metric_name) +
    scale_y_continuous(breaks = ybreaks) + 
    scale_x_continuous(breaks = xbreaks)

  if (direct_labels) {

    p <- p + geom_dl(
      aes(y = mean_metric, label = model_label),
      method = list(dl.trans(x = x + 0.2), "last.points", "bumpup",  cex = 0.8)) +
      scale_x_continuous(expand = expansion(mult = c(0.06, 0.4)), breaks = xbreaks)

  }

  p
}


dat_tidy <- read_csv(results_file) |>
  mutate(model = str_remove(model, "results/arx_estim/")) |>
  select(model, dataset, n_obs, npars_x, rho, rep, variable, rmse, mae, energy_score, pars, zero_padding)



# metric comparison

(dat_tidy |>
  filter(variable == "postR2", pars == "beta", npars_x == 400) |>
  pivot_longer(cols = c("rmse", "energy_score", "mae")) |>
  ggplot(aes(y = value, x = model, group = rep)) +
   #  stat_summary() +
   geom_point() +
   geom_line() +
  facet_wrap(zero_padding~name, scales = "free_y") + ggtitle("beta")) /
(dat_tidy |>
  filter(variable == "postR2", pars == "phi", npars_x == 400) |>
  pivot_longer(cols = c("rmse", "energy_score", "mae")) |>
  ggplot(aes(y = value, x = model, group = rep)) +
   #  stat_summary() +
   geom_point() +
   geom_line() + 
  facet_wrap(zero_padding~name, scales = "free_y") + ggtitle("phi"))



rmse_120_beta_rho0_padded <- metric_plot(
  dat_tidy |> filter(dataset == "cement", model != "gaussian_arx"),
  nobs = 120,
  zero_padding = TRUE,
  rho = 0,
  "rmse",
  "RMSE",
  "beta",
  xbreaks = c(20, 100, 200, 400)
) +
  papertheme + ggtitle(NULL, subtitle = "$\\beta$") +
theme(axis.title.x = element_blank()) +
coord_capped_cart(bottom = "right", gap = 0)

rmse_120_beta_rho05_padded <- metric_plot(
  dat_tidy |> filter(dataset == "cement", model != "gaussian_arx"),
  nobs = 120,
  zero_padding = TRUE,
  rho = 0.5,
  "rmse",
  "RMSE",
  "beta",
  xbreaks = c(20, 100, 200, 400)
) +
  papertheme +
coord_capped_cart(bottom = "right", gap = 0)

gc()

rmse_120_beta_rho09_padded <- metric_plot(
  dat_tidy |> filter(dataset == "cement", model != "gaussian_arx"),
  nobs = 120,
  zero_padding = TRUE,
  rho = 0.9,
  metric = "rmse",
  metric_name = "RMSE",
  pars = "beta",
  direct_labels = TRUE,
  xbreaks = c(20, 100, 200, 400)
) + papertheme +
theme(axis.title.x = element_blank()) + 
  coord_capped_cart(bottom = "right", gap = 0)


rmse_120_beta_padded <- rmse_120_beta_rho0_padded + rmse_120_beta_rho05_padded + rmse_120_beta_rho09_padded + plot_layout(axis_titles = "collect", design = ("AAAABBBBCCCCC")) & labs(x = "No. of exogen.\npredictors")

gc()

rmse_120_phi_rho0_padded <- metric_plot(
  dat_tidy |> filter(dataset == "cement", model != "gaussian_arx"),
  nobs = 120,
  zero_padding = TRUE,
  rho = 0,
  "rmse",
  "RMSE",
  "phi",
  xbreaks = c(20, 100, 200, 400)
) +
  papertheme + ggtitle(NULL, subtitle = "$\\phi$") +
  theme(axis.title.x = element_blank()) + 
coord_capped_cart(bottom = "right", gap = 0)

gc()

rmse_120_phi_rho05_padded <- metric_plot(
  dat_tidy |> filter(dataset == "cement", model != "gaussian_arx"),
  nobs = 120,
  zero_padding = TRUE,
  rho = 0.5,
  "rmse",
  "RMSE",
  "phi",
  xbreaks = c(20, 100, 200, 400)
) +
  papertheme +
coord_capped_cart(bottom = "right", gap = 0)

gc()

rmse_120_phi_rho09_padded <- metric_plot(
  dat_tidy |> filter(dataset == "cement", model != "gaussian_arx"),
  nobs = 120,
  zero_padding = TRUE,
  rho = 0.9,
  metric = "rmse",
  metric_name = "RMSE",
  pars = "phi",
  direct_labels = TRUE,
  xbreaks = c(20, 100, 200, 400)
) + papertheme +
theme(axis.title.x = element_blank()) + 
coord_capped_cart(bottom = "right", gap = 0)

rmse_120_phi_padded <- rmse_120_phi_rho0_padded + rmse_120_phi_rho05_padded + rmse_120_phi_rho09_padded + plot_layout(axis_titles = "collect", design = ("AAAABBBBCCCCC")) & labs(x = "No. of exogen.\npredictors")


ggsave(paste0(summary_path, "arx_estim_rmse_beta.pdf"), rmse_120_beta_padded, width = 5.5, height = 2.3)

bayesflow::save_tikz_plot(rmse_120_beta_padded, paste0(summary_path, "arx_estim_rmse_beta.tex"), width = 5.5, height = 2.3)



ggsave(paste0(summary_path, "arx_estim_rmse_phi.pdf"), rmse_120_phi_padded, width = 5.5, height = 2.3)

bayesflow::save_tikz_plot(rmse_120_phi_padded, paste0(summary_path, "arx_estim_rmse_phi.tex"), width = 5.5, height = 2.3)



metric_plot_indiv <- function(dat, nobs, zero_padding, metric, metric_name, pars) {

  dat |>
    filter(
      variable == "phi[1]",
      nobs == !!nobs,
      zero_padding == !!zero_padding,
      pars == !!pars
    ) |>
    mutate(model_label = if_else(rho == 0.9, model_labels[model], NA)) |>
    ggplot(aes(y = .data[[metric]], x = npars_x, color = model, group = model)) +
    geom_point(position = position_dodge(width = 15)) +
    geom_dl(aes(label = model_label), method = list(dl.trans(x = x + 0.2), "last.points", "bumpup", cex = 0.8)) +
    facet_rep_wrap(~rho, ncol = 3, labeller = as_labeller(dataset_labeller), repeat.tick.labels = TRUE) +
    scale_color_manual(
      name = "Prior",
      labels = model_labels,
      guide = "none",
      values = colours
    ) +
    ylab(metric_name) +
    xlab("No. of exogenous predictors") +
    scale_x_continuous(breaks = c(20, 100, 200, 400), expand = expansion(mult = c(0.06, 0.2))) +
    scale_y_log10()

}


rmse_120_beta_rho0_unpadded <- metric_plot(
  dat_tidy |> filter(dataset == "cement", model != "gaussian_arx"),
  nobs = 120,
  zero_padding = FALSE,
  rho = 0,
  "rmse",
  "RMSE",
  "beta",
  xbreaks = c(20, 100, 200, 400)
) +
  papertheme + ggtitle(NULL, subtitle = "$\\beta$") +
theme(axis.title.x = element_blank()) +
coord_capped_cart(bottom = "right", gap = 0)

rmse_120_beta_rho05_unpadded <- metric_plot(
  dat_tidy |> filter(dataset == "cement", model != "gaussian_arx"),
  nobs = 120,
  zero_padding = FALSE,
  rho = 0.5,
  "rmse",
  "RMSE",
  "beta",
  xbreaks = c(20, 100, 200, 400)
) +
  papertheme +
coord_capped_cart(bottom = "right", gap = 0)

gc()

rmse_120_beta_rho09_unpadded <- metric_plot(
  dat_tidy |> filter(dataset == "cement", model != "gaussian_arx", zero_padding == FALSE),
  nobs = 120,
  zero_padding = FALSE,
  rho = 0.9,
  metric = "rmse",
  metric_name = "RMSE",
  pars = "beta",
  direct_labels = TRUE,
  xbreaks = c(20, 100, 200, 400)
) + papertheme +
theme(axis.title.x = element_blank()) + 
  coord_capped_cart(bottom = "right", gap = 0)


rmse_120_beta_unpadded <- rmse_120_beta_rho0_unpadded + rmse_120_beta_rho05_unpadded + rmse_120_beta_rho09_unpadded + plot_layout(axis_titles = "collect", design = ("AAAABBBBCCCCC")) & labs(x = "No. of exogen.\npredictors")

gc()

ggsave(paste0(summary_path, "arx_estim_rmse_beta_nonzero.pdf"), rmse_120_beta_unpadded, width = 5.5, height = 2.3)

bayesflow::save_tikz_plot(rmse_120_beta_unpadded, paste0(summary_path, "arx_estim_rmse_beta_nonzero.tex"), width = 5.5, height = 2.3)






# energy

energy_score_120_beta_rho0_padded <- metric_plot(
  dat_tidy |> filter(dataset == "cement", model != "gaussian_arx"),
  nobs = 120,
  zero_padding = TRUE,
  rho = 0,
  "energy_score",
  "ENERGY_SCORE",
  "beta",
  xbreaks = c(20, 100, 200, 400)
) +
  papertheme + ggtitle(NULL, subtitle = "$\\beta$") +
theme(axis.title.x = element_blank()) +
coord_capped_cart(bottom = "right", gap = 0)

energy_score_120_beta_rho05_padded <- metric_plot(
  dat_tidy |> filter(dataset == "cement", model != "gaussian_arx"),
  nobs = 120,
  zero_padding = TRUE,
  rho = 0.5,
  "energy_score",
  "ENERGY_SCORE",
  "beta",
  xbreaks = c(20, 100, 200, 400)
) +
  papertheme +
coord_capped_cart(bottom = "right", gap = 0)

gc()

energy_score_120_beta_rho09_padded <- metric_plot(
  dat_tidy |> filter(dataset == "cement", model != "gaussian_arx"),
  nobs = 120,
  zero_padding = TRUE,
  rho = 0.9,
  metric = "energy_score",
  metric_name = "ENERGY_SCORE",
  pars = "beta",
  direct_labels = TRUE,
  xbreaks = c(20, 100, 200, 400)
) + papertheme +
theme(axis.title.x = element_blank()) + 
  coord_capped_cart(bottom = "right", gap = 0)


energy_score_120_beta_padded <- energy_score_120_beta_rho0_padded + energy_score_120_beta_rho05_padded + energy_score_120_beta_rho09_padded + plot_layout(axis_titles = "collect", design = ("AAAABBBBCCCCC")) & labs(x = "No. of exogen.\npredictors")

gc()

energy_score_120_phi_rho0_padded <- metric_plot(
  dat_tidy |> filter(dataset == "cement", model != "gaussian_arx"),
  nobs = 120,
  zero_padding = TRUE,
  rho = 0,
  "energy_score",
  "ENERGY_SCORE",
  "phi",
  xbreaks = c(20, 100, 200, 400)
) +
  papertheme + ggtitle(NULL, subtitle = "$\\phi$") +
  theme(axis.title.x = element_blank()) + 
coord_capped_cart(bottom = "right", gap = 0)

gc()

energy_score_120_phi_rho05_padded <- metric_plot(
  dat_tidy |> filter(dataset == "cement", model != "gaussian_arx"),
  nobs = 120,
  zero_padding = TRUE,
  rho = 0.5,
  "energy_score",
  "ENERGY_SCORE",
  "phi",
  xbreaks = c(20, 100, 200, 400)
) +
  papertheme +
coord_capped_cart(bottom = "right", gap = 0)

gc()

energy_score_120_phi_rho09_padded <- metric_plot(
  dat_tidy |> filter(dataset == "cement", model != "gaussian_arx"),
  nobs = 120,
  zero_padding = TRUE,
  rho = 0.9,
  metric = "energy_score",
  metric_name = "ENERGY_SCORE",
  pars = "phi",
  direct_labels = TRUE,
  xbreaks = c(20, 100, 200, 400)
) + papertheme +
theme(axis.title.x = element_blank()) + 
coord_capped_cart(bottom = "right", gap = 0)

energy_score_120_phi_padded <- energy_score_120_phi_rho0_padded + energy_score_120_phi_rho05_padded + energy_score_120_phi_rho09_padded + plot_layout(axis_titles = "collect", design = ("AAAABBBBCCCCC")) & labs(x = "No. of exogen.\npredictors")


ggsave(paste0(summary_path, "arx_estim_energy_score_beta.pdf"), energy_score_120_beta_padded, width = 5.5, height = 2.3)

bayesflow::save_tikz_plot(energy_score_120_beta_padded, paste0(summary_path, "arx_estim_energy_score_beta.tex"), width = 5.5, height = 2.3)



ggsave(paste0(summary_path, "arx_estim_energy_score_phi.pdf"), energy_score_120_phi_padded, width = 5.5, height = 2.3)

bayesflow::save_tikz_plot(energy_score_120_phi_padded, paste0(summary_path, "arx_estim_energy_score_phi.tex"), width = 5.5, height = 2.3)



energy_score_120_beta_rho0_unpadded <- metric_plot(
  dat_tidy |> filter(dataset == "cement", model != "gaussian_arx"),
  nobs = 120,
  zero_padding = FALSE,
  rho = 0,
  "energy_score",
  "ENERGY_SCORE",
  "beta",
  xbreaks = c(20, 100, 200, 400)
) +
  papertheme + ggtitle(NULL, subtitle = "$\\beta$") +
theme(axis.title.x = element_blank()) +
coord_capped_cart(bottom = "right", gap = 0)

energy_score_120_beta_rho05_unpadded <- metric_plot(
  dat_tidy |> filter(dataset == "cement", model != "gaussian_arx"),
  nobs = 120,
  zero_padding = FALSE,
  rho = 0.5,
  "energy_score",
  "ENERGY_SCORE",
  "beta",
  xbreaks = c(20, 100, 200, 400)
) +
  papertheme +
coord_capped_cart(bottom = "right", gap = 0)

gc()

energy_score_120_beta_rho09_unpadded <- metric_plot(
  dat_tidy |> filter(dataset == "cement", model != "gaussian_arx", zero_padding == FALSE),
  nobs = 120,
  zero_padding = FALSE,
  rho = 0.9,
  metric = "energy_score",
  metric_name = "ENERGY_SCORE",
  pars = "beta",
  direct_labels = TRUE,
  xbreaks = c(20, 100, 200, 400)
) + papertheme +
theme(axis.title.x = element_blank()) + 
  coord_capped_cart(bottom = "right", gap = 0)


energy_score_120_beta_unpadded <- energy_score_120_beta_rho0_unpadded + energy_score_120_beta_rho05_unpadded + energy_score_120_beta_rho09_unpadded + plot_layout(axis_titles = "collect", design = ("AAAABBBBCCCCC")) & labs(x = "No. of exogen.\npredictors")

gc()

ggsave(paste0(summary_path, "arx_estim_energy_score_beta_nonzero.pdf"), energy_score_120_beta_unpadded, width = 5.5, height = 2.3)

bayesflow::save_tikz_plot(energy_score_120_beta_unpadded, paste0(summary_path, "arx_estim_energy_score_beta_nonzero.tex"), width = 5.5, height = 2.3)



# MAE

mae_120_beta_rho0_padded <- metric_plot(
  dat_tidy |> filter(dataset == "cement", model != "gaussian_arx"),
  nobs = 120,
  zero_padding = TRUE,
  rho = 0,
  "mae",
  "MAE",
  "beta",
  xbreaks = c(20, 100, 200, 400)
) +
  papertheme + ggtitle(NULL, subtitle = "$\\beta$") +
theme(axis.title.x = element_blank()) +
coord_capped_cart(bottom = "right", gap = 0)

mae_120_beta_rho05_padded <- metric_plot(
  dat_tidy |> filter(dataset == "cement", model != "gaussian_arx"),
  nobs = 120,
  zero_padding = TRUE,
  rho = 0.5,
  "mae",
  "MAE",
  "beta",
  xbreaks = c(20, 100, 200, 400)
) +
  papertheme +
coord_capped_cart(bottom = "right", gap = 0)

gc()

mae_120_beta_rho09_padded <- metric_plot(
  dat_tidy |> filter(dataset == "cement", model != "gaussian_arx"),
  nobs = 120,
  zero_padding = TRUE,
  rho = 0.9,
  metric = "mae",
  metric_name = "MAE",
  pars = "beta",
  direct_labels = TRUE,
  xbreaks = c(20, 100, 200, 400)
) + papertheme +
theme(axis.title.x = element_blank()) + 
  coord_capped_cart(bottom = "right", gap = 0)


mae_120_beta_padded <- mae_120_beta_rho0_padded + mae_120_beta_rho05_padded + mae_120_beta_rho09_padded + plot_layout(axis_titles = "collect", design = ("AAAABBBBCCCCC")) & labs(x = "No. of exogen.\npredictors")

gc()

mae_120_phi_rho0_padded <- metric_plot(
  dat_tidy |> filter(dataset == "cement", model != "gaussian_arx"),
  nobs = 120,
  zero_padding = TRUE,
  rho = 0,
  "mae",
  "MAE",
  "phi",
  xbreaks = c(20, 100, 200, 400)
) +
  papertheme + ggtitle(NULL, subtitle = "$\\phi$") +
  theme(axis.title.x = element_blank()) + 
coord_capped_cart(bottom = "right", gap = 0)

gc()

mae_120_phi_rho05_padded <- metric_plot(
  dat_tidy |> filter(dataset == "cement", model != "gaussian_arx"),
  nobs = 120,
  zero_padding = TRUE,
  rho = 0.5,
  "mae",
  "MAE",
  "phi",
  xbreaks = c(20, 100, 200, 400)
) +
  papertheme +
coord_capped_cart(bottom = "right", gap = 0)

gc()

mae_120_phi_rho09_padded <- metric_plot(
  dat_tidy |> filter(dataset == "cement", model != "gaussian_arx"),
  nobs = 120,
  zero_padding = TRUE,
  rho = 0.9,
  metric = "mae",
  metric_name = "MAE",
  pars = "phi",
  direct_labels = TRUE,
  xbreaks = c(20, 100, 200, 400)
) + papertheme +
theme(axis.title.x = element_blank()) + 
coord_capped_cart(bottom = "right", gap = 0)

mae_120_phi_padded <- mae_120_phi_rho0_padded + mae_120_phi_rho05_padded + mae_120_phi_rho09_padded + plot_layout(axis_titles = "collect", design = ("AAAABBBBCCCCC")) & labs(x = "No. of exogen.\npredictors")


ggsave(paste0(summary_path, "arx_estim_mae_beta.pdf"), mae_120_beta_padded, width = 5.5, height = 2.3)

bayesflow::save_tikz_plot(mae_120_beta_padded, paste0(summary_path, "arx_estim_mae_beta.tex"), width = 5.5, height = 2.3)



ggsave(paste0(summary_path, "arx_estim_mae_phi.pdf"), mae_120_phi_padded, width = 5.5, height = 2.3)

bayesflow::save_tikz_plot(mae_120_phi_padded, paste0(summary_path, "arx_estim_mae_phi.tex"), width = 5.5, height = 2.3)



mae_120_beta_rho0_unpadded <- metric_plot(
  dat_tidy |> filter(dataset == "cement", model != "gaussian_arx"),
  nobs = 120,
  zero_padding = FALSE,
  rho = 0,
  "mae",
  "MAE",
  "beta",
  xbreaks = c(20, 100, 200, 400)
) +
  papertheme + ggtitle(NULL, subtitle = "$\\beta$") +
theme(axis.title.x = element_blank()) +
coord_capped_cart(bottom = "right", gap = 0)

mae_120_beta_rho05_unpadded <- metric_plot(
  dat_tidy |> filter(dataset == "cement", model != "gaussian_arx"),
  nobs = 120,
  zero_padding = FALSE,
  rho = 0.5,
  "mae",
  "MAE",
  "beta",
  xbreaks = c(20, 100, 200, 400)
) +
  papertheme +
coord_capped_cart(bottom = "right", gap = 0)

gc()

mae_120_beta_rho09_unpadded <- metric_plot(
  dat_tidy |> filter(dataset == "cement", model != "gaussian_arx", zero_padding == FALSE),
  nobs = 120,
  zero_padding = FALSE,
  rho = 0.9,
  metric = "mae",
  metric_name = "MAE",
  pars = "beta",
  direct_labels = TRUE,
  xbreaks = c(20, 100, 200, 400)
) + papertheme +
theme(axis.title.x = element_blank()) + 
  coord_capped_cart(bottom = "right", gap = 0)


mae_120_beta_unpadded <- mae_120_beta_rho0_unpadded + mae_120_beta_rho05_unpadded + mae_120_beta_rho09_unpadded + plot_layout(axis_titles = "collect", design = ("AAAABBBBCCCCC")) & labs(x = "No. of exogen.\npredictors")

gc()

ggsave(paste0(summary_path, "arx_estim_mae_beta_nonzero.pdf"), mae_120_beta_unpadded, width = 5.5, height = 2.3)

bayesflow::save_tikz_plot(mae_120_beta_unpadded, paste0(summary_path, "arx_estim_mae_beta_nonzero.tex"), width = 5.5, height = 2.3)
