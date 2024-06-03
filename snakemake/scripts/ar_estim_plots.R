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
  results_file <- "snakemake/summary/ar_estim_results.csv"

 }

colours <- c("black", "#4477AA", "#EE6677", "#228833", "grey", "yellow", "#CCBB44", "#AA3377")

dataset_labeller <- c(
  ar8_delayed = "Delayed",
  ar8_minnesota = "Minnesota",
  ar8_damposc = "Dampened oscillations",
  ar8_delayed_small = "Delayed AR(8), small R2",
  ar8_minnesota_small = "Minnesota AR(8), small R2",
  ar8_damposc_small = "Dampened oscillations AR(8), small R2",
  whitenoise = "White noise"
)

model_labels <- c(
  arr2_flat = "ARR2 (simplex, flat)",
  arr2_gamma_ncp_flat = "ARR2 (flat)",
  arr2_gamma_ncp_flat0.5 = "ARR2 (flat, 0.5)",
  arr2_gamma_ncp_flat_heavytail = "ARR2 (flat-heavytail)",
  arr2_gamma_ncp_minnesota_heavytail = "ARR2 (minn-heavytail)",
  arr2_minnesota = "ARR2 (simplex, minnesota)",
  arr2_gamma_ncp_minnesota = "ARR2 (Minn.)",
  gaussian = "Gauss.",
  minnesota =  "Minn.",
  rhs = "RHS",
  l1ball = "L1 ball",
  arr2_dc_flat = "ARR2 (dc, flat)",
  arr2_dc_minnesota = "ARR2 (dc, minnesota)"
)


papertheme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
    axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black")
  )

r2_plot <- function(dat, nobs, breaks, dodge_width = 0, ...) {
  dat |>
    filter(
      variable == "postR2",
      nobs == !!nobs,
      zero_padding == FALSE,
      pars == "postR2"
    ) |>
    group_by(npars, model, dataset) |>
    summarise(
      mean_post_R2 = mean(mean),
      q5_post_R2 = mean(q5),
      q95_post_R2 = mean(q95),
      mean_true_R2 = mean(true_par),
      q5_true_R2 = quantile(true_par, probs = 0.05),
      q95_true_R2 = quantile(true_par, probs = 0.95)
    ) |>
    mutate(model_label = model_labels[model]) |>
    ggplot(aes(x = npars, colour = model, group = model)) +
    geom_dl(aes(y = mean_post_R2, label = model_label), method = list(dl.trans(x = x + 0.2), "last.points", "bumpup", cex = 0.8)) +
    geom_pointrange(
      aes(y = mean_post_R2,
          ymin = q5_post_R2,
          ymax = q95_post_R2
          ), position = position_dodge2(width = dodge_width)) +
    geom_line(aes(y = mean_true_R2), colour = "black") +
    geom_line(aes(y = q5_true_R2), colour = "black", linetype = "dashed") +
    geom_line(aes(y = q95_true_R2), colour = "black", linetype = "dashed") +
    scale_color_manual(
      name = "Prior",
      labels = model_labels,
      guide = "none",
      values = colours
    ) +
    facet_rep_wrap(~dataset, ncol = 2, labeller = as_labeller(dataset_labeller), repeat.tick.labels = TRUE) +
    scale_x_continuous(expand = expansion(mult = c(0.06, 0.4)), breaks = breaks) +
    xlab("AR order")+
    ylab("Posterior $R^2$")

}


metric_plot <- function(dat, nobs, zero_padding, dataset_name, metric, metric_name, xbreaks = waiver(), ybreaks = waiver(), ylimits = waiver(), dodge_width = 0, errorbar_width = 0, direct_labels = TRUE, pars = "phi") {

  p <- dat |>
    filter(
      variable == "phi[1]",
      nobs == !!nobs,
      zero_padding == !!zero_padding,
      pars == !!pars,
      dataset == dataset_name
    ) |>
    group_by(npars, model, dataset) |>
    summarise(mean_metric = mean(.data[[metric]]), stdev = sd(.data[[metric]]), nreps = n()) |>
    mutate(model_label = model_labels[model]) |>
    ggplot(aes(x = as.integer(npars), color = model, group = model)) +
    geom_pointrange(aes(
      y = mean_metric,
      ymin = mean_metric - stdev / sqrt(nreps),
      ymax = mean_metric + stdev / sqrt(nreps),
      ),
      position = position_dodge(width = dodge_width),
      size = 0.001
      ) +
    geom_errorbar(aes(
      ymin = mean_metric - stdev / sqrt(nreps),
      ymax = mean_metric + stdev / sqrt(nreps),
      width = errorbar_width
    ),
    position = position_dodge(width = dodge_width),
    ) +
    geom_line(
      aes(
        y = mean_metric),
      position = position_dodge(width = dodge_width)
    ) +
    facet_rep_wrap(~dataset, ncol = 3, scales = "free_y", labeller = as_labeller(dataset_labeller), repeat.tick.labels = TRUE) +
    xlab("Model AR order") +
    ylab(metric_name) +
    scale_color_manual(
      name = "Prior",
      labels = model_labels,
      guide = "none",
      values = colours
    ) +
#    scale_y_log10(breaks = ybreaks, limits = ylimits) +
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
  mutate(model = str_remove(model, "results/ar_estim/")) |>
  filter(dataset != "whitenoise",
         !(model %in% c("arr2_gamma_ncp_flat_highR2", "arr2_gamma_ncp_flat_flatR2")))


minnesota_ar8_rmse_padded <- metric_plot(dat_tidy, dataset_name = "ar8_minnesota", direct_labels = FALSE, nobs = 120, zero_padding = TRUE, metric = "rmse", metric_name = "RMSE", xbreaks = seq(10, 60, 10), ybreaks =  seq(0.02, 0.3, length.out = 5), dodge_width = 1) + papertheme + coord_capped_cart(bottom = "right", gap = 0) +
  theme(
    axis.title.x = element_blank(),
    )

delayed_ar8_rmse_padded <- metric_plot(dat_tidy, dataset_name = "ar8_delayed", direct_labels = FALSE, nobs = 120, zero_padding = TRUE, metric = "rmse", metric_name = "RMSE", xbreaks = seq(10, 60, 10),  ybreaks = seq(0.02, 0.3, length.out = 5), dodge_width = 1) + papertheme + coord_capped_cart(bottom = "right", gap = 0) + theme(
  axis.title.y = element_blank()
)

damposc_ar8_rmse_padded <- metric_plot(dat_tidy, dataset_name = "ar8_damposc", direct_labels = TRUE, nobs = 120, zero_padding = TRUE, metric = "rmse", metric_name = "RMSE", xbreaks = seq(10, 60, 10),  ybreaks = seq(0.02, 0.3, length.out = 5), dodge_width = 1) + papertheme + coord_capped_cart(bottom = "right", gap = 0) +
  theme(
    axis.title.x = element_blank(),
        axis.title.y = element_blank()
        )

(rmse_120_padded <- minnesota_ar8_rmse_padded + delayed_ar8_rmse_padded + damposc_ar8_rmse_padded + plot_layout(axis_titles = "collect", design = ("AAAABBBBCCCCC")) & scale_y_continuous(limits = c(0, 0.32), breaks = c(0, 0.1, 0.2, 0.3)))


ggsave(paste0(summary_path, "ar_estim_padded_120_extras.pdf"), rmse_120_padded, width = 5.5, height = 2.1)
bayesflow::save_tikz_plot(rmse_120_padded, paste0(summary_path, "ar_estim_padded_120_extras.tex"), width = 5.5, height = 2.1)
