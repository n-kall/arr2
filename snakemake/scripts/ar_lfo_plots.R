library(tidyverse)
library(patchwork)
library(vroom)
library(directlabels)
library(lemon)

if (exists("snakemake")) {
  results_path <- snakemake@config[["results_path"]]
  summary_path <- snakemake@config[["summary_path"]]
  results_files <- snakemake@input[["results"]]
} else {
  summary_file <- "snakemake/summary/ar_lfo.csv"
  summary_path <- "snakemake/summary/"
}

colours <- c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#AA3377")

dataset_labeller <- c(
  ar8_delayed = "Delayed",
  ar8_minnesota = "Minnesota",
  ar8_damposc = "Dampened oscillations",
  whitenoise = "White noise"
)

model_labels <- c(
  arr2_flat = "ARR2 (simplex, flat)",
  arr2_gamma_ncp_flat = "ARR2 (flat)",
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

metric_plot <- function(dat, nobs, metric, metric_name, dataset_name, dodge_width = 0, direct_labels = TRUE, errorbar_width = 0, distance = 0, xbreaks = waiver(), ybreaks = waiver()) {

  p <- dat |>
    mutate(npars = as.numeric(npars) * as.numeric(nobs)/24) |>
    filter(
      nobs == !!nobs,
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

dat_tidy <- read_csv(summary_file)

dat_tidy <- dat_tidy |>
  separate(dataset, into = c("dataset", "ncol", "nobs"), "-") |>
  mutate(mlpd = elpd / (as.numeric(nobs)/2))


minnesota_mlpd <- metric_plot(
  dat_tidy,
  nobs = 120,
  metric = "mlpd",
  metric_name = "MLPD",
  xbreaks = seq(10, 60, by = 10),
  dodge_width = 1,
  dataset_name = "ar8_minnesota",
  direct_labels = FALSE
) +
  papertheme + coord_capped_cart(bottom = "right", gap = 0) +
theme(axis.title.x = element_blank())


delayed_mlpd <- metric_plot(
  dat_tidy,
  nobs = 120,
  metric = "mlpd",
  metric_name = "MLPD",
  xbreaks = seq(10, 60, by = 10),
  dodge_width = 1,
  dataset_name = "ar8_delayed",
  direct_labels = FALSE
) +
  papertheme + coord_capped_cart(bottom = "right", gap = 0) +
theme(
   axis.title.y = element_blank()
)

damposc_mlpd <- metric_plot(
  dat_tidy,
  nobs = 120,
  metric = "mlpd",
  metric_name = "MLPD",
  xbreaks = seq(10, 60, by = 10),
  dodge_width = 1,
  dataset_name = "ar8_damposc"
) +
  papertheme + coord_capped_cart(bottom = "right", gap = 0) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()
        )

(mlpd_120 <- minnesota_mlpd + delayed_mlpd + damposc_mlpd + plot_layout(design = "AAAABBBBCCCCC") & coord_cartesian(ylim = c(-2, -1.4)))


bayesflow::save_tikz_plot(mlpd_120, paste0(summary_path, "ar_mlpd_120.tex"), width = 5.5, height = 2)

ggsave(paste0(summary_path, "ar_mlpd_120.pdf"), mlpd_120, width = 5.5, height = 2)
