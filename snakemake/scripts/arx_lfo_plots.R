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
  results_path <- "snakemake/results/arx_lfo/"
  results_files <- list.files(results_path, pattern = "*.csv", full.names = TRUE)
  summary_path <- "snakemake/summary/"
}

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
  arxr2_gamma_minnesota = "ARR2 (gamma, minnesota)",
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


metric_plot <- function(dat, metric, metric_name, direct_labels = FALSE) {

  p <- dat |>
    group_by(npars_x, model, dataset, rho) |>
    summarise(mean_metric = mean(.data[[metric]]), stdev = sd(.data[[metric]]), nreps = n())  |>
    mutate(model_label = if_else(rho == 0.9, model_labels[model], NA)) |>
    ggplot(aes(x = as.integer(npars_x), color = model, group = model)) +
    geom_pointrange(aes(
      y = mean_metric,
      ymin = mean_metric - stdev /sqrt(nreps),
      ymax = mean_metric + stdev /sqrt(nreps)
    ),
    position = position_dodge(width = 15),
    size = 0.001
    ) +
    geom_line(
    aes(
      y = mean_metric),
    position = position_dodge(width = 15)
    ) +
#  facet_rep_wrap(~rho, ncol = 3, labeller = as_labeller(dataset_labeller), repeat.tick.labels = TRUE) +
    xlab("Number of exogenous predictors") +
    ylab(metric_name) +
  scale_color_manual(
    name = "Prior",
    labels = model_labels,
    guide = "none",
    values = colours
  ) +
    scale_x_continuous(breaks = c(20, 100, 200, 400), expand = expansion(mult = c(0.06, 0.2)))
  if (direct_labels) {

    p <- p + geom_dl(aes(y = mean_metric, label = model_label), method = list(dl.trans(x = x + 0.2), "last.points", "bumpup", cex = 0.8)) +
      scale_x_continuous(expand = expansion(mult = c(0.06, 0.4)), breaks = c(20, 100, 200, 400))
  }
  p
}

dat <- vroom_chunked(results_files, id = "filename")

dat_tidy <- dat |>
  mutate(
    filename = str_remove(
      filename,
      paste0(results_path, "/"))
  ) |>
  separate(filename, into = c("model", "dataset", "rep"), sep = "__") |>
  mutate(rep = str_remove(rep, ".csv")) |>
  separate(dataset, into = c("dataset", "n", "n_obs", "r", "rho", "mcol", "npars_x"), sep = "-") |>
  filter(model != "gaussian_arx", dataset == "cement") |>
  mutate(mlpd = elpd / (as.numeric(n_obs)/2))



mlpd_120_rho0 <- metric_plot(
  dat_tidy |> filter(rho == 0), "mlpd", "MLPD") + papertheme + coord_capped_cart(bottom = "right", gap = 0) +
  theme(axis.title.x = element_blank())

mlpd_120_rho05 <- metric_plot(
  dat_tidy |> filter(rho == 0.5), "mlpd", "MLPD") + papertheme + coord_capped_cart(bottom = "right", gap = 0) +
  theme(axis.title.y = element_blank()
        )

mlpd_120_rho09 <- metric_plot(
  dat_tidy |> filter(rho == 0.9), "mlpd", "MLPD", direct_labels = TRUE) + papertheme + coord_capped_cart(bottom = "right", gap = 0) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())


mlpd_120 <- mlpd_120_rho0 + mlpd_120_rho05 + mlpd_120_rho09 + plot_layout(axis_titles = "collect", design = ("AAAAABBBBBCCCCCC")) & labs(x = "No. of exogen.\npredictors")# & scale_y_continuous(limits = c(-1.025*2, -0.7*2))




rhat_120_rho0 <- metric_plot(
  dat_tidy |> filter(rho == 0), "rhat", "RHAT") + papertheme + coord_capped_cart(bottom = "right", gap = 0) +
  theme(axis.title.x = element_blank())

rhat_120_rho05 <- metric_plot(
  dat_tidy |> filter(rho == 0.5), "rhat", "RHAT") + papertheme + coord_capped_cart(bottom = "right", gap = 0) +
  theme(axis.title.y = element_blank()
        )

rhat_120_rho09 <- metric_plot(
  dat_tidy |> filter(rho == 0.9), "rhat", "RHAT", direct_labels = TRUE) + papertheme + coord_capped_cart(bottom = "right", gap = 0) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())


ggsave(paste0(summary_path, "arx_mlpd_120.pdf"), mlpd_120, width = 5.5, height = 2.2)

bayesflow::save_tikz_plot(mlpd_120, paste0(summary_path, "arx_mlpd_120.tex"), width = 5.5, height = 2.2)
