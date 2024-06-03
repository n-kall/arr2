library(tidyverse)
library(vroom)

if (exists("snakemake")) {
  results_path <- snakemake@config[["results_path"]]
  summary_path <- snakemake@config[["summary_path"]]
  results_files <- snakemake@input[["results"]]
  functions <- snakemake@input[["functions"]]
} else {
  results_path <- "snakemake/results/arx_estim"
  results_files <- list.files(results_path, pattern = "*.csv", full.names = TRUE)
  summary_path <- "snakemake/summary/"
}

source(functions)

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
  mutate(
    n_obs = as.numeric(n_obs),
    rho = as.numeric(rho),
    npars_x = as.numeric(npars_x)
  )

write_csv(dat_tidy, paste0(summary_path, "arx_estim_results.csv"))
