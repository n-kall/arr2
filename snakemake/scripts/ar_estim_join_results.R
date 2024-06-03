library(tidyverse)
library(vroom)

if (exists("snakemake")) {
  results_path <- snakemake@config[["results_path"]]
  summary_path <- snakemake@config[["summary_path"]]
  results_files <- snakemake@input[["results"]]
  source(snakemake@input[["functions"]])
} else {
  results_path <- "snakemake/results/ar_estim/"
  results_files <- list.files(results_path, pattern = "*.csv", full.names = TRUE)
  summary_path <- "snakemake/summary/"
  source("snakemake/scripts/functions.R")
}

dat <- vroom_chunked(results_files, id = "filename")

dat_tidy <- dat |>
  mutate(
    filename = str_remove(
      filename,
      paste0(results_path, "/ar/")
    )
  ) |>
  separate(filename, into = c("model", "npars", "dataset", "rep"), sep = "__") |>
  mutate(
    rep = str_remove(rep, ".csv")
  ) |>
  separate(dataset, into = c("dataset", "remove", "nobs"), sep = "-") |>
  mutate(npars = as.numeric(npars) * as.numeric(nobs)/24)

write_csv(dat_tidy, paste0(summary_path, "ar_estim_results.csv"))
