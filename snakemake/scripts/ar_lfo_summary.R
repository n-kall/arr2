library(tidyverse)
library(vroom)

if (exists("snakemake")) {
  output_file <- snakemake@output[["summary_file"]]
  results_path <- snakemake@config[["results_path"]]
  summary_path <- snakemake@config[["summary_path"]]
  results_files <- snakemake@input[["results_files"]]
  functions <- snakemake@input[["functions"]]
} else {
  results_path <- "snakemake/results/ar_lfo"
  results_files <- list.files(results_path, pattern = "*.csv", full.names = TRUE)
  summary_path <- "snakemake/summary/"
  output_file <- "snakemake/summary/ar_lfo.csv"
}

source(functions)

dat <- vroom_chunked(results_files, id = "filename", delim = ",")

dat_tidy <- dat |>
  mutate(
    filename = str_remove(
      filename,
      paste0(results_path, "/"))
  ) |>
  separate(filename, into = c("model", "npars", "dataset", "rep"), sep = "__") |>
  mutate(rep = str_remove(rep, ".csv"))

write_csv(dat_tidy, output_file)
