library(tidyverse)
library(vroom)

if (exists("snakemake")) {
  output_file <- snakemake@output[[1]]
  results_path <- snakemake@config[["results_path"]]
  summary_path <- snakemake@config[["summary_path"]]
  results_files <- snakemake@input[["results"]]
  functions <- snakemake@input[["functions"]]
} else {
  results_path <- "snakemake/results/inflation"
  results_files <- list.files(results_path, pattern = "*.csv", full.names = TRUE)
  summary_path <- "snakemake/summary/"
}

source(functions)


dat <- vroom_chunked(results_files, id = "filename")

dat_tidy <- dat |>
  mutate(
    filename = str_remove(
      filename,
      paste0(results_path, "inflation/"))
  ) |>
  separate(filename, into = c("model", "dataset", "fold"), sep = "__") |>
  mutate(fold = str_remove(fold, ".csv"))

write_csv(dat_tidy, output_file)
