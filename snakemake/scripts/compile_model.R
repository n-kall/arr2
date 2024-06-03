library(cmdstanr)
library(stringr)

stanfile <- snakemake@input[[1]]
binaryfile <- snakemake@output[[1]]

m <- cmdstan_model(stan_file = stanfile, exe_file = binaryfile)
