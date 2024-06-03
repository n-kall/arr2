library(bayesflow) 
library(stringr)

functions <- snakemake@input[["functions"]]
SEED <- as.numeric(snakemake@config[["seed"]])

wildcards <- snakemake@wildcards[names(snakemake@wildcards) != ""]
dgp_parameters <- lapply(wildcards[names(wildcards) != "dgp"], as.numeric)
dgp_script <- snakemake@input[["script"]]
output_file <- snakemake@output[["dataset_file"]]

source(functions)
source(dgp_script)

datasets <- do.call(dgp, args = dgp_parameters)

saveRDS(datasets, output_file)
