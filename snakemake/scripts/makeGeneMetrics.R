source(here::here("functions.R"))
library(here)

file_paths <- list.files(path = here::here("snakemake/Rdatas/residuals/"), 
                         pattern = "\\.rds", full.names = TRUE)
file_names <-  gsub(pattern = "\\.rds$", replacement = "", x = basename(file_paths))

library(plyr)
library(doMC)
registerDoMC(64)
data_list <- llply(file_paths, readRDS, .parallel = TRUE)

names(data_list) <- file_names  

# Process row wise metrics
library(dplyr)
library(Rfast)
library(matrixStats)

gene_metric_dfs <- list(
    means = calculate_row_wise_metric(data_list,rowMeans),
    var = calculate_row_wise_metric(data_list,rowVars),
    sd = calculate_row_wise_metric(data_list,rowSds),
    cv = calculate_row_wise_metric(data_list,rowcvs)
)

#saveRDS(gene_metric_dfs,  file = snakemake@output[[1]])
file.remove(here::here("snakemake/Rdatas/gene_metrics.RDS"))
saveRDS(gene_metric_dfs,  file = here::here("snakemake/Rdatas/gene_metrics.RDS"))
