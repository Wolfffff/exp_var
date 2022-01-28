my_logfile = snakemake@log[["log"]]
snakemake@source("logger.R")
log4r_info("Starting.")
print = log4r_info
log4r_info("Loading packages") 

source(here::here("functions.R"))
library(here)

file_paths <- list.files(path = here::here("snakemake/Rdatas/residuals/"), 
                         pattern = "\\.rds", full.names = TRUE)
file_names <-  gsub(pattern = "\\.rds$", replacement = "", x = basename(file_paths))

library(plyr)
library(doMC)
registerDoMC(64)

print("Reading residual files")
data_list <- llply(file_paths, readRDS, .parallel = TRUE)

names(data_list) <- file_names  
print(paste("Included", length(file_names), "studies:", paste(file_names, collapse = ", ")))

# Process row wise metrics
library(dplyr)
library(Rfast)
library(matrixStats)

max_missingness = 0.50
print(paste("Calculating gene level metrics with a max missingness of:", max_missingness))
gene_metric_dfs <- list(
    means = calculate_row_wise_metric_sparse(data_list, rowMeans, max_missingness),
    var = calculate_row_wise_metric_sparse(data_list, rowVars, max_missingness),
    sd = calculate_row_wise_metric_sparse(data_list, rowSds, max_missingness),
    cv = calculate_row_wise_metric_sparse(data_list, rowcvs, max_missingness)
)
print(paste("Kept", nrow(gene_metric_dfs[[1]]),"genes"))

#saveRDS(gene_metric_dfs,  file = snakemake@output[[1]])

file.remove(here::here("snakemake/Rdatas/gene_metrics.RDS"))
print("Saving data")
saveRDS(gene_metric_dfs,  file = here::here("snakemake/Rdatas/gene_metrics.RDS"))
print("Done!")
