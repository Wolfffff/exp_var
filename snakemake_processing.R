# rc3 input 
source(here::here("functions.R"))
library(here)

file_paths <- list.files(path = here("snakemake/Rdatas/residuals/"), pattern = "\\.rds", full.names = TRUE)
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

row_mean_dfs <- list(
    means = calculate_row_wise_metric(data_list,rowMeans),
    var = calculate_row_wise_metric(data_list,rowVars),
    sd = calculate_row_wise_metric(data_list,rowSds),
    cv = calculate_row_wise_metric(data_list,rowcvs)
)
