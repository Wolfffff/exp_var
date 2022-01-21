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

results_list_rc3

# Process row wise metrics
library(dplyr)
library(Rfast)
library(matrixStats)


row_means_df = calculate_row_wise_metric(data_list,rowMeans)
row_var_df = calculate_row_wise_metric(data_list,rowVars)
row_sd_df = calculate_row_wise_metric(data_list,rowSds)
row_cv_df = calculate_row_wise_metric(data_list,rowcvs)

