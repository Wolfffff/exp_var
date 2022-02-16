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
    mean = calculate_row_wise_metric_sparse(data_list, rowMeans, max_missingness),
    var = calculate_row_wise_metric_sparse(data_list, rowVars, max_missingness),
    sd = calculate_row_wise_metric_sparse(data_list, rowSds, max_missingness),
    cv = calculate_row_wise_metric_sparse(data_list, rowcvs, max_missingness)
)
print(paste("Kept", nrow(gene_metric_dfs[[1]]),"genes"))

#saveRDS(gene_metric_dfs,  file = snakemake@output[[1]])

file.remove(here::here("snakemake/Rdatas/gene_metrics.RDS"))
print("Saving data")
saveRDS(gene_metric_dfs,  file = here::here("snakemake/Rdatas/gene_metrics.RDS"))


file_paths <- list.files(path = here::here("snakemake/Rdatas/networkStats/"), 
                         pattern = "\\.csv", full.names = TRUE)
file_names <-  gsub(pattern = "\\.csv$", replacement = "", x = basename(file_paths))4)

print("Reading network stats files")
networkStats_list <- llply(file_paths, read.csv, .parallel = TRUE)

x = networkStats_list[[1]]
connectivity_list = llply(networkStats_list, 
                          function(x) dplyr::select(x, Gene, WeightedDegree_fdr_1e.2))

print(paste("Calculating gene level mean connectivity unsing FDR of 1e-2 and Spearman correlations."))
connectivity = do.call(rbind, connectivity_list) %>%
    dplyr::filter(Gene %in% gene_metric_dfs$mean$Genes) %>%
    group_by(Gene) %>%
    summarise(mean = mean(WeightedDegree_fdr_1e.2, na.rm = T),
              median = median(WeightedDegree_fdr_1e.2, na.rm = T))
attributes(connectivity) = list("fdr" = 1e-2), "correlation" = "spearman"

print("Saving data")
saveRDS(connectivity,  file = here::here("snakemake/Rdatas/gene_connectivity.RDS"))
print("Done!")
