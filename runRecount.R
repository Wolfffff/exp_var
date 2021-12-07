library(ExpressionAtlas)
library(plyr)
library(tidyverse)
library(limma)
library(sva)
library(edgeR)
library(ggplot2)
library(janitor)
library(foreach)
library(doParallel)
library(biomaRt)
library(ggfortify)
library(patchwork)
library(cowplot)
library(clusterProfiler)
library(AnnotationDbi)
#pak::pkg_install(c("DESeq2", "vsn"))
library(org.Hs.eg.db)
library(rrcov)
library(DESeq2)
library(vsn)
library(viridis)
source("functions.R")

# Set timeout to avoid failure when trying to download GTEx or other large datasets
options(timeout = 1800)

experimental_metadata_ea <- read_csv("expression_atlas_metadata.csv")
experimental_metadata_ea <- experimental_metadata_ea[-1,] # Remove the top row (GTEx) to lower the time cost when testing

# Load metadata from recount3
experimental_metadata_rc3 <- read_csv("recount3_metadata.csv")


plots_dir <- "data/"

# Setup filters for removing samples accordingly to the metadata
feature_vec <- list()
feature_vec[["disease"]] <- c("normal", "control", "", NA, 
                              "non inflammatory bowel disease control")
feature_vec[["treatment"]] <- c("normal", "control", "", NA)

# This is commented out as we're sticking with Expression Atlas for now.
library(recount3)

# Move cache to deal with quota issues
cache = recount3_cache(cache_dir="cache")
human_projects <- available_projects(bfc=cache)

# exp_data_rc3 = list()
# # recount_cache_rm()
# for (id in experimental_metadata_rc3$id) {
#   # Load the project
#   proj_info <- subset(
#     human_projects,
#     project == id & project_type == "data_sources"
#   )
#   # Tape to deal with acquisition issues
#   rse <- create_rse(proj_info, bfc=cache)
#   if (proj_info$file_source == "gtex") {
#     metadata_df <- colData(rse)[,grepl("gtex",colnames(colData(rse)),fixed=TRUE)]
#   } else {
#     #Could probably just use expand_sra_attributes
#     metadata_df <- convert_metadata_to_df_rc3(colData(rse)$sra.sample_attributes)
#   }
  
#   # Crude way to set NA to ""
#   metadata_df@listData <- lapply(metadata_df@listData, function(x) {
#     x[is.na(x)] = ""
#     x
#   })

#   rse@colData <- metadata_df
#   exp_data_rc3[[id]] <- rse
# }
# saveRDS(exp_data_rc3, file = "cache/recount3_data.RDS")
exp_data_rc3 = readRDS("cache/recount3_data.RDS")

source("./main_processing_loop.R")
parallel = FALSE
if(.Platform$OS.type == "unix") {
  parallel = TRUE
  library(doMC)
  registerDoMC(32)
} 

results_list_rc3 <- llply(names(exp_data_rc3), 
                      main_loop,  
                      exp_data = exp_data_rc3, 
                      experimental_metadata = experimental_metadata_rc3, 
                      feature_vec = feature_vec,
                      assay_name = "raw_counts", 
                      .parallel = parallel)  

names(results_list_rc3) <- names(exp_data_rc3)
                    
for(dset_name in names(results_list_rc3)) 
  save_plot(filename=paste0(plots_dir, dset_name,"_pca.png"), results_list_rc3[[dset_name]]$plotPanel, 
            base_height = 6, base_asp = 1.2, ncol = 2, nrow = 2)

save(results_list_rc3, file = "results_list_rc3.RData")