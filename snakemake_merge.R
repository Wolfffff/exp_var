library(here)

file_paths <- list.files(path = here("snakemake/Rdatas/residuals/"), pattern = "\\.rds", full.names = TRUE)

file_names <-  gsub(pattern = "\\.rds$", replacement = "", x = basename(file_paths))

data_list <- lapply(file_paths, readRDS)

names(data_list) <- file_names  