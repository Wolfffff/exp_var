save.image(snakemake@log[["env"]])

my_logfile = snakemake@log[["log"]]
snakemake@source("logger.R")
log4r_info("Starting.")
print = log4r_info
log4r_info("Loading packages") 

source(here::here("functions.R"))
library(here)

print("Reading residual RDS object...")
x <- readRDS(snakemake@input[[1]])
 
print("Extracting residuals...")
sbm_output = as_tibble(x$residuals_noOut) %>% 
    dplyr::mutate(Gene = rownames(x$residuals_noOut)) %>% 
    dplyr::select(Gene, everything())

print("Setting rownames...")
x = as.data.frame(sbm_output)
rownames(x) = remove_id_ver(x$Gene)
x$Gene = NULL

print("Writing data...")
write.table(x, file = snakemake@output[[1]], sep = "\t")
print("Done!")
