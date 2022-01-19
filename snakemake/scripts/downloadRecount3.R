my_logfile = snakemake@log[[1]]
snakemake@source("logger.R")
log4r_info("Starting.")

log4r_info("Loading packages") 
source("../functions.R")

# Load metadata from recount3
log4r_info("Reading metadata")
experimental_metadata_rc3 <- read_csv(snakemake@input[["metadata"]])

log4r_info("Downloading cache")
human_projects <- available_projects(bfc = cache)
id = snakemake@wildcards[["id"]]

log4r_info("Downloading data")
exp_data_rc3 = downloadRecount3(id)

log4r_info("Saving data")
saveRDS(exp_data_rc3, file = snakemake@output[[1]])

log4r_info("Done!")
