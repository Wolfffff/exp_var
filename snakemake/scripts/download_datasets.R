
options(recount3_url = "http://duffel.rail.bio/recount3")
if (snakemake@params[["source"]]=="recount3") {

    my_logfile = snakemake@log[[1]]
    snakemake@source("logger.R")
    log4r_info("Starting.")

    log4r_info("Loading packages")
    source(here::here("functions.R"))

    id = snakemake@wildcards[["id"]]

    log4r_info("Downloading data")
    exp_data_rc3 = downloadRecount3(id)

    log4r_info("Reading metadata")
    experimental_metadata_rc3 <- read_csv(snakemake@input[["metadata"]])

    log4r_info("Saving data")
    saveRDS(exp_data_rc3, file = snakemake@output[[1]])

    log4r_info("Done!")
} else if (snakemake@params[["source"]]=="EA") {
    my_logfile = snakemake@log[[1]]
    snakemake@source("logger.R")
    log4r_info("Starting.")

    # TODO: move this
    library(ExpressionAtlas)

    log4r_info("Loading packages")
    source(here::here("functions.R"))


    id = snakemake@wildcards[["id"]]
    rnaseq_exps <- getAtlasData(id)

    all_exps <- rnaseq_exps
    exps <- names(rnaseq_exps@listData)

    exp_data_ea <- lapply(exps, FUN = function(x) {
        return(all_exps[[x]]$rnaseq)
    })
    names(exp_data_ea) <- exps
    saveRDS(exp_data_ea[[1]], file = snakemake@output[[1]])
    log4r_info("Done!")
}

