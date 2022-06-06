library(lubridate)
save.image(snakemake@log[["env"]])
# setwd(here::here("snakemake"))
# load("Rdatas/env/preProcess/BONE_MARROW.Rdata")

my_logfile = snakemake@log[["log"]]
snakemake@source("logger.R")
log4r_info("Starting.")
print = log4r_info
log4r_info("Loading packages")
source("../functions.R")
source("scripts/getMTgenes.R")

log4r_info("Reading metadata")
experimental_metadata <- read_csv(snakemake@input[["metadata"]])

dset_name = snakemake@wildcards[["id"]]

log4r_info("Reading data")
dset <- readRDS(snakemake@input[["data"]])

assay_name = snakemake@params[["assay_name"]]
# Catch empty columns_to_ignore
columns_to_ignore = tryCatch(
expr = {
    unlist(strsplit(experimental_metadata[experimental_metadata$id == dset_name, ]$columns_to_ignore, split = ";"))
},
error = function(e) {
    columns_to_ignore <- c("")
    columns_to_ignore
}
)
columns_to_ignore = c(columns_to_ignore, crap_cols)

metadata <- colData(dset)
metadata <- data.frame(metadata)
counts <- assays(dset)[[assay_name]]
log4r_info("Switch case")
switch(dset_name,
SKIN = {
    selection = metadata[["gtex.smtsd"]] != "Cells - Cultured fibroblasts"
    metadata = metadata[selection,]
    counts   =   counts[, selection]
},
ESOPHAGUS = {
    selection = metadata[["gtex.smtsd"]] == "Esophagus - Mucosa"
    metadata = metadata[ selection,]
    counts   =   counts[, selection]
}
)

filtered_data <- make_filtered_data(counts, metadata, feature_vec)
metadata <- filtered_data$metadata
counts <- filtered_data$counts
n_samples <- dim(metadata)[1]
metadata$sample_id <- rownames(metadata)

countdata.list <- DGEList(counts = counts, samples = metadata, genes = rownames(dset))
rep_names = c("technical_replicate_group", "wells.replicate", "individual", "sample")
if (any(names(metadata) %in% rep_names)) {
    log4r_info("Summing technical replicates...")
    rep_col = names(metadata)[names(metadata) %in% rep_names]
    if(length(rep_col) > 1){
        rep_col = rep_col[which.min(sapply(rep_col, function(x) length(levels(factor(metadata[,x])))))]
        log4r_info(paste0("More than one rep column. Choosing ", rep_col))
    }

    countdata.list <- sumTechReps(countdata.list, metadata[,rep_col])
}

# TODO: Remove calcNormFactors and change up naming
log4r_info("Calculating normalization factors...")
countdata.norm <- calcNormFactors(countdata.list)

log4r_info("Trimming...")
# cutoff <- inv_log2_plus05(1)

# Removing genes where max expression is less than cpm
drop <- which(apply(cpm(countdata.norm), 1, max) < 1)
countdata.norm <- countdata.norm[-drop, ]

# Removing genes with low average expression
drop <- which(apply(cpm(countdata.norm), 1, mean) < 5)
countdata.norm <- countdata.norm[-drop, ]

# Removing mitochondrial genes
drop <- which(remove_id_ver(countdata.norm$genes[,1]) %in% mt_gene_ids)
countdata.norm <- countdata.norm[-drop, ]

# Remove top 3 genes from BLOOD, hemoglobins
if(dset_name == "BLOOD"){
    bigones = sort(apply(countdata.norm$counts, 1, max), decreasing = T)
    remove_genes = which(rownames(countdata.norm) %in% names(bigones)[1:3])
    print(paste("Removed genes:", paste(countdata.norm$genes[remove_genes,], collapse = " ")))
    countdata.norm  =   countdata.norm[-remove_genes,]
}
# Remove top 1 genes from STOMACH, pepsin
if(dset_name %in% c("STOMACH")){
    bigones = sort(apply(countdata.norm$counts, 1, max), decreasing = T)
    remove_genes = which(rownames(countdata.norm) %in% names(bigones)[1])
    print(paste("Removed genes:", paste(countdata.norm$genes[remove_genes,], collapse = " ")))
    countdata.norm  =   countdata.norm[-remove_genes,]
}

# TODO: recalculate normalization factors
# Recalc with updated metadata
values_count <- sapply(lapply(countdata.norm$samples, unique), length)
countdata.norm$samples <- countdata.norm$samples[, names(countdata.norm$samples[, values_count > 1])]

log4r_info("Cleaning up redundant features and large factors...")
results =  remove_redundant_features(countdata.norm$samples)
countdata.norm$samples <- results$metadata

countdata.norm$samples  <- remove_large_factors(
    countdata.norm$samples,
    columns_to_ignore
)$metadata
results <- select_meta(countdata.norm$samples)
countdata.norm$samples = results$metadata
## TODO: whatever this is
# metadata= data.frame(rep_col=rep_col,columns_to_ignore = paste(columns_to_ignore, collapse = ";"))
# write.table(metadata, file = "metadata.csv", append=TRUE,col.names=!file.exists("metadata.csv"), sep = ",")

log4r_info("Saving data...")
saveRDS(list(data = countdata.norm, columns_to_ignore = columns_to_ignore),
        file = snakemake@output[[1]])
log4r_info("Done!")
