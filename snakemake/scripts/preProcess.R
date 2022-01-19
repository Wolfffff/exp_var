save.image(snakemake@log[["env"]])

my_logfile = snakemake@log[["log"]]
snakemake@source("logger.R")
log4r_info("Starting.")
print = log4r_info
log4r_info("Loading packages") 
source("../functions.R")

log4r_info("Reading metadata")
experimental_metadata <- read_csv(snakemake@input[["metadata"]])

feature_vec <- list()
feature_vec[["disease"]] <- c("normal", "control", "", NA,
                              "non inflammatory bowel disease control")
feature_vec[["treatment"]] <- c("normal", "control", "", NA)

dset_name = snakemake@wildcards[["id"]]

log4r_info("Reading data")
dset <- readRDS(snakemake@input[["data"]])

assay_name = snakemake@config[["assay_name"]]

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

switch(dset_name,
SKIN = {
    remove_cell_culture = metadata[["gtex.smtsd"]] != "Cells - Cultured fibroblasts"
    metadata = metadata[ remove_cell_culture,]
    counts   =   counts[, remove_cell_culture]
},
ESOPHAGUS = {
    remove_cell_culture = metadata[["gtex.smtsd"]] == "Esophagus - Mucosa"
    metadata = metadata[ remove_cell_culture,]
    counts   =   counts[, remove_cell_culture]
}
)

log4r_info(paste0("Unfiltered count dimensions: ", dim(counts)[1], " x ", dim(counts)[2]))
log4r_info(paste0("Unfiltered metadata dimensions: ", dim(metadata)[1], " x ", dim(metadata)[2]))

filtered_data <- make_filtered_data(counts, metadata, feature_vec)
metadata <- filtered_data$metadata
counts <- filtered_data$counts
n_samples <- dim(metadata)[1]
metadata$sample_id <- rownames(metadata)

log4r_info("Normalizing and estimating mean-variance weights...")
countdata.list <- DGEList(counts = counts, samples = metadata, genes = rownames(dset))
rep_names = c("technical_replicate_group", "wells.replicate", "individual")
if (any(names(metadata) %in% rep_names)) {
log4r_info("Summing technical replicates...")
rep_col = names(metadata)[names(metadata) %in% rep_names]
countdata.list <- sumTechReps(countdata.list, metadata[[rep_col]])
}

log4r_info("Calculating normalization factors...")
countdata.norm <- calcNormFactors(countdata.list)

log4r_info("Trimming...")
cutoff <- inv_log2_plus05(1)
drop <- which(apply(cpm(countdata.norm), 1, max) < cutoff)
countdata.norm <- countdata.norm[-drop, ]

# recalc with updated metadata
values_count <- sapply(lapply(countdata.norm$samples, unique), length)
countdata.norm$samples <- countdata.norm$samples[, names(countdata.norm$samples[, values_count > 1])]

log4r_info("Making design matrix...")
countdata.norm$samples  <- remove_redundant_features(countdata.norm$samples )

# Removing columns with a crazy number of levels that mess everything up.
# (this is why we have random effects by the way)
countdata.norm$samples  <- remove_large_factors(
    countdata.norm$samples,
    columns_to_ignore
)
countdata.norm$samples <- select_meta(countdata.norm$samples)

# Switch to DESeq2

# Remove top 3 genes from BLOOD
if(dset_name == "BLOOD"){
    bigones = sort(apply(countdata.norm$counts, 1, max), decreasing = T)
    remove_genes = which(rownames(countdata.norm) %in% names(bigones)[1:3])
    countdata.norm  =   countdata.norm[-remove_genes,]
}
# Remove top 2 genes from LIVER
if(dset_name == "LIVER"){
    bigones = sort(apply(countdata.norm$counts, 1, max), decreasing = T)
    remove_genes = which(rownames(countdata.norm) %in% names(bigones)[1:2])
    countdata.norm  =   countdata.norm[-remove_genes,]
}
# Remove top 1 genes from COLON and STOMACH
if(dset_name %in% c("COLON", "STOMACH")){
    bigones = sort(apply(countdata.norm$counts, 1, max), decreasing = T)
    remove_genes = which(rownames(countdata.norm) %in% names(bigones)[1])
    countdata.norm  =   countdata.norm[-remove_genes,]
}

log4r_info("Saving data...")
saveRDS(list(data = countdata.norm, columns_to_ignore = columns_to_ignore), 
        file = snakemake@output[[1]])
log4r_info("Done!")