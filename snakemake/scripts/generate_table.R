
library(lubridate)

my_logfile = snakemake@log[["log"]]
snakemake@source("logger.R")
log4r_info("Starting.")
print = log4r_info
log4r_info("Loading packages")

source("../functions.R")
source("scripts/getMTgenes.R")

#provide column names
output_raw_metadata_df <- data.frame(matrix(ncol = 17, nrow = 0))
colnames(output_raw_metadata_df) <- c('id','datetime','group','tissue','raw_counts_genes','raw_counts_individuals','raw_metadata_individuals', 'raw_metadata_features', 'filtered_counts_genes', ' filtered_counts_individuals','filtered_metadata_individuals', 'filtered_metadata_features', 'to_ignore_manually_curated','to_ignore_redundant_features', 'to_ignore_large_factors', 'replicate_col')


for (idx in 1:length(snakemake@input[["raw"]])){
    raw_counts = readRDS(snakemake@input[["raw"]][idx])
    raw_metadata =  read_csv(snakemake@input[["raw_metadata"]][idx])
    dset_name = snakemake@params[["dset_names"]][idx]
    assay_name = snakemake@params[["assay_name"]][idx]

    columns_to_ignore = tryCatch(
    expr = {
            unlist(strsplit(raw_metadata[raw_metadata$id == dset_name, ]$columns_to_ignore, split = ";"))
        },
        error = function(e) {
            columns_to_ignore <- c("")
            columns_to_ignore
        }
    )
    columns_to_ignore = c(columns_to_ignore, crap_cols)
    columns_to_ignore_scsv = paste(columns_to_ignore, collapse = ";")

    metadata <- colData(raw_counts)
    metadata <- data.frame(metadata)
    counts <- assays(raw_counts)[[assay_name]]
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
    raw_meta = metadata
    filtered_data <- make_filtered_data(counts, metadata, feature_vec)

    metadata <- filtered_data$metadata
    counts <- filtered_data$counts
    n_samples <- dim(metadata)[1]
    metadata$sample_id <- rownames(metadata)

    countdata.list <- DGEList(counts = counts, samples = metadata, genes = rownames(raw_counts))
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
    else{
       rep_col = ""
        log4r_info("No technical replicate group found")
    }


    # Residual
# log4r_info("Calculating normalization factors...")
countdata.norm <- countdata.list#calcNormFactors()

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
results =  remove_redundant_features(countdata.norm$samples )
countdata.norm$samples <- results$metadata
redundant_features = results$redundant_features



results <- remove_large_factors(
    countdata.norm$samples,
    columns_to_ignore
)

countdata.norm$samples  <- results$metadata
large_factors = results$large_factors

results =  select_meta(countdata.norm$samples)
countdata.norm$samples <- results$metadata
filtered = results$filtered



log4r_info("Writing output!")
md = data.frame(id = dset_name, group=raw_metadata[raw_metadata$id == dset_name, ]$group, tissue=raw_metadata[raw_metadata$id == dset_name, ]$tissue,
     datetime = now(), raw_counts_genes = dim(raw_counts)[1], raw_counts_individuals= dim(raw_counts)[2],
    raw_metadata_individuals = dim(raw_meta)[1], raw_metadata_features = dim(raw_meta)[2],  filtered_counts_genes = dim(countdata.norm$counts)[1],
    filtered_counts_individuals = dim(countdata.norm$counts)[2], filtered_metadata_individuals = dim(countdata.norm$samples)[1],
    filtered_metadata_features = dim(countdata.norm$samples)[2],
    to_ignore_manually_curated = columns_to_ignore_scsv,
    to_ignore_redundant_features = paste(redundant_features, collapse=";"),
    to_ignore_large_factors =  paste(large_factors, collapse=";"),
    final_features_removed = paste(filtered,collapse=";") ,replicate_col=rep_col)

    log4r_info(dim(raw_metadata))
# output_raw_metadata_df <- rbind(output_raw_metadata_df, raw_metadata)
    write.table(md,  file = "raw_metadata_df.csv",col.names=!file.exists("raw_metadata_df.csv"), sep = ",",append=file.exists("raw_metadata_df.csv"), row.names=F)
}

# write.table(output_raw_metadata_df,  file = "raw_metadata_df.csv", append=TRUE,col.names=!file.exists("raw_metadata_df.csv"), sep = ",") # Save metadata to file

