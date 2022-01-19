my_logfile = snakemake@log[[1]]
snakemake@source("logger.R")
log4r_info("Starting.")
print = log4r_info
log4r_info("Loading packages") 
source("../functions.R")

dset_name = snakemake@wildcards[["id"]]

log4r_info("Reading data")
dset <- readRDS(snakemake@input[["data"]])

countdata.norm = dset$data
columns_to_ignore = dset$columns_to_ignore

log4r_info("Making design matrix")
design <- make_design_matrix(countdata.norm$samples, columns_to_ignore)
print(paste("Design matrix size:", paste(dim(design), collapse = " x ")))

pca_on_raw <- pca_plot(countdata.norm$counts, color = rep("1", ncol(countdata.norm$counts)))

# Switch to DESeq2

print(paste0("Filtered count dimensions: ", 
             dim(countdata.norm$counts)[1], " x ", dim(countdata.norm$counts)[2]))
countdata_resids <- DESeq2_vst_lm(countdata.norm, design = design, label = dset_name)#removeBatchEffect(countdata.voom, covariates = design)

png(snakemake@output[["mean_var"]])
    meanSdPlot(countdata_resids)
dev.off()

# Dropping outlier
countdata.norm_noOut <- countdata.norm

rpca_resid <- PcaGrid(t(countdata_resids), 20, crit.pca.distances = 0.99)
countdata.norm_noOut$counts <- countdata.norm_noOut$counts[, rpca_resid@flag]
countdata.norm_noOut$samples <- countdata.norm_noOut$samples[rpca_resid@flag, ,drop = FALSE]

# PCA plot with Batch effects (this plot happens here to make use of the outlier tags from the robust PCA)
pca_on_resids <- pca_plot(countdata_resids, color = !rpca_resid@flag)
#scree_on_resids <- scree_plot(countdata_resids)

print(paste0(
"Filtered metadata dimensions: ",
paste(dim(countdata.norm_noOut$samples), collapse = " x ")
))

design_noOut <- make_design_matrix(countdata.norm_noOut$samples, columns_to_ignore)
countdata_resids_noOut <- DESeq2_vst_lm(countdata.norm_noOut, design = design_noOut)
pca_on_resids_noOut <- pca_plot(countdata_resids_noOut, color = rep("1", ncol(countdata_resids_noOut)))

# Batch effects With PC1

PCs <- pca(countdata_resids_noOut)
countdata.norm_noOut$samples$PC1 <- PCs$pc[1, ]
design_with_pc1 <- make_design_matrix(countdata.norm_noOut$samples, columns_to_ignore)

countdata_resids_with_pc1 <- DESeq2_vst_lm(countdata.norm_noOut, design = design_with_pc1)

# PCA plot with Batch effects and PC1
pca_on_resids_with_pc1 <- pca_plot(countdata_resids_with_pc1, color = rep("1", ncol(countdata_resids_with_pc1)))

print("Writing figures")

plt <- plot_grid(nrow = 2, scale = 0.9,
                pca_on_raw + ggtitle("Uncorrected"),
                pca_on_resids + ggtitle("Known effects"),
                pca_on_resids_noOut + ggtitle("Known effects, no outliers"),
                pca_on_resids_with_pc1 + ggtitle("Known effects and PC1")
) 
save_plot(filename = snakemake@output[["pca_panel"]], plt,
          base_height = 6, base_asp = 1.2, ncol = 2, nrow = 2)

plot_list = list(uncorrected = pca_on_raw, 
                 batch = pca_on_resids, 
                 clean = pca_on_resids_noOut, 
                 pc1 = pca_on_resids_with_pc1)
saveRDS(plot_list, snakemake@output[["all_plots"]])


print("Saving data")
output = list(
    name = dset_name,
    n_samples = dim(countdata.norm_noOut$samples)[1],
    normcounts_raw = countdata.norm,
    normcounts_noOut = countdata.norm_noOut,
    residuals_noOut = countdata_resids_noOut,
    residuals_raw = countdata_resids,
    residuals_pc1 = countdata_resids_with_pc1,
    metadata = countdata.norm_noOut$samples,
    design = design_with_pc1
)
saveRDS(output, file = snakemake@output[["residuals"]])

log4r_info("Done!")