main_count_processing <- function(dset_name,
                                  exp_data,
                                  experimental_metadata,
                                  feature_vec,
                                  assay_name) {
  print(dset_name)
  dset <- exp_data[[dset_name]]

  # Catch empty columns_to_ignore
  columns_to_ignore = tryCatch(
    expr = {
      unlist(strsplit(experimental_metadata[experimental_metadata$id == dset_name, ]$columns_to_ignore, split = ";"))
    },
    error = function(e) {
      columns_to_ignore <- c("")
    }
  )
  columns_to_ignore = c(columns_to_ignore, crap_cols)

  metadata <- colData(dset)
  metadata <- data.frame(metadata)

  counts <- assays(dset)[[assay_name]]
  print(paste0("Unfiltered count dimensions: ", dim(counts)[1], " x ", dim(counts)[2]))
  print(paste0("Unfiltered metadata dimensions: ", dim(metadata)[1], " x ", dim(metadata)[2]))

  filtered_data <- make_filtered_data(counts, metadata, feature_vec)
  metadata <- filtered_data$metadata
  counts <- filtered_data$counts
  n_samples <- dim(metadata)[1]
  metadata$sample_id <- rownames(metadata)

  print("Normalizing and estimating mean-variance weights...")
  countdata.list <- DGEList(counts = counts, samples = metadata, genes = rownames(dset))
  rep_names = c("technical_replicate_group", "wells.replicate")
  if (any(names(metadata) %in% rep_names)) {
    print("Summing technical replicates...")
    rep_col = names(metadata)[names(metadata) %in% rep_names]
    countdata.list <- sumTechReps(countdata.list, metadata[[rep_col]])
  }
  countdata.list$samples <- remove_redundant_features(countdata.list$samples)

  # Removing columns with a crazy number of levels that mess everything up.
  # (this is why we have random effects by the way)
  countdata.list$samples <- remove_large_factors(
    countdata.list$samples,
    columns_to_ignore
  )
  print("Calculating normalization factors...")
  countdata.norm <- calcNormFactors(countdata.list)

  print("Trimming...")
  cutoff <- 1
  drop <- which(apply(cpm(countdata.norm), 1, max) < cutoff)
  countdata.norm <- countdata.norm[-drop, ]

  # recalc with updated metadata
  values_count <- sapply(lapply(countdata.norm$samples, unique), length)
  countdata.norm$samples <- countdata.norm$samples[, names(countdata.norm$samples[, values_count > 1])]

  print("Making design matrix...")
  design <- make_design_matrix(countdata.norm$samples, columns_to_ignore)
  print(paste("Design matrix size:", paste(dim(design), collapse = " x ")))
  print("Voom!")

  jpeg(paste0(plots_dir, dset_name, "_voom.jpg"))
  countdata.voom <- voom(countdata.norm, design = design, plot = T)
  dev.off()

  fit <- lmFit(countdata.voom, design)
  ebfit <- eBayes(fit)

  print("Corrected residual plot")
  jpeg(paste0(plots_dir, dset_name, "_eBayes_SA.jpg"))
  plotSA(ebfit, main="Final model: Mean-variance trend", ylab = "Sqrt( standard deviation )")
  dev.off()
  
  # Raw PCA plot
  design_intercept = model.matrix(~1, data = countdata.norm$samples)
  raw_voomed_resid = removeBatchEffect(countdata.voom, covariates = design_intercept)
  pca_on_voom <- pca_plot(raw_voomed_resid, color = rep(1, ncol(raw_voomed_resid)))
  screen_on_voom <- scree_plot(raw_voomed_resid)

  # Null model for SVA
  # mod0 = model.matrix(~1,data=countdata.list$samples)
  # We need to rebuild and ignore lib.size!

  # svobj = sva(countdata.norm$counts,mod0 = mod0,mod = design)

  # Batch effects
  # countdata_resids <- removeBatchEffect(cpm(countdata.voom, log=TRUE, prior.count=3), covariates = design)#removeBatchEffect(countdata.voom, covariates = design)
  countdata_resids <- ebfit#removeBatchEffect(ebfit, covariates = design)
  rownames(countdata_resids) <- countdata.voom$genes[, 1]

  # PCA plot with Batch effects
  pca_on_resids <- pca_plot(countdata_resids)
  scree_on_resids <- scree_plot(countdata_resids)
  
  print("Corrected residual plot")


  # Null model for SVA
  # mod0 = model.matrix(~1,data=countdata.list$samples)
  # We need to rebuild and ignore lib.size!

  # svobj = sva(countdata.norm$counts,mod0 = mod0,mod = design)



  # Dropping outlier
  rpca_resid <- PcaGrid(t(countdata_resids), 20, crit.pca.distances = 0.99)
  countdata.norm_noOut <- countdata.norm
  countdata.norm_noOut$counts <- countdata.norm_noOut$counts[, rpca_resid@flag]
  countdata.norm_noOut$samples <- countdata.norm_noOut$samples[rpca_resid@flag, ]

  # PCA plot with Batch effects
  pca_on_resids <- pca_plot(countdata_resids, color = !rpca_resid@flag)
  #scree_on_resids <- scree_plot(countdata_resids)


  print(paste0(
    "Filtered metadata dimensions: ",
    paste(dim(countdata.norm_noOut$samples), collapse = " x ")
  ))

  design_noOut <- make_design_matrix(countdata.norm_noOut$samples, columns_to_ignore)
  countdata.voom_noOut <- voom(countdata.norm_noOut, design = design_noOut)
  countdata_resids_noOut <- removeBatchEffect(countdata.voom_noOut, covariates = design_noOut)
  rownames(countdata_resids_noOut) <- countdata.voom_noOut$genes[, 1]
  pca_on_resids_noOut <- pca_plot(countdata_resids_noOut, color = rep(1, ncol(countdata_resids_noOut)))

  # Batch effects With PC1

  PCs <- pca(pca_on_resids_noOut)
  countdata.norm_noOut$samples$PC1 <- PCs$pc[1, ]
  design_with_pc1 <- make_design_matrix(countdata.norm_noOut$samples, columns_to_ignore)

  countdata_resids_with_pc1 <- removeBatchEffect(countdata.voom_noOut, covariates = design_with_pc1)
  rownames(countdata_resids_with_pc1) <- countdata.voom_noOut$genes[, 1]

  # PCA plot with Batch effects and PC1
  pca_on_resids_with_pc1 <- pca_plot(countdata_resids_with_pc1)

  # print("Writing figures")

  plt <- plot_grid(
    nrow = 2, scale = 0.9,
    pca_on_voom + ggtitle("Uncorrected"),
    pca_on_resids + ggtitle("Known effects"),
    pca_on_resids_noOut + ggtitle("Known effects, no outliers"),
    pca_on_resids_with_pc1 + ggtitle("Known effects and PC1")
  )

  # print("Appending results and metadata to lists")
  list(
    n_samples = dim(countdata.norm_noOut$samples)[1],
    raw_voom = countdata.voom,
    voom_noOut = countdata.voom_noOut,
    plotPanel = plt,
    plot_list = list(uncorrected = pca_on_voom, batch = pca_on_resids, clean = pca_on_resids_noOut, pc1 = pca_on_resids_with_pc1),
    residuals = countdata_resids_noOut,
    raw_residuals = countdata_resids,
    pc1_residuals = countdata_resids_with_pc1,
    metadata = countdata.norm_noOut$samples
  )
}