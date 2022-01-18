#pak::pkg_install(c("DESeq2", "vsn"))
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
library(org.Hs.eg.db)
library(rrcov)
library(DESeq2)
library(vsn)
library(viridis)
library(recount3)
library(tictoc)

# Set timeout to avoid failure when trying to download GTEx or other large datasets
options(timeout = 1800)

# Setting some awful metadata fields
crap_cols = c("alias", "Alias", "Broker.name", "broker.name", "Description", "Title", "ENA.checklist", 
              "ENA.FIRST.PUBLIC", "ENA.LAST.UPDATE", "isolate", "INSDC.center.alias", "sample_id",
              "INSDC.center.name", "INSDC.first.public", "INSDC.last.update", "INSDC.status", 
              "Sample.Name", "SRA.accession", "title", "gtex.smrin", "rownames", "tcga.xml_month_of_form_completion", 
              "tcga.xml_year_of_form_completion", "tcga.xml_year_of_initial_pathologic_diagnosis", 
              "tcga.xml_initial_pathologic_diagnosis_method", "tcga.cgc_case_year_of_diagnosis", 
              "tcga.gdc_metadata_files.file_size.analysis", "tcga.xml_breast_carcinoma_surgical_procedure_name", 
              "tcga.xml_day_of_form_completion", "tcga.cgc_sample_shortest_dimension", "tcga.xml_stage_event_system_version",
              "tcga.gdc_cases.samples.portions.analytes.concentration", "tcga.cgc_case_histological_diagnosis",
              "tcga.gdc_cases.samples.portions.analytes.aliquots.concentration")


# Setup filters for removing samples accordingly to the metadata
feature_vec <- list()
feature_vec[["disease"]] <- c("normal", "control", "", NA, 
                            "non inflammatory bowel disease control")
feature_vec[["treatment"]] <- c("normal", "control", "", NA)
feature_vec[["tcga.cgc_sample_sample_type"]] <- c("Solid Tissue Normal")
feature_vec[["diagnosis"]] <- c("Control")
feature_vec[["Healthy"]] <- c("Healthy")


downloadRecount3 <- function(id){
  # Load the project
  print(id)
  proj_info <- subset(
    human_projects,
    project == id & project_type == "data_sources"
  )
  print(proj_info)
  # Tape to deal with acquisition issues
  cache_path = paste("cache", id, sep = "/")
  dir.create(file.path(cache_path))
  local_cache <- recount3_cache(cache_dir = cache_path)
  rse <- create_rse(proj_info, bfc=local_cache)
  if (proj_info$file_source == "gtex") {
    metadata_df <- colData(rse)[,grepl("gtex",colnames(colData(rse)),fixed=TRUE)]
  } else if(proj_info$file_source == "tcga") {
    metadata_df <- colData(rse)[,grepl("tcga",colnames(colData(rse)),fixed=TRUE)]
  } else {
    #Could probably just use expand_sra_attributes
    metadata_df <- convert_metadata_to_df_rc3(colData(rse)$sra.sample_attributes, rse)
  }
  single_mask = sapply(metadata_df , function(x) length(table(x)) > 1) 
  metadata_df <- metadata_df[, single_mask, drop = FALSE]
  
  # Crude way to set NA to ""
  metadata_df@listData <- lapply(metadata_df@listData, function(x) {
    x[is.na(x)] = ""
    x
  })

  rse@colData <- metadata_df
  rse
}

select_meta = function(metadata){
    crap_cols = c(crap_cols, grep("^tcga..*_percent", names(metadata), value = TRUE))
    crap_cols = c(crap_cols, grep("^tcga..*analytes", names(metadata), value = TRUE))
    crap_cols = c(crap_cols, grep("^tcga..*xml_has_", names(metadata), value = TRUE))
    crap_cols = c(crap_cols, grep("^tcga..*_dimension$", names(metadata), value = TRUE))
    crap_cols = c(crap_cols, grep("^tcga..*tumor", names(metadata), value = TRUE))
    crap_cols = c(crap_cols, grep("^tcga..*cancer", names(metadata), value = TRUE))
    crap_cols = c(crap_cols, grep("^tcga..*pathologic", names(metadata), value = TRUE))
    crap_cols = c(crap_cols, grep("^tcga..*_collection_indicator$", names(metadata), value = TRUE))

    mask_no_levels = sapply(metadata, function(x) length(unique(x)) <= 1 | length(unique(x)) >= 20)
    mask_crap_cols = names(metadata) %in% crap_cols
    mask_dominant_level = laply(metadata, function(x) max(table(x)) > length(x)*(3/4))
    mask_small_levels = laply(metadata, function(x) sum(table(x) <= 2) >= 3)
    mask_high_missing = laply(metadata, function(x) sum(x == "")/length(x)) > 0.1
    mask = rowSums(cbind(mask_no_levels, mask_crap_cols, mask_high_missing, mask_dominant_level, mask_small_levels)) >= 1
    sel_metadata = metadata[,!mask, drop = FALSE]
    sel_metadata
}

pca <- function(x, space = c("rows", "columns"),
                center = TRUE, scale = FALSE) {
  space <- match.arg(space)
  if (space == "columns") {
    x <- t(x)
  }
  x <- t(scale(t(x), center = center, scale = scale))
  x <- x / sqrt(nrow(x) - 1)
  s <- svd(x)
  loading <- s$u
  colnames(loading) <- paste0("Loading", 1:ncol(loading))
  rownames(loading) <- rownames(x)
  pc <- diag(s$d) %*% t(s$v)
  rownames(pc) <- paste0("PC", 1:nrow(pc))
  colnames(pc) <- colnames(x)
  pve <- s$d^2 / sum(s$d^2)
  if (space == "columns") {
    pc <- t(pc)
    loading <- t(loading)
  }
  return(list(pc = pc, loading = loading, pve = pve))
}

make_filtered_data = function(counts, metadata, feature_vec){
  #print(colnames(metadata))
  filtered_metadata <- metadata
  filtered_counts <- as_data_frame(counts)
  for (column in names(feature_vec)) {
    if (column %in% colnames(filtered_metadata)){
      control_names_tb = table(filtered_metadata[,column], useNA = "ifany")
      control_names_tb = control_names_tb[names(control_names_tb) %in% feature_vec[[column]]]
      if(is.na(names(control_names_tb))){
        filtered_metadata[,column][is.na(filtered_metadata[,column])] = "control"
        names(control_names_tb) = "control"
      }
      control_name = names(control_names_tb)[control_names_tb == max(control_names_tb)]
      filtered_counts <- filtered_counts[,filtered_metadata[,column] == control_name]
      filtered_metadata <- filtered_metadata[filtered_metadata[,column] == control_name,]
    }
  }
  return(list(counts = filtered_counts, 
              metadata = filtered_metadata))
}

get_control_cols <- function(metadata, columns_to_ignore){
  names(metadata)[!(names(metadata) %in% c(columns_to_ignore,"lib.size","norm.factors"))]
}

remove_large_factors = function(metadata, columns_to_ignore){
  n_samples = nrow(metadata)
  cols_to_control = get_control_cols(metadata, columns_to_ignore)
  for(col in cols_to_control){
    n_coef = length(unique(metadata[[col]])) - 1
    #print(paste(col, n_coef))
    #cut_off = n_samples/5
    cut_off = 20
    if(n_coef > cut_off | n_coef <= 0){
      print(paste("Dropping column", col, "from metadata because it has", n_coef, "levels. Cut off is", cut_off, "."))
      metadata = metadata[,names(metadata) != col, drop=FALSE]
      cols_to_control = cols_to_control[cols_to_control!=col]
     }
  }
  return(metadata)
}

remove_redundant_features <- function(metadata){  
    # Removing redundant features - see https://stackoverflow.com/questions/38222318/how-to-remove-duplicate-columns-content-in-data-table-r

  features_pair <- combn(names(metadata), 2, simplify = F) # list all column pairs
  to_remove <- c() # init a vector to store duplicates
  for (pair in features_pair) { # put the pairs for testing into temp objects
    f1 <- pair[1]
    f2 <- pair[2]
    if (!(f1 %in% to_remove) & !(f2 %in% to_remove)) {
      if (all(as.numeric(as.factor(metadata[[f1]])) == as.numeric(as.factor(metadata[[f2]])))) {
        cat(f1, "and", f2, "are equal.\n")
        to_remove <- c(to_remove, f2) # build the list of duplicates
      }
    }
  }

  metadata <- metadata[,!(names(metadata) %in% as.list(to_remove))]
  metadata
}

pca_plot = function(resids, color = NULL){
  pca_resid = pca(resids)
  results = t(pca_resid$pc)
  mean = data.frame(matrix(colMeans(results), 
                    ncol = dim(results)[2]))
  colnames(mean) = colnames(results)
  if(is.null(color)){
    rpca_resid = PcaGrid(t(resids), 10, crit.pca.distances = 0.99)
    results = data.frame(results, outlier = !rpca_resid@flag)
  } else {
    results = data.frame(results, outlier = color)
  }
  colors = c("black",  "tomato3", viridis(length(unique(color))))
  plot <-  ggplot(results, aes(PC1, PC2, color = outlier)) + 
    geom_point() + 
    stat_ellipse(level = 0.99, color = "black") + 
    stat_ellipse(level = 0.99, type = "norm", color = "black", linetype = 2) + 
    geom_point(data= mean, color = "tomato3", size = 7) +
    coord_fixed(ratio=1) +
    labs(x = paste0("PC1 (", round(pca_resid$pve[1]*100, 2), "%)"),
         y = paste0("PC2 (", round(pca_resid$pve[2]*100, 2), "%)")) +
    theme_cowplot() + 
    scale_color_manual(values = colors[1:length(unique(color))]) +
    theme(legend.position = "none")
  plot
}


scree_plot <- function(resids){
  var_explained <- pca(resids)$pve
  scree_on_resids  <- qplot(c(1:10), var_explained[1:10]) +
      geom_line() +
      xlab("Principal Component") +
      ylab("Variance Explained") +
      ggtitle("Scree Plot") +
      ylim(0, 1) + 
      coord_fixed(ratio=1) +
      theme_cowplot()
}

vector_cor = function(x, y) x%*%y/sqrt(x%*%x * y%*%y)

make_design_matrix = function(metadata, columns_to_ignore){
  cols_to_control = get_control_cols(metadata, columns_to_ignore)
  if(length(cols_to_control) == 0) return(model.matrix(~1,data = metadata))
  b <- paste0(" ", cols_to_control, collapse=" +")
  model <- as.formula(paste0("~ 1 +",b))
  print(paste0("~ 1 +",b))
  design <- model.matrix(model, data = metadata)
  design = design[, qr(design)$pivot[seq_len(qr(design)$rank)], drop = FALSE]
}

map_to_cols_rc3 <- function(s) {
  s %>%
    str_split("\\|", simplify = TRUE) %>%
    str_split(";;", simplify = TRUE) %>%
    as_tibble() %>%
    t() %>%
    row_to_names(row_number = 1) %>%
    repair_names()
}

convert_metadata_to_df_rc3 <- function(sample_attributes, rse_local) {
  meta <- lapply(sample_attributes, FUN = map_to_cols_rc3)
  names(meta) <- rownames(colData(rse_local))
  meta_df <- dplyr::bind_rows(meta, .id = "rownames")
  meta_df <- DataFrame(meta_df)
  rownames(meta_df) <- meta_df$rownames
  return(meta_df)
}

DESeq2_vst_lm <- function(countdata_norm, design=NULL, label=NULL){
  counts <- countdata_norm$counts
  metadata <- countdata_norm$samples
  if(is.null(design)){
      design = matrix(1, nrow = ncol(counts), ncol = 1)
  }
   # Switch to DESeq2
  dds  <- DESeqDataSetFromMatrix(countData = counts, 
                                 colData = metadata, design = design)
  # Use vst wrapper for varianceStabilizingTransformation
  vsd <- vst(dds, blind = FALSE)
  if(!is.null(label)){
    jpeg(paste0(plots_dir, label, "_meanSd_vst.jpg"))
    meanSdPlot(assay(vsd))
    dev.off()
  }
  countdata_resids <- removeBatchEffect(assay(vsd), covariates = design)#removeBatchEffect(countdata.voom, covariates = design)
  
  rownames(countdata_resids) <- countdata_norm$genes[, 1]
  return(countdata_resids)
}

inv_log2_plus05 = function(x){
  return(2**x - 0.5)
}

cpm_lm <- function(countdata_norm, design=NULL){
  counts <- countdata_norm$counts
  metadata <- countdata_norm$samples
  if(is.null(design)){
      design = matrix(1, nrow = ncol(counts), ncol = 1)
  }
  cpm_counts  <- cpm(counts, log = TRUE, prior.count = 5)
  countdata_resids <- removeBatchEffect(cpm_counts, covariates = design)
  rownames(countdata_resids) <- countdata_norm$genes[, 1]
  return(countdata_resids)
}
