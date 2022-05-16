#pak::pkg_install(c("DESeq2", "vsn"))
library(ExpressionAtlas)
library(plyr)
library(tidyverse)
library(limma)
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
library(iMKT)
library(wesanderson)


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
feature_vec[["source_name"]] <- c("normal_skin","healthy children without bacterial colonization")
feature_vec[["condition"]] <- c("Control","")
feature_vec[["PTSD"]] <- c("No")
feature_vec[["ptsd"]] <- c("Never")
feature_vec[["infections.agent"]] = c("n/a","")
feature_vec[["disease.state"]] = c("healthy")
feature_vec[["cancer.type"]] = c("HC")
feature_vec[["group"]] = c("Control", NA)
feature_vec[["time"]] = c('Convalescent')


downloadRecount3 <- function(id){
  # Load the project

  # Tape to deal with acquisition issues
  cache_path = paste("cache", id, sep = "/")
  print(cache_path)
  if(!dir.exists(cache_path))
    dir.create(cache_path, showWarnings = FALSE)
  local_cache <- recount3_cache(cache_dir = cache_path)
  human_projects <- available_projects(bfc = local_cache)
  proj_info <- subset(
    human_projects,
    project == id & project_type == "data_sources"
  )
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
  filtered_metadata <- metadata
  filtered_counts <- as_data_frame(counts)
  for (column in names(feature_vec)) {
    if (column %in% colnames(filtered_metadata)){
      control_names_tb = table(filtered_metadata[,column], useNA = "ifany")
      control_names_tb = control_names_tb[names(control_names_tb) %in% feature_vec[[column]]]
      if(length(control_names_tb) == 0){
        next
      }
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
  list(metadata = metadata, redundant_features = to_remove)
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
  model <- as.formula(paste0("~ 0 +",b))
  print(paste0("~ 0 +",b))
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

DESeq2_vst_lm <- function(countdata_norm, design=NULL, label=NULL, plot_file=NULL){
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
    if(is.null(plot_file)) plot_file = paste0(plots_dir, label, "_meanSd_vst.png")
    png(plot_file)
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

calculate_row_wise_metric <- function(results_list,f){
  summarized_list = vector("list", length = length(results_list))
  names(summarized_list) = names(results_list)
  for (dset_name in names(results_list)) {
    exprDf = results_list[[dset_name]]
    gene_vars = f(exprDf$residuals_noOut)
    summarized_list[[dset_name]] = data.frame(
        Genes = row.names(exprDf$residuals_noOut),
        var = gene_vars) %>%
        tidyr::separate("Genes", c("Genes", NA))
  }
  summarized_df = purrr::reduce(summarized_list, inner_join, by = "Genes")
  colnames(summarized_df)[-1] = names(results_list)
  return(summarized_df)
}

calculate_row_wise_metric_sparse <- function(results_list, f, max_missingness = 0.50){
  summarized_list = vector("list", length = length(results_list))
  names(summarized_list) = names(results_list)
  n_dsets = length(results_list)
  for (dset_name in names(results_list)) {
    exprDf = results_list[[dset_name]]
    gene_vars = f(exprDf$residuals_noOut)
    summarized_list[[dset_name]] = data.frame(
        Genes = row.names(exprDf$residuals_noOut),
        var = gene_vars) %>%
        tidyr::separate("Genes", c("Genes", NA))
  }
  summarized_df = purrr::reduce(summarized_list, full_join, by = "Genes")
  colnames(summarized_df)[-1] = names(results_list)
  x = summarized_df[1,-1]
  rowmask = apply(summarized_df[,-1], 1, function(x) sum(is.na(x))/n_dsets < max_missingness)
  summarized_df = summarized_df[rowmask,]
  return(summarized_df)
}

PopHumanAnalysis_modified <- function(genes=c("gene1","gene2","..."), pops=c("pop1","pop2","..."), cutoffs=c(0,0.05,0.1), recomb=TRUE/FALSE, bins=0, test=c("standardMKT","DGRP","FWW","asymptoticMKT","iMKT"), xlow=0, xhigh=1, plot=FALSE) {

  ## Get PopHuman data
  if (exists("PopHumanData") == TRUE) {
    data <- get("PopHumanData")
  } else {
    loadPopHuman()
    data <- get("PopHumanData") }

  ## Check input variables
  ## Numer of arguments
  if (nargs() < 3 && nargs()) {
    stop("You must specify 3 arguments at least: genes, pops, recomb (T/F).\nIf test = asymptoticMKT or test = iMKT, you must specify xlow and xhigh values.") }

  ## Argument genes
  if (length(genes) == 0 || genes == "" || !is.character(genes)) {
    stop("You must specify at least one gene.") }
  if (!all(genes %in% data$ID) == TRUE) {
    difGenes <- setdiff(genes, data$ID)
    difGenes <- paste(difGenes, collapse=", ")
    stopMssg <- paste0("MKT data is not available for the requested gene(s).\nRemember to use gene IDs from Ensembl (ENSG...).\nThe genes that caused the error are: ", difGenes, ".")
    stop(stopMssg) }

  ## Argument pops
  if (length(pops) == 0 || pops == "" || !is.character(pops)) {
    stop("You must specify at least one pop.") }
  if (!all(pops %in% data$pop) == TRUE) {
    correctPops <- c("ACB","ASW","BEB","CDX","CEU","CHB","CHS","CLM","ESN","FIN","GBR","GIH","GWD","IBS","ITU","JPT","KHV","LWK","MSL","MXL","PEL","PJL","PUR","STU","TSI","YRI")
    difPops <- setdiff(pops, correctPops)
    difPops <- paste(difPops, collapse=", ")
    stopMssg <- paste0("MKT data is not available for the sequested pops(s).\nSelect among the following pops:\nACB, ASW, BEB, CDX, CEU, CHB, CHS, CLM, ESN, FIN, GBR, GIH, GWD, IBS, ITU, JPT, KHV, LWK, MSL, MXL, PEL, PJL, PUR, STU, TSI, YRI!.\nThe pops that caused the error are: ", difPops, ".")
    stop(stopMssg) }

  ## Argument recomb
  if (recomb != TRUE && recomb != FALSE) {
    stop("Parameter recomb must be TRUE or FALSE.") }

  ## Argument bins
  if (recomb == TRUE) {
    if (!is.numeric(bins) || bins == 0  || bins == 1) {
      stop("If recomb = TRUE, you must specify the number of bins to use (> 1).") }
    if (bins > round(length(genes)/2)) {
      stop("Parameter bins > (genes/2). At least 2 genes for each bin are required.") }
  }
  if (recomb == FALSE && bins != 0) {
    warning("Parameter bins not used! (recomb=F selected)") }

  ## Argument test and xlow + xhigh (when necessary)
  if(missing(test)) {
    test <- "standardMKT"
  }
  else if (test != "standardMKT" && test != "DGRP" && test != "FWW" && test != "asymptoticMKT" && test != "iMKT") {
    stop("Parameter test must be one of the following: standardMKT, DGRP, FWW, asymptoticMKT, iMKT")
  }
  if (length(test) > 1) {
    stop("Select only one of the following tests to perform: standardMKT, DGRP, FWW, asymptoticMKT, iMKT") }
  if ((test == "standardMKT" || test == "DGRP" || test == "FWW") && (xlow != 0 || xhigh != 1)) {
    warningMssgTest <- paste0("Parameters xlow and xhigh not used! (test = ",test," selected)")
    warning(warningMssgTest) }

  ## Arguments xlow, xhigh features (numeric, bounds...) checked in checkInput()

  ## Perform subset
  subsetGenes <- data[(data$ID %in% genes & data$pop %in% pops), ]
  subsetGenes$ID <- as.factor(subsetGenes$ID)
  subsetGenes <- droplevels(subsetGenes)

  ## If recomb analysis is selected
  if (recomb == TRUE) {

    ## Declare output list (each element 1 pop)
    outputList <- list()

    for (k in levels(subsetGenes$pop)) {
      print(paste0("pop = ", k))

      ## Declare bins output list (each element 1 bin)
      outputListBins <- list()

      x <- subsetGenes[subsetGenes$pop == k, ]
      x <- x[order(x$recomb), ]

      ## create bins
      binsize <- round(nrow(x)/bins) ## Number of genes for each bin
      count <- 1
      x$Group <- ""
      dat <- x[FALSE, ] ## Create df with colnames

      for (i in 0:nrow(x)) {
        if (i%%binsize == 0) { ## Only if reminder of division = 0 (equally sized bins)
          i1 <- i + binsize
          if (i == 0) {
            g1 <- x[i:binsize,]
            group <- count
            g1$Group <- group
            dat[i:binsize,] <- g1
            count <- count+1 }
          else if (i1 <= nrow(x)) {
            ii <- i+1
            g1 <- x[ii:i1,]
            group <- count
            g1$Group <- group
            dat[ii:i1,] <- g1
            count <- count+1 }
        }
      }
      dat$Group <- as.factor(dat$Group)

      ## Iterate through each recomb bin
      for (j in levels(dat$Group)) {
        print(paste0("Recombination bin = ", j))
        x1 <- dat[dat$Group == j, ]

        ## Recomb stats from bin j
        numGenes <- nrow(x1)
        minRecomb <- min(x1$recomb, na.rm=T)
        medianRecomb <- median(x1$recomb, na.rm=T)
        meanRecomb <- mean(x1$recomb, na.rm=T)
        maxRecomb <- max(x1$recomb, na.rm=T)
        recStats <- cbind(numGenes,minRecomb,medianRecomb,meanRecomb,maxRecomb)
        recStats <- as.data.frame(recStats)
        recStats <- list("Recombination bin Summary"=recStats)

        ## Set counters to 0
        Pi <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
        P0 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
        f <- seq(0.025,0.975,0.05)
        mi <- 0; m0 <- 0
        Di <- 0; D0 <- 0

        ## Group genes
        x1 <- droplevels(x1)
        for (l in levels(x1$ID)) {
          x2 <- x1[x1$ID == l, ]

          ## DAF
          x2$DAF0f <- as.character(x2$DAF0f); x2$DAF4f <- as.character(x2$DAF4f)
          daf0f <- unlist(strsplit(x2$DAF0f, split=";"))
          daf4f <- unlist(strsplit(x2$DAF4f, split=";"))
          daf0f <- as.numeric(daf0f); daf4f <- as.numeric(daf4f)
          Pi <- Pi + daf0f; P0 <- P0 + daf4f

          ## Divergence
          mi <- mi + x2$mi; m0 <- m0 + x2$m0
          Di <- Di + x2$di; D0 <- D0 + x2$d0
        }

        ## Proper formats
        daf <- cbind(f, Pi, P0); daf <- as.data.frame(daf)
        names(daf) <- c("daf","Pi","P0")
        div <- cbind(mi, Di, m0, D0); div <- as.data.frame(div)
        names(div) <- c("mi","Di","m0","D0")

        ## Check data inside each test!

        ## Transform daf20 to daf10 (faster fitting) for asymptoticMKT and iMKT
        if (nrow(daf) == 20) {
          daf1 <- daf
          daf1$daf10 <- sort(rep(seq(0.05,0.95,0.1),2)) ## Add column with the daf10 frequencies
          daf1 <- daf1[c("daf10","Pi","P0")] ## Keep new frequencies, Pi and P0
          daf1 <- aggregate(. ~ daf10, data=daf1, FUN=sum)  ## Sum Pi and P0 two by two based on daf
          colnames(daf1)<-c("daf","Pi","P0") ## Final daf columns name
        }

        ## Perform test
        if(test == "standardMKT") {
          output <- standardMKT(daf, div)
          output <- c(output, recStats) } ## Add recomb summary for bin j
        else if(test == "DGRP" && plot == FALSE) {
          output <- DGRP(daf, div, listCutoffs=cutoffs)
          output <- c(output, recStats) }
        else if(test == "DGRP" && plot == TRUE) {
          output <- DGRP(daf, div, listCutoffs=cutoffs, plot=TRUE)
          output <- c(output, recStats) }
        else if(test == "FWW" && plot == FALSE) {
          output <- FWW(daf, div, listCutoffs=cutoffs)
          output <- c(output, recStats) }
        else if(test == "FWW" && plot == TRUE) {
          output <- FWW(daf, div, listCutoffs=cutoffs, plot=TRUE)
          output <- c(output, recStats) }
        else if(test == "asymptoticMKT") {
          output <- asymptoticMKT(daf1, div, xlow, xhigh)
          output <- c(output, recStats) }
        else if(test == "iMKT" && plot == FALSE) {
          output <- iMKT(daf1, div, xlow, xhigh)
          output <- c(output, recStats) }
        else if(test == "iMKT" && plot == TRUE) {
          output <- iMKT(daf1, div, xlow, xhigh, plot=TRUE)
          output <- c(output, recStats) }

        ## Fill list with each bin
        outputListBins[[paste("Recombination bin = ",j)]] <- output
      }

      ## Fill list with each pop
      outputList[[paste("pop = ",k)]] <- outputListBins
    }

    ## Warning if some genes are lost. Bins must be equally sized.
    if (nrow(dat) != length(genes)) {
      missingGenes <- round(length(genes) - nrow(dat))
      genesNames <- as.vector(tail(x, missingGenes)$ID)
      genesNames <- paste(genesNames, collapse=", ")
      warningMssg <- paste0("The ",missingGenes," gene(s) with highest recombination rate estimates (", genesNames, ") was/were excluded from the analysis in order to get equally sized bins.\n")
      warning(warningMssg) }

    ## Return output
    cat("\n")
    return(outputList)
  }

  ## If NO recombination analysis selected
  else if (recomb == FALSE) {

    ## Declare output list (each element 1 pop)
    outputList <- list()

    for (i in levels(as.factor(subsetGenes$pop))) {
      print(paste0("pop = ", i))
      x <- subsetGenes[subsetGenes$pop == i, ]

      ## Set counters to 0
      Pi <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
      P0 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
      f <- seq(0.025,0.975,0.05)
      mi <- 0; m0 <- 0
      Di <- 0; D0 <- 0

      ## Group genes
      for (j in levels(x$ID)) {
        x1 <- x[x$ID == j, ]

        ## DAF
        x1$DAF0f <- as.character(x1$DAF0f); x1$DAF4f <- as.character(x1$DAF4f)
        daf0f <- unlist(strsplit(x1$DAF0f, split=";"))
        daf4f <- unlist(strsplit(x1$DAF4f, split=";"))
        daf0f <- as.numeric(daf0f); daf4f <- as.numeric(daf4f)
        Pi <- Pi + daf0f; P0 <- P0 + daf4f

        ## Divergence
        mi <- mi + x1$mi; m0 <- m0 + x1$m0
        Di <- Di + x1$di; D0 <- D0 + x1$d0
      }

      ## Proper formats
      daf <- cbind(f, Pi, P0); daf <- as.data.frame(daf)
      names(daf) <- c("daf","Pi","P0")
      div <- cbind(mi, Di, m0, D0); div <- as.data.frame(div)
      names(div) <- c("mi","Di","m0","D0")

      ## Check data inside each test!

      ## Transform daf20 to daf10 (faster fitting) for asymptoticMKT and iMKT
      if (nrow(daf) == 20) {
        daf1 <- daf
        daf1$daf10 <- sort(rep(seq(0.05,0.95,0.1),2)) ## Add column with the daf10 frequencies
        daf1 <- daf1[c("daf10","Pi","P0")] ## Keep new frequencies, Pi and P0
        daf1 <- aggregate(. ~ daf10, data=daf1, FUN=sum)  ## Sum Pi and P0 two by two based on daf
        colnames(daf1)<-c("daf","Pi","P0") ## Final daf columns name
      }

      ## Perform test
      if(test == "standardMKT") {
        output <- standardMKT(daf, div) }
      else if(test == "DGRP" && plot == FALSE) {
        output <- DGRP(daf, div, listCutoffs=cutoffs) }
      else if(test == "DGRP" && plot == TRUE) {
        output <- DGRP(daf, div, listCutoffs=cutoffs, plot=TRUE) }
      else if(test == "FWW" && plot == FALSE) {
        output <- FWW(daf, div, listCutoffs=cutoffs) }
      else if(test == "FWW" && plot == TRUE) {
        output <- FWW(daf, div, listCutoffs=cutoffs, plot=TRUE) }
      else if(test == "asymptoticMKT") {
        output <- asymptoticMKT(daf1, div, xlow, xhigh) }
      else if(test == "iMKT" && plot == FALSE) {
        output <- iMKT(daf1, div, xlow, xhigh) }
      else if(test == "iMKT" && plot == TRUE) {
        output <- iMKT(daf1, div, xlow, xhigh, plot=TRUE) }

      ## Fill list with each pop
      outputList[[paste("pop = ",i)]] <- output
    }

    ## Return output
    cat("\n")
    return(outputList)
  }
}

remove_id_ver <- function(x){
  dot_split = function(x)strsplit(x, "\\.")[[1]][1]
  sapply(x, dot_split, USE.NAMES = FALSE)
}

quantile_violin_plot = function(x, y, ntiles = 10){
    df = data.frame(x=x,y=y)
    df = df %>% mutate(quantilegroup = ntile(x, ntiles))
    ggplot(df, aes(x = factor(quantilegroup), group = quantilegroup, y)) +
      geom_violin(fill = wes_palette("Royal2")[5], alpha = 0.5) +
      scale_x_discrete(labels=c(1:ntiles))
}

