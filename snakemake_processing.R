source(here::here("functions.R"))
library(here)

file_paths <- list.files(path = here::here("snakemake/Rdatas/residuals/"), 
                         pattern = "\\.rds", full.names = TRUE)
file_names <-  gsub(pattern = "\\.rds$", replacement = "", x = basename(file_paths))

library(plyr)
library(doMC)
registerDoMC(64)
data_list <- llply(file_paths, readRDS, .parallel = TRUE)

names(data_list) <- file_names  

# Process row wise metrics
library(dplyr)
library(Rfast)
library(matrixStats)

gene_metric_dfs <- list(
    means = calculate_row_wise_metric(data_list,rowMeans),
    var = calculate_row_wise_metric(data_list,rowVars),
    sd = calculate_row_wise_metric(data_list,rowSds),
    cv = calculate_row_wise_metric(data_list,rowcvs)
)

saveRDS(gene_metric_dfs,  file = here::here("snakemake/Rdatas/mean_variance.RDS"))

pak::pkg_install("sergihervas/iMKT")
library(iMKT)
loadPopHuman()

# Divergence from chimpanzees -- where available.
imkt_results = NULL
for (gene in unique(PopHumanData$ID)){
  columns_to_ignore = tryCatch(
    expr = {
      results = PopHumanAnalysis_modified(genes=gene, pops=c("CEU"), recomb=FALSE, test="standardMKT", plot=TRUE)
      results_flat = unlist(r$`pop =  CEU`)
      results_flat["gene"] <- gene
      names(results_flat) <- make.names(names(results_flat), unique=TRUE)
      results_df = t(as.data.frame(results_flat))
      rownames(results_df) = gene
      imkt_results <- rbind(imkt_results, results_df)
    },
    error = function(e) {
        print(paste0("Skipped gene: ", gene))
        print(e)
    }
  )
}
write.csv(imkt_results, file="imkt_results.csv")


