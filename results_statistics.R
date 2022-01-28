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

rank_df = read.csv(here::here("data/pca_ranks.csv"), header = TRUE)[, -1]


# %%
organism = "org.Hs.eg.db"

library(organism, character.only = TRUE)
library(clusterProfiler)
library(plyr)
library(dplyr)

rank_df = read.csv(here::here("data/pca_ranks.csv"), header = TRUE)[, -1]

top_quantiles = list()
metric_list = c("means","var","sd","cv")
for(metric in metric_list){
    cutoff = quantile(rank_df[[metric]], .95)
    subset = rank_df[rank_df[[metric]] >= cutoff,]$gene
    top_quantiles[[metric]] = subset
}
term2gene_df = ptwas_table[, c("Trait","Gene")]
ptwas_table_merged = merge(term2gene_df,ptwas_traits, by.x = "Trait", by.y = "ID", all.x = TRUE)
# term2gene_df = data.frame(disease="1",

ego = enrichGO(gene  = rank_df$gene,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05,
                         readable      = TRUE)
library(enrichplot)
png(here::here("eGO_universe.png"), height = 3840, width = 2160)
barplot(ego, showCategory=100) 
dev.off()