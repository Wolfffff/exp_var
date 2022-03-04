# %%
source(here::here("functions.R"))
library(here)
library(sjmisc)
library(pryr)

# %%

# %%
# Read in the data for each chunk in parallel if needed. 
file_paths <- list.files(path = here::here("snakemake/Rdatas/residuals/"), 
                         pattern = "\\.rds", full.names = TRUE)
file_names <-  gsub(pattern = "\\.rds$", replacement = "", x = basename(file_paths))

library(plyr)
library(doMC)
registerDoMC(64)
data_list <- llply(file_paths, readRDS, .parallel = TRUE)

names(data_list) <- file_names

residuals_noout <- list()
for (name in names(data_list)){
    residuals_noout[[name]] <- as.data.frame(data_list[[name]]$residuals_noOut) %>% rotate_df()
}

# %%
# %%


residuals_df <-bind_rows(residuals_noout,.id = 'source')
dim(residuals_df)
n_dsets = length(data_list)
source = residuals_df$source

residuals_array = as.matrix(residuals_df[,-1])
residuals_array <- Rfast::transpose(residuals_array)


residuals_df <- residuals_df[1:10, 1:10] %>% rotate_df()
residuals_df$genes = rownames(residuals_df)
summarized_df <- residuals_df
# summarized_df = purrr::reduce(summarized_list, dplyr::full_join, by = "Genes")
# colnames(summarized_df)[-1] = names(results_list)
# x = summarized_df[1,-(1:2)]  
max_missingness =0.5
rowmask = apply(summarized_df[,-c(1:2)], 1, function(x) sum(is.na(x))/n_dsets < max_missingness) 
summarized_df = summarized_df[rowmask,]
# %%

# %%
# GO analysis on ea
organism = "org.Hs.eg.db"
metrics = c("mean", "sd")

library(organism, character.only = TRUE)
library(clusterProfiler)
library(plyr)
library(dplyr)

rank_df = read.csv(here::here("data/pca_ranks.csv"), header = TRUE)[, -1]

upper_quantiles = list()
lower_quantiles = list()
for(metric in c("mean","sd")){
    cutoff = quantile(rank_df[[metric]], .95)
    subset = rank_df[rank_df[[metric]] >= cutoff,]$Gene
    upper_quantiles[[metric]] = subset
    cutoff = abs(quantile(-rank_df[[metric]], .95))
    lower_quantiles[[metric]] = rank_df[rank_df[[metric]] <= quantile(rank_df[[metric]], .05),]$Gene
}
# %%

# %%
# GO analysis against our list of genes
metric = "sd"
library(enrichplot)

local_go_upper = enrichGO(gene  = upper_quantiles[[metric]],
                          universe = rank_df$Gene,
                          OrgDb         = org.Hs.eg.db,
                          keyType       = 'ENSEMBL',
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05,
                          readable      = TRUE)
png(here::here("data/plots/local_go_upper.png"), height = 3840, width = 2160)
barplot(local_go_upper, showCategory=100) 
dev.off()

local_go_lower = enrichGO(gene  = lower_quantiles[[metric]],
                          universe = rank_df$Gene,
                          OrgDb         = org.Hs.eg.db,
                          keyType       = 'ENSEMBL',
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05,
                          readable      = TRUE)
png(here::here("data/plots/local_go_lower.png"), height = 3840, width = 2160)
barplot(local_go_lower , showCategory=100) 
dev.off()
# %%




# %%
file_paths <- list.files(path = here::here("snakemake/Rdatas/residuals/"), 
                         pattern = "\\.rds", full.names = TRUE)
file_names <-  gsub(pattern = "\\.rds$", replacement = "", x = basename(file_paths))

library(plyr)
library(doMC)
registerDoMC(64)
data_list <- llply(file_paths, readRDS, .parallel = TRUE)

names(data_list) <- file_names

# %%

# %%
library("biomaRt")
# listMarts()
mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")



# %%

metrics <- readRDS("snakemake/Rdatas/gene_metrics.RDS")
