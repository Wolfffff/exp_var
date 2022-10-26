library(ggthemes)
library(ggrepel)
library(here)
library(sjmisc)
library(pryr)
library(GOxploreR)
library(moments)
library(data.table)
library(enrichplot)
library("AnnotationDbi")
library("org.Hs.eg.db")

# Global -> Global

# %%
lcl_table <- read.table(here::here("data/10Jun21_eQTL_SNPs_sharing_EDITED.txt"), header=T, sep="\t")

library(stringr)
uniq_lcl_eqtls = unique(lcl_table$gene)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

rank_df = read.csv(here::here("data/pca_ranks.csv"), header = TRUE)[, -1]
# %%


# %%
library("org.Hs.eg.db") # remember to install it if you don't have it already
rank_df$Gene <- mapIds(org.Hs.eg.db, keys = rank_df$Gene, keytype = "ENSEMBL", column="SYMBOL")
# %%
organism = "org.Hs.eg.db"
metrics = c("mean", "sd")
#
df <- data.frame(tissue=character(), number_of_genes_tissue = integer(), number_of_genes_topbottom = integer(), total_count = integer(), upper_count = integer(), lower_count=integer(), stringsAsFactors=F)
tail_size = 0.05

upper_quantiles = list()
lower_quantiles = list()
metric = "sd"
cutoff = quantile(rank_df[[metric]],1 - tail_size)
subset = rank_df[rank_df[[metric]] >= cutoff,]$Gene
upper_quantiles[[metric]] = subset
cutoff = abs(quantile(rank_df[[metric]], tail_size))
lower_quantiles[[metric]] = rank_df[rank_df[[metric]] <= cutoff,]$Gene

metric = "sd"

rank_df_with_disease = rank_df %>% mutate(sd_disease = if_else(Gene %in% uniq_lcl_eqtls, 1, 0))
lower_with_disease = as.data.frame(lower_quantiles) %>% mutate(sd_disease = if_else(sd%in% uniq_lcl_eqtls, 1, 0))
upper_with_disease = as.data.frame(upper_quantiles) %>% mutate(sd_disease = if_else(sd%in% uniq_lcl_eqtls, 1, 0))

tmp_df <- data.frame(tissue="cross_tissue", number_of_genes_tissue = dim(as.data.frame(rank_df))[1], number_of_genes_topbottom = dim(as.data.frame(upper_quantiles))[1], total_count = length(uniq_lcl_eqtls), upper_count =sum(upper_with_disease$sd_disease), lower_count=sum(lower_with_disease$sd_disease), stringsAsFactors=F)
df = tmp_df
outputdir = here::here("data/annotation/lcl")

results_df = df
results_df$upper_or =(results_df$upper_count / (results_df$total_count - results_df$upper_count))/(5/95)
results_df$lower_or =(results_df$lower_count / (results_df$total_count - results_df$lower_count))/(5/95)

results_df$phyper_upper = phyper(results_df$upper_count, results_df$number_of_genes_topbottom, (results_df$number_of_genes_tissue - results_df$number_of_genes_topbottom), results_df$total_count)
results_df$phyper_lower = phyper(results_df$lower_count, results_df$number_of_genes_topbottom, (results_df$number_of_genes_tissue - results_df$number_of_genes_topbottom), results_df$total_count)
results_df$phyper_upper_bh_adj = p.adjust(results_df$phyper_upper, method = "BH")
results_df$phyper_lower_bh_adj = p.adjust(results_df$phyper_lower, method = "BH")


results_df[results_df$phyper_upper_bh_adj < 0.05, ] %>% dplyr::select(tissue, phyper_upper_bh_adj, phyper_upper_bh_adj)
results_df[results_df$phyper_lower_bh_adj < 0.05,] %>% dplyr::select(tissue, phyper_upper_bh_adj, phyper_lower_bh_adj)
