# Script to calculate enrichment of differentially expressed LCLs across environments in high and low expression genes


# %%
library(here)
library(tidyverse)
library(stringr)
library(AnnotationDbi)
# %%


# %%
organism <- "org.Hs.eg.db"
library(organism, character.only = TRUE)
# %%


# %%
# Load the data and convert Ensembl IDs to gene symbols for matching

# LCL enrichment table from https://github.com/AmandaJLea/LCLs_gene_exp/blob/main/main_results/10Jun21_eQTL_SNPs_sharing_EDITED.txt
lcl_table <- read.table(here::here("data/10Jun21_eQTL_SNPs_sharing_EDITED.txt"), header = T, sep = "\t")
uniq_lcl_eqtls <- unique(lcl_table$gene)
rank_df <- read.csv(here::here("data/pca_ranks.csv"), header = TRUE)[, -1]
rank_df$Gene <- mapIds(org.Hs.eg.db,
    keys = rank_df$Gene,
    keytype = "ENSEMBL", column = "SYMBOL"
)
# %%


# %%
# Select metric and take top and bottom n% of genes

metric <- "sd"
tail_size <- 0.05

upper_quantiles <- list()
lower_quantiles <- list()

cutoff <- quantile(rank_df[[metric]], 1 - tail_size)
subset <- rank_df[rank_df[[metric]] >= cutoff, ]$Gene
upper_quantiles[[metric]] <- subset

cutoff <- abs(quantile(rank_df[[metric]], tail_size))
lower_quantiles[[metric]] <- rank_df[rank_df[[metric]] <= cutoff, ]$Gene
# %%


# %%
# Find overlap between LCL eQTLs and top/bottom n% of genes, calculate odds ratio, and perform hypergeometric test

# Create dummy var noting overlap with top/bottom n% of genes
rank_df_with_dummyvar <- rank_df %>% mutate(sd_dummyvar = if_else(Gene %in% uniq_lcl_eqtls, 1, 0))
lower_with_dummyvar <- as.data.frame(lower_quantiles) %>% mutate(sd_dummyvar = if_else(sd %in% uniq_lcl_eqtls, 1, 0))
upper_with_dummyvar <- as.data.frame(upper_quantiles) %>% mutate(sd_dummyvar = if_else(sd %in% uniq_lcl_eqtls, 1, 0))

results_df <- data.frame(tissue = "cross_tissue", number_of_genes_tissue = dim(as.data.frame(rank_df))[1], number_of_genes_topbottom = dim(as.data.frame(upper_quantiles))[1], total_count = length(uniq_lcl_eqtls), upper_count = sum(upper_with_dummyvar$sd_dummyvar), lower_count = sum(lower_with_dummyvar$sd_dummyvar), stringsAsFactors = F)

# Calculate odds ratios
results_df$upper_or <- (results_df$upper_count / (results_df$total_count - results_df$upper_count)) / (int(100*tail_size) / (100 - int(100*tail_size)))
results_df$lower_or <- (results_df$lower_count / (results_df$total_count - results_df$lower_count)) / (int(100*tail_size) / (100 - int(100*tail_size)))

# Perform hypergeometric test and update results_df
results_df$phyper_upper <- phyper(results_df$upper_count, results_df$number_of_genes_topbottom, (results_df$number_of_genes_tissue - results_df$number_of_genes_topbottom), results_df$total_count)
results_df$phyper_lower <- phyper(results_df$lower_count, results_df$number_of_genes_topbottom, (results_df$number_of_genes_tissue - results_df$number_of_genes_topbottom), results_df$total_count)
results_df$phyper_upper_bh_adj <- p.adjust(results_df$phyper_upper, method = "BH")
results_df$phyper_lower_bh_adj <- p.adjust(results_df$phyper_lower, method = "BH")

write_csv(data.frame(results_df), here::here("data/annotation/lcl_output_table.csv"))
# %%
