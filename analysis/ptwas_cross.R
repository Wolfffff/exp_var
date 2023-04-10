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
ptwas_table <- read.table(here::here("data/annotation/ptwas_table3_sig_genetrait_interactions.txt"), header=T, sep=" ")
ptwas_traits <- read.csv(here::here("data/annotation/ptwas_traits.csv"))
ptwas_metadata <- read.csv(here::here("data/annotation/ptwas_metadata.csv"))

library(stringr)
ptwas_table$Gene =  str_split_fixed(ptwas_table$Gene,'\\.',Inf)[,1]
uniq_disease_linked_genes = unique(ptwas_table$Gene)

organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

rank_df = read.csv(here::here("data/pca_ranks.csv"), header = TRUE)[, -1]
# %%


# %%
library(plyr)
library(dplyr)

term2gene_df = ptwas_table[, c("Trait","Gene")]
ptwas_table_merged = merge(term2gene_df,ptwas_traits, by.x = "Trait", by.y = "ID", all.x = TRUE)
# %%

# %%
organism = "org.Hs.eg.db"
metrics = c("mean", "sd")

df <- data.frame(tissue=character(), ptwas_category = character(), number_of_genes_tissue = integer(), number_of_genes_topbottom = integer(), category_count=integer(), total_count = integer(), upper_count = integer(), lower_count=integer(), stringsAsFactors=F)
tail_size = 0.05

upper_quantiles = list()
lower_quantiles = list()
for(metric in c("mean","sd")){
    cutoff = quantile(rank_df[[metric]],1 - tail_size)
    subset = rank_df[rank_df[[metric]] >= cutoff,]$Gene
    upper_quantiles[[metric]] = subset
    cutoff = abs(quantile(rank_df[[metric]], tail_size))
    lower_quantiles[[metric]] = rank_df[rank_df[[metric]] <= cutoff,]$Gene
}
# %%

# %%
# GO analysis against our list of genes
metric = "sd"

# outputdir = here::here("data/annotation/per_tissue_ptwas", current_tissue)
# dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)
for (cat in unique(ptwas_table_merged$Category)){
    if(is.na(cat)){
        next
    }
    cat_df = ptwas_table_merged[ptwas_table_merged$Category == cat,]
    rank_df_with_disease = rank_df %>% mutate(sd_disease = if_else(Gene %in% cat_df$Gene, 1, 0))
    lower_with_disease = as.data.frame(lower_quantiles) %>% mutate(sd_disease = if_else(sd%in% cat_df$Gene, 1, 0))
    upper_with_disease = as.data.frame(upper_quantiles) %>% mutate(sd_disease = if_else(sd%in% cat_df$Gene, 1, 0))
    # print()
    tmp_df <- data.frame(tissue="cross_tissue" , ptwas_category = cat, number_of_genes_tissue = dim(as.data.frame(rank_df))[1], number_of_genes_topbottom = dim(as.data.frame(upper_quantiles))[1],
    category_count=length(unique(cat_df$Gene)), total_count = sum(rank_df_with_disease$sd_disease), upper_count =sum(upper_with_disease$sd_disease), lower_count=sum(lower_with_disease$sd_disease), stringsAsFactors=F)
    df = rbind(df,tmp_df)
}
outputdir = here::here("data/annotation/cross_tissue_ptwas")


# %%


# %%
results_df = df
results_df$upper_or =(results_df$upper_count / (results_df$total_count - results_df$upper_count))/(5/95)
results_df$lower_or =(results_df$lower_count / (results_df$total_count - results_df$lower_count))/(5/95)

results_df$phyper_upper = phyper(results_df$upper_count, results_df$number_of_genes_topbottom, (results_df$number_of_genes_tissue - results_df$number_of_genes_topbottom), results_df$total_count)
results_df$phyper_lower = phyper(results_df$lower_count, results_df$number_of_genes_topbottom, (results_df$number_of_genes_tissue - results_df$number_of_genes_topbottom), results_df$total_count)
results_df$phyper_upper_bh_adj = p.adjust(results_df$phyper_upper, method = "BH")
results_df$phyper_lower_bh_adj = p.adjust(results_df$phyper_lower, method = "BH")


results_df[results_df$phyper_upper_bh_adj < 0.05, ] %>% dplyr::select(tissue, ptwas_category, phyper_upper_bh_adj, phyper_upper_bh_adj)

results_df[results_df$phyper_lower_bh_adj < 0.05,] %>% dplyr::select(tissue, ptwas_category, phyper_upper_bh_adj, phyper_lower_bh_adj)
# %%
