source(here::here("functions.R"))

metric_df = readRDS(here::here("snakemake/Rdatas/gene_metrics.RDS"))
ea_df = read.csv(here::here("snakemake/metadata/EA_metadata.csv"), header=T)
rc3_df = read.csv(here::here("snakemake/metadata/recount3_metadata.csv"), header=T)
metadata_df = bind_rows(ea_df,rc3_df)
rank_df = read.csv(here::here("data/pca_ranks.csv"), header = TRUE)[, -1]



quantile_violin_plot(rank_df$sd,rank_df$mean) + ylab("SD Rank") + xlab("metric quantile") + geom_boxplot(width=0.1)# + theme_minimal() + theme(legend.position = "none")
ggsave("example.jpg")


# %%
# %%
ptwas_table <- read.table(here::here("data/annotation/ptwas_table3_sig_genetrait_interactions.txt"), header=T, sep=" ")
ptwas_traits <- read.csv(here::here("data/annotation/ptwas_traits.csv"))
ptwas_metadata <- read.csv(here::here("data/annotation/ptwas_metadata.csv"))

library(stringr)
ptwas_table$Gene =  str_split_fixed(ptwas_table$Gene,'\\.',Inf)[,1]
uniq_disease_linked_genes = unique(ptwas_table$Gene)

organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
# %%

# %%
library(clusterProfiler)
library(plyr)
library(dplyr)

rank_df = read.csv(here::here("data/pca_ranks.csv"), header = TRUE)[, -1]
term2gene_df = ptwas_table[, c("Trait","Gene")]
ptwas_table_merged = merge(term2gene_df,ptwas_traits, by.x = "Trait", by.y = "ID", all.x = TRUE)
# %%

rank_df = dplyr::rename(rank_df, Gene = gene)
cat_df = ptwas_table_merged
rank_df_with_disease = rank_df %>% mutate(disease = if_else(Gene %in% cat_df$Gene, 1, 0))

lapply(log_reg_results, print)

library(ggpubr)
library(ggplot2)
library(purrr)
# Boxplot ranks
for (cat in unique(ptwas_table_merged$Category)){
    if(is.na(cat)){
        next
    }
    cat_df = ptwas_table_merged#[ptwas_table_merged$Category == cat,]
    rank_df_with_disease = rank_df %>% mutate(disease = if_else(Gene %in% cat_df$Gene, 1, 0))
    quantile_violin_plot( rank_df_with_disease$disease,rank_df_with_disease$sd,ntiles=100) + ylab("SD Rank") + xlab("metric quantile") + geom_boxplot(width=0.1) + stat_summary(fun = "mean", geom = "point", color = "red")
    ggsave(here::here(paste0("data/plots/violin_plots/",cat, ".jpg")), width = 18, height = 6, units = "in", dpi = 300)
}
rank_df_with_disease = rank_df %>% mutate(disease = if_else(Gene %in% ptwas_table_merged$Gene, 1, 0))
quantile_violin_plot( rank_df_with_disease$disease,rank_df_with_disease$sd,ntiles=100) + ylab("SD Rank") + xlab("metric quantile") + geom_boxplot(width=0.1) + stat_summary(fun = "mean", geom = "point", color = "red")
ggsave(here::here("data/plots/violin_plots/ptwas_table_merged.jpg"), width = 18, height = 6, units = "in", dpi = 300)

# %%


# %%
library(psych)
library(corrplot)

rank_df = read.csv(here::here("data/pca_ranks.csv"), header = TRUE)[, -1]
rank_df = dplyr::rename(rank_df, gene = Gene)
pi_vals = read.csv(here::here("data/annotation/pi_ceu_results.csv"))

merged = merge(rank_df, pi_vals, by.x = "gene", by.y = "gene")
merged_filtered = merged[,which(names(merged) %in% c("gene","mean","sd", "pi"))]
quantile_violin_plot( merged_filtered$pi  ,merged_filtered$sd,ntiles=100) + ylab("SD Rank") + xlab("metric quantile") + geom_boxplot(width=0.1) + stat_summary(fun = "mean", geom = "point", color = "red")
ggsave(here::here("data/plots/violin_plots/pi.jpg"), width = 18, height = 6, units = "in", dpi = 300)

# %%