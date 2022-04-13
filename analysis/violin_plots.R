# %%
source(here::here("functions.R"))

metric_df = readRDS(here::here("snakemake/Rdatas/gene_metrics.RDS"))
ea_df = read.csv(here::here("snakemake/metadata/EA_metadata.csv"),
                 header=T, comment.char = "#")
rc3_df = read.csv(here::here("snakemake/metadata/recount3_metadata.csv"),
                  header=T, comment.char = "#")
metadata_df = bind_rows(ea_df,rc3_df)
rank_df = read.csv(here::here("data/pca_ranks.csv"), header = TRUE)[, -1]



p = quantile_violin_plot(rank_df$sd,rank_df$mean) + 
    ylab("SD Rank") + 
    xlab("metric quantile") + 
    geom_boxplot(width=0.1) + theme_classic() + theme(legend.position = "none")
save_plot("test.png", p)


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

# %%
library(psych)
library(corrplot)

rank_df = read.csv(here::here("data/pca_ranks.csv"), header = TRUE)[, -1]
rank_df = dplyr::rename(rank_df, gene = Gene)
pi_vals = read.csv(here::here("data/annotation/pi_ceu_results.csv"))

merged = merge(rank_df, pi_vals, by.x = "gene", by.y = "gene")
merged_filtered = merged[,which(names(merged) %in% c("gene","mean","sd", "pi"))]
p = quantile_violin_plot( merged_filtered$pi  ,merged_filtered$sd,ntiles=10) + 
    ylab("SD Rank") + xlab("pi value quantile") + geom_boxplot(width=0.1) + 
    stat_summary(fun = "mean", geom = "point", color = "red")
save_plot("test.png", p)
save_plot(here::here("data/plots/violin_plots/pi.jpg"), p, base_height = 6, base_asp = 1.) 

library(ppcor)
merged = merged[complete.cases(merged),]
pcor(merged[,c("sd", "pi", "mean")], method = "spearman") 
pcor(rank_df[,c("sd", "mean_connectivity", "mean")], method = "spearman") 


# %%

# %%
p = quantile_violin_plot(rank_df$median_connectivity, rank_df$sd, ntiles=10) + ylab("SD Rank") + xlab("Median connectivity quantile") + geom_boxplot(width=0.1) + stat_summary(fun = "mean", geom = "point", color = "red")
save_plot(here::here("data/plots/violin_plots/connectivity.jpg"), p, base_height = 6, base_asp = 1.) 


# %%