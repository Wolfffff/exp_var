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
upper_quantiles = list()
lower_quantiles = list()
for(metric in c("mean","sd")){
    cutoff = quantile(rank_df[[metric]], .95)
    subset = rank_df[rank_df[[metric]] >= cutoff,]$gene
    upper_quantiles[[metric]] = subset
    cutoff = abs(quantile(-rank_df[[metric]], .95))
    lower_quantiles[[metric]] = rank_df[rank_df[[metric]] <= quantile(rank_df[[metric]], .05),]$gene
}
term2gene_df = ptwas_table[, c("Trait","Gene")]
ptwas_table_merged = merge(term2gene_df,ptwas_traits, by.x = "Trait", by.y = "ID", all.x = TRUE)
# %%


# %%
rank_df = dplyr::rename(rank_df, Gene = gene)

library(cowplot)

log_reg_results = list()
for (cat in unique(ptwas_table_merged$Category)){
    cat_df = ptwas_table_merged[ptwas_table_merged$Category == cat,]
    rank_df_with_disease = rank_df %>% mutate(disease = if_else(Gene %in% cat_df$Gene, 1, 0))
    print(table(rank_df_with_disease$disease))
    log_reg_results[[as.character(cat)]] = glm(disease ~ mean + sd , data = rank_df_with_disease, family = "binomial") %>% summary
}

lapply(log_reg_results, print)

library(ggpubr)
library(ggplot2)
library(purrr)
# Boxplot ranks
t_tests = list()
t_tests_plots = list()
for (cat in unique(ptwas_table_merged$Category)){
    if(is.na(cat)){
        next
    }
    cat_df = ptwas_table_merged[ptwas_table_merged$Category == cat,]
    rank_df_with_disease = rank_df %>% mutate(disease = if_else(Gene %in% cat_df$Gene, 1, 0))
    t_tests_plots[[cat]] <- ggplot(rank_df_with_disease, aes(group=disease, y=sd)) + 
    geom_boxplot() + ggtitle(paste0(cat," linked genes from PTWAS"))
    t_tests[[cat]] = t.test(rank_df_with_disease[rank_df_with_disease$disease == 0,]$sd,rank_df_with_disease[rank_df_with_disease$disease == 1,]$sd,)
}

library(cowplot)
png(here::here("data/plots/disease_boxplots.png"), height = 3840, width = 3840)
ggarrange(plotlist=t_tests_plots)
dev.off()

# %%