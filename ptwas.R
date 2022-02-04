ptwas_table <- read.table("data/annotation/ptwas_table3_sig_genetrait_interactions.txt", header=T, sep=" ") 
ptwas_traits <- read.csv("data/annotation/ptwas_traits.csv")
ptwas_metadata<- read.csv("data/annotation/ptwas_metadata.csv")

library(stringr)
ptwas_table$Gene =  str_split_fixed(ptwas_table$Gene,'\\.',Inf)[,1]
uniq_disease_linked_genes = unique(ptwas_table$Gene)


organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

library(clusterProfiler)
library(plyr)
library(dplyr)

rank_df = read.csv(here::here("data/pca_ranks.csv"), header = TRUE)[, -1]
upper_quantiles = list()
lower_quantiles = list()
for(metric in c("means","sd")){
    cutoff = quantile(rank_df[[metric]], .95)
    subset = rank_df[rank_df[[metric]] >= cutoff,]$gene
    upper_quantiles[[metric]] = subset
    cutoff = abs(quantile(-rank_df[[metric]], .95))
    lower_quantiles[[metric]] = rank_df[rank_df[[metric]] <= quantile(rank_df[[metric]], .05),]$gene
}
term2gene_df = ptwas_table[, c("Trait","Gene")]
ptwas_table_merged = merge(term2gene_df,ptwas_traits, by.x = "Trait", by.y = "ID", all.x = TRUE)
# term2gene_df = data.frame(disease="1",
                        #   Gene=uniq_disease_linked_genes)


# newtable <- merge(table1,table2, by  = "pid") 
upper_go_enrichment = enrichGO(gene  = upper_quantiles[["sd"]],
                         universe      = rank_df$gene,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05,
                         readable      = TRUE)
gse = upper_go_enrichment
library(enrichplot)
png(here::here("plot.png"), height = 2160, width = 2160)
barplot(gse, showCategory=20) 
dev.off()

lower_go_enrichment = enrichGO(gene  = lower_quantiles[["sd"]],
                         universe      = rank_df$gene,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05,
                         readable      = TRUE)
gse = lower_go_enrichment
library(enrichplot)
png(here::here("plot.png"), height = 2160, width = 2160)
barplot(gse, showCategory=20) 
dev.off()


# %%
rank_df = dplyr::rename(rank_df, Gene = gene)

head(rank_df_with_disease)
library(cowplot)
plot = ggplot(rank_df_with_disease, aes(means, sd, groups = disease, color = disease)) + geom_point()
save_plot("test.png", plot)

log_reg_results = list()
for (cat in unique(ptwas_table_merged$Category)){
    cat_df = ptwas_table_merged[ptwas_table_merged$Category == cat,]
    rank_df_with_disease = rank_df %>% mutate(disease = if_else(Gene %in% cat_df$Gene, 1, 0))
    print(table(rank_df_with_disease$disease))
    log_reg_results[[as.character(cat)]] = glm(disease ~ means + sd , data = rank_df_with_disease, family = "binomial") %>% summary
}

lapply(log_reg_results, print)

# Boxplot ranks
t_tests = list()
for (cat in unique(ptwas_table_merged$Category)){
    if(is.na(cat)){
        next
    }
    cat_df = ptwas_table_merged[ptwas_table_merged$Category == cat,]
    rank_df_with_disease = rank_df %>% mutate(disease = if_else(Gene %in% cat_df$Gene, 1, 0))
    ggplot(rank_df_with_disease, aes(group=disease, y=sd)) + 
    geom_boxplot()
    ggsave(paste0(cat,".png"))
    t_tests[[cat]] = t.test(rank_df_with_disease[rank_df_with_disease$disease == 0,]$sd,rank_df_with_disease[rank_df_with_disease$disease == 1,]$sd,)

}
