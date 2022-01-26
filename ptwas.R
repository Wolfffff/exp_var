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
top_quantiles = list()
for(metric in c("means","var","sd","cv")){
    cutoff = quantile(rank_df[[metric]], .95)
    subset = rank_df[rank_df[[metric]] >= cutoff,]$gene
    top_quantiles[[metric]] = subset
}
term2gene_df = ptwas_table[, c("Trait","Gene")]
ptwas_table_merged = merge(term2gene_df,ptwas_traits, by.x = "Trait", by.y = "ID", all.x = TRUE)
# term2gene_df = data.frame(disease="1",
                        #   Gene=uniq_disease_linked_genes)


newtable <- merge(table1,table2, by  = "pid") 
GO_enrichment = enrichGO(gene  = top_quantiles[["var"]],
                         universe      = rank_df$gene,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05,
                         readable      = TRUE)
gse = GO_enrichment

library(enrichplot)
png(here::here("plot.png"), height = 2160, width = 2160)
barplot(gse, showCategory=20) 
dev.off()

require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

emapplot(gse, showCategory = 10)

cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)

ridgeplot(gse) + labs(x = "enrichment distribution")


# %%
# rank_df = dplyr::rename(rank_df, Gene = gene)
rank_df_with_disease = left_join(rank_df, ptwas_table_merged, by = "Gene") %>%
    # mutate(disease = ifelse(is.na(disease), 0, 1)) %>%
    relocate(Gene)

head(rank_df_with_disease)

library(cowplot)
plot = ggplot(rank_df_with_disease, aes(means, sd, groups = disease, color = disease)) + geom_point()
save_plot("test.png", plot)

log_reg_results = list()
for (cat in unique(rank_df_with_disease$Category)){
    cat_df = rank_df_with_disease[rank_df_with_disease$disease == cat,]
    log_reg_results[[cat]] = glm(means + sd, data = rank_df_with_disease, family = "binomial") %>% summary 
}
glm(dummy() ~ m          
