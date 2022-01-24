ptwas_table <- read.table("data/annotation/ptwas_table3_sig_genetrait_interactions.txt", header=T, sep=" ") 
ptwas_traits <- read.csv("data/annotation/ptwas_traits.csv")
ptwas_metadata<- read.csv("data/annotation/ptwas_metadata.csv")

library(stringr)
ptwas_table$Gene =  str_split_fixed(ptwas_table$Gene,'\\.',Inf)[,1]
uniq_disease_linked_genes = unique(ptwas_table$Gene)


# Analysis
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

library(clusterProfiler)
library(plyr)
library(dplyr)

rank_df = read.csv(here::here("data/pca_ranks.csv"), header = TRUE)[, -1]
for(metric in c("means","var","sd","cv")){
    cutoff = quantile(rank_df[[metric]], .95)
    subset = rank_df[rank_df[[metric]] >= cutoff]
    gene_ids = subset$gene
}
term2gene_df = ptwas_table[, c("Trait","Gene")]
term2gene_df = data.frame(disease="1",
                          Gene=uniq_disease_linked_genes)

ego <- enricher(gene = gene_ids,
                universe = unique(c(uniq_disease_linked_genes,rank_df$gene)), 
                pvalueCutoff = 10, pAdjustMethod="none",
                TERM2GENE = term2gene_df)

gse <- enricher()

require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

emapplot(gse, showCategory = 10)

cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)

ridgeplot(gse) + labs(x = "enrichment distribution")

rank_df = dplyr::rename(rank_df, Gene = gene)
rank_df_with_disease = left_join(rank_df, term2gene_df, by = "Gene") %>%
    mutate(disease = ifelse(is.na(disease), 0, 1)) %>%
    relocate(Gene)

head(rank_df_with_disease)

library(cowplot)
plot = ggplot(rank_df_with_disease, aes(means, sd, groups = disease, color = disease)) + geom_point()
save_plot("test.png", plot)
glm(disease ~ means + sd, data = rank_df_with_disease, family = "binomial") %>% summary           
