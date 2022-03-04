# %%
source(here::here("functions.R"))
library(here)
library(sjmisc)
library(pryr)
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

pw_upper <- pairwise_termsim(local_go_upper) 
p1 <- emapplot(pw_upper, showCategory = 20) + theme_cowplot() + 
theme(legend.position = "none", panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.title.x = element_blank(), 
      axis.title.y = element_blank(), 
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()) + ggtitle("High variation")
p1$data$color = mean(p1$data$color)
local_go_lower = enrichGO(gene  = lower_quantiles[[metric]],
                          universe = rank_df$Gene,
                          OrgDb         = org.Hs.eg.db,
                          keyType       = 'ENSEMBL',
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05,
                          readable      = TRUE)
pw_lower <- pairwise_termsim(local_go_lower) 
p2 <- emapplot(pw_lower, showCategory = 20) + theme_cowplot() + 
theme(legend.position = "none", panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.title.x = element_blank(), 
      axis.title.y = element_blank(), 
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()) + ggtitle("Low variation")
p2$data$color = mean(p1$data$color)
p2$data$name[2] = "RNA splicing, via transesterification\n reactions with bulged adenosine as nucleophile"

p12 = p1 + p2

save_plot(here::here("data/plots/local_go_lowerUpper.png"),p12, base_height = 8, base_asp = 1.6, ncol = 2)
# %%





