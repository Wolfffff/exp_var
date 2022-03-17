# Script for making superheat correlation plots -- 

source(here::here("functions.R"))

metric_df = readRDS(here::here("snakemake/Rdatas/gene_metrics.RDS"))
ea_df = read.csv(here::here("snakemake/metadata/EA_metadata.csv"), header=T)
rc3_df = read.csv(here::here("snakemake/metadata/recount3_metadata.csv"), header=T)
metadata_df = bind_rows(ea_df,rc3_df)

library(corrplot)
library(vegan)
library(ape)
library(Hmisc)
library(superheat)

rank_list = list()
metric_cor_list = list()
metric = "sd"
print(metric)
rank_mat = ldply(metric_df[[metric]][,-1], rank)
rank_mat = t(rank_mat[,-1])
rownames(rank_mat) = rownames(metric_df[[metric]][,-1])
colnames(rank_mat) = colnames(metric_df[[metric]][,-1])


mat = as.matrix(metric_df[[metric]][,-1])
ord1 = match(colnames(mat),metadata_df$id)
ord2 = clusters = metadata_df$group[match(colnames(mat),metadata_df$id)]
leveled = factor(ord2,levels = c("GTEx", "TCGA", "Other - Expression Atlas", "Other - recount3"))
ord3 = sort.int(leveled,index.return=T)
ord3_ix = ord3$ix
ord3_x = ord3$`x`



mat = mat[,ord3_ix]
M = rcorr(mat, type = "spearman")$r


png("superheat.png", height = 2160, width = 2160)
superheat(M, membership.cols = ord3_x,#, row.dendrogram=TRUE, col.dendrogram=TRUE,
heat.lim = c(-1, 1),# X.text = round(as.matrix(M4), 2),
X.text.size = 3, grid.hline = FALSE,
# bottom.label.text.angle = 90,
          grid.vline = FALSE, legend.width = 11, legend.text.size = 36,legend.height = 0.5,
          legend.vspace = -0.2) + theme_minimal()
          
dev.off()