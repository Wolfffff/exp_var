# Script for making pretty correlation matrices

source(here::here("functions.R"))

metric_df = readRDS(here::here("snakemake/Rdatas/gene_metrics.RDS"))
ea_df = read.csv(here::here("snakemake/metadata/EA_metadata.csv"), header=T)
rc3_df = read.csv(here::here("snakemake/metadata/recount3_metadata.csv"), header=T)
metadata_df = bind_rows(ea_df,rc3_df)

library(corrplot)
library(vegan)
library(ape)
library(Hmisc)

rank_list = list()
metric_cor_list = list()
metric = "sd"
print(metric)
rank_mat = ldply(metric_df[[metric]][,-1], rank)
rank_mat = t(rank_mat[,-1])
rownames(rank_mat) = rownames(metric_df[[metric]][,-1])
colnames(rank_mat) = colnames(metric_df[[metric]][,-1])

M = rcorr(as.matrix(metric_df[[metric]][,-1]), type = "spearman")$r
ord = corrMatOrder(M, order = 'AOE')
M2 = M[ord, ord]
png(here::here(paste0(metric,"_corr_plot_unordered.png")), height = 6080, width = 6080)
corrplot.mixed(M2, upper = "ellipse")
dev.off()

M2
library(superheat)
M.size <- scale(M) + 2


# png("superheat.png", height = 2160, width = 2160)
# superheat(M,#, row.dendrogram=TRUE, col.dendrogram=TRUE,
# heat.lim = c(-1, 1), X.text = round(as.matrix(M), 1),
# X.text.size = M.size, grid.hline = FALSE,
#           grid.vline = FALSE)
# dev.off()


metadata_df = metadata_df[match(rownames(M2),metadata_df$id),]
M3 = M2[match(metadata_df$id,rownames(M2)),match(metadata_df$id,rownames(M2))]

corr_mat_df <- as.data.frame(M3)
clusters = metadata_df$group[match(rownames(M3),metadata_df$id)]
M4 = M3[order(clusters),order(clusters)]
clusters = metadata_df$group[match(rownames(M4),metadata_df$id)]

png("superheat.png", height = 2160, width = 2160)
superheat(M4, membership.cols = clusters, order.cols = order(clusters), order.rows = order(clusters),#, row.dendrogram=TRUE, col.dendrogram=TRUE,
heat.lim = c(-1, 1), X.text = round(as.matrix(M4), 2),
X.text.size = 3, grid.hline = FALSE,
# bottom.label.text.angle = 90,
          grid.vline = FALSE, legend.width = 11, legend.text.size = 36,legend.height = 0.5,
          legend.vspace = -.2) + theme_minimal()
          
dev.off()



# %%
gene_metrics = readRDS(here::here("snakemake/Rdatas/gene_metrics.RDS"))
