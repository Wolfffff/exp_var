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
library(ggthemes)
library(patchwork)

rank_list = list()
metric_cor_list = list()
metric = "sd"
print(metric)
rank_mat = ldply(metric_df[[metric]][,-1], rank, na.last = "keep")
rank_mat = t(rank_mat[,-1])
rownames(rank_mat) = metric_df[[metric]][, 1]
colnames(rank_mat) = colnames(metric_df[[metric]][,-1])

metric_matrix = as.matrix(metric_df[[metric]][,-1])
mat = metric_matrix[complete.cases(metric_matrix), ]

ord1 = match(colnames(mat),metadata_df$id)
ord2 = clusters = metadata_df$group[match(colnames(mat),metadata_df$id)]
leveled = factor(ord2,levels = c("GTEx", "TCGA", "Other - Expression Atlas", "Other - recount3"))
ord3 = sort.int(leveled,index.return=T)
ord3_ix = ord3$ix
ord3_x = ord3$`x`

mat = mat[,ord3_ix]
M = rcorr(mat, type = "spearman")$r

melted_cormat <- reshape2::melt(M)
#plot heatmap
{
heatmap = ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile() +
    scale_fill_viridis_c(alpha = 1)  + theme_tufte() +
    labs(y = "", x = "", fill = "") +
    scale_x_discrete(breaks = c(12, 30, 47), label = c("gTEX", "TCGA", "Other")) +
    theme(legend.position = "bottom", legend.key.width= unit(4, 'cm')) 
b_size_df = data.frame(start = c(0, 24, 35), end = c(24, 35, 60))  + .5
heatmap = heatmap + geom_rect(data = b_size_df, color = "#36454F", alpha = 0, size = 1,
                                aes(x = NULL, y = NULL, fill = NULL, xmin=start, xmax=end,
                                    ymin=start, ymax=end))
}
save_plot("test.png", heatmap, base_height = 6, base_asp = 1.4)

save_plot(here::here("data/plots/SpearmanCorrelations/sd_corr_plot_heatmap.png"), heatmap, base_height = 8, base_asp = 1)

histogram = data.frame(Correlations = M[lower.tri(M)]) %>%
    ggplot(aes(Correlations)) + geom_histogram(bins = 100) +
    theme_tufte() + labs(x = "Spearman correlation across studies", y = "Counts") 
save_plot("data/plots/SpearmanCorrelations/sd_corr_plot_histogram.png", plot = histogram, base_height = 4, base_asp = 1.5)

PCoA = readRDS(here::here("snakemake/Rdatas/plots/PCoA_plot_sd.RDS"))
density = readRDS(here::here("snakemake/Rdatas/plots/density_plot.RDS"))
density[["scaled"]] = density[["scaled"]] + 
                      theme_tufte() + theme(legend.position = "none") +
                      ggtitle("D.") + 
                      theme(plot.title = element_text(size = 30)) +
                      theme(axis.title = element_text(size = 18),
                            axis.text = element_text(size = 12)) 
density[["unscaled"]] = density[["unscaled"]] + 
                        theme_tufte() + theme(legend.position = "none")
density = density[["scaled"]] +
          inset_element(density[["unscaled"]], 
                        0.33, 0.33, 1, 1)
layout <- 
"AAACC
AAACC
AAADD
BBBDD"
panel = heatmap + ggtitle("A.") + theme(plot.title = element_text(size = 30),
                                        axis.title = element_text(size = 18),
                                        axis.text = element_text(size = 12)) +
        histogram + ggtitle("B.") + theme(plot.title = element_text(size = 30),
                                          axis.title = element_text(size = 18),
                                          axis.text = element_text(size = 12)) +
        PCoA + theme_tufte() + theme(legend.position = c(0.25, 0.14), 
                     legend.title = element_blank(), 
                     legend.background = element_blank(), 
                     legend.box.background = element_rect(colour = "black"),
                     legend.text = element_text(size = 12),
                     legend.margin =margin(r=1.5,l=1.5,t=0.,b=0.),
                     axis.title = element_text(size = 18),
                     axis.text = element_text(size = 12)) + 
                     ggtitle("C.") + theme(plot.title = element_text(size = 30)) +
        density  +
        plot_layout(design = layout)
save_plot(here::here("data/plots/fig1_panel.png"), panel, base_height = 7, base_asp = 1.4, ncol = 2, nrow = 2)

