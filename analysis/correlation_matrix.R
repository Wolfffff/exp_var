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
library(ggrepel)

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
study_label = rep(" ", 60)
study_label[seq(2, 60, 2)] = seq(2, 60, 2)
heatmap = ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile() +
    scale_fill_viridis_c(alpha = 1)  + theme_tufte() +
    labs(y = "    GTEx          TCGA            Misc.", 
         x = "   GTEx                TCGA                  Misc.", fill = "") +
    scale_x_discrete(breaks = c(12, 30, 47), label = c("gTEX", "TCGA", "Misc")) +
    scale_y_discrete(label = study_label) +
    theme(legend.position = "bottom", 
          legend.key.width= unit(4, 'cm'), axis.text.y = element_text(size = 15)) 
b_size_df = data.frame(start = c(0, 24, 35), end = c(24, 35, 57))  + .5
heatmap = heatmap + geom_rect(data = b_size_df, color = "black", alpha = 0, size = 2,
                                aes(x = NULL, y = NULL, fill = NULL, xmin=start, xmax=end,
                                    ymin=start, ymax=end))
}
save_plot("test.png", heatmap, base_height = 6, base_asp = 1.4)

save_plot(here::here("data/plots/SpearmanCorrelations/sd_corr_plot_heatmap.png"), heatmap, base_height = 8, base_asp = 1)
{
histogram = data.frame(Correlations = M[lower.tri(M)]) %>%
    ggplot(aes(Correlations)) + geom_histogram(bins = 100) +
    theme_tufte() + labs(x = "Spearman correlation across studies", y = "Counts") 
save_plot("data/plots/SpearmanCorrelations/sd_corr_plot_histogram.png", plot = histogram, base_height = 4, base_asp = 1.5)

pcoa_df = readRDS(here::here("snakemake/Rdatas/PCoA_df_sd.RDS"))
pcoa_df$source = gsub("Other", "Misc.", pcoa_df$source)
pallet = wes_palette("Royal1", 4)
pallet[3] = wes_palette("Rushmore1")[3]
PCoA = ggplot(pcoa_df, aes(Axis.1, Axis.2, label = study, color = source)) + 
    geom_point() + geom_text_repel(max.overlaps = 10, show.legend = FALSE, size = 6) + 
    coord_fixed() +
        scale_color_manual(values = pallet) + labs(x = "PCoA axis 1", y = "PCoA axis 2", color = "Study\nsource") + theme_tufte() + 
        #scale_x_continuous(limits = c(-.45, .4)) + scale_y_continuous(limits = c(-.5, .4)) +
               theme(legend.position = "bottom", 
                     legend.title = element_blank(), 
                     legend.background = element_blank(), 
                     legend.text = element_text(size = 23),
                     legend.margin =margin(r=1.5,l=1.5,t=0.,b=0.),
                     axis.title = element_text(size = 28),
                     axis.text = element_text(size = 18)) + 
                     ggtitle("C.") + theme(plot.title = element_text(size = 30)) + 
                     guides(color=guide_legend(ncol=2, override.aes = list(size=3))) 
save_plot("test.png", PCoA, base_height = 7, base_asp = 1.4)

density = readRDS(here::here("snakemake/Rdatas/plots/density_plot.RDS"))
density[["scaled"]] = density[["scaled"]] + 
                      theme_tufte() + theme(legend.position = "none") +
                      ggtitle("D.") + 
                      theme(plot.title = element_text(size = 30)) +
                      theme(axis.title = element_text(size = 28),
                            axis.text = element_text(size = 18)) 
density[["unscaled"]] = density[["unscaled"]] + 
                        theme_tufte() + theme(legend.position = "none") +
                        theme(axis.title.y = element_blank(), axis.title = element_text(size = 28)) 
density = density[["scaled"]] +
          inset_element(density[["unscaled"]], 
                        0.33, 0.33, 1, 1)
{layout <- 
"AAACCC
AAACCC
AAACCC
BBBDDD
BBBDDD"}
panel = heatmap + ggtitle("A.") + theme(plot.title = element_text(size = 30),
                                        axis.title = element_text(size = 30),
                                        axis.text = element_text(size = 15),
                                        legend.text = element_text(size = 18)) +
        histogram + ggtitle("B.") + theme(plot.title = element_text(size = 30),
                                          axis.title = element_text(size = 28),
                                          axis.text = element_text(size = 18)) +
        PCoA + 
        density  +
        plot_layout(design = layout)
save_plot(here::here("test.png"), panel, base_height = 7.5, base_asp = 1.2, ncol = 2, nrow = 2)
save_plot(here::here("data/plots/fig1.png"), panel, base_height = 7.5, base_asp = 1.2, ncol = 2, nrow = 2)
}







