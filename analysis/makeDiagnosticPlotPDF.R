library(tidyverse)
library(png)
library(cowplot)
library(patchwork)
library(ggthemes)
library(grid)

ea_df = read.csv(here::here("snakemake/metadata/EA_metadata.csv"),
                 header=T, comment.char = "#")
rc3_df = read.csv(here::here("snakemake/metadata/recount3_metadata.csv"),
                  header=T, comment.char = "#")
metadata_df = bind_rows(ea_df,rc3_df)

metadata_df <- arrange(metadata_df, id)

ord2 = metadata_df$group
leveled = factor(ord2,levels = c("GTEx", "TCGA", "Other - Expression Atlas", "Other - recount3"))
ord3 = sort.int(leveled,index.return=T)
ord3_ix = ord3$ix
ids = metadata_df$id[ord3_ix]

pca_plots = lapply(ids[1:2], function(x) 
    readRDS(here::here(paste0("snakemake/Rdatas/plots/all_pca_", x, ".rds"))))

mean_var_plots = lapply(ids[1:2], function(x) 
    readPNG(here::here(paste0("data/plots/meanVar/meanVar_", x, "_residual.png")), native = TRUE))

save_plot("test.png", rasterGrob(mean_var_plots[[1]]))
layout = 
"AB
CD"
p = with(pca_plots[[1]],
    uncorrected + ggtitle("Uncorrected") + 
    batch + ggtitle("Known batch effects controlled") + 
    clean + ggtitle("Batch effects controlled + outliers removed") + 
    rasterGrob(mean_var_plots[[1]]) + 
    plot_layout(design = layout) +
    plot_annotation(title = ids[1], theme = theme_tufte()))

save_plot("test.png", p, base_height = 6, base_asp = 1.2, ncol = 3)
