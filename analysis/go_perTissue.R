# %%
source(here::here("functions.R"))
library(ggthemes)
library(ggrepel)
library(here)
library(sjmisc)
library(pryr)
library(GOxploreR)
library(moments)
library(data.table)
library(enrichplot)
library("AnnotationDbi")
library("org.Hs.eg.db")
# %%

# %%
# GO analysis on ea
organism = "org.Hs.eg.db"
metrics = c("mean", "sd")

tissues = tools::file_path_sans_ext(dir(here::here("data/pca_ranks_tissue")))
tissue_rank_files = dir(here::here("data/pca_ranks_tissue"), full.names = T)
rank_df_list = lapply(tissue_rank_files, \(x) read.csv(x, header = TRUE))
names(rank_df_list) = tissues


runGOtissue = function(current_tissue = NULL){
    rank_df = rank_df_list[[current_tissue]]
    tail_size = 0.05

    upper_quantiles = list()
    lower_quantiles = list()
    for(metric in c("mean","sd")){
        cutoff = quantile(rank_df[[metric]],1 - tail_size)
        subset = rank_df[rank_df[[metric]] >= cutoff,]$Gene
        upper_quantiles[[metric]] = subset
        cutoff = abs(quantile(rank_df[[metric]], tail_size))
        lower_quantiles[[metric]] = rank_df[rank_df[[metric]] <= cutoff,]$Gene
    }
    # %%

    # %%
    # GO analysis against our list of genes
    metric = "sd"

    outputdir = here::here("data/annotation/per_tissue_go", current_tissue)
    dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)

    local_go_upper = enrichGO(gene  = upper_quantiles[[metric]],
                            universe = rank_df$Gene,
                            OrgDb         = org.Hs.eg.db,
                            keyType       = 'ENSEMBL',
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.01,
                            qvalueCutoff  = 0.01,
                            readable      = TRUE)
    local_go_upper_table = local_go_upper@result |> 
    filter(p.adjust < 0.01, Count > 4) |> 
    mutate(varRank = "high")

    write_csv(local_go_upper_table, file.path(outputdir, "go_upper_quantile.csv"))
    pw_upper <- pairwise_termsim(local_go_upper) 
    pw_upper <- simplify(pw_upper, cutoff=0.7, by="p.adjust", select_fun=min)
    plot_upper <- emapplot(pw_upper, showCategory = 10, cex_label_category = 1.2) + 
        ggtitle("A. High variation") +
            theme_tufte() + 
            theme(legend.position = "none") + 
            theme(plot.title = element_text(size=28),
                    axis.title = element_blank(),
                    axis.text = element_blank(),
                    axis.ticks = element_blank()) 
    plot_upper$data$color = mean(plot_upper$data$color)
    local_go_lower = enrichGO(gene  = lower_quantiles[[metric]],
                            universe = rank_df$Gene,
                            OrgDb         = org.Hs.eg.db,
                            keyType       = 'ENSEMBL',
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.01,
                            qvalueCutoff  = 0.01,
                            readable      = TRUE)
    local_go_lower_table = local_go_lower@result |> 
    filter(p.adjust < 0.01, Count > 4) |> 
    mutate(varRank = "low")
    write_csv(local_go_lower_table, file.path(outputdir, "go_lower_quantile.csv"))
    write_csv(rbind(local_go_upper_table, local_go_lower_table), 
                    file.path(outputdir, "go_upperLower_quantile.csv"))
    pw_lower <- pairwise_termsim(local_go_lower) 
    pw_lower <- simplify(pw_lower, cutoff=0.7, by="p.adjust", select_fun=min)
    plot_lower <- emapplot(pw_lower, showCategory = 10, cex_label_category = 1.2) + 
        ggtitle("B. Low variation") +            
            theme_tufte() + 
            theme(legend.position = "none") + 
            theme(plot.title = element_text(size=28),
                    axis.title = element_blank(),
                    axis.text = element_blank(),
                    axis.ticks = element_blank()) 
    plot_lower$data$color = mean(plot_lower$data$color)
    #p2$data$name[3] = "RNA splicing, via transesterification\n reactions with bulged adenosine as nucleophile"

    p12 = plot_upper + plot_lower + plot_layout(ncol=2)

    plot_dir = here::here("data/plots/per_tissue_go")
    dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
    save_plot(file.path(plot_dir, paste0(current_tissue, ".png")), p12, base_height = 6, base_asp = 2)
    return(0)
}

library(doMC)
registerDoMC(length(tissues))

iteration <- function(idx){
  tryCatch(
    runGOtissue(idx)
    ,error = function(e) print(paste('error',idx))
    )
}

llply(tissues, iteration, .parallel = TRUE)
# %%
