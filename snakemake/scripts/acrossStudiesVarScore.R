my_logfile = snakemake@log[["log"]]
snakemake@source("logger.R")
log4r_info("Starting.")
print = log4r_info
log4r_info("Loading packages") 
source(here::here("functions.R"))

metric_df = readRDS(here::here("snakemake/Rdatas/gene_metrics.RDS"))
connectivity_df = readRDS(here::here("snakemake/Rdatas/gene_connectivity.RDS"))


#pak::pkg_install(c("corrplot", "vegan", "ape", "Hmisc", "ggrepel", "wesanderson", "missMDA"))
library(missMDA)
library(corrplot)
library(vegan)
library(ape)
library(Hmisc)
library(ggrepel)
library(wesanderson)

ea_df = read.csv(here::here("snakemake/metadata/EA_metadata.csv"), header=T, comment.char = "#")
rc3_df = read.csv(here::here("snakemake/metadata/recount3_metadata.csv"), header=T, comment.char = "#")
metadata_df = bind_rows(ea_df,rc3_df)

rank_list = list()
metric_cor_list = list()
for (metric in c("mean", "sd")){
    print(metric)
    rank_mat = ldply(metric_df[[metric]][,-1], rank, na.last = "keep")
    rank_mat = t(rank_mat[,-1])
    rownames(rank_mat) = metric_df[[metric]][, 1]
    colnames(rank_mat) = colnames(metric_df[[metric]][,-1])

    metric_matrix = as.matrix(metric_df[[metric]][,-1])
    metric_matrix = metric_matrix[complete.cases(metric_matrix), ]
    metric_cor = rcorr(metric_matrix, type = "spearman")$r

    res <- pcoa(abs(1 - metric_cor))

    pcoa_df = data.frame(res$vectors[, 1:2],
               study = rownames(data.frame(res$vectors[, 1:2])))    
    pcoa_df$source = metadata_df[match(pcoa_df$study, metadata_df$id), "group"]

    pallet = wes_palette("Royal1", 4)
    pallet[3] = wes_palette("Rushmore1")[3]
    pcoa_plot = ggplot(pcoa_df, aes(Axis.1, Axis.2, label = study, color = source)) + 
        geom_point() + geom_text_repel(max.overlaps = 15, show.legend = FALSE) + coord_fixed() +
         scale_color_manual(values = pallet) + labs(x = "PCoA Axis 1", y = "PCoA Axis 2", color = "Study\nsource") + theme_cowplot() +
        theme(legend.position = "top")
    saveRDS(pcoa_plot, here::here(paste0("snakemake/Rdatas/plots/PCoA_plot_", metric, ".RDS")))
    save_plot(here::here(paste0("data/plots/SpearmanCorrelations/",metric,"_PCoA_plot.png")), pcoa_plot, base_height = 7, base_asp = 2)

    eig = eigen(metric_cor)
    if(all(eig$vectors[,1] < 0)){
        eig$vectors[,1] <- -eig$vectors[,1]
    }

    nb <- estim_ncpPCA(rank_mat,method.cv = "gcv",  ncp.min = 0, ncp.max = 10, verbose = FALSE)
    imputed_rank = imputePCA(rank_mat, ncp = nb$ncp, verbose = FALSE)

    PC_scores = as.matrix(imputed_rank$completeObs) %*% eig$vectors
    vars = apply(PC_scores, 2, function(x) var(x))
    print(paste(round(vars/sum(vars) * 100, 1)[1:5], collapse = "% "))
    gene_rank = rank(PC_scores[,1], na.last = "keep")
    names(gene_rank) = metric_df[[metric]][,1]
    rank_list[[metric]] = gene_rank
    metric_cor_list[[metric]] = metric_cor
}
saveRDS(PC_scores, file=snakemake@output[[3]])

rank_df = data.frame(bind_cols(rank_list))
rank_df$Gene = metric_df[[1]][,1]
names(connectivity_df) = c("Gene", "mean_connectivity", "median_connectivity")
rank_df = inner_join(rank_df, connectivity_df, by = "Gene") %>% dplyr::relocate(Gene) %>% as_tibble

write.csv(rank_df, file=here::here("data/pca_ranks.csv"))
saveRDS(rank_df, file=snakemake@output[[1]])
saveRDS(metric_cor_list, file=snakemake@output[[2]])
