my_logfile = snakemake@log[["log"]]
snakemake@source("logger.R")
log4r_info("Starting.")
print = log4r_info
log4r_info("Loading packages") 
source(here::here("functions.R"))

metric_df = readRDS(here::here("snakemake/Rdatas/gene_metrics.RDS"))

pak::pkg_install(c("corrplot", "vegan", "ape", "Hmisc"))
library(corrplot)
library(vegan)
library(ape)
library(Hmisc)

rank_list = list()
metric_cor_list = list()
for (metric in c("means", "var", "sd", "cv")){
    print(metric)
    rank_mat = ldply(metric_df[[metric]][,-1], rank)
    rank_mat = t(rank_mat[,-1])
    rownames(rank_mat) = rownames(metric_df[[metric]][,-1])
    colnames(rank_mat) = colnames(metric_df[[metric]][,-1])

    metric_cor = rcorr(as.matrix(metric_df[[metric]][,-1]), type = "spearman")$r
    png(here::here(paste0("data/plots/SpearmanCorrelations/",metric,"_corr_plot.png")), height = 6080, width = 6080)
    corrplot.mixed(metric_cor, upper = "ellipse")
    dev.off()

    res <- pcoa(abs(1 - metric_cor))
    png(here::here(paste0("data/plots/SpearmanCorrelations/",metric,"_PCoA_plot.png")), height = 1080, width = 1080)
    biplot(res)
    dev.off()

    eig = eigen(metric_cor)
    PC_scores = as.matrix(rank_mat) %*% eig$vectors
    print(paste(round(eig$values/sum(eig$values) * 100, 1)[1:5], collapse = "% "))
    gene_rank = rank(PC_scores[,1])
    names(gene_rank) = metric_df[[metric]][,1]
    rank_list[[metric]] = gene_rank
    metric_cor_list[[metric]] = metric_cor
}

rank_df = data.frame(bind_cols(rank_list))
rank_df$gene = metric_df[[1]][,1]
write.csv(rank_df, file=here::here("data/pca_ranks.csv"))
saveRDS(rank_df, file=snakemake@output[[1]])
saveRDS(metric_cor_list, file=snakemake@output[[2]])
