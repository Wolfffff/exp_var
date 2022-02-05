source(here::here("functions.R"))


metric_df = readRDS(here::here("snakemake/Rdatas/gene_metrics.RDS"))
rank_list = list()
pak::pkg_install("corrplot")
library(corrplot)
library(vegan)
library(ape)

for (metric in c("mean", "sd")){
    print(metric)
    rank_mat = ldply(metric_df[[metric]][,-1], rank)
    rank_mat = t(rank_mat[,-1])
    rownames(rank_mat) = rownames(metric_df[[metric]][,-1])
    colnames(rank_mat) = colnames(metric_df[[metric]][,-1])

    metric_cor = cor(metric_df[[metric]][,-1], method = "s")
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
}

rank_df = data.frame(bind_cols(rank_list))
rank_df$gene = metric_df[[1]][,1]
write.csv(rank_df, file=here::here("data/pca_ranks.csv"))




#  Modularity
pak::pkg_install("diogro/evolqg")
devtools::install_github("diogro/yamda-r", subdir = "package")
library(evolqg)
library(yamdar)

groups = c(rep("EA", 8), rep("gTEX", 17),  rep("tcga", 12), rep("misc", 4))
names(groups) =  colnames(metric_cor)
toHypotMatrix(groups)

TestModularity(var_cor, toHypotMatrix(groups), permutations = 10000)
