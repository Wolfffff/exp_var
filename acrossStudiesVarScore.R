source(here::here("functions.R"))


metric_df = readRDS(here::here("snakemake/Rdatas/gene_metrics.RDS"))
rank_list = list()
pak::pkg_install("corrplot")
library(corrplot)

for (metric in c("means", "var", "sd", "cv")){
    print(metric)
    rank_mat = ldply(metric_df[[metric]][,-1], rank)
    rank_mat = t(rank_mat[,-1])
    rownames(rank_mat) = rownames(metric_df[[metric]][,-1])
    colnames(rank_mat) = colnames(metric_df[[metric]][,-1])

    metric_cor = cor(metric_df[[metric]][,-1], method = "s")
    png(here::here(paste0("data/plots/SpearmanCorrelations/",metric,"_corr_plot.png")), height = 6080, width = 6080)
    corrplot.mixed(metric_cor, upper = "ellipse")
    dev.off()

    PC_scores = as.matrix(rank_mat) %*% eigen(metric_cor)$vectors
    
    gene_rank = rank(PC_scores[,1])
    names(gene_rank) = metric_df[[metric]][,1]
    rank_list[[metric]] = gene_rank
}

rank_df = data.frame(bind_cols(rank_list))
rank_df$gene = metric_df[[1]][,1]
write.csv(rank_df, file=here::here("data/pca_ranks.csv"))

res <- pcoa(abs(1 - metric_cor))
res$values
png("../test.png")
biplot(res)
dev.off()

#  Modularity
pak::pkg_install("diogro/evolqg")
devtools::install_github("diogro/yamda-r", subdir = "package")
library(evolqg)
library(yamdar)

groups = c(rep("EA", 8), rep("gTEX", 17),  rep("tcga", 12), rep("misc", 4))
names(groups) =  colnames(var_cor)
toHypotMatrix(groups)

TestModularity(var_cor, toHypotMatrix(groups), permutations = 10000)
