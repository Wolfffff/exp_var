source("functions.R")
mean_var = readRDS("Rdatas/mean_variance.RDS")

rank_mat = ldply(mean_var$var$merged[,-1], rank)
rank_mat = t(rank_mat[,-1])
colnames(rank_mat) = colnames(mean_var$var$merged)

pak::pkg_install("corrplot")
library(corrplot)
var_cor = cor(mean_var$var$merged[,-1], method = "s")
png("data/var_corr_plot.png", height = 2080, width = 2080)
corrplot.mixed(var_cor, upper = "ellipse")
dev.off()

all(cor(rank_mat) == var_cor)

PC_scores = as.matrix(rank_mat) %*% eigen(var_cor)$vectors
rank(PC_scores[,1])
gene_var_rank = rank(PC_scores[,1])
names(gene_var_rank) = (mean_var$var$merged[,1])
saveRDS(gene_var_rank, file = "./Rdatas/gene_var_rank.RDS")
data.frame(gene = names(gene_var_rank), rank = gene_var_rank, row.names = NULL) %>%
    write_csv("Rdatas/gene_var_rank.csv")
which(gene_var_rank == 2)
