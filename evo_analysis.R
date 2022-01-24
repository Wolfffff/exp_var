metric_df = readRDS(here::here("snakemake/Rdatas/mean_variance.RDS"))
rank_df = read.csv(here::here("data/pca_ranks.csv"), header = TRUE)[, -1]
imkt_df = read.csv(here::here("data/annotation/imkt_results.csv"))
imkt_df$gene = str_split_fixed(imkt_df$gene,'\\.',Inf)[,1]
imkt_df = imkt_df[,-which(names(imkt_df) %in% c("X"))]

merged = merge(rank_df, imkt_df, by.x = "gene", by.y = "gene")

p_format <- function(x, ndp=5)
{
  out <- format(round(as.numeric(x),ndp),ns=ndp,scientific=T,just="none")
}



library(psych)
merged_filtered = merged[,which(names(merged) %in% c("gene","means", "var", "sd", "cv", "alpha.symbol","Divergence.metrics.omega"))]
cor_mat = cor(merged_filtered[,-1],method="spearman")

# Worth looking into -- FDR corrected here.
cor_test_mat = corr.test(merged_filtered[,-1],method="spearman",adjust="fdr")
png(here::here("data/plots/SpearmanCorrelations/corr_plot_with_pvals.png"), height = 2160, width = 2160)
corrplot(cor_mat,method='ellipse')
mtext("Spearman Correlation Plot - Numbers are p-vals", at=3.5, line=-0.5, cex=4)
pos <- expand.grid(1:ncol(cor_test_mat$p), ncol(cor_test_mat$p):1)
text(pos, p_format(cor_test_mat$p))
dev.off()