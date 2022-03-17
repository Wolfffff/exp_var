# File for evo analysis including mapping PopHuman
# %%
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
library(corrplot)
merged_filtered = merged[,which(names(merged) %in% c("gene","mean", "sd", "alpha.symbol","Divergence.metrics.omega"))]
cor_mat = cor(merged_filtered[,-1],method="spearman")

# Worth looking into -- FDR corrected here.
cor_test_mat = corr.test(merged_filtered[,-1],method="spearman",adjust="fdr")
png(here::here("data/plots/SpearmanCorrelations/corr_plot_with_pvals.png"), height = 2160, width = 2160)
corrplot(cor_mat,method='ellipse')
mtext("Spearman Correlation Plot - Numbers are p-vals", at=3.5, line=-0.5, cex=4)
pos <- expand.grid(1:ncol(cor_test_mat$p), ncol(cor_test_mat$p):1)
text(pos, p_format(cor_test_mat$p))
dev.off()

# %%

# %%
# pi vals

library(GenomicFeatures)
library(rtracklayer)
# BiocManager::install("GenomicFeatures")

gtf <- makeTxDbFromGFF(here::here("data/annotation/Homo_sapiens.GRCh37.87.gtf")) #change me!
gene_annotations <- genes(gtf)

pi_ceu = import.bw(here::here("data/annotation/Pi_CEU_10kb.bw"))
seqlevels(pi_ceu)<- sub('chr','',seqlevels(pi_ceu))

genes_with_pi <- mergeByOverlaps(gene_annotations,pi_ceu)
results_df = data.frame(gene = character(),pi = numeric())
for(gene in unique(genes_with_pi$gene_id)){
    results = c(gene, mean(genes_with_pi[genes_with_pi$gene_id == gene,]$score))
    results_df[nrow(results_df)+1,] = results
}

write.csv(results_df, here::here("data/annotation/pi_ceu_results.csv"),row.names = F)

# %% 

# %% 

library(psych)
library(corrplot)

rank_df = read.csv(here::here("data/pca_ranks.csv"), header = TRUE)[, -1]
pi_vals = read.csv(here::here("data/annotation/pi_ceu_results.csv"))

merged = merge(rank_df, pi_vals, by.x = "gene", by.y = "gene")
merged_filtered = merged[,which(names(merged) %in% c("gene","mean","sd", "pi"))]
cor_mat = cor(merged_filtered[,-1],method="spearman")
cor_test_mat = corr.test(merged_filtered[,-1],method="spearman",adjust="fdr")
png(here::here("data/plots/SpearmanCorrelations/pi_corr_plot_with_pvals.png"), height = 2160, width = 2160)
corrplot(cor_mat,method='ellipse')
mtext("Spearman Correlation Plot - Numbers are p-vals", at=3.5, line=-0.5, cex=4)
pos <- expand.grid(1:ncol(cor_test_mat$p), ncol(cor_test_mat$p):1)
text(pos, p_format(cor_test_mat$p))
dev.off()
# %% 

# %%
rank_df[order(rank_df$sd,decreasing=TRUE),][1:10,]
low <- grep("ENSG00000136709",rownames(brain$residuals_noOut))
high <- grep("ENSG00000179388",rownames(brain$residuals_noOut))
sd(brain$residuals_noOut[high,])
apply(brain$residuals_noOut,1,sd) %>% summary
#%%