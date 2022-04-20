source(here::here("functions.R"))

metric_df = readRDS(here::here("snakemake/Rdatas/gene_metrics.RDS"))
connectivity_df = readRDS(here::here("snakemake/Rdatas/gene_connectivity.RDS"))

#pak::pkg_install(c("corrplot", "vegan", "ape", "Hmisc", "ggrepel", "wesanderson"))
library(corrplot)
library(vegan)
library(ape)
library(Hmisc)
library(ggrepel)
library(wesanderson)

ea_df = read.csv(here::here("snakemake/metadata/EA_metadata.csv"), header=T, comment.char = "#")
rc3_df = read.csv(here::here("snakemake/metadata/recount3_metadata.csv"), header=T, comment.char = "#")
metadata_df = bind_rows(ea_df,rc3_df)

tissue_counts = table(metadata_df$tissue)
tissues = unique(metadata_df$tissue)

for(t in tissues){
    rank_list = list()
    metric_cor_list = list()
    studies = metadata_df[metadata_df$tissue == t, "id"]
    if(tissue_counts[t] > 1){
        for (metric in c("mean", "sd")){
                rank_mat = ldply(metric_df[[metric]][, studies], rank, na.last = "keep")
                rank_mat = t(rank_mat[,-1])
                rownames(rank_mat) = metric_df[[metric]][, 1]
                colnames(rank_mat) = colnames(metric_df[[metric]][, studies])

                metric_matrix = as.matrix(metric_df[[metric]][, studies])
                metric_matrix = metric_matrix[complete.cases(metric_matrix), ]
                metric_cor = rcorr(metric_matrix, type = "spearman")$r

                eig = eigen(metric_cor)
                if(all(eig$vectors[,1] < 0)){
                    eig$vectors[,1] <- -eig$vectors[,1]
                }
                
                PC_scores = as.matrix(rank_mat) %*% eig$vectors

                gene_rank = rank(PC_scores[,1], na.last = "keep")
                names(gene_rank) = metric_df[[metric]][,1]
            } 
            # else{
            #     metric_df[[metric]][, studies] 
            #     gene_rank = rank(metric_df[[metric]][, studies], na.last = "keep")
            # }
            rank_list[[metric]] = gene_rank
        
        rank_df = data.frame(bind_cols(rank_list))
        rank_df$Gene = metric_df[[1]][,1]
        rank_df = relocate(rank_df, Gene)
        rank_df = rank_df[complete.cases(rank_df),]
        print(c(t, tissue_counts[t], nrow(rank_df)))
        write.csv(rank_df, file=here::here(paste0("data/pca_ranks_tissue/", t, ".csv")), row.names = FALSE)
    }
}
