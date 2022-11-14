# Script to calculate enrichment of differentially expressed LCLs across environments in high and low expression genes


# %%
library(here)
library(tidyverse)
library(stringr)
library(AnnotationDbi)
# %%


# %%

organism <- "org.Hs.eg.db"
library(organism, character.only = TRUE)

# %%


# %%
lcl_count_of_genes_per_treatment = list(
    ACRYL = 797,
    BAFF = 195,
    BPA = 717,
    ETOH = 14,
    FSL1 = 102,
    GARD = 700,
    IFNG = 532,
    # IGF = 0,
    PFOA = 801,
    DEX = 3214,
    TUNIC = 60
)
# %%

# %%
# Load the data and convert Ensembl IDs to gene symbols for matching

# LCL enrichment table from https://github.com/AmandaJLea/LCLs_gene_exp/blob/main/main_results/10Jun21_eQTL_SNPs_sharing_EDITED.txt
lcl_table <- read.table(here::here("data/10Jun21_eQTL_SNPs_sharing_EDITED.txt"), header = T, sep = "\t")
uniq_lcl_eqtls <- unique(lcl_table$gene)
rank_df <- read.csv(here::here("data/pca_ranks.csv"), header = TRUE)[, -1]

# %%
run_hypergeometric <- function(fdr){
# %%
library(qvalue)
files <- list.files(path="/Genomics/ayroleslab2/scott/exp_var/data/annotation/Lea_LCL_results", pattern="*.txt", full.names=TRUE, recursive=FALSE)
list_of_dfs <- lapply(files, function(x) {
    df <- read_tsv(x) # load file
    if("P.Value" %in% colnames(df)) {
        df <- df[order(df$P.Value),] # sort by p-value
        # df$qvalue <- qvalue(df$P.Value)$qvalues#p.adjust(df$P.Value, method = "BH")
        # df <- df[df$qvalue < fdr,]
        name = str_split(x, "/")[[1]][9]
        name = str_split(name, "\\.")[[1]][1]
        name = str_split(name, "_")[[1]][4]
        df <- df[1:lcl_count_of_genes_per_treatment[[name]],]

    } else {
        df$gene <- mapIds(org.Hs.eg.db,
        keys = df$gene,
        keytype = "SYMBOL", column = "ENSEMBL"
        )
    }
    df
})
names(list_of_dfs) <- gsub(".txt", "", basename(files))

n_env_responsive_genes <- lapply(list_of_dfs, FUN=function(x) {
    nrow(x)
})
n_env_responsive_genes[['10Jun21_eQTL_SNPs_sharing_EDITED']] <- unique(list_of_dfs[['10Jun21_eQTL_SNPs_sharing_EDITED']])$gene %>% unique %>% length
# %%

# list_of_dfs[["23May22_limma_treatment_ACRYL"]]


# %%
# Select metric and take top and bottom n% of genes

metric <- "sd"
tail_size <- 0.05

upper_quantiles <- list()
lower_quantiles <- list()

cutoff <- quantile(rank_df[[metric]], 1 - tail_size)
subset <- rank_df[rank_df[[metric]] >= cutoff, ]$Gene
upper_quantiles[[metric]] <- subset

cutoff <- abs(quantile(rank_df[[metric]], tail_size))
lower_quantiles[[metric]] <- rank_df[rank_df[[metric]] <= cutoff, ]$Gene
# %%


# %%
# Find overlap between LCL eQTLs and top/bottom n% of genes, calculate odds ratio, and perform hypergeometric test
df <- data.frame(source=character(),
                 ptwas_category = character(),
                 number_of_genes_tissue = integer(),
                 number_of_genes_topbottom = integer(),
                 category_count=integer(),
                 total_count = integer(),
                 upper_count = integer(),
                 lower_count=integer(),
                 num_env_responsive = integer(),
                 fdr=numeric(),
                 stringsAsFactors=F)

# Create dummy var noting overlap with top/bottom n% of genes
for (name in names(list_of_dfs)){
    uniq_lcl_eqtls = unique(list_of_dfs[[name]]$gene)
    rank_df_with_dummyvar <- rank_df %>% mutate(sd_dummyvar = if_else(Gene %in% uniq_lcl_eqtls, 1, 0))
    lower_with_dummyvar <- as.data.frame(lower_quantiles) %>% mutate(sd_dummyvar = if_else(sd %in% uniq_lcl_eqtls, 1, 0))
    upper_with_dummyvar <- as.data.frame(upper_quantiles) %>% mutate(sd_dummyvar = if_else(sd %in% uniq_lcl_eqtls, 1, 0))

    results_df <- data.frame(source = name, number_of_genes_tissue = dim(as.data.frame(rank_df))[1], number_of_genes_topbottom = dim(as.data.frame(upper_quantiles))[1],
                total_count = length(uniq_lcl_eqtls), upper_count = sum(upper_with_dummyvar$sd_dummyvar), lower_count = sum(lower_with_dummyvar$sd_dummyvar), n_env_responsive = n_env_responsive_genes[[name]],fdr=fdr, stringsAsFactors = F)

    # Calculate odds ratios
    results_df$upper_or <- (results_df$upper_count / (results_df$total_count - results_df$upper_count)) / ((100*tail_size) / (100 - (100*tail_size)))
    results_df$lower_or <- (results_df$lower_count / (results_df$total_count - results_df$lower_count)) / ((100*tail_size) / (100 - (100*tail_size)))

    # Perform hypergeometric test
    q = results_df$upper_count
    m = results_df$number_of_genes_topbottom
    n = (results_df$number_of_genes_tissue - results_df$number_of_genes_topbottom)
    k = n_env_responsive_genes[[name]]# Number of environmentally responsive genes

    # Perform hypergeometric test and update results_df
    results_df$phyper_upper <- phyper(q, m, n, k, lower.tail = FALSE)

    q = results_df$lower_count
    results_df$phyper_lower <- phyper(q, m, n, k, lower.tail = FALSE)

    results_df$phyper_upper_bh_adj <- p.adjust(results_df$phyper_upper, method = "BH")
    results_df$phyper_lower_bh_adj <- p.adjust(results_df$phyper_lower, method = "BH")
    df <- rbind(df, results_df)
}

df

# %%
}
output <- run_hypergeometric(0)
write_csv(output, here::here("data/Lea_LCL_results_hypergeometric_test.csv"))
# write_csv(data.frame(output), here::here(paste0("data/annotation/lcl_output_table_fdr",gsub('\\.', '', toString(fdr)),".csv")))
