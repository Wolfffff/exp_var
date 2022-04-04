# %%
source(here::here("functions.R"))
library(here)
library(sjmisc)
library(pryr)
library(GOxploreR)
# %%

# %%
# GO analysis on ea
organism = "org.Hs.eg.db"
metrics = c("mean", "sd")

rank_df = read.csv(here::here("data/pca_ranks.csv"), header = TRUE)[, -1]

go_terms <- Level2GOTermBP(level = 1, organism = "Human")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")


# %%

# %%
call_go <- function(id) {
  print(id)
  gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'go_id', "name_1006", "namespace_1003"),
                filters = 'go', values = id, mart = ensembl)
  return(gene.data)
}


library(GO.db)
library(biomaRt)
get_go_genes <- function(go_id) {
  print(go_id)
  go_ids <- c(GOBPOFFSPRING[[go_id]], go_id)
  all <- lapply(go_ids, call_go)
  return(all)
}

go_gene_groups <- lapply(go_terms, get_go_genes)
names(go_gene_groups) <- go_terms

go_gene_groups_flat <- lapply(go_gene_groups, rbindlist)
save.image(file='go_traversal_1.RData')
go_gene_overlapping <- list()
for (term in go_terms) {
  go_gene_overlap <- go_gene_groups_flat[[term]][go_gene_groups_flat[[term]]$ensembl_gene_id %in% rank_df[, 1],]
  go_gene_overlap <- go_gene_overlap[!duplicated(go_gene_overlap$ensembl_gene_id),]
  go_gene_overlapping[[term]] <- go_gene_overlap
}
save.image(file='go_traversal_2.RData')

# %%

lapply(go_gene_overlapping,dim)
