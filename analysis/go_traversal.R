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
  gene.data = mget(c(id),org.Hs.egGO2ALLEGS)
  entrezgene = unlist(gene.data)
  genes <- getBM(filters="entrezgene_id", attributes=c("ensembl_gene_id","entrezgene_id"), values=entrezgene, mart=ensembl)
  return(genes)
}


library(GO.db)
library(biomaRt)

go_gene_groups <- lapply(go_terms, call_go)
names(go_gene_groups) <- go_terms
go_gene_overlapping <- list()
for (term in go_terms) {
  go_gene_overlap <- go_gene_groups[[term]][go_gene_groups[[term]]$ensembl_gene_id %in% rank_df$Gene,]
  go_gene_overlap <- go_gene_overlap[!duplicated(go_gene_overlap$ensembl_gene_id),]
  go_gene_overlapping[[term]] <- go_gene_overlap
}

lapply(go_gene_overlapping,dim)
save.image(file='go_traversal.RData')

# %%





