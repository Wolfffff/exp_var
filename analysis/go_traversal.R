# %%
source(here::here("functions.R"))
library(here)
library(sjmisc)
library(pryr)
library(GOxploreR)
# %%

# %%
# GO analysis on ea
  getAllBPChildren <- function(goids)
    {
      ans <- unique(unlist(mget(goids, GOBPCHILDREN), use.names=FALSE))
      ans <- ans[!is.na(ans)]
    }
organism = "org.Hs.eg.db"
metrics = c("mean", "sd")

rank_df = read.csv(here::here("data/pca_ranks.csv"), header = TRUE)[, -1]

# go_terms <- Level2GOTermBP(level = 1, organism = "Human")
go_terms <-  getAllBPChildren("GO:0008150")
# level2_BP_terms <- getAllBPChildren(level1_BP_terms)  # 256 terms
# level3_BP_terms <- getAllBPChildren(level2_BP_terms)  # 3059 terms
# level4_BP_terms <- getAllBPChildren(level3_BP_terms)  # 9135 terms
# level5_BP_terms <- getAllBPChildren(level4_BP_terms)  # 15023 terms
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")


# %%

# %%
call_go <- function(id) {
  tryCatch({
      gene.data = mget(c(id),org.Hs.egGO2ALLEGS)
      entrezgene = unlist(gene.data)
      genes <- getBM(filters="entrezgene_id", attributes=c("ensembl_gene_id","entrezgene_id"), values=entrezgene, mart=ensembl)
      genes
  }, error = function(e) {
    print(e)
    NULL
  })
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
lapply(go_gene_groups,dim)
names(go_gene_groups)
save.image(file='go_traversal.RData')


# %%


# %%
GO <- as.list(GOTERM)
my.term <- GO$[[term]]@Term
# %%




