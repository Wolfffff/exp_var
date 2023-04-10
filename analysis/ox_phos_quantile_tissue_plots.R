# %%
source(here::here("functions.R"))
library(ggthemes)
library(ggrepel)
library(here)
library(sjmisc)
library(pryr)
library(GOxploreR)
library(moments)
library(data.table)
library("AnnotationDbi")
library("org.Hs.eg.db")
GO <- as.list(GOTERM)


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

# go_terms <- Level2GOTermBP(level = 1, organism = "Human")
level1_BP_terms <- getAllBPChildren("GO:0008150")
level2_BP_terms <- getAllBPChildren(level1_BP_terms)  # 256 terms
level3_BP_terms <- getAllBPChildren(level2_BP_terms)  # 3059 terms
level4_BP_terms <- getAllBPChildren(level3_BP_terms)  # 9135 terms
level5_BP_terms <- getAllBPChildren(level4_BP_terms)  # 15023 terms
levels_1_5 = list(level1_BP_terms, level2_BP_terms, level3_BP_terms, level4_BP_terms, level5_BP_terms)

level = 3
if(level == 1) go_terms <- level1_BP_terms
if(level > 1)  go_terms <- setdiff(levels_1_5[[level]], do.call(c, levels_1_5[1:(level - 1)]))

# ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

# %%

# %%
call_go <- function(id) {
  tryCatch({
      gene.data = mget(c(id),org.Hs.egGO2ALLEGS)
      entrezgene = unlist(gene.data)
      genes = data.frame(entrez=entrezgene)
      genes$ensembl = mapIds(org.Hs.eg.db,
                       keys=genes$entrez, 
                        column="ENSEMBL",
                        keytype="ENTREZID",
                        multiVals="first")
      genes
  }, error = function(e) {
    print(e)
    NULL
  })
}


library(GO.db)
library(biomaRt)

# go_gene_groups <- lapply(go_terms, call_go)
# names(go_gene_groups) <- go_terms
# go_gene_overlapping <- list()
# for (term in go_terms) {
#   go_gene_overlap <- go_gene_groups[[term]][go_gene_groups[[term]]$ensembl %in% rank_df$Gene,]
#   go_gene_overlap <- go_gene_overlap[!duplicated(go_gene_overlap$ensembl),]
#   term <- GO[[term]]@Term
#   go_gene_overlapping[[term]] <- go_gene_overlap
# }

# ldply(go_gene_overlapping,dim) %>% filter(V1 > 20)
# names(go_gene_groups)
# save.image(file='go_traversal.RData')
load(here::here('go_traversal.RData'))


# %%

# %%
library(stringr)
metrics = c("mean", "sd")

tissues = tools::file_path_sans_ext(dir(here::here("data/pca_ranks_tissue")))
tissue_rank_files = dir(here::here("data/pca_ranks_tissue"), full.names = T)
rank_df_list = lapply(tissue_rank_files, \(x) read.csv(x, header = TRUE))
names(rank_df_list) = tissues


runGOtissue = function(current_tissue = NULL){
    
  rank_df = rank_df_list[[current_tissue]]
  complete_list <- ldply(go_gene_overlapping,rbind)
  list_of_go_genes = unique(complete_list$ensembl)
  rank_df_overlap = rank_df[rank_df$Gene %in% list_of_go_genes,]
  rank_df_overlap <- rank_df_overlap %>% arrange(sd)
  rank_df_overlap$sd <- seq(1,nrow(rank_df_overlap))

  # Rerank based only on genes represented in
  rank_df <- rank_df_overlap

  n_classes = 10
  quantile_lower = as.vector(quantile(rank_df$sd, seq(0, 1, length.out = n_classes + 1))[1:n_classes])

  classifyQuantileFast = function(x){
    laply(x, function(rank) paste0("quantile_", str_pad(sum(rank >= quantile_lower), 2, pad = 0)))
  }
  rank_class_df = tibble(gene = rank_df$Gene, quantile = classifyQuantileFast(rank_df$sd))

  shannon <- function(x) -sum(((x <- na.omit(x[x!=0]))/sum(x)) * log(x/sum(x)))

  termTable <- function(x){
    out = as.data.frame(matrix(0, ncol = n_classes, nrow = 1))
    names(out) = paste0("quantile_", str_pad(1:n_classes, 2, pad = 0))
    counts = table(rank_class_df$quantile[match(x$ensembl, rank_class_df$gene)])
    out = out + counts[names(out)]
    out[is.na(out)] = 0
    out
  } 
  chiSqTest = function(x){
    obs = termTable(x) 
    chisq.test(obs)$p.value
  }

  x = go_gene_overlapping[[1]]
  shannonGOterm = function(x){
    tx = termTable(x)
    shannon(tx/sum(tx))
  }
  skewnessGOterm = function(x){
    tx = termTable(x)
    vec = transpose(tx)$V1
    vec_repped <- rep(seq(1, length(vec)), vec)
    skewness(vec_repped)
  }

  mask = filter(ldply(go_gene_overlapping,dim), V1 > 20)$`.id`
  goTerm_shannon = ldply(go_gene_overlapping[mask], shannonGOterm)

  sig_terms_df = ldply(go_gene_overlapping[mask], 
                      function(x) c(p.value = chiSqTest(x), 
                                    N = nrow(x), 
                                    H = shannonGOterm(x),
                                    Skew = skewnessGOterm(x))) %>% 
    mutate(p.adjusted = p.adjust(p.value),
          significant = p.adjusted < 0.01) %>% 
    arrange(Skew) %>% 
    filter(N > 100) %>%
    as.tibble
    n_terms = nrow(sig_terms_df)
  df1 = ldply(go_gene_overlapping, termTable) %>% melt %>% mutate(class = "Low-variance bias")
  df1 = df1[df1$.id == "oxidative phosphorylation",]
  p1 = ggplot(df1, aes(x=.id, y=value, fill=variable)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) + 
    scale_fill_viridis_d(option="inferno", labels = 1:10) + ggtitle(str_to_title(current_tissue)) +
    theme_tufte() + xlab("") + ylab("Counts") + labs(fill = "Decile") +
                            theme(plot.title = element_text(size = 50),
                                  legend.title = element_text(size = 50),
                                  legend.text = element_text(size = 50),
                                  axis.title = element_text(size = 50),
                                  axis.text = element_text(size = 30),
                                  axis.text.x = element_text(size=45,angle = 0, hjust = 0.5),
                                  strip.text.x = element_text(size = 32))
  ggsave(here::here(paste0("plots/ox_phos_sd_dist_",current_tissue,".png")), p1, width = 20, height = 10)
}

library(doMC)
registerDoMC(length(tissues))
llply(tissues, runGOtissue, .parallel = TRUE)