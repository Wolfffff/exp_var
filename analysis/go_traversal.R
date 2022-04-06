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
level1_BP_terms <- getAllBPChildren("GO:0008150")
level2_BP_terms <- getAllBPChildren(level1_BP_terms)  # 256 terms
level3_BP_terms <- getAllBPChildren(level2_BP_terms)  # 3059 terms
level4_BP_terms <- getAllBPChildren(level3_BP_terms)  # 9135 terms
level5_BP_terms <- getAllBPChildren(level4_BP_terms)  # 15023 terms
levels_1_5 = list(level1_BP_terms, level2_BP_terms, level3_BP_terms, level4_BP_terms, level5_BP_terms)

level = 3
if(level == 1) go_terms <- level1_BP_terms
if(level > 1)  go_terms <- setdiff(levels_1_5[[level]], do.call(c, levels_1_5[1:(level - 1)]))

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

ldply(go_gene_overlapping,dim) %>% filter(V1 > 20)
qldply(go_gene_groups,dim)
names(go_gene_groups)
save.image(file='go_traversal.RData')


# %%


# %%
GO <- as.list(GOTERM)
my.term <- GO$[[term]]@Term
# %%

# %%
library(stringr)

rank_df = read.csv("data/pca_ranks.csv")
n_classes = 4
quantile_lower = as.vector(quantile(rank_df$sd, seq(0, 1, length.out = n_classes + 1))[1:n_classes])

classifyQuantileFast = function(x){
  laply(x, function(rank) paste0("quantile_", str_pad(sum(rank >= quantile_lower), 2, pad = 0)))
}
rank_class_df = tibble(gene = rank_df$Gene, 
                       quantile = classifyQuantileFast(rank_df$sd))

shannon <- function(x) -sum(((x <- na.omit(x[x!=0]))/sum(x)) * log(x/sum(x)))

termTable <- function(x){
  out = as.data.frame(matrix(0, ncol = n_classes, nrow = 1))
  names(out) = paste0("quantile_", str_pad(1:n_classes, 2, pad = 0))
  counts = table(rank_class_df$quantile[match(x$ensembl_gene_id, rank_class_df$gene)])
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
mask = filter(ldply(go_gene_overlapping,dim), V1 > 20)$`.id`
goTerm_shannon = ldply(go_gene_overlapping[mask], shannonGOterm)
png("test.png")
hist(goTerm_shannon$V1)
dev.off()

ldply(go_gene_overlapping[single_quant[,1]], termTable)
sig_terms_df = ldply(go_gene_overlapping[mask], 
                     function(x) c(p.value = chiSqTest(x), 
                                   N = nrow(x), 
                                   H = shannonGOterm(x))) %>% 
  mutate(p.adjusted = p.adjust(p.value),
         significant = p.adjusted < 0.01) %>% 
  arrange(H) %>% 
  filter(N > 100) %>%
  as.tibble


png("test.png")
plot(sort(goTerm_shannon$V1))
dev.off()

n_plots = 30
library(reshape2)
df2 = ldply(go_gene_overlapping[sig_terms_df$.id[1:n_plots]], termTable) %>% melt %>% mutate(class = "top")
n_terms = nrow(sig_terms_df)
df1 = ldply(go_gene_overlapping[sig_terms_df$.id[(n_terms-n_plots):n_terms]], termTable) %>% melt %>% mutate(class = "bottom")
p = ggplot(data=rbind(df1, df2), aes(x=.id, y=value, fill=variable)) +
geom_bar(stat="identity", color="black", position=position_dodge()) + facet_wrap(~class, ncol = 1, scale="free") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                          legend.position = "none") 
save_plot("test.png", p, base_height = 6)

p = ggplot(sig_terms_df, aes(N, H, color = significant)) + geom_point()
save_plot("test.png", p, base_height = 6)

# %%