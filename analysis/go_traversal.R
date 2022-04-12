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

go_gene_groups <- lapply(go_terms, call_go)
names(go_gene_groups) <- go_terms
go_gene_overlapping <- list()
for (term in go_terms) {
  go_gene_overlap <- go_gene_groups[[term]][go_gene_groups[[term]]$ensembl %in% rank_df$Gene,]
  go_gene_overlap <- go_gene_overlap[!duplicated(go_gene_overlap$ensembl),]
  term <- GO[[term]]@Term
  go_gene_overlapping[[term]] <- go_gene_overlap
}

ldply(go_gene_overlapping,dim) %>% filter(V1 > 20)
names(go_gene_groups)
save.image(file='go_traversal.RData')



# %%

# %%
library(stringr)

rank_df = read.csv("data/pca_ranks.csv")
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

ldply(go_gene_overlapping[single_quant[,1]], termTable)
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

p = ggplot(sig_terms_df, aes(H, Skew)) + geom_point()
save_plot("test.png", p, base_height = 5)

p = ggplot(sig_terms_df, aes(H, Skew)) + geom_point() + geom_label_repel(label = sig_terms_df$.id, max.overlaps = 27) +
  theme_tufte() +
  xlab("Shannon Entropy") + ylab("Skewness") +
                           theme(plot.title = element_text(size = 30),
                                          axis.title = element_text(size = 18),
                                          axis.text = element_text(size = 18),
                                          strip.text.x = element_text(size = 32)) +
                            geom_hline(yintercept=0, linetype = "dashed") +
                            annotate("text", x = 1.5, y = 1, label = 'bold("Low variation bias")',parse=TRUE) +
                            annotate("text", x = 1.5, y = -1, label = 'bold("High variation bias")', parse = TRUE)
save_plot("test.png", p, base_width = 6.5*2, base_height = 11*0.25*2)
save_plot(here::here("data/plots/GOterm_entropy_by_skewness.png"), p, base_width = 6.5*2, base_height = 11*0.25*2)


p = ggplot(sig_terms_df, aes(x=Skew)) + geom_histogram()
save_plot("test.png", p, base_height = 5)


png("test.png")
plot(sort(goTerm_shannon$V1))
dev.off()

n_plots = 5
library(reshape2)
df2 = ldply(go_gene_overlapping[sig_terms_df$.id[1:n_plots]], termTable) %>% melt %>% mutate(class = "High variation bias")
n_terms = nrow(sig_terms_df)
df1 = ldply(go_gene_overlapping[sig_terms_df$.id[(n_terms-n_plots -4):n_terms]], termTable) %>% melt %>% mutate(class = "Low variation bias")
df1$.id = vector <- sub("^(\\S+) (\\S+) ", "\\1 \\2\n", df1$.id)
df2$.id = vector <- sub("^(\\S+) (\\S+) ", "\\1 \\2\n", df2$.id)
p1 = ggplot(df1, aes(x=.id, y=value, fill=variable)) +
geom_bar(stat="identity", color="black", position=position_dodge()) + facet_wrap(~class, ncol = 1, scale="free") +
  scale_fill_viridis_d(option="inferno", labels = 1:10) + 
  theme_tufte() + theme(axis.text.x = element_text(size=18,angle = 0, hjust = 0.5)) + xlab("") + ylab("Counts") + labs(fill = "Decile") +
                          theme(plot.title = element_text(size = 30),
                                legend.title = element_text(size = 25),
                                legend.text = element_text(size = 15),
                                axis.title = element_text(size = 18),
                                axis.text = element_text(size = 24),
                                strip.text.x = element_text(size = 32))
p2 = ggplot(df2, aes(x=.id, y=value, fill=variable)) +
geom_bar(stat="identity", color="black", position=position_dodge()) + facet_wrap(~class, ncol = 1, scale="free") +
  scale_fill_viridis_d(option="inferno", labels = 1:10) + 
  theme_tufte() + theme(axis.text.x = element_text(size=24,angle = 0, hjust = 0.5)) + xlab("") + ylab("Counts") + labs(fill = "Decile") +
                          theme(plot.title = element_text(size = 30),
                                legend.title = element_text(size = 25),
                                legend.text = element_text(size = 15),
                                axis.title = element_text(size = 18),
                                axis.text = element_text(size = 24),
                                strip.text.x = element_text(size = 32))
save_plot("test.png", p1+p2 + plot_layout(ncol=1, guides = "collect"), base_height = 6.5*2.2,base_width=13*2.2)
save_plot(here::here("data/plots/GOterm_decile_barplot.png"), p1+p2 + plot_layout(ncol=1), base_height = 6.5*2.2,base_width=13*2.2)

# %%

# %%

# %%

# %%
df = data.frame(ensembl=list_of_go_genes)
tx = termTable(df)
df = transpose(data.frame(tx))
df$V2 <- factor(seq_along(df$V1), levels = rownames(df))
p = ggplot(df,aes(x=V2,y=V1))+geom_bar(stat="identity")
save_plot("test.png", p, base_height = 6)
# %%




# %%
# Lets explore the genes that don't appear in GO terms!
rank_df_excl_overlap = rank_df[!(rank_df$Gene %in% list_of_go_genes),]
df = data.frame(ensembl=rank_df_excl_overlap$Gene)
tx = termTable(df)
df = transpose(data.frame(tx))
df$V2 <- factor(seq_along(df$V1), levels = rownames(df))
p = ggplot(df,aes(x=V2,y=V1))+geom_bar(stat="identity")
save_plot("test.png", , base_height = 6)
# %%