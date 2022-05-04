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
load('go_traversal.RData')


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

entropy_by_skewness = ggplot(sig_terms_df, aes(H, Skew)) + 
  geom_point() + 
  geom_label_repel(label = sig_terms_df$.id, max.overlaps = 27) +
  theme_tufte() +
  xlab("Shannon Entropy") + ylab("Skewness") +
                           theme(plot.title = element_text(size = 30),
                                          axis.title = element_text(size = 18),
                                          axis.text = element_text(size = 18),
                                          strip.text.x = element_text(size = 32)) +
                            geom_hline(yintercept=0, linetype = "dashed") +
                            annotate("text", x = 1.5, y = 1, label = 'bold("Low variation bias")',parse=TRUE) +
                            annotate("text", x = 1.5, y = -1, label = 'bold("High variation bias")', parse = TRUE)
save_plot("test.png", entropy_by_skewness, base_width = 6.5*2, base_height = 11*0.25*2)
save_plot(here::here("data/plots/GOterm_entropy_by_skewness.png"), entropy_by_skewness, base_width = 6.5*2, base_height = 11*0.25*2)


p = ggplot(sig_terms_df, aes(x=Skew)) + geom_histogram()
save_plot("test.png", p, base_height = 5)


png("test.png")
plot(sort(goTerm_shannon$V1))
dev.off()

n_plots = 5
library(reshape2)
df2 = ldply(go_gene_overlapping[sig_terms_df$.id[1:n_plots]], termTable) %>% melt %>% mutate(class = "High variation bias")
n_terms = nrow(sig_terms_df)
df1 = ldply(go_gene_overlapping[sig_terms_df$.id[(n_terms-n_plots+1):n_terms]], termTable) %>% melt %>% mutate(class = "Low variation bias")
df1$.id = vector <- sub("^(\\S+) (\\S+) ", "\\1 \\2\n", df1$.id)
df2$.id = vector <- sub("^(\\S+) (\\S+) ", "\\1 \\2\n", df2$.id)
p1 = ggplot(df1, aes(x=.id, y=value, fill=variable)) +
geom_bar(stat="identity", color="black", position=position_dodge()) + 
  scale_fill_viridis_d(option="inferno", labels = 1:10) + ggtitle("A. Low variation bias") +
  theme_tufte() + xlab("") + ylab("Counts") + labs(fill = "Decile") +
                          theme(plot.title = element_text(size = 50),
                                legend.title = element_text(size = 35),
                                legend.text = element_text(size = 32),
                                axis.title = element_text(size = 32),
                                axis.text = element_text(size = 28),
                                axis.text.x = element_text(size=40,angle = 0, hjust = 0.5),
                                strip.text.x = element_text(size = 32))
p2 = ggplot(df2, aes(x=.id, y=value, fill=variable)) +
geom_bar(stat="identity", color="black", position=position_dodge()) + 
  scale_fill_viridis_d(option="inferno", labels = 1:10) + ggtitle("B. High variation bias") +
  theme_tufte() + xlab("") + ylab("Counts") + labs(fill = "Decile") +
                          theme(plot.title = element_text(size = 50),
                                legend.title = element_text(size = 35),
                                legend.text = element_text(size = 32),
                                axis.title = element_text(size = 32),
                                axis.text = element_text(size = 28),
                                axis.text.x = element_text(size=40,angle = 0, hjust = 0.5),
                                strip.text.x = element_text(size = 32))
p = p1+ p2 + plot_layout(ncol=1, guides = "collect")
save_plot("test.png", p, base_height = 6.5*2.2,base_width=13*2.2)
save_plot(here::here("data/plots/GOterm_decile_barplot.png"), p, base_height = 6.5*2.2,base_width=13*2.2)
high_term = p2
low_term = p1
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

# %%
source(here::here("functions.R"))
library(here)
library(sjmisc)
library(pryr)
# %%

# %%
# GO analysis on ea
organism = "org.Hs.eg.db"
metrics = c("mean", "sd")

rank_df = read.csv(here::here("data/pca_ranks.csv"), header = TRUE)[, -1]

tail_size = 0.05

upper_quantiles = list()
lower_quantiles = list()
for(metric in c("mean","sd")){
    cutoff = quantile(rank_df[[metric]],1 - tail_size)
    subset = rank_df[rank_df[[metric]] >= cutoff,]$Gene
    upper_quantiles[[metric]] = subset
    cutoff = abs(quantile(rank_df[[metric]], tail_size))
    lower_quantiles[[metric]] = rank_df[rank_df[[metric]] <= cutoff,]$Gene
}
# %%

# %%
# GO analysis against our list of genes
metric = "sd"
library(enrichplot)

local_go_upper = enrichGO(gene  = upper_quantiles[[metric]],
                          universe = rank_df$Gene,
                          OrgDb         = org.Hs.eg.db,
                          keyType       = 'ENSEMBL',
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.001,
                          qvalueCutoff  = 0.001,
                          readable      = TRUE)
write_csv(local_go_upper@result |> filter(p.adjust < 0.01, Count > 4), here::here("data/annotation/go_upper_quantile.csv"))
pw_upper <- pairwise_termsim(local_go_upper) 
pw_upper <- simplify(pw_upper, cutoff=0.7, by="p.adjust", select_fun=min)
plot_upper <- emapplot(pw_upper, showCategory = 10, cex_label_category = 1.2) + 
      ggtitle("A. High variation") +
          theme_tufte() + 
          theme(legend.position = "none") + 
          theme(plot.title = element_text(size=28),
                axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank()) 
plot_upper$data$color = mean(plot_upper$data$color)
local_go_lower = enrichGO(gene  = lower_quantiles[[metric]],
                          universe = rank_df$Gene,
                          OrgDb         = org.Hs.eg.db,
                          keyType       = 'ENSEMBL',
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.001,
                          qvalueCutoff  = 0.001,
                          readable      = TRUE)
write_csv(local_go_lower@result |> filter(p.adjust < 0.01, Count > 4), here::here("data/annotation/go_lower_quantile.csv"))
pw_lower <- pairwise_termsim(local_go_lower) 
pw_lower <- simplify(pw_lower, cutoff=0.7, by="p.adjust", select_fun=min)
plot_lower <- emapplot(pw_lower, showCategory = 10, cex_label_category = 1.2) + 
    ggtitle("B. Low variation") +            
          theme_tufte() + 
          theme(legend.position = "none") + 
          theme(plot.title = element_text(size=28),
                axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank()) 
plot_lower$data$color = mean(plot_lower$data$color)
#p2$data$name[3] = "RNA splicing, via transesterification\n reactions with bulged adenosine as nucleophile"

p12 = plot_upper + plot_lower + plot_layout(ncol=2)

save_plot(here::here("data/plots/local_go_lowerUpper.png"), p12, base_height = 6, base_asp = 2)
# %%


# %%
{layout <- 
"AABB
CCCC
DDDD"}
fig3 =  plot_upper + 
        plot_lower + 
        entropy_by_skewness + ggtitle("C. Entropy by Skewness") + 
        (p1+ p2 + plot_layout(ncol=1, guides = "collect")) + plot_layout(design = layout)

save_plot("test.png", fig3, base_height = 10*2.2,base_width=10*2.2)
# %%