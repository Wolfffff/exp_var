source(here::here("functions.R"))
library(here)

file_paths <- list.files(path = here::here("snakemake/Rdatas/residuals/"), 
                         pattern = "\\.rds", full.names = TRUE)
file_names <-  gsub(pattern = "\\.rds$", replacement = "", x = basename(file_paths))

library(plyr)
library(doMC)
registerDoMC(64)
data_list <- llply(file_paths, readRDS, .parallel = TRUE)

names(data_list) <- file_names

rank_df = read.csv(here::here("data/pca_ranks.csv"), header = TRUE)[, -1]


# %%
organism = "org.Hs.eg.db"

library(organism, character.only = TRUE)
library(clusterProfiler)
library(plyr)
library(dplyr)

rank_df = read.csv(here::here("data/pca_ranks.csv"), header = TRUE)[, -1]

top_quantiles = list()
metric_list = c("means","var","sd","cv")
for(metric in metric_list){
    cutoff = quantile(rank_df[[metric]], .95)
    subset = rank_df[rank_df[[metric]] >= cutoff,]$gene
    top_quantiles[[metric]] = subset
}
term2gene_df = ptwas_table[, c("Trait","Gene")]
ptwas_table_merged = merge(term2gene_df,ptwas_traits, by.x = "Trait", by.y = "ID", all.x = TRUE)
# term2gene_df = data.frame(disease="1",

ego = enrichGO(gene  = rank_df$gene,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05,
                         readable      = TRUE)
library(enrichplot)
png(here::here("data/plots/eGO_universe.png"), height = 3840, width = 2160)
barplot(ego, showCategory=100) 
dev.off()

png(here::here("data/plots/cv_sd.png"), height = 2160, width = 2160)
ggplot(rank_df, aes(x = cv, y = sd)) + geom_point()
dev.off()

library(ape)

gff <- read.gff(here::here("data/annotation/Homo_sapiens.GRCh37.87.gff3"), na.strings = c(".", "?"), GFF3 = TRUE)

source("https://bioconductor.org/biocLite.R")
biocLite("Homo.sapiens")
library(Homo.sapiens)

# ensembl=useMart("ensembl")
library("biomaRt")
# listMarts()
mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
elocs <- getBM(attributes=getAtt,filters="hgnc_symbol",value="ZMYM4",m
art=ensembl)

library(gread)
gff_file <- file.path(here::here("data/annotation/Homo_sapiens.GRCh37.87.gff3"))
gtf <- read_gff(gff_file)
ans <- construct_introns(gtf, update=TRUE)[] # default
# same as above, but returns only constructed introns
introns <- construct_introns(gtf, update=FALSE)

library('GenomicFeatures')


library(gread)
seqs = biomaRt::getSequence(id = Ensembl_IDs, 
           type="ensembl_gene_id",
           seqType = "gene_exon", 
           mart = mart)

# %%
# 
library(GenomicFeatures)
library(rtracklayer)
# BiocManager::install("GenomicFeatures")

gtf <- makeTxDbFromGFF(here::here("data/annotation/Homo_sapiens.GRCh37.87.gtf")) #change me!
gene_annotations <- genes(gtf)
# exons <- exonsBy(gtf, by="gene")

# #make introns
# exons <- reduce(exons)
# exons <- exons[sapply(exons, length) > 1]

# introns <- lapply(exons, function(x) {
#     #Make a "gene" GRange object
#     gr = GRanges(seqnames=seqnames(x)[1], ranges=IRanges(start=min(start(x)),
#         end=max(end(x))), 
#         strand=strand(x)[1])
#     db = disjoin(c(x, gr))
#     ints = db[countOverlaps(db, x) == 0]
#     #Add an ID
#     if(as.character(strand(ints)[1]) == "-") {
#         ints$exon_id = c(length(ints):1)
#     } else {
#         ints$exon_id = c(1:length(ints))
#     }
#     ints
# })
# introns <- GRangesList(introns)

# %%