# %%
# 
library(ape)

gff <- read.gff(here::here("data/annotation/Homo_sapiens.GRCh37.87.gff3"), na.strings = c(".", "?"), GFF3 = TRUE)

source("https://bioconductor.org/biocLite.R")
biocLite("Homo.sapiens")
library(Homo.sapiens)

# ensembl=useMart("ensembl")
library("biomaRt")
# listMarts()

mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")

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

# %%
# 
library(GenomicFeatures)
library(rtracklayer)
# BiocManager::install("GenomicFeatures")

gtf <- makeTxDbFromGFF(here::here("data/annotation/Homo_sapiens.GRCh37.87.gtf")) #change me!
gene_annotations <- genes(gtf)
exons <- exonsBy(gtf, by="gene")

#make introns
exons <- reduce(exons)
exons <- exons[sapply(exons, length) > 1]

introns <- lapply(exons, function(x) {
    #Make a "gene" GRange object
    gr = GRanges(seqnames=seqnames(x)[1], ranges=IRanges(start=min(start(x)),
        end=max(end(x))), 
        strand=strand(x)[1])
    db = disjoin(c(x, gr))
    ints = db[countOverlaps(db, x) == 0]
    #Add an ID
    if(as.character(strand(ints)[1]) == "-") {
        ints$exon_id = c(length(ints):1)
    } else {
        ints$exon_id = c(1:length(ints))
    }
    ints
})
introns <- GRangesList(introns)

# %%