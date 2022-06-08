library(tidyverse)
library(data.table)
metadata_df <- read_csv("data/raw_metadata_df.csv")
metadata_df$bibkey <- NA
metadata_df$citation <- NA
metadata_df[which(metadata_df$group == "GTEx"),]$citation = "GTEx"
metadata_df[which(metadata_df$group == "TCGA"),]$citation = "TCGA"
metadata_df[which(metadata_df$group %in% c("Other - recount3", "Other - Expression Atlas")),]$group = "Misc"
metadata_df <- metadata_df[order(metadata_df$id),]


id_citation_map <- read_csv("data/id_citation_mapping.csv")

library(stringr)


map <- id_citation_map$citation
lookup <- setNames(map, id_citation_map$id)


metadata_df[metadata_df$group == "GTEx",]$bibkey <- lookup["GTEx"]
metadata_df[metadata_df$group == "TCGA",]$bibkey <- lookup["TCGA"]
metadata_df[metadata_df$group == "Misc",]$bibkey <- lookup[metadata_df[metadata_df$group == "Misc",]$id]

metadata_df$citation <- paste0("@",metadata_df$bibkey)

metadata_df$id <- paste(metadata_df$id, paste0(" (", str_to_title(metadata_df$tissue),")"), sep="")

library(bibtex)
library(bib2df)
bib <- bib2df("data/references.bib")
bib <- bib[!duplicated(bib$BIBTEXKEY),]
bib$BIBTEXKEY <- as.character(bib$BIBTEXKEY)
metadata_df$bibkey <- as.character(metadata_df$bibkey)
bib = bib[which(bib$BIBTEXKEY %in% metadata_df$bibkey),]
author_lastname <- unlist(lapply(bib$AUTHOR, function(x) {
  strsplit(x,",")[[1]][1]
}))

years <- unlist(lapply(bib$YEAR, function(x) {
  strsplit(x,",")[[1]][1]
}))

bib$verbose_citation <- paste0(author_lastname, " et al., ", years)
metadata_df <- merge(metadata_df, bib, by.x="bibkey", by.y = "BIBTEXKEY")
metadata_df$citation <- paste0(metadata_df$verbose_citation, " - ", metadata_df$citation)
grouped_df <- data.table::data.table(metadata_df)[,lapply(.SD, paste, collapse = ", "),'group']


grouped_df = grouped_df[c(2,3,1),]

grouped_df [1]$citation <-paste0("The GTEx Consortium, 2020 - @", lookup["GTEx"])
grouped_df [2]$citation <- paste0("The Cancer Genome Atlas Research Network et al., 2013 - @", lookup["TCGA"])
grouped_df <- grouped_df[1:2,c("id","citation")]
grouped_df <- rbind(grouped_df, metadata_df[metadata_df$group == "Misc",c("id", "citation")])

library(kableExtra)
tbl = knitr::kable(grouped_df,"pipe")
writeLines(tbl, 'data/T1.md'')
