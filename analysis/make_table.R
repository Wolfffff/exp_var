library(tidyverse)
library(data.table)
metadata_df <- read_csv("data/raw_metadata_df.csv")
metadata_df$citation <- NA
metadata_df[which(metadata_df$group == "GTEx"),]$citation = "GTEx"
metadata_df[which(metadata_df$group == "TCGA"),]$citation = "TCGA"
metadata_df[which(metadata_df$group %in% c("Other - recount3", "Other - Expression Atlas")),]$group = "Misc"


id_citation_map <- read_csv("data/id_citation_mapping.csv")
map <- id_citation_map$citation
lookup <- setNames(map, id_citation_map$id)
metadata_df[metadata_df$group == "Misc",]$citation <- lookup[unique(metadata_df[metadata_df$group == "Misc",]$id)]
grouped_df = data.table::data.table(metadata_df)[,lapply(.SD, paste, collapse = ","),'group']
grouped_df [1]$citation <- "GTEx"
grouped_df [2]$citation <- "TCGA"
grouped_df[,c("id", "citation")]
