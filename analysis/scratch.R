source(here::here("functions.R"))
library(here)

ea_df = read.csv(here::here("snakemake/metadata/EA_metadata.csv"), header=T, comment.char = "#")
rc3_df = read.csv(here::here("snakemake/metadata/recount3_metadata.csv"), header=T, comment.char = "#")
metadata_df = bind_rows(ea_df,rc3_df)

md = read.csv(here::here("data/raw_metadata_df.csv"), header=T, comment.char = "#")
sum(md$filtered_metadata_individuals)

pca_ranks = read.csv(here::here("data/pca_ranks.csv"), header=T, comment.char = "#")
num_genes = length(pca_ranks)

sample_metadata = read.csv(here::here("data/sample_metadata.csv"), header=T, comment.char = "#")
# sample_metadata[sample_metadata$filtered_metadata_individuals==min(sample_metadata$filtered_counts_individuals),]
range(sample_metadata$filtered_counts_individuals)
mean(sample_metadata$filtered_counts_individuals)
median(sample_metadata$filtered_counts_individuals)

