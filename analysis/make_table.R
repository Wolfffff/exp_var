df1<-read.csv("snakemake/uf_metadata.csv",head=TRUE,row.names=NULL)
df2<-read.csv("snakemake/f_metadata.csv",head=TRUE, row.names=NULL)

merged_df <- merge(df1,df2,by="name",all.x=T, all.y=T)


filtered_df = data.frame(study = merged_df$name, raw_gene_count = merged_df$count_rows.x, raw_sample_count = merged_df$metadata_rows.x,
 filtered_gene_count = merged_df$count_rows.y, filtered_sample_count = merged_df$metadata_rows.y)

ea_ids = read.csv("snakemake/metadata/EA_metadata.csv")
rc3_ids = read.csv("snakemake/metadata/recount3_metadata.csv")

filtered_df$source = ifelse(filtered_df$study %in% ea_ids, "Expression Atlas", "recount3")

filtered_df <- filtered_df[order(filtered_df$source),]


