source(here::here("functions.R"))

data = readRDS(here::here("snakemake/Rdatas/raw/SRP162654.rds"))

metadata <- colData(data)
mask_no_levels = sapply(metadata, function(x) length(unique(x)) <= 1)
mask_crap_cols = names(metadata) %in% crap_cols
mask = unlist(Map(`|`, mask_no_levels, mask_crap_cols))
names(metadata)[!mask]
head(metadata[,!mask]) 

col = metadata[["donor.id"]]
length(col)
head(col)
length(unique(col))
table(col)
col == ""