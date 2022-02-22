# Script for analysing the drivers of correlation strucuture

source(here::here("functions.R"))

metric_df = readRDS(here::here("snakemake/Rdatas/gene_metrics.RDS"))
ea_df = read.csv(here::here("snakemake/metadata/EA_metadata.csv"), header=T)
rc3_df = read.csv(here::here("snakemake/metadata/recount3_metadata.csv"), header=T)
metadata_df = bind_rows(ea_df,rc3_df)

corr_mat = readRDS(here::here("snakemake/Rdatas/gene_var_matrices.RDS"))$sd

ids = rownames(corr_mat)
(z1 = outer(ids, ids, paste, sep = ":"))
grp = metadata_df[match(ids, metadata_df$id), "group"]
(z2 = outer(grp, grp, paste, sep = ":"))

lt = function(x) x[lower.tri(x)]
sort_label = function(x){
    sapply(x, function(x){
    strs = str_split(x, ":")[[1]]
    paste(str_sort(strs), collapse = ":")
    })
}
corr_df = tibble(corr = lt(corr_mat), 
                 pair = sort_label(lt(z1)), 
                 source = sort_label(lt(z2))) %>%
    mutate(pair2 = pair) %>%
    separate(pair2, c("Study1", "Study2"), sep = ":")



