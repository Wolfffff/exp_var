names(exp_data_rc3)

metadata <- colData(exp_data_rc3[["ERP107748"]])
mask_no_levels = sapply(metadata, function(x) length(unique(x)) > 1)
names(metadata)[mask_no_levels]
head(metadata[,mask_no_levels])

col = metadata[["title"]]
length(col)
head(col)
length(unique(col))
table(col) %>% head
col == ""
