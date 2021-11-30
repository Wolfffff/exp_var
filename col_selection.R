names(exp_data_rc3)

crap_cols = c("alias", "Alias", "Broker.name", "broker.name", "Description", "Title", "ENA.checklist", "ENA.FIRST.PUBLIC", "ENA.LAST.UPDATE", "isolate", "INSDC.center.alias", 
"INSDC.center.name", "INSDC.first.public", "INSDC.last.update", "INSDC.status", "Sample.Name", "SRA.accession", "title")

metadata <- colData(exp_data_rc3[["ERP017126"]])
mask_no_levels = sapply(metadata, function(x) length(unique(x)) <= 1)
mask_crap_cols = names(metadata) %in% crap_cols
mask = unlist(Map(`|`, mask_no_levels, mask_crap_cols))
names(metadata)[!mask]
head(metadata[,!mask])

col = metadata[["cell.gating"]]
length(col)
head(col)
length(unique(col))
table(col) %>% head
col == ""
