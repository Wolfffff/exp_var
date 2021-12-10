names(exp_data_rc3)

crap_cols = c("alias", "Alias", "Broker.name", "broker.name", "Description", "Title", "ENA.checklist", 
              "ENA.FIRST.PUBLIC", "ENA.LAST.UPDATE", "isolate", "INSDC.center.alias", 
              "INSDC.center.name", "INSDC.first.public", "INSDC.last.update", "INSDC.status", 
              "Sample.Name", "SRA.accession", "title", "gtex.smrin")#, 

             #"gtex.smubrid")

dset_name = "ESOPHAGUS"
exp_data = exp_data_rc3
metadata <- colData(exp_data_rc3[["ESOPHAGUS"]])
mask_no_levels = sapply(metadata, function(x) length(unique(x)) <= 1 | length(unique(x)) > 100)
mask_crap_cols = names(metadata) %in% crap_cols
mask = unlist(Map(`|`, mask_no_levels, mask_crap_cols))
names(metadata)[!mask]
head(metadata[,!mask])

sapply(metadata[,!mask, drop = FALSE], function(x) length(unique(x)))
#sapply(metadata[,mask], function(x) length(unique(x)))
nrow(metadata)
table(metadata[,!mask])
col = metadata[["gtex.smtsd"]]
length(col)
head(col)
length(unique(col))
table(col)
col == ""