names(exp_data_rc3)

crap_cols = c("alias", "Alias", "Broker.name", "broker.name", "Description", "Title", "ENA.checklist", 
              "ENA.FIRST.PUBLIC", "ENA.LAST.UPDATE", "isolate", "INSDC.center.alias", 
              "INSDC.center.name", "INSDC.first.public", "INSDC.last.update", "INSDC.status", 
              "Sample.Name", "SRA.accession", "title", "gtex.smrin", "rownames", "tcga.xml_month_of_form_completion", 
              "tcga.xml_year_of_form_completion", "tcga.xml_year_of_initial_pathologic_diagnosis", 
              "tcga.xml_initial_pathologic_diagnosis_method", "tcga.cgc_case_year_of_diagnosis", 
              "tcga.gdc_metadata_files.file_size.analysis", "tcga.xml_breast_carcinoma_surgical_procedure_name", 
              "tcga.xml_day_of_form_completion", "tcga.cgc_sample_shortest_dimension", "tcga.xml_stage_event_system_version",
              "tcga.gdc_cases.samples.portions.analytes.concentration", "tcga.gdc_cases.samples.portions.analytes.aliquots.concentration")
crap_cols = c(crap_cols, grep("^tcga..*_percent", names(metadata), value = TRUE))
crap_cols = c(crap_cols, grep("^tcga..*analytes", names(metadata), value = TRUE))
crap_cols = c(crap_cols, grep("^tcga..*xml_has_", names(metadata), value = TRUE))
crap_cols = c(crap_cols, grep("^tcga..*_dimension$", names(metadata), value = TRUE))
crap_cols = c(crap_cols, grep("^tcga..*tumor", names(metadata), value = TRUE))
crap_cols = c(crap_cols, grep("^tcga..*cancer", names(metadata), value = TRUE))
crap_cols = c(crap_cols, grep("^tcga..*pathologic", names(metadata), value = TRUE))
crap_cols = c(crap_cols, grep("^tcga..*_collection_indicator$", names(metadata), value = TRUE))


feature_vec <- list()
feature_vec[["disease"]] <- c("normal", "control", "", NA, 
                              "non inflammatory bowel disease control")
feature_vec[["treatment"]] <- c("normal", "control", "", NA)
feature_vec[["tcga.cgc_sample_sample_type"]] <- c("Solid Tissue Normal")
feature_vec[["diagnosis"]] <- c("Control")
feature_vec[["Healthy"]] <- c("Healthy")

metadata <- colData(exp_data_rc3[[20]])

for (column in names(feature_vec)) {
    if (column %in% colnames(metadata)){
      control_names_tb = table(metadata[,column], useNA = "ifany")
      control_names_tb = control_names_tb[names(control_names_tb) %in% feature_vec[[column]]]
      if(is.na(names(control_names_tb))){
        metadata[,column][is.na(metadata[,column])] = "control"
        names(control_names_tb) = "control"
      }
      control_name = names(control_names_tb)[control_names_tb == max(control_names_tb)]
      metadata <- metadata[metadata[,column] == control_name,]
    }
  }
 

mask_no_levels = sapply(metadata, function(x) length(unique(x)) <= 1 | length(unique(x)) >= 20)
mask_crap_cols = names(metadata) %in% crap_cols
mask_dominant_level = laply(metadata, function(x) max(table(x)) > length(x)*(3/4))
mask_small_levels = laply(metadata, function(x) sum(table(x) <= 2) >= 3)
mask_high_missing = laply(metadata, function(x) sum(x == "")/length(x)) > 0.1
mask = rowSums(cbind(mask_no_levels, mask_crap_cols, mask_high_missing, mask_dominant_level, mask_small_levels)) >= 1
sel_metada = metadata[,!mask, drop = FALSE]

names(sel_metada)
head(sel_metada)

design = make_design_matrix(sel_metada, crap_cols)
dim(design)
dim(sel_metada)

ldply(sel_metada, function(x) length(unique(x)))
llply(sel_metada, function(x) table(x))

#sapply(metadata[,mask], function(x) length(unique(x)))
nrow(sel_metada)
table(sel_metada)
col = sel_metada[["tcga.gdc_cases.samples.portions.analytes.a260_a280_ratio"]]
sum(col == "")/length(col)
length(col)
head(col)
length(unique(col))
table(col)
col == ""
