names(exp_data_rc3)

source("functions.R")

filter_meta = function(metadata){
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
    metadata
}


samples = ldply(names(exp_data_rc3)[18:30], function(x) {print(x); nrow(filter_meta(colData(exp_data_rc3[[x]])))})
samples$dataset = names(exp_data_rc3)[18:30]
metadata <- colData(exp_data_rc3[["KIRC"]])
sel_metada = select_meta(filter_meta(metadata))
nrow(sel_metada)

names(metadata)
head(sel_metada)

design = make_design_matrix(sel_metada, crap_cols)
dim(design)
dim(sel_metada)

ldply(sel_metada, function(x) length(unique(x)))
llply(sel_metada, function(x) table(x))

#sapply(metadata[,mask], function(x) length(unique(x)))
nrow(sel_metada)
table(sel_metada)
col = metadata[["tcga.cgc_sample_sample_type"]]
names(dplyr::select(as.data.frame(metadata), contains("type")))
sum(col == "")/length(col)
length(col)
head(col)
length(unique(col))
table(col)
col == ""
