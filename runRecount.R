source("functions.R")

experimental_metadata_ea <- read_csv("expression_atlas_metadata.csv")
# Remove the top row (GTEx) to lower the time cost when testing
experimental_metadata_ea <- experimental_metadata_ea[-1, ]

# Load metadata from recount3
experimental_metadata_rc3 <- read_csv("recount3_metadata.csv")

plots_dir <- "data/"

# Setup filters for removing samples accordingly to the metadata
feature_vec <- list()
feature_vec[["disease"]] <- c("normal", "control", "", NA,
                              "non inflammatory bowel disease control")
feature_vec[["treatment"]] <- c("normal", "control", "", NA)

# Move cache to deal with quota issues
cache <- recount3_cache(cache_dir = "cache")
human_projects <- available_projects(bfc = cache)

source("./main_processing_loop.R")
parallel <- FALSE
if (.Platform$OS.type == "unix") {
  parallel <- TRUE
  library(doMC)
  registerDoMC(64)
}

pull_data <- TRUE
if (pull_data) {
  tic()
  exp_data_rc3 = llply(experimental_metadata_rc3$id, downloadRecount3, .parallel = parallel)
  names(exp_data_rc3) = experimental_metadata_rc3$id
  toc()
  #saveRDS(exp_data_rc3, file = "cache/recount3_data.RDS")
}
#exp_data_rc3 <- readRDS("cache/recount3_data.RDS")

results_list_rc3 <- llply(names(exp_data_rc3),
                      main_loop,
                      exp_data = exp_data_rc3,
                      experimental_metadata = experimental_metadata_rc3,
                      feature_vec = feature_vec,
                      assay_name = "raw_counts",
                      .parallel = parallel)
names(results_list_rc3) <- names(exp_data_rc3)

names(results_list_rc3)[(is.na(results_list_rc3))]
# re_run = "SRP216404"
# dset_name = re_run
# results_list_rc3[[re_run]] <- main_loop(dset_name = re_run, exp_data = exp_data_rc3,
#                       experimental_metadata = experimental_metadata_rc3,
#                       feature_vec = feature_vec,
#                       assay_name = "raw_counts")

save(results_list_rc3, file = "Rdatas/results_list_rc3.RData")
#load("Rdatas/results_list_rc3.RData")
for (dset_name in names(results_list_rc3)) {
    if(!is.na(results_list_rc3[[dset_name]]))
        save_plot(filename = paste0(plots_dir, dset_name, "_pca.png"),
                  results_list_rc3[[dset_name]]$plotPanel,
                  base_height = 6, base_asp = 1.2, ncol = 2, nrow = 2)
}

variance_df_list = vector("list", length = length(results_list_rc3))
names(variance_df_list) = names(results_list_rc3)
for (dset_name in names(results_list_rc3)) {
  print(dset_name)
  exprDf = results_list_rc3[[dset_name]]
  gene_vars = rowVars(exprDf$residuals_noOut)
  variance_df_list[[dset_name]] = data.frame(
      Genes = names(gene_vars),
      vars = gene_vars)
  }
variance_df = reduce(variance_df_list, inner_join, by = "Genes")
colnames(variance_df)[-1] = names(results_list_rc3)

pak::pkg_install("corrplot")
library(corrplot)
var_cor = cor(variance_df[,-1], method = "s")
png("data/var_corr_plot_rc3.png", height = 1080, width = 1080)
corrplot.mixed(var_cor, upper = "ellipse")
dev.off()
