source("functions.R")

load("Rdatas/results_list_ea.RData")
load("Rdatas/results_list_rc3.RData")

mean_df_list = vector("list", length = length(results_list_rc3))
names(mean_df_list) = names(results_list_rc3)
for (dset_name in names(results_list_rc3)) {
  print(dset_name)
  exprDf = results_list_rc3[[dset_name]]
  gene_vars = rowMeans(exprDf$residuals_noOut)
  mean_df_list[[dset_name]] = data.frame(
      Genes = names(gene_vars),
      vars = gene_vars)
  }
mean_df_rc3 = reduce(mean_df_list, inner_join, by = "Genes")
mean_df_rc3 = tidyr::separate(mean_df_rc3, "Genes", c("Genes", NA))
colnames(mean_df_rc3)[-1] = names(results_list_rc3)

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
variance_df_rc3 = reduce(variance_df_list, inner_join, by = "Genes")
variance_df_rc3 = tidyr::separate(variance_df_rc3, "Genes", c("Genes", NA))
colnames(variance_df_rc3)[-1] = names(results_list_rc3)

mean_df_list = vector("list", length = length(results_list_ea))
names(mean_df_list) = names(results_list_ea)
for (dset_name in names(results_list_ea)) {
  print(dset_name)
  exprDf = results_list_ea[[dset_name]]
  gene_vars = rowMeans(exprDf$residuals_noOut)
  mean_df_list[[dset_name]] = data.frame(
      Genes = names(gene_vars),
      vars = gene_vars)
  }
mean_df_ea = reduce(mean_df_list, inner_join, by = "Genes")
colnames(mean_df_ea)[-1] = names(results_list_ea)

variance_df_list = vector("list", length = length(results_list_ea))
names(variance_df_list) = names(results_list_ea)
for (dset_name in names(results_list_ea)) {
  print(dset_name)
  exprDf = results_list_ea[[dset_name]]
  gene_vars = rowVars(exprDf$residuals_noOut)
  variance_df_list[[dset_name]] = data.frame(
      Genes = names(gene_vars),
      vars = gene_vars)
  }
variance_df_ea = reduce(variance_df_list, inner_join, by = "Genes")
colnames(variance_df_ea)[-1] = names(results_list_ea)

merged_variance_df = inner_join(variance_df_ea, variance_df_rc3, by = "Genes")
merged_mean_df = inner_join(mean_df_ea, mean_df_rc3, by = "Genes")

pak::pkg_install("corrplot")
library(corrplot)
var_cor = cor(merged_variance_df[,-1], method = "s")
png("data/var_corr_plot.png", height = 2080, width = 2080)
corrplot.mixed(var_cor, upper = "ellipse")
dev.off()

mean_var = list(mean = list(ea = mean_df_ea, rc3 = mean_df_rc3, 
                            merged = merged_mean_df),
                var  = list(ea = variance_df_ea, rc3 = variance_df_rc3, 
                            merged = merged_variance_df))
saveRDS(mean_var, file = "Rdatas/mean_variance.RDS")
