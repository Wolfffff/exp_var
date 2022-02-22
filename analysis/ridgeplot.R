
# %%
library(here)
library(sjmisc)
library(pryr)
source(here::here("functions.R"))
library(here)

file_paths <- list.files(path = here::here("snakemake/Rdatas/residuals/"), 
                         pattern = "\\.rds", full.names = TRUE)
file_names <-  gsub(pattern = "\\.rds$", replacement = "", x = basename(file_paths))

library(plyr)
library(doMC)
registerDoMC(64)

print("Reading residual files")
data_list <- llply(file_paths, readRDS, .parallel = TRUE)

names(data_list) <- file_names  
print(paste("Included", length(file_names), "studies:", paste(file_names, collapse = ", ")))

# Process row wise metrics
library(dplyr)
library(Rfast)
library(matrixStats)

max_missingness = 0.50
print(paste("Calculating gene level metrics with a max missingness of:", max_missingness))
f =  rowMeans
results_list = data_list
summarized_list = vector("list", length = length(results_list))
names(summarized_list) = names(results_list)
n_dsets = length(results_list)
for (dset_name in names(results_list)) {
    exprDf = results_list[[dset_name]]
    gene_vars = f(exprDf$residuals_noOut)
    summarized_list[[dset_name]] = data.frame(
        Genes = row.names(exprDf$residuals_noOut),
        var = gene_vars) %>% 
        tidyr::separate("Genes", c("Genes", NA))
}
summarized_df = purrr::reduce(summarized_list, full_join, by = "Genes")
colnames(summarized_df)[-1] = names(results_list)
x = summarized_df[1,-1]  
rowmask = apply(summarized_df[,-1], 1, function(x) sum(is.na(x))/n_dsets < max_missingness) 
genemask = summarized_df$Genes[rowmask]
# %%

residuals_noout <- list()
for (name in names(data_list)){
    df <- as.data.frame(data_list[[name]]$residuals_noOut) %>% rotate_df()
    df$source = name
    colnames(df) = remove_id_ver(colnames(df))
    # Median center
    mean = mean(unlist(c(df[,!(names(df) %in% c("source"))])))
    df[,!(names(df) %in% c("source"))] = df[,!(names(df) %in% c("source"))] - mean
    residuals_noout[[name]] <- df
}

# %%

residuals_df <-bind_rows(residuals_noout)
residuals_df_filtered = residuals_df[,(colnames(residuals_df) %in% c(genemask,"source"))]

residuals_df_long <- residuals_df_filtered %>% tibble::rownames_to_column(var = "Row") %>% pivot_longer(cols = starts_with("E"), names_to = "gene",values_to = "value")

# %% 
#%%
ggplot(residuals_df_long, aes(x = value, y = source, fill = stat(quantile))) +
  stat_density_ridges(quantile_lines = TRUE,
                      calc_ecdf = TRUE,
                      geom = "density_ridges_gradient",
                      quantiles = c(0.05, 0.95)) +
  scale_fill_manual(name = "Prob.", values = c("#E2FFF2", "white", "#B0E0E6"),
                    labels = c("(0, 5%]", "(5%, 95%]", "(95%, 1]")) + theme_minimal()
ggsave("ridgeplot_raw_residuals_mean_centered", width = 8, height = 6, units = "in", dpi = 300)
# %%

#%%
metric_df = readRDS(here::here("snakemake/Rdatas/gene_metrics.RDS"))
metric_df$sd[,-1] = metric_df$sd[,-1] %>% mutate_all(scale, scale=FALSE)
metric_df_sd_long = metric_df$sd %>% pivot_longer(cols = -Genes, names_to = "study",values_to = "value")
ggplot(metric_df_sd_long, aes(x = value, y = study, fill = stat(quantile))) +
  stat_density_ridges(quantile_lines = TRUE,
                      calc_ecdf = TRUE,
                      geom = "density_ridges_gradient",
                      quantiles = c(0.025, 0.8)) +
  scale_fill_manual(name = "Prob.", values = c("#E2FFF2", "white", "#B0E0E6"),
                    labels = c("(0, 5%]", "(5%, 95%]", "(95%, 1]")) + theme_minimal()
ggsave("ridgeplot_sd_by_study_mean_centered.png", width = 24, height = 12, units = "in", dpi = 300)
#%%


#%%
metric_df = readRDS(here::here("snakemake/Rdatas/gene_metrics.RDS"))
metric_df$sd[,-1] = metric_df$sd[,-1] %>% mutate_all(scale, scale = FALSE)
metric_df_sd_long = metric_df$sd %>% pivot_longer(cols = -Genes, names_to = "study",values_to = "value")
ggplot(metric_df_sd_long, aes(x = value,group= study,color=study)) + #, fill = stat(quantile))) +
  geom_density()
ggsave("stat_density_sd_by_study_mean_centered.png", width = 24, height = 12, units = "in", dpi = 300)
#%%

# %%
# Fit chisq in loop to find best fit df -- there must be a better way
# https://stackoverflow.com/questions/28922782/r-fitting-chi-squared-distribution-with-large-x-range 
library(MASS)

remove_na_fit_dist = function(x){
    x<-x[!is.na(x)]
    max_df <- 100 # max degree of freedom to test (here from 1 to 100)
    chi_df_disp <- rep(NA,max_df)

    # loop across degree of freedom
    for (i in 1:max_df) {
    chi_adjusted <- (x/mean(x))*i # Adjust the chi-sq distribution so that the mean matches the tested degree of freedom 
    chi_fit <- fitdistr(chi_adjusted,"chi-squared",start=list(df=i),method="BFGS") ## Fitting
    chi_df_disp[i] <- chi_fit$estimate/i # This is going to give you the dispersion between the fitted df and the tested df
    }

    # Find the value with the smallest dispersion (i.e. the best match between the estimated df and the tested df)
    real_df <- which.min(abs(chi_df_disp-1))
    real_df # print the real degree of freedom after correction
}

res = apply(metric_df$sd[,-1], 2, remove_na_fit_dist)
# %%

# write.csv(x= metric_df$sd, here::here("sd_data.csv"))

# %%
n_samples = lapply(residuals_noout, nrow)
cor(res,unlist(n_samples)
# %%


max()