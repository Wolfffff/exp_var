
# %%
library(here)
library(sjmisc)
library(pryr)
source(here::here("functions.R"))
library(here)

file_paths <- list.files(path = here::here("snakemake/Rdatas/residuals/"),
                         pattern = "\\.rds",
                         full.names = TRUE)
file_names <- gsub(pattern = "\\.rds$", replacement = "", x = basename(file_paths))

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

max_missingness <- 0.5
print(paste("Calculating gene level metrics with a max missingness of:", max_missingness))
f <- rowMeans
results_list <- data_list
summarized_list <- vector("list", length = length(results_list))
names(summarized_list) <- names(results_list)
n_dsets <- length(results_list)
for (dset_name in names(results_list)) {
  exprDf <- results_list[[dset_name]]
  gene_vars <- f(exprDf$residuals_noOut)
  summarized_list[[dset_name]] <- data.frame(Genes = row.names(exprDf$residuals_noOut),
                                             var = gene_vars) %>%
    tidyr::separate("Genes", c("Genes", NA))
}
summarized_df <- purrr::reduce(summarized_list, full_join, by = "Genes")
colnames(summarized_df)[-1] <- names(results_list)
x <- summarized_df[1, -1]
rowmask <- apply(summarized_df[, -1], 1, function(x) sum(is.na(x))/n_dsets < max_missingness)
# %%









# %%
source(here::here("functions.R"))
library(here)
library(sjmisc)
library(pryr)

# %%

# %% Read in the data for each chunk in parallel if needed.
file_paths <- list.files(path = here::here("snakemake/Rdatas/residuals/"),
                         pattern = "\\.rds",
                         full.names = TRUE)
file_names <- gsub(pattern = "\\.rds$", replacement = "", x = basename(file_paths))

library(plyr)
library(doMC)
registerDoMC(64)
data_list <- llply(file_paths, readRDS, .parallel = TRUE)

names(data_list) <- file_names

r

# %% %%


residuals_df <- bind_rows(residuals_noout)
dim(residuals_df)
n_dsets <- length(data_list)
# residuals_df <- residuals_df %>% rotate_df() residuals_df$genes = rownames(residuals_df)
# summarized_df <- residuals_df # summarized_df = purrr::reduce(summarized_list, dplyr::full_join,
# by = 'Genes') # colnames(summarized_df)[-1] = names(results_list) # x = summarized_df[1,-(1:2)]
# max_missingness =0.5 rowmask = apply(summarized_df[,-c(1:2)], 1, function(x)
# sum(is.na(x))/n_dsets < max_missingness) summarized_df = summarized_df[rowmask,] %%

residuals_df_long <- residuals_df %>%
  tibble::rownames_to_column(var = "Row") %>%
  pivot_longer(cols = starts_with("E"), names_to = "gene", values_to = "value")

# %%


# %%
residuals_df_filtered <- residuals_df[, which(colnames(df) %in% summarized_df$Genes[rowmask])]
residuals_df_long <- residuals_df_filtered %>%
  tibble::rownames_to_column(var = "Row") %>%
  pivot_longer(cols = starts_with("E"), names_to = "gene", values_to = "value")

# %% %%
library(ggridges)
library(ggplot2)
ggplot(residuals_df_long,
       aes(x = value, y = source, fill = source)) + geom_density_ridges() + theme_ridges() +
  theme(legend.position = "none")
ggsave("test.png")
# %% %%
