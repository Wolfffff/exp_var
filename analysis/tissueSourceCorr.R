# Script for analysing the drivers of correlation strucuture

# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# cmdstanr::install_cmdstan()
# pak::pkg_install(c("coda", "mvtnorm", "devtools", "loo", "dagitty", "shape", 
#                    "bayesplot", "patchwork", "sjmisc"))
# pak::pkg_install("rmcelreath/rethinking")

source(here::here("functions.R"))


library(rethinking)
library(cmdstanr)
library(bayesplot)
library(patchwork)


metric_df = readRDS(here::here("snakemake/Rdatas/gene_metrics.RDS"))
ea_df = read.csv(here::here("snakemake/metadata/EA_metadata.csv"),
                 header=T, comment.char = "#")
rc3_df = read.csv(here::here("snakemake/metadata/recount3_metadata.csv"),
                  header=T, comment.char = "#")
metadata_df = bind_rows(ea_df, rc3_df) %>%
    mutate(group = gsub("Other - Expression Atlas", "Other", group),
           group = gsub("Other - recount3", "Other", group))

file_paths <- list.files(path = here::here("snakemake/Rdatas/residuals/"), 
                         pattern = "\\.rds", full.names = TRUE)
file_names <-  gsub(pattern = "\\.rds$", replacement = "", x = basename(file_paths))

library(plyr)
library(doMC)
registerDoMC(64)

print("Reading residual files")
data_list <- llply(file_paths, readRDS, .parallel = TRUE)
names(data_list) <- file_names

corr_mat = readRDS(here::here("snakemake/Rdatas/gene_var_matrices.RDS"))$sd

ids = rownames(corr_mat)
(z1 = outer(ids, ids, paste, sep = ":"))
grp = metadata_df[match(ids, metadata_df$id), "group"]
(z2 = outer(grp, grp, paste, sep = ":"))
tissue = metadata_df[match(ids, metadata_df$id), "tissue"]
(z3 = outer(tissue, tissue, paste, sep = ":"))
tissue = metadata_df[match(ids, metadata_df$id), "tissue"]
(z4 = outer(tissue, tissue, `==`))
n_samples = sapply(data_list, function(x) ncol(x$residuals_noOut))[ids]
(z5 = outer(n_samples, n_samples, function(x, y) sqrt(x * y)))

lt = function(x) x[lower.tri(x)]
sort_label = function(x){
    sapply(x, function(x){
    strs = str_split(x, ":")[[1]]
    paste(str_sort(strs), collapse = ":")
    })
}
corr_df = tibble(corr = atanh(lt(corr_mat)),
                 pair = sort_label(lt(z1)),
                 source = sort_label(lt(z2)),
                 tissue = sort_label(lt(z3)),
                 bool_tissue = lt(z4),
                 sample_size = scale(lt(z5))[,1]) %>%
    mutate(pair2 = pair) %>%
    separate(pair2, c("Study1", "Study2"), sep = ":") %>%
    mutate(source = factor(source),
           s1 = match(Study1, ids),
           s2 = match(Study2, ids),
           n_source = as.numeric(source),
           n_bool_tissue = as.numeric(bool_tissue) + 1)


quantile(lt(corr_mat), c(0.025, 0.975))
quantile(tanh(mean(corr_df$corr) + rnorm(100000, 0, 0.25)), c(0.025, 0.975))

save_plot("test.png", ggplot(corr_df, aes(sample_size, corr)) + geom_point())
rethinking_data = dplyr::select(corr_df,
                                corr, s1, s2,
                                n_source, n_bool_tissue,
                                sample_size) %>%
                                as.list()
rethinking_data$index = 1:60
fit_stan <- ulam(
    alist(
        corr ~ normal(mu, sigma),
        mu <- a + as[s1] + as[s2] + b[n_source] + c[n_bool_tissue] + d*sample_size,
        as[index] ~ normal(0, 0.25),
        b[n_source] ~ normal(0, 0.25),
        c[n_bool_tissue] ~ normal(0, 0.25),
        d ~ normal(0, 0.5),
        a ~ normal(0, 1),
        sigma ~ exponential(1)
    ), data = rethinking_data, chains = 1, cores = 1, iter = 2)
mod <- cmdstan_model(cmdstanr_model_write(rethinking::stancode(fit_stan)))
fit <- mod$sample(
  data = rethinking_data,
  chains = 8,
  parallel_chains = 8,
  iter_warmup = 2000,
  iter_sampling = 2000,
  adapt_delta = 0.99, max_treedepth = 14
)
fit$summary() %>% as.data.frame()

fit$summary() %>%
    dplyr::filter(grepl('b', variable)) %>%
    mutate(source = levels(corr_df$source)) %>%
    relocate(source)

fit$summary() %>%
    dplyr::filter(grepl('c|d', variable))

fit$summary() %>%
    dplyr::filter(grepl('as', variable)) %>%
    mutate(id = ids) %>%
    relocate(id) %>%
    print(n = 60)

p_tissue = mcmc_intervals(fit$draws("c")) +
    scale_y_discrete(labels = c("Different tissue", "Same tissue")) +
    ggtitle("B. Tissue congruence effect")

p_source = mcmc_intervals(fit$draws("b")) +
    scale_y_discrete(labels = levels(corr_df$source)) +
    ggtitle("C. Pairwise source effect")

p_study = mcmc_intervals(fit$draws("as")) +
    scale_y_discrete(labels = paste0(paste(ids, grp, sep = " ("), ")")) +
    ggtitle("A. Pairwise random effect - Study (Source)")

p_model = p_study + (p_tissue / p_source)  +
  plot_annotation(title = "Modeling the driver of across-study variance correlations",
    caption = "Linear effect model coefficients with Fisher z-transformed spearman correlations as the response. \n A: The pairwise random effect captures the non-independence of the correlation values and estimates the contribution \n of each study to the correlation. For example: comparisons involving bone_marrow tend to be lower than the others. \n B and C: Fixed effect estimates: correlations among studies that use the same tissue are higher, and correlations \ninvolving studies in the \"Other\" category (non gTEX and TCGA) tend to be lower.")
save_plot(here::here("data/plots/correlationModeling.png"),
          p_model, base_height = 6, base_asp = 1, ncol = 2, nrow = 2)


