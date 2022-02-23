# Script for analysing the drivers of correlation strucuture

install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()

pak::pkg_install(c("coda","mvtnorm","devtools","loo","dagitty","shape"))
pak::pkg_install("rmcelreath/rethinking")

library(rethinking)
library(cmdstanr)

source(here::here("functions.R"))

metric_df = readRDS(here::here("snakemake/Rdatas/gene_metrics.RDS"))
ea_df = read.csv(here::here("snakemake/metadata/EA_metadata.csv"), header=T)
rc3_df = read.csv(here::here("snakemake/metadata/recount3_metadata.csv"), header=T)
metadata_df = bind_rows(ea_df,rc3_df)
table(metadata_df$tissue)
corr_mat = readRDS(here::here("snakemake/Rdatas/gene_var_matrices.RDS"))$sd

ids = rownames(corr_mat)
(z1 = outer(ids, ids, paste, sep = ":"))
grp = metadata_df[match(ids, metadata_df$id), "group"]
(z2 = outer(grp, grp, paste, sep = ":"))

lt = function(x) x[lower.tri(x)]
sort_label = function(x){
    sapply(x, function(x){
    strs = str_split(x, ":")[[1]]
    paste(str_sort(strs), collapse = ":")
    })
}
corr_df = tibble(corr = 2*atanh(lt(corr_mat)), 
                 pair = sort_label(lt(z1)), 
                 source = sort_label(lt(z2))) %>%
    mutate(pair2 = pair) %>%
    separate(pair2, c("Study1", "Study2"), sep = ":") %>%
    mutate(source = factor(source),
           s1 = match(corr_df$Study1, ids),
           s2 = match(corr_df$Study2, ids),
           n_source = as.numeric(source))


rethinking_data = dplyr::select(corr_df, corr, s1, s2, n_source) %>% as.list()
rethinking_data$index = 1:60
fit_stan <- ulam(
    alist(
        corr ~ normal( mu , sigma ),
        mu <- a + as[s1] + as[s2] + b[n_source],
        as[index] ~ normal(0, 1),
        b[n_source] ~ normal(0, 1),
        a ~ normal( 0 , 1 ),
        sigma ~ exponential( 1 )
    ), data = rethinking_data, chains = 1, cores = 1, iter = 1)
mod <- cmdstan_model( cmdstanr_model_write(rethinking::stancode(fit_stan)) )
fit <- mod$sample(
  data = rethinking_data, 
  seed = 123, 
  chains = 4, 
  parallel_chains = 4,
  adapt_delta = 0.99, max_treedepth = 15
)
fit$summary() %>% as.data.frame()

fit$summary() %>% 
    dplyr::filter(grepl('b', variable)) %>% 
    mutate(source = levels(corr_df$source)) %>% 
    relocate(source)

fit$summary() %>% 
    dplyr::filter(grepl('as', variable)) %>% 
    mutate(id = ids) %>% 
    relocate(id) %>% 
    print(n=60)

levels(corr_df$source)
ids
