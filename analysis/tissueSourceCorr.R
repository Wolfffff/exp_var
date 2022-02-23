# Script for analysing the drivers of correlation strucuture

install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()

pak::pkg_install(c("coda","mvtnorm","devtools","loo","dagitty","shape"))
pak::pkg_install("rmcelreath/rethinking")

library(rethinking)
library(cmdstanr)

source(here::here("functions.R"))

metric_df = readRDS(here::here("snakemake/Rdatas/gene_metrics.RDS"))
ea_df = read.csv(here::here("snakemake/metadata/EA_metadata.csv"), header=T, comment.char = "#")
rc3_df = read.csv(here::here("snakemake/metadata/recount3_metadata.csv"), header=T, comment.char = "#")
metadata_df = bind_rows(ea_df,rc3_df)
corr_mat = readRDS(here::here("snakemake/Rdatas/gene_var_matrices.RDS"))$sd

ids = rownames(corr_mat)
(z1 = outer(ids, ids, paste, sep = ":"))
grp = metadata_df[match(ids, metadata_df$id), "group"]
(z2 = outer(grp, grp, paste, sep = ":"))
tissue = metadata_df[match(ids, metadata_df$id), "tissue"]
(z3 = outer(tissue, tissue, paste, sep = ":"))
tissue = metadata_df[match(ids, metadata_df$id), "tissue"]
(z4 = outer(tissue, tissue, `==`))

lt = function(x) x[lower.tri(x)]
sort_label = function(x){
    sapply(x, function(x){
    strs = str_split(x, ":")[[1]]
    paste(str_sort(strs), collapse = ":")
    })
}
corr_df = tibble(corr = 2*atanh(lt(corr_mat)), 
                 pair = sort_label(lt(z1)), 
                 source = sort_label(lt(z2)),
                 tissue = sort_label(lt(z3)),
                 bool_tissue = lt(z4)) %>%
    mutate(pair2 = pair) %>%
    separate(pair2, c("Study1", "Study2"), sep = ":") %>%
    mutate(source = factor(source),
           s1 = match(corr_df$Study1, ids),
           s2 = match(corr_df$Study2, ids),
           n_source = as.numeric(source), 
           n_bool_tissue = as.numeric(bool_tissue)+1)


rethinking_data = dplyr::select(corr_df, corr, s1, s2, n_source, n_bool_tissue) %>% as.list()
rethinking_data$index = 1:60
fit_stan <- ulam(
    alist(
        corr ~ normal( mu , sigma ),
        mu <- a + as[s1] + as[s2] + b[n_source] + c[n_bool_tissue],
        as[index] ~ normal(0, .3),
        b[n_source] ~ normal(0, .3),
        c[n_bool_tissue] ~ normal(0, .3),
        a ~ normal( 0 , 1 ),
        sigma ~ exponential( 1 )
    ), data = rethinking_data, chains = 1, cores = 1, iter = 2)
mod <- cmdstan_model( cmdstanr_model_write(rethinking::stancode(fit_stan)) )
fit <- mod$sample(
  data = rethinking_data, 
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
    dplyr::filter(grepl('c', variable))

fit$summary() %>% 
    dplyr::filter(grepl('as', variable)) %>% 
    mutate(id = ids) %>% 
    relocate(id) %>% 
    print(n=60)

library(bayesplot)
mcmc_hist(fit$draws("mu"))


p_tissue = mcmc_intervals(fit$draws("c")) + 
    scale_y_discrete(labels = c("Different tissue", "Same tissue"))

p_source = mcmc_intervals(fit$draws("b")) + 
    scale_y_discrete(labels = levels(corr_df$source))

p_study = mcmc_intervals(fit$draws("as"))


library(patchwork)
p_model = p_study + (p_tissue/p_source)
save_plot("test.png", p_model, base_height = 7, base_asp = 1, ncol = 2, nrow = 2)
