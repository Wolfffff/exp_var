
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
# %%

# %%
library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

metric_df = readRDS(here::here("snakemake/Rdatas/gene_metrics.RDS"))
metric_df_copy = metric_df

metric_df$sd[,-1] = metric_df$sd[,-1] %>% mutate_all(scale, scale = FALSE)
# quantile normalize
# metric_df_copy$sd[,-1] = as.data.frame(normalize.quantiles(as.matrix(metric_df_copy$sd[,-1])))
metric_df_copy$sd[,-1] = metric_df_copy$sd[,-1] %>% mutate_all(scale, scale = TRUE)

metric_df_sd_long = metric_df$sd %>% pivot_longer(cols = -Genes, names_to = "study",values_to = "value")
metric_df_sd_long_scaled = metric_df_copy$sd %>% pivot_longer(cols = -Genes, names_to = "study",values_to = "value")

rank_df = read.csv(here::here("data/pca_ranks.csv"), header = TRUE)[, -1]
high_gene = rank_df[which(rank_df$sd >= sort(rank_df$sd, decreasing=TRUE)[1]),]$Gene
low_gene = rank_df[which(rank_df$sd <= sort(rank_df$sd)[1]),]$Gene

print(low_gene)
print(high_gene)

high_df = metric_df_sd_long[metric_df_sd_long$Genes %in% high_gene,]
low_df = metric_df_sd_long[metric_df_sd_long$Genes %in% low_gene,]

high_df_scaled = metric_df_sd_long_scaled[metric_df_sd_long_scaled$Genes %in% high_gene,]
low_df_scaled = metric_df_sd_long_scaled[metric_df_sd_long_scaled$Genes %in% low_gene,]

upper_col = wes_palette("Royal1")[1]
lower_col = wes_palette("Royal1")[2]
unscaled <- ggplot(metric_df_sd_long, aes(x = value,group= study,color=study)) + scale_color_manual(values=col_vector) +#, fill = stat(quantile))) +
  geom_density() + geom_rug(data = low_df,aes(x=value), col = lower_col) + geom_rug(data = high_df,aes(x=value), color = upper_col) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position = "none") +
  xlab("Standard deviation") + ylab("Density") #+ ggtitle("Density Plot Standard Deviation by Study")
# ggsave("stat_density_sd_by_study_mean_centered_rug.png", width = 24, height = 12, units = "in", dpi = 300)


scaled <- ggplot(metric_df_sd_long_scaled, aes(x = value,group= study,color=study)) +
  scale_color_manual(values=col_vector) +
  geom_density() +
  geom_rug(data = low_df_scaled,aes(x=value), col = lower_col) +
  geom_rug(data = high_df_scaled,aes(x=value), color = upper_col) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position = "none") +
  xlab("Z-normalized standard deviation") + ylab("Density")# + ggtitle("Density Plot of Z-score Normalized Standard Deviation by Study")
#%%

# %%
library(wesanderson)
library(patchwork)
scaled = scaled #+ xlim(-3,10)
unscaled = unscaled + xlim(-1,2)
density_plot = scaled + inset_element(unscaled, 0.33, 0.33, 1, 1)
density_plot
ggsave(here::here("data/plots/sd_dist.png"), width = 12, height = 6, units = "in", dpi = 300)

saveRDS(list(scaled = scaled, unscaled = unscaled), here::here("snakemake/Rdatas/plots/density_plot.RDS"))
# %%


# %%
metric_df_sd_long_scaled
names(metric_df_sd_long_scaled)[3] <- "scaled_sd"
rank_df = read.csv(here::here("data/pca_ranks.csv"), header = TRUE)[, -1]
rank_df = rename(rank_df, Genes = Gene)
scaled_sd_ranked = inner_join(metric_df_sd_long_scaled, dplyr::select(rank_df, Genes, sd), by= "Genes")
library(ggthemes)
p = ggplot(scaled_sd_ranked, aes(sd, scaled_sd, color = study, group = study)) + geom_smooth() + theme_tufte() + theme(legend.position = "none")
save_plot("example.png", p, base_height = 7, base_asp = 2)
# %%
