source(here::here("functions.R"))

metric_df = readRDS(here::here("snakemake/Rdatas/gene_metrics.RDS"))
ea_df = read.csv(here::here("snakemake/metadata/EA_metadata.csv"), header=T)
rc3_df = read.csv(here::here("snakemake/metadata/recount3_metadata.csv"), header=T)
metadata_df = bind_rows(ea_df,rc3_df)
rank_df = read.csv(here::here("data/pca_ranks.csv"), header = TRUE)[, -1]
library(dplyr)

quantile_violin_plot = function(x,y){
    df = data.frame(x=x,y=y)
    df = df %>% mutate(quantilegroup = ntile(x, 10))
    ggplot(df, aes(x=factor(quantilegroup),group=quantilegroup,y)) + geom_violin(fill = wes_palette("Royal2")[5]) +
    scale_x_discrete(labels=c(1:10))
}
library(wesanderson)

quantile_violin_plot(rank_df$sd,rank_df$mean) + ylab("SD Rank") + xlab("metric quantile") + geom_boxplot(width=0.1)# + theme_minimal() + theme(legend.position = "none")
ggsave("example.jpg")