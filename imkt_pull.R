# Get PopHuman data
pak::pkg_install("sergihervas/iMKT")
library(iMKT)
loadPopHuman()

# Divergence from chimpanzees -- where available.
imkt_results = NULL
for (gene in unique(PopHumanData$ID)){
  columns_to_ignore = tryCatch(
    expr = {
      results = PopHumanAnalysis_modified(genes=gene, pops=c("CEU"), recomb=FALSE, test="standardMKT", plot=TRUE)
      results_flat = unlist(results$`pop =  CEU`)
      results_flat["gene"] <- gene
      names(results_flat) <- make.names(names(results_flat), unique=TRUE)
      results_df = t(as.data.frame(results_flat))
      rownames(results_df) = gene
      imkt_results <- rbind(imkt_results, results_df)
    },
    error = function(e) {
        print(paste0("Skipped gene: ", gene))
        print(e)
    }
  )
}
write.csv(imkt_results, file="imkt_results.csv")

