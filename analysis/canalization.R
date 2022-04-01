A = c(1, 1, 1, 1)
B = rep(0, 4)
P = c(0.5, 0.5, 0.5, 0.5) + c(0, 0, 1, 1)

Norm <- function(x) sqrt(sum(x^2))

getDist = function(P, A, B = rep(0, length(A))){
    pa = P - A
    ba = B - A
    t = (pa %*% ba) / (ba %*% ba)
    d = Norm(pa - c(t) * ba)
    d
}
getDist(P, A)

corrMat = readRDS(here::here("snakemake/Rdatas/gene_var_matrices.RDS"))$sd
PCscores = readRDS(here::here("snakemake/Rdatas/gene_rank_pca_scores.RDS"))
ranks = readRDS(here::here("snakemake/Rdatas/gene_var_rank.RDS"))

str(corrMat)
str(PCscores)
PCscores = scale(PCscores, scale = FALSE)
PC1 = -eigen(corrMat)$vectors[,1]
ISO_dist = apply(PCscores[,1:10], 1, getDist, rep(1, 60))
summary(ISO_dist)
png("test.png")
plot(ranks$sd, ISO_dist)
dev.off()

png("test.png")
plot(abs(PCscores[,2]), ISO_dist)
dev.off()

png("test.png")
plot(PCscores[,1:2])
dev.off()

cor(PCscores[,1:2], method = "spearman")
str(PCscores)
