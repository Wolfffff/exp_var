## run the wilcoxon rank sum test (between the genes in the highest and lowest 5% of the variance/mean rank metrics) to associate ranks with the chromHMM annotations ##


#################################
## arguments to pass to script ##

args = commandArgs( trailingOnly = TRUE)

## args[ 1] is distance around gene used in kb
## args[ 2] is the crosse-tissue variance and mean rank metric file
## args[ 3] is the directory where all the tissue-level variance and mean rank metrics files are


##############################
## set up workspace (local) ##

library(tidyverse)
library(ggthemes)
library(gridExtra)
library(RColorBrewer)
library(ggrepel)

load( "~/Documents/AyrolesLab/varGeneExp/Figures/pcor_rankVchromHMM10_forHeatmapPlot_FigS1.RData")


####################################
## run the Wilcoxon rank sum test ##

wilcoxRes.var <- NULL
wilcoxRes.mean <- NULL

thisVarrank <- read.csv( "~/Documents/AyrolesLab/varGeneExp/annotationData/pca_ranks.csv", header = T) %>% dplyr::select( Gene, mean, sd) ## cross-tissue ranks
colnames( thisVarrank)[ 3] <- c( "varrank")

thisVarrank <- thisVarrank %>% mutate( varrankBin = ntile( thisVarrank$varrank, n = 20))
thisVarrank <- thisVarrank %>% mutate( meanrankBin = ntile( thisVarrank$mean, n = 20))

thisChromHMM <- get( "chromHMMcov.genome")

finalDF <- merge( thisVarrank, thisChromHMM, by = 1)

for ( k in 6:ncol( finalDF)) {
  
  thisCat <- colnames( finalDF)[ k] ## annotation
  
  tmpDF.tissue <- finalDF[ , c( 4, k)] ## varrank
  tmpDF.tissue <- tmpDF.tissue[ complete.cases( tmpDF.tissue), ]
  
  medianL <- median( tmpDF.tissue[ which( tmpDF.tissue$varrankBin == 1), 2]) ## lowest 5%
  medianH <- median( tmpDF.tissue[ which( tmpDF.tissue$varrankBin == 20), 2]) ## highest 5%
  SEML <- sd( tmpDF.tissue[ which( tmpDF.tissue$varrankBin == 1), 2]) / sqrt( length( which( tmpDF.tissue$varrankBin == 1))) ## lowest 5%
  SEMH <- sd( tmpDF.tissue[ which( tmpDF.tissue$varrankBin == 20), 2]) / sqrt( length( which( tmpDF.tissue$varrankBin == 20))) ## highest 5%
  SDL <- sd( tmpDF.tissue[ which( tmpDF.tissue$varrankBin == 1), 2]) ## lowest 5%
  SDH <- sd( tmpDF.tissue[ which( tmpDF.tissue$varrankBin == 20), 2]) ## highest 5%
  
  wilcoxRes.var <- rbind( wilcoxRes.var, data.frame( category = thisCat, 
                                                     medianL = medianL, 
                                                     medianH = medianH, 
                                                     SEML = SEML,
                                                     SEMH = SEMH,
                                                     SDL = SDL,
                                                     SDH = SDH,
                                                     medianDiffHsubtractL = medianH - medianL, 
                                                     wilcoxP = wilcox.test( tmpDF.tissue[ which( tmpDF.tissue$varrankBin == 20), 2], 
                                                                            tmpDF.tissue[ which( tmpDF.tissue$varrankBin == 1), 2])$p.value
                                                     )
                          )
  
  tmpDF.tissue <- finalDF[ , c( 5, k)] ## meanrank
  tmpDF.tissue <- tmpDF.tissue[ complete.cases( tmpDF.tissue), ]
  
  medianL <- median( tmpDF.tissue[ which( tmpDF.tissue$meanrankBin == 1), 2]) ## lowest 5%
  medianH <- median( tmpDF.tissue[ which( tmpDF.tissue$meanrankBin == 20), 2]) ## highest 5%
  SEML <- sd( tmpDF.tissue[ which( tmpDF.tissue$meanrankBin == 1), 2]) / sqrt( length( which( tmpDF.tissue$meanrankBin == 1))) ## lowest 5%
  SEMH <- sd( tmpDF.tissue[ which( tmpDF.tissue$meanrankBin == 20), 2]) / sqrt( length( which( tmpDF.tissue$meanrankBin == 20))) ## highest 5%
  SDL <- sd( tmpDF.tissue[ which( tmpDF.tissue$meanrankBin == 1), 2]) ## lowest 5%
  SDH <- sd( tmpDF.tissue[ which( tmpDF.tissue$meanrankBin == 20), 2]) ## highest 5%
  
  wilcoxRes.mean <- rbind( wilcoxRes.mean, data.frame( category = thisCat, 
                                                     medianL = medianL, 
                                                     medianH = medianH, 
                                                     SEML = SEML,
                                                     SEMH = SEMH,
                                                     SDL = SDL,
                                                     SDH = SDH,
                                                     medianDiffHsubtractL = medianH - medianL, 
                                                     wilcoxP = wilcox.test( tmpDF.tissue[ which( tmpDF.tissue$meanrankBin == 20), 2], 
                                                                            tmpDF.tissue[ which( tmpDF.tissue$meanrankBin == 1), 2])$p.value
                                                     )
                          )
  
  
}

wilcoxRes.var <- wilcoxRes.var %>% mutate( padj = p.adjust( wilcoxP, method = "BH")) %>% mutate( sig = NA)
wilcoxRes.var$sig[ which( wilcoxRes.var$padj < 0.05)] <- c( "yes")
wilcoxRes.var$sig[ which( wilcoxRes.var$padj >= 0.05)] <- c( "no")

wilcoxRes.mean <- wilcoxRes.mean %>% mutate( padj = p.adjust( wilcoxP, method = "BH")) %>% mutate( sig = NA)
wilcoxRes.mean$sig[ which( wilcoxRes.mean$padj < 0.05)] <- c( "yes")
wilcoxRes.mean$sig[ which( wilcoxRes.mean$padj >= 0.05)] <- c( "no")


####################
## plot lineplots ##

## variance

# newDF <- data.frame( cbind( wilcoxRes.var[ , c( 1, 2, 4, 9)], type = "low")) %>% rename( median = medianL, SEM = SEML)
# newDF <- rbind( newDF, data.frame( cbind( wilcoxRes.var[ , c( 1, 3, 5, 9)], type = "high")) %>% rename( median = medianH, SEM = SEMH))
newDF <- data.frame( cbind( wilcoxRes.var[ , c( 1, 2, 6, 11)], type = "low")) %>% rename( median = medianL, SD = SDL)
newDF <- rbind( newDF, data.frame( cbind( wilcoxRes.var[ , c( 1, 3, 7, 11)], type = "high")) %>% rename( median = medianH, SD = SDH))
newDF$type <- factor( newDF$type, levels = c( "low", "high"))
newDF$type <- factor( newDF$type, labels = c( "Bottom 5% of Genes\nRanked by Variance", "Top 5% of Genes\nRanked by Variance"))

# get color palette and make non-significant comparisons black
colors <- brewer.pal( 10, "Paired")
colors <- colorRampPalette( colors)( 13)

tmpDF.col <- data.frame( color = colors, category = levels( as.factor( newDF$category)))

tmpDF.col.var <- unique( merge( tmpDF.col, newDF[ , c( "category", "sig")] ))
tmpDF.col.var$color <- as.character( tmpDF.col.var$color)
tmpDF.col.var$color[ which( tmpDF.col.var$sig == "no")] <- c( "#000000")

library(ggthemes)
svg( file = paste0( "top_bottom_5pVarrank_fxnlGen_lineplot_10kb.svg"), width = 12.5, height = 8.5)
ggplot( newDF, aes( x = type, y = median, group = category, color = category)) +
  geom_line( size = 2, alpha = 0.6, position = position_dodge( width = 0.25)) +
  geom_point( size = 3, alpha = 0.6, position = position_dodge( width = 0.25)) +
  geom_pointrange( aes( ymin = median - SD, ymax = median + SD), position = position_dodge( width = 0.25)) +
  geom_label_repel(size = 6,label = newDF$category, position = position_dodge( width = 0.25), family = "serif",max.overlaps = Inf) +
  ylab( "Median proportion of gene region in indicated annotation") +
  xlab("") + 
  scale_color_manual( values = tmpDF.col.var$color) + 
  theme_tufte( base_size = 22) +
  theme(axis.line = element_line("black"), legend.position = "none",
  axis.title = element_text(size = 22),axis.text=element_text(size=22,color="black"))
  # theme( panel.grid.major = element_line( color = "grey90"), legend.position = "none")
dev.off()

gg1 <- ggplot( newDF, aes( x = type, y = median, group = category, color = category)) + geom_line( size = 2, alpha = 0.6, position = position_dodge( width = 0.25)) + geom_point( size = 3, alpha = 0.6, position = position_dodge( width = 0.25)) + geom_pointrange( aes( ymin = median - SD, ymax = median + SD), position = position_dodge( width = 0.25)) + geom_label_repel( label = newDF$category, position = position_dodge( width = 0.25), family = "serif") + theme_bw() + theme( legend.position = "none") + ylab( "Median proportion of gene region in indicated annotation") + xlab( "") + scale_color_manual( values = tmpDF.col.var$color) + theme_tufte( base_size = 20) + theme( panel.grid.major = element_line( color = "grey90"), legend.position = "none")

## mean rank

# newDF <- data.frame( cbind( wilcoxRes.mean[ , c( 1, 2, 4, 9)], type = "low")) %>% rename( median = medianL, SEM = SEML)
# newDF <- rbind( newDF, data.frame( cbind( wilcoxRes.mean[ , c( 1, 3, 5, 9)], type = "high")) %>% rename( median = medianH, SEM = SEMH))
newDF <- data.frame( cbind( wilcoxRes.mean[ , c( 1, 2, 6, 11)], type = "low")) %>% rename( median = medianL, SD = SDL)
newDF <- rbind( newDF, data.frame( cbind( wilcoxRes.mean[ , c( 1, 3, 7, 11)], type = "high")) %>% rename( median = medianH, SD = SDH))
newDF$type <- factor( newDF$type, levels = c( "low", "high"))
newDF$type <- factor( newDF$type, labels = c( "Bottom 5%", "Top 5%"))

# get color palette and make non-significant comparisons black
colors <- brewer.pal( 10, "Paired")
colors <- colorRampPalette( colors)( 13)

tmpDF.col <- data.frame( color = colors, category = levels( as.factor( newDF$category)))

tmpDF.col.var <- unique( merge( tmpDF.col, newDF[ , c( "category", "sig")] ))
tmpDF.col.var$color <- as.character( tmpDF.col.var$color)
tmpDF.col.var$color[ which( tmpDF.col.var$sig == "no")] <- c( "#000000")


pdf( file = paste0( "top_bottom_5pMeanrank_fxnlGen_lineplot_10kb.pdf"), width = 8, height = 9)
ggplot( newDF, aes( x = type, y = median, group = category, color = category)) + geom_line( size = 2, alpha = 0.6, position = position_dodge( width = 0.25)) + geom_point( size = 3, alpha = 0.6, position = position_dodge( width = 0.25)) + geom_pointrange( aes( ymin = median - SD, ymax = median + SD), position = position_dodge( width = 0.25)) + geom_label_repel( label = newDF$category, position = position_dodge( width = 0.25), family = "serif") + theme_bw() + theme( legend.position = "none") + ylab( "Median proportion of gene region in indicated annotation") + xlab( "") + scale_color_manual( values = tmpDF.col.var$color) + theme_tufte( base_size = 20) + theme( panel.grid.major = element_line( color = "grey90"), legend.position = "none")
dev.off()


gg2 <- ggplot( newDF, aes( x = type, y = median, group = category, color = category)) + geom_line( size = 2, alpha = 0.6, position = position_dodge( width = 0.25)) + geom_point( size = 3, alpha = 0.6, position = position_dodge( width = 0.25)) + geom_pointrange( aes( ymin = median - SD, ymax = median + SD), position = position_dodge( width = 0.25)) + geom_label_repel( label = newDF$category, position = position_dodge( width = 0.25), family = "serif") + theme_bw() + theme( legend.position = "none") + ylab( "Median proportion of gene region in indicated annotation") + xlab( "") + scale_color_manual( values = tmpDF.col.var$color) + theme_tufte( base_size = 20) + theme( panel.grid.major = element_line( color = "grey90"), legend.position = "none")

pdf( file = paste0( "top_bottom_5pVarMeanrank_fxnlGen_lineplots_10kb.pdf"), width = 16, height = 9)
grid.arrange( gg1, gg2, ncol = 2)
dev.off()

