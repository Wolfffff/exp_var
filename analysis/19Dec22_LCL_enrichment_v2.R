###########
# ENV RESPONSIVE GENES
###########

setwd('/Users/alea/Dropbox/Amanda_files/ayroles_lab/cell_stimulation_exp/Sep2022_high_var_genes')

# limma
coef=read.delim('31Mar21_DE_coef.txt')
se=read.delim('31Mar21_DE_SE.txt')
pval=read.delim('31Mar21_DE_pval.txt')

# mashR
lfsr=read.delim('31Mar21_LFSR_DE_mashr.txt')
pm=read.delim('31Mar21_pm_DE_mashr.txt')

# assess sharing
lfsr$sig<-apply(lfsr,1,function(x) length(which(x<0.1)))
lfsr$id<-1:dim(lfsr)[1]
genes=read.delim("~/Dropbox/Amanda_files/ayroles_lab/cell_stimulation_exp/Feb2021_NovaSeq_FINAL/31Mar21_all_runs_voom_norm_geneIDs.txt")
lfsr$gene<-genes$V1

lfsr_sig1<-subset(lfsr,sig>0 )

lfsr_sig1$within2_v2<-NA
for (i in 1:dim(lfsr_sig1)[1]){
	z<-min(lfsr_sig1[i,1:11])[1]
	focal<-pm[lfsr_sig1$id[i],which(lfsr_sig1[i,1:11] == z)][1]
	not_focal<-pm[lfsr_sig1$id[i],which(lfsr_sig1[i,1:11] != z)]
	not_focal2<-log2(not_focal/as.numeric(focal))
	lfsr_sig1$within2_v2[i]<- length(which(not_focal2 > -1 & not_focal2<1))
}

lfsr2<-unique(merge(lfsr_sig1[,c('within2_v2','gene')],lfsr,by='gene',all.y=T))
lfsr2$within2_v3<-lfsr2$within2_v2+1
lfsr2$within2_v3[which(is.na(lfsr2$within2_v2))]<-0

pval_FDR<-pval
for (i in 1:11){
pval_FDR[,i]<-qvalue(pval[,i])$qvalues
}

###########
# CROSS TISSUE GENES
# are high vs low variance genes more likely to be responsive to a particular env?
###########

# provides the ranks (colname sd for variance ranks) and genes tested! We've just been taking the top and bottom 5%. Note that low rank implies low variance (i.e. rank 1 is the lowest variance gene).

data=read.csv('pca_ranks.csv')
quant<-quantile(data$sd,seq(0,1,0.05))
data$cat<-'neither'
data$cat[which(data$sd<quant[2])]<-'low'
data$cat[which(data$sd>quant[20])]<-'high'

# low vs others
pval_FDR$gene<-genes$V1
both2=merge(data,pval_FDR,by.x='Gene',by.y='gene')

gene_count<-matrix(ncol=3,nrow=11)
for (i in c(8:14,16:18)){
both2$fdr_cat<-'2_NS'
both2$fdr_cat[(which(both2[,i]<0.1))]<-'1_S'
gene_count[i-7,]<- table(both2$fdr_cat,both2$cat)[1,] }

both2=subset(both2,cat!='high')

output<-matrix(ncol=4,nrow=11)
for (i in c(8:14,16:18)){
both2$fdr_cat<-'2_NS'
both2$fdr_cat[(which(both2[,i]<0.1))]<-'1_S'
output[i-7,1]<-fisher.test(table(both2$fdr_cat,both2$cat))$estimate
output[i-7,2]<-fisher.test(table(both2$fdr_cat,both2$cat))$p.value
output[i-7,3]<-fisher.test(table(both2$fdr_cat,both2$cat))$conf.int[1]
output[i-7,4]<-fisher.test(table(both2$fdr_cat,both2$cat))$conf.int[2]
}

output<-as.data.frame(output)
output$treatment<-names(both2)[8:18]
output$fdr<-p.adjust(output$V2,method='BH')
write.table(output,'low_vs_others.txt',row.names=F,sep='\t')

# high vs others
pval_FDR$gene<-genes$V1
both2=merge(data,pval_FDR,by.x='Gene',by.y='gene')
both2=subset(both2,cat!='low')

output<-matrix(ncol=4,nrow=11)
for (i in c(8:14,16:18)){
both2$fdr_cat<-'2_NS'
both2$fdr_cat[(which(both2[,i]<0.1))]<-'1_S'
output[i-7,1]<-fisher.test(table(both2$fdr_cat,both2$cat))$estimate
output[i-7,2]<-fisher.test(table(both2$fdr_cat,both2$cat))$p.value
output[i-7,3]<-fisher.test(table(both2$fdr_cat,both2$cat))$conf.int[1]
output[i-7,4]<-fisher.test(table(both2$fdr_cat,both2$cat))$conf.int[2]
}

output<-as.data.frame(output)
output$treatment<-names(both2)[8:18]
output$fdr<-p.adjust(output$V2,method='BH')
write.table(output,'high_vs_others.txt',row.names=F,sep='\t')

###########
# CROSS TISSUE GENES
# additional analyses
###########

# high vs others

both=merge(data,lfsr2,by.x='Gene',by.y='gene')
both$env_responsive<-'2_No'
both$env_responsive[which(both$sig>0)]<-'1_Yes'

tmp<-subset(both,env_responsive=='1_Yes')
boxplot(tmp$within2_v3~tmp$cat)
summary(glm(tmp$within2_v3~tmp$cat),family='Poisson')
aggregate(tmp$within2_v3~tmp$cat,FUN=median)

# are high vs low variance genes more likely to be responsive to ANY env? 
both=subset(both,cat!='low')
fisher.test(table(both$env_responsive,both$cat),alternative='greater')

# low vs others

# are high vs low variance genes more likely to be responsive to ANY env? 
both=merge(data,lfsr2,by.x='Gene',by.y='gene')
both$env_responsive<-'2_No'
both$env_responsive[which(both$sig>0)]<-'1_Yes'
both=subset(both,cat!='high')
fisher.test(table(both$env_responsive,both$cat),alternative='greater')

###########
# CROSS TISSUE GENES
# plot
###########

plot=read.delim('SI_table2_lcl_enrichment_AJL.txt')
plot=subset(plot,DE.genes.that.overlap.with.focal.set>4)

library(ggplot2)
ggplot(plot, aes(x=Treatment, y=log2(Odds.Ratio))) +geom_hline(yintercept=0,lty=2)+ geom_point(size=2)+theme_bw(13)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ geom_errorbar(aes(ymin=log2(X95..confidence.interval), ymax=log2(X95..confidence.interval.1)), width=.1, position=position_dodge(0.1))+facet_wrap(~Focal.set,scales = "free_x")+ylab('Log2 odds ratio')+xlab('Environment')
    
###########
# BLOOD GENES
# are high vs low variance genes more likely to be responsive to a particular env?
###########

# provides the ranks (colname sd for variance ranks) and genes tested! We've just been taking the top and bottom 5%. Note that low rank implies low variance (i.e. rank 1 is the lowest variance gene).

data=read.csv('blood.csv')
quant<-quantile(data$sd,seq(0,1,0.05))
data$cat<-'neither'
data$cat[which(data$sd<quant[2])]<-'low'
data$cat[which(data$sd>quant[20])]<-'high'

# low vs others
pval_FDR$gene<-genes$V1
both2=merge(data,pval_FDR,by.x='Gene',by.y='gene')
both2=subset(both2,cat!='high')

output<-matrix(ncol=4,nrow=11)
for (i in c(5:11,13:15)){
both2$fdr_cat<-'2_NS'
both2$fdr_cat[(which(both2[,i]<0.1))]<-'1_S'
output[i-4,1]<-fisher.test(table(both2$fdr_cat,both2$cat))$estimate
output[i-4,2]<-fisher.test(table(both2$fdr_cat,both2$cat))$p.value
output[i-4,3]<-fisher.test(table(both2$fdr_cat,both2$cat))$conf.int[1]
output[i-4,4]<-fisher.test(table(both2$fdr_cat,both2$cat))$conf.int[2]
}

output<-as.data.frame(output)
output$treatment<-names(both2)[5:15]
output$fdr<-p.adjust(output$V2,method='BH')
write.table(output,'low_vs_others_blood.txt',row.names=F,sep='\t')

# high vs others
pval_FDR$gene<-genes$V1
both2=merge(data,pval_FDR,by.x='Gene',by.y='gene')
both2=subset(both2,cat!='low')

output<-matrix(ncol=4,nrow=11)
for (i in c(5:11,13:15)){
both2$fdr_cat<-'2_NS'
both2$fdr_cat[(which(both2[,i]<0.1))]<-'1_S'
output[i-4,1]<-fisher.test(table(both2$fdr_cat,both2$cat))$estimate
output[i-4,2]<-fisher.test(table(both2$fdr_cat,both2$cat))$p.value
output[i-4,3]<-fisher.test(table(both2$fdr_cat,both2$cat))$conf.int[1]
output[i-4,4]<-fisher.test(table(both2$fdr_cat,both2$cat))$conf.int[2]
}

output<-as.data.frame(output)
output$treatment<-names(both2)[5:15]
output$fdr<-p.adjust(output$V2,method='BH')
write.table(output,'high_vs_others_blood.txt',row.names=F,sep='\t')









