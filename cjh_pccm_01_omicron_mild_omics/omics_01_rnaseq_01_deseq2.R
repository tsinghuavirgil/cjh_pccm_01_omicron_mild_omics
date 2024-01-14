rm(list = ls())
######
setwd('/xxx/')
getwd()
library(BiocManager)
library(colorspace)
library(dplyr)
library(plyr)
library(DESeq2)
library(BiocParallel)
library(RColorBrewer)
library(gplots)
library(amap)
library(gmodels)
library(ggpubr)
library(ggplot2)
library(ggthemes)
library(limma)
library(wesanderson)
library(corrplot)
library(EnhancedVolcano)
########################
count_data = read.table('All_Sample_Count.txt', header=T, row.names='Geneid',
                        com='', quote='', check.names=F, sep='\t')

countdata <- as.matrix(count_data)
head(countdata)
countdata <- countdata[rowSums(countdata)>46,]
countdata
dim(countdata)

meta_data <- read.table('meta_data_02.txt', header=T, row.names=1, com='',
                        quote='', check.names=F, sep='\t', colClasses='factor')
meta_data
all(rownames(meta_data) %in% colnames(countdata))
countdata <- countdata[, rownames(meta_data)]
all(rownames(meta_data) == colnames(countdata))

#########
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = countdata,
                                            colData = meta_data,
                                            design= ~ name + time)
dds <- DESeq(ddsFullCountTable)  

normalized_counts <- counts(dds, normalized=TRUE)
head(normalized_counts)

normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
rld <- rlog(dds, blind=FALSE)
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pearson_cor <- as.matrix(cor(rlogMat, method="pearson"))
hc <- hcluster(t(rlogMat), method="pearson")
fig_rnaseq_01 <- plot(hc,  cex = .5)
print(fig_rnaseq_01)
dev.new()
heatmap.2(pearson_cor, Rowv=as.dendrogram(hc), symm=T, trace="none",
          col=hmcol, margins=c(11,11), main="The pearson correlation of each sample")
dev.off()
#########
plotPCA(rld, 'time')
plotPCA(rld, 'name')
rld_mat <- assay(rld)
pca_rld <- prcomp(rld_mat)
pca_data <- data.frame(pca_rld$rotation)
merge_pca <- cbind(pca_data, meta_data)
#############################################################
fig_pca <- ggscatter(merge_pca, x = 'PC1', y = 'PC2',
                     color = 'time', #shape = 'name_c',
                     ellipse = TRUE, ellipse.type = "convex",
                     palette = c('#30A9DE', '#f26d5b', '#f9a11b', '#9DC8C8', '#666699'),
) + theme_base()
ggsave(fig_pca_id, file='rnaseq_03_pca_06_sample.pdf', width = 7, height = 4)
#########################################################################################

sample_day0 = 'day0'
sample_day1 = 'day1'
sample_day3 = 'day3'
sample_day5 = 'day5'
sample_day7 = 'day7'
###################
deg_day0 <- counts(dds, normalized=TRUE)[, colData(dds)$time == sample_day0]

if (is.vector(deg_day0)){
  basemean_day0 <- as.data.frame(deg_day0)
} else {
  basemean_day0 <- as.data.frame(rowMeans(deg_day0))
}
colnames(basemean_day0) <- sample_day0
head(basemean_day0)

deg_day1 <- counts(dds, normalized=TRUE)[, colData(dds)$time == sample_day1]

if (is.vector(deg_day1)){
  basemean_day1 <- as.data.frame(deg_day1)
} else {
  basemean_day1 <- as.data.frame(rowMeans(deg_day1))
}
colnames(basemean_day1) <- sample_day1
head(basemean_day1)

deg_day3 <- counts(dds, normalized=TRUE)[, colData(dds)$time == sample_day3]

if (is.vector(deg_day3)){
  basemean_day3 <- as.data.frame(deg_day3)
} else {
  basemean_day3 <- as.data.frame(rowMeans(deg_day3))
}
colnames(basemean_day3) <- sample_day3
head(basemean_day3)

deg_day5 <- counts(dds, normalized=TRUE)[, colData(dds)$time == sample_day5]

if (is.vector(deg_day5)){
  basemean_day5 <- as.data.frame(deg_day5)
} else {
  basemean_day5 <- as.data.frame(rowMeans(deg_day5))
}
colnames(basemean_day5) <- sample_day5
head(basemean_day5)

deg_day7 <- counts(dds, normalized=TRUE)[, colData(dds)$time == sample_day7]

if (is.vector(deg_day7)){
  basemean_day7 <- as.data.frame(deg_day7)
} else {
  basemean_day7 <- as.data.frame(rowMeans(deg_day7))
}
colnames(basemean_day7) <- sample_day7
head(basemean_day7)

day1_vs_day0 <- c('time', sample_day1, sample_day0)
res_day1_vs_day0 <- results(dds,  contrast = day1_vs_day0)
summary(res_day1_vs_day0)
sum(res_day1_vs_day0$padj < 0.05, na.rm=TRUE)
res_day1_vs_day0_2 <- cbind(basemean_day0, basemean_day1, as.data.frame(res_day1_vs_day0))
head(res_day1_vs_day0_2)
summary(res_day1_vs_day0_2)

res_day1_vs_day0_2$padj[is.na(res_day1_vs_day0_2$padj)] <- 1
res_day1_vs_day0_2 <- res_day1_vs_day0_2[order(res_day1_vs_day0_2$padj),]
head(res_day1_vs_day0_2)

write.csv(res_day1_vs_day0_2, file='volcano_01_day1vs0.csv', quote=F)

sig_day1_vs_day0 <- res_day1_vs_day0_2[which(res_day1_vs_day0_2$padj<0.05
                                             & abs(res_day1_vs_day0_2$log2FoldChange)>1),]
sig_day1_vs_day0[which(sig_day1_vs_day0$log2FoldChange>0),'up_down']<-'up'
sig_day1_vs_day0[which(sig_day1_vs_day0$log2FoldChange<0),'up_down']<-'down'
head(sig_day1_vs_day0)
dim(sig_day1_vs_day0)
write.csv(sig_day1_vs_day0, file='deseq2_01_day1_vs_day0_02_degs_new.csv', quote=F)

day3_vs_day0 <- c('time', sample_day3, sample_day0)
res_day3_vs_day0 <- results(dds,  contrast = day3_vs_day0)
summary(res_day3_vs_day0)

res_day3_vs_day0_2 <- cbind(basemean_day0, basemean_day3, as.data.frame(res_day3_vs_day0))
head(res_day3_vs_day0_2)
summary(res_day3_vs_day0_2)

res_day3_vs_day0_2$padj[is.na(res_day3_vs_day0_2$padj)] <- 1

res_day3_vs_day0_2 <- res_day3_vs_day0_2[order(res_day3_vs_day0_2$padj),]
head(res_day3_vs_day0_2)

write.csv(res_day3_vs_day0_2, file='volcano_02_day3vs0.csv', quote=F)

sig_day3_vs_day0 <- res_day3_vs_day0_2[which(res_day3_vs_day0_2$padj<0.05
                                             & abs(res_day3_vs_day0_2$log2FoldChange)>1),]
sig_day3_vs_day0[which(sig_day3_vs_day0$log2FoldChange>0),'up_down']<-'up'
sig_day3_vs_day0[which(sig_day3_vs_day0$log2FoldChange<0),'up_down']<-'down'
head(sig_day3_vs_day0)
write.csv(sig_day3_vs_day0, file='02_day3_vs_day0_02_degs_new.csv', quote=F)

day5_vs_day0 <- c('time', sample_day5, sample_day0)
res_day5_vs_day0 <- results(dds,  contrast = day5_vs_day0)
summary(res_day5_vs_day0)

res_day5_vs_day0_2 <- cbind(basemean_day0, basemean_day5, as.data.frame(res_day5_vs_day0))
head(res_day5_vs_day0_2)
summary(res_day5_vs_day0_2)

res_day5_vs_day0_2$padj[is.na(res_day5_vs_day0_2$padj)] <- 1

res_day5_vs_day0_2 <- res_day5_vs_day0_2[order(res_day5_vs_day0_2$padj),]
head(res_day5_vs_day0_2)

write.csv(res_day5_vs_day0_2, file='volcano_03_day5vs0.csv', quote=F)

sig_day5_vs_day0 <- res_day5_vs_day0_2[which(res_day5_vs_day0_2$padj<0.05
                                             & abs(res_day5_vs_day0_2$log2FoldChange)>1),]
sig_day5_vs_day0[which(sig_day5_vs_day0$log2FoldChange>0),'up_down']<-'up'
sig_day5_vs_day0[which(sig_day5_vs_day0$log2FoldChange<0),'up_down']<-'down'
head(sig_day5_vs_day0)
write.csv(sig_day5_vs_day0, file='03_day5_vs_day0_02_degs_new.csv', quote=F)

day7_vs_day0 <- c('time', sample_day7, sample_day0)
res_day7_vs_day0 <- results(dds,  contrast = day7_vs_day0)
summary(res_day7_vs_day0)

res_day7_vs_day0_2 <- cbind(basemean_day0, basemean_day7, as.data.frame(res_day7_vs_day0))
head(res_day7_vs_day0_2)
summary(res_day7_vs_day0_2)

res_day7_vs_day0_2$padj[is.na(res_day7_vs_day0_2$padj)] <- 1

res_day7_vs_day0_2 <- res_day7_vs_day0_2[order(res_day7_vs_day0_2$padj),]
head(res_day7_vs_day0_2)

write.csv(res_day7_vs_day0_2, file='volcano_04_day7vs0.csv', quote=F)

sig_day7_vs_day0 <- res_day7_vs_day0_2[which(res_day7_vs_day0_2$padj<0.05
                                             & abs(res_day7_vs_day0_2$log2FoldChange)>1),]
sig_day7_vs_day0[which(sig_day7_vs_day0$log2FoldChange>0),'up_down']<-'up'
sig_day7_vs_day0[which(sig_day7_vs_day0$log2FoldChange<0),'up_down']<-'down'
head(sig_day7_vs_day0)
write.csv(sig_day7_vs_day0, file='04_day7_vs_day0_02_degs_new.csv', quote=F)


day3_vs_day1 <- c('time', sample_day3, sample_day1)
res_day3_vs_day1 <- results(dds,  contrast = day3_vs_day1)
summary(res_day3_vs_day1)

res_day3_vs_day1_2 <- cbind(basemean_day1, basemean_day3, as.data.frame(res_day3_vs_day1))
head(res_day3_vs_day1_2)
summary(res_day3_vs_day1_2)

res_day3_vs_day1_2$padj[is.na(res_day3_vs_day1_2$padj)] <- 1

res_day3_vs_day1_2 <- res_day3_vs_day1_2[order(res_day3_vs_day1_2$padj),]
dim(res_day3_vs_day1_2)

write.csv(res_day3_vs_day1_2, file='deseq2_05_day3_vs_day1_01_res.csv', quote=F)

sig_day3_vs_day1 <- res_day3_vs_day1_2[which(res_day3_vs_day1_2$padj<0.05
                                             & abs(res_day3_vs_day1_2$log2FoldChange)>1),]
sig_day3_vs_day1[which(sig_day3_vs_day1$log2FoldChange>0),'up_down']<-'up'
sig_day3_vs_day1[which(sig_day3_vs_day1$log2FoldChange<0),'up_down']<-'down'
dim(sig_day3_vs_day1)
write.csv(sig_day3_vs_day1, file='deseq2_05_day3_vs_day1_02_degs_new.csv', quote=F)


day5_vs_day1 <- c('time', sample_day5, sample_day1)
res_day5_vs_day1 <- results(dds,  contrast = day5_vs_day1)
summary(res_day5_vs_day1)

res_day5_vs_day1_2 <- cbind(basemean_day1, basemean_day5, as.data.frame(res_day5_vs_day1))
head(res_day5_vs_day1_2)
summary(res_day5_vs_day1_2)

res_day5_vs_day1_2$padj[is.na(res_day5_vs_day1_2$padj)] <- 1

res_day5_vs_day1_2 <- res_day5_vs_day1_2[order(res_day5_vs_day1_2$padj),]
head(res_day5_vs_day1_2)

write.csv(res_day5_vs_day1_2, file='deseq2_06_day5_vs_day1_01_res.csv', quote=F)

sig_day5_vs_day1 <- res_day5_vs_day1_2[which(res_day5_vs_day1_2$padj<0.05
                                             & abs(res_day5_vs_day1_2$log2FoldChange)>1),]
sig_day5_vs_day1[which(sig_day5_vs_day1$log2FoldChange>0),'up_down']<-'up'
sig_day5_vs_day1[which(sig_day5_vs_day1$log2FoldChange<0),'up_down']<-'down'
head(sig_day5_vs_day1)
write.csv(sig_day5_vs_day1, file='deseq2_06_day5_vs_day1_02_degs.csv', quote=F)


day5_vs_day3 <- c('time', sample_day5, sample_day3)
res_day5_vs_day3 <- results(dds,  contrast = day5_vs_day3)
summary(res_day5_vs_day3)

res_day5_vs_day3_2 <- cbind(basemean_day3, basemean_day5, as.data.frame(res_day5_vs_day3))
head(res_day5_vs_day3_2)
summary(res_day5_vs_day3_2)

res_day5_vs_day3_2$padj[is.na(res_day5_vs_day3_2$padj)] <- 1

res_day5_vs_day3_2 <- res_day5_vs_day3_2[order(res_day5_vs_day3_2$padj),]
head(res_day5_vs_day3_2)

write.csv(res_day5_vs_day3_2, file='deseq2_07_day5_vs_day3_01_res.csv', quote=F)

sig_day5_vs_day3 <- res_day5_vs_day3_2[which(res_day5_vs_day3_2$padj<0.05
                                             & abs(res_day5_vs_day3_2$log2FoldChange)>1),]
sig_day5_vs_day3[which(sig_day5_vs_day3$log2FoldChange>0),'up_down']<-'up'
sig_day5_vs_day3[which(sig_day5_vs_day3$log2FoldChange<0),'up_down']<-'down'
head(sig_day5_vs_day3)
write.csv(sig_day5_vs_day3, file='deseq2_07_day5_vs_day3_02_degs.csv', quote=F)

