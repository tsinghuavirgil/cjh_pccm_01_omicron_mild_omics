rm(list = ls())
setwd('/xxx/')
getwd()
library(reshape2)
library(ggplot2)
#library(limma)
library(pheatmap)
library(Hmisc)
##############
data<-read.csv('unique_01.txt',header=T, row.names = 1, sep = '\t', check.names=F)#
data[1]
colnames(data)
rownames(data)
class(data)

################
#annotation_col<-read.csv('heatmap_02_nametime.csv',header=T,row.names=1)
#annotation_col <- as.data.frame(annotation_col)
#annotation_col = annotation_col
pheatmap(log2(data + 1), scale='row', cluster_rows=T, cluster_cols=F,
         #annotation_col = annotation_col,
         cellwidth = 5, cellheight = 5, border_color='white',
         color = colorRampPalette(c('#7f9eb2','white', '#881600'))(100),
         fontsize_row = 6, fontsize_col = 6)
