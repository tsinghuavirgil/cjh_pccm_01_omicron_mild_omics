rm(list=ls())
setwd('/xxx/')
getwd()
options(scipen=1000000)
#######
library('dplyr')
library('ggplot2')
library('gplots')
library(BiocParallel)
library(RColorBrewer)
suppressPackageStartupMessages(library(vioplot))
suppressPackageStartupMessages(library(mclust))

scatter_colbrewer <- c('#36A0D0', '#E56B59', '#EF9D25', '#9AC3C4', '#636393'
)#
plot(1: 5,rep(1, 5), col = scatter_colbrewer, pch = 16, cex = 2)
brewer_pal<-colorRampPalette(scatter_colbrewer)
###########
fig_violin <- as.data.frame(read.csv('score_ifn_stimulated_genes.csv', header=T, row.names= 1))
head(fig_violin)
class(fig_violin)
score_color = length(unique(fig_violin$time))
###########
score_fig <- ggplot(fig_violin,aes(x=time,
                                   y=log2(score_ifn_stimulated_genes + 1)))+
  theme_classic()
score_fig <- score_fig + geom_violin(aes(fill = time),trim = FALSE)
scale_fill_manual(values = brewer_pal(score_color)) +
  geom_boxplot(width=.05)
print(score_fig)
ggsave(score_fig, file='violin_ifn_stimulated_genes.pdf', width = 6, height = 4)
