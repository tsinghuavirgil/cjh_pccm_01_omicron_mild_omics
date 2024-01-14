rm(list=ls())
setwd('/xxx/')
getwd()
################
library('dplyr')
library('ggplot2')
library('gplots')
suppressPackageStartupMessages(library(mclust))

###########
#read data#
###########
countdata <- as.matrix(read.table('state_01_nametime_02_all_05_tpm_sort.txt', header=T, row.names= 1))
coldata <- read.table('state_02_meta.txt', header=T, row.names= 1)
all(rownames(coldata) %in% colnames(countdata))
countdata <- countdata[, rownames(coldata)]
all(rownames(coldata) == colnames(countdata))
class(countdata)
head(countdata)
###############
bulkseq_set <- read.table('state_03_gene_sets.txt', sep = '\t', header = TRUE)
head(bulkseq_set)
dim(bulkseq_set)
table(bulkseq_set)
reads_single_phase = countdata
#######################

#interferon_stimulated_genes
ifn_stimulated_genes = bulkseq_set$gene_symbol[which(bulkseq_set$state %in% c('interferon_stimulated_genes'))]
ifn_stimulated_genes_single_phase = as.matrix(reads_single_phase[rownames(reads_single_phase) %in% (ifn_stimulated_genes) ,])
ifn_stimulated_genes_combined_matrix = rbind(ifn_stimulated_genes_single_phase,average=apply(ifn_stimulated_genes_single_phase,2,mean))
ifn_stimulated_genes_cor_matrix = cor(t(ifn_stimulated_genes_combined_matrix))
ifn_stimulated_genes_cor_vector = ifn_stimulated_genes_cor_matrix[,dim(ifn_stimulated_genes_cor_matrix)[1]]
ifn_stimulated_genes_single_phase_restricted = ifn_stimulated_genes_single_phase[rownames(ifn_stimulated_genes_single_phase) %in% names(ifn_stimulated_genes_cor_vector[ifn_stimulated_genes_cor_vector >= 0.1]),]
ifn_stimulated_genes_score = apply(ifn_stimulated_genes_single_phase_restricted,2,mean)
score_mtx = cbind(coldata,ifn_stimulated_genes_score)
write.csv(score_mtx, 'score_ifn_stimulated_genes.csv')
