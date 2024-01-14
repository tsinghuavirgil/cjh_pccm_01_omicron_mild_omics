suppressMessages({
library(scRepertoire)
library(Seurat)
library(ggsci)
library(tidyverse)
library(ggraph)
library(circlize)
library(scales)
library(venn)
library(GSVA)
library(ComplexHeatmap)
source('utils.R')})

obj <- readRDS("Tcell_combined.rds")
obj$patient <- ifelse(obj$stim%in%c("DYM_1_5","DYM_2_5","DYM_3_5","DYM_4_5","DYM_5_5"), "DYM",
                         ifelse(obj$stim%in%c("JQW_1_5","JQW_2_5","JQW_3_5","JQW_4_5","JQW_5_5"), "JQW",
                         ifelse(obj$stim%in%c("MHL_1_5","MHL_2_5","MHL_3_5","MHL_4_5","MHL_5_5"), "MHL",
                         ifelse(obj$stim%in%c("YMW_1_5","YMW_2_5","YMW_3_5","YMW_4_5","YMW_5_5"), "YMW",
                         ifelse(obj$stim%in%c("SCC_1_5","SCC_2_5","SCC_3_5","SCC_4_5","SCC_5_5"), "SCC", "LB"
                            )))))
obj$stage <- ifelse(obj$stim%in%c("DYM_1_5","JQW_1_5","MHL_1_5","YMW_1_5","SCC_1_5"), "before",
                         ifelse(obj$stim%in%c("DYM_2_5","JQW_2_5","MHL_2_5","YMW_2_5","SCC_2_5"), "day1",
                         ifelse(obj$stim%in%c("DYM_3_5","JQW_3_5","MHL_3_5","YMW_3_5","SCC_3_5","LB_3_5"), "day3",
                         ifelse(obj$stim%in%c("DYM_4_5","JQW_4_5","MHL_4_5","YMW_4_5","SCC_4_5","LB_4_5"), "day5", "day7"
                            ))))
obj@meta.data$seurat_clusters <- factor(obj@meta.data$seurat_clusters, levels = sort(as.numeric(as.character(unique(obj$seurat_clusters)))))

anno <- RenameIdents(obj,
`0` = "Effect CD8T",
`1` = "Th1|Th17 CD4T",
`2` = "Effect CD8T",
`3` = "Naive CD4T",
`4` = "Naive CD4T",
`5` = "Memory CD4T",
`6` = "Naive CD8T",
`7` = "Naive CD4T",
`8` = "MAIT",
`9` = "γδT",
`10` = "Effect Memory CD8T",
`11` = "NK",
`12` = "TREG CD4T",
`13` = "Th2",
`14` = "NKT")
obj$celltype <- Idents(anno)

### 合并多个样本的filtered_contig_annotations.csv
sample="DYM_1_TCRvsDYM_2_TCRvsDYM_3_TCRvsDYM_4_TCRvsDYM_5_TCRvsJQW_1_TCRvsJQW_2_TCRvsJQW_3_TCRvsJQW_4_TCRvsJQW_5_TCRvsMHL_1_TCRvsMHL_2_TCRvsMHL_3_TCRvsMHL_4_TCRvsMHL_5_TCRvsYMW_1_TCRvsYMW_2_TCRvsYMW_3_TCRvsYMW_4_TCRvsYMW_5_TCRvsSCC_1_TCRvsSCC_2_TCRvsSCC_3_TCRvsSCC_4_TCRvsSCC_5_TCRvsLB_3_TCRvsLB_4_TCRvsLB_5_TCR"
newSample="DYM_1_5vsDYM_2_5vsDYM_3_5vsDYM_4_5vsDYM_5_5vsJQW_1_5vsJQW_2_5vsJQW_3_5vsJQW_4_5vsJQW_5_5vsMHL_1_5vsMHL_2_5vsMHL_3_5vsMHL_4_5vsMHL_5_5vsYMW_1_5vsYMW_2_5vsYMW_3_5vsYMW_4_5vsYMW_5_5vsSCC_1_5vsSCC_2_5vsSCC_3_5vsSCC_4_5vsSCC_5_5vsLB_3_5vsLB_4_5vsLB_5_5"
group="DYMvsDYMvsDYMvsDYMvsDYMvsJQWvsJQWvsJQWvsJQWvsJQWvsMHLvsMHLvsMHLvsMHLvsMHLvsYMWvsYMWvsYMWvsYMWvsYMWvsSCCvsSCCvsSCCvsSCCvsSCCvsLBvsLBvsLB"
path="./data/TCR/st/"

#########################################1.数据处理
samples <- strsplit(sample,'vs')[[1]]
newSample <- strsplit(newSample,'vs')[[1]]
group <- strsplit(group, 'vs')[[1]]

contig_list <- list()
for (each in samples){
contig_obj <- read.csv(paste0(path,'/', each,"/filtered_contig_annotations.csv"))
contig_list[[each]] <- contig_obj
    }

numsap = 1
for(i in seq_along(contig_list)){
    
if('-1'%in%contig_list[[i]]$barcode[1]){
    contig_list[[i]]$barcode <- str_replace_all(contig_list[[i]]$barcode, '-1',paste0('-',numsap))
}else{
    contig_list[[i]]$barcode <- paste0(contig_list[[i]]$barcode,'-',numsap)
}
numsap=numsap+1
}

combined <- combineTCR(contig_list,
                samples = newSample,
                ID = group,
                cells = "T-AB",
                removeNA = FALSE,
                removeMulti = FALSE,
                filterMulti = FALSE)

for(i in seq_along(combined)){
combined[[i]] <- stripBarcode(combined[[i]] , column = 1, connector = "_", num_connects = 5)

#count <- as.data.frame(t(data.frame(strsplit(combined[[i]][, 1], 
#        paste("['", "_", "']", sep = "")), stringsAsFactors = FALSE)), 
#        stringsAsFactors = FALSE)[5:6]
#combined[[i]][, 1] <- paste(count[, 1], count[, 2], sep = "_")
   
}

seurat <- combineExpression(combined, obj,
                  cloneCall="strict",
                  proportion = FALSE,
                  cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

##TCR和转录组overlap
all_combined<- do.call("rbind", combined)
vdj_b <- all_combined$barcode
sc_b <- colnames(obj)
options(repr.plot.height=7, repr.plot.width=7)
venn(
  list(TCR = vdj_b, SC = sc_b), 
  zcolor = "style",      # add color
  opacity = 0.25, 
  ellipse = F)

seurat@meta.data$TCR <- ifelse(rownames(seurat@meta.data) %in% rownames(na.omit(seurat@meta.data)),"TCR","NA")
DimPlot(seurat, group.by = "TCR") +
    scale_color_manual(values = c("grey","red")) + 
  theme(plot.title = element_blank())

## Clonal vs Non-Clonal
seurat@meta.data$Clonal <- ifelse(seurat@meta.data$cloneType=="Single (0 < X <= 1)","Non-Clonal","Clonal")
DimPlot(seurat, group.by = "Clonal") +
    scale_color_manual(values = c("red","grey")) + 
  theme(plot.title = element_blank())

## clonal水平稳态
slot(seurat, "meta.data")$cloneType <- factor(slot(seurat, "meta.data")$cloneType, 
                levels = c("Hyperexpanded (100 < X <= 500)", 
                           "Large (20 < X <= 100)", 
                            "Medium (5 < X <= 20)", 
                            "Small (1 < X <= 5)", 
                            "Single (0 < X <= 1)", NA))
options(repr.plot.height=7, repr.plot.width=9)
DimPlot(seurat, group.by = "cloneType") +
    scale_color_manual(values = pal_npg("nrc")(9), na.value="grey") + 
  theme(plot.title = element_blank())

occupiedscRepertoire(seurat, 
                     x.axis = "celltype",
                     label=F)+
theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1)) + scale_fill_npg()

occupiedscRepertoire(seurat, 
                     x.axis = "celltype",
                     label=F,
                     proportion = T)+
theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1)) + scale_fill_npg()

occupiedscRepertoire(seurat, 
                     x.axis = "stim",
                     label=F,facet.by = "patient")+
theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1)) + scale_fill_npg()+
facet_wrap(~patient, ncol = 6, scales = "free_x")

occupiedscRepertoire(seurat, 
                     x.axis = "stim",
                     label=F,
                     proportion = T,facet.by = "patient")+
theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1)) + scale_fill_npg()+
facet_wrap(~patient, ncol = 6, scales = "free_x")

occupiedscRepertoire(seurat, 
                     x.axis = "stim",
                     label=F,facet.by = "stage")+
theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1)) + scale_fill_npg()+
facet_wrap(~stage, ncol = 6, scales = "free_x")

occupiedscRepertoire(seurat, 
                     x.axis = "stim",
                     label=F,
                     proportion = T,facet.by = "stage")+
theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1)) + scale_fill_npg()+
facet_wrap(~stage, ncol = 6, scales = "free_x")

### 每种细胞类型的克隆水平
for (each in unique(seurat$celltype)){
p1 <- occupiedscRepertoire(subset(seurat, subset=celltype==each), 
                     x.axis = "stim",
                     label=F,facet.by = "patient")+
theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1)) + scale_fill_npg()+
facet_wrap(~patient, ncol = 6, scales = "free_x")+ ggtitle(each)

p2 <- occupiedscRepertoire(subset(seurat, subset=celltype==each), 
                     x.axis = "stim",
                     label=F,
                     proportion = T,facet.by = "patient")+
theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1)) + scale_fill_npg()+
facet_wrap(~patient, ncol = 6, scales = "free_x")
}

## clonotype共享
compareClonotypes(seurat, numbers = 5, 
                  cloneCall="strict", graph = "alluvial", split.by = "stim") +
            theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))+  NoLegend()


## 克隆多样性
clonalDiversity(combined, 
                cloneCall = "strict", 
                group.by = "sample", 
                x.axis = "ID", 
                n.boots = 100)+ 
        scale_color_manual(values = MYCOLOR)
