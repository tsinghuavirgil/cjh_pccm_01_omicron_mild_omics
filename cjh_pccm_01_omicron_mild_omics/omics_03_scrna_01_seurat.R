library(Seurat)
library(dplyr)
library(stringr)
library(argparse)
library(reshape2)
library(tidyr)
library(ggplot2)
library(scales)
library(cowplot)
library(purrr)
library(harmony)

source('/PROJ/k8s/DATA/nas-9b08bbd0-5dcd-4f3a-bdc8-6f0c22252535/pipeline/seekonetools/src/utils/utils.R')


parser = ArgumentParser()
parser$add_argument("--path", help="/path/to/samples gene bar")
parser$add_argument("--compare", help="samplevssamplevssample,using vs as the split")
parser$add_argument("--species", help="GRCh38 or mm10")
parser$add_argument("--outdir", help="outdir of project")
parser$add_argument("--Nfeatures",help="nfeatures of FindVariableFeatures()", default='2000' )
parser$add_argument("--prefix", help="prefix of results")
parser$add_argument("--resolution", help="resolution for cluster",default='0.6')
parser$add_argument("--dims", help="dims_use for cluster",default='20')
parser$add_argument("--min_cells", help="genes expressed in min cells",default='3')
parser$add_argument("--x_low_cutoff", help="x_low_cutoff for search high variable genes",default='0.0125')
parser$add_argument("--x_high_cutoff", help="x_high_cutoff for search high variable genes",default='3')
parser$add_argument("--y_cutoff", help="y_cutoff for search high variable genes",default='0.5')
args <- parser$parse_args()


str(args)

path=args$path
compare=args$compare
####gene_path=args$gene_path
species = args$species
outdir=args$outdir
prefix=args$prefix
resolution=args$resolution
dims=args$dims
min_cells=args$min_cells
x_low_cutoff=args$x_low_cutoff
x_high_cutoff=args$x_high_cutoff
y_cutoff=args$y_cutoff
Nfeatures = args$Nfeatures

resolution <- as.numeric(resolution)
min_cells <- as.numeric(min_cells)
Nfeatures<-as.numeric(Nfeatures)
dims<- as.numeric(dims)
x_low_cutoff <- as.numeric(x_low_cutoff)
x_high_cutoff <- as.numeric(x_high_cutoff)
y_cutoff <- as.numeric(y_cutoff)

if (is.null(prefix)){
	prefix=compare
}

ob.list <- list()
samples<-strsplit(compare,'vs')[[1]]

numsap=1

for (each in samples){
  	pbmc <- readRDS(paste0(path,'/',each,'/',each,'_seurat.rds'))
        if(grepl('-1',colnames(pbmc@assays$RNA@counts)[1])){
	colnames(pbmc@assays$RNA@counts) <- str_replace_all(colnames(pbmc@assays$RNA@counts), '-1',paste0('-',numsap))
        }else{
        colnames(pbmc@assays$RNA@counts) <- paste0(colnames(pbmc@assays$RNA@counts),'-',numsap)
        }
	ob <- CreateSeuratObject(counts =pbmc@assays$RNA@counts,project =each,min.cells = min_cells)
	ob$stim <-each
#	ob <- NormalizeData(ob)
#	ob <- FindVariableFeatures(ob,  selection.method = "vst",nfeatures = Nfeatures)
	numsap=numsap+1
	ob.list[[each]] <- ob
}


obj.integrated<-merge(ob.list[[1]],ob.list[[2]])

for (i in (3:length(ob.list))){
        obj.integrated<-merge(obj.integrated,ob.list[[i]])
}


obj.integrated <- NormalizeData(obj.integrated, normalization.method = "LogNormalize", scale.factor = 10000)
obj.integrated <- FindVariableFeatures(obj.integrated, selection.method = "vst", nfeatures = 2000)

obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
obj.integrated <- RunPCA(obj.integrated, npcs = dims, verbose = FALSE)

obj.integrated <- obj.integrated %>% RunHarmony("stim", plot_convergence = TRUE)

obj.integrated <- RunUMAP(obj.integrated, reduction = "harmony", dims = 1:dims)
obj.integrated <- RunTSNE(obj.integrated, reduction = "harmony", dims = 1:dims)
obj.integrated <- FindNeighbors(obj.integrated,reduction = "harmony", dims = 1:dims)
obj.integrated <- FindClusters(obj.integrated, resolution = resolution)


##umap
umap<-obj.integrated@reductions$umap@cell.embeddings
umap<-cbind(rownames(umap),umap)
colnames(umap)[1]<-'Barcode'

tsne<-obj.integrated@reductions$tsne@cell.embeddings
tsne<-cbind(rownames(tsne),tsne)
colnames(tsne)[1]<-'Barcode'

###cluster
clus<-as.data.frame(Idents(obj.integrated))
clus<-cbind(rownames(clus),clus)
colnames(clus)<-c('Barcode','Cluster')


sample_number <- length(samples)
num <- sample_number + 1

b <- mutate(clus,Sample=as.numeric(str_split(clus$Barcode,'-',simplify = TRUE)[,2]))

c <- b %>%
  group_by(Cluster) %>%
  count(Sample) %>%
  spread(key=Sample,value=n)
colnames(c)[2:num] <- samples

#calculate the relative percent of each sample in each cluster
df <- c[,-1]
rowsum <- rowSums(df,na.rm=T)
df_r <- df/rowsum

#print out the relative table with heads
df_r <- cbind(as.numeric(rownames(df_r))-1, df_r)
colnames(df_r)[1] <- "Cluster"
colnames(df_r)[2:num] <- samples



colour1 <-MYCOLOR[1:sample_number]
##colour1 <- hue_pal()(sample_number)

td <- gather(df_r,key="Cluster Name",value="Cells Ratio",-Cluster)
td[,1] <- factor(td[,1], levels = sort(as.numeric(df_r$Cluster)))
td[,2] <- factor(td[,2], levels = samples)


plt<- ggplot(td,aes(x=td[,1],y=td[,3],fill=td[,2]))+
  geom_bar(position = 'stack',stat="identity")+
  labs(x="Cluster Name",y="Cells Ratio")+
  theme(panel.background=element_rect(fill='transparent', color='black'),
        legend.key=element_rect(fill='transparent', color='transparent'),axis.text = element_text(color="black"))+
  scale_y_continuous(expand=c(0.001,0.001))+
  scale_fill_manual(values=colour1)+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1,ncol=1,title = 'Sample'))

#clusters in sample percent
df_c <- t(t(df)/rowSums(t(df),na.rm=T))
row.names(df_c) <-c(1:dim(df_c)[1])
df_c <- as.data.frame(cbind((as.numeric(rownames(df_c))-1), df_c))
colnames(df_c)[1] <- "Cluster"

cluster_number <- length(row.names(df_c))

#colour2 <- hue_pal()(cluster_number)
colour2 <-MYCOLOR[1:cluster_number]


td_c <- gather(df_c,key="Sample Name",value="Cells Ratio",-Cluster)
td_c[,1] <- factor(td_c[,1], levels = sort(as.numeric(df_c$Cluster)))
td_c[,2] <- factor(td_c[,2], levels = samples)

plt_c<- ggplot(td_c,aes(x=td_c[,2],y=td_c[,3],fill=td_c[,1]))+
  geom_bar(position = 'stack',stat="identity")+
  labs(x="Sample Name",y="Cells Ratio")+
  theme(panel.background=element_rect(fill='transparent', color='black'),
        legend.key=element_rect(fill='transparent', color='transparent'),axis.text = element_text(color="black"))+
  scale_y_continuous(expand=c(0.001,0.001))+
  scale_fill_manual(values=colour2)+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1,ncol=1,title = 'Cluster'))


####Find All Markers
obj.markers <- FindAllMarkers(obj.integrated, only.pos = F, min.pct = 0.25,logfc.threshold=0.25)

obj.markers <- obj.markers[,c(7,6,2,1,5,3,4)]

cluster <- unique(obj.markers$cluster)
for (i in cluster) {
  data <- filter(obj.markers,cluster == i)
  data1 <- filter(obj.markers,cluster == i, p_val_adj < 0.05, avg_log2FC > 0)
}

top9 <- obj.markers %>% group_by(cluster) %>% top_n(n = 9, wt = avg_log2FC)
obj.integrated <- ScaleData(obj.integrated, features=top9$gene, assay='RNA')

p2 <- refineDoHeatmap(obj.integrated, features = unique(top9$gene), label=FALSE, assay='RNA')


l0 <- group_split(top9)
names(l0) <- group_keys(top9)[['cluster']]
for (x in names(l0)){
    tmp0 <- arrangeTop9(
                FeaturePlot(obj.integrated, features=l0[[x]][['gene']], reduction='umap', pt.size=0.8, combine=FALSE)
            )
    tmp1 <- arrangeTop9(
                VlnPlot(obj.integrated, features=l0[[x]][['gene']], pt.size=0, cols=MYCOLOR[1:length(levels(Idents(obj.integrated)))], combine=FALSE),
                legend.position='none')
}










