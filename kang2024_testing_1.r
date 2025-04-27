# setwd('')
dir.create('results')
options(stringsAsFactors = F,check.bounds = F)
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(ggsci)

#
# dir_name <- list.dirs("GSE193337_RAW/", full.names = FALSE, recursive = FALSE)
# dir_name
# datalist <- list()
# for (i in seq_along(dir_name)){
#   dir_10x <- paste0("GSE193337_RAW/", dir_name[i])
#   my_data <- Read10X(data.dir = dir_10x)
#   colnames(my_data) <- paste0(dir_name[i], colnames(my_data))
#   datalist[[i]] <- CreateSeuratObject(counts = my_data, project = dir_name[i],
#                                 min.cells = 3, min.features = 250)
#   datalist[[i]]$Samples <- dir_name[i]
#   datalist[[i]]$type <- substr(dir_name[i], 1, 1)
# }
# names(datalist) <- dir_name
# for (i in seq_along(datalist)){
#     sce <- datalist[[i]]
#     sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
#     sce[["percent.Ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")
#     datalist[[i]] <- sce
#     rm(sce)
# }
# sce <- merge(datalist[[1]], y=datalist[2:length(datalist)])
# raw_cell=sce@meta.data
# raw_count <- table(raw_cell$Samples)
# raw_count
# sum(raw_count) # 12554

# pearplot_befor<-VlnPlot(sce, group.by ='Samples', features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"), pt.size = 0, ncol = 4)
# pearplot_befor
# # ggsave('results/pearplot_befor.pdf',pearplot_befor,height = 5,width = 15)
# ggsave('results/pearplot_befor.jpg',pearplot_befor,height = 5,width = 15,dpi = 300)

# sample_color<-pal_nejm(alpha = 0.5)(8)[1:8]
# sample_color
# Feature_ber1<-FeatureScatter(sce,feature1 = 'nFeature_RNA', feature2 = 'nCount_RNA', group.by = 'Samples', cols = sample_color)
# Feature_ber2<-FeatureScatter(sce,feature1 = 'percent.mt', feature2 = 'nCount_RNA', group.by = 'Samples', cols = sample_color)
# Feature_ber3<-FeatureScatter(sce,feature1 = 'percent.mt', feature2 = 'nFeature_RNA', group.by = 'Samples', cols = sample_color)
# Feature_ber1=Feature_ber1+theme(legend.position = 'none')
# Feature_ber2=Feature_ber2+theme(legend.position = 'none')
# Feature_ber<-ggarrange(Feature_ber1,Feature_ber2,Feature_ber3,ncol = 3,nrow = 1,widths =c(1,1,1.2))
# # ggsave('results/Feature_cor.pdf',Feature_ber,height = 5,width = 17)
# ggsave('results/Feature_cor.jpg',Feature_ber,height = 5,width = 17,dpi = 300)

# datalist <- lapply(X = datalist, FUN = function(x) {
# x<-subset(x,subset = nFeature_RNA < 5000 &
# percent.mt < 15)
# })
# sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
# clean_cell=sce@meta.data
# clean_count <- table(clean_cell$Samples)
# clean_count
# sum(clean_count) # 6410
# pearplot_after <- VlnPlot(sce,group.by ='Samples', features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"), pt.size = 0, ncol = 4)
# pearplot_after
# # ggsave('results/pearplot_after.pdf',pearplot_after,height = 5,width = 15)
# ggsave('results/pearplot_after.jpg',pearplot_after,height = 5,width = 15,dpi = 300)
# save(datalist, file = 'datalist.RData')

load('datalist.RData')

sce <- merge(datalist[[1]], y=datalist[2:length(datalist)])
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000, mean.cutoff = c(0.0125,3), dispersion.cutoff = c(1.5,Inf))
sce <- ScaleData(sce, features = rownames(sce))
sce <- RunPCA(sce, features = VariableFeatures(sce))
elbowplot <- ElbowPlot(sce, ndims=50, reduction="pca")
# elbowplot
ggsave('results/elbowplot.pdf',elbowplot,height = 5,width = 5)

Dims <- 30
sce <- RunUMAP(sce, dims = 1:Dims, reduction = "pca")
raw.umap<-DimPlot(sce,group.by='Samples', reduction="umap", label = "T", pt.size = 0.2, label.size = 0) + ggtitle('')
# raw.umap
ggsave('results/raw.umap.pdf', raw.umap, height = 7, width = 7)
library(clustree)
sce <- FindNeighbors(sce, dims = 1:Dims)
sce <- FindClusters(
object = sce, resolution = c(seq(.1, 1, .1))
)
colnames(sce@meta.data)
clustree(sce@meta.data, prefix = "RNA_snn_res.")
pdf('results/clust.snn_res.pdf',he=15,wi=15)
clustree(sce@meta.data, prefix = "RNA_snn_res.")
dev.off()
Resolution <- 0.8
sce <- FindNeighbors(object = sce, dims = 1:Dims)
sce <- FindClusters(object = sce, resolution = Resolution)
DefaultAssay(sce) <- "RNA" #fibroblast:9,11,21
VlnPlot(sce,features = c('CD68','CD86','CD163'),pt.size = 0,group.by = 'seurat_clusters',ncol = 2)
library(randomcoloR)
allcolour <- c(pal_npg(alpha = 0.8)(9), pal_igv(alpha = 0.8)(9), pal_jama(alpha = 0.8)(7), pal_jco(alpha = 0.8)(9), pal_nejm(alpha = 0.8)(8))
# length(table(sce@active.ident))
# 28
mycolor1 <- allcolour[seq_along(table(sce$seurat_clusters))]
figs2b<-FeaturePlot(sce, features = c('CD68','CD86','CD163'), pt.size = 0.3,reduction = 'umap',ncol = 2)
figs2a<-DimPlot(sce,cols =mycolor1 ,group.by = 'seurat_clusters', reduction="umap",
label = "T", pt.size = 0.3,
label.size = 5) +
theme(axis.line = element_line(size=0.1, colour = "black"), #axis.text = element_blank(), #axis.title = element_blank(), axis.ticks = element_blank()
) +ggtitle('')
figs2ab<-ggarrange(figs2a,figs2b,nrow = 1,ncol = 2,widths = c(1,1),labels = c('A','B'))
# figs2ab
table(sce$seurat_clusters)
#
# save(sce, file = 'sce1.RData')
load('sce1.RData')
Idents(sce) <- "seurat_clusters"
sce <- subset(sce, idents =c())
Resolution <- 0.1
DefaultAssay(sce) <- "RNA"
sce <- FindNeighbors(object = sce, dims = 1:Dims)
sce <- FindClusters(object = sce, resolution = Resolution)
DefaultAssay(sce) <- "RNA"

VlnPlot(sce, features = c('CD68','CD86','CD163'), pt.size = 0, group.by = 'seurat_clusters', ncol = 2)
sce <- RunUMAP(sce, dims = 1:Dims, reduction = "pca", perplexity = 30, max_iter = 1000)
figs2c <- DimPlot(sce, cols=mycolor1, group.by = 'seurat_clusters', reduction="umap",
label = "T", pt.size = 0.2, label.size = 5) +
theme(axis.line = element_line(size=0.1, colour = "black"), axis.ticks = element_blank()
) + ggtitle('')
figs2d <- FeaturePlot(sce, features = c('CD68','CD86','CD163'), pt.size = 0.1, reduction = 'umap', ncol = 2)
fig2cd <- ggarrange(figs2c, figs2d, nrow = 1, ncol = 2, widths = c(1,1), labels = c('C','D'))
# fig2cd
figs2 <- ggarrange(figs2ab, fig2cd, nrow = 2, ncol = 1)
ggsave("results/FigS2.pdf", figs2, height = 15, width = 15)

Logfc <- 0.5
Minpct <- 0.35
DefaultAssay(sce) <- "RNA"
Idents(sce) <- "seurat_clusters"
sce.markers <- FindAllMarkers(object = sce, logfc.threshold = Logfc, min.pct = Minpct, only.pos = T)
sce.markers["pct.diff"] <- sce.markers$pct.1 - sce.markers$pct.2
sce.markers <- sce.markers[sce.markers$p_val_adj < 0.05,]
# sce.markers
# length(unique(sce.markers$gene))
head(sce.markers)

write.table(sce.markers,'results/scRNA_marker_gene.txt',quote = F,row.names = F,sep='\t')
Top5 <- sce.markers %>% group_by(cluster) %>% slice_max(n =5, order_by = avg_log2FC)
Top5 <- intersect(unique(Top5$gene),rownames(sce@assays$RNA@meta.features))
sc_marker_dotplot <- DotPlot(object = sce, features = Top5,cols=c("blue", "red"),scale = T)+
RotatedAxis()+ ggtitle("Top 5 Marker Genes")+
theme(plot.title = element_text(hjust = 0.5)) +xlab('')
sc_marker_dotplot
ggsave('results/sc_marker_dotplot.pdf',sc_marker_dotplot,height = 7,width = 9)