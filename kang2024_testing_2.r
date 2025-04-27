options(stringsAsFactors = FALSE, check.bounds = FALSE)
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
library(clustree)
library(randomcoloR)

load("datalist.RData")

sce <- merge(datalist[[1]], y = datalist[2:length(datalist)])

print("Running 'NormalizeData'")
sce <- NormalizeData(sce, normalization.method = "LogNormalize",
  scale.factor = 10000
)

print("Running 'FindVariableFeatures'")
sce <- FindVariableFeatures(sce, selection.method = "vst",
  nfeatures = 2000,
  mean.cutoff = c(0.0125, 3),
  dispersion.cutoff = c(1.5, Inf)
)

print("Running 'SketchData'")
sce <- SketchData(
  object = sce,
  ncells = 10000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

print("Running 'DefaultAssay'")
DefaultAssay(sce) <- "sketch"

print("Running 'ScaleData'")
sce <- ScaleData(sce)  # , features = rownames(sce))

print("Running 'RunPCA'")
sce <- RunPCA(sce)  # , features = VariableFeatures(sce))

elbowplot <- ElbowPlot(sce, ndims = 50, reduction = "pca")
ggsave("results/elbowplot.pdf", elbowplot, height = 5, width = 5)

print("Running 'RunUMAP'")
dims <- 30
sce <- RunUMAP(sce, dims = 1:dims, reduction = "pca")
raw_umap <- DimPlot(sce, group.by = "Samples",
                    reduction = "umap",
                    label = "T",
                    pt.size = 0.2,
                    label.size = 0) + ggtitle("")
ggsave("results/raw.umap.pdf", raw_umap, height = 7, width = 7)

print("Running 'FindNeighbors'")
sce <- FindNeighbors(sce, dims = 1:dims)

print("Running 'FindClusters'")
sce <- FindClusters(
  object = sce, resolution = c(seq(.1, 1, .1))
)

# colnames(sce@meta.data)
# clustree(sce@meta.data, prefix = "RNA_snn_res.")
# pdf("results/clust.snn_res.pdf", he = 15, wi = 15)
# clustree(sce@meta.data, prefix = "RNA_snn_res.")

resolution <- 0.8

print("Running 'FindNeighbors'")
sce <- FindNeighbors(object = sce, dims = 1:dims)

print("Running 'FindClusters'")
sce <- FindClusters(object = sce, resolution = resolution)

# print("Running 'DefaultAssay(sce) <- RNA'")
# DefaultAssay(sce) <- "RNA"  # fibroblast: 9, 11, 21

print("Running 'VlnPlot'")
VlnPlot(sce, features = c("CD68", "CD86", "CD163"),
        pt.size = 0, group.by = "seurat_clusters", ncol = 2)

allcolour <- c(
  pal_npg(alpha = 0.8)(9),
  pal_igv(alpha = 0.8)(9),
  pal_jama(alpha = 0.8)(7),
  pal_jco(alpha = 0.8)(9),
  pal_nejm(alpha = 0.8)(8)
)

length(table(sce@active.ident))  # 28
mycolor1 <- allcolour[seq_along(table(sce$seurat_clusters))]

figs2b <- FeaturePlot(sce, features = c("CD68", "CD86", "CD163"),
                      pt.size = 0.3, reduction = "umap", ncol = 2)

figs2a <- DimPlot(sce, cols = mycolor1,
                  group.by = "seurat_clusters", reduction = "umap",
                  label = "T", pt.size = 0.3,
                  label.size = 5) +
  theme(axis.line = element_line(size = 0.1, colour = "black"),
  ) + ggtitle("")

figs2ab <- ggarrange(figs2a, figs2b, nrow = 1,
                     ncol = 2, widths = c(1, 1), labels = c("A", "B"))

table(sce$seurat_clusters)

# save(sce, file = 'sce1.RData')
# load('sce1.RData')

Idents(sce) <- "seurat_clusters"
sce <- subset(sce, idents = c())
resolution <- 0.1

# DefaultAssay(sce) <- "RNA"

print("Running 'FindNeighbors' with DefaultAssay <- sketch")
sce <- FindNeighbors(object = sce, dims = 1:dims)

print("Running 'FindClusters' with DefaultAssay <- sketch")
sce <- FindClusters(object = sce, resolution = resolution)

# DefaultAssay(sce) <- "RNA"

VlnPlot(sce, features = c("CD68", "CD86", "CD163"),
  pt.size = 0, group.by = "seurat_clusters", ncol = 2
)

sce <- RunUMAP(sce, dims = 1:dims, reduction = "pca",
               perplexity = 30, max_iter = 1000)

figs2c <- DimPlot(sce, cols = mycolor1,
  group.by = "seurat_clusters", reduction = "umap",
  label = "T", pt.size = 0.2, label.size = 5
) + theme(
  axis.line = element_line(size = 0.1, colour = "black"),
  axis.ticks = element_blank()
) + ggtitle("")

figs2d <- FeaturePlot(sce, features = c("CD68", "CD86", "CD163"),
                      pt.size = 0.1, reduction = "umap", ncol = 2)

fig2cd <- ggarrange(figs2c, figs2d,
                    nrow = 1, ncol = 2, widths = c(1, 1), labels = c("C", "D"))

# fig2cd
figs2 <- ggarrange(figs2ab, fig2cd, nrow = 2, ncol = 1)
ggsave("results/FigS2.pdf", figs2, height = 15, width = 15)

logfc <- 0.5  # original is 0.5
minpct <- 0.35  # original is 0.35
DefaultAssay(sce) <- "RNA"
Idents(sce) <- "seurat_clusters"

sce <- JoinLayers(sce)

print("Running 'FindAllMarkers")
sce_markers <- FindAllMarkers(object = sce, logfc.threshold = logfc,
                              min.pct = minpct, only.pos = TRUE)

sce_markers["pct.diff"] <- sce_markers$pct.1 - sce_markers$pct.2
sce_markers <- sce_markers[sce_markers$p_val_adj < 0.05, ]

length(unique(sce_markers$gene))
head(sce_markers)

write.table(sce_markers, "results/scRNA_marker_gene.txt",
            quote = FALSE, row.names = FALSE, sep = "\t")

top5 <- sce_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  ungroup()

top5 <- intersect(unique(top5$gene), rownames(sce@assays$RNA))

sc_marker_dotplot <- DotPlot(
  object = sce, features = top5, cols = c("blue", "red"), scale = TRUE
) + RotatedAxis() + ggtitle("Top 5 Marker Genes") +
  theme(plot.title = element_text(hjust = 0.5)) + xlab("")

ggsave("results/sc_marker_dotplot.pdf",
       sc_marker_dotplot, height = 7, width = 9)