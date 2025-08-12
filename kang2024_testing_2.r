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

# print("Running 'SketchData'")
# sce <- SketchData(
#   object = sce,
#   ncells = 10000,
#   method = "LeverageScore",
#   sketched.assay = "sketch"
# )

# print("Running 'DefaultAssay'")
# DefaultAssay(sce) <- "sketch"

print("Running 'ScaleData'")
sce <- ScaleData(sce, features = rownames(sce))

print("Running 'RunPCA'")
sce <- RunPCA(sce, features = VariableFeatures(sce))

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

pdf("results/clust.snn_res.pdf", he = 15, wi = 15)
clustree(sce@meta.data, prefix = "RNA_snn_res.")
dev.off()

resolution <- 0.8  # according to page 3/13, resolution = 0.1.

print("Running 'FindNeighbors'")
sce <- FindNeighbors(object = sce, dims = 1:dims)

print("Running 'FindClusters'")
sce <- FindClusters(object = sce, resolution = resolution)

print("Running 'DefaultAssay(sce) <- RNA'")
DefaultAssay(sce) <- "RNA"  # fibroblast: 9, 11, 21

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

save(sce, file = "sce1.RData")
load("sce1.RData")

Idents(sce) <- "seurat_clusters"
sce <- subset(sce, idents = c())
resolution <- 0.1

DefaultAssay(sce) <- "RNA"

# print("Running 'FindNeighbors' with DefaultAssay <- sketch")
sce <- FindNeighbors(object = sce, dims = 1:dims)

# print("Running 'FindClusters' with DefaultAssay <- sketch")
sce <- FindClusters(object = sce, resolution = resolution)

DefaultAssay(sce) <- "RNA"

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

print("Running 'FindAllMarkers'")
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

top5 <- intersect(unique(top5$gene), rownames(sce@assays$RNA@features))

sc_marker_dotplot <- DotPlot(
  object = sce, features = top5, cols = c("blue", "red"), scale = TRUE
) + RotatedAxis() + ggtitle("Top 5 Marker Genes") +
  theme(plot.title = element_text(hjust = 0.5)) + xlab("")

save(sc_marker_dotplot, file = "sc_marker_dotplot.RData")
ggsave("results/sc_marker_dotplot.pdf",
       sc_marker_dotplot, height = 7, width = 9)

###
bubble_df <- as.matrix(sce[["RNA"]]$data[top5, ])
bubble_df <- t(bubble_df)
bubble_df <- as.data.frame(scale(bubble_df))
bubble_df$CB <- rownames(bubble_df)
bubble_df <- merge(
  bubble_df,
  data.frame(
    CB = rownames(sce@meta.data),
    celltype = sce@meta.data$seurat_clusters
  ),
  by = "CB"
)
bubble_df$CB <- NULL
celltype_v <- c()
gene_v <- c()
mean_v <- c()
ratio_v <- c()
for (i in unique(bubble_df$celltype)) {
  bubble_df_small <- bubble_df %>% filter(celltype == i)
  for (j in top5) {
    exp_mean <- mean(bubble_df_small[, j])
    exp_ratio <- sum(
      bubble_df_small[, j] > min(bubble_df_small[, j])
    ) / length(bubble_df_small[, j])
    celltype_v <- append(celltype_v, i)
    gene_v <- append(gene_v, j)
    mean_v <- append(mean_v, exp_mean)
    ratio_v <- append(ratio_v, exp_ratio)
  }
}

plotdf <- data.frame(
  celltype = celltype_v, gene = gene_v, exp = mean_v, ratio = ratio_v
)

plotdf$celltype <- factor(
  plotdf$celltype, levels = unique(as.character(sce_markers$cluster))
)
plotdf$gene <- factor(plotdf$gene, levels = rev(as.character(top5)))
plotdf$exp <- ifelse(plotdf$exp > 3, 3, plotdf$exp)
sc_marker_dotplot1 <- plotdf %>%
  ggplot(
    aes(x = celltype, y = gene, size = ratio, color = exp)
  ) + geom_point() + scale_x_discrete("") + scale_y_discrete("") +
  scale_color_gradientn(
    colours = rev(c("#FFD92F", "#FEE391", brewer.pal(11, "Spectral")[7:11]))
  ) + scale_size_continuous(limits = c(0, 1)) + theme_bw() + theme(
  axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
)

ggsave(
  "results/sc_marker_dotplot1.pdf", sc_marker_dotplot1, height = 7, width = 9
)

mycolor <- ggsci::pal_jama()(9)
fig1a <- DimPlot(sce, dims = c(1, 2), group.by = "Samples", reduction = "umap",
  label = "F", pt.size = 0.5,
  label.size = 5
) +  theme(
  axis.line = element_line(size = 0.1, colour = "black"),
  # axis.text = element_blank(),
  # axis.title = element_blank(),
  # axis.ticks = element_blank(),
) + ggtitle("") + guides(colour = guide_legend(ncol = 1))
save(fig1a, file = "fig1a.RData")

fig1b <- DimPlot(
  sce, cols = mycolor, group.by = "seurat_clusters",
  reduction = "umap", split.by = "type",
  label = "F", pt.size = 0.5, label.size = 5
) + theme(axis.line = element_line(size = 0.1, colour = "black"),
  # axis.text = element_blank(),
  # axis.title = element_blank(),
  # axis.ticks = element_blank(),
) + ggtitle("")
save(fig1b, file = "fig1b.RData")

Idents(sce) <- "seurat_clusters"
library("ggplot2")

sample_clust <- as.matrix(table(sce$Samples, sce$seurat_clusters))
sample_clust <- apply(sample_clust, 1, function(x){return(x/sum(x))})
sample_clust <- reshape2::melt(sample_clust)
colnames(sample_clust) <- c("cluster", "Samples", "proportion")
sample_clust$cluster <- paste0("CAF_", sample_clust$cluster)
write.table(
  sample_clust, "results/sample_clust1.txt",
  quote = FALSE, row.names = TRUE, sep = "\t"
)
clust_freq <- as.data.frame(table(sce$Samples))

colnames(clust_freq) <- c("Samples", "cell_num")

clust_freq <- clust_freq[order(clust_freq$cell_num, decreasing = TRUE), ]
clust_freq$Samples <- factor(clust_freq$Samples, levels = clust_freq$Samples)
sample_clust$Samples <- factor(
  sample_clust$Samples, levels = clust_freq$Samples
)

fig1e1 <- ggplot(
  sample_clust, aes(x = Samples, y = proportion, fill = cluster)
) + geom_bar(stat = "identity", position = "fill") +
  ggtitle("") + scale_fill_manual(values = mycolor) +
  theme_bw() + theme(
  axis.ticks.length = unit(0.1, "cm"),
  legend.position = "left"
) + xlab("") +
  coord_flip() + scale_y_continuous(expand = expand_scale(mult = c(0, 0))
)

sample_color <- pal_nejm(alpha = 0.5)(8)[1:8]
fig1e2 <- ggplot(clust_freq, aes(x = Samples, y = cell_num, fill = Samples)) +
  geom_bar(stat = "identity") + ggtitle("") +
  theme_bw() + scale_fill_manual(values = sample_color) +
  theme(
    axis.ticks.length = unit(0, "cm"),
    axis.text.y = element_blank(),
    axis.title.y = element_blank()
  ) + coord_flip() + scale_y_continuous(
  expand = expand_scale(mult = c(0, 0))
) + ylim(0, max(clust_freq$cell_num) + 10)

fig1e3 <- ggpubr::ggarrange(
  fig1e1,
  fig1e2,
  nrow = 1,
  ncol = 2,
  widths = c(2, 1)
)
# save fig1e3 object for later use
save(fig1e3, file = "fig1e3.RData")

library(clusterProfiler)
library(org.Hs.eg.db)
ids <- bitr(
  sce_markers$gene, "SYMBOL", "ENTREZID", "org.Hs.eg.db"
) ## 将 SYMBOL 转成ENTREZID
sce_markers2 <- merge(sce_markers, ids, by.x = "gene", by.y = "SYMBOL")
gcsample <- split(sce_markers2$ENTREZID, sce_markers2$cluster)

## KEGG
sce_markers2_enrich_res <- compareCluster(
  gcsample, fun = "enrichKEGG",
  organism = "hsa", pvalueCutoff = 0.05
)

fig1f <- dotplot(sce_markers2_enrich_res) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
  axis.text.y = element_text(size = 10)
)
save(fig1f, file = "fig1f.RData")

save(sce, file = "sce.RData")