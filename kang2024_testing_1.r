dir.create("results")
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

dir_name <- list.dirs("GSE193337_RAW/", full.names = FALSE, recursive = FALSE)

datalist <- list()
for (i in seq_along(dir_name)){
  dir_10x <- paste0("GSE193337_RAW/", dir_name[i])
  my_data <- Read10X(data.dir = dir_10x)
  colnames(my_data) <- paste0(dir_name[i], colnames(my_data))
  datalist[[i]] <- CreateSeuratObject(
    counts = my_data,
    project = dir_name[i],
    min.cells = 3, min.features = 250
  )
  datalist[[i]]$Samples <- dir_name[i]
  datalist[[i]]$type <- substr(dir_name[i], 1, 1)
}
names(datalist) <- dir_name
for (i in seq_along(datalist)){
  sce <- datalist[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  sce[["percent.Ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")
  datalist[[i]] <- sce
  # rm(sce)
}
sce <- merge(datalist[[1]], y = datalist[2:length(datalist)])
raw_cell <- sce@meta.data
raw_count <- table(raw_cell$Samples)

sum(raw_count) # 12554

pearplot_befor <- VlnPlot(
  sce,
  group.by = "Samples",
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.Ribo"),
  pt.size = 0,
  ncol = 4
)
ggsave("results/pearplot_befor.pdf", pearplot_befor, height = 5, width = 15)
# ggsave("results/pearplot_befor.jpg", pearplot_befor, height = 5, width = 15, dpi = 300)

sample_color <- pal_nejm(alpha = 0.5)(8)[1:8]

Feature_ber1 <- FeatureScatter(
  sce,
  feature1 = "nFeature_RNA",
  feature2 = "nCount_RNA",
  group.by = "Samples",
  cols = sample_color
)

Feature_ber2 <- FeatureScatter(
  sce,
  feature1 = "percent.mt",
  feature2 = "nCount_RNA",
  group.by = "Samples",
  cols = sample_color
)

Feature_ber3 <- FeatureScatter(
  sce,
  feature1 = "percent.mt",
  feature2 = "nFeature_RNA",
  group.by = "Samples",
  cols = sample_color
)

Feature_ber1 <- Feature_ber1 + theme(legend.position = "none")
Feature_ber2 <- Feature_ber2 + theme(legend.position = "none")
Feature_ber <- ggarrange(
  Feature_ber1,
  Feature_ber2,
  Feature_ber3,
  ncol = 3,
  nrow = 1,
  widths = c(1, 1, 1.2)
)
ggsave("results/Feature_cor.pdf", Feature_ber, height = 5, width = 17)
# ggsave("results/Feature_cor.jpg", Feature_ber, height = 5, width = 17, dpi = 300)

datalist <- lapply(X = datalist, FUN = function(x) {
  x <- subset(x, subset = nFeature_RNA < 5000 & percent.mt < 15)
})

sce <- merge(datalist[[1]], y = datalist[2:length(datalist)])

clean_cell <- sce@meta.data

clean_count <- table(clean_cell$Samples)

# display clean_count value with explanatory text
print(paste0("sum of clean_count (expected 6410): ", sum(clean_count)))

pearplot_after <- VlnPlot(
  sce,
  group.by = "Samples",
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.Ribo"),
  pt.size = 0,
  ncol = 4
)
ggsave("results/pearplot_after.pdf", pearplot_after, height = 5, width = 15)
# ggsave("results/pearplot_after.jpg", pearplot_after, height = 5, width = 15, dpi = 300)

save(datalist, file = "datalist.RData")
