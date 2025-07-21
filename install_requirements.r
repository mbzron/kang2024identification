# install requirements for single-cell seq analysis (split over multiple lines for readability)
install.packages(c(
  "Seurat",
  "dplyr",
  "ggplot2",
  "magrittr",
  "gtools",
  "stringr",
  "Matrix",
  "tidyverse",
  "patchwork",
  "data.table",
  "RColorBrewer",
  "ggpubr",
  "ggsci",
  "clustree",
  "randomcoloR",
  "BiocManager",
  "devtools",
  "GMD",
  "pheatmap"
), dependencies = TRUE)

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("GSEABase")
BiocManager::install("GSVA")

library(devtools)
devtools::install_github("navinlabcode/copykat")
devtools::install_github("cran/GMD")
