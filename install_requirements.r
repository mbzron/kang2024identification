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
), dependencies = TRUE)

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
