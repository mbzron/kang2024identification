load("sce.RData")
library(copykat)
library(Seurat)
library(ggplot2)
library(ggpubr)

keep <- c(
  "sce",
  "fig1e3",
  "fig1f",
  "fig1a",
  "fig1b",
  "sc_marker_dotplot"
)
variables_to_remove <- c()
variables_to_remove <- setdiff(ls(), keep)
rm(list = variables_to_remove)

stopifnot("sce" %in% ls())

# source("copykat.r")

# copykat_test <- copykat(
#   rawmat = sce@assays$RNA$counts,
#   id_type = "S",
#   cell_line = "no",
#   ngene_chr = 5,
#   win_size = 25,
#   KS.cut = 0.15,
#   sam_name = "LUAD",
#   distance = "euclidean",
#   n.cores = 6
# )

if (!file.exists("copykat.test.RData")) {
  source("copykat_breakdown.r")
  copykat_test <- copykat_simplified(
    n_cores = 6, subset_percentage = 0.25, plot_heatmap = FALSE
  )
  save(copykat_test, file = "copykat.test.RData")
}

copykat_test <- read.delim(
  "LUAD_copykat_prediction.txt", sep = "\t", header = TRUE
)
head(copykat_test)
table(copykat_test$copykat.pred)

rownames(copykat_test) <- copykat_test$cell.names
copykat_test <- copykat_test[rownames(sce@meta.data), ]
sce <- AddMetaData(sce, copykat_test$copykat.pred, col.name = "copykat.pred")
sce$copykat.pred[is.na(sce$copykat.pred)] <- "Unknown"
table(sce$copykat.pred)

sce$copykat.pred <- ifelse(
  sce$copykat.pred == "aneuploid",
  "malignant",
  "no_malignant"
)

save(sce, file = "sce.RData")

fig1h <- DimPlot(
  sce,
  cols = c("red", "blue"),
  group.by = "copykat.pred",
  reduction = "umap",
  label = "F",
  pt.size = 0.5,
  label.size = 5
) + theme(
  axis.line = element_line(
    size = 0.1,
    colour = "black"
  )
)

fig1ef <- ggarrange(fig1e3, fig1f, fig1h, labels = c("D", "E", "F"),
                    nrow = 1, ncol = 3, widths = c(1.2, 1.3, 1))
fig1ab <- ggarrange(fig1a, fig1b, nrow = 1, ncol = 2, labels = c("A", "B"),
                    widths = c(1, 1.5))
fig1 <- ggarrange(fig1ab, sc_marker_dotplot, fig1ef, labels = c("", "C", ""),
                  nrow = 3, ncol = 1, heights = c(2, 1, 1))

ggsave(filename = "results/Fig1.pdf", plot = fig1, height = 15, width = 18)
# ggsave(filename = "results/Fig1.jpg", plot = fig1, height = 15, width = 18)