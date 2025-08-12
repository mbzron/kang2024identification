load("sce.RData")
library(copykat)
library(Seurat)
library(ggplot2)
library(ggpubr)

source("logging.r")

use_existing_files <- TRUE
data_subset_proportion <- 0.1
n_cores <- 6

log_message("Loading existing basa.rds ...")
basa <- readRDS("basa.rds")

WNS <- basa$WNS
preN <- basa$preN
CL <- basa$cl
basel <- basa$basel

log_message("step 6: Loading existing Aj.rds ...")
Aj <- readRDS("Aj.rds")

###############################################
###############################################

use_existing_copykat_test <- FALSE
if (!file.exists("copykat.test.RData") || !use_existing_copykat_test) {
  source("copykat_breakdown.r")
  copykat_test <- copykat_simplified(
    Aj = Aj,
    n_cores = n_cores,
    WNS = WNS,
    subset_percentage = data_subset_proportion, plot_heatmap = FALSE
  )
  save(copykat_test, file = "copykat.test.RData")
} else {
  log_message("Loading existing copykat.test.RData to variable copykat_test...")
  load("copykat.test.RData")
}

copykat_test <- read.delim(
  "LUAD_copykat_prediction.txt", sep = "\t", header = TRUE
)
# head(copykat_test)
# table(copykat_test$copykat.pred)

rownames(copykat_test) <- copykat_test$cell.names

# print(nrow(copykat_test))
# print(nrow(sce@meta.data))
# print(any(rownames(sce@meta.data) %in% rownames(copykat_test)))

copykat_test <- copykat_test[rownames(sce@meta.data), ]
sce <- AddMetaData(sce, copykat_test$copykat.pred, col.name = "copykat.pred")
sce$copykat.pred[is.na(sce$copykat.pred)] <- "Unknown"
# table(sce$copykat.pred)

sce$copykat.pred <- ifelse(
  sce$copykat.pred == "aneuploid",
  "malignant",
  "no_malignant"
)

log_message("saving SCE object with copykat predictions...")
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

log_message("loading precomputed figures...")
load("fig1e3.RData")
load("fig1f.RData")

fig1ef <- ggarrange(fig1e3, fig1f, fig1h, labels = c("D", "E", "F"),
                    nrow = 1, ncol = 3, widths = c(1.2, 1.3, 1))

load("fig1a.RData")
load("fig1b.RData")

fig1ab <- ggarrange(fig1a, fig1b, nrow = 1, ncol = 2, labels = c("A", "B"),
                    widths = c(1, 1.5))


load("sc_marker_dotplot.RData")
fig1 <- ggarrange(fig1ab, sc_marker_dotplot, fig1ef, labels = c("", "C", ""),
                  nrow = 3, ncol = 1, heights = c(2, 1, 1))

log_message("saving figure 1...")
ggsave(filename = "results/Fig1.pdf", plot = fig1, height = 15, width = 18)

# ggsave(filename = "results/Fig1.jpg", plot = fig1, height = 15, width = 18)