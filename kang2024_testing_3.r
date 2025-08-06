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
  "sc_marker_dotplot",
  "cyclegenes"
)
variables_to_remove <- c()
variables_to_remove <- setdiff(ls(), keep)
rm(list = variables_to_remove)

stopifnot("sce" %in% ls())

# Define log file
log_file <- "my_log.txt"

# Function to write to log
log_message <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  full_message <- paste0("[", timestamp, "] ", message, "\n")
  print(full_message)
  cat(full_message, file = log_file, append = TRUE)
}

n_cores <- 6
use_existing_files <- TRUE
data_subset_proportion <- 0.1
id.type <- "S"
LOW.DR <- 0.05
UP.DR <- 0.1

log_message("step1: read and filter data ...")

rawmat <- sce@assays$RNA$counts

der <- apply(rawmat, 1, function(x) (sum(x > 0))) / ncol(rawmat)

if (sum(der > LOW.DR) >= 1) {
  rawmat <- rawmat[which(der > LOW.DR), ]
  log_message(paste(nrow(rawmat), " genes past LOW.DR filtering", sep = ""))
}

WNS1 <- "data quality is ok"

if (nrow(rawmat) < 7000) {
  WNS1 <- "low data quality"
  UP.DR <- LOW.DR
  log_message("WARNING: low data quality; assigned LOW.DR to UP.DR...")
}

log_message("step 2: annotations gene coordinates ...")

anno.mat <- annotateGenes.hg20(mat = rawmat, ID.type = id.type)
anno.mat <- anno.mat[order(anno.mat$abspos, decreasing = FALSE), ]
HLAs <- anno.mat$hgnc_symbol[grep("^HLA-", anno.mat$hgnc_symbol)]
toRev <- which(anno.mat$hgnc_symbol %in% c(as.vector(cyclegenes[[1]]), HLAs))

rawmat3 <- data.matrix(anno.mat[, 8:ncol(anno.mat)])
norm.mat <- log(sqrt(rawmat3) + sqrt(rawmat3 + 1))
norm.mat <- apply(norm.mat, 2, function(x) (x <- x - mean(x)))
colnames(norm.mat) <- colnames(rawmat3)

log_message("step 3: smoothing data with dlm ...")
dlm.sm <- function(c) {
  model <- dlm::dlmModPoly(order = 1, dV = 0.16, dW = 0.001)
  x <- dlm::dlmSmooth(norm.mat[, c], model)$s
  x <- x[2:length(x)]
  x <- x - mean(x)
}

test.mc <- parallel::mclapply(
  seq_len(ncol(norm.mat)), dlm.sm, mc.cores = n_cores
)
norm.mat.smooth <- matrix(unlist(test.mc), ncol = ncol(norm.mat), byrow = FALSE)
colnames(norm.mat.smooth) <- colnames(norm.mat)

log_message("step 4: measuring baselines ...")

if (use_existing_files && file.exists("basa.rds")) {
  log_message("Loading existing basa.rds ...")
  basa <- readRDS("basa.rds")
} else {
  log_message("Creating new basa.rds ...")
  basa <- baseline.norm.cl(
    norm.mat.smooth = norm.mat.smooth,
    min.cells = 5,
    n.cores = n_cores
  )
  saveRDS(basa, file = "basa.rds")
}

if (!file.exists("copykat.test.RData")) {
  source("copykat_breakdown.r")
  copykat_test <- copykat_simplified(
    basa = basa,
    norm.mat.smooth = norm.mat.smooth, # not sure about this one
    rawmat3 = rawmat3, # not sure about this one
    n_cores = n_cores,
    subset_percentage = data_subset_proportion, plot_heatmap = FALSE
  )
  save(copykat_test, file = "copykat.test.RData")
} else {
  log_message("Loading existing copykat.test.RData to variable copykat_test...")
  load("copykat.test.RData")
}

# copykat_test <- read.delim(
#   "LUAD_copykat_prediction.txt", sep = "\t", header = TRUE
# )
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