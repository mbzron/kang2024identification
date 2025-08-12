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

source("logging.r")

use_existing_files <- TRUE
data_subset_proportion <- 0.1

n_cores <- 6
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

use_existing_norm.mat.smooth <- FALSE
if (use_existing_norm.mat.smooth && file.exists("norm_mat_smooth.RData")) {
  log_message("Loading existing norm.mat.smooth.RData ...")
  norm.mat.smooth <- load("norm_mat_smooth.RData")
} else {
  log_message("Creating new norm.mat.smooth.RData with parallel::mclapply...")
  test.mc <- parallel::mclapply(
    seq_len(ncol(norm.mat)), dlm.sm, mc.cores = n_cores
  )
  norm.mat.smooth <- matrix(
    unlist(test.mc), ncol = ncol(norm.mat), byrow = FALSE
  )
  colnames(norm.mat.smooth) <- colnames(norm.mat)
  save(norm.mat.smooth, file = "norm_mat_smooth.RData")
}

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

WNS <- basa$WNS
preN <- basa$preN
CL <- basa$cl
basel <- basa$basel

###############################################
###############################################

if (WNS == "unclassified.prediction") {
  Tc <- colnames(rawmat)[
    which(
      as.numeric(
        apply(
          rawmat[which(rownames(rawmat) %in% c("PTPRC", "LYZ", "PECAM1")), ],
          2, mean
        )
      ) > 1
    )
  ]

  preN <- intersect(Tc, colnames(norm.mat.smooth))

  if (length(preN) > 5) {
    print("start manual mode")
    WNS <- paste(
      "copykat failed in locating normal cells; manual adjust performed with",
      length(preN),
      " immune cells",
      sep = ""
    )
    print(paste0("WNS: ", WNS))
    basel <- apply(
      norm.mat.smooth[, which(colnames(norm.mat.smooth) %in% preN)],
      1, mean
    )
  } else {
    print("running baseline.GMM")
    basa <- baseline.GMM(
      CNA.mat = norm.mat.smooth,
      max.normal = 5,
      mu.cut = 0.05,
      Nfraq.cut = 0.99,
      RE.before = basa,
      n.cores = n_cores
    )
    basel <- basa$basel
    WNS <- basa$WNS
    preN <- basa$preN
  }
}

rm(basa)

norm.mat.relat <- norm.mat.smooth - basel

DR2 <- apply(rawmat3, 1, function(x) (sum(x > 0))) / ncol(rawmat3)
norm.mat.relat <- norm.mat.relat[which(DR2 >= UP.DR), ]
anno.mat2 <- anno.mat[which(DR2 >= UP.DR), ]

CL <- CL[which(names(CL) %in% colnames(norm.mat.relat))]
CL <- CL[order(match(names(CL), colnames(norm.mat.relat)))]

# # possible alernative syntax:
# common <- intersect(colnames(norm.mat.relat), names(CL))
# CL <- CL[match(common, names(CL))]

# log_message("step 5: segmentation...")

# sam_name <- "LUAD"
# sample.name <- paste(sam_name, "_copykat_", sep = "")
# luad_filename <- paste(sample.name, "CNA_raw_results.txt", sep = "")

# if (use_existing_files && file.exists(luad_filename)) {
#   log_message(paste("Loading existing", luad_filename, "..."))
#   RNA.copycat <- read.delim(luad_filename, sep = "\t", header = TRUE)
# } else {
#   log_message(paste("Creating new", luad_filename, "..."))

#   win_size <- 25
#   KS.cut = 0.1
#   results <- CNA.MCMC(
#     clu = CL,
#     fttmat = norm.mat.relat,
#     bins = win_size,
#     cut.cor = KS.cut,
#     n.cores = n_cores
#   )

#   if (length(results$breaks) < 25) {
#     log_message("too few breakpoints detected; decreased KS.cut to 50%")
#     results <- CNA.MCMC(
#       clu = CL,
#       fttmat = norm.mat.relat,
#       bins = win_size,
#       cut.cor = 0.5 * KS.cut,
#       n.cores = n_cores
#     )
#   }

#   if (length(results$breaks) < 25) {
#     log_message("too few breakpoints detected; decreased KS.cut to 75%")
#     results <- CNA.MCMC(
#       clu = CL,
#       fttmat = norm.mat.relat,
#       bins = win_size,
#       cut.cor = 0.5 * 0.5 * KS.cut,
#       n.cores = n_cores
#     )
#   }

#   if (length(results$breaks) < 25) {
#     log_message(
#       stop(
#         "too few breakpoints detected; 
#         please check your data quality or increase KS.cut"
#       )
#     )
#   }

#   colnames(results$logCNA) <- colnames(norm.mat.relat)
#   results.com <- apply(results$logCNA, 2, function(x) (x <- x - mean(x)))
#   RNA.copycat <- cbind(anno.mat2[, 1:7], results.com)

#   write.table(
#     RNA.copycat,
#     luad_filename,
#     sep = "\t",
#     row.names = FALSE,
#     quote = FALSE
#   )

# }

# if (!file.exists("Aj.rds")) {
#   log_message("step 6: convert to genomic bins...")
#   Aj <- copykat::convert.all.bins.hg20(
#     DNA.mat = copykat::DNA.hg20,
#     RNA.mat = RNA.copycat,
#     n.cores = n_cores
#   )
#   saveRDS(Aj, file = "Aj.rds")
# } else {
#   log_message("step 6: Loading existing Aj.rds ...")
#   Aj <- readRDS("Aj.rds")
# }

# rm(RNA.copycat)  # this isn't needed anymore

# ###############################################
# ###############################################



# use_existing_copykat_test <- FALSE
# if (!file.exists("copykat.test.RData") || !use_existing_copykat_test) {
#   source("copykat_breakdown.r")
#   copykat_test <- copykat_simplified(
#     Aj = Aj,
#     n_cores = n_cores,
#     WNS = WNS,
#     subset_percentage = data_subset_proportion, plot_heatmap = FALSE
#   )
#   save(copykat_test, file = "copykat.test.RData")
# } else {
#   log_message("Loading existing copykat.test.RData to variable copykat_test...")
#   load("copykat.test.RData")
# }

# copykat_test <- read.delim(
#   "LUAD_copykat_prediction.txt", sep = "\t", header = TRUE
# )
# head(copykat_test)
# table(copykat_test$copykat.pred)

# rownames(copykat_test) <- copykat_test$cell.names

# print(nrow(copykat_test))
# print(nrow(sce@meta.data))
# print(any(rownames(sce@meta.data) %in% rownames(copykat_test)))

# copykat_test <- copykat_test[rownames(sce@meta.data), ]
# sce <- AddMetaData(sce, copykat_test$copykat.pred, col.name = "copykat.pred")
# sce$copykat.pred[is.na(sce$copykat.pred)] <- "Unknown"
# # table(sce$copykat.pred)

# sce$copykat.pred <- ifelse(
#   sce$copykat.pred == "aneuploid",
#   "malignant",
#   "no_malignant"
# )

# log_message("saving SCE object with copykat predictions...")
# save(sce, file = "sce.RData")

# fig1h <- DimPlot(
#   sce,
#   cols = c("red", "blue"),
#   group.by = "copykat.pred",
#   reduction = "umap",
#   label = "F",
#   pt.size = 0.5,
#   label.size = 5
# ) + theme(
#   axis.line = element_line(
#     size = 0.1,
#     colour = "black"
#   )
# )

# log_message("loading precomputed figures...")
# load("fig1e3.RData")
# load("fig1f.RData")

# fig1ef <- ggarrange(fig1e3, fig1f, fig1h, labels = c("D", "E", "F"),
#                     nrow = 1, ncol = 3, widths = c(1.2, 1.3, 1))

# load("fig1a.RData")
# load("fig1b.RData")

# fig1ab <- ggarrange(fig1a, fig1b, nrow = 1, ncol = 2, labels = c("A", "B"),
#                     widths = c(1, 1.5))


# load("sc_marker_dotplot.RData")
# fig1 <- ggarrange(fig1ab, sc_marker_dotplot, fig1ef, labels = c("", "C", ""),
#                   nrow = 3, ncol = 1, heights = c(2, 1, 1))

# log_message("saving figure 1...")
# ggsave(filename = "results/Fig1.pdf", plot = fig1, height = 15, width = 18)

# # ggsave(filename = "results/Fig1.jpg", plot = fig1, height = 15, width = 18)