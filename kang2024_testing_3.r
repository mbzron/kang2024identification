load("sce.RData")
library(copykat)

copykat <- function(
  rawmat = rawdata, id.type = "S", cell.line = "no",
  ngene.chr = 0, LOW.DR = 0.05, UP.DR = 0.1,
  win.size = 25, norm.cell.names = "", KS.cut = 0.1,
  sam.name = "", distance = "euclidean", n.cores = 1
) {
  start_time <- Sys.time()
  set.seed(1)
  sample.name <- paste(sam.name, "_copykat_", sep = "")
  print("running copykat v1.0.4")
  print("step1: read and filter data ...")
  print(paste(nrow(rawmat), " genes, ", ncol(rawmat), " cells in raw data",
              sep = ""))

  der <- apply(rawmat, 1, function(x) sum(x > 0)) / ncol(rawmat)
  if (sum(der > LOW.DR) >= 1) {
    rawmat <- rawmat[which(der > LOW.DR), ]
    print(paste(nrow(rawmat), " genes past LOW.DR filtering", sep = ""))
  }

  wns1 <- "data quality is ok"
  if (nrow(rawmat) < 7000) {
    wns1 <- "low data quality"
    UP.DR <- LOW.DR
    print("WARNING: low data quality; assigned LOW.DR to UP.DR...")
  }
  print(wns1)

  print("step 2: annotations gene coordinates ...")
  anno_mat <- annotateGenes.hg20(mat = rawmat, ID.type = id.type)
  anno_mat <- anno_mat[order(anno_mat$abspos, decreasing = FALSE), ]

  # HLAs <- anno_mat$hgnc_symbol[grep("^HLA-", anno_mat$hgnc_symbol)]
  # Filtering code if necessary ...

  rawmat3 <- data.matrix(anno_mat[, 8:ncol(anno_mat)])
  norm_mat <- log(sqrt(rawmat3) + sqrt(rawmat3 + 1))
  norm_mat <- apply(norm_mat, 2, function(x) x - mean(x))
  colnames(norm_mat) <- colnames(rawmat3)

  print("step 3: smoothing data with dlm ...")
  dlm.sm <- function(c) {
    model <- dlm::dlmModPoly(order = 1, dV = 0.16, dW = 0.001)
    x <- dlm::dlmSmooth(norm_mat[, c], model)$s
    x <- x[-1] - mean(x)
    return(x)
  }

  test_mc <- parallel::mclapply(1:ncol(norm_mat), dlm.sm, mc.cores = n.cores)
  norm_mat_smooth <- matrix(unlist(test_mc), ncol = ncol(norm_mat), byrow = FALSE)
  colnames(norm_mat_smooth) <- colnames(norm_mat)

  print("step 5: segmentation...")
  results <- CNA.MCMC(clu = CL, fttmat = norm.mat.relat, bins = win.size,
                      cut.cor = KS.cut, n.cores = n.cores)

  if (length(results$breaks) < 25) {
    print("too few breakpoints detected; decreased KS.cut to 50%")
    results <- CNA.MCMC(clu = CL, fttmat = norm.mat.relat, bins = win.size,
                        cut.cor = 0.5 * KS.cut, n.cores = n.cores)
  }

  if (length(results$breaks) < 25) {
    print("too few breakpoints detected; decreased KS.cut to 75%")
    results <- CNA.MCMC(clu = CL, fttmat = norm.mat.relat, bins = win.size,
                        cut.cor = 0.5 * 0.5 * KS.cut, n.cores = n.cores)
  }

  # Additional steps with proper indentation and styling
  end_time <- Sys.time()
  print(end_time - start_time)
  return(results)
}

copykat_test <- copykat(rawmat = sce@assays$RNA@counts,
                        id.type = "S", cell.line = "no", ngene.chr = 5, 
                        win.size = 25, KS.cut = 0.15, 
                        sam.name = "LUAD", distance = "euclidean", n.cores = 1)

save(copykat_test, file = "copykat.test.RData")
copykat_test <- read.delim(
  "LUAD_copykat_prediction.txt", sep = "\t", header = TRUE
)
head(copykat_test)
table(copykat_test$copykat.pred)

rownames(copykat_test) <- copykat_test$cell.names
copykat_test <- copykat_test[rownames(sce@meta.data),]
sce <- AddMetaData(sce, copykat_test$copykat.pred, col.name = "copykat.pred")
sce$copykat.pred[is.na(sce$copykat.pred)] <- "Unknown"
table(sce$copykat.pred)

sce$copykat.pred <- ifelse(sce$copykat.pred == "aneuploid", "malignant",
                           "no_malignant")
save(sce, file = "sce.RData")

fig1h <- DimPlot(sce, cols = c("red", "blue"), group.by = "copykat.pred",
                 reduction = "umap", label = "F", pt.size = 0.5,
                 label.size = 5) +
  theme(
    axis.line = element_line(size = 0.1, colour = "black")
  )

fig1ef <- ggarrange(fig1e3, fig1f, fig1h, labels = c("D", "E", "F"),
                    nrow = 1, ncol = 3, widths = c(1.2, 1.3, 1))
fig1ab <- ggarrange(fig1a, fig1b, nrow = 1, ncol = 2, labels = c("A", "B"),
                    widths = c(1, 1.5))
fig1 <- ggarrange(fig1ab, sc_marker_dotplot, fig1ef, labels = c("", "C", ""),
                  nrow = 3, ncol = 1, heights = c(2, 1, 1))

ggsave(filename = "results/Fig1.pdf", plot = fig1, height = 15, width = 18)
ggsave(filename = "results/Fig1.jpg", plot = fig1, height = 15, width = 18)