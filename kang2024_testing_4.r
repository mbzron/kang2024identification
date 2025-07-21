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

# load('../01.scRNA/sce.RData')
load('sce.RData')

ssGSEAScore_by_genes <- function(gene.exp, genes) {

  gs <- GSEABase::GeneSet(
    setName = "GeneSet",
    setIdentifier = paste0("101"),
    geneIds = unique(genes),
    GSEABase::SymbolIdentifier()
  )

  gsc <- GSEABase::GeneSetCollection(list(gs))
  fl <- tempfile()
  GSEABase::toGmt(gsc, fl)
  cgeneset <- GSEABase::getGmt(fl)

  ssGSEA.geneset <- GSVA::gsva(
    as.matrix(gene.exp),
    cgeneset,
    method = "ssgsea",
    min.sz = 1,
    max.sz = 5000,
    verbose = TRUE
  )

}

pathway.score <- function(exp, gene) {

  filename <- "pathway_scores.txt"

  if (file.exists(filename)) {

    # load pathway scores from file
    print(paste0("Loading pathway scores from ", filename))
    pathway_score <- read.table(
      filename,       # The name of your file
      sep = "\t",     # Tab-separated file
      header = TRUE,  # Assuming the first row is headers/column names
      row.names = 1,  # First column contains row names
      quote = ""      # Specify quote as "" if no quotes were used
    )

  } else {

    pathway_score <- data.frame()

    for (i in unique(gene[, 2])){

      # display progress
      print(paste0("Processing pathway: ", i))

      gene_set <- gene[gene[, 2] == i, 1]
      score <- ssGSEAScore_by_genes(exp, gene_set)
      rownames(score) <- i
      pathway_score <- rbind.data.frame(pathway_score, score)

    }

    # write pathway_scores to file
    write.table(
      pathway_score,
      filename,
      quote = FALSE,
      row.names = TRUE,
      sep = "\t"
    )

    print(paste0("Pathway scores saved to ", filename))

  }

  t(pathway_score)

}

#_pmid_29625050
tumor.pathway <- read.delim(
  "pmid_29625050_pathway.txt", sep = "\t", header = TRUE
)
head(tumor.pathway)

tumor.pathway <- tumor.pathway[, c("Gene", "OG.TSG")]

#每一个细胞计算得分
tumor.pathway.score <- pathway.score(
  exp = as.matrix(sce@assays$RNA$counts), gene = tumor.pathway
)
head(tumor.pathway.score)

tumor.pathway.score.group <- merge(
  data.frame(cell.names = rownames(sce@meta.data), sce@meta.data),
  data.frame(cell.names = rownames(tumor.pathway.score), tumor.pathway.score),
  by = "cell.names"
)

rownames(tumor.pathway.score.group) <- tumor.pathway.score.group$cell.names
head(tumor.pathway.score.group)

tumor.pathway.score.group <- tumor.pathway.score.group[, -1]

copykat.test <- data.frame(
  cell.names = rownames(sce@meta.data),
  copykat.pred = sce@meta.data$copykat.pred
)
head(copykat.test)
table(copykat.test$copykat.pred)

tumor.score.copy <- cbind.data.frame(
  tumor.pathway.score.group[copykat.test$cell.names, ],
  copykat.pred = copykat.test$copykat.pred
)

head(tumor.score.copy)
table(tumor.score.copy$copykat.pred)
tumor.score.copy$seurat_clusters <- paste0('CAF_', tumor.score.copy$seurat_clusters)

library(pheatmap)
head(tumor.score.copy)
table(tumor.score.copy$copykat.pred)

print(colnames(tumor.score.copy))

tumor.score.copy <- tumor.score.copy[,
  c("Samples",
    "seurat_clusters",
    "Cell.Cycle",
    "HIPPO",
    "MYC",
    "NOTCH",
    "NRF1",
    "PI3K",
    "TGF.Beta",
    "RAS",
    "TP53",
    "WNT",
    "copykat.pred"
  )
]

# print(colnames(tumor.score.copy)[9])
# colnames(tumor.score.copy)[9] <- "TGF-Beta"
mat <- tumor.score.copy[, as.character(unique(tumor.pathway$OG.TSG))]

anno_col <- tumor.score.copy[, c("seurat_clusters", "copykat.pred")]

anno_col <- anno_col[order(anno_col$copykat.pred, anno_col$seurat_clusters), ]

pdf("results/Fig2a.pdf", height = 9, width = 12)
pheatmap::pheatmap(
  t(mat[rownames(anno_col), ]),
  scale = "row",
  show_colnames = FALSE,
  annotation_col = anno_col,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  annotation_names_row = FALSE,
  annotation_colors = list(
    copykat.pred = c("malignant" = "red", "no_malignant" = "blue")
  ),
  breaks = unique(c(seq(-2, 2, length = 100)))
)
dev.off()

write.table(
  tumor.score.copy,
  "results/tumor.score.copy.txt",
  quote = FALSE,
  row.names = TRUE,
  sep = "\t"
)

clust.malig <- table(
  tumor.score.copy$copykat.pred,
  tumor.score.copy$seurat_clusters
)

write.table(
  clust.malig,
  "results/clust.malig.txt",
  quote = FALSE,
  sep = "\t",
  row.names = TRUE
)

library(ggplot2)

plotiBarplot <- function(
  dat,
  palette,
  ist = FALSE,
  margin = TRUE,
  lineCol = "black",
  legTitle = "Group",
  showValue = FALSE,
  showLine = TRUE) {

  xlb <- ""
  ylb <- ""
  lineW <- 0.5
  xangle <- 0
  isAuto <- TRUE
  #library(tidyverse)
  #library(reshape2)
  #library(optparse)

  if (ist) {
    dat <- t(dat)
  }
  lbc <- colnames(dat)
  lbr <- row.names(dat)
  bk_dat <- dat
  if (margin) {
    dat <- dat %*% diag(1 / c(apply(t(dat), 1, sum)))
  }
  row.names(dat) <- paste0("R", seq_len(nrow(dat)))
  colnames(dat) <- paste0("C", seq_len(ncol(dat)))
  row.names(bk_dat) <- paste0("R", seq_len(nrow(bk_dat)))
  colnames(bk_dat) <- paste0("C", seq_len(ncol(bk_dat)))
  # df <- cbind(bg = paste0("R", seq_len(nrow(dat)), dat)
  # colnames(df) <- c("bg", paste0("C", seq_len(ncol(dat))))

  tp.dat <- as.data.frame(cbind(bg = row.names(dat), dat))
  tp.dat[, 1] <- as.character(tp.dat[, 1])
  for (i in 2:ncol(tp.dat)) {
    tp.dat[, i] <- as.numeric(as.character(tp.dat[, i]))
  }

  mt.df <- reshape2::melt(tp.dat)
  colnames(mt.df) <- c("bg", "variable", "value")
  pg <- ggplot(mt.df, aes(x = variable, y = value, fill = bg)) +
    geom_bar(stat = "identity", width = lineW, col = lineCol)
  if (showLine) {
    for (i in 2:(ncol(tp.dat) - 1)) {
      tmp <- tp.dat[order(tp.dat[, 1], decreasing = TRUE), ]
      tmp[, i] <- base::cumsum(tmp[, i])
      tmp[, i + 1] <- base::cumsum(tmp[, i + 1])
      colnames(tmp)[c(i, i + 1)] <- c("STY", "ED")
      tmp1 <- cbind(
        tmp,
        STX = rep(i - 1 + lineW / 2, nrow(tmp)),
        EDX = rep(i - lineW / 2, nrow(tmp))
      )

      pg <- pg + geom_segment(
        data = tmp1,
        aes(x = STX, xend = EDX, y = STY, yend = ED)
      )
    }
  }
  if (showValue) {
    pg <- pg + geom_text(
      data = mt.df,
      aes(label = sprintf("%0.2f", round(value, digits = 2))),
      position = position_stack(vjust = 0.5)
    )
  }
  pg <- pg + scale_x_discrete(
    breaks = paste0("C", seq_len(ncol(dat))),
    label = lbc
  )
  pg <- pg + labs(x = xlb, y = ylb) + theme(legend.position = "bottom")
  # pg <- pg + scale_fill_discrete(breaks = paste0("R", seq_len(dat)), label = lbr, name=legTitle)
  pg <- pg + scale_fill_manual(
    breaks = paste0("R", seq_len(nrow(dat))),
    label = lbr,
    name = legTitle,
    values = palette
  )
  if (xangle > 0) {
    pg <- pg + theme(
      axis.text.x = element_text(angle = xangle, hjust = 1),
      legend.position = "bottom"
    )
  }
  g.tb <- matrix(0, nrow = ncol(dat), ncol = ncol(dat))
  for (i in seq_len(ncol(dat))) {
    for (j in seq_len(ncol(dat))) {
      if (i != j) {
        g.tb[i, j] <- round(-log10((chisq.test(bk_dat[, c(i, j)])$p.value)), 2)
      }
    }
  }
  colnames(g.tb) <- lbc
  row.names(g.tb) <- lbc
  g.tb <- reshape2::melt(g.tb)
  colnames(g.tb) <- c("A1", "A2", "A3")
  g.tb$A4 <- paste0(g.tb[, 3], ifelse(g.tb[, 3] > -log10(0.05), "(*)", ""))
  stable.p <- ggplot(g.tb, aes(A1, A2)) +
  geom_tile(aes(fill = A3), colour = "white") +
  xlab("") +
  ylab("") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  geom_text(aes(x = A1, y = A2, label = A4)) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )
  stable.p <- stable.p + ggtitle("-log10(anova p value)")
  if (isAuto) {
    g1 <- ggpubr::ggarrange(
      stable.p,
      pg,
      ncol = 1,
      nrow = 2,
      heights = c(0.5, 1),
      align = "hv"
    )
    return(g1)
  } else {
    return(list(Bar = pg, Table = stable.p))
  }
}

#cellcolor <- ggsci::pal_jama()(9)
fig2b <- plotiBarplot(
  dat = clust.malig,
  palette = c("red", "blue"),
  ist = FALSE,
  margin = TRUE,
  lineCol = "black",
  legTitle = "Predict",
  showValue = TRUE,
  showLine = TRUE
)

ggsave("results/Fig2b.pdf", fig2b, height = 7, width = 10)
head(tumor.pathway.score.group)
table(tumor.score.copy$seurat_clusters)

#C0
tumor.pathway.score.group1 <- tumor.score.copy[
  tumor.score.copy$seurat_clusters == "CAF_0",
]

library(ggpubr)
library(reshape2)

# 箱线图
Muti_Boxplot <- function(
  dat,
  group,
  group_cols,
  leg,
  test_method = 'wilcox.test',
  ylabs
) {
  dat1 <- reshape2::melt(cbind.data.frame(dat, group))
  p <- ggboxplot(
    dat1,
    x = "variable",
    y = "value",
    fill = "group",
    color = "black",
    palette = group_cols,
    ylab = ylabs,
    xlab = "",
    add = "boxplot"
  ) +
  stat_compare_means(
    mapping = aes(group = group),
    data = dat1,
    method = test_method,
    symnum.args = list(
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols = c("***", "**", "*", "ns")
    ), label = "p.signif"
  ) + theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) + labs(fill = leg)
  return(p)
}

fig2d <- Muti_Boxplot(
  dat = tumor.pathway.score.group1[,
    as.character(unique(tumor.pathway$OG.TSG))
  ],
  group = tumor.pathway.score.group1$copykat.pred,
  group_cols = ggsci::pal_lancet()(9)[c(2, 1)],
  test_method = "wilcox.test",
  leg = "CAF_0", ylab = "GSVA Score"
)
ggsave("results/Fig2c.pdf", fig2d, height = 5, width = 9)

#C1
tumor.pathway.score.group2 <- tumor.score.copy[
  tumor.score.copy$seurat_clusters == "CAF_1",
]
fig2e <- Muti_Boxplot(
  dat = tumor.pathway.score.group2[,
    as.character(unique(tumor.pathway$OG.TSG))
  ],
  group = tumor.pathway.score.group2$copykat.pred,
  group_cols = ggsci::pal_lancet()(9)[c(2,1)],
  test_method = "wilcox.test",
  leg = "CAF_1",
  ylab = "GSVA Score"
)
ggsave("results/Fig2d.pdf", fig2e, height = 5, width = 9)

#C2
tumor.pathway.score.group3 <- tumor.score.copy[
  tumor.score.copy$seurat_clusters == "CAF_2",
]

fig2f <- Muti_Boxplot(
  dat = tumor.pathway.score.group3[,
    as.character(unique(tumor.pathway$OG.TSG))
  ],
  group = tumor.pathway.score.group3$copykat.pred,
  group_cols = ggsci::pal_lancet()(9)[c(2,1)],
  test_method = "wilcox.test",
  leg = "CAF_2",
  ylab = "GSVA Score"
)

ggsave("results/Fig2e.pdf", fig2f, height = 5, width = 9)

#C3
tumor.pathway.score.group4 <- tumor.score.copy[
  tumor.score.copy$seurat_clusters == "CAF_3",
]
fig2g <- Muti_Boxplot(
  dat = tumor.pathway.score.group4[,
    as.character(unique(tumor.pathway$OG.TSG))
  ],
  group = tumor.pathway.score.group4$copykat.pred,
  group_cols = ggsci::pal_lancet()(9)[c(2, 1)],
  test_method = "wilcox.test",
  leg = "CAF_3", ylab = "GSVA Score"
)

ggsave("results/Fig2f.pdf", fig2g, height = 5, width = 9)

#C4
tumor.pathway.score.group5 <- tumor.score.copy[
  tumor.score.copy$seurat_clusters == "CAF_4",
]
fig2h <- Muti_Boxplot(
  dat = tumor.pathway.score.group5[,
    as.character(unique(tumor.pathway$OG.TSG))
  ],
  group = tumor.pathway.score.group5$copykat.pred,
  group_cols = ggsci::pal_lancet()(9)[c(2, 1)],
  test_method = "wilcox.test",
  leg = "CAF_4", ylab = "GSVA Score"
)

ggsave("results/Fig2g.pdf", fig2h, height = 5, width = 9)
save.image("all.RData")