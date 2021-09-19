## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(coreheat)

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("coreheat")

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("vfey/coreheat")

## ---- out.width='85%', fig.width=6, fig.height=6, fig.align='center'----------
set.seed(1234)
mat <- matrix(c(rnorm(100, mean = 1), rnorm(100, mean = -1)), nrow = 20)
rownames(mat) <- paste0("gene-", 1:20)
colnames(mat) <- paste0(c("A", "B"), rep(1:5, 2))
cormap2(mat, convert = FALSE, main="Random matrix")

## ---- out.width='85%', fig.width=6, fig.height=6, fig.align='center'----------
## Read data and prepare input data frame
fl <- system.file("extdata", "PrCaTCGASample.txt", package = "coreheat", mustWork = TRUE)
dat0 <- read.delim(fl, stringsAsFactors=FALSE)
dat1 <- data.frame(dat0[, grep("TCGA", names(dat0))], row.names=dat0$ensembl_gene_id)
cormap2(dat1, main="TCGA data frame + ID conversion")

## ---- out.width='85%', fig.width=6, fig.height=6, fig.align='center'----------
dat2 <- as.matrix(dat0[, grep("TCGA", names(dat0))])
sym <- dat0$hgnc_symbol
cormap2(dat1, convert=FALSE, lab=sym, genes2highl=c("GNAS","NCOR1","AR", "ATM"),
main="TCGA matrix + custom labels")

## ---- out.width='85%', fig.width=10, fig.height=10,  fig.align='center'-------
expr <- Biobase::ExpressionSet(as.matrix(dat1))
cormap2(expr, add.sig=TRUE, autoadj=FALSE, cex=1, main="TCGA ExpressionSet object + ID conversion")

## ---- out.width='85%', fig.width=6, fig.height=6, fig.align='center'----------
fl <- system.file("extdata", "PrCaTCGASample_500.txt", package = "coreheat", mustWork = TRUE)
dat0_500 <- read.delim(fl, stringsAsFactors=FALSE)
dat1_500 <- data.frame(dat0_500[, grep("TCGA", names(dat0_500))], row.names=dat0_500$ensembl_gene_id)
cormap2(dat1_500, biomart = TRUE, main="TCGA 500 genes")

## ---- out.width='85%', fig.width=7, fig.height=7, fig.align='center'----------
cm <- cormap2(dat1_500, cut.thr = 1.5, cut.size = 40, biomart = TRUE,
              genes2highl = c("RPLP1", "NOP53", "RPL13", "RPS17", "RPL11", "RPS24", "RPL30", "EEF1D", "RPL27",
                              "SNHG29", "HSPA1A", "HSPA1B", "DNAJB1", "BHLHE40", "NR4A1", "JUND", "HLA-A", "HLA-B",
                              "BTF3", "EEF1A1P9"),
              main="TCGA 500 genes", verbose = TRUE)

## ---- out.width='85%', fig.width=6, fig.height=6, fig.align='center'----------
cm2 <- cormap2(cormat = cm, cut.thr = 1.2, cut.size = 20, biomart = TRUE,
               genes2highl = c("RPLP1", "NOP53", "RPL13", "RPS17", "RPL11", "RPS24", "RPL30", "EEF1D", "RPL27",
                               "SNHG29", "HSPA1A", "HSPA1B", "DNAJB1", "BHLHE40", "NR4A1", "JUND", "HLA-A", "HLA-B",
                               "BTF3", "EEF1A1P9"),
               main="TCGA 500 genes", cex = 0.9, verbose = TRUE)

## ---- out.width='85%', fig.width=6, fig.height=6, fig.align='center'----------
cm3 <- cormap2(cormat = cm, cor.thr = 0.6, cor.window = 5, cor.cluster = 1, biomart = TRUE,
               genes2highl = c("RPLP1", "NOP53", "RPL13", "RPS17", "RPL11", "RPS24", "RPL30", "EEF1D", "RPL27",
                               "SNHG29", "HSPA1A", "HSPA1B", "DNAJB1", "BHLHE40", "NR4A1", "JUND", "HLA-A", "HLA-B",
                               "BTF3", "EEF1A1P9"),
               add.sig = TRUE, cex = 2, main="TCGA 500 genes", verbose = TRUE)

## ---- out.width='85%', fig.width=6, fig.height=6, fig.align='center'----------
cm3 <- cormap2(cormat = cm, cor.thr = 0.6, cor.mar = 0.3, biomart = TRUE,
               genes2highl = c("RPLP1", "NOP53", "RPL13", "RPS17", "RPL11", "RPS24", "RPL30", "EEF1D", "RPL27",
                               "SNHG29", "HSPA1A", "HSPA1B", "DNAJB1", "BHLHE40", "NR4A1", "JUND", "HLA-A",
                               "HLA-B"),
               main="TCGA 500 genes", cex = 1, verbose = T)

## ---- out.width='85%', fig.width=8, fig.height=8, fig.align='center'----------
cm3 <- cormap2(cormat = cm, cor.thr = -0.2, cor.mar = 0.2, biomart = F,
               genes2highl = c("RPLP1", "NOP53", "RPL13", "RPS17", "RPL11", "RPS24", "RPL30", "EEF1D", "RPL27",
                               "SNHG29", "HSPA1A", "HSPA1B", "DNAJB1", "BHLHE40", "NR4A1", "JUND", "HLA-A",
                               "HLA-B"),
               add.sig = TRUE, cex = 2, main="TCGA 500 genes")

## ---- out.width='85%', fig.width=8, fig.height=8, fig.align='center'----------
cormap_filt(dat1_500, main = "TCGA 500 genes", cut.thr = 21, p.thr = 0.01, cex.filt = 0.8)

