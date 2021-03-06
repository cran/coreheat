% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cormap2.R
\name{cormap_filt}
\alias{cormap_filt}
\title{Automatically split clusters based on noise level and hierarchy}
\usage{
cormap_filt(
  x,
  na.frac = 0.1,
  method = "ward.D",
  do.abs = TRUE,
  main = "correlation map",
  postfix = NULL,
  p.thr = 0.01,
  cex = 0.2,
  cex.clust = cex,
  cex.filt = cex,
  cut.thr = NULL,
  cor.thr = NULL,
  cor.cluster = 1,
  cor.window = NULL,
  do.plots = c("dend", "full.heat", "filt.heat"),
  genes2highl = NULL,
  order.list = TRUE,
  convert = TRUE,
  biomart = FALSE,
  add.sig = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{x}{(\code{ExpressionSet}, \code{data.frame} or \code{numeric}). A numeric data frame, matrix or an ExpressionSet object.}

\item{na.frac}{(\code{numeric}). Fraction of missing values allowed per row of the input matrix. Defaults to \code{0.1} which
means LESS than 10 per cent of the values in one row are allowed to be NAs.}

\item{method}{(\code{character}). The agglomeration method used for clustering. See help for \code{\link[stats]{hclust}}.
Defaults to "ward.D".}

\item{do.abs}{(\code{logical}). Should the distances for clustering be calculated based on the absolute correlation values?
In other words, should the sign of the correlation be ignored in favor of its strength?}

\item{main}{(\code{character}). The main title of the plot. Defaults to \code{""}.}

\item{postfix}{(\code{character} of \code{logical}). A plot sub-title. Will be printed below the main title. Defaults to \code{NULL}.}

\item{p.thr}{(\code{numeric}). P-Value threshold for filtering sub-clusterd with significant correlations. Defaults to \code{0.01}.}

\item{cex}{(\code{numeric}). Font size for the heatmap of the unfiltered correlation matrix. Defaults to \code{0.2}.}

\item{cex.clust}{(\code{numeric}). Font size for the dendrogram plot of the unfiltered correlation matrix clusters.
Defaults to \code{cex}.}

\item{cex.filt}{(\code{numeric}). Font size for the heatmap of the filtered correlation matrix. Defaults to \code{cex}.}

\item{cut.thr}{(\code{numeric}). Threshold at which dendrogram branches are to be cut. Passed on to argument \code{h} in
\code{\link[stats]{cut.dendrogram}}. Defaults to \code{NULL} meaning no cutting.}

\item{cor.thr}{(\code{numeric}). Correlation threshold to filter the correlation matrix for plotting. Defaults to \code{NULL} meaning
no filtering. Note that this value will be applied to margin \option{cor.mar} of the values per row.}

\item{cor.cluster}{(\code{numeric}). The correlation cluster along the diagonal 'line' in the heatmap that should be
zoomed into. A sliding window of size \code{cor.window} will be moved along the diagonal of the correlation
matrix to find the cluster with the most corelation values meeting \code{core.thr}. Defaults to \code{1}.}

\item{cor.window}{(\code{numeric}). The size of the sliding window (see \code{cor.cluster}). Defaults to \code{NULL}.
Note that this works only for positive correlations.}

\item{do.plots}{(\code{character}). The plots to be produced. A character vector containing one or more of "dend"
to produce the dendrogram plot, "full.heat" to produce the heatmap of the unfiltered correlation matrix, and
"filt.heat" to produce the heatmap of the filtered correlation matrix. Defaults to all three plots.}

\item{genes2highl}{(\code{character}). Vector of gene symbols (or whatever labels are used) to be highlighted.
If not \code{NULL} will draw a semi-transparent rectangle around the labels and rows or columns in the heatmap
labels.}

\item{order.list}{(\code{logical}). Should the order of the correlation matrix, i.e. the 'list' of labels be reversed?
Meaningful if the order of input variables should be preserved because \code{\link[graphics]{image}} turns the input
matrix. Defaults to \code{TRUE}.}

\item{convert}{(\code{logical}). Should an attempt be made to convert IDs provided as row names of the input or in \option{lab}?
Defaults to \code{TRUE}. Conversion will be done using BioMart or an annotation package, depending on \option{biomart}.}

\item{biomart}{(\code{logical}). Should BioMart (or an annotation package) be used to convert IDs? If \code{TRUE}
the \code{todisp2} function in package \code{convertid} attempts to access the BioMart API to convert ENSG IDs to Gene Symbols
Defaults to \code{FALSE} which will use the traditional \code{AnnotationDbi} Bimap interface.}

\item{add.sig}{(\code{logical}). Should significance asterisks be drawn? If \code{TRUE} P-Values for correlation significance
are calculated and encoded as asterisks. See 'Details'.}

\item{verbose}{(\code{logical}). Should verbose output be written to the console? Defaults to \code{FALSE}.}
}
\value{
A \code{list}. If the dendrogram is being cut, i.e., \code{cut.thr} is not \code{NULL}, a list of
\tabular{ll}{
\tab clusters: the list of cluster labels from \code{lower} component of the \code{cut.dendrogram} output which
    is list with the branches obtained from cutting the tree\cr
\tab filt: the index of the cluster labels passing the signal metrics threshold\cr
\tab filt_cluster: the list of the filtered cluster labels\cr
\tab h: the cut threshold\cr
\tab p.thr: the P-Value threshold for filtering sub-clusters\cr
\tab metric: the signal metrics for all sub-clusters\cr
\tab cormat: the clustered (ordered) correlation matrix\cr
\tab hclust: a list of hierarchical clustering metrics (output of \code{\link[stats]{hclust}})\cr
\tab pvalues: the correlation P-Value matrix\cr
}

If no tree cutting is applied, a list of
\tabular{ll}{
\tab cormat: the clustered (ordered) correlation matrix\cr
\tab hclust: a list of hierarchical clustering metrics (output of \code{\link[stats]{hclust}})\cr
\tab pvalues: the correlation P-Value matrix\cr
}
}
\description{
\command{cormat_filt} splits (cuts) the dendrogram at a given threshold dividing it into larger or
    smaller "sub-clusters". Correlation P-Values (see \code{\link[coreheat]{eset_cor}}) are converted to represent
    significance as a sub-cluster-wise signal metric used for filtering. Optionally, up to 3 plots are produced,
    the third one being a filtered heatmap based on significance and three height cutting.
}
\details{
P-Values are calculated from the t-test value of the correlation coefficient: \eqn{t = r x sqrt(n-2) / sqrt(1-r^2)},
    where r is the correlation coefficient, n is the number of samples with no missing values for each gene (row-wise
    \code{ncol(eset)} minus the number of columns that have an NA). P-Values are then calculated using \code{\link[stats]{pt}} and
    corrected account for the two-tailed nature of the test, i.e., the possibility of positive as well as negative correlation.
    The approach to calculate correlation significance was adopted from \emph{Miles, J., & Banyard, P. (2007)} on
    "Calculating the exact significance of a Pearson correlation in MS Excel".

    To obtain a suitable metric for isolating significant sub-clusters, P-Values are represented as \eqn{-log10(median(pval))}
    where \code{pval} is the \emph{median of the parallel maximum of all P-Values belonging to the sub-cluster and
    \code{1e-38}} to avoid values of zero (0).
}
