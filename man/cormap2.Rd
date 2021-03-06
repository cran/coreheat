% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cormap2.R
\name{cormap2}
\alias{cormap2}
\title{Draw correlation maps from large datasets.}
\usage{
cormap2(
  x,
  cormat = NULL,
  lab = NULL,
  convert = TRUE,
  biomart = FALSE,
  cluster_correlations = TRUE,
  main = "",
  postfix = NULL,
  cex = NULL,
  na.frac = 0.1,
  cor.cluster = 1,
  cor.window = NULL,
  cor.thr = NULL,
  cor.mar = 0.5,
  cut.thr = NULL,
  cut.size = 5,
  autoadj = TRUE,
  labelheight = NULL,
  labelwidth = NULL,
  add.sig = FALSE,
  genes2highl = NULL,
  order.list = TRUE,
  doPlot = TRUE,
  updateProgress = NULL,
  verbose = FALSE
)
}
\arguments{
\item{x}{(\code{ExpressionSet}, \code{data.frame} or \code{numeric}). A numeric data frame, matrix or an ExpressionSet object.}

\item{cormat}{(\code{numeric}). A correlation matrix. If this not \code{NULL} then \option{x} is ignored. Defaults to \code{NULL}.}

\item{lab}{(\code{character}). Optional row/column labels for the heatmap. Defaults to NULL meaning the row names of the input data
are used. Note that the order of the labels must match the order of the row names of the input data!}

\item{convert}{(\code{logical}). Should an attempt be made to convert IDs provided as row names of the input or in \option{lab}?
Defaults to \code{TRUE}. Conversion will be done using BioMart or an annotation package, depending on \option{biomart}.}

\item{biomart}{(\code{logical}). Should BioMart (or an annotation package) be used to convert IDs? If \code{TRUE}
the \code{todisp2} function in package \code{convertid} attempts to access the BioMart API to convert ENSG IDs to Gene Symbols
Defaults to \code{FALSE} which will use the traditional \code{AnnotationDbi} Bimap interface.}

\item{cluster_correlations}{(\code{logical}). Should the correlation matrix be clustered before plotting? Defaults to \code{TRUE}.}

\item{main}{(\code{character}). The main title of the plot. Defaults to \code{""}.}

\item{postfix}{(\code{character} of \code{logical}). A plot sub-title. Will be printed below the main title. Defaults to \code{NULL}.}

\item{cex}{(\code{numeric}). Font size. Defaults to \code{0.5} if \code{autoadj} is \code{FALSE}. See 'Details'.}

\item{na.frac}{(\code{numeric}). Fraction of missing values allowed per row of the input matrix. Defaults to \code{0.1} which
means LESS than 10 per cent of the values in one row are allowed to be NAs.}

\item{cor.cluster}{(\code{numeric}). The correlation cluster along the diagonal 'line' in the heatmap that should be
zoomed into. A sliding window of size \code{cor.window} will be moved along the diagonal of the correlation
matrix to find the cluster with the most corelation values meeting \code{core.thr}. Defaults to \code{1}.}

\item{cor.window}{(\code{numeric}). The size of the sliding window (see \code{cor.cluster}). Defaults to \code{NULL}.
Note that this works only for positive correlations.}

\item{cor.thr}{(\code{numeric}). Correlation threshold to filter the correlation matrix for plotting. Defaults to \code{NULL} meaning
no filtering. Note that this value will be applied to margin \option{cor.mar} of the values per row.}

\item{cor.mar}{(\code{numeric}). Margin of the values per row of the correlation matrix the \option{cor.thr} filter needs to
meet. Defaults to \code{0.5} meaning at least 50 per cent of the values in a row need to meet the threshold in order to keep the row.}

\item{cut.thr}{(\code{numeric}). Threshold at which dendrogram branches are to be cut. Passed on to argument \code{cutHeight} in
\code{\link[WGCNA]{cutreeStatic}}. Defaults to \code{NULL} meaning no cutting.}

\item{cut.size}{(\code{numeric}). Minimum number of objects on a dendrogram branch considered a cluster. Passed on to argument
\code{minSize} in \code{\link[WGCNA]{cutreeStatic}}. Defaults to \code{5}.}

\item{autoadj}{(\code{logical}). Should plot measures be adjusted automatically? Defaults to \code{TRUE}.}

\item{labelheight}{(\code{numeric} or \code{lcm(numeric)}). Relative or absolute height (using \code{\link[graphics]{lcm}},
see \code{\link[graphics]{layout}}) of the labels. Defaults to \code{0.2} if \code{autoadj} is \code{FALSE}. See 'Details'.}

\item{labelwidth}{(\code{numeric} or \code{lcm(numeric)}). Relative or absolute width (using \code{\link[graphics]{lcm}},
see \code{\link[graphics]{layout}}) of the labels. Defaults to \code{0.2} if \code{autoadj} is \code{FALSE}. See 'Details'.}

\item{add.sig}{(\code{logical}). Should significance asterisks be drawn? If \code{TRUE} P-Values for correlation significance
are calculated and encoded as asterisks. See 'Details'.}

\item{genes2highl}{(\code{character}). Vector of gene symbols (or whatever labels are used) to be highlighted.
If not \code{NULL} will draw a semi-transparent rectangle around the labels and rows or columns in the heatmap
labels.}

\item{order.list}{(\code{logical}). Should the order of the correlation matrix, i.e. the 'list' of labels be reversed?
Meaningful if the order of input variables should be preserved because \code{\link[graphics]{image}} turns the input
matrix. Defaults to \code{TRUE}.}

\item{doPlot}{(\code{logical}). Draw the plot? Defaults to \code{TRUE}.}

\item{updateProgress}{(\code{function}). Function for updating a progress bar in a Shiny web application. This was added here
for the \strong{BioCPR} application.}

\item{verbose}{(\code{logical}). Should verbose output be written to the console? Defaults to \code{FALSE}.}
}
\value{
Invisibly returns the correlation matrix, though the function is mainly called for its side-effect of producing
    a heatmap (if \code{doPlot = TRUE} which is the default).
}
\description{
\command{cormap2()} generates pair-wise correlations from an input ExpressionSet object, a \code{data.frame} or a
    numerical \code{matrix}. With the default options it also produces a heatmap.
}
\details{
P-Values are calculated from the t-test value of the correlation coefficient: \eqn{t = r x sqrt(n-2) / sqrt(1-r^2)},
    where r is the correlation coefficient, n is the number of samples with no missing values for each gene (row-wise
    \code{ncol(eset)} minus the number of columns that have an NA). P-Values are the calculated using \code{\link[stats]{pt}} and
    corrected account for the two-tailed nature of the test, i.e., the possibility of positive as well as negative correlation.
    The approach to calculate correlation significance was adopted from \emph{Miles, J., & Banyard, P. (2007)} on
    "Calculating the exact significance of a Pearson correlation in MS Excel".

    The asterisks encode significance as follows:
    \tabular{ll}{
    \tab P < 0.05:  *\cr
    \tab P < 0.01:  **\cr
    \tab P < 0.001: ***\cr
    }

    The label measures (\code{labelheight}, \code{labelwidth} and \code{cex}) are adjusted automatically by default
    with argument \code{autoadj=TRUE} and have default values which are hard coded into the helper function
    \code{heatmap.cor}. The values calculated by the helper function \code{plotAdjust} can be overridden by setting
    any of those arguments to a valid \code{numeric} or \code{lcm(numeric)} value.
}
\examples{
# 1. Generate a random 10x10 matrix with two distinct sets and plot it with
# default settings without ID conversion since the IDs are made up:
set.seed(1234)
mat <- matrix(c(rnorm(100, mean = 1), rnorm(100, mean = -1)), nrow = 20)
rownames(mat) <- paste0("gene-", 1:20)
colnames(mat) <- paste0(c("A", "B"), rep(1:5, 2))
cormap2(mat, convert=FALSE, main="Random matrix")

# 2. Use a real-world dataset from TCGA (see README file in inst/extdata directory).
# Package 'convertid' is used to convert Ensembl Gene IDs to HGNC Symbols
## Read data and prepare input data frame
fl <- system.file("extdata", "PrCaTCGASample.txt", package = "coreheat", mustWork = TRUE)
dat0 <- read.delim(fl, stringsAsFactors=FALSE)
dat1 <- data.frame(dat0[, grep("TCGA", names(dat0))], row.names=dat0$ensembl_gene_id)
cormap2(dat1, main="TCGA data frame + ID conversion")

# 3. Use separately supplied IDs with a matrix created from the data frame of the
# previous example and highlight genes of interest
dat2 <- as.matrix(dat0[, grep("TCGA", names(dat0))])
sym <- dat0$hgnc_symbol
cormap2(dat1, convert=FALSE, lab=sym, genes2highl=c("GNAS","NCOR1","AR", "ATM"),
main="TCGA matrix + custom labels")

# 4. Use an ExpressionSet object and add significance asterisks
## For simplicity reasons we create the ExpressionSet from a matrix created
## from the data frame in the second example
expr <- Biobase::ExpressionSet(as.matrix(dat1))
cormap2(expr, add.sig=TRUE, main="TCGA ExpressionSet object + ID conversion")

# More examples can be found in the vignette.
}
\references{
Miles, J., & Banyard, P. (2007). \emph{Understanding and using statistics in psychology: A practical introduction.}
    Sage Publications Ltd. \url{https://psycnet.apa.org/record/2007-06525-000}.
}
\seealso{
\code{\link[stats]{pt}}

\code{\link[base]{tcrossprod}}
}
