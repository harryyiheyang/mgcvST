\name{mgcvST.wgcna}
\alias{mgcvST.wgcna}
\alias{print.mgcvST_wgcna}
\title{Identify co-expression modules within explicitly selected gene blocks}
\usage{
mgcvST.wgcna(fitmgcvST, indices, group = NULL,
  wgcna.para = NULL, verbose = FALSE)

\method{print}{mgcvST_wgcna}(x, ...)
}
\arguments{
\item{fitmgcvST}{A compact fit returned by \code{mgcvST.estimate()} or
\code{inlaST.estimate()}.}
\item{indices}{Required gene IDs, integer positions in
\code{fitmgcvST$feature_id}, or a named list of such vectors. Each block must
contain at least two distinct genes. A vector defines the block named
\code{selected}. Gene order is preserved.}
\item{group}{Spatial score-component names. \code{NULL} uses the sole available
component; it is an error when the fit contains multiple score components.}
\item{wgcna.para}{\code{NULL}, an empty list, or named partial overrides of
\code{networkType = "signed"}, \code{power = 6}, \code{TOMType = "signed"},
\code{hclustMethod = "average"}, \code{minClusterSize = 20L}, and
\code{deepSplit = 1L}. Unspecified settings keep these defaults; unknown
names are errors. Blocks smaller than \code{minClusterSize} retain their
matrices and receive grey (zero) labels, without changing the requested size.}
\item{verbose}{Whether to display compact progress messages.}
\item{x}{An \code{mgcvST_wgcna} result.}
\item{...}{Unused.}
}
\value{
An \code{mgcvST_wgcna} object with \code{modules} (component, feature ID,
integer module, color), named \code{networks} (feature IDs, covariance,
correlation, adjacency, TOM, tree, labels, modules, coordinate count \code{q},
and status), \code{score} (aligned \code{A}, group names and widths),
\code{settings}, and \code{timing}. Zero labels mean unassigned (grey).
Module labels are local to each input block. The print method returns
\code{x} invisibly.
}
\description{
Constructs the original score covariance \code{crossprod(A) / nrow(A)} from
an \code{mgcvST.estimate()} or \code{inlaST.estimate()} fit, converts it to
correlation, and performs WGCNA. \code{indices} selects genes, for example
the genes in a connected component identified after \code{mgcvST.test()}.
}
\details{
A named list analyzes several gene blocks separately. The function never
pools blocks, selects significant pairs, or reruns estimation or testing.
All pairs within each selected block enter its covariance.

Shared SPDE factors are constructed once. All backends compute only the score
vectors needed here; pair-test calibration matrices are unnecessary. The
result retains these vectors and the matrices needed for downstream analysis
without modifying the supplied fit. No centering, ridge, or covariance
projection is applied.

\code{group} selects spatial score components, not connected gene blocks.
With multiple components, select their names explicitly. The selected
coordinate groups are concatenated and divided by their total coordinate
count. For a legacy single-SPDE compact fit this group is named \code{global}.

The WGCNA, dynamicTreeCut, and fastcluster packages are required only when calling this
function. Biological annotations and module eigengenes are not inferred
from the covariance input.
}
\examples{
\dontrun{
W <- mgcvST.wgcna(fit, indices = genes)
W <- mgcvST.wgcna(fit,
  indices = list(block1 = genes1, block2 = genes2),
  wgcna.para = list(power = 6))
W$modules
W$networks$block1$covariance
}
}
\seealso{\code{\link{mgcvST.estimate}}, \code{\link{inlaST.estimate}},
\code{\link{mgcvST.test}}}
