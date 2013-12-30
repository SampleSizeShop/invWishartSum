% File src/library/base/man/load.Rd
\name{approximateInverseWishart}
\alias{approximateInverseWishart}
\title{Obtain a single inverse Wishart distribution which approximates the distribution of the sum of inverse Wisharts}
\description{
  The distribution of a sum of arbitrary inverse Wishart matrices is not known.  This function
  calculates a single inverse Wishart distribution which approximates the distribution of the sum.
}
\usage{
approximateInverseWishart(invWishartList, method="trace", scaleMatrixList=NULL, cell=NULL)
}
\arguments{
  \item{invWishartList}{a list of inverse Wishart distributions of equal dimension.}
  \item{method}{the approximation method.  Supported values are "trace", "logDeterminant", or "cellVariance."}
  \item{scaleMatrixList}{an optional list of matrices by which to scale each inverse Wishart.}
}
\seealso{
  \code{\link{wishart}, \link{inverseWishart}}.
}
\examples{
## coming soon
}
