\name{approximateInverseWishartScaled.trace}
\alias{approximateInverseWishartScaled.trace}
\title{approximateInverseWishartScaled.trace}
\usage{
approximateInverseWishartScaled.trace(invWishartList, scaleMatrixList)
}
\arguments{
  \item{invWishartList}{the list of inverse Wishart
  matrices in the sum}

  \item{scaleMatrixList}{the list of matrices which will be
  pre- and post-multiplied onto each corresponding inverse
  Wishart to build the quadratic forms}
}
\value{
the approximating inverse Wishart object
}
\description{
Approximate the distribution of a sum of quadratic forms in
inverse Wishart matrices with a single inverse Wishart. The
two-moment approximation matches the expectation of the sum
and the variance of the trace of the sum.
}
\references{
Kreidler, S. M., Muller, K. E., & Glueck, D. H. An
Approximation to the Distribution of the Sum of Inverse
Wishart Matrices, In review.
}
\seealso{
\code{\link{inverseWishart}}
}

