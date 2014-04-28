\name{approximateInverseWishart.trace}
\alias{approximateInverseWishart.trace}
\title{approximateInverseWishart.trace}
\usage{
approximateInverseWishart.trace(invWishartList)
}
\arguments{
  \item{invWishartList}{the list of inverse Wishart
  matrices in the sum}
}
\value{
the approximating inverse Wishart object
}
\description{
Approximate the distribution of a sum of inverse Wishart
matrices with a single inverse Wishart.
}
\details{
The two-moment approximation matches the expectation of the
sum and the variance of the trace of the sum.
}
\references{
Kreidler, S. M., Muller, K. E., & Glueck, D. H. An
Approximation to the Distribution of the Sum of Inverse
Wishart Matrices, In review.
}
\seealso{
\code{\link{inverseWishart}}
}

