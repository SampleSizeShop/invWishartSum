\name{approximateInverseWishart.logDeterminant}
\alias{approximateInverseWishart.logDeterminant}
\title{approximateInverseWishart.logDeterminant}
\usage{
approximateInverseWishart.logDeterminant(invWishartList)
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
matrices with a single inverse Wishart.The two-moment
approximation matches the expectation of the sum and the
expectation of the log determinant of the sum. The solution
to the resulting equations is obtained by numerical
optimization.
}
\references{
Based on the approach of: Granstrom, K., & Orguner, U.
(2012). On the reduction of Gaussian inverse Wishart
mixtures. In 2012 15th International Conference on
Information Fusion (FUSION) (pp. 2162-2169).
}
\seealso{
\code{\link{inverseWishart}}
}

