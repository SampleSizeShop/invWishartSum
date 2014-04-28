\name{operand.logDeterminant}
\alias{operand.logDeterminant}
\title{operand.logDeterminant}
\usage{
operand.logDeterminant(N, dim, weightedPrecisionSum, expectationSum)
}
\arguments{
  \item{N}{the current degrees of freedom}

  \item{dim}{the dimension of the inverse Wishart matrices
  in the sum}

  \item{weightedPrecisionSum}{weighted sum of precision
  matrices}

  \item{expectationSum}{sum of the expected values of each
  inverse Wishart}
}
\value{
will return 0 or near 0 at optimal value of N
}
\description{
Equation used for optimization when matching the
expectation and log determinant to approximate the
distribution of a sum of inverse Wishart matrices.
Optimization occurs across the degrees of freedom, N.
}
\note{
Based on the approach of: Granstrom, K., & Orguner, U.
(2012). On the reduction of Gaussian inverse Wishart
mixtures. In 2012 15th International Conference on
Information Fusion (FUSION) (pp. 2162-2169).
}
\keyword{internal}

