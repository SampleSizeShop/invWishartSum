\name{operand.trace}
\alias{operand.trace}
\title{operand.trace}
\usage{
operand.trace(N, dim, g1, g2, g3, g4)
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
expectation and variance of the trace to approximate the
distribution of a sum of inverse Wishart matrices.
Optimization occurs across the degrees of freedom, N.
}
\note{
THIS FUNCTION IS OBSOLETE SINCE THIS APPROACH HAS A
CLOSED-FORM SOLUTION
}
\keyword{internal}

