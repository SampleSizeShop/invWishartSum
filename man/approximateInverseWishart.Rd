\name{approximateInverseWishart}
\alias{approximateInverseWishart}
\title{approximateInverseWishart}
\usage{
approximateInverseWishart(invWishartList, method = "trace",
  scaleMatrixList = NULL, cell = NULL)
}
\arguments{
  \item{invWishartList}{the list of inverse Wishart
  matrices in the sum}

  \item{method}{the moment-matching approach to use.  Valid
  values are \code{trace}, \code{log-determinant}, and
  \code{cell}}

  \item{scaleMatrixList}{optional, list of scale matrices
  to form a sum of quadratic forms. There must be one scale
  matrix for each inverse Wishart in the
  \code{invWishartList}.}

  \item{cell}{optional, a vector containing the row and
  column of the cell used for variance matching when using
  \code{method='cell'}}
}
\value{
the approximating inverse Wishart object
}
\description{
Approximate the distribution of a sum of inverse Wishart
matrices or a sum of quadratic forms in inverse Wishart
matrices with a single inverse Wishart.
}
\details{
Three approximation methods are currently supported
\enumerate{ \item{\code{trace}: }{The default method which
matches the expectation of the sum and the variance of the
trace of the sum} \item{\code{logDeterminant: }}{A method
which matches the expectation of the sum and the
expectation of the log determinant of the sum}
\item{\code{cell}: }{A method which matches the expectation
of the sum and the variance of a specified cell of the sum}
}

For quadratic forms, a list of scale matrices should be
specified which will be pre- and post- multiplied onto each
inverse Wishart.  For example, for the inputs
\code{invWishartList = c(X1, X2)} and \code{scaleMatrixList
= c(A,B)}, the function will calculate the approximate
distribution of \code{A'X1A + B'X2B}
}
\references{
The \code{trace} method implements the approach of:\cr
Kreidler, S. M., Muller, K. E., & Glueck, D. H. An
Approximation to the Distribution of the Sum of Inverse
Wishart Matrices, In review.

The \code{logDeterminant} method is based on the approach
of:\cr Granstrom, K., & Orguner, U. (2012). On the
reduction of Gaussian inverse Wishart mixtures. In 2012
15th International Conference on Information Fusion
(FUSION) (pp. 2162-2169).
}
\seealso{
\code{\link{inverseWishart}}
}

