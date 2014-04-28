\name{approximateInverseWishart.cell}
\alias{approximateInverseWishart.cell}
\title{approximateInverseWishart.cell}
\usage{
approximateInverseWishart.cell(invWishartList, row = 1, column = 1)
}
\arguments{
  \item{invWishartList}{the list of inverse Wishart
  matrices in the sum}

  \item{row}{the row of the cell used in the variance
  match}

  \item{column}{the column of the cell used in the variance
  match}
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
sum and the variance of the specified cell of the sum.
}
\note{
This approach has lower accuracy that the
\code{\link{approximateInverseWishart.trace}} function and
is not recommended.
}
\seealso{
\code{\link{inverseWishart}}
}

