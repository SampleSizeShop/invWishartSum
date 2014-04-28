\docType{methods}
\name{invert}
\alias{invert}
\alias{invert,inverseWishart,inverseWishart-method}
\alias{invert,inverseWishart-method}
\alias{invert,wishart,wishart-method}
\alias{invert,wishart-method}
\title{invert}
\usage{
invert(object)

\S4method{invert}{wishart}(object)

\S4method{invert}{inverseWishart}(object)
}
\arguments{
  \item{a}{\code{\link{wishart}} or
  \code{\link{inverseWishart}} distribution object}
}
\value{
the distribution object for the inverse of the specified
Wishart or inverse Wishart
}
\description{
Form the distribution of the inverse of the specified
Wishart or inverse Wishart.
}
\note{
Based on Theorem 3.4.1, p 111 from

Gupta, A. K., & Nagar, D. K. (2000). Matrix variate
distributions. Boca Raton, FL: Chapman & Hall.
}
\seealso{
\code{\link{wishart}} and \code{\link{inverseWishart}}
}

