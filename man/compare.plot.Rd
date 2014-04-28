\name{compare.plot}
\alias{compare.plot}
\title{compare.plot}
\usage{
compare.plot(sumReplicates, approxReplicates, ylim = c(0, 1), xlim = c(0,
  1), col = c("red", "black"))
}
\arguments{
  \item{sumReplicates}{a sample from the distribution of a
  sum of inverse Wishart matrices}

  \item{approxReplicates}{a sample from a single inverse
  Wishart which approximates the distribution of the sum of
  inverse Wishart matrices.}

  \item{ylim}{optional Y-axis plot limits}

  \item{xlim}{optional X-axis plot limits}

  \item{col}{optional list of colors for the empirical and
  approximate densities}
}
\value{
plot of univariate densities
}
\description{
Creates a comparison plot showing the univariate density of
each cell of the empirical sum of inverse Wishart matrices
and a sample from a single inverse Wishart which
approximates the distribution of the sum of inverse Wishart
matrices.
}

