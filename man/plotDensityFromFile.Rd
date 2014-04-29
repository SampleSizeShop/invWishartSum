\name{plotDensityFromFile}
\alias{plotDensityFromFile}
\title{plotDensityFromFile}
\usage{
plotDensityFromFile(sumFile, approxFile, filename, ylim = c(0, 1),
  xlim = c(0, 1), col = c("blue", "black"))
}
\arguments{
  \item{sumFile}{file containing a sample from the
  distribution of a sum of inverse Wishart matrices}

  \item{approxFile}{file containing a sample from a single
  inverse Wishart which approximates the distribution of
  the sum of inverse Wishart matrices.}

  \item{filename}{filename for the output pdf plot}

  \item{ylim}{optional Y-axis plot limits}

  \item{xlim}{optional X-axis plot limits}

  \item{col}{optional list of colors for the empirical and
  approximate densities}
}
\value{
plot of univariate densities as a pdf
}
\description{
Loads samples from files for the empirical sum of inverse
Wishart matrices and a single inverse Wishart which
approximates the distribution of the sum of inverse Wishart
matrices. Creates a comparison plot showing the univariate
density of each cell.
}

