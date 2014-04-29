\name{plotCovarianceFromFile}
\alias{plotCovarianceFromFile}
\title{plotCovarianceFromFile}
\usage{
plotCovarianceFromFile(sumFile, approxFile, cell1 = c(1, 1), cell2 = c(1,
  2), height = 4, width = 4, filename, lims = c(-2, 2, -2, 2),
  col = c("blue", "black"))
}
\arguments{
  \item{sumFile}{file containing a sample from the
  distribution of a sum of inverse Wishart matrices}

  \item{approxFile}{file containing a sample from a single
  inverse Wishart which approximates the distribution of
  the sum of inverse Wishart matrices.}

  \item{cell1}{the row and column of the first cell in the
  bivariate density}

  \item{cell2}{the row and column of the second cell in the
  bivariate density}

  \item{height}{height of the output pdf figure}

  \item{width}{width of the output pdf figure}

  \item{filename}{filename for the output pdf plot}

  \item{lims}{plot limits}

  \item{col}{optional list of colors for the empirical and
  approximate densities}
}
\value{
comparison plot of bivariate densities
}
\description{
Loads samples from files for the empirical sum of inverse
Wishart matrices and a single inverse Wishart which
approximates the distribution of the sum of inverse Wishart
matrices. Creates a comparison plot showing the bivariate
density of two cells.
}

