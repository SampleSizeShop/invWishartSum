\name{compare.covarDiff}
\alias{compare.covarDiff}
\title{compare.covarDiff}
\usage{
compare.covarDiff(empiricalReps, approxReps, cell1 = c(1, 2), cell2 = c(1,
  3), style = "contour", lims = c(-2, 2, -2, 2), col = c("red", "black"))
}
\arguments{
  \item{empiricalReps}{a sample from the distribution of a
  sum of inverse Wishart matrices}

  \item{approxReps}{a sample from a single inverse Wishart
  which approximates the distribution of the sum of inverse
  Wishart matrices.}

  \item{cell1}{the row and column of the first cell in the
  bivariate density}

  \item{cell2}{the row and column of the second cell in the
  bivariate density}

  \item{style}{the type of 3D plot, either
  \code{'contour'}, \code{'image'}, or \code{'persp'}}

  \item{ylim}{optional Y-axis plot limits}

  \item{xlim}{optional X-axis plot limits}

  \item{col}{optional list of colors for the empirical and
  approximate densities}
}
\value{
comparison plot of the difference between the bivariate
densities
}
\description{
Creates a comparison plot showing the difference in the
bivariate densities of two cells of the empirical sum of
inverse Wishart matrices and a sample from a single inverse
Wishart which approximates the distribution of the sum of
inverse Wishart matrices.
}

