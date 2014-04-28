\name{plotWishartElement}
\alias{plotWishartElement}
\title{plotWishartElement}
\usage{
plotWishartElement(empirical, approximate, row, column, col = c("blue",
  "black"), lty = c(1, 1), ylim = c(0, 1), xlim = c(0, 1))
}
\arguments{
  \item{empirical}{a sample from the distribution of a sum
  of inverse Wishart matrices}

  \item{approximate}{a sample from a single inverse Wishart
  which approximates the distribution of the sum of inverse
  Wishart matrices.}

  \item{row}{cell row}

  \item{column}{cell column}

  \item{col}{optional list of colors for the empirical and
  approximate densities}

  \item{lty}{optional list of line styles for the empirical
  and approximate densities}

  \item{ylim}{optional Y-axis plot limits}

  \item{xlim}{optional X-axis plot limits}
}
\value{
plot of empirical and approximate denisities
}
\description{
Plot the density of a single cell of the empirical inverse
Wishart sum and the approximating inverse Wishart density
}

