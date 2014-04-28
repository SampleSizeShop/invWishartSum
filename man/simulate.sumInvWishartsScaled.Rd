\name{simulate.sumInvWishartsScaled}
\alias{simulate.sumInvWishartsScaled}
\title{simulate.sumInvWishartsScaled}
\usage{
\method{simulate}{sumInvWishartsScaled}(invWishartList, scaleMatrixList,
  replicates = 1000)
}
\arguments{
  \item{invWishartList}{list of inverse Wisharts}

  \item{scaleMatrixList}{list of matrices to pre- and
  post-multiply onto the corresponding inverse Wishart to
  build a quadratic form}

  \item{replicates}{the number of samples to generate}
}
\value{
A sample of matrices from the distribution of the sum of
the specified quadratic forms
}
\description{
Generate a sample from the distribution of the sum of
quadratic forms in inverse Wishart matrices
}

