\name{generateEmpiricalSumData}
\alias{generateEmpiricalSumData}
\title{generateEmpiricalSumData}
\usage{
generateEmpiricalSumData(invWishartList, scaleMatrixList, dfScaleList,
  replicates = 1000, outputDir = ".", filenamePrefix = "sumReplicates")
}
\arguments{
  \item{invWishartList}{list of inverse Wisharts}

  \item{scaleMatrixList}{list of matrices to pre- and
  post-multiply onto the corresponding inverse Wishart to
  build a quadratic form}

  \item{dfScaleList}{list of scale factors for the degrees
  of freedom for each inverse Wishart in the sum}

  \item{replicates}{the number of samples to generate}

  \item{outputDir}{the output directory for the .Rdata
  files}

  \item{filenamePrefix}{}
}
\value{
Produces one Rdata file for each df scale factor
}
\description{
Generate .Rdata files containing an empirical sample of the
sum of inverse Wisharts, or quadratic forms in inverse
Wisharts
}
\details{
Filenanes for the .Rdata files are automatically generated
and will have the form:\cr
\code{<outputDir>/<filenamePrefix>_dfScale<dfScale-value>.Rdata}
}

