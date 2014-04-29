\name{getEnergyDistanceForFilelist}
\alias{getEnergyDistanceForFilelist}
\title{getEnergyDistanceForFilelist}
\usage{
getEnergyDistanceForFilelist(sumFileList, approxFileList)
}
\arguments{
  \item{sumFileList}{list of files containing the empirical
  samples}

  \item{approxFileList}{list of files containing the
  samples for the approximating inverse Wishart}
}
\value{
the list of energy distance values
}
\description{
Calculate the energy distance between samples contained in
the specified .Rdata files containing samples for the sum
of inverse Wisharts and samples for the approximating
inverse Wishart.
}
\note{
This function assumes that the filename lists are in the
same order to pair the samples correctly
}

