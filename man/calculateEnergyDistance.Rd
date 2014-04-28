\name{calculateEnergyDistance}
\alias{calculateEnergyDistance}
\title{calculateEnergyDistance}
\usage{
calculateEnergyDistance(sample1, sample2)
}
\arguments{
  \item{sample1}{a list of matrix samples from the first
  matrix variate distribution}

  \item{sample1}{a list of matrix samples from the second
  matrix variate distribution}
}
\value{
the energy distance
}
\description{
Estimate the energy distance between two samples from
matrix variate distributions of equal dimension. Assumes
equal sample sizes for each distribution.
}
\note{
This function is computationally expensive and may take
several minutes or longer to run, depending on the size of
each sample.
}
\references{
Szekely, G. J., & Rizzo, M. L. (2004). Testing for Equal
Distributions in High Dimension. InterStat, 5.
}

