\name{runSimulationStudy}
\alias{runSimulationStudy}
\title{runSimulationStudy}
\usage{
runSimulationStudy(study.seed = 1812, study.data.dir = ".",
  study.figures.dir = ".")
}
\arguments{
  \item{study.seed}{the random number seed (defaults to
  1812)}

  \item{study.data.dir}{the directory into which data files
  are written (defaults to current working directory)}

  \item{study.figures.dir}{the directory into which pdf
  figures are written (defaults to the current working
  directory)}
}
\description{
This function reproduces the simulation study results for
the manuscript:\cr Kreidler, S. M., Muller, K. E., &
Glueck, D. H. An Approximation to the Distribution of the
Sum of Inverse Wishart Matrices, In review.
}
\note{
The energy distance calculations will take several minutes
to a few hours to run depending on processor speed
}

