The invWishartSum Package
=========================

The invWishartSum package for R (>3.0.0) implements an approximation to the
distribution of a sum of inverse Wishart random matrices.  

The package provides companion code for the manuscript:

Kreidler, S. M., Muller, K. E., & Glueck, D. H. An Approximation to 
the Distribution of the Sum of Inverse Wishart Matrices, In review.

### Instructions for replicating the manuscript results 

The results in the above manuscript were produced using R version 3.0.0. To reproduce the results,
perform the following steps:

* Install R version 3.0.x or higher by following the instructions at http://www.r-project.org
* From the R environment, install and load the "devtools" package
```R
> install.packages("devtools")
> library(devtools)
```
* Install the "invWishartSum" package directly from Github.com
```R
> install_github(repo="invWishartSum", user="SampleSizeShop", ref="develop")
```
* Load the library
```R
> library(invWishartSum)
```
* Run the simulation study to produce the manuscript results (takes about 1 hour to run). You may specify output directories for data files (study.data.dir) and figures (study.figures.dir). If omitted, both default to the current working directories.
```R
> runSimulationStudy(study.data.dir="myDataDir", study.figures.dir="myFiguresDir")
```

