#####################################################################
#
# Copyright (C) 2014 Sarah M. Kreidler
#
# This file is part of the R package invWishartSum, which provides
# functions to calculate the approximate distribution of a sum
# of inverse central Wishart matrices and certain quadratic forms
# in inverse central Wishart matrices.
# 
# invWishartSum is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# invWishartSum is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with invWishartSum.  If not, see <http://www.gnu.org/licenses/>.
#
#####################################################################

#' generateEmpiricalSumData
#'
#' Generate .Rdata files containing an empirical sample of the sum
#' of inverse Wisharts, or quadratic forms in inverse Wisharts
#' 
#' Filenanes for the .Rdata files are automatically generated and will have the
#' form:\cr
#' \code{<outputDir>/<filenamePrefix>_dfScale<dfScale-value>.Rdata}
#'
#' @param invWishartList list of inverse Wisharts
#' @param scaleMatrixList list of matrices to pre- and post-multiply onto the corresponding
#' inverse Wishart to build a quadratic form
#' @param dfScaleList list of scale factors for the degrees of freedom for each 
#' inverse Wishart in the sum
#' @param replicates the number of samples to generate
#' @param outputDir the output directory for the .Rdata files
#' @param filenamePrefix 
#' @return Produces one Rdata file for each df scale factor
#'
generateEmpiricalSumData <- function(invWishartList, scaleMatrixList, dfScaleList,
                                     replicates=1000,
                                     outputDir=".", filenamePrefix="sumReplicates") {
  filenames=rep("", length(dfScaleList))
  idx = 1
  for(scale in dfScaleList) {

    # scale the degrees of freedom
    dfScaleInvWishartList = lapply(invWishartList, function(invWishart, scale) {
      return (new("inverseWishart", df=(invWishart@df * scale), 
                  precision=invWishart@precision))
    }, scale=scale)
    class(dfScaleInvWishartList) = "inverseWishart"
    
    # get replicates of the sum of the specified inverse Wisharts
    if (!is.null(scaleMatrixList)) {
      # get replicates of the sum of the specified inverse Wisharts
      sumReplicates = simulate.sumInvWishartsScaled(dfScaleInvWishartList, 
                                                    scaleMatrixList,
                                                    replicates=replicates)
    } else {
      # get replicates of the sum of the specified inverse Wisharts
      sumReplicates = simulate.sumInvWisharts(dfScaleInvWishartList, 
                                              replicates=replicates)
    }
    
    # save as an .Rdata file
    filename = paste(c(outputDir, "/", filenamePrefix, "_dfScale", scale, ".Rdata"), collapse="")
    save(sumReplicates, file=filename)
    
    filenames[idx] = filename
    idx = idx + 1
  }
  
  return(filenames);
}

#
# Generate .Rdata files containing a sample of the approximating Wishart
#
#' generateApproximationData
#'
#' Generate .Rdata files containing an empirical sample of a single inverse Wishart
#' used to approximate the distribution of a sum of inverse Wisharts or quadratic forms
#' in inverse Wisharts.  The approximating Wishart is calculating from the input list
#' of inverse Wishart matrices.
#' 
#' Filenanes for the .Rdata files are automatically generated and will have the
#' form:\cr
#' \code{<outputDir>/<filenamePrefix>_dfScale<dfScale-value>.Rdata}
#'
#' @param invWishartList list of inverse Wisharts
#' @param scaleMatrixList list of matrices to pre- and post-multiply onto the corresponding
#' inverse Wishart to build a quadratic form
#' @param dfScaleList list of scale factors for the degrees of freedom for each 
#' inverse Wishart in the sum
#' @param replicates the number of samples to generate
#' @param outputDir the output directory for the .Rdata files
#' @param filenamePrefix 
#' @seealso \code{\link{approximateInverseWishart}}
#' @return Produces one Rdata file for each df scale factor
#'
generateApproximationData <- function(invWishartList, scaleMatrixList, dfScaleList,
                                      replicates=1000,
                                      outputDir=".", filenamePrefix="approxReplicates") {
  filenames=rep("", length(dfScaleList))
  idx = 1
  
  for(scale in dfScaleList) {
    
    # scale the degrees of freedom
    dfScaleInvWishartList = lapply(invWishartList, function(invWishart, scale) {
      return (new("inverseWishart", df=(invWishart@df * scale), 
                  precision=invWishart@precision))
    }, scale=scale)
    class(dfScaleInvWishartList) = "inverseWishart"
    
    # calculate the approximating inverse Wishart
    approxInvWishart = approximateInverseWishart(dfScaleInvWishartList, 
                                                 scaleMatrixList,
                                                 method="trace")
    
    # get replicates of the approximating inv Wishart
    approxReplicates = rInverseWishart(approxInvWishart, n=replicates)
    
    # save as an .Rdata file
    filename = paste(c(outputDir, "/", filenamePrefix, "_dfScale", scale, ".Rdata"), collapse="")
    save(approxReplicates, file=filename)
    
    filenames[idx] = filename
    idx = idx + 1
  }
  
  return(filenames)
}

#' getEnergyDistanceForFilelist
#'
#' Calculate the energy distance between samples contained in the specified
#' .Rdata files containing samples for the sum of inverse Wisharts and
#' samples for the approximating inverse Wishart.  
#'
#' @param sumFileList list of files containing the empirical samples
#' @param approxFileList list of files containing the samples for the 
#' approximating inverse Wishart
#' @return the list of energy distance values
#' 
#' @note This function assumes that the filename lists are in the same
#' order to pair the samples correctly
#' 
getEnergyDistanceForFilelist <- function(sumFileList, approxFileList) {
  if (length(sumFileList) != length(approxFileList)) {
    stop("the number of sum sample files must match the number of approx sample files")
  }
  
  # containing for energy distance values
  edistList = rep(0,length(sumFileList))
  
  for(i in 1:length(sumFileList)) {
    print(paste(c("Calculating energy distance for samples in files [", sumFileList[i], 
          "] and [", approxFileList[i], "]"), collapse=""))
    
    # load the sum sample from the Rdata file - sets a variable called sumReplicates
    load(sumFileList[i])
    # load the approx sample from the Rdata file - sets a variable called approxReplicates
    load(approxFileList[i])
    
    # calculate the energy distance

    edistList[i] = calculateEnergyDistance(sumReplicates, approxReplicates)
    
  }
  
  return(edistList)
  
}

#' plotDensityFromFile
#' 
#' Loads samples from files for the empirical sum of inverse Wishart matrices and a single 
#' inverse Wishart which approximates the distribution of the sum of inverse Wishart matrices.
#' Creates a comparison plot showing the univariate density of each cell.
#' 
#' @param sumFile file containing a sample from the distribution of a sum of inverse Wishart matrices
#' @param approxFile file containing a sample from a single inverse Wishart which approximates the distribution
#' of the sum of inverse Wishart matrices.
#' @param filename filename for the output pdf plot
#' @param ylim optional Y-axis plot limits 
#' @param xlim optional X-axis plot limits 
#' @param col optional list of colors for the empirical and approximate densities
#' @return plot of univariate densities as a pdf
#'
plotDensityFromFile <- function(sumFile, approxFile, filename,
                                ylim=c(0,1), xlim=c(0,1), col=c("blue", "black")) {
  print(paste(c("Plotting univariate densities for files [", sumFile, 
                "] and [", approxFile, "]"), collapse=""))
  
  # load the sum sample from the Rdata file - sets a variable called sumReplicates
  load(sumFile)
  # load the approx sample from the Rdata file - sets a variable called approxReplicates
  load(approxFile)
  # plot the densities
  pdf(filename)
  compare.plot(sumReplicates, approxReplicates, ylim=ylim, xlim=xlim, col=col)
  dev.off()
}

#' plotCovarianceFromFile
#' 
#' Loads samples from files for the empirical sum of inverse Wishart matrices and a single 
#' inverse Wishart which approximates the distribution of the sum of inverse Wishart matrices.
#' Creates a comparison plot showing the bivariate density of two cells.
#' 
#' @param sumFile file containing a sample from the distribution of a sum of inverse Wishart matrices
#' @param approxFile file containing a sample from a single inverse Wishart which approximates the distribution
#' of the sum of inverse Wishart matrices.
#' @param cell1 the row and column of the first cell in the bivariate density
#' @param cell2 the row and column of the second cell in the bivariate density
#' @param height height of the output pdf figure
#' @param width width of the output pdf figure
#' @param filename filename for the output pdf plot
#' @param lims plot limits 
#' @param col optional list of colors for the empirical and approximate densities
#' @return comparison plot of bivariate densities
#'
plotCovarianceFromFile <- function(sumFile, approxFile, cell1=c(1,1), cell2=c(1,2),
                                   height=4, width=4,
                                   filename, lims=c(-2,2,-2,2),
                                   col=c("blue", "black")) {
  print(paste(c("Plotting bivariate densities for files [", sumFile, 
                "] and [", approxFile, "]"), collapse=""))
  
  # load the sum sample from the Rdata file - sets a variable called sumReplicates
  load(sumFile)
  # load the approx sample from the Rdata file - sets a variable called approxReplicates
  load(approxFile)
  # plot the densities
  pdf(filename, height=height, width=width)
  compare.covarPlot(sumReplicates, approxReplicates, lims=lims, col=col, style="contour")
  
  dev.off()
}


#' runSimulationStudy
#' 
#' This function reproduces the simulation study results for the manuscript:\cr
#' Kreidler, S. M., Muller, K. E., & Glueck, D. H. An Approximation to 
#' the Distribution of the Sum of Inverse Wishart Matrices, In review.
#'
#' @param study.seed the random number seed (defaults to 1812)
#' @param study.data.dir the directory into which data files are written (defaults to
#' current working directory)
#' @param study.figures.dir the directory into which pdf figures are written (defaults
#' to the current working directory)
#' 
#' @note 
#' The energy distance calculations will take several minutes to a few hours 
#' to run depending on processor speed
#'
runSimulationStudy <- function(study.seed=1812, study.data.dir=".", study.figures.dir=".") {
  
  ##### Define Test Cases #####
  #
  # Test Case 1: sum of inverse ChiSquare variables (i.e. 1x1 inverse Wishart)
  #
  case1.invChiSqList = c(
    invert(new("wishart", df=6, covariance=matrix(c(1)))),
    invert(new("wishart", df=8, covariance=matrix(c(1)))),
    invert(new("wishart", df=9, covariance=matrix(c(1))))
  )
  
  # 
  # Test case 2: sum of 3x3 inverse Wishart matrices 
  #
  case2.invWishartList = c(
    invert(new("wishart", df=6, covariance=matrix(c(1,0.3,0.3,0.3,1,0.3,0.3,0.3,1), nrow=3))),
    invert(new("wishart", df=8, covariance=matrix(c(1,0.3,0.3,0.3,3,0.3,0.3,0.3,7), nrow=3))),
    invert(new("wishart", df=9, covariance=matrix(c(5,2,1,2,2,2,1,2,8), nrow=3)))
  )
  
  #
  # Test case 3: sum of quadratic forms sharing the same base covariance
  #
  case3.scaleMatrixList = list(
    matrix(c(1,0,0,0,0,0,0,0,1,0,0,0), nrow=2, byrow=TRUE),
    matrix(c(1,0,0,0,0,0,
             0,1,0,0,0,0,
             0,0,1,0,0,0), nrow=3, byrow=TRUE),
    cbind(matrix(rep(0,9), nrow=3, byrow=TRUE), diag(3))
  )
  sigma = matrix(c(1,0.3,0.4,0.3,1,0.3,0.4,0.3,1), nrow=3, byrow=3) 
  case3.invWishartList = c(
    invert(new("wishart", df=6, covariance=diag(3)[c(1,3),] %*% sigma %*% t(diag(3)[c(1,3),]))),
    invert(new("wishart", df=8, covariance=sigma)),
    invert(new("wishart", df=9, covariance=sigma))
  )
  
  ##### Generate samples for each test case and each scale factor #####
  set.seed(demo.seed)
  dfScaleList = 2^(0:5)
  replicates=1000
  # 
  # Generate empirical and approximate data for test case 1 (sum of chi-squared)
  #
  case1.approxFiles = generateApproximationData(case1.invChiSqList, NULL, dfScaleList=dfScaleList, 
                                                replicates=replicates, filenamePrefix="case1ApproxReplicates", 
                                                outputDir=demo.data.dir)
  case1.sumFiles = generateEmpiricalSumData(case1.invChiSqList, NULL, dfScaleList=dfScaleList, 
                                            replicates=replicates, filenamePrefix="case1SumReplicates", 
                                            outputDir=demo.data.dir)
  # 
  # Generate empirical and approximate data for test case 2 (sum of Wisharts)
  #
  case2.approxFiles = generateApproximationData(case2.invWishartList, NULL, dfScaleList=dfScaleList, 
                                                replicates=replicates, filenamePrefix="case2ApproxReplicates", 
                                                outputDir=demo.data.dir)
  case2.sumFiles = generateEmpiricalSumData(case2.invWishartList, NULL, dfScaleList=dfScaleList, 
                                            replicates=replicates, filenamePrefix="case2SumReplicates", 
                                            outputDir=demo.data.dir)
  # 
  # Generate empirical and approximate data for test case 3 (sum of quadratic forms of Wisharts)
  #
  case3.approxFiles = generateApproximationData(case3.invWishartList, case3.scaleMatrixList, 
                                                dfScaleList=dfScaleList, replicates=replicates, 
                                                filenamePrefix="case3ApproxReplicates", outputDir=demo.data.dir)
  case3.sumFiles = generateEmpiricalSumData(case3.invWishartList, case3.scaleMatrixList, 
                                            dfScaleList=dfScaleList, replicates=replicates,
                                            filenamePrefix="case3SumReplicates", outputDir=demo.data.dir)
  
  
  ##### Calculate the energy distance for each test case and scale factor #####
  #
  # Note - the energy distance is computationally expensive and may take several
  # minutes to run
  #
  case1.edist = getEnergyDistanceForFilelist(case1.sumFiles, case1.approxFiles)
  case2.edist = getEnergyDistanceForFilelist(case2.sumFiles, case2.approxFiles)
  case3.edist = getEnergyDistanceForFilelist(case3.sumFiles, case3.approxFiles)
  
  # Build a data frame with the energy distance results
  edist = data.frame(dfScale=dfScaleList,
                     chiSQEdist=case1.edist,
                     invWishartEdist=case2.edist,
                     singularInvWishartEdist=case3.edist)
  # write the energy distance table to disk
  write.csv(edist, paste(c(demo.data.dir, "energyDistance.csv"), collapse="/"), row.names=FALSE)
  
  # plot the energy distance by df scale factor
  edist = read.csv(paste(c(demo.data.dir, "energyDistance.csv"), collapse="/"))
  pdf(file=paste(c(demo.figures.dir, "energyDistanceByScaleFactor.pdf"), collapse="/"), family="Times")
  plot(edist$dfScale, edist$chiSQEdist, "l", cex.lab=1.5, ylim=c(0,10),
       xlab=expression(nu), ylab="Energy Distance", las=1)
  points(edist$dfScale, edist$chiSQEdist, pch=0)
  lines(edist$dfScale, edist$invWishartEdist, lty=2)
  points(edist$dfScale, edist$invWishartEdist, pch=1)
  lines(edist$dfScale, edist$singularInvWishartEdist, lty=3)
  points(edist$dfScale, edist$singularInvWishartEdist, pch=2)
  legend("topright", c("Inverse chi-squared sum", "Inverse Wishart sum", "Quadratic form sum"), 
         lty = c(1,2,3), cex=1.5, pch=c(0,1,2))
  dev.off()
  
  
  ##### Plot the univariate densities for the sum 
  ##### and the approximating inverse Wishart 
  #####
  idx=which(grepl("dfScale4.Rdata", case1.sumFiles))
  plotDensityFromFile(case1.sumFiles[idx], case1.approxFiles[idx], 
                      ylim=c(0,40), xlim=c(0,0.25), 
                      filename=paste(c(demo.figures.dir, "invChiSqDensity.pdf"), collapse="/"))
  plotDensityFromFile(case2.sumFiles[idx], case2.approxFiles[idx], 
                      ylim=c(0,80), xlim=c(-0.1,0.15), 
                      filename=paste(c(demo.figures.dir, "invWishartDensity.pdf"), collapse="/"))
  plotDensityFromFile(case3.sumFiles[idx], case3.approxFiles[idx], 
                      ylim=c(0,100), xlim=c(-0.1,0.15), 
                      filename=paste(c(demo.figures.dir, "quadraticFormDensity.pdf"), collapse="/"))
  
  ##### Plot selected covariances for the p > 1 cases #####
  idx=which(grepl("dfScale4.Rdata", case2.sumFiles))
  plotCovarianceFromFile(case2.sumFiles[idx], case2.approxFiles[idx], cell1=c(1,1),
                         cell2=c(1,2), height=4, width=6,
                         lims=c(-0.04,0,-0.025,0.01), 
                         filename=paste(c(demo.figures.dir, "invWishartCovar.pdf"), collapse="/"))
  plotCovarianceFromFile(case3.sumFiles[idx], case3.approxFiles[idx], cell1=c(1,1),
                         cell2=c(1,2), height=4, width=6,
                         lims=c(-0.03,0.03,-0.06,0), 
                         filename=paste(c(demo.figures.dir, "quadraticFormCovar.pdf"), collapse="/"))
}

