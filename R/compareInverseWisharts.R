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


#
# Functions to compare an approximate inverse Wishart to an empirical
# sum of inverse Wishart matrices
#
#
library(MASS)
#
# Simulate the sum of inverse Wishart matrices with zero-padding scale matrices
# 
simulate.sumInvWishartsScaled = function(invWishartList, scaleMatrixList, replicates=1000) {
  
  # get the Wisharts corresponding to the inverse Wisharts
  wishartList = lapply(invWishartList, invert)
  
  # simulate the Wishart replicates
  wishartReplicates = lapply(wishartList, function(obj, replicates) {
    replicateArray = rWishart(replicates, obj@df, obj@covariance)
    # convert the weird array from rWishart to a list of matrices
    wishartReplicates = list()
    for(i in 1:replicates) {
      wishartReplicates[[i]] = as.matrix(replicateArray[,,i])
    }
    return(wishartReplicates)
  }, replicates=replicates)
  
  # invert the Wisharts and scale them
  invWishartReplicates = lapply(1:length(invWishartList), function(i) {
    scaleMatrix = scaleMatrixList[[i]]
    replicateList = wishartReplicates[[i]]
    return (lapply(replicateList, function(replicate) { 
      t(scaleMatrix) %*% solve(replicate) %*% scaleMatrix
      }))
  })
  
  # sum them and return the result
  sumReplicates = lapply(1:replicates, function(i) {
    return (Reduce('+',lapply(invWishartReplicates, '[[', i)))
  })  
  
  return(sumReplicates)
}

#
# Simulate the sum of inverse Wishart matrices
# 
simulate.sumInvWisharts = function(invWishartList, replicates=1000) {
  
  # get the Wisharts corresponding to the inverse Wisharts
  wishartList = lapply(invWishartList, invert)
  
  # simulate the Wishart replicates
  wishartReplicates = lapply(wishartList, function(obj, replicates) {
    replicateArray = rWishart(replicates, obj@df, obj@covariance)
    # convert the weird array from rWishart to a list of matrices
    wishartReplicates = list()
    for(i in 1:replicates) {
      wishartReplicates[[i]] = as.matrix(replicateArray[,,i])
    }
    return(wishartReplicates)
  }, replicates)
  
  # invert the Wisharts
  invWishartReplicates = lapply(wishartReplicates, function(obj) {
    return (lapply(obj, solve))
  })
  
  # sum them and return the result
  sumReplicates = lapply(1:replicates, function(i) {
    return (Reduce('+',lapply(invWishartReplicates, '[[', i)))
  })  
  
  return(sumReplicates)
}

#
# Get a sample of inverse wisharts
#
rInverseWishart = function(invWishart, n=1000) {
  if (class(invWishart) != "inverseWishart") {
    stop("input is not an inverse Wishart object"); 
  }
  if (n <= 0) {
    stop("Invalid number of replicates")
  }
  
  # get the corresponding Wishart
  wishart = invert(invWishart)
  # generate replicates from said Wishart
  wishartArray = rWishart(n, wishart@df, wishart@covariance)
  # make a list rather than the array thing that comes back from the 
  # rWishart function
  wishartReplicates = lapply(1:n, function(i) { return(as.matrix(wishartArray[,,i])) })

  # invert those crazy wishart matrices and return
  return(lapply(wishartReplicates, solve))  
}

##### Functions for visually comparing densities #####

#
# Plot the density of a single cell of the empirical inverse Wishart
# sum and the approximating inverse Wishart density
#
plotWishartElement = function(empirical, approximate, row, column, 
                              col=c("black", "black"),
                              lty=c(1,1), ylim=c(0,1), xlim=c(0,1)) {
  empElt = sapply(1:length(empirical),function(i) {return (empirical[[i]][row,column])})  
  approxElt = sapply(1:length(approximate),function(i) {return (approximate[[i]][row,column])})
  
  if (sum(as.numeric(empElt!=0)) == 0) {
    # for singular Wisharts there may be all 0's in some cells
    plot(density(approxElt), main="", 
         ylim=ylim,
         xlim=xlim,
         xlab="",ylab="",xaxt='n',yaxt='n',
         col=col[2], lty=lty[2])
  } else {
    plot(density(empElt), main="", 
         ylim=ylim,
         xlim=xlim,
         xlab="",ylab="",xaxt='n',yaxt='n',
         col=col[1], lty=lty[1])
    lines(density(approxElt), col=col[2], lty=lty[2])
  }
  
}

#
# Simulate the specified number of replicates for the
# sum of inverse Wishart matrices and the approximating inverse Wishart.
# Then plot the result
#
compare.plot = function(sumReplicates, approxReplicates, ylim=c(0,1), xlim=c(0,1), 
                        col=c("red", "black")) {
  
  # plot the replicates
  dim = nrow(sumReplicates[[1]])
  
  #
  # Plot the density of each cell 
  #
  par(mfrow=c(dim,dim),mar=c(0,0,0,0), oma=c(1,1,1,1))
  for(r in 1:dim) {
    for(c in 1:dim) {
      plotWishartElement(
        empirical=sumReplicates, approximate=approxReplicates, 
        row=r, column=c, ylim=ylim, xlim=xlim, col=col)
    }
  }
  
}

#
# Produce a comparison plot of the covariance components
#
compare.covarPlot = function(empiricalReps, approxReps, cell1=c(1,2), cell2=c(1,3),
                             style="contour", lims=c(-2,2,-2,2), col=c("red", "black")) {
  empirical.kde <- kde2d(
    sapply(empiricalReps, function(x, cell) { 
      return(x[cell[1], cell[2]])
    }, cell=cell1),
    sapply(empiricalReps, function(x, cell) { 
      return(x[cell[1], cell[2]])
    }, cell=cell2), 
    lims=lims, n=50)
  
  approx.kde <- kde2d(
    sapply(approxReps, function(x, cell) { 
      return(x[cell[1], cell[2]])
    }, cell=cell1),
    sapply(approxReps, function(x, cell) { 
      return(x[cell[1], cell[2]])
    }, cell=cell2), 
    lims=lims, n=50)
  
  par(mfrow=c(1,2))
  if (style == "image") {
    image(empirical.kde, col = colorRampPalette(c("blue", "green", "yellow", "white"))(256))
    image(approx.kde, col = colorRampPalette(c("blue", "green", "yellow", "white"))(256))
  } else if (style == "persp") {
    persp(empirical.kde, phi = 45, theta = 30, border=col[1], zlim=c(0,5))
    persp(approx.kde, phi = 45, theta = 30, border=col[2], zlim=c(0,5))
  } else {
    contour(empirical.kde, col=col[1], asp=1)
    contour(approx.kde, col=col[2], asp=1)
  }
  
}

#
# Plot the difference of the covariances
#
compare.covarDiff = function(empiricalReps, approxReps, cell1=c(1,2), cell2=c(1,3),
                          style="contour", lims=c(-2,2,-2,2), col=c("red", "black")) {
  
  empirical.kde <- kde2d(
    sapply(empiricalReps, function(x, cell) { 
      return(x[cell[1], cell[2]])
    }, cell=cell1),
    sapply(empiricalReps, function(x, cell) { 
      return(x[cell[1], cell[2]])
    }, cell=cell2), 
    lims=lims, n=100)
  
  approx.kde <- kde2d(
    sapply(approxReps, function(x, cell) { 
      return(x[cell[1], cell[2]])
    }, cell=cell1),
    sapply(approxReps, function(x, cell) { 
      return(x[cell[1], cell[2]])
    }, cell=cell2), 
    lims=lims, n=100)
  
  diff.kde = empirical.kde
  diff.kde[[3]] = abs(empirical.kde[[3]] - approx.kde[[3]]) 
    
  if (style == "image") {
    image(diff.kde, col = colorRampPalette(c("blue", "green", "yellow", "white"))(256))
  } else if (style == "persp") {
    persp(empirical.kde, phi = 45, theta = 30, border=col[1], zlim=c(0,20))
    persp(approx.kde, phi = 45, theta = 30, border=col[2], zlim=c(0,20))
    persp(diff.kde, phi = 45, theta = 30, border="green", zlim=c(0,20))
  } else {
    contour(diff.kde, col=col[1])
  }
  
}

#
# Calculate the euclidean distance between two matrices
#
euclideanMatrixDistance <- function(m1, m2) {
  sqrt(sum(apply(m1 - m2, c(1,2), function(x) { return(x^2)})))   
}

#
# Estimate the energy distance between two inverse Wishart samples
#
# Assumes equal sample sizes
#
calculateEnergyDistance = function(sample1, sample2) {
  dim = nrow(sample1[[1]])
  n = length(sample1)
  
  # sum of distances between samples 1 and 2
  s1s1Dist = 0
  s2s2Dist = 0
  s1s2Dist = 0
  # calculate the sums
  for(i in 1:n) {
    for(j in 1:n) {
      s1s2Dist <- s1s2Dist + euclideanMatrixDistance(sample1[[i]], sample2[[j]])
      s1s1Dist <- s1s1Dist + euclideanMatrixDistance(sample1[[i]], sample1[[j]])
      s2s2Dist <- s2s2Dist + euclideanMatrixDistance(sample2[[i]], sample2[[j]])
    }
  }
  
  # calculate the energy distance
  eDist = (
    (n*n/(n+n)) *
      (2/(n*n) * (s1s2Dist) -
         1/n^2 * (s1s1Dist) -
         1/n^2 * (s2s2Dist)
        )
  )
  
  return(eDist)
  
}


