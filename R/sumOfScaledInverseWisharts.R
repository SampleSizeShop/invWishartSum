#
# Functions to calculate an approximation to the distribution of
# the sum of inverse Wishart matrices.
#
# Also handles pseudo-inverse Wishart matrices
#
# Provides three methods for approximating the distribution of the sum
#
# trace - matches the expectation of the sum and the expectation of the trace
# (from Kreidler et al. 2014)
#
# cellVariance - matches the expectation of the sum and the variance of a single
#  cell. (from Kreidler et al. 2014)
#
# logDeterminant - matches the expectation of the sum and the expectation of the
#  log determinant of the sum (from Granstrom et al. 2012)
#
#


#
# Calculate the degrees of freedom and scale matrix of the
# inverse Wishart which best approximates the distribution
# of sum of the specified list of Wishart matrices
#
# This function implements the method described by Granstrom and Orguner (2012)
# which matches the expectation and the expectation of the log determinant.
# This approach requires numerical optimization
#
#
approximateInverseWishart.logDeterminant = function(invWishartList) {
  
}

#
#
#
#
approximateInverseWishart.trace = function(invWishartList) {
  
}

#
# Approximates the distribution of the sum of the inverse Wisharts
# by matching:
# - The expectation of the sum
# - The variance of the specified cell of the sum
#
approximateInverseWishart.cell = function(invWishartList, row=1, column=1) {

  # sum of cell / (df - dim - 1)
  sum1 = sum(sapply(invWishartList, function(obj, row, column) {
    dim = nrow(obj@precision)
    return (obj@precision[row, column]/(obj@df - dim - 1))
  }, row=row, column=column))
  
  # sum of cell^2 / ((df - dim - 1)^2*(df - dim - 3))
  sum2 = sum(sapply(invWishartList, function(obj, row, column) {
    dim = nrow(obj@precision)
    return ((obj@precision[row, column])^2/((obj@df - dim - 1)^2*(obj@df - dim - 3)))
  }, row=row, column=column))
  
  # get the dimension of the inverse Wisharts
  dim = nrow(invWishartList[[1]]@precision)
  
  # calculate the degrees of freedom for the approximate Wishart
  dfStar = (dim + 3 + (sum1^2 * (1/sum2))) 

  # calculate the precision matrix for the approximate Wishart
  precisionStar = ((2 + (sum1^2*(1/sum2)))
                   *Reduce("+",lapply(invWishartList, function(obj) {
                     dim = nrow(obj@precision)
                     return((1/(obj@df - dim - 1))*obj@precision)
                   })))
    
  return(new("inverseWishart", df=dfStar, precision=precisionStar))  
}

#
# Calculate an inverse Wishart or pseudo inverse Wishart which
# best approximates the distribution of the sum of the specified
# inverse Wishart matrices.
#
# Optionally, scale matrices can be specified which will be pre and post
# multiplied onto each Wishart.  For example, for the inputs:
# wishartList = (X1, X2)
# scaleMatrixList = (A,B)
#
# this function will calculate the approximate distribution of
#  A'X1A + B'X2B
#
# When the method="cellVariance" is used, the user must specify the cell
# which will be matched
#
approximateInverseWishart = function(invWishartList, scaleMatrixList=NULL, method="trace", cell=NULL) {
  # TODO: check class of invWishartList
  
  # scale the inverse Wisharts if needed
  scaledInvWishartList = invWishartList
  if (!is.null(scaleMatrixList)) {
    if (length(scaleMatrixList) != length(invWishartList)) {
      stop("The scale matrix list and wishart list do not have the same length")
    }  
    scaledInvWishartList = lapply(i:length(invWishartList), function(i) {
      return (t(scaleMatrixList[i]) %*% invWishartList[i]@covariance %*% scaleMatrixList[i])
    })
  }
  
  # make sure all of the scaled Wisharts are the same size
  dimension = nrow(scaledInvWishartList[[1]]@precision)
  if (sum(sapply(scaledInvWishartList, function(obj, dimension) { 
    return(ifelse(nrow(obj@precision) != dimension, 1, 0))}, dimension)) > 0) {
    stop("All inverse wishart matrices must have equal dimension")
  }

  #
  # Call appropriate approximation function
  #
  switch(method,
         trace = {
           return(approximateInverseWishart.trace(scaledInvWishartList))
         },
         logDeterminant = {
           return(approximateInverseWishart.logDeterminant(scaledInvWishartList))
         },
         cellVariance = {
           if (is.null(cell) || length(cell) != 2) {
             stop("Invalid cell specified")
           }
           row = cell[1]
           column = cell[2]
           dim = nrow(scaledInvWishartList[[1]]@precision)
           if (row <= 0 || row > dim || column <= 0 || column > dim) {
             stop("Specified cell for variance match is out of bounds")
           }
           return(approximateInverseWishart.cell(scaledInvWishartList, row, column))
         },
         stop("Unknown approximation method (valid values are 'trace', 'logDeterminant', or 'cellVariance'")
  )       
}

