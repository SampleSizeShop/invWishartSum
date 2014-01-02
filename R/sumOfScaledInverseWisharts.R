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
# Equation to optimize for matching the expectation and log determinant
#
operand.logDeterminant = function(N, dim, weightedPrecisionSum, expectationSum) {
  return(
    log(det((N - dim - 1)*weightedPrecisionSum)) - dim*log(2) 
    - sum(sapply(1:dim, function(i) {return(digamma((N-dim-i)/2))}))
    - expectationSum
    
    )
}

# N = seq(10,10000,by=100)
# plot(N,sapply(N,function(obj) { return(operand.logDeterminant(obj,3,weightedPrecisionSum,expectationSum))}),
#      "l")

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
  # get the dimension of the inverse Wisharts
  dim = nrow(invWishartList[[1]]@precision)
  
  # weighted sum of precision matrices
  # sum of precision matrices weighted by the degrees of freedom
  weightedPrecisionSum = Reduce("+",lapply(invWishartList, function(obj, dim) {
    return((1/(obj@df - dim - 1))*obj@precision)
  }, dim=dim))
  
  # calculate the log determinant / digamm sum term
  expectationSum = sum(sapply(invWishartList, function(obj, dim) {
    log(det(obj@precision)) - dim*log(2) 
    - sum(sapply(1:dim, function(i) {return(digamma((obj@df-dim-i)/2))} ))    
  }, dim=dim))
    
  # start at the average of the df's for each wishart
  dfMean = mean(sapply(invWishartList, function(obj) { return(obj@df)}))
  # optimize the system of equations to obtain an expression for the
  # approximate degrees of freedom
  result = nlm(operand.logDeterminant, dfMean, dim=dim, 
               weightedPrecisionSum=weightedPrecisionSum, 
               expectationSum=expectationSum) 
  dfStar = result$estimate
  
  # use dfStar to calculate the precision matrix
  precisionStar = (dfStar - dim - 1) * precisionSum
  
  return(new("inverseWishart", df=dfStar, precision=precisionStar))
}

#
# Function to be optimized for the trace
#
operand.trace = function(N, dim, g1, g2, g3, g4) {
  return(
    ((2/(N - dim - 3))*g1
    + (4/((N - dim)*(N - dim - 3))*(g2+(N - dim - 1)*g3))
      - g4
    )^2)
}

#
# Calculates the variance of the trace of an inverse Wishart
#
traceVariance = function(invWishart) {
  # get the dimension of the inverse Wishart
  dim = nrow(invWishart@precision)
  # array of numbers 1 to dim
  indices = 1:dim
  
  # get the sum of the variances
  varianceSum = sum(sapply(indices, function(i, dim) {
    return (2 * invWishart@precision[i,i] / 
              ((invWishart@df - dim - 1)^2*(invWishart@df - dim - 3)))
  }, dim=dim))
  
  # generate a list of unique pairs of elements
  pairs = combn(indices, 2, simplify=FALSE)
  # get the sum of the covariances
  covarianceSum = 4 * sum(sapply(pairs, function(pair, dim) { 
    i = pair[1]
    j = pair[2]
    return(
      (invWishart@precision[i,i]*invWishart@precision[j,j] + 
         (invWishart@df - dim - 1)*invWishart@precision[i,j]^2) 
      / ((invWishart@df - dim)*(invWishart@df - dim - 1)^2*(invWishart@df - dim - 3))) 
  }, dim=dim))
  
  # return the variance of the trace
  return(varianceSum + 4 * covarianceSum)
  
}

#
# Convenience routine to calculate the trace
#
trace = function(m) {
  return (sum(diag(m)))
}

#
# Approximates the distribution of the sum of the inverse Wisharts
# by matching:
# - The expectation of the sum
# - The expectation of the trace of the sum
#
approximateInverseWishart.trace = function(invWishartList) {
  # get the dimension of the inverse Wisharts
  dim = nrow(invWishartList[[1]]@precision)
  
  # sum of precision matrices weighted by the degrees of freedom
  weightedPrecisionSum = Reduce("+",lapply(invWishartList, function(obj, dim) {
    return((1/(obj@df - dim - 1))*obj@precision)
  }, dim=dim))
  
  # sum of weighted averages squared of each diagonal element of the precision matrices
  weightedDiagElementSqSum = sum(sapply(1:dim, function(i, dim) {
    return(sum(sapply(invWishartList, function(obj, dim, i) {
      return((1/(obj@df - dim - 1))*obj@precision[i,i])
    }, dim=dim, i=i))^2)
  }, dim=dim))
  
  # generate a list of unique pairs of indices
  pairs = combn(1:dim, 2, simplify=FALSE)
  
  # sum of the weighted products of all pairs of diagonal elements of the
  # precision matrices
  weightedDiagElementProdSum = sum(sapply(pairs, function(pair, dim) {
    return(sum(sapply(invWishartList, function(obj, dim, i) {
      return((1/(obj@df - dim - 1))*obj@precision[i,i])
    }, dim=dim, i=pair[1]))
           *sum(sapply(invWishartList, function(obj, dim, i) {
             return((1/(obj@df - dim - 1))*obj@precision[i,i])
           }, dim=dim, i=pair[2])))
  }, dim=dim))
  
  # sum of weighted averages squared of each off-diagonal element
  weightedOffDiagElementSqSum = sum(sapply(pairs, function(pair, dim) {
    return(sum(sapply(invWishartList, function(obj, dim, i, j) {
      return((1/(obj@df - dim - 1))*obj@precision[i,j])
    }, dim=dim, i=pair[1], j=pair[2]))^2)
  }, dim=dim))
  
  # sum of the variances of the traces of the inverse Wisharts
  traceVarianceSum = sum(sapply(invWishartList, traceVariance)) 
    
  # start at the average of the df's for each wishart
  dfMean = mean(sapply(invWishartList, function(obj) { return(obj@df)}))
  # optimize the system of equations to obtain an expression for the
  # approximate degrees of freedom
  result = nlm(operand.trace, dfMean, dim=dim, 
               g1=weightedDiagElementSqSum,
               g2=weightedDiagElementProdSum,
               g3=weightedOffDiagElementSqSum,
               g4=traceVarianceSum)
            
  dfStar = result$estimate

  
  # use dfStar to calculate the precision matrix
  precisionStar = (dfStar - dim - 1) * weightedPrecisionSum
  
  return(new("inverseWishart", df=dfStar, precision=precisionStar))
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

