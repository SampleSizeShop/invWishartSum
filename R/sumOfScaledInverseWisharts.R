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


#' operand.logDeterminant
#' 
#' Equation used for optimization when matching the expectation and log determinant to 
#' approximate the distribution of a sum of inverse Wishart matrices.  Optimization
#' occurs across the degrees of freedom, N.
#' 
#' @param N the current degrees of freedom
#' @param dim the dimension of the inverse Wishart matrices in the sum
#' @param weightedPrecisionSum weighted sum of precision matrices
#' @param expectationSum sum of the expected values of each inverse Wishart
#' @return will return 0 or near 0 at optimal value of N
#' @keywords internal
#' @note
#' Based on the approach of:
#' Granstrom, K., & Orguner, U. (2012). On the reduction of Gaussian inverse Wishart mixtures. 
#' In 2012 15th International Conference on Information Fusion (FUSION) (pp. 2162-2169).
#'
operand.logDeterminant = function(N, dim, weightedPrecisionSum, expectationSum) {
  return(
    (log(det((N - dim - 1)*weightedPrecisionSum)) - dim*log(2) 
    - sum(sapply(1:dim, function(i) {return(digamma((N-dim-i)/2))}))
    - expectationSum
    
    )^2)
}

# N = seq(10,10000,by=100)
# plot(N,sapply(N,function(obj) { return(operand.logDeterminant(obj,3,weightedPrecisionSum,expectationSum))}),
#      "l")

#' approximateInverseWishart.logDeterminant
#' 
#' Approximate the distribution of a sum of inverse Wishart matrices with
#' a single inverse Wishart.The two-moment approximation matches the expectation
#' of the sum and the expectation of the log determinant of the sum.
#' The solution to the resulting equations is obtained by numerical optimization.
#' 
#' @param invWishartList the list of inverse Wishart matrices in the sum
#' @return the approximating inverse Wishart object
#' @seealso \code{\link{inverseWishart}}
#' @references
#' Based on the approach of:
#' Granstrom, K., & Orguner, U. (2012). On the reduction of Gaussian inverse Wishart mixtures. 
#' In 2012 15th International Conference on Information Fusion (FUSION) (pp. 2162-2169).
#'
approximateInverseWishart.logDeterminant = function(invWishartList) {
  # get the dimension of the inverse Wisharts
  dim = nrow(invWishartList[[1]]@precision)
  
  # weighted sum of precision matrices
  # sum of precision matrices weighted by the degrees of freedom
  weightedPrecisionSum = solve(Reduce("+",lapply(invWishartList, function(obj, dim) {
    return(((obj@df - dim - 1))*solve(obj@precision))
  }, dim=dim)))
  
  
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
  precisionStar = (dfStar - dim - 1) * weightedPrecisionSum
  
  return(new("inverseWishart", df=dfStar, precision=precisionStar))
}

#' operand.trace
#' 
#' Equation used for optimization when matching the expectation and variance of the trace to 
#' approximate the distribution of a sum of inverse Wishart matrices.  Optimization
#' occurs across the degrees of freedom, N. 
#' 
#' @param N the current degrees of freedom
#' @param dim the dimension of the inverse Wishart matrices in the sum
#' @param weightedPrecisionSum weighted sum of precision matrices
#' @param expectationSum sum of the expected values of each inverse Wishart
#' @return will return 0 or near 0 at optimal value of N
#' @keywords internal
#' @note
#' THIS FUNCTION IS OBSOLETE SINCE THIS APPROACH HAS A CLOSED-FORM SOLUTION
#'
operand.trace = function(N, dim, g1, g2, g3, g4) {
  return(
    ((2/(N - dim - 3))*g1
    + (4/((N - dim)*(N - dim - 3))*(g2+(N - dim - 1)*g3))
      - g4
    )^2)
}


#' traceVariance
#' 
#' Calculate the variance of the trace of an inverse Wishart matrix
#' 
#' @param invWishart the inverse Wishart object
#' @seealso \code{\link{inverseWishart}}
#' @return the variance of the trace of the specified inverse Wishart
#'
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
  pairs = ifelse(dim > 1, combn(1:dim, 2, simplify=FALSE), list(c(1,1)))
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

#' approximateInverseWishartScaled.trace
#' 
#' Approximate the distribution of a sum of quadratic forms in inverse Wishart 
#' matrices with a single inverse Wishart. The two-moment approximation matches 
#' the expectation of the sum and the variance of the trace of the sum.
#' 
#' @param invWishartList the list of inverse Wishart matrices in the sum
#' @param scaleMatrixList the list of matrices which will be pre- and post-multiplied
#' onto each corresponding inverse Wishart to build the quadratic forms
#' @return the approximating inverse Wishart object
#' @seealso \code{\link{inverseWishart}}
#' 
#' @references
#' Kreidler, S. M., Muller, K. E., & Glueck, D. H. An Approximation to 
#' the Distribution of the Sum of Inverse Wishart Matrices, In review.
#'
approximateInverseWishartScaled.trace = function(invWishartList, scaleMatrixList) {

  # scaled precision matrices
  scaledPrecisionMatrices = lapply(1:length(scaleMatrixList), function(i) {
    scaleMatrix = scaleMatrixList[[i]]
    invWishart = invWishartList[[i]]
    return (t(scaleMatrix) %*% invWishart@precision %*% scaleMatrix)
  })
  # get the dimension of the approximating inverse Wishart
  dim = nrow(scaledPrecisionMatrices[[1]])
  
  # sum of precision matrices weighted by the degrees of freedom
  weightedPrecisionSum = Reduce("+",lapply(1:length(scaleMatrixList), function(i) {
    precision = scaledPrecisionMatrices[[i]]
    invWishart = invWishartList[[i]]
    return((1/(invWishart@df - nrow(invWishart@precision) - 1)) * precision)
  }))    
  
  # sum of weighted averages squared of each diagonal element of the precision matrices
  weightedDiagElementSqSum = sum(sapply(1:dim, function(cellIdx) {
    return(sum(sapply(1:length(scaleMatrixList), function(i, cellIdx) {
      precision = scaledPrecisionMatrices[[i]]
      invWishart = invWishartList[[i]]
      return((1/(invWishart@df - nrow(invWishart@precision) - 1))*precision[cellIdx,cellIdx])
    }, cellIdx=cellIdx))^2)
  }))
  
  # generate a list of unique pairs of indices
  pairs = ifelse(dim > 1, combn(1:dim, 2, simplify=FALSE), list(c(1,1)))
  
  # sum of the weighted products of all pairs of diagonal elements of the
  # precision matrices
  weightedDiagElementProdSum = sum(sapply(pairs, function(pair) {
    return(sum(sapply(1:length(scaleMatrixList), function(i, cellIdx) {
      precision = scaledPrecisionMatrices[[i]]
      invWishart = invWishartList[[i]]
      return((1/(invWishart@df - nrow(invWishart@precision) - 1))*precision[cellIdx,cellIdx])
    }, cellIdx=pair[1]))
           * sum(sapply(1:length(scaleMatrixList), function(i, cellIdx) {
             precision = scaledPrecisionMatrices[[i]]
             invWishart = invWishartList[[i]]
             return((1/(invWishart@df - nrow(invWishart@precision) - 1))*precision[cellIdx,cellIdx])
           }, cellIdx=pair[2])))
  }))
  
  # sum of weighted averages squared of each off-diagonal element
  weightedOffDiagElementSqSum = sum(sapply(pairs, function(pair) {
    return(sum(sapply(1:length(scaleMatrixList), function(i, cellIdx1, cellIdx2) {
      precision = scaledPrecisionMatrices[[i]]
      invWishart = invWishartList[[i]]
      return((1/(invWishart@df - nrow(invWishart@precision) - 1))*precision[cellIdx1,cellIdx2])
    }, cellIdx1=pair[1], cellIdx2=pair[2]))^2)
  }))
  
  # sum of the variances of the traces of the inverse Wisharts
  traceVarianceSum = sum(sapply(invWishartList, traceVariance)) 
  
  # start at the average of the df's for each wishart
  dfMean = mean(sapply(invWishartList, function(obj) { return(obj@df)}))
  
  #
  # calculate intermediate values
  #
  b = -(2 * weightedDiagElementSqSum + 
         4 * weightedOffDiagElementSqSum +
         2 * traceVarianceSum * dim + 
         3 * traceVarianceSum)
  c = (2 * weightedDiagElementSqSum * dim -
         4 * weightedDiagElementProdSum + 
         4 * weightedOffDiagElementSqSum * dim + 
         4 * weightedOffDiagElementSqSum + traceVarianceSum * dim^2 +
         3 * traceVarianceSum * dim)
  
  tmp = (-1*b + c(-1,1)*sqrt(b^2 - 4 * traceVarianceSum * c)) /
    (2 * traceVarianceSum)

  # calculate approximate N
  dfStar = max(tmp)
  
  # use dfStar to calculate the precision matrix
  precisionStar = (dfStar - dim - 1) * weightedPrecisionSum

  return(new("inverseWishart", df=dfStar, precision=precisionStar))
}


#' approximateInverseWishart.trace
#' 
#' Approximate the distribution of a sum of inverse Wishart 
#' matrices with a single inverse Wishart.
#' 
#' The two-moment approximation matches the expectation
#' of the sum and the variance of the trace of the sum.
#' 
#' @param invWishartList the list of inverse Wishart matrices in the sum
#' @return the approximating inverse Wishart object
#' @seealso \code{\link{inverseWishart}}
#' 
#' @references
#' Kreidler, S. M., Muller, K. E., & Glueck, D. H. An Approximation to 
#' the Distribution of the Sum of Inverse Wishart Matrices, In review.
#'
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
  pairs = ifelse(dim > 1, combn(1:dim, 2, simplify=FALSE), list(c(1,1)))
  
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
  
  
  #
  # calculate intermediate values
  #
  b = (2 * weightedDiagElementSqSum + 
         4 * weightedOffDiagElementSqSum +
         2 * traceVarianceSum * dim + 
         3 * traceVarianceSum)
  c = (2 * weightedDiagElementSqSum * dim -
         4 * weightedDiagElementProdSum + 
         4 * weightedOffDiagElementSqSum * dim + 
         4 * weightedOffDiagElementSqSum + traceVarianceSum * dim^2 +
         3 * traceVarianceSum * dim)
  
  # calculate approximate N
  dfStar = max(
    (-1*b + c(-1,1)*sqrt(b^2 - 4 * traceVarianceSum * c)) /
      (-2 * traceVarianceSum)
    )
  
  # use dfStar to calculate the precision matrix
  precisionStar = (dfStar - dim - 1) * weightedPrecisionSum
  
  return(new("inverseWishart", df=dfStar, precision=precisionStar))

}

#' approximateInverseWishart.cell
#' 
#' Approximate the distribution of a sum of inverse Wishart 
#' matrices with a single inverse Wishart.
#' 
#' The two-moment approximation matches the expectation
#' of the sum and the variance of the specified cell of the sum.
#' 
#' @param invWishartList the list of inverse Wishart matrices in the sum
#' @param row the row of the cell used in the variance match
#' @param column the column of the cell used in the variance match
#' @return the approximating inverse Wishart object
#' @seealso \code{\link{inverseWishart}}
#' 
#' @note
#' This approach has lower accuracy that the \code{\link{approximateInverseWishart.trace}}
#' function and is not recommended.
#' 
#'
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

#' approximateInverseWishart
#' 
#' Approximate the distribution of a sum of inverse Wishart 
#' matrices or a sum of quadratic forms in inverse Wishart matrices
#' with a single inverse Wishart.
#'
#' Three approximation methods are currently supported
#' \enumerate{
#' \item{\code{trace}: }{The default method which matches the expectation of the sum and the 
#' variance of the trace of the sum}
#' \item{\code{logDeterminant: }}{A method which matches the expectation of the sum and the 
#' expectation of the log determinant of the sum}
#' \item{\code{cell}: 
#' }{A method which matches the expectation of the sum and the 
#' variance of a specified cell of the sum}
#' }
#' 
#' For quadratic forms, a list of scale matrices should be specified which will be pre- and post-
#' multiplied onto each inverse Wishart.  For example, for the inputs 
#' \code{invWishartList = c(X1, X2)} and \code{scaleMatrixList = c(A,B)}, 
#' the function will calculate the approximate distribution of
#'  \code{A'X1A + B'X2B}
#'
#' @param invWishartList the list of inverse Wishart matrices in the sum
#' @param method the moment-matching approach to use.  Valid values are \code{trace},
#' \code{log-determinant}, and \code{cell} 
#' @param scaleMatrixList optional, list of scale matrices to form a sum of quadratic forms.  
#' There must be one scale matrix for each inverse Wishart in the \code{invWishartList}.
#' @param cell optional, a vector containing the row and column of the cell used for variance
#' matching when using \code{method='cell'} 
#' @return the approximating inverse Wishart object
#' @seealso \code{\link{inverseWishart}}
#' 
#' @references
#' The \code{trace} method implements the approach of:\cr
#' Kreidler, S. M., Muller, K. E., & Glueck, D. H. An Approximation to 
#' the Distribution of the Sum of Inverse Wishart Matrices, In review.
#' 
#' The \code{logDeterminant} method is based on the approach of:\cr
#' Granstrom, K., & Orguner, U. (2012). On the reduction of Gaussian inverse Wishart mixtures. 
#' In 2012 15th International Conference on Information Fusion (FUSION) (pp. 2162-2169).
#'
#'
approximateInverseWishart = function(invWishartList, method="trace", 
                                     scaleMatrixList=NULL, cell=NULL) {
  # check class of invWishartList
  if (class(invWishartList) != "inverseWishart") {
    stop("input list does not contain inverse Wisharts")
  }
  
  # make sure the wisharts and scale matrices are valid
  if (!is.null(scaleMatrixList)) {
    if (length(scaleMatrixList) != length(invWishartList)) {
      stop("The number of inverse Wisharts must match the number of scale matrices")
    }
    # make sure all of the scale matrices conform with their inverse Wishart,
    # as well as each otherare the same size
    dimension = ncol(scaleMatrixList[[1]])
    if (sum(sapply(1:length(scaleMatrixList), function(i, dimension) { 
      invWishart = invWishartList[[i]]
      scaleMatrix = scaleMatrixList[[i]]
      return(ifelse(ncol(scaleMatrix) != dimension
                    || nrow(scaleMatrix) != nrow(invWishart@precision), 
                    1, 0))}, dimension)) > 0) {
      stop("Scale matrices do not conform")
    }
  } else {
    # we aren't scaling, so make sure the inverse Wisharts have the same dimension
    dimension = nrow(invWishartList[[1]]@precision)
    if (sum(sapply(invWishartList, function(obj, dimension) { 
      return(ifelse(nrow(obj@precision) != dimension, 1, 0))}, dimension)) > 0) {
      stop("All inverse Wisharts must have the same dimension")
    }
  }

  #
  # Call appropriate approximation function
  #
  if (method == "trace") {
    if (!is.null(scaleMatrixList)) {
      return(approximateInverseWishartScaled.trace(invWishartList, scaleMatrixList))
    } else {
      return(approximateInverseWishart.trace(invWishartList))
    }    
    
  } else if (method == "logDeterminant") {
    if (!is.null(scaleMatrixList)) {
      stop("scale matrices are only supported for the 'trace' method")
    }
    return(approximateInverseWishart.logDeterminant(invWishartList))
    
  } else if (method == "cellVariance") {
    if (!is.null(scaleMatrixList)) {
      stop("scale matrices are only supported for the 'trace' method")
    }
    if (is.null(cell) || length(cell) != 2) {
      stop("Invalid cell specified")
    }
    row = cell[1]
    column = cell[2]
    dim = nrow(invWishartList[[1]]@precision)
    if (row <= 0 || row > dim || column <= 0 || column > dim) {
      stop("Specified cell for variance match is out of bounds")
    }
    return(approximateInverseWishart.cell(invWishartList, row, column))
  } else {
    stop("Unknown approximation method (valid values are 'trace', 'logDeterminant', or 'cellVariance'")
  }
      
}



