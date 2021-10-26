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
# Classes describing a Wishart distribution, an inverse Wishart
# and methods for moving between the two.
#
#
#

#' isSquare 
#' 
#' check if a matrix is square
#' 
#' @param x a \code{matrix} object
#' @return TRUE if the matrix is square, and FALSE otherwise
#' @keywords internal
isSquare = function(x) {
  if (!is.matrix(x)) {
    stop("argument must be a matrix")
  }
  return(nrow(x) == ncol(x))
}

#' isPositiveDefinite
#' 
#' check if a matrix is positive definite
#'
#' @param x a \code{matrix} object
#' @param tol a value indicating the smallest integer considered positive
#' @return TRUE if the matrix is positive definite, and FALSE otherwise
#' @keywords internal
isPositiveDefinite = function(x, tol=1e-18) {
  if (!is.matrix(x) || !isSquare(x)) {
    stop("argument must be a square matrix")
  }
  eigenValue = sapply(eigen(x,only.values=TRUE)$values, 
                      function(obj) {return(ifelse(abs(obj) < tol, 0, obj))})
  if (any(eigenValue <= 0)) {
    return(FALSE)
  }
  return(TRUE)
}

#'
#' wishart
#'
#' Class representing the Wishart distribution 
#'
#' @slot df The degrees of freedom for the Wishart distribution.  Must be a positive integer.
#' @slot covariance The covariance of the Wishart distribution.  Must be a square, positive
#' definite matrix.    
#'
#' @name wishart 
#' @rdname wishart
#' @note For theoretical details, please see
#' 
#'Gupta, A. K., & Nagar, D. K. (2000). Matrix variate distributions. 
#'Boca Raton, FL: Chapman & Hall.
#' 
setClass("wishart",
         representation ( df = "numeric",
                          covariance = "matrix"
         ),
         prototype(df=10,
                   covariance=matrix(c(1,0.2,0.2,1), nrow=2)
                   ),
         validity = function(object) {
           if (object@df < 0) {
             stop("The df must be positive")
           }
           if (!isPositiveDefinite(object@covariance)) {
             stop("The covariance matrix must be positive definite.")
           }
           return(TRUE);
         }
         )

#'
#' wishart
#'
#' Class representing the inverse Wishart distribution 
#'
#' @slot df The degrees of freedom for the Wishart distribution.  Must be a positive integer.
#' @slot precision The precision matrix of the inverse Wishart distribution.  Must be a square, positive
#' definite matrix.    
#'
#' @name inverseWishart 
#' @rdname inverseWishart
#' @note For theoretical details, please see
#' 
#'Gupta, A. K., & Nagar, D. K. (2000). Matrix variate distributions. 
#'Boca Raton, FL: Chapman & Hall.
#' 
setClass("inverseWishart",
         representation ( df = "numeric",
                          precision = "matrix"
         ),
         prototype(df=10,
                   precision=matrix(c(1,0.2,0.2,1), nrow=2)
         ),
         validity = function(object) {
           if (object@df < 0) {
             stop("The df must be positive")
           }
           if (!isPositiveDefinite(object@precision)) {
             stop("The covariance matrix must be positive definite.")
           }
           return(TRUE);
         }
)

#' invert
#' 
#' Form the distribution of the inverse of the specified Wishart or inverse Wishart.
#'
#' @param a \code{\link{wishart}} or \code{\link{inverseWishart}} distribution object
#'
#' @return the distribution object for the inverse of the specified Wishart or inverse Wishart
#' 
#' @seealso \code{\link{wishart}} and \code{\link{inverseWishart}}
#' 
#' @export
#' @docType methods
#' @rdname invert-methods
#' @name invert
#' @note
#' Based on Theorem 3.4.1, p 111 from 
#' 
#' Gupta, A. K., & Nagar, D. K. (2000). Matrix variate distributions. 
#' Boca Raton, FL: Chapman & Hall.
#' 
setGeneric(
  "invert",
  function(object) {
    standardGeneric("invert")
  }
)

#' @rdname invert-methods
#' @aliases invert,wishart,wishart-method
setMethod("invert", 
          signature("wishart"),
          function(object) {
            dimension = nrow(object@covariance)
            precision = solve(object@covariance)
            df = object@df + dimension + 1
            return(new("inverseWishart", df=df, precision=precision))
          })

#' @rdname invert-methods
#' @aliases invert,inverseWishart,inverseWishart-method
setMethod("invert", 
          signature("inverseWishart"),
          function(object) {
            dimension = nrow(object@precision)
            covariance = solve(object@precision)
            df = object@df - dimension - 1
            return(new("wishart", df=df, covariance=covariance))
          })

#' Combine Values into a vector or list
#' 
#' Concatenate Wishart objects
#'
#' @rdname c-methods
#' @aliases c,wishart,wishart-method
setMethod("c", signature(x = "wishart"), function(x, ...){
  elements = list(x, ...)
  
  wishartList = list()
  for (i in 1:length(elements)){
    wishartList[[i]] = new("wishart", df = elements[[i]]@df, covariance = elements[[i]]@covariance)  
  }
  
  class(wishartList) = "wishart"
  
  wishartList 
})

#' Concatenate inverse Wishart objects
#'
#' @rdname c-methods
#' @aliases c,inverseWishart,inverseWishart-method
setMethod("c", signature(x = "inverseWishart"), function(x, ...){
  elements = list(x, ...)
  
  invWishartList = list()
  for (i in 1:length(elements)){
    invWishartList[[i]] = new("inverseWishart", 
                            df = elements[[i]]@df, precision = elements[[i]]@precision)  
  }
  
  class(invWishartList) = "inverseWishart"
  
  invWishartList 
})





