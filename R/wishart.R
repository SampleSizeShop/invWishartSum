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

#
# check if a matrix is square
#
isSquare = function(x) {
  if (!is.matrix(x)) {
    stop("argument must be a matrix")
  }
  return(nrow(x) == ncol(x))
}

#
# check if a matrix is positive definite
#
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

#
# Class representing the Wishart distribution
#
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

#
# Class representing the inverse Wishart distribution
#
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

# generic for inverting a Wishart or inverse Wishart distribution
setGeneric(
  "invert",
  function(object) {
    standardGeneric("invert")
  }
)
# 
# For a given Wishart, create the corresponding inverse Wishart
# Based on Theorem 3.4.1, p 111 from Gupta & Nagar (2000)
#
setMethod("invert", 
          signature("wishart"),
          function(object) {
            dimension = nrow(object@covariance)
            precision = solve(object@covariance)
            df = object@df + dimension + 1
            return(new("inverseWishart", df=df, precision=precision))
          })
# 
# For a given Wishart, create the corresponding inverse Wishart
# Based on Theorem 3.4.1, p 111 from Gupta & Nagar (2000)
#
setMethod("invert", 
          signature("inverseWishart"),
          function(object) {
            dimension = nrow(object@precision)
            covariance = solve(object@precision)
            df = object@df - dimension - 1
            return(new("wishart", df=df, covariance=covariance))
          })

#
# Concatenation function for Wisharts
#
setMethod("c", signature(x = "wishart"), function(x, ...){
  elements = list(x, ...)
  
  wishartList = list()
  for (i in 1:length(elements)){
    wishartList[i] = new("wishart", df = elements[[i]]@df, covariance = elements[[i]]@covariance)  
  }
  
  class(wishartList) = "wishart"
  
  wishartList 
})
#
# Concatenation function for inverse Wisharts
#
setMethod("c", signature(x = "inverseWishart"), function(x, ...){
  elements = list(x, ...)
  
  invWishartList = list()
  for (i in 1:length(elements)){
    invWishartList[i] = new("inverseWishart", 
                            df = elements[[i]]@df, precision = elements[[i]]@precision)  
  }
  
  class(invWishartList) = "inverseWishart"
  
  invWishartList 
})





