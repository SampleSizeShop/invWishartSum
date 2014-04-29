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
# Approximate the distribution of a sum of quadratic forms in
# inverse Wishart matrices
#

#
# Create a list of scale matrices for each Wishart
#
scaleMatrixList = list(
  matrix(c(1,0,0,0,0,0,0,0,1,0,0,0), nrow=2, byrow=TRUE),
  matrix(c(1,0,0,0,0,0,
           0,1,0,0,0,0,
           0,0,1,0,0,0), nrow=3, byrow=TRUE),
  cbind(matrix(rep(0,9), nrow=3, byrow=TRUE), diag(3))
)

#
# Create a list of inverse Wishart matrices with dimension
# either 3x3 or 2x2
#
sigma = matrix(c(1,0.3,0.4,0.3,1,0.3,0.4,0.3,1), nrow=3, byrow=3) 
quadFormList = c(
  invert(new("wishart", df=26, covariance=diag(3)[c(1,3),] %*% sigma %*% t(diag(3)[c(1,3),]))),
  invert(new("wishart", df=38, covariance=sigma)),
  invert(new("wishart", df=29, covariance=sigma))
)

# calculate the approximating inverse Wishart
approxInvWishart = approximateInverseWishart(quadFormList, scaleMatrixList=scaleMatrixList)

# output the approximating Wishart
approxInvWishart

