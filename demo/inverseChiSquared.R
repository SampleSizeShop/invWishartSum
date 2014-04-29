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
# Define a list of inverse Wishart distributions.  In this case,
# we define inverse Chi-square variables
#
invChiSqList = c(
  invert(new("wishart", df=23, covariance=matrix(c(1)))),
  invert(new("wishart", df=14, covariance=matrix(c(1)))),
  invert(new("wishart", df=30, covariance=matrix(c(1))))
)

# calculate the approximating inverse Wishart
approxInvWishart = approximateInverseWishart(invChiSqList)

# output the approximating Wishart
approxInvWishart

