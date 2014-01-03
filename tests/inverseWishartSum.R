#
# Test cases of a sums of inverse Wishart matrices
#
#
set.seed(2014)

# Wishart Test Case 1
# Specify our list of inverse Wisharts to sum
#
case1.invWishartList = c(
  invert(new("wishart", df=10, covariance=matrix(c(1,0.3,0.3,0.3,1,0.3,0.3,0.3,1), nrow=3))),
  invert(new("wishart", df=8, covariance=matrix(c(1,0.3,0.3,0.4,3,0.3,0.3,0.3,7), nrow=3))),
  invert(new("wishart", df=6, covariance=matrix(c(5,2,1,2,2,2,1,2,8), nrow=3)))
  )

checkApproximation(case1.invWishartList, ylim=c(0,4), xlim=c(-2,4), method="cellVariance")
checkApproximation(case1.invWishartList, ylim=c(0,4), xlim=c(-2,4), method="trace")
checkApproximation(case1.invWishartList, ylim=c(0,4), xlim=c(-2,4), method="logDeterminant")

#
# Wishart Test Case 2
#
case2.invWishartList = c(
  invert(new("wishart", df=10, 
             covariance=matrix(c(1, 0.8, 0.5, 0.8, 1, 0.8, 0.5, 0.8, 1), 
                               nrow=3, byrow=TRUE))),
  invert(new("wishart", df=7, 
             covariance=matrix(c(1, 0.2, 0.2, 0.2, 1, 0.2, 0.2, 0.2, 1), 
                               nrow=3, byrow=TRUE)))
  )
case2.xlim=c(-4,4)
case2.ylim=c(0,2)
checkApproximation(case2.invWishartList, xlim=case2.xlim, ylim=case2.ylim, method="cellVariance")
checkApproximation(case2.invWishartList, xlim=case2.xlim, ylim=case2.ylim, method="trace")
checkApproximation(case2.invWishartList, xlim=case2.xlim, ylim=case2.ylim, method="logDeterminant")

#
# Scaled, zero-padded Wisharts
#
case3.scaleMatrixList = list(
  matrix(c(1,0,0,0,0,0,0,0,1,0,0,0), nrow=2, byrow=TRUE),
  matrix(c(1,0,0,0,0,0,
           0,1,0,0,0,0,
           0,0,1,0,0,0), nrow=3, byrow=TRUE),
  cbind(matrix(rep(0,9), nrow=3, byrow=TRUE), diag(3))
  )
sigma = matrix(c(5,2,1,2,2,2,1,2,8), nrow=3) 
case3.invWishartList = c(
  invert(new("wishart", df=10, covariance=diag(3)[c(1,3),] %*% sigma %*% t(diag(3)[c(1,3),]))),
  invert(new("wishart", df=10, covariance=sigma)),
  invert(new("wishart", df=25, covariance=sigma))
)

checkApproximation(case3.invWishartList, scaleMatrixList=case3.scaleMatrixList, 
                   xlim=c(-0.25,0.25), ylim=c(0,50), method="trace")
