#
# Test cases of a sums of inverse Wishart matrices
#
#
set.seed(2014)

#
# Plot the approximation and calculate the KL-divergence
#
checkApproximation = function(invWishartList, scaleMatrixList=NULL,
                              ylim, xlim, method="trace",
                              cell=c(1,1), col=c("red", "black")) {
  print(method)

  #
  # get the approximating inverse Wishart
  #
  if (method == "trace") {
    approxInvWishart = approximateInverseWishart(invWishartList, scaleMatrixList, 
                                                 method="trace")
  } else if (method == "logDeterminant") {
    approxInvWishart = approximateInverseWishart(invWishartList, scaleMatrixList,  
                                                 method="logDeterminant")
  } else if (method == "cellVariance") {
    approxInvWishart = approximateInverseWishart(invWishartList, scaleMatrixList,  
                                                 method="cellVariance", cell=cell)
  } else {
    stop("invalid method")
  }

  compare.plot(invWishartList, approxInvWishart, replicates=10000, 
               ylim=ylim, xlim=xlim, col=col)
  
}

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
  invert(new("wishart", df=40, 
             covariance=matrix(c(1, 0.8, 0.5, 0.8, 1, 0.8, 0.5, 0.8, 1), 
                               nrow=3, byrow=TRUE))),
  invert(new("wishart", df=20, 
             covariance=matrix(c(1, 0.2, 0.2, 0.2, 1, 0.2, 0.2, 0.2, 1), 
                               nrow=3, byrow=TRUE)))
  )
checkApproximation(case2.invWishartList, xlim=c(-0.5,0.5), ylim=c(0,40), method="cellVariance")
checkApproximation(case2.invWishartList, xlim=c(-0.5,0.5), ylim=c(0,40), method="trace")
checkApproximation(case2.invWishartList, xlim=c(-0.5,0.5), ylim=c(0,40), method="logDeterminant")

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
  invert(new("wishart", df=8, covariance=sigma)),
  invert(new("wishart", df=6, covariance=sigma))
)

checkApproximation(case3.invWishartList, scaleMatrixList=case3.scaleMatrixList, 
                   xlim=c(-0.5,0.5), ylim=c(0,40), method="trace")
