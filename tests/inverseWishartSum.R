#
# Test cases of a sums of inverse Wishart matrices
#
#
set.seed(2014)

# degrees of freedom scale factors
scaleList = 1:5

#
# Test Case 1: inverse ChiSquare
#
case1.invChiSqList = c(
  invert(new("wishart", df=6, covariance=matrix(c(1)))),
  invert(new("wishart", df=8, covariance=matrix(c(1)))),
  invert(new("wishart", df=9, covariance=matrix(c(1))))
)

# build energy distance results
energyTable = data.frame(
  dfScale = scaleList,
  chiSquareEnergy = sapply(scaleList, function(scale) {
    # scale the degrees of freedom
    invWishartList = sapply(case1.invChiSqList, function(invWishart, scale) {
      return (new("inverseWishart", df=(invWishart@df * scale), 
                  precision=invWishart@precision))
    }, scale)
    # get the approximating inverse Wishart
    approxInvWishart = approximateInverseWishart(invWishartList, scaleMatrixList,
                                                 method=method)
    # get replicates of the sum of the specified inverse Wisharts
    sumReplicates = simulate.sumInvWisharts(invWishartList, replicates=replicates)
    # get replicates of the approximating inv Wishart
    approxReplicates = rInverseWishart(approxInvWishart, n=replicates)
    
    return(calculateEnergyDistance(sumReplicates, approxReplicates))
  })
)

calculateEnergyDistance(sumReplicates, approxReplicates)
checkApproximation(case1.invChiSqList, ylim=c(0,10), xlim=c(-2,4), method="trace")


# inverse Chi square

# Wishart Test Case 1
# Specify our list of inverse Wisharts to sum
#
case1.invWishartList = c(
  invert(new("wishart", df=100, covariance=matrix(c(1,0.3,0.3,0.3,1,0.3,0.3,0.3,1), nrow=3))),
  invert(new("wishart", df=80, covariance=matrix(c(1,0.3,0.3,0.4,3,0.3,0.3,0.3,7), nrow=3))),
  invert(new("wishart", df=60, covariance=matrix(c(5,2,1,2,2,2,1,2,8), nrow=3)))
  )

checkApproximation(case1.invWishartList, ylim=c(0,4), xlim=c(-2,4), method="trace")

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
checkApproximation(case2.invWishartList, xlim=case2.xlim, ylim=case2.ylim, method="trace")


sumReps = simulate.sumInvWisharts(case2.invWishartList, replicates=10000)
approxIW = approximateInverseWishart(case2.invWishartList)
approxReps = rInverseWishart(approxIW, n=10000)

compareCovariance = function(empiricalReps, approxReps, cell1=c(1,2), cell2=c(1,3),
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
  
  if (style == "image") {
    image(empirical.kde, col = colorRampPalette(c("blue", "green", "yellow", "white"))(256))
    image(approx.kde, col = colorRampPalette(c("blue", "green", "yellow", "white"))(256))
  } else if (style == "persp") {
    persp(empirical.kde, phi = 45, theta = 30, border=col[1], zlim=c(0,5))
    persp(approx.kde, phi = 45, theta = 30, border=col[2], zlim=c(0,5))
  } else {
    contour(empirical.kde, col=col[1])
    contour(approx.kde, col=col[2])
  }

  
  #image(empirical.kde)
 # persp(bivn.kde, phi = 45, theta = 30)
  
}

covarianceDiff = function(empiricalReps, approxReps, cell1=c(1,2), cell2=c(1,3),
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
  
  #diff.kde <- kde2d(repDiffCell1, repDiffCell2, 
  #  lims=lims, n=100)
  
  if (style == "image") {
    image(diff.kde, col = colorRampPalette(c("blue", "green", "yellow", "white"))(256))
  } else if (style == "persp") {
    persp(empirical.kde, phi = 45, theta = 30, border=col[1], zlim=c(0,20))
    persp(approx.kde, phi = 45, theta = 30, border=col[2], zlim=c(0,20))
    persp(diff.kde, phi = 45, theta = 30, border="green", zlim=c(0,20))
  } else {
    contour(diff.kde, col=col[1])
  }
  
  
  #image(empirical.kde)
  # persp(bivn.kde, phi = 45, theta = 30)
  
}

# cells = (c(1,1))
# combos = combn(1:9, 2)
# 
# par(mfrow=c(1,3))
# compareCovariance(sumReps, approxReps, cell1=c(1,2), cell2=c(1,3), style="persp")
# compareCovariance(sumReps, approxReps, cell1=c(2,2), cell2=c(2,1), style="image",
#                   lims=c(-2,3,-2,2))
# compareCovariance(sumReps, approxReps, cell1=c(3,3), cell2=c(2,2), style="image",
#                   lims=c(-1,2,0,3))
# compareCovariance(sumReps, approxReps, cell1=c(1,3), cell2=c(3,1), style="image",
#                   lims=c(-2,3,-2,2))
# 
# covarianceDiff(sumReps, approxReps, cell1=c(1,3), cell2=c(3,1), style="persp",
#                lims=c(-2,3,-2,2))

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
sigma = matrix(c(1,0.3,0.4,0.3,1,0.3,0.4,0.3,1), nrow=3, byrow=3) 
case3.invWishartList = c(
  invert(new("wishart", df=10, covariance=diag(3)[c(1,3),] %*% sigma %*% t(diag(3)[c(1,3),]))),
  invert(new("wishart", df=10, covariance=sigma)),
  invert(new("wishart", df=20, covariance=sigma))
)

checkApproximation(case3.invWishartList, scaleMatrixList=case3.scaleMatrixList, 
                   xlim=c(-0.25,0.75), ylim=c(0,30), method="trace")
