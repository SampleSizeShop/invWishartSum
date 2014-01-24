#
# Test cases of a sums of inverse Wishart matrices
#
#
set.seed(2014)

#
# Test Case 1: sum of inverse ChiSquare variables (i.e. 1x1 inverse Wishart)
#
case1.invChiSqList = c(
  invert(new("wishart", df=6, covariance=matrix(c(1)))),
  invert(new("wishart", df=8, covariance=matrix(c(1)))),
  invert(new("wishart", df=9, covariance=matrix(c(1))))
)

# 
# Test case 2: sum of 3x3 inverse Wishart matrices 
#
case2.invWishartList = c(
  invert(new("wishart", df=6, covariance=matrix(c(1,0.3,0.3,0.3,1,0.3,0.3,0.3,1), nrow=3))),
  invert(new("wishart", df=8, covariance=matrix(c(1,0.3,0.3,0.3,3,0.3,0.3,0.3,7), nrow=3))),
  invert(new("wishart", df=9, covariance=matrix(c(5,2,1,2,2,2,1,2,8), nrow=3)))
)

#
# Test case 3: sum of scaled, zero-padded Wisharts sharing the same base covariance
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
  invert(new("wishart", df=6, covariance=diag(3)[c(1,3),] %*% sigma %*% t(diag(3)[c(1,3),]))),
  invert(new("wishart", df=8, covariance=sigma)),
  invert(new("wishart", df=9, covariance=sigma))
)

# degrees of freedom scale factors
scaleList = 2^(0:5)
# test cases
testCases.dist = list(case1.invChiSqList, case2.invWishartList, case3.invWishartList)
testCases.scale = list(NULL, NULL, case3.scaleMatrixList)

#
# Calculate the energy distance for each level
#
replicates=100
edist = data.frame(scale=scaleList)
for(case in 1:length(testCases.dist)) {
  print(paste(c("Processing case", case), collapse=""))
  # get current case
  case.dist = testCases.dist[[case]]
  case.scale = testCases.scale[[case]]
  # allocate an array to store the energy distances
  case.edist = rep(0,length(scaleList))
  
  for(scaleIdx in 1:length(scaleList)) {
    print(paste(c("Scale: ", scaleList[scaleIdx]), collapse=""))
    # scale the degrees of freedom
    invWishartList = lapply(case.dist, function(invWishart, scale) {
      return (new("inverseWishart", df=(invWishart@df * scale), 
                  precision=invWishart@precision))
    }, scale=scaleList[scaleIdx])
    class(invWishartList) = "inverseWishart"
    
    # calculate the approximating inverse Wishart
    approxInvWishart = approximateInverseWishart(invWishartList, case.scale,
                                                 method="trace")
    
    # get replicates of the sum of the specified inverse Wisharts
    if (!is.null(case.scale)) {
      # get replicates of the sum of the specified inverse Wisharts
      sumReplicates = simulate.sumInvWishartsScaled(invWishartList, case.scale,
                                                    replicates=replicates)
    } else {
      # get replicates of the sum of the specified inverse Wisharts
      sumReplicates = simulate.sumInvWisharts(invWishartList, replicates=replicates)
    } 
    
    # get replicates of the approximating inv Wishart
    approxReplicates = rInverseWishart(approxInvWishart, n=replicates)
    
    case.edist[scaleIdx] = calculateEnergyDistance(sumReplicates, approxReplicates)
    
  }

  edist <- cbind(edist, case.edist)
  names(edist)[case+1] = paste("case",case,"_edist",collapse="",sep="")
  
}





# inverse Chi square



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
