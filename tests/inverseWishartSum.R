#
# Test cases of a sums of inverse Wishart matrices
#
#
source("../R/wishart.R")
source("../R/sumOfScaledInverseWisharts.R")
source("../R/compareInverseWisharts.R")

#
# Evaluate the accuracy of the approximation
# for the input inverse Wishart list
#
# Calculates the approximating inverse Wishart for the
# list of Wishart matrices, generates replicates for each
# scale factor, calculates the energy distance between the empirical
# and approximate distributions, and produces plots for scale=4
#
# RETURNS:
#  list of energy distances by scale factor
#
evaluateApproximation <- function(invWishartList, scaleMatrixList,
                                  ylim=NULL, xlim=NULL, 
                                  legendX=0, legendY=1, legendCex=1,
                                  replicates=1000,
                                  scaleList = c(1),
                                  plotScaleValue=1,
                                  plotFilename=NULL,
                                  rDataFilename=NULL) {
  # containing for energy distance values
  case.edist = rep(0,length(scaleList))
  
  for(scaleIdx in 1:length(scaleList)) {
    print(paste(c("Scale: ", scaleList[scaleIdx]), collapse=""))
    # scale the degrees of freedom
    scaledInvWishartList = lapply(invWishartList, function(invWishart, scale) {
      return (new("inverseWishart", df=(invWishart@df * scale), 
                  precision=invWishart@precision))
    }, scale=scaleList[scaleIdx])
    class(scaledInvWishartList) = "inverseWishart"
    
    # calculate the approximating inverse Wishart
    approxInvWishart = approximateInverseWishart(scaledInvWishartList, 
                                                 scaleMatrixList,
                                                 method="trace")
    
    # get replicates of the sum of the specified inverse Wisharts
    if (!is.null(scaleMatrixList)) {
      # get replicates of the sum of the specified inverse Wisharts
      sumReplicates = simulate.sumInvWishartsScaled(scaledInvWishartList, 
                                                    scaleMatrixList,
                                                    replicates=replicates)
    } else {
      # get replicates of the sum of the specified inverse Wisharts
      sumReplicates = simulate.sumInvWisharts(scaledInvWishartList, 
                                              replicates=replicates)
    } 
    
    # get replicates of the approximating inv Wishart
    approxReplicates = rInverseWishart(approxInvWishart, n=replicates)
    
    if (rDataFilename != NULL) {
      save(sumReplicates, approxReplicates, 
           filename=paste(c(rDataFilename, "_scale", scaleList[scaleIdx], ".Rdata"), collapse=""))
    }
    
    # calculate the energy distance
    case.edist[scaleIdx] = calculateEnergyDistance(sumReplicates, approxReplicates)
    
    # create a plot for scale factor 4
    if (scaleList[scaleIdx] == plotScaleValue) {
      compare.plot(sumReplicates, approxReplicates, 
                   filename=plotFilename, 
                   ylim=ylim, xlim=xlim,
                   legendX=legendX, legendY=legendY,
                   legendCex=legendCex,
                   col=c("red", "black"))
      
    }
  }
  
  return(case.edist)
  
}


##### Define Test Cases #####
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

##### Evaluate the accuracy of each test case #####
set.seed(1066)
scaleList = 2^(0:5)
replicates=1000

# Check test case 1: inverse ChiSquare variables
case1.edist = evaluateApproximation(case1.invChiSqList, NULL,
                                  ylim=c(0,30), xlim=c(0,0.25), 
                                  legendX=0.05, legendY=-1.6, legendCex=1.5,
                                  replicates=replicates,
                                  scaleList = scaleList,
                                  plotScaleValue=4,
                                  plotFilename="../figures/case1InvChiSquareElementDensity.png")

# Check test case 2: inverse Wishart matrices
case2.edist =  evaluateApproximation(case2.invWishartList, NULL,
                                     ylim=c(0,80), xlim=c(-0.1, 0.2), 
                                     legendX=-0.5, legendY=-3, legendCex=1.5,
                                     replicates=replicates,
                                     scaleList = scaleList,
                                     plotScaleValue=4,
                                     plotFilename="../figures/case2InvWishartElementDensity.png")

# Check test case 3: singular inverse Wisharts
case3.edist =  evaluateApproximation(case3.invWishartList, case3.scaleMatrixList,
                                     ylim=c(0,120), xlim=c(-0.1,0.25), 
                                     legendX=-0.5, legendY=-2, legendCex=1.5,
                                     replicates=replicates,
                                     scaleList = scaleList,
                                     plotScaleValue=4,
                                     plotFilename="../figures/case3SingularInvWishartElementDensity.png")

#
# Build a data frame with the energy distance results
#
edist = data.frame(dfScale=scaleList,
                   chiSQEdist=case1.edist,
                   invWishartEdist=case2.edist,
                   singularInvWishartEdist=case3.edist)
# write the energy table to disk
write.csv(edist, "../data/energyDistance.csv", row.names=FALSE)



