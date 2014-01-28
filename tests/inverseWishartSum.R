#
# Test cases of a sums of inverse Wishart matrices
#
#
source("../R/wishart.R")
source("../R/sumOfScaledInverseWisharts.R")
source("../R/compareInverseWisharts.R")

#
# Generate .Rdata files containing the empirical sample of the sum
# of inverse Wisharts
#
#
generateEmpiricalSumData <- function(invWishartList, scaleMatrixList, dfScaleList,
                                     replicates=1000,
                                     outputDir=".", filenamePrefix="sumReplicates") {
  filenames=rep("", length(dfScaleList))
  idx = 1
  for(scale in dfScaleList) {

    # scale the degrees of freedom
    dfScaleInvWishartList = lapply(invWishartList, function(invWishart, scale) {
      return (new("inverseWishart", df=(invWishart@df * scale), 
                  precision=invWishart@precision))
    }, scale=scale)
    class(dfScaleInvWishartList) = "inverseWishart"
    
    # get replicates of the sum of the specified inverse Wisharts
    if (!is.null(scaleMatrixList)) {
      # get replicates of the sum of the specified inverse Wisharts
      sumReplicates = simulate.sumInvWishartsScaled(dfScaleInvWishartList, 
                                                    scaleMatrixList,
                                                    replicates=replicates)
    } else {
      # get replicates of the sum of the specified inverse Wisharts
      sumReplicates = simulate.sumInvWisharts(dfScaleInvWishartList, 
                                              replicates=replicates)
    }
    
    # save as an .Rdata file
    filename = paste(c(outputDir, filenamePrefix, "_dfScale", scale, ".Rdata"), collapse="")
    save(sumReplicates, file=filename)
    
    filenames[idx] = filename
    idx = idx + 1
  }
  
  return(filenames);
}

#
# Generate .Rdata files containing a sample of the approximating Wishart
#
generateApproximationData <- function(invWishartList, scaleMatrixList, dfScaleList,
                                      replicates=1000,
                                      outputDir=".", filenamePrefix="approxReplicates") {
  filenames=rep("", length(dfScaleList))
  idx = 1
  
  for(scale in dfScaleList) {
    
    # scale the degrees of freedom
    dfScaleInvWishartList = lapply(invWishartList, function(invWishart, scale) {
      return (new("inverseWishart", df=(invWishart@df * scale), 
                  precision=invWishart@precision))
    }, scale=scale)
    class(dfScaleInvWishartList) = "inverseWishart"
    
    # calculate the approximating inverse Wishart
    approxInvWishart = approximateInverseWishart(dfScaleInvWishartList, 
                                                 scaleMatrixList,
                                                 method="trace")
    
    # get replicates of the approximating inv Wishart
    approxReplicates = rInverseWishart(approxInvWishart, n=replicates)
    
    # save as an .Rdata file
    filename = paste(c(outputDir, filenamePrefix, "_dfScale", scale, ".Rdata"), collapse="")
    save(approxReplicates, file=filename)
    
    filenames[idx] = filename
    idx = idx + 1
  }
  
  return(filenames)
}

#
# Calculate the energy distance between samples contained in the specified 
# .Rdata files
#
# Calculates the energy distance between each pair of sum replicates and 
# approximate replicates.  Assumes that the filename lists are in the same
# order to pair the samples correctly
#
getEnergyDistanceForFilelist <- function(sumFileList, approxFileList) {
  if (length(sumFileList) != length(approxFileList)) {
    stop("the number of sum sample files must match the number of approx sample files")
  }
  
  # containing for energy distance values
  edistList = rep(0,length(sumFileList))
  
  for(i in 1:length(sumFileList)) {
    print(paste(c("Calculating energy distance for samples in files [", sumFileList[i], 
          "] and [", approxFileList[i], "]"), collapse=""))
    
    # load the sum sample from the Rdata file - sets a variable called sumReplicates
    load(sumFileList[i])
    # load the approx sample from the Rdata file - sets a variable called approxReplicates
    load(approxFileList[i])
    
    # calculate the energy distance

    edistList[i] = calculateEnergyDistance(sumReplicates, approxReplicates)
    
  }
  
  return(edistList)
  
}

#
# Plot the univariate densities for the sum and the approximation
#
plotDensityFromFile <- function(sumFile, approxFile, filename,
                                ylim=c(0,1), xlim=c(0,1), col=c("red", "black")) {
  print(paste(c("Plotting univariate densities for files [", sumFile, 
                "] and [", approxFile, "]"), collapse=""))
  
  # load the sum sample from the Rdata file - sets a variable called sumReplicates
  load(sumFile)
  # load the approx sample from the Rdata file - sets a variable called approxReplicates
  load(approxFile)
  # plot the densities
  pdf(paste(c("../../text/Figures/", filename), collapse=""))
  compare.plot(sumReplicates, approxReplicates, ylim=ylim, xlim=xlim, col=col)
  dev.off()
}

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

##### Generate samples for each test case and each scale factor #####
set.seed(1812)
outDataDir = "../data/"
dfScaleList = 2^(0:5)
replicates=1000
# test case 1
case1.approxFiles = generateApproximationData(case1.invChiSqList, NULL, dfScaleList=dfScaleList, 
                                             replicates=replicates, filenamePrefix="case1ApproxReplicates", 
                                             outputDir=outDataDir)
case1.sumFiles = generateEmpiricalSumData(case1.invChiSqList, NULL, dfScaleList=dfScaleList, 
                                          replicates=replicates, filenamePrefix="case1SumReplicates", 
                                          outputDir=outDataDir)
# test case 2
case2.approxFiles = generateApproximationData(case2.invWishartList, NULL, dfScaleList=dfScaleList, 
                                              replicates=replicates, filenamePrefix="case2ApproxReplicates", 
                                              outputDir=outDataDir)
case2.sumFiles = generateEmpiricalSumData(case2.invWishartList, NULL, dfScaleList=dfScaleList, 
                                          replicates=replicates, filenamePrefix="case2SumReplicates", 
                                          outputDir=outDataDir)
# test case 3
case3.approxFiles = generateApproximationData(case3.invWishartList, case3.scaleMatrixList, 
                                              dfScaleList=dfScaleList, replicates=replicates, 
                                              filenamePrefix="case3ApproxReplicates", outputDir=outDataDir)
case3.sumFiles = generateEmpiricalSumData(case3.invWishartList, case3.scaleMatrixList, 
                                          dfScaleList=dfScaleList, replicates=replicates,
                                          filenamePrefix="case3SumReplicates", outputDir=outDataDir)


##### Calculate the energy distance for each test case and scale factor #####
case1.edist = getEnergyDistanceForFilelist(case1.sumFiles, case1.approxFiles)
case2.edist = getEnergyDistanceForFilelist(case2.sumFiles, case2.approxFiles)
case3.edist = getEnergyDistanceForFilelist(case3.sumFiles, case3.approxFiles)

# Build a data frame with the energy distance results
edist = data.frame(dfScale=dfScaleList,
                   chiSQEdist=case1.edist,
                   invWishartEdist=case2.edist,
                   singularInvWishartEdist=case3.edist)
# write the energy table to disk
write.csv(edist, "../data/energyDistance.csv", row.names=FALSE)

# plot the energy distance by df scale factor
pdf(file="../../text/Figures/energyDistanceByScaleFactor.pdf")
plot(edist$dfScale, edist$chiSQEdist, "l", cex.lab=1.5, ylim=c(0,10),
     xlab="Degrees of Freedom Scale Factor", ylab="Energy Distance")
lines(edist$dfScale, edist$invWishartEdist, lty=2)
lines(edist$dfScale, edist$singularInvWishartEdist, lty=3)
legend("topright", c("Chi-Square", "Inverse Wishart", "Singular inverse Wishart"), 
       lty = c(1,2,3), cex=1.5)
dev.off()

##### Plot the univariate densities for the sum 
##### and the approximating inverse Wishart 
#####
idx = which(dfScaleList==4)
plotDensityFromFile(case1.sumFiles[idx], case1.approxFiles[idx], 
                    ylim=c(0,40), xlim=c(0,0.25), filename="invChiSqDensity.pdf")
plotDensityFromFile(case2.sumFiles[idx], case2.approxFiles[idx], 
                    ylim=c(0,80), xlim=c(-0.1,0.15), filename="invWishartDensity.pdf")
plotDensityFromFile(case3.sumFiles[idx], case3.approxFiles[idx], 
                    ylim=c(0,100), xlim=c(-0.1,0.15), filename="singularInvWishartDensity.pdf")


##### Plot selected covariances for the p > 1 cases #####

##### Evaluate the accuracy of each test case #####

# 
# # Check test case 1: inverse ChiSquare variables
# case1.edist = evaluateApproximation(case1.invChiSqList, NULL,
#                                   ylim=c(0,30), xlim=c(0,0.25), 
#                                   legendX=0.05, legendY=-1.6, legendCex=1.5,
#                                   replicates=replicates,
#                                   scaleList = scaleList,
#                                   plotScaleValue=4,
#                                   plotFilename="../figures/case1InvChiSquareElementDensity.png")
# 
# # Check test case 2: inverse Wishart matrices
# case2.edist =  evaluateApproximation(case2.invWishartList, NULL,
#                                      ylim=c(0,80), xlim=c(-0.1, 0.2), 
#                                      legendX=-0.5, legendY=-3, legendCex=1.5,
#                                      replicates=replicates,
#                                      scaleList = scaleList,
#                                      plotScaleValue=4,
#                                      plotFilename="../figures/case2InvWishartElementDensity.png")
# 
# # Check test case 3: singular inverse Wisharts
# case3.edist =  evaluateApproximation(case3.invWishartList, case3.scaleMatrixList,
#                                      ylim=c(0,120), xlim=c(-0.1,0.25), 
#                                      legendX=-0.5, legendY=-2, legendCex=1.5,
#                                      replicates=replicates,
#                                      scaleList = scaleList,
#                                      plotScaleValue=4,
#                                      plotFilename="../figures/case3SingularInvWishartElementDensity.png")





