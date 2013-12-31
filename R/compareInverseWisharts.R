#
# Functions to compare an approximate inverse Wishart to an empirical
# sum of inverse Wishart matrices
#
#

#
# Simulate the sum of inverse Wishart matrices
# 
simulate.sumInvWisharts = function(invWishartList, replicates=1000) {
  
  # get the Wisharts corresponding to the inverse Wisharts
  wishartList = lapply(invWishartList, invert)
  
  # simulate the Wishart replicates
  wishartReplicates = lapply(wishartList, function(obj, replicates) {
    replicateArray = rWishart(replicates, obj@df, obj@covariance)
    # convert the weird array from rWishart to a list of matrices
    wishartReplicates = list()
    for(i in 1:replicates) {
      wishartReplicates[i] = as.matrix(replicateArray[,,i])
    }
    return(wishartReplicates)
  }, replicates)
  
  # invert the Wisharts
  invWishartReplicates = lapply(wishartReplicates, function(obj) {
    return lapply(obj, invert)
  })
  
  # sum them and return the result
  
}

##### Functions for visually comparing densities #####

#
# Plot the density of a single cell of the empirical inverse Wishart
# sum and the approximating inverse Wishart density
#
plotWishartElement = function(empirical, approximate, row, column, 
                              col=c("black", "black"),
                              lty=c(3,1)) {
  empElt = sapply(1:length(empirical),function(i) {return (empirical[[i]][row,column])})
  approxElt = sapply(1:length(approximate),function(i) {return (approximate[[i]][row,column])})
  plot(density(empirical), main="", ylim=c(0,20),
       xlab="",ylab="",xaxt='n',yaxt='n',
       col=col[1], lty=lty[1])
  lines(density(wSum11), col=col[2], lty=lty[2])
  
}

#
# Simulate the specified number of replicates for the
# sum of inverse Wishart matrices and the approximating inverse Wishart.
# Then plot the result
#
compare.plot = function(invWishartList, approximatingInvWishart, replicates=1000,
                        filename=NULL) {
  
}
png(file="invWishartPrelimTest.png",
    width=1000, height=1000, pointsize=18)
par(mfrow=c(3,3),mar=c(0,0,0,0), oma=c(4,1,1,1))
plotWishartElement(1,1)
plotWishartElement(1,2)
plotWishartElement(1,3)
plotWishartElement(2,1)
plotWishartElement(2,2)
par(xpd=NA)
legend(0,-22,c("Empirical", "Approximate"),lty=c(1,1),
       lwd=c(3,4),horiz=TRUE,cex=2,
       col=c("purple","green"))
par(xpd=FALSE)
plotWishartElement(2,3)
plotWishartElement(3,1)
plotWishartElement(3,2)
plotWishartElement(3,3)


dev.off()



##### Calculate the Kullback-Leibler divergence #####