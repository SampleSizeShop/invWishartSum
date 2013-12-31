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
      wishartReplicates[[i]] = as.matrix(replicateArray[,,i])
    }
    return(wishartReplicates)
  }, replicates)
  
  # invert the Wisharts
  invWishartReplicates = lapply(wishartReplicates, function(obj) {
    return (lapply(obj, solve))
  })
  
  # sum them and return the result
  sumReplicates = lapply(1:replicates, function(i) {
    return (Reduce('+',lapply(invWishartReplicates, '[[', i)))
  })  
  
  return(sumReplicates)
}

#
# Get a sample of inverse wisharts
#
rInverseWishart = function(invWishart, n=1000) {
  if (class(invWishart) != "inverseWishart") {
    stop("input is not an inverse Wishart object"); 
  }
  if (n <= 0) {
    stop("Invalid number of replicates")
  }
  
  # get the corresponding Wishart
  wishart = invert(invWishart)
  # generate replicates from said Wishart
  wishartArray = rWishart(n, wishart@df, wishart@covariance)
  # make a list rather than the array thing that comes back from the 
  # rWishart function
  wishartReplicates = lapply(1:n, function(i) { return(as.matrix(wishartArray[,,i])) })

  # invert those crazy wishart matrices and return
  return(lapply(wishartReplicates, solve))  
}

##### Functions for visually comparing densities #####

#
# Plot the density of a single cell of the empirical inverse Wishart
# sum and the approximating inverse Wishart density
#
plotWishartElement = function(empirical, approximate, row, column, 
                              col=c("black", "black"),
                              lty=c(1,1), ylim=c(0,1), xlim=c(0,1)) {
  empElt = sapply(1:length(empirical),function(i) {return (empirical[[i]][row,column])})
  approxElt = sapply(1:length(approximate),function(i) {return (approximate[[i]][row,column])})
  plot(density(empElt), main="", 
       ylim=ylim,
       xlim=xlim,
       xlab="",ylab="",xaxt='n',yaxt='n',
       col=col[1], lty=lty[1])
  lines(density(approxElt), col=col[2], lty=lty[2])
  
}

#
# Simulate the specified number of replicates for the
# sum of inverse Wishart matrices and the approximating inverse Wishart.
# Then plot the result
#
compare.plot = function(invWishartList, approximatingInvWishart, replicates=1000,
                        filename=NULL, ylim=c(0,1), xlim=c(0,1),
                        col=c("black", "black")) {
  # get replicates of the sum of the specified inverse Wisharts
  sumReplicates = simulate.sumInvWisharts(invWishartList, replicates=replicates)
  
  # get replicates of the approximating inv Wishart
  approxReplicates = rInverseWishart(approximatingInvWishart, n=replicates)
  
  # plot the replicates
  dim = nrow(approximatingInvWishart@precision)
  # write to a file if specified
  if (!is.null(filename)) {
    png(file=filename,
        width=1000, height=1000, pointsize=18)
  }
  #
  # Plot the density of each cell 
  #
  par(mfrow=c(dim,dim),mar=c(0,0,0,0), oma=c(4,1,1,1))
  for(r in 1:dim) {
    for(c in 1:dim) {
      plotWishartElement(
        empirical=sumReplicates, approximate=approxReplicates, 
        row=r, column=c, ylim=ylim, xlim=xlim, col=col)
    }
  }
  # close output device if writing to a file
  if (!is.null(filename)) {
    dev.off()
  }
  
}




##### Calculate the Kullback-Leibler divergence #####