# function to calculate the adjusted rand statistic
# x and y are vectors containing the two partitions to be compared
randError <- function(x, y) {
  # Get crosstabs
  ctab <- table(x,y);
  
  # Calculate 4 intermediary sums
  cellsum <- sum(ctab*(ctab-1)/2)
  totsum <- sum(ctab)*(sum(ctab)-1)/2
  
  # Use matrix multiplication to get row and column marginal sums
  rows <- ctab %*% rep(1,ncol(ctab))
  rowsum <- sum(rows*(rows-1)/2)
  cols <- rep(1,nrow(ctab)) %*% ctab
  colsum <- sum(cols*(cols-1)/2)
  
  # Put them together
  adj.rand <- (cellsum - (rowsum*colsum/totsum))/(.5*(rowsum +colsum)-(rowsum*colsum/totsum))
  return (adj.rand)
}

