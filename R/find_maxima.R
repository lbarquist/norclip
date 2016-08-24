#' Find local maxima in the empirical distribution of a vector of values
#'
#'
#' @param data A vector containing numeric data points
#' @param n Number of local maxima to find in the empirical distribution
#' @param plot Plot density and maxima?
#'
#' @return find_maxima Returns a vector containing n local maxima ordered by
#' decreasing density.
#'
#' @examples
#'
#' @export
#'
#'

find_maxima <- function(data, n=2, plot=TRUE){
  maxima <- NULL
  density <- density(data)
  y <- density$y
  for ( i in 2:(length(y)-1) ){
    if ( (y[i] > y[i-1]) & (y[i] > y[i+1]) ) {
      maxima <- c(maxima,i)
    }
  }
  if ( length(maxima) == 0 ) {
    stop('Distribution appears monotonic')
  }

  if ( length(maxima) < n ) {
    stop(paste('Only',length(maxima), "maxima found!", sep=" "))
  }

  positions <- density$x[maxima]
  densities <- density$y[maxima]

  order <- order(densities, decreasing=TRUE)

  selected_positions <- positions[order[1:n]]

  if(plot){
    hist(data, breaks="FD", freq=FALSE, main=paste("First",n, "local maxima", sep= " "), col="lightgrey")
    lines(density, lwd=2)
    l_ply(selected_positions, function(x){abline(v=x, col="red", lwd=2)})
  }

  return(selected_positions)
}
