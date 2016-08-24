#' A function for finding local minima in the empirical distribution of a vector
#' of values within a certain range
#'
#'
#' @param data A vector containing numeric data points
#' @param n Number of local minima to find in the empirical distribution
#' @param range Range in data within which to find minima
#' @param plot Plot density and found minima?
#'
#' @return find_minima Return a vector containing n local minima ordered by
#' decreasing density.
#'
#' @examples
#'
#' @export
#'
#'

find_minima_in_range <- function(data, n=1, range=c(min(data), max(data)), plot=TRUE){
  minima <- NULL
  density <- density(data)

  density_range <- density$x > range[1] & density$x < range[2]

  y <- density$y[density_range]
  spacer <- length(which(density$x <= range[1]))

  for ( i in 2:(length(y)-1) ){
    if ( (y[i] < y[i-1]) & (y[i] < y[i+1]) ) {
      minima <- c(minima,i + spacer)
    }
  }
  if ( length(minima) == 0 ) {
    stop('Distribution appears monotonic')
  }

  if ( length(minima) < n ) {
    stop(paste('Only',length(minima), "minima found!", sep=" "))
  }

  positions <- density$x[minima]
  densities <- density$y[minima]

  order <- order(densities, decreasing=FALSE)

  selected_positions <- positions[order[1:n]]

  if(plot){
    hist(data, breaks="FD", freq=FALSE, main=paste("First",n, "local minima", sep= " "), col="lightgrey")
    lines(density, lwd=2)
    l_ply(selected_positions, function(x){abline(v=x, col="red", lwd=2)})
  }

  return(selected_positions)
}
