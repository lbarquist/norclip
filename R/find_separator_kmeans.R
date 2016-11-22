#' Find separation point between two clusters in a 1D vector
#'
#'
#'
#' @param data A 1-dimensional vector
#' @param plot Plot histogram?
#'
#' @return Returns value for separation point between clusters
#'
#' @details This function computes scale factors for CLIP-seq data, under the
#' assumption that it follows a bi- or multi-modal distribution. For each
#' matched experiment, low count reads are first filtered using an elliptical
#' filter based on the read count standard deviation. Log2 ratios of read counts
#' per nucleotide are calculated, and a local minima is identified between the
#' two highest density maxima. Positions with a log2 ratio less than this minima
#' are then used to fit a scale factor for each experiment. Optionally using the
#' \code{crossnormalize} parameter, all libraries can be scaled based on scale
#' factors computed between control libraries to produce consistent values for
#' hypothesis testing. Note positions with less than \code{bg_cut} reads in
#' any library will be excluded.
#'
#'
#' @examples
#'
#' @seealso \code{\link{norclip}}, \code{\link{clipScaleFactors}}
#'
#' @import Ckmeans.1d.dp
#' @export
#'
find_separator_kmeans <- function(data, plot=TRUE){
  clustering <- Ckmeans.1d.dp(data, 2)
  foreground <- which(clustering$centers == max(clustering$centers))
  fg_indices <- which(clustering$cluster == foreground)
  seperator <- min(data[fg_indices])

  if(plot){
    hist(data, breaks="FD", main=paste("K-means clustering breakpoint"), col="lightgrey")
    abline(v=seperator, col="red")
  }

  return(seperator)
}
