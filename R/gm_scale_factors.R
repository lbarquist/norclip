#' A function for determining RNA-seq scale factors using the geometric mean
#' across samples, as introduced by Anders and Huber, Genome Biology, 2010.
#'
#'
#' @param matrix A matrix containing columns corresponding to samples, and rows
#' to per feature read counts
#'
#' @return gm_scale_factor Return an array of scale factors.
#'
#' @examples
#'
#' @export
#'
#'

gm_scale_factors <- function(matrix, locfunc=median){
  log_mean <- rowMeans(log(matrix))
  nz <- which(log_mean != 0 & !is.infinite(log_mean))
  SF <- aaply(matrix[nz,],2,function(x){exp(locfunc(log(x) - log_mean[nz]))})

  return(SF)
}
