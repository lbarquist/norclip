#' Construct elliptical filters for matched CLIP-seq experiments
#'
#'
#' @param exp A (possibly IRanges rle) vector containing read-counts per
#' nucleotide in a crosslinked sample
#' @param ctrl A (possibly IRanges rle) vector containing read-counts per
#' nucleotide in a non-crosslinked control
#' @param sdn The number of standard deviations to use in constructing the
#' elliptical filter
#'
#' @return filter_elliptical Return a (possibly IRanges rle) elliptical filter
#' mask.
#'
#'
#' @examples
#'
#' @export
#'
#'

filter_elliptical <- function(exp, ctrl, sdn=8){
  a <- (mean(exp[exp > 0]) + sdn*sd(exp[exp > 0]))^2
  b <- (mean(ctrl[ctrl > 0]) + sdn*sd(ctrl[ctrl > 0]))^2
  return(((exp)^2)/a + ((ctrl)^2)/b > 1)
}
