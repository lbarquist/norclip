#' A function for importing a wig file and transforming it into a vector
#' containing per nucleotide counts.
#'
#'
#' @param forward_path Path to forward strand wig file.
#' @param reverse_path Path to reverse strand wig file.
#'
#' @return wig2rle Return an IRanges rle coverage vector.
#'
#' @examples
#'
#' @export
#'
#'

wig2rle <- function(forward_path, reverse_path){
  x <- import(forward_path)
  x <- coverage(x, weight="score")

  y <- import(reverse_path)
  y <- coverage(y, weight="score")

  r <- unlist(c(x,  abs(y)), use.names=F)
  return(unlist(r))
}
