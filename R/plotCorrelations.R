#' Plot clustered correlation heatmap for a list of RLE coverage plots
#'
#'
#'
#' @param wigs A list of IRanges rle coverage vectors.
#' @param data_table Vector metadata in data frame format, see package vignette
#' for details.
#' @param method Method for cor function, defaults to "pearson"
#' @param log Log transform read counts?
#' @param heatscale Color range to use for heatscale.
#'
#' @return Returns a correlation matrix.
#'
#' @details This function calculates correlation values between a set of
#' (possibly IRanges RLE) coverage plots, plotting them using the
#' \code{\link{ggheat}} function. Additionally, the correlation matrix is
#' returned.
#'
#'
#' @examples
#'
#' @seealso \code{\link{ggheat}}, \code{\link{runDiagnostics}}
#'
#' @export


plotCorrelations <- function(wigs, data_table, method="pearson", log=F,
                             heatscale= c(low='darkblue',high='lightblue')){
  colnames(data_table) <- c("identifier","type","direction","file")
  nz_pos <- plyr::llply(wigs, function(this_rle){return(which(this_rle != 0))})
  nz <- Reduce(union, nz_pos)

  vec_list <- plyr::llply(wigs, function(this_rle){as.vector(this_rle[nz],
                                                             mode="numeric")})
  nz_mat <- do.call(rbind, vec_list)
  these_names <- c()
  for(index in 1:length(data_table[,1])){
    these_names <- c(these_names,
                     paste(data_table$identifier[index], data_table$type[index],
                           data_table$direction[index], sep="."))
  }
  if(log){
    cormat <- cor(t(log(nz_mat, base=2)), method="pearson")
  } else {
    cormat <- cor(t(nz_mat), method="pearson")
  }
  rownames(cormat) <- these_names
  colnames(cormat) <- these_names
  print(ggheat(cormat, clustering="both", labRow=T, labCol=T,
               heatscale=heatscale, border=T))
  return(cormat)
}
