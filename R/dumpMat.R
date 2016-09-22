#' Calculate scale factors for CLIP-seq data
#'
#'
#'
#' @param wigs A list of IRanges rle coverage vectors.
#' @param data_table Vector metadata in data frame format, see package vignette
#' for details.
#' @param sdn Number of standard deviations for elliptical filtering, see
#' \code{\link{filter_elliptical}} for details.
#' @param crossnormalize Crossnormalize the control libraries?
#' @param plot Produce diagnostic plots?
#' @param bg_cut Cut-off for control library normalization.
#'
#' @return Returns an array of scale factors for the provided libraries.
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
#' hypothesis testing. Note positions that with less than \code{bg_cut} reads in
#' any library will be excluded.
#'
#'
#' @examples
#'
#' @seealso \code{\link{norclip}}, \code{\link{loadData}},
#' \code{\link{gm_scale_factors}}, \code{\link{filter_elliptical}}
#' \code{\link{clip_scat}}
#'
#' @export


plotCorrelations <- function(wigs, data_table){
  colnames(data_table) <- c("identifier","type","direction","file")
  nz_pos <- plyr::llply(wigs, function(this_rle){return(which(this_rle != 0))})
  nz <- Reduce(intersect, nz_pos)

  vec_list <- plyr::llply(wigs, function(this_rle){as.vector(this_rle[nz],
                                                             mode="numeric")})
  nz_mat <- do.call(rbind, vec_list)
  these_names <- c()
  for(index in 1:length(data_table[,1])){
    these_names <- c(these_names,
                     paste(data_table$identifier[index], data_table$type[index],
                           data_table$direction[index], sep="."))
  }
  cormat <- cor(t(nz_mat))
  rownames(cormat) <- these_names
  colnames(cormat) <- these_names
  ggheat(cormat, clustering="both", labRow=T, labCol=T)
}
