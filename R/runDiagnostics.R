#' A function to generate diagnostic plots for CLIP-seq data.
#'
#' This function produces various diagnostic plots to help the user understand
#' if their data is suitable for \code{\link{norclip}} analysis, and if the
#' default filtering values are appropriate. These include unfiltered and
#' filtered 2D density plots comparing crosslinked libraries to uncrosslinked
#' control libraries, and 1D histograms of the result log2 ratios. See vignette
#' for advice on interpreting these plots.
#'
#'
#' @param wigs A list of IRanges rle coverage vectors.
#' @param data_table Vector metadata in data frame format, see package vignette
#' for details.
#' @param sdn Number of standard deviations for elliptical filtering, see
#' \code{\link{filter_elliptical}} for details.
#' @param colramp Color scheme for 2D density plots. Defaults to mimic the color
#' scheme of
#'
#' @return No explicit return, will produce 3 plots per experiment group:
#' an unfiltered 2D density plot, an elliptically filtered 2D density plot,
#' and a density plot for the log2 ratio of experiment to control library read
#' counts.
#'
#' @examples
#'
#' @seealso \code{\link{norclip}}, \code{\link{loadData}},
#' \code{\link{filter_elliptical}}, \code{\link{clip_scat}}

runDiagnostics <- function(wigs, data_table, sdn=8, colramp=rev(rainbow(10, end = 4/6))){
  message("Running CLIP diagnostics")
  #data_table <-read.table(file, stringsAsFactors = F)
  colnames(data_table) <- c("identifier","type","direction","file")
  #wigs <- loadData(data_table)

  uids <- as.vector(unique(data_table$identifier), mode="list")

  factors <- l_ply(uids, function(this_id){
    efi <- which(data_table$identifier == this_id & data_table$type == "E" &
                   data_table$direction =="F")
    eri <- which(data_table$identifier == this_id & data_table$type == "E" &
                   data_table$direction =="R")
    cfi <- which(data_table$identifier == this_id & data_table$type == "C" &
                   data_table$direction =="F")
    cri <- which(data_table$identifier == this_id & data_table$type == "C" &
                   data_table$direction =="R")

    erle <- c(wigs[[efi]], wigs[[eri]])
    crle <- c(wigs[[cfi]], wigs[[cri]])

    #filt <- filter_elliptical(erle, crle, sdn=0)
    #evec <- as.vector(erle[filt], mode="integer")
    #cvec <- as.vector(crle[filt], mode="integer")



#    smoothScatter(evec, cvec, nbin=5000, nrpoints=0,colramp=colramp, pch=19, cex=.4, xlab="XL", ylab="NXL", main=this_id)

    filt_0 <- which(erle > 0 | crle > 0)

    clip_scat(as.vector(erle[filt_0], mode="integer"), as.vector(crle[filt_0], mode="integer"), elliptical=0, main=paste("Unfiltered 2D density plot for", this_id))

    filt <- filter_elliptical(erle, crle, sdn=sdn)

    evec <- as.vector(erle[filt], mode="integer")
    cvec <- as.vector(crle[filt], mode="integer")
    rm(erle, crle)

    clip_scat(evec, cvec, elliptical = 0, main=paste("Filtered 2D density plot for", this_id))

    ratio <- log2(evec/cvec)

    hist(ratio, breaks="FD", freq=FALSE, main=paste("Ratio density for ", this_id, sep= " "), col="lightgrey")
  }, .progress="text")



}
