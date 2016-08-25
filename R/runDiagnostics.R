#' A function
#'
#'
#' @param forward_path .
#' @param reverse_path Path to reverse strand wig file.
#'
#' @return loadData Return a list of IRanges rle coverage vectors.
#'
#' @examples
#'
#'

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
