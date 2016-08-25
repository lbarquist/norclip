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
#' @export

naiveScaleFactors <- function(wigs, data_table, sdn=8, bg_cut=5, plot=F){
  colnames(data_table) <- c("identifier","type","direction","file")
  uids <- as.vector(unique(data_table$identifier), mode="list")

  vecs <- llply(uids, function(this_id){
    efi <- which(data_table$identifier == this_id & data_table$type == "E" &
                   data_table$direction =="F")
    eri <- which(data_table$identifier == this_id & data_table$type == "E" &
                   data_table$direction =="R")
    erle <- c(wigs[[efi]], wigs[[eri]])

    cfi <- which(data_table$identifier == this_id & data_table$type == "C" &
                   data_table$direction =="F")
    cri <- which(data_table$identifier == this_id & data_table$type == "C" &
                   data_table$direction =="R")
    crle <- c(wigs[[cfi]], wigs[[cri]])
    return(list(exp=erle,ctrl=crle))
  })

  names(vecs) <- unlist(uids)

  #return(vecs)

  nz_pos_exp <- llply(vecs, function(this_uid){
    return(which(this_uid$exp > bg_cut))
  })

  nz_exp <- Reduce(intersect, nz_pos_exp)
  rm(nz_pos_exp)

  nz_pos_ctrl <- llply(vecs, function(this_uid){
    return(which(this_uid$ctrl > bg_cut))
  })

  nz_ctrl <- Reduce(intersect, nz_pos_ctrl)
  rm(nz_pos_ctrl)

  nz <- union(nz_exp, nz_ctrl)

  message(paste(length(nz)," background positions used for naive normalization", sep=""))

  vecs <- unlist(vecs)

  nz_ar <- laply(vecs, function(this_rle){
    return(as.vector(this_rle[nz], mode="integer"))
  })

  sfs <- gm_scale_factors(t(nz_ar))
  sfs <- matrix(sfs, ncol=2, byrow=T)
  rownames(sfs) <- unlist(uids)
  colnames(sfs) <- c("E", "C")

  if(plot){
    l_ply(uids, function(this_id){
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

      message(paste(this_id,sfs[this_id,"E"], sfs[this_id,"C"]))
      clip_scat(erle, crle, sf_exp=sfs[this_id,"E"], sf_ctrl=sfs[this_id, "C"], xyline=T, elliptical=sdn, main=paste("Naive 2D density plot for", this_id))
    })
  }

  return(sfs)
}
