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

clipScaleFactors <- function(wigs, data_table, sdn=8, crossnormalize=T, plot=T, bg_cut=5){
  colnames(data_table) <- c("identifier","type","direction","file")
  uids <- as.vector(unique(data_table$identifier), mode="list")

  message("Calculating size factors")

  fg_sfs <- laply(uids, function(this_id){
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

    filt <- filter_elliptical(erle, crle, sdn=sdn)

    evec <- as.vector(erle[filt], mode="integer")
    cvec <- as.vector(crle[filt], mode="integer")
    rm(erle, crle)

    ratio <- log2(evec / cvec)

    maxima <- find_maxima(ratio, plot=plot)

    minima <- find_minima_in_range(ratio, range=sort(maxima), plot=plot)

    scale_indices <- which(ratio < minima)

    sf <- median(evec[scale_indices] / cvec[scale_indices])

    if(plot){
      nf <- median(evec / cvec)
      plot(sort(evec[scale_indices] / cvec[scale_indices]), ylab="XL+ to XL- ratio", xlab="Sort Index", pch=20)
      abline(h=nf, col="red")
      abline(h=sf, col="blue")
      mtext(paste("naive scale factor:", nf), side=3, adj=0, line=0, col="red")
      mtext(paste("norclip scale factor:",sf), side=3, adj=0, line=1, col="blue")
      mtext(this_id, side=3, adj=0, line=2)
    }

    return(sf)
  }, .progress="text")

  names(fg_sfs) <- unlist(uids)
  #return(fg_sfs)



  if(crossnormalize){
    bg_vecs <- llply(uids, function(this_id){
      cfi <- which(data_table$identifier == this_id & data_table$type == "C" &
                     data_table$direction =="F")
      cri <- which(data_table$identifier == this_id & data_table$type == "C" &
                     data_table$direction =="R")
      crle <- c(wigs[[cfi]], wigs[[cri]])
      return(crle)
    })

    #return(bg_vecs)

    nz_pos <- llply(bg_vecs, function(this_rle){
      return(which(this_rle > bg_cut))
    })

    nz <- Reduce(intersect, nz_pos)
    rm(nz_pos)
    message(paste(length(nz)," background positions used for crossnormalization", sep=""))

    nz_ar <- laply(bg_vecs, function(this_rle){
      return(as.vector(this_rle[nz], mode="integer"))
    })

    bg_sfs <- gm_scale_factors(t(nz_ar))
    names(bg_sfs) <- uids
  } else {
    bg_sfs <- rep(1, times=length(uids))
  }

  names(bg_sfs) <- unlist(uids)

  sfs <- aaply(seq(from=1,to=length(unlist(uids))), 1, function(this_ind){
    return(cbind(bg_sfs[this_ind] * fg_sfs[this_ind], bg_sfs[this_ind]))
  })
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
      clip_scat(erle, crle, sf_exp=sfs[this_id,"E"], sf_ctrl=sfs[this_id, "C"], xyline=T, elliptical=sdn, main=paste("2D density plot for", this_id))
    })
  }

  return(sfs)
}
