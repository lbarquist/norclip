#' Load CLIP-seq data for analysis with \code{norclip}
#'
#' @param data_table A data frame describing the experiments.
#'
#' @return Returns a list of IRanges rle coverage vectors.
#'
#' @details This function loads data to be used for analysis in various
#' functions provided by the \code{norclip} package. The format of the input
#' data frame should be \code{identifier}, \code{type}, \code{direction},
#' \code{file}. See the vignette for details.
#'
#' @examples
#'
#' @export

loadData <- function(data_table){
  stopifnot(is.data.frame(data_table), ncol(data_table) == 4)
  colnames(data_table) <- c("identifier","type","direction","file")

  #sanity check
  plyr::l_ply(unique(data_table$identifier), function(x){
    tmp <- data_table[data_table$identifier == x,]
    if(length(tmp[, "identifier"]) != 4){
      stop("Number of entries for identifer ",x," not equal to four.")
    }
    if(length(tmp[tmp$type == "E" & tmp$direction == "F",1]) != 1){
      stop("Incorrect number of forward experimental files for ", x)
    }
    if(length(tmp[tmp$type == "E" & tmp$direction == "R",1]) != 1){
      stop("Incorrect number of reverse experimental files for ", x)
    }
    if(length(tmp[tmp$type == "C" & tmp$direction == "F",1]) != 1){
      stop("Incorrect number of forward control files for ", x)
    }
    if(length(tmp[tmp$type == "C" & tmp$direction == "R",1]) != 1){
      stop("Incorrect number of reverse control files for ", x)
    }
  })

  message("Assuming ", length(unique(data_table$identifier)), " experiments")

  message("Reading wiggle files")

  wigs <- plyr::llply(data_table$file, function(x){
    this_wig <- import(as.character(x))
    this_wig <- abs(coverage(this_wig, weight="score"))
    return(this_wig)
  }, .progress="text")

  #check replicons identical, get lengths
  replicons <- names(wigs[[1]])
  plyr::l_ply(wigs, function(x){
    if(length(setdiff(replicons, names(x))) > 0){
      stop("Replicons differ between experiments")
    }
  })

  #implement sampling strategy for large genomes?
  max_len <- as.list(rep(0, length(replicons)))
  max_len <- setNames(max_len, replicons)
  for (name in replicons){
    for (i in 1:length(wigs)){
      if(length(wigs[[i]][[name]]) > max_len[[name]]){
        max_len[[name]] <- length(wigs[[i]][[name]])
      }
    }
  }

  #adjust lengths
  for (name in replicons) {
    for (i in 1:length(wigs)){
      diff <- max_len[[name]] - length(wigs[[i]][[name]])
      if(diff > 0){
        wigs[[i]][[name]] <- c(wigs[[i]][[name]], rep(0, times=diff))
      }
    }
  }

  #convert all replicons to single RLE
  wigs <- plyr::llply(wigs, function(x){
    return(unlist(x[order(replicons)]))
  })

  return(wigs)
}
